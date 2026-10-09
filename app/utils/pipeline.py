import time
import xml.etree.ElementTree as ET
import numpy as np
import streamlit as st
import pandas as pd
from Bio import Entrez
import urllib.parse
import os
import urllib.parse
import gseapy as gp
import io, zipfile
import re
import urllib.parse

GRAPHQL_URL = "https://dgidb.org/api/graphql"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/"

# Methods
import requests
from pathlib import Path


def download_mesh_xml(save_dir="data/mesh", year=2026):
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    url = f"https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc{year}.xml"
    file_path = Path(save_dir) / f"desc{year}.xml"

    if not file_path.exists():
        print(f"Downloading {file_path.name}...")
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(file_path, "wb") as f:
                for chunk in r.iter_content(1024):
                    f.write(chunk)
            print(f"Downloaded file: {file_path}")
        else:
            raise ConnectionError(f"Could not download: {r.status_code}")
    else:
        print(f"File exists: {file_path}")

    return str(file_path)

def load_mesh_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    records = []
    for descriptor in root.findall(".//DescriptorRecord"):
        mesh_id = descriptor.findtext("DescriptorUI")
        name = descriptor.findtext("DescriptorName/String")
        # Some descriptors have multiple TreeNumbers
        tree_numbers = [tn.text for tn in descriptor.findall(".//TreeNumberList/TreeNumber")]
        if not tree_numbers:
            tree_numbers = np.nan  # For later filtering
        records.append({
            "MeSH_ID": mesh_id,
            "Name": name,
            "TreeNumber": tree_numbers
        })
    return pd.DataFrame(records)

from concurrent.futures import ThreadPoolExecutor

@st.cache_data(show_spinner=False)
def cached_open_targets(gene_tuple):
    genes = list(gene_tuple)
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(find_possible_target_of_drugs, genes))
    return [r for r in results if r is not None]

@st.cache_data(show_spinner=False)
def cached_dgidb(gene_list_tuple):
    return get_drug_targets_dgidb_graphql(list(gene_list_tuple))

@st.cache_data(show_spinner=False)
def get_mesh_df():
    path = download_mesh_xml()
    df = load_mesh_xml(path)
    df = df[df["TreeNumber"].apply(
        lambda x: any(tn.startswith("C") for tn in x) if isinstance(x, list) else False
    )]
    return df

def get_disease_name(mesh_id):
    """
    Returns the corresponding name of a disease given a MESH id provided by the user. 
    """
    try:
        search_handle = Entrez.esearch(db="mesh", term=mesh_id)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        if not search_record['IdList']:
            return None, None
        uid = search_record['IdList'][0]
        summary_handle = Entrez.esummary(db="mesh", id=uid)
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()
        disease_url = f"https://www.ncbi.nlm.nih.gov/mesh/{uid}"
        return summary_record[0]['DS_MeshTerms'][0], disease_url
    except:
        return None, None

def generate_expression_atlas_link(disease_name):
    """
    Generates link to Expression Atlas for a given disease to extract gene information. 
    """
    encoded_disease = urllib.parse.quote(disease_name)
    base_url = (
        "https://www.ebi.ac.uk/gxa/search?geneQuery=%5B%5D"
        "&species=Homo%20sapiens"
        f"&conditionQuery=[%7B%22value%22%3A%22{encoded_disease}%22%7D]"
        "&ds=%7B%22kingdom%22%3A%5B%22animals%22%5D%2C%22regulation%22%3A%5B%22UP%22%5D%7D"
        "&bs=%7B%22homo%20sapiens%22%3A%5B%22ORGANISM_PART%22%5D%7D"
        "#differential"
    )
    return base_url


@st.cache_data(show_spinner=False)
def get_gene_names_batch(ensembl_ids):
    url = "https://rest.ensembl.org/lookup/id"
    headers = {"Content-Type": "application/json"}

    payload = {"ids": list(ensembl_ids)}
    response = requests.post(url, json=payload, headers=headers)
    if response.status_code != 200:
        st.warning(f"Ensembl request failed with status {response.status_code}")
        return {}
    data = response.json()
    gene_map = {}
    for gene_id, info in data.items():
        if info is None:
            gene_map[gene_id] = "Not Found"
        else:
            gene_map[gene_id] = info.get("display_name", "Not Found")

    return gene_map

def fetch_gene_names(df):
    """
    Fetches gene names from ensembl and groups them according to log_2 fold change.
    """
    required_columns = ["Gene", "log_2 fold change", "Adjusted p-value"]
    for col in required_columns:
        if col not in df.columns:
            st.error(f"Missing required column: '{col}'")
            st.stop()

    gene_ids = tuple(
        str(gene_id)
        for gene_id in df["Gene"].dropna().unique()
    )

    gene_map = get_gene_names_batch(gene_ids)

    df["Gene Name"] = df["Gene"].map(gene_map)
    df_filtered = df[df["Gene Name"] != "Not Found"]

    grouped = df_filtered.groupby(
        ["Gene", "Gene Name"], as_index=False
    ).agg({
        "log_2 fold change": "sum",
        "Adjusted p-value": "sum"
    })

    return grouped


def find_possible_target_of_drugs(ensembl_id):
    """
    Given an Ensembl Gene ID, it asks the OpenTragetS API to check if it's a drug target.
    Returns the gene symbol, whether it's a known drug target, and associated approved drugs.
    """
    
    #query string to get general information about AR and genetic constraint and tractability assessments 
    query_string = """
      query target($ensemblId: String!){
        target(ensemblId: $ensemblId){
          id
          approvedSymbol
          biotype
          geneticConstraint {
            constraintType
            exp
            obs
            score
            oe
            oeLower
            oeUpper
          }
          tractability {
            label
            modality
            value
          }
        }
      }"""

    # Set variables object of arguments to be passed to endpoint
    variables = {"ensemblId": ensembl_id}

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Perform POST request and check status code of response
    try:
        r = requests.post(base_url, json={"query": query_string, "variables": variables})
        if r.status_code != 200:
            return None
        data = r.json()['data']['target']
        return {
            "Gene Symbol": data.get("approvedSymbol", ""),
            "Ensembl ID": data.get("id", ""),
            "Name": data.get("approvedName", ""),
            "Biotype": data.get("biotype", ""),
            "Tractability": [
                t["label"] for t in data.get("tractability", []) if t["value"]
            ]
        }
    except:
        return None


def analyze_pathways(df, number, retries=3, delay=5):
    """
    Analyzes n (given number) pathways in which the genes interact.
    Retries Enrichr API call up to 'retries' times in case of timeout.
    """
    # Prepare gene list
    df_genes = df["Gene Name"].dropna().astype(str).str.strip().str.upper().unique().tolist()

    # Enrichment with retry
    for attempt in range(retries):
        try:
            enr = gp.enrichr(
                gene_list=df_genes,
                gene_sets="Reactome_2022",
                organism="hsapiens",
                outdir=None
            )
            if enr.results.empty:
                return None
            top_pathways = enr.results.copy()
            break  # success, exit retry loop
        except gp.enrichr.EnrichrAPIError as e:
            st.warning(f"Enrichr API error (attempt {attempt+1}/{retries}): {e}")
            if attempt < retries - 1:
                st.info(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                st.error("Enrichr server is not responding. Try again later.")
                return None

    # Create Reactome links
    def make_reactome_link(term):
        match = re.search(r'(R-HSA-\d+)', term)
        if match:
            rid = match.group(1)
            return f'<a href="https://reactome.org/PathwayBrowser/#/{rid}" target="_blank">{term}</a>'
        return term

    top_pathways["Reactome Link"] = top_pathways["Term"].apply(make_reactome_link)

    # Parse overlap info
    top_pathways[["Input Genes", "Pathway Genes"]] = (
        top_pathways["Overlap"].str.split("/", expand=True).astype(int)
    )
    top_pathways["Input %"] = top_pathways["Input Genes"] / top_pathways["Pathway Genes"] * 100
    top_pathways = top_pathways.sort_values("Adjusted P-value", ascending=True).head(number)

    # Sum log2fc for overlapping genes per pathway
    df["Gene Name"] = df["Gene Name"].astype(str).str.strip().str.upper()
    sum_fc = []
    for _, row in top_pathways.iterrows():
        genes = [g.strip().upper() for g in row["Genes"].split(";")]
        overlap_df = df[df["Gene Name"].isin(genes)]
        sum_fc.append(overlap_df["log_2 fold change"].sum())
    top_pathways["Sum log2fc"] = sum_fc

    return top_pathways

def get_overlapping_genes(df, selected_pathway_row):
    """
    Get overlapped genes as a way to identify most important genes in a given specific pathway
    """
    # Normalize input genes
    df["Gene Name"] = df["Gene Name"].astype(str).str.strip().str.upper()

    # Get genes from selected pathway
    pathway_genes = [g.strip().upper() for g in selected_pathway_row["Genes"].split(";")]

    # Get overlapping genes from input
    overlap_df = df[df["Gene Name"].isin(pathway_genes)].copy()
    overlap_df["abs_fc"] = overlap_df["log_2 fold change"].abs()
    overlap_df = overlap_df.sort_values("abs_fc", ascending=False)

    return overlap_df



def get_drug_targets_dgidb_graphql(gene_names):

    """
    Queries the DGIdb GraphQL API for drug–gene interaction data for a list of gene names.

    For each gene, the function retrieves associated drugs, interaction types, 
    directionality, interaction scores, sources, and PMIDs (if available), 
    and compiles the results into a pandas DataFrame.
    """

    all_results = []
    for gene_name in gene_names:
        graphql_query = f"""
        {{
          genes(names: ["{gene_name}"]) {{
            nodes {{
              interactions {{
                drug {{
                  name
                  conceptId
                }}
                interactionScore
                interactionTypes {{
                  type
                  directionality
                }}
                interactionAttributes {{
                  name
                  value
                }}
                publications {{
                  pmid
                }}
                sources {{
                  sourceDbName
                }}
              }}
            }}
          }}
        }}
        """
        response = requests.post(GRAPHQL_URL, json={"query": graphql_query})
        if response.status_code == 200:
            try:
                data = response.json()
                if 'data' in data and 'genes' in data['data'] and len(data['data']['genes']['nodes']) > 0:
                    interactions = data['data']['genes']['nodes'][0].get('interactions', [])
                    for interaction in interactions:
                        drug_name = interaction['drug'].get('name', 'Unknown')
                        score = interaction.get('interactionScore', 'N/A')
                        types = interaction.get('interactionTypes', [])
                        interaction_type = types[0].get('type', 'N/A') if types else 'N/A'
                        directionality = types[0].get('directionality', 'N/A') if types else 'N/A'
                        sources = interaction.get('sources', [])
                        source = sources[0]['sourceDbName'] if sources else 'N/A'
                        pmids = interaction.get('publications', [])
                        pmid = pmids[0]['pmid'] if pmids else 'N/A'

                        all_results.append({
                            'Gene': gene_name,
                            'Drug': drug_name,
                            'Interaction Type': interaction_type,
                            'Directionality': directionality,
                            'Source': source,
                            'PMID': pmid,
                            'Interaction Score': score
                        })
            except:
                continue

    return pd.DataFrame(all_results)


def drug_with_links(df):
    """
    Returns the df with links to the resources
    """
    df_with_links = df.copy()
    df_with_links["Gene"] = df_with_links["Gene"].apply(
            lambda gene: f'<a href="https://dgidb.org/results?searchType=gene&searchTerms={urllib.parse.quote(gene)}" target="_blank">{gene}</a>'
            if gene else ""
    )
    df_with_links["Drug"] = df_with_links["Drug"].apply(
        lambda drug: f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(drug)}" target="_blank">{drug}</a>'
    )
    df_with_links["PMID"] = df_with_links["PMID"].apply(
        lambda pmid: f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" target="_blank">{pmid}</a>'
        if pmid else "NaN"
    )
    return df_with_links

def normalize_disease_name(name: str) -> str:
    """
    Normalizes the disease name to match those in Expression Atlas
    """
    # If there's a comma, move the part after the comma to the front
    if "," in name:
        parts = [part.strip() for part in name.split(",")]
        # e.g., ["Muscular Dystrophy", "Duchenne"] → "Duchenne Muscular Dystrophy"
        return f"{parts[1]} {parts[0]}"
    return name


def add_links_to_final_table(df):
    """
    Returns the final table df with links to Ensembl, DGIDB, PubMed and Reactome
    """
    df = df.copy()
    
    # Gene → Ensembl
    if "Ensembl ID" in df.columns:
        df["Ensembl ID"] = df["Ensembl ID"].apply(
            lambda gene_id: f'<a href="https://www.ensembl.org/search/results?query={urllib.parse.quote(str(gene_id))}" target="_blank">{gene_id}</a>'
        )
    
    # Drug → DGIdb
    if "Drug" in df.columns:
        def drug_links(drugs_str):
            if pd.isna(drugs_str):
                return ""
            drugs = [d.strip() for d in drugs_str.split(";")]
            return "; ".join([f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(d)}" target="_blank">{d}</a>' for d in drugs])
        df["Drug"] = df["Drug"].apply(drug_links)
    
    # PMID → PubMed
    if "PMID" in df.columns:
        def pmid_links(pmids_str):
            if pd.isna(pmids_str):
                return ""
            pmids = [p.strip() for p in str(pmids_str).split(";")]
            return "; ".join([f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}" target="_blank">{p}</a>' for p in pmids])
        df["PMID"] = df["PMID"].apply(pmid_links)
    
    # Pathways → Reactome
    if "Pathways" in df.columns:
        def pathway_links(pathways_str):
            if pd.isna(pathways_str):
                return ""
            pathways = [p.strip() for p in pathways_str.split(";")]
            linked_pathways = []
            for p in pathways:
                match = re.search(r"(R-HSA-\d+)", p)  # Extract Reactome ID
                if match:
                    reactome_id = match.group(1)
                    linked_pathways.append(f'<a href="https://reactome.org/content/detail/{reactome_id}" target="_blank">{p}</a>')
                else:
                    linked_pathways.append(p)
            return "; ".join(linked_pathways)
        df["Pathways"] = df["Pathways"].apply(pathway_links)
    
    return df

def save_pathway_csvs(df_selected, top_pathways):
    """
    Creates csv files for each of the top N pathways, with info about genes
    """
    csv_files = {}

    for _, row in top_pathways.iterrows():
        pathway_name = row["Term"]
        pathway_genes = get_overlapping_genes(df_selected, row)

        if not pathway_genes.empty:
            safe_name = pathway_name.replace("/", "_").replace(" ", "_")

            buf = io.StringIO()
            pathway_genes.to_csv(buf, index=False)
            csv_files[f"{safe_name}_genes.csv"] = buf.getvalue().encode("utf-8")

    return csv_files

def save_drug_csvs(df_selected, top_pathways):
    """
    Creates csv files for each of the top N pathways with info about genes-drugs
    """
    csv_files = {}

    for _, row in top_pathways.iterrows():
        pathway_name = row["Term"]
        pathway_genes = get_overlapping_genes(df_selected, row)

        drug_df = get_drug_targets_dgidb_graphql(pathway_genes["Gene Name"].tolist())
        if not drug_df.empty:
            drug_df = drug_with_links(drug_df)
            safe_name = pathway_name.replace("/", "_").replace(" ", "_")

            # Write to memory instead of disk
            buf = io.StringIO()
            drug_df.to_csv(buf, index=False)
            csv_files[f"{safe_name}_drugs.csv"] = buf.getvalue().encode("utf-8")

    return csv_files

@st.cache_data(show_spinner=False)
def get_ensembl_ids_batch(gene_names):
    """Map human gene symbols to Ensembl gene IDs."""

    gene_names = tuple(
        dict.fromkeys(
            str(name).strip()
            for name in gene_names
            if name is not None and str(name).strip()
        )
    )

    if not gene_names:
        return {}

    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json"}

    response = requests.post(
        url,
        json={"symbols": list(gene_names)},
        headers=headers,
        timeout=60,
    )
    response.raise_for_status()

    results = response.json()

    return {
        symbol: results[symbol]["id"]
        for symbol in gene_names
        if results.get(symbol) and results[symbol].get("id")
    }


def prepare_gene_csv(df_raw):
    """Prepare an uploaded gene CSV for the existing DisGenix pipeline."""

    df = df_raw.copy()

    # Ensure the expected columns exist.
    if "Gene Name" not in df.columns:
        raise ValueError("The CSV must contain a 'Gene Name' column.")

    if "log_2 fold change" not in df.columns:
        df["log_2 fold change"] = np.nan

    # Clean gene symbols.
    df["Gene Name"] = df["Gene Name"].astype("string").str.strip()

    df = df.dropna(subset=["Gene Name"]).copy()
    df = df[
        (df["Gene Name"] != "")
        & (df["Gene Name"].str.lower() != "nan")
    ].copy()

    if df.empty:
        st.error("No valid gene names were found in the uploaded CSV.")
        st.stop()

    # Ensure fold changes are numeric, preserving missing values.
    df["log_2 fold change"] = pd.to_numeric(
        df["log_2 fold change"],
        errors="coerce",
    )

    # Map gene symbols to Ensembl IDs.
    gene_names = tuple(df["Gene Name"].dropna().unique())
    gene_mapping = get_ensembl_ids_batch(gene_names)

    df["Gene"] = df["Gene Name"].map(gene_mapping)

    # Report genes that could not be mapped.
    unmapped = df.loc[df["Gene"].isna(), "Gene Name"].unique()

    if len(unmapped) > 0:
        st.warning(
            f"{len(unmapped)} gene(s) could not be mapped to Ensembl IDs "
            "and will be excluded from the downstream analysis: "
            + ", ".join(map(str, unmapped[:20]))
            + (" ..." if len(unmapped) > 20 else "")
        )

    df = df.dropna(subset=["Gene"]).copy()

    if df.empty:
        st.error(
            "None of the uploaded genes could be mapped to Ensembl IDs. "
            "Please check that the CSV contains valid human gene symbols."
        )
        st.stop()

    # Combine duplicate gene symbols, if present.
    df = (
        df.groupby(["Gene", "Gene Name"], as_index=False)
        ["log_2 fold change"]
        .mean()
    )

    # Match the structure expected by the existing pipeline.
    df = df[["Gene", "Gene Name", "log_2 fold change"]]

    return df.reset_index(drop=True)