import plotly.express as px
import datetime
from collections import defaultdict
import streamlit.components.v1 as components
from utils.pipeline import *
from utils.utils import estimate_table_height

st.set_page_config(
    page_title="Home - DisGenix CNB",
    page_icon="../assets/disgenixLogo.png",
)

components.html(
    """
    <script src="https://cdn.counter.dev/script.js" 
            data-id="7573e275-619f-450a-b1c1-a2a3b586a601" 
            data-utcoffset="1"></script>
    """,
    height=0,  # no visible element needed
)

st.markdown("""
    <style>
    .fixed-logo {
        position: fixed;
        top: 0;
        left: 0;
        padding: 5px;
        z-index: 9999;
        background: transparent;
    }
    html, body, .stApp {
        width: 100%;
        margin: 0;
        padding: 0;
        overflow-x: hidden;
        font-size: 18px;
    }
    .block-container {
        max-width: 80% !important;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    .main {
        background: linear-gradient(to bottom right, #e3f2fd, #fce4ec);
        font-family: 'Arial', sans-serif;
        padding: 10px;
    }
    .stButton > button {
        background-color: #6a1b9a;
        color: white;
        font-weight: bold;
        border-radius: 2px;
    }
    .stDownloadButton > button {
        background-color: #2e7d32;
        color: white;
        font-weight: bold;
        border-radius: 2px;
    }
    </style>
""", unsafe_allow_html=True)

# Image set with st.image()
col1, col2 = st.columns([1,5],gap=None, vertical_alignment='bottom')
with col1:
    st.image("../assets/disgenixLogo.png", width=200)
with col2:
    st.image("../assets/CNB_2025.png", width=200)

st.logo(image="../assets/disgenixLogo2.png", size="large")
# Title and description
st.title("DisGenix")
st.markdown(
    "<p style='font-size:18px; font-weight:bold;'>Integrative analysis of upregulated disease genes, pathways, and drug interactions.</p>",
    unsafe_allow_html=True
)

# Button to go to demo page
if st.button("Pre-computed example of MeSH ID: D003110"):
    # Redirect to the demo
    st.switch_page("pages/2_Demo.py")

# Style and layout (tooltip)
st.markdown("""
<style>        
.email-label-wrapper {
  display: flex;
  align-items: center;
  gap: 6px;
  font-weight: 500;
  font-size: 16px;
  margin-bottom: -10px;
}

.tooltip {
  position: relative;
  display: inline-block;
  font-size: 14px;
  vertical-align: middle;
  cursor: pointer;
  color: #1E90FF;
}

.tooltip .tooltiptext {
  visibility: hidden;
  width: 280px;
  background-color: #555;
  color: #fff;
  text-align: left;
  border-radius: 6px;
  padding: 6px;
  position: absolute;
  z-index: 1;
  bottom: 125%;  /* aparece encima del icono */
  left: 50%;
  margin-left: -140px;
  opacity: 0;
  transition: opacity 0.3s;
}

.tooltip .tooltiptext::after {
  content: "";
  position: absolute;
  top: 100%;
  left: 50%;
  margin-left: -5px;
  border-width: 5px;
  border-style: solid;
  border-color: #555 transparent transparent transparent;
}

.tooltip:hover .tooltiptext {
  visibility: visible;
  opacity: 1;
}
</style>
""", unsafe_allow_html=True)

#Step 1: Text input for email
st.write('Enter user e-mail ℹ️: ')
email = st.text_input("User e-mail:", key="user_email", label_visibility="collapsed")
Entrez.email = email

options = ['Disease name', 'Disease MeSH ID']
searchBy = st.pills('Input: ', options, selection_mode="single", default=None)

#Step 2: Retrieve MeSH ID
spinner = st.spinner

if searchBy:
    normalized_disease = None
    if searchBy == options[1]:
        input = st.text_input("🔍 Enter MeSH ID (e.g., D003920 for Diabetes Mellitus):")
        if input:
            mesh_id = input
            with st.spinner("Fetching disease information..."):
                disease, disease_url = get_disease_name(mesh_id)

            if not disease:
                st.error("No valid disease found. "
                         "Please enter a recognized MeSH ID.")
                st.stop()

            if disease:
                normalized_disease = normalize_disease_name(disease)
                st.success(f"🎯 Disease **[{normalized_disease}]({disease_url})** identified")
    elif searchBy == options[0]:
        with st.spinner("Downloading MeSH data..."):
            try:
                mesh_df = get_mesh_df()
            except Exception as e:
                st.error(f"Problem downloading MeSH data: {e}")
                st.stop()

        input_text = st.text_input("🔍 Enter disease name (type at least 3 letters):")

        normalized_disease = None
        filtered_names = []
        if input_text and len(input_text) >= 3:
            filtered_names = mesh_df[mesh_df["Name"].str.contains(input_text, case=False, na=False)]["Name"].tolist()

        if filtered_names:
            selected_disease = st.selectbox("Matching diseases:", filtered_names)
            if selected_disease:
                mesh_id = mesh_df.loc[mesh_df["Name"] == selected_disease, "MeSH_ID"].values[0]
                disease_url = f"https://www.ncbi.nlm.nih.gov/mesh/{mesh_id}"
                normalized_disease = selected_disease
                st.success(f"🎯 Disease **[{normalized_disease}]({disease_url})** identified")

    uploaded_file = None
    if normalized_disease:
        # Step 3: File upload only after disease name is given
        url = generate_expression_atlas_link(disease_name=normalized_disease)
        uploaded_file = st.file_uploader(
            f"📁 Upload Differential Expression File (TSV from Expression Atlas: [link]({url}))",
            type=["tsv"]
        )


    if uploaded_file:
        df_raw = pd.read_csv(uploaded_file, sep="\t")
        st.write("Uploaded correctly")

        with st.spinner("Mapping Ensembl IDs to gene names..."):
            try:
                df_selected = fetch_gene_names(df_raw)
            except ValueError as e:
                st.error(f"Problem with uploaded file: {e}")
                st.stop()

        # Step 4: Obtain Ensembl ID with links for the genes
        df_selected_with_links = df_selected.copy()
        df_selected_with_links["Gene"] = df_selected_with_links["Gene"].apply(
            lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
        )
        df_selected_with_links = df_selected_with_links.sort_values(by="log_2 fold change", ascending=False)
        #Rename columns to keep consistency
        df_selected_with_links_newNames = df_selected_with_links.rename(columns={
            "Gene":"Ensembl ID",
            "Gene Name":"Gene"
        })


        # Title
        st.markdown("## Gene Table with Links to Ensembl")
        # Check if the DataFrame is empty
        if df_selected_with_links_newNames.empty:
            st.warning("The gene table is empty. Stopping the program.")
            st.stop()

        # Conversion of dataframe to HTML
        html_table = df_selected_with_links_newNames.to_html(escape=False, index=False, table_id="geneTable")

        # HTML code with CSS fixed in black and white
        html_code = f"""
        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

        <style>
        /* Apply style to whole table */
        .dataTables_wrapper, .dataTables_wrapper * {{
            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
            color: #31333f !important;
            font-size: 15px !important;
        }}

        /* Table */
        #geneTable, #geneTable thead, #geneTable tbody, #geneTable tr, #geneTable td, #geneTable th {{
            background-color: white !important;
            color: #31333f !important;
            border-color: #ccc !important;
        }}

        /* Headers */
        #geneTable thead th {{
            font-weight: 600 !important;
        }}

        /* Search inputs */
        #geneTable input {{
            background-color: white !important;
            color: #31333f !important;
            font-size: 14px !important;
            border: 1px solid #ccc !important;
        }}

        /* Entry selections */
        .dataTables_length select {{
            background-color: white !important;
            color: #31333f !important;
            border: 1px solid #ccc !important;
            font-size: 14px !important;
        }}

        /* Results info */
        .dataTables_info {{
            color: #31333f !important;
        }}

        /* Pages */
        .dataTables_paginate a {{
            color: #31333f !important;
            font-size: 14px !important;
            font-weight: normal !important;
            padding: 6px 12px;
            border-radius: 6px;
            text-decoration: none;
            margin: 0 2px;
        }}

        .dataTables_paginate a.current {{
            background-color: #f0f0f0 !important;
            font-weight: bold !important;
            border: 1px solid #aaa !important;
        }}

        .dataTables_paginate a:hover {{
            background-color: #e3e3e3 !important;
        }}
        </style>

        <script>
        $(document).ready(function() {{
            var table = $('#geneTable').DataTable({{
                scrollCollapse: true,
                paging: true,
                orderCellsTop: true,
                fixedHeader: true
            }});

            // Add search inputs to each column
            $('#geneTable thead tr').clone(true).appendTo('#geneTable thead');
            $('#geneTable thead tr:eq(1) th').each(function(i) {{
                var title = $(this).text();
                $(this).html('<input type="text" placeholder="Search ' + title + '" />');

                $('input', this).on('keyup change', function () {{
                    if (table.column(i).search() !== this.value) {{
                        table.column(i).search(this.value).draw();
                    }}
                }});
            }});
        }});
        </script>

        <div style="overflow-x:auto; margin-bottom: 0px;"">
        {html_table}
        </div>
        """

        height = estimate_table_height(df_selected)


        # Show table
        components.html(html_code, height=height, scrolling=True)

        # Download button
        st.download_button(
            "📥 Download Genes CSV",
            df_selected_with_links.to_csv(index=False),
            "genes.csv",
            "text/csv")

        df_selected["Gene Name_raw"] = df_selected["Gene Name"]

        gname_fc = df_selected.groupby("Gene Name_raw", as_index=False)["log_2 fold change"].sum()
        gname_fc = gname_fc.sort_values("log_2 fold change", ascending=False)

        # Graph for genes and their log2 fold change
        fig = px.bar(
            gname_fc,
            x="Gene Name_raw",
            y="log_2 fold change",
            color="log_2 fold change",
            color_continuous_scale="sunset",
            title="Sum of log₂ fold change per gene",
            labels={"Gene Name_raw": "Gene Name", "log_2 fold change": "Sum of log₂ fold change"},
            height=500,
            width=2000
        )

        st.plotly_chart(fig, width="stretch")


        # Step 5: Search biotype and tractability for the genes in Open Targets
        with st.spinner("Checking Open Targets..."):
            openTargets_results = cached_open_targets(tuple(df_selected["Gene"].unique()))
            openTargets_df = pd.DataFrame(openTargets_results)

        if not openTargets_df.empty:
            if "Name" in openTargets_df.columns:
                openTargets_df = openTargets_df.drop(columns=["Name"])

                openTargets_df = openTargets_df.sort_values(by="Gene Symbol")

                def format_tractability(tags):
                    if not tags:
                        return "NaN"
                    return " ".join(
                        f'<span style="background-color:#d1c4e9; color:#4a148c; padding:4px 8px; border-radius:10px; margin:2px; display:inline-block;">{tag}</span>'
                        for tag in tags
                    )

                openTargets_df["Tractability"] = openTargets_df["Tractability"].apply(format_tractability)
                openTargets_df["Gene Symbol"] = openTargets_df["Gene Symbol"].apply(
                    lambda x: f'<span title="{x}">{x[:10]}...</span>' if len(x) > 10 else x
                )

                openTargets_df["Ensembl ID"] = openTargets_df["Ensembl ID"].apply(
                    lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                )

                openTargets_df_newNames = openTargets_df.rename(columns={
                    "Gene Symbol":"Gene"
                })

                st.markdown("# Genes with Tractability (Open Targets)")
                # Add a helpful link to the overview page
                st.markdown(
                    """
                    ℹ️ For details on how tractability is defined, see the 
                    [Open Targets Tractability Overview](https://platform-docs.opentargets.org/target/tractability).
                    """
                )

                html_ot_table = openTargets_df_newNames.to_html(escape=False, index=False, table_id="openTargetsTable")

                html_ot_scroll = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                    /* Apply style to whole table */
                    .dataTables_wrapper, .dataTables_wrapper * {{
                        font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                        color: #31333f !important;
                        font-size: 15px !important;
                    }}

                    /* Table */
                    #openTargetsTable, #openTargetsTable thead, #openTargetsTable tbody, #openTargetsTable tr, #openTargetsTable td, #openTargetsTable th {{
                        background-color: white !important;
                        color: #31333f !important;
                        border-color: #ccc !important;
                    }}

                    /* Headers */
                    #openTargetsTable thead th {{
                        font-weight: 600 !important;
                    }}

                    /* Search inputs */
                    #openTargetsTable input {{
                        background-color: white !important;
                        color: #31333f !important;
                        font-size: 14px !important;
                        border: 1px solid #ccc !important;
                    }}

                    /* Entry selections */
                    .dataTables_length select {{
                        background-color: white !important;
                        color: #31333f !important;
                        border: 1px solid #ccc !important;
                        font-size: 14px !important;
                    }}

                    /* Results info */
                    .dataTables_info {{
                        color: #31333f !important;
                    }}

                    /* Pages */
                    .dataTables_paginate a {{
                        color: #31333f !important;
                        font-size: 14px !important;
                        font-weight: normal !important;
                        padding: 6px 12px;
                        border-radius: 6px;
                        text-decoration: none;
                        margin: 0 2px;
                    }}

                    .dataTables_paginate a.current {{
                        background-color: #f0f0f0 !important;
                        font-weight: bold !important;
                        border: 1px solid #aaa !important;
                    }}

                    .dataTables_paginate a:hover {{
                        background-color: #e3e3e3 !important;
                    }}
                    </style>

                    <script>
                    $(document).ready(function() {{
                        var table = $('#openTargetsTable').DataTable({{
                            scrollCollapse: true,
                            paging: true,
                            orderCellsTop: true,
                            fixedHeader: true
                        }});

                        // Add search inputs to each column
                        $('#openTargetsTable thead tr').clone(true).appendTo('#openTargetsTable thead');
                        $('#openTargetsTable thead tr:eq(1) th').each(function(i) {{
                            var title = $(this).text();
                            $(this).html('<input type="text" placeholder="Search ' + title + '" />');

                            $('input', this).on('keyup change', function () {{
                                if (table.column(i).search() !== this.value) {{
                                    table.column(i).search(this.value).draw();
                                }}
                            }});
                        }});
                    }});
                    </script>
                    
                    <div style="overflow-x:auto">
                        {html_ot_table}
                    </div>
                """

                height = estimate_table_height(openTargets_df)
                components.html(html_ot_scroll, height=height+200, scrolling=True)

                st.download_button(
                    "📥 Download results from Open Targets",
                    openTargets_df.to_csv(index=False),
                    "openTargets_genes.csv",
                    "text/csv"
                )

                openTargets_df["Tractability_raw"] = openTargets_df["Tractability"]
                openTargets_df["Biotype_raw"] = openTargets_df["Biotype"]

                tract_tags = openTargets_df["Tractability_raw"].explode()
                bio_tags = openTargets_df["Biotype_raw"].explode()

                tract_counts = tract_tags.value_counts().reset_index()
                tract_counts.columns = ["Tractability", "Count"]
                bio_counts = bio_tags.value_counts().reset_index()
                bio_counts.columns = ["Biotype", "Count"]

                # Tractability graph
                fig = px.bar(
                    tract_counts,
                    x="Tractability",
                    y="Count",
                    color="Count",
                    color_continuous_scale="burg",
                    title="Count of genes by Tractability",
                    labels={"Tractability": "Type of tractability", "Count": "Number of genes"},
                    height=700,
                    width= 1000
                )
                fig.update_xaxes(showticklabels=False)
                st.plotly_chart(fig, width="stretch")

                # Biotype graph
                fig = px.bar(
                    bio_counts,
                    x="Biotype",
                    y="Count",
                    color="Count",
                    color_continuous_scale="blues",
                    title="Count of genes by Biotype",
                    labels={"Biotype": "Biotype", "Count": "Number of genes"},
                    height=500,
                    width= 1000
                )

                st.plotly_chart(fig, width="stretch")

            else:
                st.warning("No results found in Open Targets.")


        with st.spinner("Performing pathway analysis..."):

            # Step 6: Search pathways for the genes
            st.markdown("# Top Enriched Reactome Pathways")
            # Add a helpful link to the overview page
            st.markdown(
                """
                ℹ️ For details on how the columns are defined, see the 
                [Enrichr Help Center](https://maayanlab.cloud/Enrichr/help#background).
                """
            )
            st.markdown(
                """
                ℹ️ **Custom column in this table**  
                - **Input %**: Percentage of the pathway’s genes present in your input list. 
                """
            )

            number_pathways = st.text_input("🔍 Enter number of pathways to retrieve", value="10")
            if number_pathways:
                try:
                    number_pathways = int(number_pathways)
                    top_pathways = analyze_pathways(df_selected, number_pathways)
                    top_pathways["-log10(Adj P)"] = -np.log10(top_pathways["Adjusted P-value"])
                    top_pathways = top_pathways.sort_values("Adjusted P-value", ascending=True)

                    if top_pathways is not None:
                        html_pathway_table = top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)", "Overlap", "Input %", "Sum log2fc"]].to_html(escape=False, index=False, table_id="topPathwayTable")

                        html_code_top_pathways = f"""
                        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                        <style>
                            body, .dataTables_wrapper, .dataTables_wrapper * {{
                                font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                                color: #31333f !important;
                                font-size: 15px !important;
                            }}

                            #topPathwayTable {{
                                width: 100% !important;
                            }}

                            #topPathwayTable, #topPathwayTable thead, #topPathwayTable tbody, #topPathwayTable tr, #topPathwayTable td, #topPathwayTable th {{
                                background-color: white !important;
                                color: #31333f !important;
                                border-color: #ccc !important;
                            }}

                            #topPathwayTable thead th {{
                                font-weight: 600 !important;
                            }}

                            #topPathwayTable input {{
                                width: 100%;
                                box-sizing: border-box;
                                background-color: white !important;
                                color: #31333f !important;
                                font-size: 14px !important;
                                border: 1px solid #ccc !important;
                            }}

                            .dataTables_length select {{
                                background-color: white !important;
                                color: #31333f !important;
                                border: 1px solid #ccc !important;
                                font-size: 14px !important;
                            }}

                            .dataTables_info {{
                                color: #31333f !important;
                            }}

                            .dataTables_paginate a {{
                                color: #31333f !important;
                                font-size: 14px !important;
                                font-weight: normal !important;
                                padding: 6px 12px;
                                border-radius: 6px;
                                text-decoration: none;
                                margin: 0 2px;
                            }}

                            .dataTables_paginate a.current {{
                                background-color: #f0f0f0 !important;
                                font-weight: bold !important;
                                border: 1px solid #aaa !important;
                            }}

                            .dataTables_paginate a:hover {{
                                background-color: #e3e3e3 !important;
                            }}
                        </style>

                        <script>
                            $(document).ready(function() {{
                                var table = $('#topPathwayTable').DataTable({{
                                    scrollCollapse: true,
                                    paging: true,
                                    orderCellsTop: true,
                                    fixedHeader: false,
                                    autoWidth: false
                                }});

                                $('#topPathwayTable thead tr').clone(true).appendTo('#topPathwayTable thead');
                                $('#topPathwayTable thead tr:eq(1) th').each(function(i) {{
                                    var title = $(this).text();
                                    $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                                }});

                                $('#topPathwayTable thead').on('keyup change', 'input', function () {{
                                    let i = $(this).parent().index();
                                    table.column(i).search(this.value).draw();
                                }});
                            }});
                        </script>

                        <div style="overflow-x:auto">
                            {html_pathway_table}
                        </div>
                        """
                        height = estimate_table_height(top_pathways)
                        components.html(html_code_top_pathways, height=height, scrolling=True)

                        st.download_button(
                            "📥 Download Top N pathways CSV",
                            top_pathways.to_csv(index=False),
                            "topPathways.csv",
                            "text/csv")

                        st.markdown("# Important genes in pathway: ")
                        selected_pathway = st.selectbox(
                            "🔍 Select a pathway to see the top genes",
                            top_pathways["Term"].tolist(),
                            key="pathway"
                        )

                        # Get full row of selected pathway
                        selected_pathway_row = top_pathways[top_pathways["Term"] == selected_pathway].iloc[0]

                        # Get overlapping genes for that pathway
                        pathway_genes = get_overlapping_genes(df_selected, selected_pathway_row)

                        # Display pathway name
                        st.subheader(f"{selected_pathway}")

                        # Turn Gene IDs into Ensembl search links
                        pathway_genes["Gene"] = pathway_genes["Gene"].apply(
                            lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                        )

                        pathway_genes = pathway_genes.drop(columns=["Gene Name_raw", "abs_fc"])

                        pathway_genes_newNames = pathway_genes.rename(columns={
                            "Gene":"Ensembl ID",
                            "Gene Name":"Gene"
                        })

                        # Render HTML table with scroll
                        html_genespathway_table = pathway_genes_newNames.to_html(escape=False, index=False, table_id="pathwayGeneTable")


                        # HTML and DataTables JS
                        html_code_pathway = f"""
                        <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                        <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                        <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                        <style>
                            .dataTables_wrapper, .dataTables_wrapper * {{
                                font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                                color: #31333f !important;
                                font-size: 15px !important;
                            }}

                            #pathwayGeneTable {{
                                width: 100% !important;
                            }}

                            #pathwayGeneTable, #pathwayGeneTable thead, #pathwayGeneTable tbody, #pathwayGeneTable tr, #pathwayGeneTable td, #pathwayGeneTable th {{
                                background-color: white !important;
                                color: #31333f !important;
                                border-color: #ccc !important;
                            }}

                            #pathwayGeneTable thead th {{
                                font-weight: 600 !important;
                            }}

                            #pathwayGeneTable input {{
                                width: 100%;
                                box-sizing: border-box;
                                background-color: white !important;
                                color: #31333f !important;
                                font-size: 14px !important;
                                border: 1px solid #ccc !important;
                            }}

                            .dataTables_length select {{
                                background-color: white !important;
                                color: #31333f !important;
                                border: 1px solid #ccc !important;
                                font-size: 14px !important;
                            }}

                            .dataTables_info {{
                                color: #31333f !important;
                            }}

                            .dataTables_paginate a {{
                                color: #31333f !important;
                                font-size: 14px !important;
                                font-weight: normal !important;
                                padding: 6px 12px;
                                border-radius: 6px;
                                text-decoration: none;
                                margin: 0 2px;
                            }}

                            .dataTables_paginate a.current {{
                                background-color: #f0f0f0 !important;
                                font-weight: bold !important;
                                border: 1px solid #aaa !important;
                            }}

                            .dataTables_paginate a:hover {{
                                background-color: #e3e3e3 !important;
                            }}
                        </style>

                        <script>
                            $(document).ready(function() {{
                                var table = $('#pathwayGeneTable').DataTable({{
                                    scrollCollapse: true,
                                    paging: true,
                                    orderCellsTop: true,
                                    fixedHeader: false,
                                    autoWidth: false
                                }});

                                $('#pathwayGeneTable thead tr').clone(true).appendTo('#pathwayGeneTable thead');
                                $('#pathwayGeneTable thead tr:eq(1) th').each(function(i) {{
                                    var title = $(this).text();
                                    $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                                }});

                                $('#pathwayGeneTable thead').on('keyup change', 'input', function () {{
                                    let i = $(this).parent().index();
                                    table.column(i).search(this.value).draw();
                                }});
                            }});
                        </script>

                        <div style="overflow-x:auto">
                            {html_genespathway_table}
                        </div>
                        """

                        height = estimate_table_height(pathway_genes)
                        components.html(html_code_pathway, height=height, scrolling=True)

                        st.download_button(
                            "📥 Download Important Genes CSV",
                            pathway_genes.to_csv(index=False),
                            "genes_pathway.csv",
                            "text/csv"
                        )

                except ValueError:
                    st.error("Please enter a valid integer.")

                st.markdown("# Drug-Gene Interactions")
                with st.spinner("Searching drug-gene interactions..."):

                    selected_pathway = st.selectbox(
                        "🔍 Select a pathway to see drugs",
                    top_pathways["Term"].tolist(),
                    key="drugs"
                    )
                    selected_pathway_row = top_pathways[top_pathways["Term"] == selected_pathway].iloc[0]

                    pathway_genes = get_overlapping_genes(df_selected, selected_pathway_row)
                    drug_df = cached_dgidb(tuple(sorted(pathway_genes["Gene Name"].unique())))

                if not drug_df.empty:

                    drug_df_with_links = drug_with_links(drug_df)

                    drug_df_with_links = drug_df_with_links.replace(r'^\s*$', np.nan, regex=True)

                    # HTML table
                    html_drug_table = drug_df_with_links.to_html(escape=False, index=False, table_id="drugTable")

                    html_code_drug = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                        .dataTables_wrapper, .dataTables_wrapper * {{
                            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                            color: #31333f !important;
                            font-size: 15px !important;
                        }}

                        #drugTable {{
                            width: 100% !important;
                        }}

                        #drugTable, #drugTable thead, #drugTable tbody, #drugTable tr, #drugTable td, #drugTable th {{
                            background-color: white !important;
                            color: #31333f !important;
                            border-color: #ccc !important;
                        }}

                        #drugTable thead th {{
                            font-weight: 600 !important;
                        }}

                        #drugTable input {{
                            width: 100%;
                            box-sizing: border-box;
                            background-color: white !important;
                            color: #31333f !important;
                            font-size: 14px !important;
                            border: 1px solid #ccc !important;
                        }}

                        .dataTables_length select {{
                            background-color: white !important;
                            color: #31333f !important;
                            border: 1px solid #ccc !important;
                            font-size: 14px !important;
                        }}

                        .dataTables_info {{
                            color: #31333f !important;
                        }}

                        .dataTables_paginate a {{
                            color: #31333f !important;
                            font-size: 14px !important;
                            font-weight: normal !important;
                            padding: 6px 12px;
                            border-radius: 6px;
                            text-decoration: none;
                            margin: 0 2px;
                        }}

                        .dataTables_paginate a.current {{
                            background-color: #f0f0f0 !important;
                            font-weight: bold !important;
                            border: 1px solid #aaa !important;
                        }}

                        .dataTables_paginate a:hover {{
                            background-color: #e3e3e3 !important;
                        }}
                    </style>

                    <script>
                        $(document).ready(function() {{
                            var table = $('#drugTable').DataTable({{
                                scrollCollapse: true,
                                paging: true,
                                orderCellsTop: true,
                                fixedHeader: false,
                                autoWidth: false
                            }});

                            $('#drugTable thead tr').clone(true).appendTo('#drugTable thead');
                            $('#drugTable thead tr:eq(1) th').each(function(i) {{
                                var title = $(this).text();
                                $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                            }});

                            $('#drugTable thead').on('keyup change', 'input', function () {{
                                let i = $(this).parent().index();
                                table.column(i).search(this.value).draw();
                            }});
                        }});
                    </script>

                    <div style="overflow-x:auto">
                        {html_drug_table}
                    </div>
                    """
                    height = estimate_table_height(drug_df)
                    components.html(html_code_drug, height=height, scrolling=True)


                    st.download_button("📥 Download Drug Interactions CSV", drug_df.to_csv(index=False), "drug_interactions.csv", "text/csv")

                    # Gene selection and plot
                    unique_genes = drug_df_with_links['Gene'].apply(lambda x: re.search(r'>(.*?)<', x).group(1)).unique()
                    selected_gene = st.selectbox("Select a gene to view its drug interaction scores", unique_genes)

                    # Filter by selected gene
                    gene_df = drug_df_with_links[drug_df['Gene'] == selected_gene]
                    gene_df = gene_df.sort_values(by="Interaction Score", ascending=False)

                    # Plot
                    fig = px.bar(
                        gene_df,
                        x='Drug',
                        y='Interaction Score',
                        hover_data=['Interaction Type', 'Source', 'PMID'],
                        title=f"Drug Interactions for {selected_gene}",
                        labels={'Interaction Score': 'Interaction Score', 'Drug': 'Drug'},
                    )
                    fig.update_layout(xaxis_tickangle=-45)

                    st.plotly_chart(fig)

                    if 'Interaction Type' in drug_df.columns:
                        st.markdown("# Distribution of Interaction Types")

                        # Count each type
                        interaction_counts = drug_df['Interaction Type'].value_counts().reset_index()
                        interaction_counts.columns = ['Interaction Type', 'Count']

                        # Pie chart
                        pie_fig = px.pie(
                            interaction_counts,
                            values='Count',
                            names='Interaction Type',
                            title='Percentage of Interaction Types',
                            hole=0.4
                        )

                        # Group by type of interaction
                        interaction_grouped = (
                            drug_df.groupby('Interaction Type', group_keys=False)
                            .apply(lambda df: [f"{g} → {d}" for g, d in zip(df['Gene'], df['Drug'])], include_groups=False)
                        )



                        # Convert to DataFrame with columns according to interaction type
                        max_len = interaction_grouped.map(len).max()
                        interaction_wide = pd.DataFrame({
                            interaction_type: values + [""] * (max_len - len(values))
                            for interaction_type, values in interaction_grouped.items()
                        })

                        # Layout: 2 columns
                        col1, col2 = st.columns([1.5, 1.8])

                        with col1:
                            st.plotly_chart(pie_fig)

                        with col2:
                            st.markdown("**Gene–Drug Interactions by Interaction Type**")
                            st.dataframe(interaction_wide, width="stretch")

                    else:
                        st.info("No 'Interaction Type' column found.")

                else:
                    st.warning("No drug-gene interactions found.")

            else:
                st.warning("No enriched pathways found.")
        # -- Merge all data into a full report table --
        with st.spinner("Creating summary table..."):

            # First, clean and prepare dataframes
            merged = df_selected.copy()

            # Merge with Open Targets
            if not openTargets_df.empty:
                ot_clean = openTargets_df.copy()
                ot_clean["Tractability"] = ot_clean["Tractability"].str.replace(r'<.*?>', '', regex=True)
                merged = pd.merge(merged, ot_clean.drop(columns=["Biotype_raw"]), how="left", left_on="Gene Name", right_on="Gene Symbol")

                merged["abs_fc"] = merged["log_2 fold change"].abs()

                # Drug–gene interactions for all genes
                all_drug_df = cached_dgidb(tuple(sorted(merged["Gene Name"].unique())))

            if not all_drug_df.empty:
                drug_clean = all_drug_df.copy()
                drug_summary = drug_clean.groupby("Gene").agg({
                    "Drug": lambda x: "; ".join(sorted(set(x))),
                    "Interaction Type": lambda x: "; ".join(sorted(set(x))),
                    "PMID": lambda x: "; ".join(sorted(set(str(p) for p in x if p != 'N/A'))),
                    "Interaction Score": "mean"
                }).reset_index()

                merged = pd.merge(merged, drug_summary, how="left", left_on="Gene Name", right_on="Gene")

            if not top_pathways.empty:
                gene_to_pathways = defaultdict(list)
                for _, row in top_pathways.iterrows():
                    pathway_name = row["Term"]
                    genes_in_pathway = [g.strip().upper() for g in row["Genes"].split(";")]
                    for gene in genes_in_pathway:
                        gene_to_pathways[gene].append(pathway_name)

                merged["Gene Name"] = merged["Gene Name"].astype(str).str.strip().str.upper()
                merged["Pathways"] = merged["Gene Name"].apply(lambda g: "; ".join(gene_to_pathways.get(g, [])))

                drug_summary = drug_clean.groupby("Gene").agg({
                        "Drug": lambda x: "; ".join(sorted(set(x))),
                        "Interaction Type": lambda x: "; ".join(sorted(set(x))),
                        "PMID": lambda x: "; ".join(sorted(set(str(p) for p in x if p != 'N/A'))),
                        "Interaction Score": "mean"
                        }).reset_index()
                merged = pd.merge(merged, drug_summary, how="left", left_on="Gene Name", right_on="Gene")
                # Drop the right-side 'Gene' from drug_summary to not be redundant
                merged = merged.drop(columns=["Gene_y"], errors="ignore")

                # Rename Gene_x column
                if "Gene_x" in merged.columns:
                        merged = merged.rename(columns={"Gene_x": "Gene"})

                # Remove suffixes from Drug/Interaction columns when duplicated
                for col in merged.columns:
                        if col.endswith("_x") and col[:-2] + "_y" in merged.columns:
                            # Keep _x version, drop _y
                            merged = merged.drop(columns=[col[:-2] + "_y"])
                            merged = merged.rename(columns={col: col[:-2]})


                if not top_pathways.empty:
                    pathway_genes_unique = pathway_genes[["Gene Name"]].drop_duplicates()
                    pathway_genes_unique["Pathway Match"] = True
                    merged = pd.merge(
                            merged,
                            pathway_genes[["Gene Name", "abs_fc"]],
                            how="left",
                            on="Gene Name"
                        )


                #Delete duplicated columns
                merged = merged.drop(columns=["Gene_y", "Gene Symbol", "Tractability_raw", "Gene Name_raw_x", "Gene Name_raw_y", "abs_fc_y", "log_2 fold change_y"], errors="ignore")
                merged = merged.rename(columns={"Gene_x": "Gene", "abs_fc_x" : "abs_fc"})
                dupes = merged.columns[merged.columns.duplicated()].tolist()
                merged = merged.loc[:, ~merged.columns.duplicated(keep="first")]

                #Add pathways in which gene interacts
                gene_to_pathways = defaultdict(list)

                if not top_pathways.empty:
                    for _, row in top_pathways.iterrows():
                        pathway_name = row["Term"]
                        genes_in_pathway = [g.strip().upper() for g in row["Genes"].split(";")]
                        for gene in genes_in_pathway:
                            gene_to_pathways[gene].append(pathway_name)

                merged["Gene Name"] = merged["Gene Name"].astype(str).str.strip().str.upper()
                merged["Pathways"] = merged["Gene Name"].apply(lambda g: "; ".join(gene_to_pathways.get(g, [])))
                merged = merged.drop(columns=["Gene Name_raw", "abs_fc", "Ensembl ID"])

                merged["Gene"] = merged["Gene"].apply(
                        lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                    )

                merged_newNames = merged.rename(columns={
                    "Gene":"Ensembl ID",
                    "Gene Name":"Gene"
                })

                # Apply to your merged table
                merged_with_links = add_links_to_final_table(merged_newNames)
                merged_with_links = merged_with_links.replace(r'^\s*$', np.nan, regex=True)


                # Show and download button
                st.markdown("## 📦 Download Full Results Table")

                html_final_table = merged_with_links.to_html(escape=False, index=False, table_id="finalTable")

                html_code_final = f"""
                    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/jquery.dataTables.min.css">
                    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
                    <script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>

                    <style>
                        .dataTables_wrapper, .dataTables_wrapper * {{
                            font-family: "Segoe UI", "Helvetica", "Arial", sans-serif !important;
                            color: #31333f !important;
                            font-size: 12px !important;
                        }}

                        #finalTable {{
                            width: 100% !important;
                        }}

                        #finalTable, #finalTable thead, #finalTable tbody, #finalTable tr, #finalTable td, #finalTable th {{
                            background-color: white !important;
                            color: #31333f !important;
                            border-color: #ccc !important;
                        }}

                        #finalTable thead th {{
                            font-weight: 600 !important;
                        }}

                        #finalTable input {{
                            width: 100%;
                            box-sizing: border-box;
                            background-color: white !important;
                            color: #31333f !important;
                            font-size: 14px !important;
                            border: 1px solid #ccc !important;
                        }}

                        .dataTables_length select {{
                            background-color: white !important;
                            color: #31333f !important;
                            border: 1px solid #ccc !important;
                            font-size: 14px !important;
                        }}

                        .dataTables_info {{
                            color: #31333f !important;
                        }}

                        .dataTables_paginate a {{
                            color: #31333f !important;
                            font-size: 14px !important;
                            font-weight: normal !important;
                            padding: 6px 12px;
                            border-radius: 6px;
                            text-decoration: none;
                            margin: 0 2px;
                        }}

                        .dataTables_paginate a.current {{
                            background-color: #f0f0f0 !important;
                            font-weight: bold !important;
                            border: 1px solid #aaa !important;
                        }}

                        .dataTables_paginate a:hover {{
                            background-color: #e3e3e3 !important;
                        }}
                    </style>

                    <script>
                        $(document).ready(function() {{
                            var table = $('#finalTable').DataTable({{
                                scrollCollapse: true,
                                paging: true,
                                orderCellsTop: true,
                                fixedHeader: false,
                                autoWidth: false
                            }});

                            $('#finalTable thead tr').clone(true).appendTo('#finalTable thead');
                            $('#finalTable thead tr:eq(1) th').each(function(i) {{
                                var title = $(this).text();
                                $(this).html('<input type="text" placeholder="Search ' + title + '" />');
                            }});

                            $('#finalTable thead').on('keyup change', 'input', function () {{
                                let i = $(this).parent().index();
                                table.column(i).search(this.value).draw();
                            }});
                        }});
                    </script>

                    <div style="overflow-x:auto">
                        {html_final_table}
                    </div>
                    """
                height = estimate_table_height(merged)
                components.html(html_code_final, height=height, scrolling=True)

                st.download_button(
                        "📥 Download Full Combined CSV",
                        merged_with_links.to_csv(index=False),
                        "full_results_table.csv",
                        "text/csv"
                )



                # --- ZIP with chosen tables (as they are) ---
                # Build the list of available dataframes (desc, filename, df)
                candidates = []
                try: candidates.append(("Genes with links",     "genes_with_links.csv", df_selected_with_links))
                except NameError: pass
                try: candidates.append(("Open Targets",         "openTargets_df.csv",   openTargets_df))
                except NameError: pass
                try: candidates.append(("Top pathways",         "top_pathways.csv",     top_pathways))
                except NameError: pass
                try: candidates.append(("Pathway genes",        "pathway_genes.csv",    pathway_genes))
                except NameError: pass
                try: candidates.append(("Drug interactions",    "drug_df.csv",          drug_df))
                except NameError: pass

                if candidates:
                    st.markdown("### Select tables to include in ZIP")

                    # Table headers
                    h1, h2 = st.columns([0.12, 0.88])
                    h1.markdown("** **")
                    h2.markdown("**Dataframe**")

                    # Individual checkboxes
                    for desc, fname, df in candidates:
                        c1, c2 = st.columns([0.12, 0.88])
                        c1.checkbox(
                            label=f"Include {desc}",
                            key=f"_dl_{fname}",
                            value=True,  # selected by default
                            label_visibility="collapsed"
                        )
                        c2.write(f"{desc} ({len(df)} rows)")

                    # Gather selected
                    selected = [(fname, df) for _, fname, df in candidates if st.session_state.get(f"_dl_{fname}", False)]

                # Create columns for the two download buttons
                col1, col2 = st.columns(2)

                # --- Column 1: Download selected tables ---
                with col1:
                    if not selected:
                        st.warning("Select at least one table to enable the download.")
                    else:
                        zip_buf = io.BytesIO()
                        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                            for fname, df in selected:
                                zf.writestr(fname, df.to_csv(index=False))
                        zip_buf.seek(0)

                        st.download_button(
                            label="📥 Download selected (.zip)",
                            data=zip_buf.getvalue(),
                            file_name=f"selected_tables_{datetime.date.today().isoformat()}.zip",
                            mime="application/zip",
                            key="zip_download_btn"
                        )

                # --- Column 2: Download all tables ---
                with col2:
                    # Define folders containing your CSVs
                    genes_folder = "all_pathway_genes_csvs"
                    drugs_folder = "all_drug_csvs"

                    # Optional extra CSVs
                    extra_files = [
                        "genes.csv",
                        "openTargets_genes.csv",
                        "topPathways.csv",
                        "genes_pathway.csv",
                        "drug_interactions.csv",
                        "full_results_table.csv"
                    ]

                    def count_csvs(folder):
                        return len([f for f in os.listdir(folder) if f.endswith(".csv")])

                    st.markdown(
                        "⚠️ **Note:** Clicking the button below will start the process of generating and packaging all tables "
                        "(it will generate all the tables for each pathway). This may take a few minutes depending on your system and file sizes."
                    )

                    if st.button("Generate and Download Tables"):
                        with st.spinner("Generating all CSVs..."):
                            drug_csvs = save_drug_csvs(df_selected, top_pathways)
                            pathway_csvs = save_pathway_csvs(df_selected, top_pathways)

                            # Build ZIP in memory
                            zip_buf = io.BytesIO()
                            with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
                                # Add pathway CSVs
                                for fname, data in pathway_csvs.items():
                                    zf.writestr(fname, data)
                                # Add drug CSVs
                                for fname, data in drug_csvs.items():
                                    zf.writestr(fname, data)
                                # Add any extra files
                                for f in extra_files:
                                    if os.path.exists(f):
                                        zf.write(f, arcname=os.path.basename(f))

                            zip_buf.seek(0)

                        st.download_button(
                            label="📥 Download All Tables",
                            data=zip_buf.getvalue(),
                            file_name="all_tables.zip",
                            mime="application/zip"
                        )

st.markdown("""
    <style>
    footer {
        visibility: hidden;
    }
    .footer {
        position: fixed;
        left: 0;
        bottom: 0;
        width: 100%;
        background-color: rgba(240,240,240,0.7);
        text-align: center;
        color: gray;
        font-size: 0.9em;
        padding: 8px 0;
    }
    </style>

    <div class="footer">
        © 2025 CNB – DisGenix 🧬
    </div>
""", unsafe_allow_html=True)