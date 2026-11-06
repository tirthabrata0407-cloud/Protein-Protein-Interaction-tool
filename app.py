
# ============================================================
# PROTEIN-PROTEIN INTERACTION ANALYZER
# Streamlit Web Application
#
# Input  : .pdb file
# Output : Excel (.xlsx) report
#
# Sheet 1: Interacting Residues (3.5 Angstrom)
# Sheet 2: Salt Bridges (4.0 Angstrom)
# ============================================================


# STEP 1: IMPORT LIBRARIES

import streamlit as st
import pandas as pd
import numpy as np
import os
import io

from Bio.PDB import PDBParser, is_aa


# ============================================================
# STEP 2: STREAMLIT PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="Protein PPI Analyzer",
    page_icon="🧬",
    layout="wide"
)


# ============================================================
# STEP 3: APPLICATION HEADER
# ============================================================

st.title("🧬 Protein-Protein Interaction Analyzer")

st.markdown(
    """
    Upload a protein structure in **`.pdb` format** to identify
    interactions between two protein chains.

    The application identifies:
    - Interacting residues and atom pairs within **3.5 Å**
    - Potential salt bridges within **4.0 Å**
    - Detailed interaction distances and residue information
    - Downloadable Excel report
    """
)

st.divider()


# ============================================================
# STEP 4: UPLOAD PDB FILE
# ============================================================

st.subheader("1. Upload Protein Structure")

uploaded_file = st.file_uploader(
    "Upload your protein PDB file",
    type=["pdb"],
    accept_multiple_files=False
)


if uploaded_file is not None:

    pdb_file = uploaded_file.name

    # Strict lowercase extension validation

    if not pdb_file.endswith(".pdb"):

        st.error(
            "Invalid file format! Please upload a file "
            "with .pdb extension in lowercase."
        )

        st.stop()

    st.success(f"Successfully uploaded: {pdb_file}")

    # Read uploaded PDB content

    pdb_content = uploaded_file.getvalue()

    # ========================================================
    # STEP 5: PARSE PDB STRUCTURE
    # ========================================================

    try:

        parser = PDBParser(QUIET=True)

        pdb_text = pdb_content.decode(
            "utf-8",
            errors="replace"
        )

        structure = parser.get_structure(
            "Protein",
            io.StringIO(pdb_text)
        )

        model = structure[0]

        # Identify protein chains

        chains = {}

        for chain in model:

            residues = [
                residue for residue in chain
                if is_aa(residue, standard=False)
            ]

            if residues:

                chains[chain.id] = residues


        if len(chains) < 2:

            st.error(
                "At least two protein chains are required "
                "for PPI analysis."
            )

            st.stop()


        # ====================================================
        # STEP 6: DISPLAY AVAILABLE CHAINS
        # ====================================================

        st.subheader("2. Protein Chain Information")

        chain_information = []

        for chain_id, residues in chains.items():

            chain_information.append({

                "Chain ID": chain_id,

                "Number of Residues": len(residues)

            })


        chain_df = pd.DataFrame(chain_information)

        st.dataframe(
            chain_df,
            use_container_width=True,
            hide_index=True
        )


        # ====================================================
        # STEP 7: SELECT TWO CHAINS
        # ====================================================

        st.subheader("3. Select Chains for Interaction Analysis")

        chain_ids = list(chains.keys())

        col1, col2 = st.columns(2)

        with col1:

            chain1_id = st.selectbox(
                "Select First Protein Chain",
                options=chain_ids,
                index=0
            )

        with col2:

            chain2_id = st.selectbox(
                "Select Second Protein Chain",
                options=chain_ids,
                index=1 if len(chain_ids) > 1 else 0
            )


        if chain1_id == chain2_id:

            st.warning(
                "Please select two different protein chains."
            )

            st.stop()


        # ====================================================
        # STEP 8: DEFINE DISTANCE PARAMETERS
        # ====================================================

        INTERACTION_CUTOFF = 3.5

        SALT_BRIDGE_CUTOFF = 4.0

        residues_1 = chains[chain1_id]

        residues_2 = chains[chain2_id]


        # ====================================================
        # STEP 9: HELPER FUNCTIONS
        # ====================================================

        def get_heavy_atoms(residue):

            return [
                atom for atom in residue
                if atom.element.upper() not in ["H", "D"]
            ]


        def residue_info(residue):

            return {

                "Residue": residue.resname,

                "Residue_Number": residue.id[1],

                "Insertion_Code": residue.id[2].strip(),

                "Residue_ID": (
                    str(residue.id[1]) + residue.id[2].strip()
                )

            }


        # ====================================================
        # STEP 10: PPI INTERACTION ANALYSIS
        # ====================================================

        def analyze_ppi(residues_1, residues_2):

            interaction_data = []

            for res1 in residues_1:

                atoms1 = get_heavy_atoms(res1)

                for res2 in residues_2:

                    atoms2 = get_heavy_atoms(res2)

                    for atom1 in atoms1:

                        for atom2 in atoms2:

                            distance = atom1 - atom2

                            if distance <= INTERACTION_CUTOFF:

                                r1 = residue_info(res1)

                                r2 = residue_info(res2)

                                interaction_data.append({

                                    "Chain_1": chain1_id,

                                    "Residue_1": r1["Residue"],

                                    "Residue_1_Number":
                                        r1["Residue_Number"],

                                    "Insertion_Code_1":
                                        r1["Insertion_Code"],

                                    "Residue_ID_1": r1["Residue_ID"],

                                    "Atom_1": atom1.name,

                                    "Chain_2": chain2_id,

                                    "Residue_2": r2["Residue"],

                                    "Residue_2_Number":
                                        r2["Residue_Number"],

                                    "Insertion_Code_2":
                                        r2["Insertion_Code"],

                                    "Residue_ID_2": r2["Residue_ID"],

                                    "Atom_2": atom2.name,

                                    "Distance_Angstrom":
                                        round(distance, 3)

                                })


            interaction_columns = [

                "Chain_1",

                "Residue_1",

                "Residue_1_Number",

                "Insertion_Code_1",

                "Residue_ID_1",

                "Atom_1",

                "Chain_2",

                "Residue_2",

                "Residue_2_Number",

                "Insertion_Code_2",

                "Residue_ID_2",

                "Atom_2",

                "Distance_Angstrom"

            ]


            interaction_df = pd.DataFrame(

                interaction_data,

                columns=interaction_columns

            )


            if not interaction_df.empty:

                interaction_df = interaction_df.sort_values(

                    by="Distance_Angstrom"

                )


            return interaction_df


        # ====================================================
        # STEP 11: SALT BRIDGE ANALYSIS
        # ====================================================

        def analyze_salt_bridges(residues_1, residues_2):

            positive_atoms = {

                "ARG": ["NE", "NH1", "NH2"],

                "LYS": ["NZ"],

                "HIS": ["ND1", "NE2"]

            }


            negative_atoms = {

                "ASP": ["OD1", "OD2"],

                "GLU": ["OE1", "OE2"]

            }


            salt_bridge_data = []


            def get_charged_atoms(residue, atom_dictionary):

                atom_names = atom_dictionary.get(
                    residue.resname,
                    []
                )

                return [

                    atom for atom in residue

                    if atom.name in atom_names

                ]


            for res1 in residues_1:

                for res2 in residues_2:

                    charged_pairs = []


                    # Chain 1 positive, Chain 2 negative

                    if (

                        res1.resname in positive_atoms

                        and res2.resname in negative_atoms

                    ):

                        charged_pairs.append((

                            res1,

                            res2,

                            positive_atoms,

                            negative_atoms

                        ))


                    # Chain 1 negative, Chain 2 positive

                    if (

                        res1.resname in negative_atoms

                        and res2.resname in positive_atoms

                    ):

                        charged_pairs.append((

                            res2,

                            res1,

                            positive_atoms,

                            negative_atoms

                        ))


                    for (

                        positive_res,

                        negative_res,

                        pos_dict,

                        neg_dict

                    ) in charged_pairs:


                        positive_list = get_charged_atoms(

                            positive_res,

                            pos_dict

                        )


                        negative_list = get_charged_atoms(

                            negative_res,

                            neg_dict

                        )


                        for pos_atom in positive_list:

                            for neg_atom in negative_list:

                                distance = pos_atom - neg_atom


                                if distance <= SALT_BRIDGE_CUTOFF:

                                    p = residue_info(positive_res)

                                    n = residue_info(negative_res)


                                    salt_bridge_data.append({

                                        "Positive_Chain":
                                            positive_res.get_parent().id,

                                        "Positive_Residue":
                                            p["Residue"],

                                        "Positive_Residue_Number":
                                            p["Residue_Number"],

                                        "Positive_Atom":
                                            pos_atom.name,

                                        "Negative_Chain":
                                            negative_res.get_parent().id,

                                        "Negative_Residue":
                                            n["Residue"],

                                        "Negative_Residue_Number":
                                            n["Residue_Number"],

                                        "Negative_Atom":
                                            neg_atom.name,

                                        "Distance_Angstrom":
                                            round(distance, 3),

                                        "Interaction_Type": (

                                            "Potential Histidine-dependent salt bridge"

                                            if "HIS" in [

                                                positive_res.resname,

                                                negative_res.resname

                                            ]

                                            else "Potential salt bridge"

                                        )

                                    })


            salt_bridge_columns = [

                "Positive_Chain",

                "Positive_Residue",

                "Positive_Residue_Number",

                "Positive_Atom",

                "Negative_Chain",

                "Negative_Residue",

                "Negative_Residue_Number",

                "Negative_Atom",

                "Distance_Angstrom",

                "Interaction_Type"

            ]


            salt_bridge_df = pd.DataFrame(

                salt_bridge_data,

                columns=salt_bridge_columns

            )


            if not salt_bridge_df.empty:

                salt_bridge_df = salt_bridge_df.sort_values(

                    by="Distance_Angstrom"

                )


            return salt_bridge_df


        # ====================================================
        # STEP 12: RUN ANALYSIS
        # ====================================================

        st.divider()

        st.subheader("4. Run Interaction Analysis")

        if st.button(
            "🔬 Analyze Protein Interactions",
            type="primary",
            use_container_width=True
        ):

            with st.spinner(
                "Analyzing protein-protein interactions..."
            ):

                interaction_df = analyze_ppi(
                    residues_1,
                    residues_2
                )

                salt_bridge_df = analyze_salt_bridges(
                    residues_1,
                    residues_2
                )


            st.success("Interaction analysis completed successfully!")


            # =================================================
            # STEP 13: DISPLAY ANALYSIS SUMMARY
            # =================================================

            st.subheader("5. Analysis Summary")

            col1, col2, col3 = st.columns(3)

            col1.metric(
                "Atom-Pair Contacts",
                len(interaction_df)
            )

            col2.metric(
                "Potential Salt Bridge Atom Pairs",
                len(salt_bridge_df)
            )

            col3.metric(
                "Selected Chains",
                f"{chain1_id} – {chain2_id}"
            )


            # =================================================
            # STEP 14: DISPLAY INTERACTION RESULTS
            # =================================================

            st.divider()

            st.subheader(
                "6. Interacting Residues (Distance ≤ 3.5 Å)"
            )

            if not interaction_df.empty:

                st.dataframe(
                    interaction_df,
                    use_container_width=True,
                    hide_index=True
                )

            else:

                st.info(
                    "No atom-pair contacts were detected "
                    "within 3.5 Å."
                )


            st.subheader(
                "7. Salt Bridges (Distance ≤ 4.0 Å)"
            )

            if not salt_bridge_df.empty:

                st.dataframe(
                    salt_bridge_df,
                    use_container_width=True,
                    hide_index=True
                )

            else:

                st.info(
                    "No potential salt bridges were detected "
                    "within 4.0 Å."
                )


            # =================================================
            # STEP 15: GENERATE EXCEL REPORT
            # =================================================

            output_name = (

                os.path.splitext(pdb_file)[0]

                + "_PPI_Analysis.xlsx"

            )


            excel_buffer = io.BytesIO()


            with pd.ExcelWriter(

                excel_buffer,

                engine="openpyxl"

            ) as writer:


                interaction_df.to_excel(

                    writer,

                    sheet_name="Interacting_Residues",

                    index=False

                )


                salt_bridge_df.to_excel(

                    writer,

                    sheet_name="Salt_Bridges",

                    index=False

                )


                # Excel formatting

                for sheet in writer.sheets.values():

                    sheet.freeze_panes = "A2"

                    sheet.auto_filter.ref = sheet.dimensions


                    for column_cells in sheet.columns:

                        max_length = max(

                            len(str(cell.value))

                            if cell.value is not None

                            else 0

                            for cell in column_cells

                        )


                        sheet.column_dimensions[

                            column_cells[0].column_letter

                        ].width = min(max_length + 3, 35)


            excel_buffer.seek(0)


            # =================================================
            # STEP 16: DOWNLOAD EXCEL REPORT
            # =================================================

            st.divider()

            st.subheader("8. Download Excel Report")

            st.download_button(

                label="📥 Download PPI Analysis Report (.xlsx)",

                data=excel_buffer.getvalue(),

                file_name=output_name,

                mime=(
                    "application/vnd.openxmlformats-officedocument."
                    "spreadsheetml.sheet"
                ),

                type="primary",

                use_container_width=True

            )


            st.success(
                f"Excel report ready: {output_name}"
            )


    except Exception as e:

        st.error(
            f"An error occurred during PDB analysis: {str(e)}"
        )


else:

    st.info(
        "Please upload a .pdb file to begin protein interaction analysis."
    )


# ============================================================
# END OF APPLICATION
# ============================================================
