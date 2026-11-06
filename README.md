# Protein-Protein Interaction Analyzer

A Streamlit-based web application for analyzing protein-protein interactions from a user-uploaded PDB structure file.

The app reads a `.pdb` file, identifies protein chains, allows the user to select two chains for comparison, and reports:

- atom-pair contacts within 3.5 Å
- potential salt bridges within 4.0 Å
- residue and atom-level interaction details
- an Excel report containing both result tables

## Features

- Upload a protein structure in `.pdb` format
- Detects amino-acid-containing chains automatically
- Lets users select any two chains for interaction analysis
- Calculates interchain atomic distances using heavy atoms only
- Identifies potential salt bridges involving:
  - positively charged residues: ARG, LYS, HIS
  - negatively charged residues: ASP, GLU
- Displays summary metrics and detailed data tables in the app
- Exports the results as an `.xlsx` file with two sheets:
  - `Interacting_Residues`
  - `Salt_Bridges`

## Requirements

Install the required Python packages:

```bash
pip install streamlit pandas numpy biopython openpyxl
```

## Run the Application

From the project directory, run:

```bash
streamlit run app.py
```

Then open the local URL shown in the terminal in your browser.

## How It Works

1. Upload a valid `.pdb` file.
2. The app parses the structure using Biopython.
3. It identifies amino-acid chains in the structure.
4. You select two chains for interaction analysis.
5. The app calculates distances between heavy atoms from each chain.
6. It reports all atom pairs within 3.5 Å and potential salt bridges within 4.0 Å.
7. You can download a detailed Excel report.

## Output

The generated Excel workbook includes:

### Interacting_Residues
Contains the residue-level atom interactions between the selected chains, including:

- chain IDs
- residue names and numbers
- atom names
- distance in Å

### Salt_Bridges
Contains potential salt-bridge interactions, including:

- positive and negative chain IDs
- residue names and residue numbers
- involved atom names
- distance in Å
- interaction type

## Notes

- The uploaded file must have a lowercase `.pdb` extension.
- At least two protein chains are required for analysis.
- Hydrogen atoms are excluded from the interaction calculations.
- This tool is intended for protein structure analysis and is best suited for static PDB coordinate inspection.

## Example Workflow

```bash
git clone https://github.com/tirthabrata0407-cloud/Protein-Protein-Interaction-tool.git
cd Protein-Protein-Interaction-tool
pip install -r requirements.txt
streamlit run app.py
```

If a `requirements.txt` file is not present in the repository, install the libraries listed in the Requirements section manually.

## License

This project is provided as-is for educational and research-oriented structural biology analysis.
