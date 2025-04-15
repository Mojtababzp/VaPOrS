# VaPOrS: Vapor Pressure & Enthalpy of Vaporization Calculation from SMILES  

A Python-based tool that takes a list of compounds in SMILES format, detects the presence and quantity of 30 structural groups based on the SIMPOL method and calculates their vapor pressure and enthalpy of vaporization.

---

## **Installation & Setup**  

1. Clone the Repository
First, download the repository to your local machine:  
```sh
git clone https://github.com/Mojtababzp/VaPOrS.git  
cd VaPOrS

2. Install Dependencies
Ensure Python (≥3.8) is installed, then install the required libraries:
```sh
pip install numpy scipy

3. Run the Code
- Option A: Use Jupyter Notebook
To launch Jupyter Notebook and open VaPOrS.ipynb, run:
```sh
jupyter notebook  
Then navigate to VaPOrS.ipynb and execute all cells.
- Option B: Use the Python Script
To run the code directly via terminal using VaPOrS.py, simply run:
```sh
python VaPOrS.py

Make sure SMILES.txt is present in the directory before execution.

Input File (SMILES.txt)
The input file already exists in the repository with several SMILES in each row as examples.
Edit this file to insert your own SMILES notations (one SMILES per line, no headers).

Output Files
After execution, the following output files will be created or updated in the same directory:
output.csv: Tabular output containing structural groups, vapor pressure (Pa), and enthalpy of vaporization (kJ/mol).
output.txt: Text output summarizing the same information in a readable format.
These files already exist in the directory as an example output for guidance.

Citation
If you use this tool in your research, please cite:
Mojtaba Bezaatpour. (2025). Mojtababzp/VaPOrS: VaPOrS v1.0.1 (v.1.0.1). Zenodo. https://doi.org/10.5281/zenodo.15222175


