# VaPOrS: Vapor Pressure & Enthalpy of Vaporization Calculation from SMILES  

**VaPOrS** is a Jupyter Notebook tool that:  
- **Identifies and quantifies 30 functional groups** in organic compounds from their **SMILES notation**.  
- **Calculates vapor pressure and enthalpy of vaporization** using the **SIMPOL method**.  
- **Outputs results in `output.csv` and `output.txt`**.  

---

## **Installation & Setup**  

1. clone the repository
First, download the repository to your local machine:  
```sh
git clone https://github.com/Mojtababzp/VaPOrS.git  
cd VaPOrS
2. Install Dependencies
Ensure Python (≥3.8) is installed, then install the required libraries:

sh
pip install numpy scipy
 
Alternatively, if running inside Jupyter Notebook, execute:
!pip install numpy scipy

3. Open and Run the Jupyter Notebook
To launch Jupyter Notebook and open VaPOrS.ipynb, run:
jupyter notebook  
Then navigate to VaPOrS.ipynb and execute all cells.

Usage
Input File (SMILES.txt)
The input file already exists in the repository with several SMILES in each row as an example.
Users must edit it  by inserting the SMILES of interest (one SMILES notation per line, without headers).
After running the notebook, the following files will be created (or updated) in the same directory, containing the employed SMILES with the occurrence of functional groups, vapor Pressure (Pa) and Enthalpy of Vaporization (kJ/mol) at the given temperature (298.15 by default) and temperature-based vapor Pressure (log P in atm) and Enthalpy of Vaporization (Delta_h in kJ/mol) for them in each row:
output.csv
output.txt
These files already exist in the directory as the outputs of the example for guidance.
