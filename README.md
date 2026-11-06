# ML-Guided Virtual Screening of Phytochemicals Against SARS-CoV-2 Mpro

*EthnoHack 2026 — Track 2: Natural Products & Drug Discovery*

**Project overview**

This project shows an integrated computational pipeline which combines **Machine Learning (ML)** and **Molecular Docking** to rapidly screen traditional plant compounds against the SARS-CoV-2 Main Protease (Mpro). 

The main purpose of this project was to assess which phytochemicals from the USDA Dr. Duke's database are most likely to inhibit the SARS-CoV-2 Main Protease (Mpro) virus, and then validating the top predictions of ML through docking using AutoDock Vina.


**The Problem Addressed**
covid-19 is a potential viral disease that spreads rapidly so new antiviral compounds must be identified which are effective against this disease to prevent its spread. However, traditional drug discovery from medicinal plants is time consuming and required high cost. Testing thousands of compounds in a wet lab takes months. So our question was: *Can computer be trained  to identify Mpro inhibitors and prioritize candidates for wet-lab testing?*

 **Methodology & Workflow**
1. **Data Cleaning:** we targeted USDA Dr. Duke's Phytochemical Database consisting of 400 plant compounds. compounds lacking chemical structures were removed, leaving **293 valid compounds**.
2. **Featurization:** **RDKit** was used to Convert molecular SMILES strings to 1024-bit Morgan Fingerprints.
3. **Machine Learning:**  **Random Forest Classifier** (Scikit-Learn) model was trained using known Mpro inhibitors and random molecules to predict inhibition probabilities.
4. **Molecular Docking:** top 10 ML hits were Docked against the Mpro crystal structure (**PDB: 6LU7**) using **AutoDock Vina** (via SwissDock).  we used Quercetin  as a positive control.

**Key Results**
| Compound | ML Probability | Docking Score (kcal/mol) |
| :--- | :--- | :--- |
| **Apigenin** (ML #1) | 70% | -5.377 |
| **Quercetin** (Control) | 59% | -5.423 |
| **Baicalin** (ML #2) | 69% | **-5.725** *(Strongest Binder)* |

*   **Apigenin** validated successfully, scoring almost identically to the known inhibitor Quercetin.
*   **Baicalin** (the ML's second choice) physically bound tighter to the viral protein than the known control.
*   **Betavulgarin** (ML #3) scored poorly (-2.239), indicating the necessity of docking to filter out ML false positives.



 **Core Code Snippets**
The complete code pipeline is available in the `EthnoHack_ML_Pipeline.ipynb` Jupyter Notebook.


1. **Generating Morgan Fingerprints**
```python
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def smiles_to_fingerprint(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is not None:
        # Generate Morgan Fingerprint (radius 2, 1024 bits)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        return np.array(fp)
    return None

print("Generating fingerprints (this might take 10-20 seconds)...")

# Apply the function to our clean dataframe
df_clean['fingerprint'] = df_clean['canonical_smiles'].apply(smiles_to_fingerprint)

# Drop any rows where RDKit failed to read the SMILES
df_clean = df_clean.dropna(subset=['fingerprint'])

print(f"Successfully generated fingerprints for {len(df_clean)} compounds!")
```

**Training the Random Forest Model**

```python
from sklearn.ensemble import RandomForestClassifier

# Generate fingerprints for the training set
def get_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))

X_train = np.array([get_fp(s) for s in train_smiles])
y_train = np.array(train_labels)

# Train the Random Forest model
print("Training the Random Forest Classifier...")
rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)
print("Model trained successfully!")
```

**Repository Contents**

EthnoHack_ML_Pipeline.ipynb: The complete Python code for data processing and ML training.

ethno_sample_400.json: The raw USDA Dr. Duke's dataset.

docking summary.xlsx: A table of ML probabilities and Docking scores.

results_*.zip: Raw output files from SwissDock for all 10 docked compounds.
*6LU7 (PDB)*: Mpro crystal structure as target of our project

**Limitations**

 All results are entirely computational. It doesnot involve any wet lab validation.

 The ML model was trained on a small dataset so it might have low accuracy.

 Toxicity and bioavailability were not evaluated.
