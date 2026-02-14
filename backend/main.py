from fastapi import FastAPI
from rdkit import Chem
from rdkit.Chem import Descriptors
from pydantic import BaseModel

app = FastAPI()

@app.get("/")
def root():
    mol = Chem.MolFromSmiles("CCO")
    return {"ok": True, "mol": mol is not None}

@app.get("/analyze")
def analyze(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    return {
        "smiles": smiles,
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": Descriptors.TPSA(mol),
        "num_atoms": mol.GetNumAtoms(),
    }

class AnalyzeRequest(BaseModel):
    smiles: str

@app.post("/analyze")
def analyze_post(req: AnalyzeRequest):
    return analyze(req.smiles)


