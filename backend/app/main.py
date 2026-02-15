from fastapi import FastAPI
from rdkit import Chem
from rdkit.Chem import Descriptors
from app.schemas import AnalyzeRequest, AnalyzeResponse, DescriptorResult
import uuid

app = FastAPI()

@app.get("/")
def root():
    return {"ok": True}

@app.post("/analyze", response_model=AnalyzeResponse)
def analyze(request: AnalyzeRequest):

    request_id = str(uuid.uuid4())

    mol = Chem.MolFromSmiles(request.smiles)

    if mol is None:
        return AnalyzeResponse(
            requestId=request_id,
            smiles=request.smiles,
            valid=False,
            error=f"Invalid SMILES: {request.smiles}"
        )

    # Salt handling
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())

    canonical = Chem.MolToSmiles(mol)

    descriptors = DescriptorResult(
        mw=round(Descriptors.MolWt(mol), 4),
        logp=round(Descriptors.MolLogP(mol), 4),
        hbd=Descriptors.NumHDonors(mol),
        hba=Descriptors.NumHAcceptors(mol),
        tpsa=round(Descriptors.TPSA(mol), 4),
        atom_count=mol.GetNumAtoms()
    )

    return AnalyzeResponse(
        requestId=request_id,
        smiles=request.smiles,
        canonicalSmiles=canonical,
        valid=True,
        descriptors=descriptors
    )
