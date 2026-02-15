from fastapi import FastAPI
from rdkit import Chem
from rdkit.Chem import Descriptors
import uuid

from app.schemas import AnalyzeRequest, AnalyzeResponse, DescriptorResult

app = FastAPI()

@app.get("/")
async def root():
    return {"ok": True}

@app.get("/health")
async def health():
    return {"status": "healthy", "apiVersion": "1.0"}

def _largest_fragment(mol: Chem.Mol) -> Chem.Mol:
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    return mol

def _calc(mol: Chem.Mol) -> DescriptorResult:
    return DescriptorResult(
        mw=round(Descriptors.MolWt(mol), 4),
        logp=round(Descriptors.MolLogP(mol), 4),
        hbd=Descriptors.NumHDonors(mol),
        hba=Descriptors.NumHAcceptors(mol),
        tpsa=round(Descriptors.TPSA(mol), 4),
        atom_count=mol.GetNumAtoms(),
    )

@app.post("/analyze", response_model=AnalyzeResponse)
async def analyze_post(req: AnalyzeRequest):
    rid = str(uuid.uuid4())

    mol = Chem.MolFromSmiles(req.smiles)
    if mol is None:
        return AnalyzeResponse(
            requestId=rid,
            smiles=req.smiles,
            valid=False,
            error=f"Invalid SMILES: {req.smiles}",
        )

    mol = _largest_fragment(mol)
    canonical = Chem.MolToSmiles(mol)
    desc = _calc(mol)

    return AnalyzeResponse(
        requestId=rid,
        smiles=req.smiles,
        canonicalSmiles=canonical,
        valid=True,
        descriptors=desc,
        error=None,
    )

@app.get("/analyze", response_model=AnalyzeResponse)
async def analyze_get(smiles: str):
    return await analyze_post(AnalyzeRequest(smiles=smiles))
