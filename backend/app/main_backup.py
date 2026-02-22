from fastapi import FastAPI, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import PlainTextResponse
import traceback
import uuid
import time

from rdkit import Chem
from rdkit.Chem import Descriptors

from app.schemas import AnalyzeRequest, AnalyzeResponse, DescriptorResult
from app.firebase_init import initialize_firebase
from app.auth import require_auth
from app.firestore_store import save_analysis, list_analyses, get_analysis

# ✅ Initialize Firebase BEFORE handling requests
initialize_firebase()

app = FastAPI()

# ✅ Show real stack trace during development
@app.exception_handler(Exception)
async def all_exception_handler(request, exc):
    return PlainTextResponse(traceback.format_exc(), status_code=500)

# ---- CORS ----
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://localhost:5174",
        "http://127.0.0.1:5173",
        "http://127.0.0.1:5174",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    return {"ok": True}

@app.get("/health")
async def health():
    return {"status": "healthy", "apiVersion": "1.0"}

@app.get("/me")
def me(uid: str = Depends(require_auth)):
    return {"uid": uid}

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
async def analyze_post(req: AnalyzeRequest, uid: str = Depends(require_auth)):
    rid = str(uuid.uuid4())

    mol = Chem.MolFromSmiles(req.smiles)
    if mol is None:
        return AnalyzeResponse(
            requestId=rid,
            smiles=req.smiles,
            valid=False,
            error=f"Invalid SMILES: {req.smiles}",
        )

    try:
        mol = _largest_fragment(mol)
        canonical = Chem.MolToSmiles(mol)
        desc = _calc(mol)

        resp = AnalyzeResponse(
            requestId=rid,
            smiles=req.smiles,
            canonicalSmiles=canonical,
            valid=True,
            descriptors=desc,
            error=None,
        )

        # 🔒 Safe Firestore save (will NOT crash endpoint)
        try:
            save_analysis(uid, {
                "requestId": resp.requestId,
                "uid": uid,
                "smiles": resp.smiles,
                "canonicalSmiles": resp.canonicalSmiles,
                "valid": resp.valid,
                "descriptors": resp.descriptors.model_dump(),
                "error": resp.error,
                "createdAt": int(time.time()),
            })
        except Exception as e:
            print("Firestore save failed:", repr(e))

        return resp

    except Exception:
        print("ANALYZE CRASH:\n", traceback.format_exc())
        raise

@app.get("/analyses")
def analyses(uid: str = Depends(require_auth)):
    return {"items": list_analyses(uid, limit=25)}

@app.get("/analyses/{request_id}")
def analyses_one(request_id: str, uid: str = Depends(require_auth)):
    doc = get_analysis(uid, request_id)
    return {"item": doc}

