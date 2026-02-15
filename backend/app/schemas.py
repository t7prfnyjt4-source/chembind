from pydantic import BaseModel, Field
from typing import Optional
import uuid

class AnalyzeRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)

class DescriptorResult(BaseModel):
    mw: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    atom_count: int

class AnalyzeResponse(BaseModel):
    apiVersion: str = "1.0"
    requestId: str
    smiles: str
    canonicalSmiles: Optional[str] = None
    valid: bool
    descriptors: Optional[DescriptorResult] = None
    error: Optional[str] = None
