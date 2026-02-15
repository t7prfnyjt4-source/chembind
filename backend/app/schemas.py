from pydantic import BaseModel, Field, field_validator
from typing import Optional
import uuid

class AnalyzeRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)

    @field_validator("smiles")
    @classmethod
    def strip_and_reject_blank(cls, v: str) -> str:
        v = v.strip()
        if v == "":
            raise ValueError("SMILES cannot be empty")
        return v

class DescriptorResult(BaseModel):
    mw: float
    logp: float
    hbd: int
    hba: int
    tpsa: float
    atom_count: int

class AnalyzeResponse(BaseModel):
    apiVersion: str = "1.0"
    requestId: str = Field(default_factory=lambda: str(uuid.uuid4()))
    smiles: str
    canonicalSmiles: Optional[str] = None
    valid: bool
    descriptors: Optional[DescriptorResult] = None
    error: Optional[str] = None
