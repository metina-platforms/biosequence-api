from pydantic import BaseModel, Field
from typing import Optional, Literal
from enum import Enum


class SequenceType(str, Enum):
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "protein"


class SequenceRequest(BaseModel):
    sequence: str = Field(..., description="The biological sequence to analyze")
    sequence_type: SequenceType = Field(default=SequenceType.DNA, description="Type of sequence")


class SequenceResponse(BaseModel):
    sequence: str = Field(..., description="The resulting sequence")


class SequenceAnalysis(BaseModel):
    sequence: str = Field(..., description="Original sequence")
    sequence_type: SequenceType = Field(..., description="Type of sequence")
    length: int = Field(..., description="Length of sequence")
    gc_content: float = Field(..., description="GC content percentage")
    molecular_weight: Optional[float] = Field(None, description="Molecular weight")
    composition: dict = Field(..., description="Nucleotide/amino acid composition")