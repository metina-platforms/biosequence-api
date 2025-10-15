from pydantic import BaseModel, Field
from typing import Optional, List


class AnnotationRequest(BaseModel):
    sequence: str = Field(..., description="The DNA/RNA sequence to analyze")
    min_length: Optional[int] = Field(None, description="Minimum length for ORF detection")
    window_size: Optional[int] = Field(None, description="Window size for sliding window analysis")


class AnnotationResponse(BaseModel):
    start: int = Field(..., description="Start position of the annotation")
    end: int = Field(..., description="End position of the annotation")
    frame: int = Field(..., description="Reading frame (1, 2, 3, -1, -2, -3)")
    sequence: str = Field(..., description="The annotated sequence")
    length: int = Field(..., description="Length of the annotation")
    annotation_type: str = Field(..., description="Type of annotation (e.g., ORF, motif)")