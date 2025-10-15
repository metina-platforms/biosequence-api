from fastapi import APIRouter, HTTPException
from app.models.annotation import AnnotationRequest, AnnotationResponse
from app.services.annotation_service import AnnotationService

router = APIRouter()
annotation_service = AnnotationService()


@router.post("/find-orfs", response_model=list[AnnotationResponse])
async def find_open_reading_frames(request: AnnotationRequest):
    """
    Find Open Reading Frames (ORFs) in a DNA sequence
    """
    try:
        orfs = annotation_service.find_orfs(request.sequence, request.min_length or 300)
        return orfs
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/find-motifs")
async def find_motifs(request: AnnotationRequest):
    """
    Find common sequence motifs in a DNA/RNA sequence
    """
    try:
        motifs = annotation_service.find_motifs(request.sequence)
        return {"motifs": motifs}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/gc-content")
async def calculate_gc_content(request: AnnotationRequest):
    """
    Calculate GC content across the sequence with sliding window
    """
    try:
        gc_content = annotation_service.calculate_gc_content_window(
            request.sequence, request.window_size or 100
        )
        return {"gc_content": gc_content}
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
