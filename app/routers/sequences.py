from fastapi import APIRouter, HTTPException
from app.models.sequence import SequenceRequest, SequenceResponse, SequenceAnalysis
from app.services.sequence_service import SequenceAnalyzer

router = APIRouter()
sequence_analyzer = SequenceAnalyzer()

@router.post("/analyze", response_model=SequenceAnalysis)
async def analyze_sequence(sequence_request: SequenceRequest):
    """
    Analyze a DNA/RNA sequence for basic properties and composition
    """
    try:
        analysis = sequence_analyzer.analyze_sequence(
            sequence_request.sequence, 
            sequence_request.sequence_type
        )
        return analysis
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/reverse-complement", response_model=SequenceResponse)
async def get_reverse_complement(sequence_request: SequenceRequest):
    """
    Get the reverse complement of a DNA sequence
    """
    try:
        result = sequence_analyzer.reverse_complement(sequence_request.sequence)
        return SequenceResponse(sequence=result)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/transcribe", response_model=SequenceResponse)
async def transcribe_dna(sequence_request: SequenceRequest):
    """
    Transcribe DNA sequence to RNA
    """
    try:
        result = sequence_analyzer.transcribe(sequence_request.sequence)
        return SequenceResponse(sequence=result)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/translate", response_model=SequenceResponse)
async def translate_sequence(sequence_request: SequenceRequest):
    """
    Translate DNA/RNA sequence to protein
    """
    try:
        result = sequence_analyzer.translate(sequence_request.sequence)
        return SequenceResponse(sequence=result)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))