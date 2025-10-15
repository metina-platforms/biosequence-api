import re
from typing import List
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from app.models.annotation import AnnotationResponse


class AnnotationService:
    
    def find_orfs(self, sequence: str, min_length: int = 300) -> List[AnnotationResponse]:
        """Find Open Reading Frames in a DNA sequence"""
        sequence = sequence.upper()
        orfs = []
        
        # Define start and stop codons
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        # Check all 6 reading frames (3 forward, 3 reverse)
        for frame in range(3):
            # Forward frames (1, 2, 3)
            orfs.extend(self._find_orfs_in_frame(
                sequence, frame, frame + 1, start_codons, stop_codons, min_length
            ))
            
            # Reverse frames (-1, -2, -3)
            reverse_seq = str(Seq(sequence).reverse_complement())
            orfs.extend(self._find_orfs_in_frame(
                reverse_seq, frame, -(frame + 1), start_codons, stop_codons, min_length
            ))
        
        return orfs
    
    def _find_orfs_in_frame(self, sequence: str, offset: int, frame: int, 
                           start_codons: List[str], stop_codons: List[str], 
                           min_length: int) -> List[AnnotationResponse]:
        """Find ORFs in a specific reading frame"""
        orfs = []
        seq_len = len(sequence)
        
        i = offset
        while i < seq_len - 2:
            # Look for start codon
            codon = sequence[i:i+3]
            if codon in start_codons:
                start_pos = i
                j = i + 3
                
                # Look for stop codon
                while j < seq_len - 2:
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf_length = j + 3 - start_pos
                        if orf_length >= min_length:
                            orfs.append(AnnotationResponse(
                                start=start_pos,
                                end=j + 3,
                                frame=frame,
                                sequence=sequence[start_pos:j+3],
                                length=orf_length,
                                annotation_type="ORF"
                            ))
                        break
                    j += 3
                i = j
            else:
                i += 3
        
        return orfs
    
    def find_motifs(self, sequence: str) -> List[dict]:
        """Find common sequence motifs"""
        sequence = sequence.upper()
        motifs = []
        
        # Common motifs to search for
        motif_patterns = {
            "TATA_box": r"TATAAA",
            "CAAT_box": r"CCAAT",
            "GC_box": r"GGGCGG",
            "Kozak_sequence": r"[AG]CC[AG]CCAUGG",
            "Poly_A_signal": r"AAUAAA"
        }
        
        for motif_name, pattern in motif_patterns.items():
            matches = list(re.finditer(pattern, sequence))
            for match in matches:
                motifs.append({
                    "motif_type": motif_name,
                    "start": match.start(),
                    "end": match.end(),
                    "sequence": match.group()
                })
        
        return motifs
    
    def calculate_gc_content_window(self, sequence: str, window_size: int = 100) -> List[dict]:
        """Calculate GC content using sliding window"""
        sequence = sequence.upper()
        gc_contents = []
        
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            window = sequence[i:i + window_size]
            if len(window) == window_size:
                gc_content = GC(Seq(window))
                gc_contents.append({
                    "start": i,
                    "end": i + window_size,
                    "gc_content": round(gc_content, 2)
                })
        
        return gc_contents