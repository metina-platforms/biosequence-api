from Bio.Seq import Seq
from Bio.SeqUtils import GC, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from app.models.sequence import SequenceAnalysis, SequenceType
from collections import Counter


class SequenceAnalyzer:

    def analyze_sequence(
        self, sequence: str, sequence_type: SequenceType
    ) -> SequenceAnalysis:
        """Analyze a biological sequence and return comprehensive analysis"""
        seq_obj = Seq(sequence.upper())

        # Calculate basic properties
        length = len(sequence)
        composition = dict(Counter(sequence.upper()))

        # Calculate GC content for DNA/RNA
        gc_content = 0.0
        mol_weight = None

        if sequence_type in [SequenceType.DNA, SequenceType.RNA]:
            gc_content = GC(seq_obj)
            mol_weight = molecular_weight(seq_obj, seq_type=sequence_type.value)
        elif sequence_type == SequenceType.PROTEIN:
            # For proteins, use ProteinAnalysis
            protein_analysis = ProteinAnalysis(sequence.upper())
            mol_weight = protein_analysis.molecular_weight()

        return SequenceAnalysis(
            sequence=sequence.upper(),
            sequence_type=sequence_type,
            length=length,
            gc_content=gc_content,
            molecular_weight=mol_weight,
            composition=composition,
        )

    def reverse_complement(self, sequence: str) -> str:
        """Get reverse complement of DNA sequence"""
        seq_obj = Seq(sequence.upper())
        return str(seq_obj.reverse_complement())

    def transcribe(self, sequence: str) -> str:
        """Transcribe DNA to RNA"""
        seq_obj = Seq(sequence.upper())
        return str(seq_obj.transcribe())

    def translate(self, sequence: str) -> str:
        """Translate DNA/RNA to protein"""
        seq_obj = Seq(sequence.upper())
        return str(seq_obj.translate())
