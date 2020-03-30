import random
from DNAToolKit import *
from utilities import coloured
import colorama
from fasta import *

# Used to colour output in terminal
colorama.init()

# Generate a random DNA sequence that is 50 nucleotides long
DNA_Sequence = ''.join([random.choice(nucleotides) for nucleotide in range(30)])
Mutated_DNA_Sequence = ''.join([random.choice(nucleotides) for nucleotide in range(30)])

# Validate and normalise the DNA sequence
DNA_Sequence = validateSequence(DNA_Sequence)

print('\n(1) DNA Sequence: ' + coloured(DNA_Sequence))
print('(2) DNA Sequence Length: ' + str(len(DNA_Sequence)))
print('(3) Nucleotide Frequency: ' + coloured(str(countNucleotideFrequency(DNA_Sequence))))
print('(4) DNA/RNA Transcription: ' + coloured(transcription(DNA_Sequence)))
print('(5) DNA Reverse Complement: ' + coloured(complement(DNA_Sequence)))
print('(6) GC Content Percentage: ' + str(calculateGCContent(DNA_Sequence)))
print('(7) GC Content Subsection k = 5: ' + str(calculateGCContentSubsection(DNA_Sequence, k=5)))
print('(8) Hamming Distance: ' + str(countPointMutations(DNA_Sequence, Mutated_DNA_Sequence)))
print('(9) Translated RNA: ' + translateDNA(DNA_Sequence, 0))
print('(10) Codons Used: ' + str(codonUsageDNA(DNA_Sequence, '*')))
print('(11) Reading Frames: ' + str(generateReadingFrames(DNA_Sequence)))
print('(12) Proteins: ' + str(convertAminoAcidsToProteins(generateReadingFrames(DNA_Sequence))))

print('\n================================ END ================================')

x = """GTGTAGCTCTCATAGTCGGTTTGCGTTGAACACTCCTGAGACCAAAGATGAATAGCCACG
AAGCCTGATTAAGCCTGTAAGGTAAAGCCGTGATAAGTTCTAACCCAGCAGAGGAACAGT
GGTACTACAGAGTCTATGATGGATGGCACGGTTAGATTTAATTCACTGTACACAATGATA
TGGCTTTGTGAAAAACTTACCTAATCCTCCGGTCCCTCATATTTTGGCGCCACGCCTAGT
TCACTATTGGCAAATGCCCTTACTGAAAGCATTGCCGGAGTGCGTGCTCTGAGAGCTCCC
GTAAAGTATCGCATAATCAACGGATACCATACCAAGGCGCCTGGCGGGTAATGCGTTGCG
TACCACAGTGCGCACCGCGCGAAGCCAAGGAACTACTAGCCTACGTCAGGGCGCCTGCTC
CGGGCAAATCACGAACGGCTGCTGTCTTACTATGTATAAATGAAAGGCATGAATGAACAA
TAAGCTTCCCTATATTCACCGTTGGACTTCGGGACATCCTAGGCAGACACTATTCTACTC
GGCTAAGCCTCTTAAGAGGTAGTAATTTAAGTACTCGTTGCAATGAAATATCACAGTAAT
ACAGATGCCACTCATTTCGTCCATCCCGAACAAGCGCCATCTGAGGGAGGTACGCCGTCG
TACGTATGTGCCTCTCTGACTAGCGCCCTAGGGACTAACGAAGAACTGGCCAGTTGCTTT
ACGGTATTGAATGGTTGACAGCCGGCTGCCCCACTTAGTATGGACTGGCATGCGCCAATT
GTGTCCCAGGCGTCCGCGTTGAACGGAATGATACGCGCGACGATCGTCGATGCCTAACTT
AATTCGGACCAACGAGGGCTGGAGTGCTAATGGGATGAGTTATTGCCGTTGATTGACCTG
TTCACGGTCCAAGAGCGTGGAGTTGTAAACTTTTTGTTGACTT""".replace('\n', '')
y = """GCACAATTAATCCGTAGGTGCCTCCCGCCAGGCATGTAAGTTCTAGAGAAGCTACAGGTAAGTGCCGCCTGGGAGCTTCCGCTACACATCTTACGAGTGG""".replace('\n', '')

print(str(locateSubseqeunce(x, y)).replace(',', '').replace('[', '').replace(']', ''))
