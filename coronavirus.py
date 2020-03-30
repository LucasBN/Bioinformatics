from coronavirus_sequence import *
import random
from DNAToolKit import *
from utilities import coloured
import colorama

CV2_Sequence = validateSequence(coronavirus)
CV1_Sequence = validateSequence(sars)

CV2_readingFrames = generateReadingFrames(CV2_Sequence)
CV2_proteins = (convertAminoAcidsToProteins(CV2_readingFrames))

CV1_readingFrames = generateReadingFrames(CV1_Sequence)
CV1_proteins = (convertAminoAcidsToProteins(CV1_readingFrames))


# Create a file containing all proteins in the DNA_Sequence
with open('proteins/CV2_proteins.txt', 'w') as f:
    for protein in CV2_proteins:
        if len(protein) > 200:
            f.write("%s\n" % protein)
        pass

# Create a file containing all proteins in the DNA_Sequence
with open('proteins/CV1_proteins.txt', 'w') as f:
    for protein in CV1_proteins:
        if len(protein) > 200:
            f.write("%s\n" % protein)
        pass
