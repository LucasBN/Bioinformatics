from coronavirus_sequence import *
import random
from DNAToolKit import *
from utilities import coloured
import colorama

DNA_Sequence = validateSequence(coronavirus)
readingFrames = generateReadingFrames(DNA_Sequence)
proteins = (convertAminoAcidsToProteins(readingFrames))


# Create a file containing all proteins in the DNA_Sequence
with open('proteins/proteins.txt', 'w') as f:
    for protein in proteins:
        f.write("%s\n" % protein)
        pass

res = (str(locateProteinMotif(proteins[293], [ 'N', '{P}', '[ST]', '{P}' ])).replace(',', '').replace('[', '').replace(']', ''))
print(res)
