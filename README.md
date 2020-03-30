# Ensures that the DNA sequence is valid and does not contain any unexpected characters
def validateSequence(sequence):

# Returns a dictionary that contains the number of time each nucleotide appears in the DNA Sequence
def countNucleotideFrequency(sequence):

# Transcribes the DNA into RNA by replacing Thymine (T) with Uracil (U)
def transcription(sequence):

# Returns the DNA Reverse Complement of a DNA sequence
def complement(sequence):

# Returns the percentage compositon by C/G of a DNA sequence
def calculateGCContent(sequence):

# Returns the percentage compositon by C/G of a DNA sequence within subsections of length k
def calculateGCContentSubsection(sequence, k=20):

# Returns the hamming distance between two DNA sequences
def countPointMutations(sequence, mutatedSequence):

# Returns the amino acid chain translated from a RNA sequence
def translateRNA(sequence, initialPosition=0):

# Returns the amino acid chain translated from a DNA sequence
def translateDNA(sequence, initialPosition=0):

# Returns the percentage compositon of each codon used in a DNA sequence
def codonUsageDNA(sequence, aminoAcid):

# Returns all 6 reading frames (in an array) of a DNA sequence
def generateReadingFrames(sequence):

# Returns a list of all possible proteins that can be made from all 6 reading frames
def convertAminoAcidsToProteins(readingFrames):

# Returns a list of all motifs within a single DNA sequence
def locateMotif(sequence, motif):

# Returns a list of all shared motifs in an array of DNA sequences (slow)
def findSharedMotifs(sequences):

# Returns True if motif is present in all DNA sequences
def motifInSequences(motif, sequences):

# Returns the longest shared motif in an array of DNA sequences
def findLongestSharedMotifs(sequences):

# Returns a list of the positions of a protein motif within a protein string
def locateProteinMotif(proteinString, motif):

# Returns a string containing the amino acid chain of a given protein (Data downloaded from https://www.uniprot.org/)
def downloadProtein(proteinID):
