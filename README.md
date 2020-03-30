## validateSequence(sequence):
Ensures that the DNA sequence is valid and does not contain any unexpected characters

## countNucleotideFrequency(sequence):
Returns a dictionary that contains the number of time each nucleotide appears in the DNA Sequence

## transcription(sequence):
Transcribes the DNA into RNA by replacing Thymine (T) with Uracil (U)

## complement(sequence):
Returns the DNA Reverse Complement of a DNA sequence

## calculateGCContent(sequence):
Returns the percentage compositon by C/G of a DNA sequence

## calculateGCContentSubsection(sequence, k=20):
Returns the percentage compositon by C/G of a DNA sequence within subsections of length k

## countPointMutations(sequence, mutatedSequence):
Returns the hamming distance between two DNA sequences

## translateRNA(sequence, initialPosition=0):
Returns the amino acid chain translated from a RNA sequence

## translateDNA(sequence, initialPosition=0):
Returns the amino acid chain translated from a DNA sequence

## codonUsageDNA(sequence, aminoAcid):
Returns the percentage compositon of each codon used in a DNA sequence

## generateReadingFrames(sequence):
Returns all 6 reading frames (in an array) of a DNA sequence

## convertAminoAcidsToProteins(readingFrames):
Returns a list of all possible proteins that can be made from all 6 reading frames

## locateMotif(sequence, motif):
Returns a list of all motifs within a single DNA sequence

## findSharedMotifs(sequences):
Returns a list of all shared motifs in an array of DNA sequences (slow)

## motifInSequences(motif, sequences):
Returns True if motif is present in all DNA sequences

## findLongestSharedMotifs(sequences):
Returns the longest shared motif in an array of DNA sequences

## locateProteinMotif(proteinString, motif):
Returns a list of the positions of a protein motif within a protein string

## downloadProtein(proteinID):
Returns a string containing the amino acid chain of a given protein (Data downloaded from https://www.uniprot.org/)
