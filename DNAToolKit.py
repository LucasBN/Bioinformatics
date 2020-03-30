import collections
from structures import *
from urllib.request import urlopen

# Ensures that the DNA sequence is valid and does not contain any unexpected characters
def validateSequence(sequence):
    for nucleotide in sequence.upper():
        if nucleotide not in nucleotides:
            return False
    return sequence.upper()

# Returns a dictionary that contains the number of time each nucleotide appears in the DNA Sequence
def countNucleotideFrequency(sequence):
    return dict(collections.Counter(sequence))

# Transcribes the DNA into RNA by replacing Thymine (T) with Uracil (U)
def transcription(sequence):
    return sequence.replace('T', 'U')

# Returns the DNA Reverse Complement of a DNA sequence
def complement(sequence):
    mapping = str.maketrans('ATCG', 'TAGC')
    return sequence.translate(mapping)[::-1]

# Returns the percentage compositon by C/G of a DNA sequence
def calculateGCContent(sequence):
    return round((sequence.count('C') + sequence.count('G')) / len(sequence) * 100)

# Returns the percentage compositon by C/G of a DNA sequence within subsections of length k
def calculateGCContentSubsection(sequence, k=20):
    result = []
    for i in range(0, len(sequence) - k + 1, k):
        subSequence = sequence[i:i+k]
        result.append(calculateGCContent(subSequence))
    return result

# Returns the hamming distance between two DNA sequences
def countPointMutations(sequence, mutatedSequence):
    c = 0
    for i in range(len(sequence)):
        if sequence[i] != mutatedSequence[i]:
            c +=1
    return c

# Returns the amino acid chain translated from a RNA sequence
def translateRNA(sequence, initialPosition=0):
    proteinString = ''
    for i in range(initialPosition, len(sequence) - 2, 3):
        proteinString += (RNACodonTable[str(sequence[i] + sequence[i+1] + sequence[i+2])])
    return proteinString

# Returns the amino acid chain translated from a DNA sequence
def translateDNA(sequence, initialPosition=0):
    proteinString = ''
    for i in range(initialPosition, len(sequence) - 2, 3):
        proteinString += (DNACodonTable[str(sequence[i] + sequence[i+1] + sequence[i+2])])
    return proteinString

# Returns the percentage compositon of each codon used in a DNA sequence
def codonUsageDNA(sequence, aminoAcid):
    tmpList = []
    for i in range(0, len(sequence) - 2, 3):
        if DNACodonTable[sequence[i:i+3]] == aminoAcid or aminoAcid == '*':
            tmpList.append(sequence[i:i+3])
    frequencyDict = dict(collections.Counter(tmpList))
    totalWight = sum(frequencyDict.values())
    for sequence in frequencyDict:
        frequencyDict[sequence] = round(frequencyDict[sequence] / totalWight, 2)
    return frequencyDict

# Returns all 6 reading frames (in an array) of a DNA sequence
def generateReadingFrames(sequence):
    readingFrames = []
    for i in range(3):
        readingFrames.append(translateDNA(sequence, i))
        readingFrames.append(translateDNA(complement(sequence), i))
    return readingFrames

# Returns a list of all possible proteins that can be made from all 6 reading frames
def convertAminoAcidsToProteins(readingFrames):
    aminoAcidChains = []
    for frame in readingFrames:
        for i in range(len(frame)):
            if frame[i] == 'M':
                for j in range(i, len(frame)):
                    if frame[j] == '_':
                        aminoAcidChains.append(frame[i:j])
                        break

    return aminoAcidChains

# Returns a list of all motifs within a single DNA sequence
def locateMotif(sequence, motif):
    results = []
    for i in range(0, len(sequence) - len(motif) + 1, 1):
        tmpStr = ""
        for j in range(len(motif)):
            tmpStr += sequence[i+j]
        if tmpStr == motif: results.append(i+1)
    return results

# Returns a list of all shared motifs in an array of DNA sequences (slow)
def findSharedMotifs(sequences):
    # Take the shortest sequence first - no motif can be longer than this sequence
    # Start with searching for motifs of length one, using each character from the shortest sequence
    # a tmp sequence variable can be used so .replace('MOTIF', '') will avoid any unnecessary repeats
    # increase the length of the motif until the motif is the size of the shortest sequence
    tmpSequences = sorted(sequences, key=len)
    possibleMotifs = []
    sharedMotifs = []
    for i in range(1, len(tmpSequences[0])+1): # length of motif
        for j in range(len(tmpSequences[0]) - i + 1): # motif start point
            if tmpSequences[0][j:j+i] != '': possibleMotifs.append(tmpSequences[0][j:j+i])
    possibleMotifs = sorted(possibleMotifs, key=len)[::-1]
    for motif in possibleMotifs:
        motifInSequences = 0
        for sequence in sequences:
            if [] == locateMotif(sequence, motif):
                break
            else:
                motifInSequences += 1
        if motifInSequences == len(sequences):
            sharedMotifs.append(motif)

    return list(dict.fromkeys(sharedMotifs))

# Returns True if motif is present in all DNA sequences
def motifInSequences(motif, sequences):
    for sequence in sequences:
        if motif not in sequence:
            return False
    return True

# Returns the longest shared motif in an array of DNA sequences
def findLongestSharedMotifs(sequences):
    # Sort sequences
    # Take the shortest sequence
    # First motif: size of entire sequence
    # Nth motif: size of entire sequence - n
    # After finding each potential motif candidate, check to see if it is in the other sequences
    # If it is, stop the program and return that motif as your longest shared motif
    # If not, carry on until program finishes
    tmpSequence = (sequences)[0]

    longest = ''

    for start in range(len(tmpSequence)):
        for end in range(len(tmpSequence), start, -1):

            if (end - start) <= len(longest):
                break

            if motifInSequences(tmpSequence[start:end], sequences):
                    longest = tmpSequence[start:end]

    return longest

# Returns a list of the positions of a protein motif within a protein string
def locateProteinMotif(proteinString, motif):

    positions = []

    for i in range(len(proteinString)-len(motif)+1):
        counter = 0
        for j in range(len(motif)):
            if motif[j][0] == '[':
                for k in range(1, len(motif[j])-1):
                    if proteinString[i+j] == motif[j][k]:
                        counter += 1
            elif motif[j][0] == '{':
                for k in range(1, len(motif[j])-1):
                    if proteinString[i+j] != motif[j][k]:
                        counter += 1
            else:
                if proteinString[i+j] == motif[j]:
                    counter += 1
        if counter == 4:
            positions.append(i+1)

    return positions

# Returns a string containing the amino acid chain of a given protein (Data downloaded from https://www.uniprot.org/)
def downloadProtein(proteinID):
    u = (f"https://www.uniprot.org/uniprot/{proteinID}.fasta")
    proteinString = ""
    with urlopen(u) as page:
      data = []
      for line in page:
          data.append(line.decode("utf-8").replace('\n', ''))
      proteinString = (''.join(data[1::]))
    return proteinString

# Returns a spliced RNA sequences with all introns removed and extrons concatenated
def spliceRNA(sequence, introns):
    finalSequence = sequence
    for intron in introns:
        finalSequence = finalSequence.replace(intron, '')
    return finalSequence

def locateSubseqeunce(sequence, subsequence):
    position = []
    counter = 0
    for i in range(len(sequence)):
        if sequence[i] == subsequence[counter]:
            counter += 1
            position.append(i+1)
        if counter == len(subsequence):
            break
    return position

# implement edit distance
# http://rosalind.info/problems/edit/
