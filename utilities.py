def coloured(sequence):
    bcolours = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0m'
    }

    tmpStr = ""

    for nucleotide in sequence:
        if nucleotide in bcolours:
            tmpStr += bcolours[nucleotide] + nucleotide
        else:
            tmpStr += bcolours['reset'] + nucleotide

    return tmpStr + '\033[0;0m'
