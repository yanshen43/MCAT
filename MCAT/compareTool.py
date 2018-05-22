from math import sqrt
import sys

# Written by Jeff Robertson (thejar@vt.edu)
# originally for CS4884 Spring 2016 Class project

# dynamic programing algorithm to compute edit distance between motifs
def compare(seq1, seq2, vals=[-1, 0, 2], pos1=None, pos2=None):
    m, n = len(seq1), len(seq2)
    table = [[0 for j in range(n+1)] for i in range(m+1)]
    def cellMax(i, j):
        up   = table[i][j-1] + vals[0]
        left = table[i-1][j] + vals[0]
        diag = table[i-1][j-1] + vals[1]
        if seq1[i-1] == seq2[j-1]:
            diag += vals[2]
        return max([up, left, diag])
    maxVal = 0
    for i in range(1, m+1):
        for j in range(1, n+1):
            table[i][j] = cellMax(i, j)
            if table[i][j] > maxVal:
                maxVal = table[i][j]
    if pos1 != None and pos2 != None:
        closest = sys.maxint
        for num1 in pos1:
            for num2 in pos2:
                closest = min(abs(num1 - num2), closest)
        if closest < 10:
            maxVal += 1/max(closest, 1)
    return maxVal/sqrt(min(len(seq1), len(seq2)))

def diff(seq1, seq2):
    tot = 0
    for c1, c2 in zip(seq1, seq2):
        tot += 1 if c1 == c2 else 0
    return tot


def test():
    # tests do not reflect normilization for length
    print compare("AAAAA", "AAAAA"), 5
    print compare("AAAAA", "AAAAC"), 4
    print compare("AAAAA", "AAAA"), 4
    print compare("AAAAA", "CCCCC"), 0
    print compare("AAAAA", "ACAGAT"), 3
    print compare("AAAAA", "ACAGATA"), 4
    print compare("AACAAA", "ACAGAT"), 4
    print compare("ACGTA", "ACAGAT"), 4
    print compare("ACGTA", "ACAGAT"), 4
    print compare("CATACGTAAG", "ACAGAT"), 4
    print compare("CATACAGATAAG", "ACAGAT"), 5
    print compare("GGAGACCCCCGCTGCGTGGAGGACAGCGCGCAGCCCCAGCGCGCGGGCCC", "CGTGGA"), 6


if __name__ == "__main__":
    test()
