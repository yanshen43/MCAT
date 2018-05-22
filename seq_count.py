import sys

with open(sys.argv[1], 'r') as f:
    seqs = 0
    for line in f:
        if line[0] == '>':
            seqs += 1
print sys.argv[1] + " contains " + str(seqs) + " sequences"
