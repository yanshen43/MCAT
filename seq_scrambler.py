from sys import argv
from random import randint
if len(argv) != 2:
    print("Invocation: %s <FASTA_FILE>"%argv[0])
    quit()
with open(argv[1], 'r') as fin:
    for line in fin:
        if len(line) < 3 or line[0] == '>':
            print line.strip()
        else:
            shuf = 0
            line = list(line.strip())
            while shuf < len(line)/2:
                # pick 2 positions and a length
                p1 = randint(0, len(line) - 1)
                p2 = randint(0, len(line) - 1)
                l = randint(5, 10)
                # if they overlap
                if min([p1, p2]) + l > max([p1, p2]) or max([p1, p2]) + l > len(line):
                    continue
                # swap them
                s1 = line[p1:p1+l]
                s2 = line[p2:p2+l]
                line[p1:p1+l] = s2
                line[p2:p2+l] = s1
                shuf += 1
            print ''.join(line)
