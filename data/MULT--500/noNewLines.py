from sys import argv
with open(argv[1], 'r') as fin, open(argv[1]+'.long', 'w') as fout:
    for line in fin:
        if len(line) < 3 or line[0] == '>':
            fout.write(line)
        else:
            fout.write(line.strip())
        if len(line) < 3:
            fout.write('\n')
