def filterTool(finame, foname, lout):
    with open(finame, 'r') as fin, open(foname, 'w') as fout:
        linelen = 0
        buf = ''
        for line in fin:
            if len(line) < 2:
                continue
            sline = line.strip()
            if sline[0] == '>':
                fout.write((buf + '\n\n' if linelen != 0 else '') + sline + '\n')
                buf = ''
            else:
                if len(sline) > linelen:
                    linelen = len(sline)
                for c in sline:
                    if c in {'A', 'C', 'G', 'T'}:
                        buf += c
                    elif c in {'a', 'c', 'g', 't'}:
                        buf += 'N'
                while len(buf) > linelen:
                    fout.write(buf[:linelen] + '\n')
                    buf = buf[linelen:]
        fout.write(buf + '\n')

if __name__ == '__main__':
    from sys import stdout
    filterTool('data/MULT--500/positive_500.marked.fasta', 'data/MULT--500/positive_500.filtered.fasta', stdout)
