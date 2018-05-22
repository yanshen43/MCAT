#! /usr/bin/python

# Written (mostly) by Jeff Robertson (thejar@vt.edu)
# originally for CS4884 Spring 2016 Class project

import subprocess
from searchTool import makePWM
import pdb

# Here be dragons! (vis code be ugly...)

# use weblogo tool to produce sequence logos from fasta files
def createSeqLogo(RES_DIR, mnum):
    print "Starting visualization"
    cmdline = "./weblogo/seqlogo -F PNG -c -f " +\
        RES_DIR + "motif" + str(mnum) + ".instances" + " > " + RES_DIR +\
        "motif_result_" + str(mnum) + ".png"
    with open(RES_DIR + "weblogo.out", 'w') as fout:
        subprocess.call([cmdline], shell=True, stdout=fout, stderr=fout)

# copy list of motif instances from motif instance files
def parseInstances(RES_DIR, mnum):
    instances = {}
    maxSeqLen = 0
    with open(RES_DIR + "motif"+str(mnum+1)+".instances", "r") as instanceFile:
        for line in instanceFile:
            if len(line) < 2:
                pass
            elif line[0] == '>':
                sequenceName = line[1:line.rfind("\t")]
                seqPos, seqLen = line.split("\t")[-1].split("|")
                maxSeqLen = max(maxSeqLen, int(seqLen))
                if not sequenceName in instances:
                    instances[sequenceName] = []
            else:
                instances[sequenceName].append(line.strip()+"\t"+seqPos+'|'+seqLen)
    return (maxSeqLen, instances)

# visualize and write results to html
def createHTML(RES_DIR, mslist, bgroundProbs):
    with open(RES_DIR+"results.html", "w") as rfile:
        rfile.write('<html><head><meta charset="ascii"></meta><link rel="stylesheet" href="http://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css"></link></head>\n<body>\n<div id="jumpTo" class="container">')
        # write "jump to" header
        for score, motif, pval, zscore in mslist:
            rfile.write('<a href="#'+motif+'">'+motif+'</a><br>')
        rfile.write('</div>')

        # write each motif 
        for mnum, (score, motif, pval, zscore) in enumerate(mslist):
            print motif, score
            ins = []
            maxLen, instances = parseInstances(RES_DIR, mnum)
            for i in instances:
                ins += [x.split("\t")[0] for x in instances[i]]
            PWMstr, logLikelihood = makePWM(len(ins[0]), ins, [bgroundProbs[x] for x in "ACGT"])
            # print header for this motif
            writeMotif(rfile, mnum, score, motif, pval, zscore, logLikelihood, PWMstr, ins)
            # print and visualize each instance
            sequences = instances.keys()
            sequences.sort()
            for sequence in sequences:
                rfile.write("<h4>"+sequence + "</h4>\n")
                for instance in instances[sequence]:
                    writeInstance(instance, maxLen, rfile)
                rfile.write("<br>\n")
            rfile.write('<a href="#jumpTo">Jump to top</a></div>\n')
            rfile.write("</div>")
        rfile.write("<br></body>\n</html>")

# this writes the header (scores, len, PWM, etc.) for a single motif
def writeMotif(rfile, mnum, score, motif, pval, zscore, logLikelihood, PWMstr, ins):
    rfile.write('<div id="'+motif+'">')
    rfile.write('<div class="container"><div style="display:inline-block;">')
    rfile.write('<h1>Motif: ' + motif + "</h1>\n")
    rfile.write("Length: " + str(len(motif))+"<br>\n")
    rfile.write("Comparison Score: {:.3f}<br>".format(score))
    rfile.write("Log Likelihood: {:.3f}<br>".format(logLikelihood))
    rfile.write("P-value: {:.3e}<br>".format(pval))
    rfile.write("Z-Score: {:.3e}<br>".format(zscore))
    rfile.write('</div><div style="display:inline-block;vertical-align:top;">')
    rfile.write('<img src="motif_result_'+str(mnum+1)+'.png"></div></div>\n')
    rfile.write('<div class="container"><h2>Position Weight Matrix</h2>\n')
    rfile.write('<pre style="margin-top: 0px; margin-bottom: 0px;">')
    rfile.write(PWMstr)
    rfile.write("</pre>\n</div>\n")
    rfile.write('<div class="container"><h2>'+str(len(ins))+" instances found</h2>\n")

# this writes a found instance (instance, vis of pos in seq, pos) of a motif
def writeInstance(instance, maxLen, rfile):
    insPos, seqLen = instance.split("\t")[-1].split('|')
    seqLen = int(seqLen)
    insPos = seqLen - int(insPos)
    # vis is a string representing the motif's position in the sequence
    vis = (seqLen)*110/maxLen * "_"
    vis = vis[:insPos*110/maxLen] + '|' + vis[insPos*110/maxLen+1:]
    vis = (maxLen - seqLen)*110/maxLen * "&nbsp" + vis
    rfile.write(instance[:instance.rfind("\t")] + "\n<br>\n")
    rfile.write('<pre style="margin-top: 0px; margin-bottom: 0px;">')
    rfile.write(vis+" + "+str(insPos)+" - "+str(seqLen - insPos)+" | "+str(seqLen))
    rfile.write("</pre>\n")

def visualizeOutput(RES_DIR, mlist, probabilities):
    for mnum, motif in enumerate(mlist):
        createSeqLogo(RES_DIR, mnum+1)
    createHTML(RES_DIR, mlist, probabilities)

def main():
    createHTML("results/", ["ACGT", "ACTG", "AGCT", "AGTC"])

if __name__ == "__main__":
    main()
