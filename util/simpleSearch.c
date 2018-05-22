/**
 * Written by Jeff Robertson (thejar@vt.edu)
 */

#include <stdio.h>
#include <malloc.h>

int main(int argc, char** argv) {
    int i;
    if (argc != 5) {
        printf("Invocation:\n%s FASTA_FILE MOTIF verbose(v) #*s\nAll command line arguments are mandatory\n", argv[0]);
        return 1;
    }
    // process commend line args
    char* fname = argv[1];
    char* motif = argv[2];
    char* verbose = argv[3];
    short vis = verbose[0] == 'v';
    char* dontCares = argv[4];
    if (vis) {
        printf("SequenceStart.MatchNum|MatchPos: Instance: *s\n");
    }
    else {
        printf("SequenceStart:MatchPos<tab>\n");
    }
    // find motif len
    int mlen = 0;
    while (motif[mlen] != 0) {
        mlen++;
    }
    // malloc and init circular buffer
    char* curr = malloc((mlen+1)*sizeof(char));
    for (i = 0; i < mlen; i++) {
        curr[i] = 'q';
    }
    curr[mlen] = 0;

    int cp = 0;
    FILE* fp = fopen(fname, "r");
    int stars = dontCares[0]-48;
    char nextChar;
    int lnum = 1; // current line num
    int last_seqName = -1; // line num of most recently seen seq name
    int seqNames[1000]; // line num of seq names (assumes <= 1000 sequences)
    int posCount = 0; // total matches found
    int posSeqCount = 0; // number of sequences containing motif
    int negSeqCount = 0; // number of sequences not containing motif
    int matchPerSeq = 0; // number of matches found in current seq
    int seqPos = 0; // position in current sequence
    while (1 == fscanf(fp, "%c", &nextChar)) {
        if (nextChar == '\n') {
            lnum++;
        }
        if (nextChar == '>') {
            // beginning of new sequence
            last_seqName = lnum;
            if (matchPerSeq == 0) {
                negSeqCount++;
            }
            matchPerSeq = 0;
            seqPos = 0;
        }
        if (nextChar != 'A' && nextChar != 'C' && nextChar != 'G' && nextChar != 'T') {
            continue;
        }

        // add new char to circular buffer
        curr[cp] = nextChar;
        cp++;
        cp %= mlen;
        seqPos++;

        // check current buffer
        int wrong = stars;
        for (i = cp; i < cp + mlen; i++) {
            if (curr[i % mlen] != motif[i - cp]) {
                wrong--;
            }
            if (wrong < 0) {
                break;
            }
        }
        if (wrong >= 0) {
            matchPerSeq++;
            if (!vis) {
                printf("%d:%d\t", last_seqName, seqPos);
            }
            else {
                printf("%d.%d|%d: %s", last_seqName, matchPerSeq, seqPos, &curr[cp]);
                curr[cp] = 0;
                printf("%s: %d\n", curr, stars - wrong);
                curr[cp] = nextChar;
            }
            posCount++;
            if (matchPerSeq == 1) {
                posSeqCount++;
            }
        }
    }
    if (vis) {
        printf("Pos Count: %d, %d seqs\nNeg Count: %d", posCount, posSeqCount, negSeqCount);
    }
    printf("\n");
    free(curr);
    fclose(fp);
    return 0;
}
