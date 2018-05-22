/**
 * Written by Jeff Robertson (thejar@vt.edu)
 */

#include <stdio.h>
#include <malloc.h>

int LEN = 9;
int size = 5*5*5*5*5*5*5*5*5;
char last;
char* fFile = "pos.fasta";
char* outfile = "counts.json";
int strict = 0;
int stars = 4;
int results = 10;
int hFilterCutoff = 1;
int sFilterCutoff = 1;

// position sequence pair represents a possible motif instance
typedef struct {
    int pos;
    int seq;
} pair;

/**
 * Map a character to its associated value
 */
int mapc(char c) {
    switch (c) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        case '*':
            return 4;
        default:
            return -1;
    }
}

/**
 * Map a motif to a numerical value
 */
int map(char * motif, int mptr) {
    int val = 0;
    int i;
    for (i = 0; i < LEN; i++) {
        val *= 5;
        val += mapc(motif[(i + mptr) % LEN]);
    }
    return val;
}

/**
 * Map a numerical value to a motif
 */
void revMap(int val, char* motif) {
    int i;
    char chars[] = {'A', 'C', 'G', 'T', '*'};
    for (i = LEN-1; i > -1; i--) {
        motif[i] = chars[val % 5];
        val /= 5;
    }
}

/**
 * Consume a line from fptr
 */
void line(FILE* fptr) {
    char c = 0;
    while (c != '\n') {
        fscanf(fptr, "%c", &c);
    }
}

/**
 * Increment the table entry associated with the given motif
 */
void inc(short* mapTable, char* motif, int mptr) {
    mapTable[map(motif, mptr)]++;
}

/**
 * Recursively add the desired number of stars to a motif
 * and increment the table entries associated with each version of the motif
 */
void addStar(short* mapTable, int nstars, int lastStar, char* motif, char* mcpy, int mptr) {
    if (nstars < 1) {
        return;
    }
    int i;
    for (i = lastStar + 1; i < LEN; i++) {
        mcpy[(i + mptr) % LEN] = '*';
        inc(mapTable, mcpy, mptr);
        addStar(mapTable, nstars - 1, i, motif, mcpy, mptr);
        mcpy[(i + mptr) % LEN] = motif[(i + mptr) % LEN];
    }
}

/**
 * Filter for homogenius regions, then call increment table for l-mer
 */
void insert(short* mapTable, char* motif, int mptr, pair p) {
    int i, j, k, l;
    // skip homogenius regions
    static pair lastHom[4];
    char c1 = motif[0], c2 = motif[LEN-1];
    int diff = 0;
    for (i = 1; i < LEN; i++) {
        if (motif[i] == '*') {
            continue;
        }
        if (motif[i] != c1 && motif[i] != c2) {
            diff++;
            if (diff > 1) {
                break;
            }
        }
    }
    if (i == LEN) {
        int indx = mapc(c1);
        pair old = lastHom[indx];
        lastHom[indx] = p;
        if (old.seq == p.seq && old.pos > p.pos - LEN) {
            return;
        }
    }
    // inc with base motif
    inc(mapTable, motif, mptr);
    char mcpy[LEN];
    for (i = 0; i < LEN; i++) {
        mcpy[i] = motif[i];
    }
    addStar(mapTable, stars, -1, motif, mcpy, mptr);
    for (i = 0; i < LEN; i++) {
        if (motif[i] != mcpy[i]) {
            printf("Error!\n%s\n", motif);
        }
    }
}

/**
 * Parse command line arguments and set respective options
 */
int parseArgs(int argc, char** argv) {
    int i, j, ff = 1, help = 0;
    if (argc > 1) {
        for (i = 1; i < argc; i++) {
            char* arg = argv[i];
            printf("%d(%s)\n", i, arg);
            if (arg[0] != '-' && ff) {
                fFile = argv[i];
                printf("setting fasta to '%s'\n", fFile);
                ff = 0;
            }
            else if (strncmp(arg, "-h", 2) == 0) {
                help = 1;
            }
            else if (strncmp(arg, "--no-filter", 11) == 0) {
                strict = 1;
            }
            else if (strncmp(arg, "-l", 2) == 0) {
                i++;
                LEN = atoi(argv[i]);
                size = 1;
                for (j = 0; j < LEN; j++) {
                    size *= 5;
                }
                printf("setting len to %d and size to %d\n", LEN, size);
            }
            else if (strncmp(arg, "-s", 2) == 0) {
                i++;
                stars = atoi(argv[i]);
                printf("setting stars to %d\n", stars);
            }
            else if (strncmp(arg, "-r", 2) == 0) {
                i++;
                results = atoi(argv[i]);
                printf("will print %d results\n", results);
            }
            else if (strncmp(arg, "-fh", 3) == 0) {
                i++;
                hFilterCutoff = atoi(argv[i]);
                printf("set homogeneous filter cutoff to %d\n", hFilterCutoff);
            }
            else if (strncmp(arg, "-fs", 3) == 0) {
                i++;
                sFilterCutoff = atoi(argv[i]);
                printf("set similarity filter cutoff to %d\n", hFilterCutoff);
            }
            else if (strncmp(arg, "-o", 2) == 0) {
                i++;
                outfile = argv[i];
                printf("set outfile to %s\n", outfile);
            }
            else {
                printf("unrecognized arg '%s'\n", arg);
            }
        }
    }
    if (argc > 13 || help || ff) {
        printf("INVOCAION:\n./hmf [-h] [-l len] [-s stars] [-fh hFilterCutoff]");
        printf(" [-fs sFilterCutoff] [-r resultsPrinted] [-o outputFile]");
        printf(" [--no-filter] fastaFile\n");
        if (help) {
            printf("\n-h\tdisplay this help messge and exit\n\n");
            printf("-l\tset the desired l-mer length (default=9)\n\n");
            printf("-s\tset the number of don't care characters (default=2)\n\n");
            printf("-fh\tset the number of characters that must differ before a\n");
            printf("\tl-mer is not considered low-complexity (default=1)\n\n");
            printf("-fs\tset the number of characters that must differ before a\n");
            printf("\tl-mer is not considered similar to a previously found l-mer (default=1)\n\n");
            printf("-r\tset the number of results to print (default=10)\n\n");
            printf("-o\tset the output file\n\n");
            printf("--no-filter\tdisable the low-complexity and similarity filters\n");
            printf("\tif the filter is on, 'h' is printed instead of a low-complexity l-mer\n");
            printf("\tand si instead of an l-mer that is similar to an already displayed l-mer\n");
            printf("\twhere i is the index of the displayed l-mer that it was similar to\n\n");
            printf("fastafile\tthe name of the fasta file this tool should be run on\n\n");
        }
        return 1;
    }
    return 0;
}

int main(int argc, char** argv) {
    int i, j, k, l;
    // parse args
    if (parseArgs(argc, argv)) {
        return 0;
    }
    
    FILE* ifile = fopen(fFile, "r");
    short* mapTable = malloc(sizeof(short) * size);
    for (i = 0; i < size; i++) {
        mapTable[i] = 0;
    }
    char nMotif[LEN + 1];
    nMotif[LEN] = 0;
    int mptr = 0;
    char c;
    pair p;
    p.pos = 0;
    p.seq = 0;
    int seqPos;
    // scan file
    while (fscanf(ifile, "%c", &c) == 1) {
        // note sequence start lines
        if (c == '>') {
            line(ifile);
            for (i = 0; i < LEN; i++) {
                fscanf(ifile, "%c", nMotif + i);
            }
            mptr = 0;
            p.seq++;
            seqPos = 0;
        }
        else if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            continue;
        }
        else {
            nMotif[mptr++] = c;
            mptr %= LEN;
            last = c;
            p.pos++;
            seqPos++;
        }
        if (seqPos >= LEN - 1) {
            // add instance to table
            insert(mapTable, nMotif, mptr, p);
        }
    }
    fclose(ifile);
    FILE* jout = fopen(outfile, "w");
    fprintf(jout, "{");
    char motif[LEN+1];
    motif[LEN] = 0;
    int counts[50000];
    for (i = 0; i < 50000; i++) {
        counts[i] = 0;
    }
    int first = 1;
    for (i = 0; i < size; i++) {
        counts[mapTable[i]]++;
        if (mapTable[i] > 0) {
            revMap(i, motif);
            if (!first) {
                fprintf(jout, ",");
            }
            fprintf(jout, "\n\"%s\":%d", motif, mapTable[i]);
            first = 0;
        }
    }
    fprintf(jout, "}\n");
    fclose(jout);
    // print top scoring motifs
    printf("max motifs are:\n");
    int prevPrinted = 0;
    int nextToPrint = 0;
    int err = 0;
    char* prevP[10];
    for (i = 0; i < 10; i++) {
       prevP[i] = malloc(LEN * sizeof(char));
    }
    
    for (i = 0; i < results; i++) {
        nextToPrint = 0;
        // find best that hasn't been printed
        for (j = 0; j < size; j++) {
            if ((i == 0 || mapTable[j] < mapTable[prevPrinted]
                       || (mapTable[j] == mapTable[prevPrinted] && j > prevPrinted))
                    && mapTable[j] > mapTable[nextToPrint]) {
                nextToPrint = j;
            }
        }
        prevPrinted = nextToPrint;
        revMap(nextToPrint, motif);
        if (!strict) {
            // filter out low complexity motifs
            int diff = 0;
            int fNonS = 0;
            while (motif[fNonS] == '*') {
                fNonS++;
            }
            for (j = fNonS+1; j < LEN; j++) {
                if (motif[j] != motif[fNonS] && motif[j] != '*') {
                    diff++;
                    if (diff > hFilterCutoff) {
                        break;
                    }
                }
            }
            if (j == LEN) {
                i--;
                fprintf(stderr, "h");
                err = 1;
                continue;
            }
            // filter out motifs similar to ones previously printed
            int same = 0;
            int sameNum;
            for (j = 0; j < i; j++) {
                for (k = 0; k < LEN; k++) {
                    same = 0;
                    for (l = 0; l < LEN; l++) {
                        if (motif[(l + k) % LEN] == prevP[j][l]
                               && prevP[j][l] != '*') {
                            same++;
                            if (same > LEN - stars - sFilterCutoff) {
                                sameNum = j;
                                break;
                            }
                        }
                    }
                    if (same > 4) {
                        break;
                    }
                }
                if (same > 4) {
                    break;
                }
            }
            if (same > 4) {
                //fprintf(stderr, "Rejecting %s because it is too similar\n", motif);
                fprintf(stderr, "s%d", sameNum);
                err = 1;
                i--;
                continue;
            }
        }
        if (err) {
            printf("\n");
            err = 0;
        }
        printf("%s(%d) with %d hits\n", motif, nextToPrint, mapTable[nextToPrint]);
        for (j = 0; j < LEN; j++) {
            prevP[i][j] = motif[j];
        }
    }
    for (i = 0; i < 10; i++) {
        free(prevP[i]);
    }

    // print counts of each
    int counti = 1, rowLen, lastCount;
    while (counti < 50000) {
        rowLen = 0;
        lastCount = counti;
        while (rowLen < 20 && counti < 50000) {
            if (counts[counti] != 0) {
                printf("%7.2d|", counti);
                rowLen++;
            }
            counti++;
        }
        printf("\n");
        for (i = lastCount; i < counti; i++) {
            if (counts[i] != 0) {
                printf("%7.2d|", counts[i]);
            }
        }
        printf("\n\n");
    }
    free(mapTable);
    return 0;
}
