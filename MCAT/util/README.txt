simpleSearch performs a naive search through a fasta file for a motif
The command line arguments specify the fasta file to search, the motif to
search for, whether or not it should output the results verbosely, and how
many *'s (don't cares) it should allow.
To indicate that the search should be verbose in its output, put 'v' (without
the quotes) as the third parameter. To not be verbose, put anything else as the
third parameter.

lmerCount counts all l-mers in a fasta file, prints the top 10 most common ones,
and prints a list of how many l-mers there were with each count (i.e. there were
x l-mers that appeared 1 time, y l-mers that appeared 2 times, z lmers that
appeared 3 times, etc...)
It also writes l-mers and counts to counts.json
