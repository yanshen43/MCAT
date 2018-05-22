all: simpleSearch lmerCount
	
simpleSearch: simpleSearch.c
	gcc -O3 -o simpleSearch simpleSearch.c

lmerCount: lmerCount.c
	gcc -O3 -o lmerCount lmerCount.c

clean:
	$(RM) simpleSearch
	$(RM) lmerCount
	$(RM) counts.json

