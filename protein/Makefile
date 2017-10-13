all: bwp-index bwp-search

bwp-index: bwp-index.c
	gcc -o bwp-index bwp-index.c

bwp-search: bwp-search.c
	gcc -o bwp-search bwp-search.c

clean:
	rm -rf bwp-search
	rm -rf bwp-index
