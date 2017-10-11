all: bwn-index bwn-search

bwn-index: bwn-index.c
	gcc -o bwn-index bwn-index.c

bwn-search: bwn-search.c
	gcc -o bwn-search bwn-search.c

clean:
	rm -rf bwn-search
	rm -rf bwn-index
