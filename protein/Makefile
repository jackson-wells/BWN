all: bwp-index bwp-search

bwp-index: bwp-index.c
	gcc -g -lm -Wall -Wextra -std=c89 -pedantic -Wmissing-prototypes -Wstrict-prototypes \
    -Wold-style-definition -ansi -o bwp-index bwp-index.c

bwp-search: bwp-search.c
	gcc -g -lm -Wall -Wextra -std=c89 -pedantic -Wmissing-prototypes -Wstrict-prototypes \
    -Wold-style-definition -ansi -o bwp-search bwp-search.c

clean:
	rm -rf bwp-search
	rm -rf bwp-index
	rm -rf core.*
	clear
