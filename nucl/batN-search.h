#ifndef BATN_SEARCH_H
#define BATN_SEARCH_H

struct index
{
	char *desc;
	int *SA;
	char *transform;
	int length;
    	int **O;
    	int C[5];
};

struct query
{
	char *fileName;
	char **name;
	int *length;
	char **sequence;
};

struct output
{
	char *sequence;
	int start;
	int end;
};

int *getLength(int seqCount);
int baseMap(char temp);
void outputToFile(struct output **out, int qsc, int isc);
struct output **search(struct query input,int qsc,struct index *interval,int isc);
void read_fasta(char *fileName, struct query *input);
struct index *getIndex(int *seqCount);
struct query manageInputs(char *argv[], int argc,int *sCount);
void handleF(struct query *input,char *fileName);
void handleS(struct query *input,char *sequence);
struct query initializeInputStruct( int sCount, int *sLength);
int *seqLength(char *fileName,int seqCount,int cCount);
int seqCount(char *fileName);
char *removePrefix(char *input);

#endif
