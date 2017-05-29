#ifndef BWN_SEARCH_H
#define BWN_SEARCH_H

struct FMidx
{
	char *desc;
	int *SA;
	char *transform;
	int length;
    	int **O;
	int **R;
    	int C[6];
};

struct input
{
	char *fileName;
	char **name;
	int *length;
	char **sequence;
};

struct output
{
	char *sequence;
	int low;
	int high;
};

int *getLength(int seqCount);
int baseMap(char temp);
void printResults(struct output **out, int qsc, int isc,struct FMidx *index);
void outputToFile(struct output **out, int qsc, int isc,struct FMidx *index,struct input query);
struct output **search(struct input query,int qsc,struct FMidx *index,int isc);
void read_fasta(char *fileName, struct input *query);
struct FMidx *getIndex(int *seqCount);
struct input manageInputs(char *argv[], int argc,int *seqCount);
void handleF(struct input *query,char *fileName);
void handleS(struct input *query,char *sequence);
struct input initializeInputStruct( int seqCount, int *seqLength);
int *getSeqLength(char *fileName,int seqCount,int charCount);
int getSeqCount(char *fileName);
char *removePrefix(char *query);
int **calculateD(struct FMidx *index,int isc,struct input query,int qsc);
struct output ***inExactSearch(struct input query,int QseqCount,struct FMidx *index,int IseqCount,int **D);
struct output *inexRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high);
char revBaseMap(int temp);

#endif
