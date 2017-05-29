#ifndef BWN_SEARCH_H
#define BWN_SEARCH_H

struct index
{
	char *desc;
	int *SA;
	char *transform;
	int length;
    	int **O;
	int **R;
    	int C[6];
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
	int low;
	int high;
};

int *getLength(int seqCount);
int baseMap(char temp);
void printResults(struct output **out, int qsc, int isc,struct index *interval);
void outputToFile(struct output **out, int qsc, int isc,struct index *interval,struct query input);
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
int **calculateD(struct index *interval,int isc,struct query input,int qsc);
struct output ***inExactSearch(struct query input,int QseqCount,struct index *interval,int IseqCount,int **D);
struct output *inexRecur(struct index interval, int *D,char *W,int i,int d, int low, int high);
char revBaseMap(int temp);

#endif
