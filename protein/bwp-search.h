#ifndef BWP_SEARCH_H
#define BWP_SEARCH_H

struct FMidx
{
	char *desc;
	int *SA;
	char *transform;
	char *reverse;
	int length;
    	int **O;
	int **R;
    	int C[20];
	char *sequence;
};

struct input
{
	char *fileName;
	char **name;
	int *length;
	char **sequence;
};

struct matches
{
	int low;
	int high;
	struct matches *next;
	int score;
	char *tb;
	int traceLength;
};

struct output 
{
	int low;
	int high;
	char *sequence;
};

struct results
{
	struct matches *match;
	int matches;
};

extern char INTERVAL_FILE[];

void getHighScore(struct matches *head, int k, int l, int score,char *traceback,int traceLength);
int getScore(int l1, int l2,int score,int p,int c);
struct matches *getUnion(struct matches *head1, struct matches *head2);
void push(struct matches** head_ref, int k,int l,int score,char *traceBack,int traceLength);
int isPresent(struct matches *head, int k, int l);
int *getLength(int seqCount);
int baseMap(char temp);
void outputToFile(struct results **out, int qsc, int isc,struct FMidx *index,struct input query);
struct output **search(struct input query,int qsc,struct FMidx *index,int isc);
void read_fasta(char *fileName, struct input *query);
struct FMidx *getIndex(void);
struct input manageInputs(char *argv[], int argc,int *seqCount);
void handleF(struct input *query,char *fileName);
void handleS(struct input *query,char *sequence);
struct input initializeInputStruct( int seqCount, int *seqLength);
int *getSeqLength(char *fileName,int seqCount);
int getCount(void);
int getSeqCount(char *fileName);
char *removePrefix(char *query);
int ***calculateD(struct FMidx *index,int isc,struct input query,int qsc);
struct results **inexactSearch(struct input query,int QseqCount,struct FMidx *index,int IseqCount,int ***D);
struct matches *inexRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high,int score,int pState,char *traceBack,int tbIdx);
char revBaseMap(int temp);
void printInResults(struct results **out,int qsc,int isc,struct FMidx *index, struct input query);
struct matches *pointToTail(struct matches *match);
struct matches *append(struct matches *match, struct matches *results);
int fileExists(char *temp);
void printKLs(struct results **out, int qsc, int isc);
void frontbacksplit(struct matches* source, struct matches** frontRef, struct matches** backRef);
struct matches* sortedmergeMatch(struct matches* a, struct matches* b);
void mergeSortMatches(struct matches** headRef);
struct matches *sortMatches(struct matches *match);
char *reverse(char *str);
char *getSequenceAlignment(int i,int j,struct FMidx *index, struct input query,struct matches *match);
#endif
