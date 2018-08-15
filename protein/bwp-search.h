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
	int keep;
};

struct output /* artifact of exact search */
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
void push(struct matches** head_ref, int k,int l,int score,char *traceBack,int traceLength, int keep);
int isPresent(struct matches *head, int k, int l);
int *getLength(int seqCount);
int baseMap(char temp);
void outputToFile(struct results **out, int qsc, int isc,struct FMidx *index,struct input query);
struct output **exactSearch(struct input query,int qsc,struct FMidx *index,int isc);
void readFasta(char *fileName, struct input *query);
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
struct results **conservedSearch(struct input query,int QseqCount,struct FMidx *index,int IseqCount,int ***D);
struct matches *conservedRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high,int score,int pState,char *traceBack,int tbIdx);
struct results **distanceSearch(struct input query,int QseqCount,struct FMidx *index,int IseqCount,int ***D);
struct matches *distanceRecur(struct FMidx index, int *D,char *W,int i,int d, int low, int high,int score,int pState,char *traceBack,int tbIdx);
char revBaseMap(int temp);
void printInResults(struct results **out,int qsc,int isc,struct FMidx *index, struct input query);
struct matches *pointToTail(struct matches *match);
struct matches *append(struct matches *match, struct matches *results);
int fileExists(char *temp);
void printKLs(struct results **out, int qsc, int isc);
void frontbacksplit(struct matches* source, struct matches** frontRef, struct matches** backRef);
struct matches* sortedmergeMatch(struct matches* a, struct matches* b);
void mergeSortMatches(struct matches** headRef);
struct matches *sortMatches(struct matches *match, struct FMidx input);
char *reverse(char *str);
char *getSequenceAlignment(int i,int j,struct FMidx *index, struct input query,struct matches *match);
int ***calculateS(struct FMidx *index,int isc, struct input query,int qsc, int *St);
struct matches *scoredRecur(struct FMidx index, int *Sp,char *W,int i, int low, int high,int score,int pState,char *traceBack,int tbIdx,int St);
struct results **scoredSearch(struct input query,int qsc,struct FMidx *index,int isc,int ***S,int *St);
void filterMatches(struct matches **head,struct FMidx input);
int roundFloat(float num);
void readSubMat(int selection);
int min(int a, int b);

#endif
