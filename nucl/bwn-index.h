#ifndef BWN_INDEX_H
#define BWN_INDEX_H

struct suffix
{
    	int pos;
    	char *string;
};

struct input
{
    	char **name;
    	int *length;
	char **sequence;
	char **reverse;
};

struct FMidx
{
    	int **O;
	int **R;
    	int C[6];
};

struct suffix **buildReverseSuffixArray(struct input query,int seqCount);
void populateReverseSuffixArray(struct suffix **m,struct input query,int seqCount);
char *reverse(char *str,int length);
int baseMap(char temp);
int mostChars(char *fileName);
void charToEnd(char *input);
void initializeSuffixArray(struct suffix **m, int *seqLength,int seqCount);
void freeSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void deleteSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void sortSuffixArray(struct suffix **m,int seqCount,int *seqLength);
char **bwt(struct suffix **m,int seqCount, int *seqLength);
struct input manageInputs(int argc, char* argv[],int *sCount);
int seqCount(char *fileName);
int *seqLength(char *fileName,int seqCount,int cCount);
struct suffix **newSuffixArray(int seqCount, int *seqLength);
struct suffix **buildSuffixArray(struct input query,int seqCount);
void printSuffixArray(struct suffix **m,int seqCount,int *seqLength);
void printBwt(struct suffix **m,char **transform,int seqCount);
int *getLengths(struct input query,int seqCount);
struct input initializeInputStruct(int sCount,int *sLength);
void read_fasta(char *fValue,struct input *query);
void handleS(struct input *query, char *sequence);
void handleF(struct input *query, char *fileName);
void populateSuffixArray(struct suffix **m,struct input query,int seqCount);
void deleteInputStruct(struct input query,int seqCount);
struct FMidx *calculateInterval(char **transform, int *seqLength,int seqCount,char **revTransform);
int *calculateO(char *sequence,int seqLength,int letterValue);
int calculateC(char *sequence, int seqLength,int letterValue);
void intervalToFile(struct FMidx *index, int seqCount,struct suffix **m,char **transform,struct input query);
void printInt(struct FMidx *index, int seqCount,int *seqLength);
void mergeSort(int low, int high,struct suffix **m,int seqNumber);
void Merge(int low,int mid, int high, struct suffix **m,int seqNumber);

#endif
