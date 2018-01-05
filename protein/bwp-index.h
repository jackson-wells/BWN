#ifndef BWP_INDEX_H
#define BWP_INDEX_H

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
    	int C[20];
};

extern char OUTPUT_FILE[]; 

struct suffix *buildSuffixArray(char *sequence);
struct suffix *populateSuffixArray(char *sequence);
char *reverse(char *str);
int baseMap(char temp);
int mostChars(char *fileName);
void charToEnd(char *input);
void initializeSuffixArray(struct suffix *m, int seqLength);
void freeSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void deleteSuffixArray(struct suffix **m, int seqCount,int *seqLength);
char **bwt(struct suffix **m,int seqCount, int *seqLength);
struct input manageInputs(int argc, char* argv[],int *seqCount);
int getSeqCount(char *fileName);
int *getLength(char *fileName,int seqCount);
struct suffix *newSuffixArray(int seqLength);
void printSuffixArray(struct suffix **m,int seqCount,int *seqLength);
void printBwt(char **transform,int seqCount);
struct input initializeInputStruct(int seqCount,int *seqLength);
void read_fasta(char *fValue,struct input *query);
void handleS(struct input *query, char *sequence);
void handleF(struct input *query, char *fileName);
void deleteInputStruct(struct input query,int seqCount);
struct FMidx *calculateInterval(char **transform, int seqCount,char **revTransform);
int *calculateO(char *sequence, int letterValue);
int calculateC(char *sequence, int letterValue);
void intervalToFile(struct FMidx *index, int seqCount,struct suffix **m,char **transform,struct input query,char **revTransform);
void mergeSort(int low, int high,struct suffix *m);
void Merge(int low,int mid, int high, struct suffix *m);
int fileExists(char *temp);
int extensionExists(char *temp);


#endif
