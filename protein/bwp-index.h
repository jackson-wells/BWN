#ifndef BWP_INDEX_H
#define BWP_INDEX_H

struct suffix
{
    	int pos;
    	char *string;
};

struct transform
{
	int *positions;
	char *string;
};

struct input
{
    	char **name;
    	int *length;
	int maxLength;
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

struct transform *calculateTransform(struct input query, int seqCount);
struct suffix *buildSuffixArray(char *sequence,int length);
int *positions(struct suffix *m, int seqLength);
struct suffix *populateSuffixArray(char *sequence, int length);
char *reverse(char *str,int length);
int baseMap(char temp);
int mostChars(char *fileName);
void charToEnd(char *input,int len);
void initializeSuffixArray(struct suffix *m, int seqLength);
void freeSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void deleteSuffixArray(struct suffix **m, int seqCount,int *seqLength);
char *bwt(struct suffix *m,int seqLength);
struct input manageInputs(int argc, char* argv[],int *seqCount);
int getSeqCount(char *fileName);
int *getLength(char *fileName,int seqCount);
struct suffix *newSuffixArray(int seqLength);
void printSuffixArray(struct suffix **m,int seqCount,int *seqLength);
void printBwt(char **transform,int seqCount);
struct input initializeInputStruct(int seqCount,int *seqLength);
void readFasta(char *fValue,struct input *query);
void handleS(struct input *query, char *sequence);
void handleF(struct input *query, char *fileName);
void deleteInputStruct(struct input query,int seqCount);
struct FMidx *calculateInterval(struct transform *transform, int seqCount,struct transform *revTransform,int *seqLength);
int *calculateO(char *sequence, int letterValue,int length);
int calculateC(char *sequence, int letterValue,int length);
void intervalToFile(struct FMidx *index, int seqCount,struct transform *transform,struct input query,struct transform *revTransform);
void mergeSort(int low, int high,struct suffix *m, int length);
void Merge(int low,int mid, int high, struct suffix *m,int seqLength);
int fileExists(char *temp);
int extensionExists(char *temp);
int getLineCount(char *fileName);
void formatFasta(char *fileName);
#endif
