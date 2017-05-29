#ifndef BWN_INDEX_H
#define BWN_INDEX_H

struct suffix
{
    	int pos;
    	char *string;
};

struct input_data
{
    	char **name;
    	int *length;
	char **sequence;
	char **reverse;
};

struct index
{
    	int **O;
	int **R;
    	int C[6];
};

struct suffix **buildReverseSuffixArray(struct input_data input,int seqCount);
void populateReverseSuffixArray(struct suffix **m,struct input_data input,int seqCount);
char *reverse(char *str,int length);
int baseMap(char temp);
int mostChars(char *fileName);
void charToEnd(char *input);
void initializeSuffixArray(struct suffix **m, int *seqLength,int seqCount);
void freeSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void deleteSuffixArray(struct suffix **m, int seqCount,int *seqLength);
void sortSuffixArray(struct suffix **m,int seqCount,int *seqLength);
char **bwt(struct suffix **m,int seqCount, int *seqLength);
struct input_data manageInputs(int argc, char* argv[],int *sCount);
int seqCount(char *fileName);
int *seqLength(char *fileName,int seqCount,int cCount);
struct suffix **newSuffixArray(int seqCount, int *seqLength);
struct suffix **buildSuffixArray(struct input_data input,int seqCount);
void printSuffixArray(struct suffix **m,int seqCount,int *seqLength);
void printBwt(struct suffix **m,char **transform,int seqCount);
int *getLengths(struct input_data input,int seqCount);
struct input_data initializeInputStruct(int sCount,int *sLength);
void read_fasta(char *fValue,struct input_data *input);
void handleS(struct input_data *input, char *sequence);
void handleF(struct input_data *input, char *fileName);
void populateSuffixArray(struct suffix **m,struct input_data input,int seqCount);
void deleteInputStruct(struct input_data input,int seqCount);
struct index *calculateInterval(char **transform, int *seqLength,int seqCount,char **revTransform);
int *calculateO(char *sequence,int seqLength,int letterValue);
int calculateC(char *sequence, int seqLength,int letterValue);
void intervalToFile(struct index *FMidx, int seqCount,struct suffix **m,char **transform,struct input_data input);
void printInt(struct index *FMidx, int seqCount,int *seqLength);
void mergeSort(int low, int high,struct suffix **m,int seqNumber);
void Merge(int low,int mid, int high, struct suffix **m,int seqNumber);

#endif
