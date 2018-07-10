#ifndef DATASET_H_
#define DATASET_H_

typedef int (*fptrAt)(void*, int, int);
typedef void (*fptrFree)(void*);

typedef struct _functions {
  // Functions
  fptrAt at;
  fptrFree free;
  
} vFunctions;

typedef struct _dataset {
  vFunctions functions;
  
  // Base variables
  int N;
  int D; // number of rows and dimension

} Dataset;

int DatasetAt(void* dataset, int n, int d) { return 0; }
 
void DatasetFree(void* dataset) {
  free(dataset);
}
 
Dataset* getDatasetInstance() {
  Dataset* dataset = (Dataset*) malloc(sizeof(Dataset));
  dataset->functions.at = DatasetAt;
  dataset->functions.free = DatasetFree;
  dataset->N = 0;
  dataset->D = 0;
  return dataset;
}

#endif /* DATASET_H_ */
