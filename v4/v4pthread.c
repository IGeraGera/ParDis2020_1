#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "mmio.h"
/* Function Delclaration */
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
);
void checkArgs(int argc, char *argv[]);
/* Global variables declaration */
uint32_t * csc_r;
uint32_t * csc_c;
uint32_t * nodes;
/* pthread parameter function */
typedef struct param{
	uint32_t indexstart;
	uint32_t indexend;
	pthread_mutex_t *mut;
}param;
/* pthread work functrion */
void *work(void *arg){
	struct param *p =  (struct param *)arg;
	for (uint32_t j=p->indexstart; j<p->indexend;j++){
		uint32_t jstart=csc_c[j];
		uint32_t jend=csc_c[j+1];
		uint32_t totval=0;
		for(uint32_t ptr=jstart ;ptr<jend; ptr++){
			uint32_t i = csc_r[ptr];
			uint32_t istart = csc_c[i];
			uint32_t iend = csc_c[i+1];
			for (uint32_t ptrj = jstart; ptrj<jend;ptrj++ ){
				for (uint32_t ptri = istart;ptri<iend;ptri++){
					if (csc_r[ptrj]<csc_r[ptri]){ break;}
					if (csc_r[ptrj]==csc_r[ptri]){
						totval++;	
						istart = ptri+1;
					}
				}
			}
		}
		pthread_mutex_lock(p->mut);
		nodes[j]=totval/2;
		pthread_mutex_unlock(p->mut);
	}
	pthread_exit(NULL);
}
int main(int argc, char *argv[]){
	/* Timer structs */
	struct timespec ts_start;
	struct timespec ts_end;
	/* Variable declaration */
	FILE *f;
	uint32_t M,N,nz,ret;
	MM_typecode matcode;
	
	/* Check if arguments are correct */ 
	checkArgs(argc,argv);
	/* If argument is 1 then open file  */		
	if ((f = fopen(argv[1],"r")) == NULL){
		printf("Problem opening the file");
		exit(1);
	}
	
	/* Read Banner*/
	if(mm_read_banner(f,&matcode) != 0){
		printf("Matrix banner read not succesful\n");
		exit(1);
	}
	/* If banner is wrong format (array,dense) exit */
       if ((mm_is_matrix(matcode)==0)||(mm_is_sparse(matcode)==0)){
       		printf("Check the format from the banner\n");
		printf("Must be matrix, sparse \n ");
		printf("Banner\n %s",mm_typecode_to_str(matcode));
       }	
	
       /* Read the matrix sizes for coordinate format*/
	if ((ret=mm_read_mtx_crd_size(f,&M,&N,&nz))!=0){
		printf("Cannot read coordinate");
		exit(0);
	}
	/* Setting number of threads from the argument */ 	
	uint32_t num_threads = atoi(argv[2]);
	/* Declaration of COO format arrays  */ 
	uint32_t *coo_r,*coo_c,*val;

    	/* Allocate memory for matrices non-zero elements *2 */
   	coo_r = (uint32_t *) malloc(2*nz * sizeof(uint32_t ));
    	coo_c = (uint32_t *) malloc(2*nz * sizeof(uint32_t ));
    	val  = (uint32_t *) malloc(2* nz * sizeof(uint32_t ));
	/* Read the  .mtx file and pass values to coo_r, coo_c, val
	 * NOTE 1: The iteration is increasing by 2
	 * NOTE 2: Only symmetrical .mtx files must be given.
	 * NOTE 3: The values must be only the lower triangle */
	for (uint32_t i=0; i<2*nz; i+=2){
		/* Scan the values from file */
		if(fscanf(f,"%d %d \n",&coo_r[i],&coo_c[i])==2){
		/* Subtract 1 to make zerobased */
		coo_r[i]--;
		coo_c[i]--;
		/* If same values on column and row ignore them. It means that
		 * there is a value on main diagonal of the matrix */
		if (coo_r[i]==coo_c[i]) continue;
		/* Put to the next position of the array the symmetrical point (i,j)->(j,i) to generate the whole matrix*/
		coo_c[i+1]=coo_r[i];
		coo_r[i+1]=coo_c[i];
		val[i]=1;
		val[i+1]=1;
		}
	}
	/* Close the file  */
	fclose(f);
	/* Declaration of CSC format arrays  */ 
	uint32_t const row = M;
	uint32_t const nnz = nz*2;
	uint32_t isOneBased =0;
    	 /*Allocate memory for csc  matrices */
  	csc_r = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
  	csc_c = (uint32_t *)malloc((row+1) * sizeof(uint32_t));
	/* Make the conversion from CSC to COO */	
	coo2csc(csc_r,csc_c,coo_r,coo_c,nnz,row,isOneBased);
	/* Allocate and init to zero nodes array */
	nodes = (uint32_t *)malloc(row * sizeof(uint32_t));
	for (uint32_t x = 0; x<row;x++) nodes[x]=0;
	/* Allocate an array of threads */
	pthread_t *threads = malloc(num_threads*sizeof(pthread_t));
	/* Alocate and initialize a mutex */
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex,NULL);
	/* Starting timing */
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	/* nodes = ((A .* (A*A)).*e)/2 */
	uint32_t interval  = row/num_threads;
	uint32_t count=0;
	uint32_t extra=row%num_threads;
	for (uint32_t j=0; j<num_threads;j++){
		struct param *arg = (struct param *)malloc(sizeof(struct param));
		arg->indexstart= (j* interval)+count;
		if (extra!=0){
			count++;
			extra--;
		}
		arg->indexend = (j*interval)+count+interval;
		//printf("%d %d\n",arg->indexstart,arg->indexend);
		arg->mut = &mutex;
		pthread_create(&threads[j],NULL,work,(void *)arg);
	}
	for (uint32_t j=0; j<num_threads;j++){
		pthread_join(threads[j],NULL);
	}
	
	/* End Timing */
	clock_gettime(CLOCK_MONOTONIC, &ts_end);
	double runtime = (ts_end.tv_sec - ts_start.tv_sec) * 1000000000 + (ts_end.tv_nsec - ts_start.tv_nsec);	
	/* Print triangles while using double type, to not lose any decimals */
	double tottriangles=0;
	for (uint32_t x=0;x<row;x++) tottriangles+=nodes[x];
	tottriangles/=3.0;
	printf("c3 %.3f\n",tottriangles);
	/* Print info  */
	printf("%s\n",argv[0]);
	printf("rows %d\n",row);
	printf("Non-zero %d\n",nnz);
	printf("runtime %f s\n" ,runtime/1000000000);
}

/* Void function that check if 1 argument was given */ 
void checkArgs(int argc, char *argv[]){
	/* Check if input arguments is more than 2 */ 
	if (argc>3){
		printf("Too many arguments where given... %d instead of 2 \n",(argc-1));
		printf("Please give in the format\n %s data num_threads\n",argv[0]);
	      	exit(1);
	/* If no argument is given exit */	
	} else if(argc==1){ 
		printf("No argument was given \n");
		exit(1);
	/* Check if only 1 arg was given */
	} else if(argc==2){
		printf("1 arg was given instead of 2 \n");
		printf("Please give in the format\n %s data num_threads\n",argv[0]);
	      	exit(1);
	}
}
/* Convert COO to CSC  */
void coo2csc(
  uint32_t       * const row,       /*!< CSC row start indices */
  uint32_t       * const col,       /*!< CSC column indices */
  uint32_t const * const row_coo,   /*!< COO row indices */
  uint32_t const * const col_coo,   /*!< COO column indices */
  uint32_t const         nnz,       /*!< Number of nonzero elements */
  uint32_t const         n,         /*!< Number of rows/columns */
  uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
) {

  // ----- cannot assume that input is already 0!
  for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


  // ----- find the correct column sizes
  for (uint32_t l = 0; l < nnz; l++)
	  col[col_coo[l] - isOneBased]++;

  // ----- cumulative sum
  for (uint32_t i = 0, cumsum = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = cumsum;
    cumsum += temp;
  }
  col[n] = nnz;
  // ----- copy the row indices to the correct place
  for (uint32_t l = 0; l < nnz; l++) {
    uint32_t col_l;
    col_l = col_coo[l] - isOneBased;

    uint32_t dst = col[col_l];
    row[dst] = row_coo[l] - isOneBased;

    col[col_l]++;
  }
  // ----- revert the column pointers
  for (uint32_t i = 0, last = 0; i < n; i++) {
    uint32_t temp = col[i];
    col[i] = last;
    last = temp;
  }

}
