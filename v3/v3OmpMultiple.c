/* ATTENTION: In the current form the input must be a coordinate  
 *
 * This is coded for the purposes of the parallel ann distibuted course.
 * In this code the triangles are detected by iterating in a CSC data form.
 * A the begining an argument in the form of .mtx (COO) is given. Then the program
 * converts the COO form onto CSC and finally it iterates through the CSC data structure
 * by following edges. 
 * */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <omp.h>
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

int main(int argc, char *argv[]){
	/* Define mutex */
	omp_lock_t mut;
	/* Timer structs */
	struct timespec ts_start;
	struct timespec ts_end;
	/* Variable declaration */
	FILE *f;
	int M,N,nz,ret;
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
	/* Declaration of COO format arrays  */ 
	uint32_t *coo_r,*coo_c,*val;
    	/* Allocate memory for coo matrices */
   	coo_r = (uint32_t *) malloc(nz * sizeof(uint32_t ));
    	coo_c = (uint32_t *) malloc(nz * sizeof(uint32_t ));
    	//val  = (uint32_t *) malloc(nz * sizeof(uint32_t ));
	/* Iterate .mtx and scan the values  */
	for (uint32_t i=0; i<nz; i++){
		fscanf(f,"%d %d \n",&coo_r[i],&coo_c[i]);
		coo_r[i]--;
		coo_c[i]--;
		//val[i]=1;
	}
	/* Close the file  */
	fclose(f);
	/* Print the .mtx file  */ 
/*    	for (int i=0; i<nz; i++)
        	fprintf(stdout, "%d %d \n", coo_r[i]+1, coo_c[i]+1); 

	printf("\n");
*/
	/* Declaration of CSC format arrays  */ 
	uint32_t const row = M;
	uint32_t const nnz = nz;
	uint32_t isOneBased =0;
    	 /*Allocate memory for csc  matrices */
  	uint32_t * csc_r = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
  	uint32_t * csc_c = (uint32_t *)malloc((row+1) * sizeof(uint32_t));
	/* Make the conversion */	
	coo2csc(csc_r,csc_c,coo_r,coo_c,nnz,row,isOneBased);
	/* Nodes Array and zero init */
	uint32_t *nodes= (uint32_t *)malloc(row * sizeof(uint32_t));
	for (uint32_t i =0; i<row;i++) nodes[i]=0;
	/* Init Mutex */
	omp_init_lock(&mut);
	/* Starting timing */
	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	/* Triangle detection --- Edge following */
	# pragma omp parallel 
	{
	# pragma omp for schedule(dynamic,300)
	for(uint32_t i=0; i<row-2;i++){
		uint32_t ptr1start=csc_c[i];
		uint32_t ptr1end=csc_c[i+1];
		# pragma omp parallel
		{
		#  pragma omp for schedule(dynamic,300)	
		for(uint32_t ptr1=ptr1start ;ptr1<ptr1end; ptr1++){
			uint32_t j = csc_r[ptr1];
			uint32_t ptr2start=csc_c[j];
			uint32_t ptr2end =csc_c[j+1];
			for (uint32_t ptr2=ptr2start ;ptr2<ptr2end; ptr2++){
				uint32_t k = csc_r[ptr2];
				for(uint32_t check = ptr1;check<ptr1end;check++){
					if (csc_r[check]==k){
						# pragma omp critical
						{
						nodes[i]++;
						nodes[j]++;
						nodes[k]++;
						}
					}
				}

			}
		}
		}
	}
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
	/* Check if input argument is 1 */ 
	if (argc>2){
		printf("Too many arguments where given... %d instead of 1 \n",(argc-1));
	      	exit(1);
	/* If no argument is given exit */	
	} else if(argc==1){ 
		printf("No argument was given \n");
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
