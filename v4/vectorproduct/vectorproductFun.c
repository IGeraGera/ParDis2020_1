#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "mmio.h"
/* C = A*v */
void CSCvectorMult(
  uint32_t const * const A_row,
  uint32_t const * const A_col,
  uint32_t const * const A_val,
  uint32_t       * const C,
  uint32_t const * const vector,
  uint32_t const         n

);
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

	/* Declaration of COO format arrays  */ 
	uint32_t *coo_r,*coo_c,*val;

    	/* Allocate memory for matrices */
   	coo_r = (uint32_t *) malloc(2*nz * sizeof(uint32_t ));
    	coo_c = (uint32_t *) malloc(2*nz * sizeof(uint32_t ));
    	val  = (uint32_t *) malloc(2* nz * sizeof(uint32_t ));
	/* Iterate .mtx and scan the values  */
	for (uint32_t i=0; i<2*nz; i+=2){
		fscanf(f,"%d %d \n",&coo_r[i],&coo_c[i]);
		coo_r[i]--;
		coo_c[i]--;
		coo_c[i+1]=coo_r[i];
		coo_r[i+1]=coo_c[i];
		val[i]=1;
		val[i+1]=1;
	}
	/* Close the file  */
	fclose(f);
	/* Declaration of CSC format arrays  */ 
	uint32_t const row = M;
	uint32_t const nnz = nz*2;
	uint32_t isOneBased =0;
    	 /*Allocate memory for csc  matrices */
  	uint32_t * csc_r = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
  	uint32_t * csc_c = (uint32_t *)malloc((row+1) * sizeof(uint32_t));
	/* Make the conversion from CSC to COO */	
	coo2csc(csc_r,csc_c,coo_r,coo_c,nnz,row,isOneBased);
	/* Allocate memory and initialize dense vector v and C*/
  	uint32_t * v = (uint32_t *)malloc((row) * sizeof(uint32_t));
  	uint32_t * C = (uint32_t *)malloc((row) * sizeof(uint32_t));
	for (int x = 0; x < row; x++) {
		v[x] = 2; 
		C[x] = 0;
	}
	CSCvectorMult(csc_r,csc_c,val,C,v,row);
}
/* C = A*v */
void CSCvectorMult(
  uint32_t const * const A_row,
  uint32_t const * const A_col,
  uint32_t const * const A_val,
  uint32_t       * const C,
  uint32_t const * const vector,
  uint32_t const         n

){
	for (uint32_t i=0; i<n;i++){
		uint32_t ptr1start=A_col[i];
		uint32_t ptr1end=A_col[i+1];
		for(uint32_t ptr1=ptr1start ;ptr1<ptr1end; ptr1++){
			uint32_t j = A_row[ptr1];
			C[i]=C[i] + A_val[ptr1] * vector[j];
		}
		printf("C[%d] = %d\n",i,C[i]);
	}
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
