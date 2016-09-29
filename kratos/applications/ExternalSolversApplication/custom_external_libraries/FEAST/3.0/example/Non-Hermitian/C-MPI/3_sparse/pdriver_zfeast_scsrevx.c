/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Expert Driver Example - CSR Storage
  !!!!!!! solving Ax=ex with A complex-symmetric (non-Hermitian)
  !!!!!!! James KEstyn, Eric Polizzi 2015
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>
#include <mpi.h>

#include "feast.h"
#include "feast_sparse.h"

int main(int argc, char **argv) {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  int loop;
  char UPLO='F'; 

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system4";
  int  N,nnz;
  double *sa;
  int *isa,*jsa;
  /*!!!!!!!!!!!!!!!!! Others */
  double *Zedge,*Zne,*Wne;
  int ccN, *Nedge, *Tedge;
  struct timeval t1, t2;
  int  i,k,err;
  int  M0,M,info;
  double Emid[2],r,epsout,trace[2];
  double *X; //! eigenvectors
  double *E,*res; //! eigenvalue+residual
/*********** MPI *****************************/
int rank,numprocs;
MPI_Init(&argc,&argv); 
MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
/*********************************************/


  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!! read input file in csr format!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  // !!!!!!!!!! form CSR arrays isa,jsa,sa 
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  sa=calloc(2*nnz,sizeof(double)); //factor 2 for complex 
  isa=calloc(N+1,sizeof(int));
  jsa=calloc(nnz,sizeof(int));

  for (i=0;i<=N;i++){
    *(isa+i)=0;
  };
  *(isa)=1;
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d%d%lf%lf\n",&i,jsa+k,sa+2*k,sa+2*k+1);
    *(isa+i)=*(isa+i)+1;
  };
  fclose(fp);
  for (i=1;i<=N;i++){
    *(isa+i)=*(isa+i)+*(isa+i-1);
  };

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if(rank==0) printf("sparse matrix -system4- size %.d\n",N);
  if(rank==0) printf("nnz %d \n",nnz);

  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in sparse format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*!!!!!!!!!!!!  FEAST */
  feastinit(feastparam);
  feastparam[0]=1;  /*change from default value */
  M0=40; // !! M0>=M

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  /*! Create Custom Contour          */
  ccN = 3;     //!! number of pieces that make up contour
  Zedge = calloc(2*ccN,sizeof(double));
  Tedge = calloc(ccN,sizeof(int));
  Nedge = calloc(ccN,sizeof(int));

  /*!!! Example contour - triangle   */
  Zedge[0] = 0.1e0;  Zedge[1] = 0.41e0;  // 1st complex #
  Zedge[2] = 4.2e0;  Zedge[3] = 0.41e0;
  Zedge[4] = 4.2e0;  Zedge[5] = -8.30e0;
  Tedge[0] = 0; Tedge[1] = 0; Tedge[2] = 0;
  Nedge[0] = 6; Nedge[1] = 6; Nedge[2] = 18;

  /*!! Note: user must specify total # of contour points and edit feastparam(8)*/
  feastparam[7]=0;
  for(i=0;i<ccN;i++) feastparam[7] = feastparam[7] + Nedge[i];
  Zne = calloc(2*feastparam[7],sizeof(double)); // Contains the complex valued contour points 
  Wne = calloc(2*feastparam[7],sizeof(double)); // Contains the complex valued integration weights

  /* !! Fill Zne/Wne */
  zfeast_customcontour(&(feastparam[7]),&ccN,Nedge,Tedge,Zedge,Zne,Wne);

  if(rank==0){
     printf("Printing Contour \n");
     for(i=0;i<feastparam[7];i++)  printf("%d: %le + i*%le\n",i,Zne[2*i],Zne[2*i+1]);
  }
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(2*M0,sizeof(double));  // eigenvalues  // factor 2 fopr complex
  res=calloc(M0,sizeof(double));// eigenvectors 
  X=calloc(2*N*M0,sizeof(double));// residual (if needed) // factor 2 for complex


  zfeast_scsrevx(&UPLO,&N,sa,isa,jsa,feastparam,&epsout,&loop,Emid,&r,&M0,E,X,&M,res,&info,Zne,Wne);

  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  if(rank==0) printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0 && rank==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace[0] = 0.0e0;
    trace[1] = 0.0e0;
    for (i=0;i<=M-1;i=i+1){
      trace[0]=trace[0]+ *(E+2*i);
      trace[1]=trace[1] + *(E+2*i+1);
    }	  
    printf("Trace %.15e + i*%.15e\n",trace[0],trace[1] );
    printf("Relative error on the Trace %.15e\n",epsout);
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e+i*%.15e %.15e\n",i,*(E+2*i),*(E+2*i+1),*(res+i));
    }
  }
  MPI_Finalize();
  return 0;
}






