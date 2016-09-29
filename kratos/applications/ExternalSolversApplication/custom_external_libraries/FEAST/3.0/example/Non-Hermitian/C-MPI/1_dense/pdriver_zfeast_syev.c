/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver Example -Dense Storage 
  !!!!!!! solving Ax=ex with A complex-symmetric (non-Hermitian)
  !!!!!!! James Kestyn, Eric Polizzi 2015
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#include <stdio.h> 
#include <stdlib.h> 
#include <sys/time.h>
#include <mpi.h> 

#include "feast.h"
#include "feast_dense.h"
int main(int argc, char **argv) {
  /*!!!!!!!!!!!!!!!!! Feast declaration variable */
  int  feastparam[64]; 
  int loop;
  char UPLO='F'; // ! 'L' or 'U' also fine

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system4";
  int  N,LDA,LDB,nnz;
  double *A,*B;
  double *Zedge,*Zne,*Wne;
  int ccN, *Nedge, *Tedge;
  double *Emid,epsout,*trace;
  double r;
  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,j,k,n2,err;
  int  M0,M,info;
  double *E,*XR,*XL; //! eigenvectors
  double *resr,*resl; //! eigenvalue+resridual
/*********** MPI *****************************/
int rank,numprocs;
MPI_Init(&argc,&argv); 
MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
/*********************************************/

  Emid = calloc(2,sizeof(double));
  trace = calloc(2,sizeof(double));
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!Read Coordinate format and convert to dense format
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  LDA=N;
  n2=2*N*N;  // factor 2 because of complex number
  A=calloc(n2,sizeof(double));
  B=calloc(n2,sizeof(double));
  for (i=0;i<=n2-1;i++){
    *(A+i)=(double) 0.0;
    *(B+i)=(double) 0.0;
  };
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d %d",&i,&j);
    err=fscanf(fp,"%lf%lf\n",A+(j-1)*2*N+2*i-2,A+(j-1)*2*N+2*i-1);
  };
  for (i=0;i<N;i++){
    B[i+N*i]=(double) 1.0;
  };
  fclose(fp);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if(rank==0) printf("dense matrix -system4- size %.d\n",N);
  
  gettimeofday(&t1,NULL);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*!!! search interval [Emin,Emax] including M eigenpairs*/
  Emid[0] = 4.0e0;
  Emid[1] = 0.0e0;
  r = 3.0e0;
  M0=40; // !! M0>=M

  feastinit(feastparam);
  feastparam[0]=1;  /*change from default value */
  
  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(2*M0,sizeof(double));  // eigenvalues
  XR=calloc(2*N*M0,sizeof(double));// right eigenvector // factor 2 because of complex number
  resr=calloc(M0,sizeof(double));// eigenvector residual 

  printf("Enter FEAST - %d\n",rank);
  /*!!!!!!!!!!!!  FEAST */
  zfeast_syev(&UPLO,&N,A,&N,feastparam,&epsout,&loop,Emid,&r,&M0,E,XR,&M,resr,&info);

  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
  info=0;
  if(rank==0) printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0 && rank==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace[0]=(double) 0.0;
    trace[1]=(double) 0.0;
    for (i=0;i<=M-1;i=i+1){
      trace[0]=trace[0]+*(E+2*i);
      trace[1]=trace[1]+*(E+2*i+1);
    }	  
    printf("TRACE %.15e + i*%.15e \n", trace[0],trace[1]);
    printf("Relative error on the Trace %.15e\n",epsout );
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e+i*%.15e \t %.15e\n",i,*(E+2*i),*(E+2*i+1),*(resr+i));
    }
  }
  MPI_Finalize();
  return 0;
}





