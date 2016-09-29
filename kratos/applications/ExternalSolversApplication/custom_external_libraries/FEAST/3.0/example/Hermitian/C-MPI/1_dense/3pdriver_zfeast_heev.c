/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! FEAST Driver dense example !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! solving Ax=ex with A complex Hermitian --- A dense matrix!!!!!!!!
  !!!!!!! Using two search intervals
  !!!!!!! by Eric Polizzi- 2009-2012!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  double epsout;
  int loop;
  char UPLO='F'; // ! 'L' or 'U' also fine

  /*!!!!!!!!!!!!!!!!! Matrix declaration variable */
  FILE *fp;
  char name[]="../../system2";
  int  N,nnz;
  double *A;

  /*!!!!!!!!!!!!!!!!! Others */
  struct timeval t1, t2;
  int  i,j,k,n2,err;
  int  M0,M,info;
  double Emin,Emax,trace;
  double *X; //! eigenvectors
  double *E,*res; //! eigenvalue+residual

/*********** MPI *****************************/
int lrank,lnumprocs,color,key;
int rank,numprocs;
MPI_Comm NEW_COMM_WORLD;
MPI_Init(&argc,&argv); 
MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 
MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
/*********************************************/
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!Read Coordinate format and convert to dense format
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  fp = fopen (name, "r");
  err=fscanf (fp, "%d%d\n",&N,&nnz);
  n2=2*N*N;  // factor 2 because of complex number
  A=calloc(n2,sizeof(double));
  for (i=0;i<=n2-1;i++){
    *(A+i)=(double) 0.0;
  };
  for (k=0;k<=nnz-1;k++){
    err=fscanf(fp,"%d %d",&i,&j);
    err=fscanf(fp,"%lf%lf\n",A+(j-1)*2*N+2*i-2,A+(j-1)*2*N+2*i-1);
  };
  fclose(fp);
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! INFORMATION ABOUT MATRIX !!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    if (rank==0){printf("dense matrix -system2- size %.d\n",N);}
  

  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!! FEAST in dense format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Definition of the two intervals 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  if (rank<=numprocs/2-1) {
    color=1;} // first interval
  else {
    color=2; //! second interval
  }

  //!!!!!!!!!!!!!!!!! create new_mpi_comm_world
 key=0;
 MPI_Comm_split(MPI_COMM_WORLD,color,key,&NEW_COMM_WORLD);
 MPI_Comm_rank(NEW_COMM_WORLD,&lrank);
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  /*!!! search interval [Emin,Emax] including M eigenpairs*/
 if (color==1) { // 1st interval
  Emin=(double) -0.35;
  Emax=(double) 0.0;
  M0=40; // !! M0>=M
  }
 else if(color==2){ // 2nd interval
  Emin=(double) 0.0;
  Emax=(double) 0.23;
  M0=40; // !! M0>=M
  }

//!!!!!!!!!!!!!!!!!! RUN INTERVALS in PARALLEL
 gettimeofday(&t1,NULL);


  /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
  E=calloc(M0,sizeof(double));  // eigenvalues
  res=calloc(M0,sizeof(double));// residual
  X=calloc(2*N*M0,sizeof(double));// eigenvectors // factor 2 because of complex number


  /*!!!!!!!!!!!!  FEAST */
  feastinit(feastparam);
  feastparam[8]=NEW_COMM_WORLD;  /*change from default value */
  zfeast_heev(&UPLO,&N,A,&N,feastparam,&epsout,&loop,&Emin,&Emax,&M0,E,X,&M,res,&info);

  gettimeofday(&t2,NULL);
  /*!!!!!!!!!! REPORT !!!!!!!!!*/
if (lrank==0) {
  printf("interval # %d\n",color);
  printf("FEAST OUTPUT INFO %d\n",info);
  if (info==0) {
    printf("*************************************************\n");
    printf("************** REPORT ***************************\n");
    printf("*************************************************\n");
    MPI_Comm_size(NEW_COMM_WORLD,&lnumprocs);
    printf("# of processors %d \n",lnumprocs);
    printf("SIMULATION TIME %f\n",(t2.tv_sec-t1.tv_sec)*1.0+(t2.tv_usec-t1.tv_usec)*0.000001);
    printf("# Search interval [Emin,Emax] %.15e %.15e\n",Emin,Emax);
    printf("# mode found/subspace %d %d \n",M,M0);
    printf("# iterations %d \n",loop);
    trace=(double) 0.0;
    for (i=0;i<=M-1;i=i+1){
      trace=trace+*(E+i);
    }	  
    printf("TRACE %.15e\n", trace);
    printf("Relative error on the Trace %.15e\n",epsout );
    printf("Eigenvalues/Residuals\n");
    for (i=0;i<=M-1;i=i+1){
      printf("   %d %.15e %.15e\n",i,*(E+i),*(res+i));
    }
  }
 }


MPI_Finalize(); /************ MPI ***************/
  return 0;
}





