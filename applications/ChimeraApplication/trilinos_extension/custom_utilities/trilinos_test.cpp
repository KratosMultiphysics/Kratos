#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_FEVector.h"
#include <Epetra_Time.h>
#include "EpetraExt_CrsMatrixIn.h"
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#define PRINT() std::cout<<"Rank "<<mpi_rank<<" : "


Epetra_FECrsMatrix* create_epetra_crsmatrix(int numProcs,
                                          int localProc,
                                          int local_n,
                                          bool callFillComplete = true,
                                          bool symmetric = true,
                                          int factor = 1);


int time_matrix_matrix_multiply(Epetra_Comm& Comm,
                                bool verbose);

int main(int argc, char *argv[])
{
    // Initializing MPI and making epetra comm
    (void)MPI_Init(&argc, &argv);
    int ierr;
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // Get current rank and a total num process
    const int mpi_rank = comm.MyPID();
    const int mpi_size = comm.NumProc();
    // const int local_elements_per_rank = 5;
    // const int num_over_lapping_elems = 0;
    // int total_num_local_elems = local_elements_per_rank
    // +num_over_lapping_elems;
    // const int total_global_rows = local_elements_per_rank * mpi_size;
    // std::vector<int> my_local_rows;
    // for(int i=0; i<local_elements_per_rank; ++i)
    // {
    //     my_local_rows.emplace_back( (mpi_rank*local_elements_per_rank)+i );
    // }
    // //defining a map as needed
    // // const Epetra_Map local_row_map_a(-1, my_local_rows.size(), my_local_rows.data(), 0, comm);
    // // const Epetra_Map row_map(-1, local_eq_ids.size(), local_eq_ids.data(), 0, BaseType::mrComm);
    // const Epetra_Map local_row_map_a(total_global_rows, total_num_local_elems, 0, comm);

    // // std::vector<int> my_local_cols;
    // // for(int i=0; i<local_elements_per_rank*mpi_size; ++i)
    // // {
    // //     my_local_cols.emplace_back( i );
    // // }
    // // const Epetra_Map local_col_map_a(-1, my_local_cols.size(), my_local_cols.data(), 0, comm);
    // comm.Barrier();
    // int nnz_per_row = 9;

    // // // Making dummy epetra matrix (square) 10x10
    // // Epetra_FECrsGraph a_graph(Copy, local_row_map_a, local_col_map_a, 10);
    // // // Fill the graph
    // // ierr = a_graph.InsertGlobalIndices(my_local_rows.size(), my_local_rows.data(), my_local_cols.size(), my_local_cols.data());
    // // ierr = a_graph.GlobalAssemble();
    // // // use the graph to make the matrix
    // // Epetra_FECrsMatrix mat_a(Copy,a_graph);
    // Epetra_FECrsMatrix mat_a(Copy, local_row_map_a, nnz_per_row);
    // if(mpi_rank == 0){
    //     int global_row = 0;
    //     int index = 0;
    //     std::vector<double> values{2.0, 5.50, 2.111};
    //     std::vector<int> col_indices{7,5,8};
    //     mat_a.ReplaceGlobalValues(8, col_indices.size(), values.data(), col_indices.data());

    //     global_row = 8;
    //     for(const auto col : col_indices){
    //         mat_a.InsertGlobalValues(global_row, 1, &values[index], &col);
    //         ++index;
    //     }
    //     index = 0;

    //     std::vector<double> values2{3.3, 1.50, 5.111,6.11,8.45};
    //     std::vector<int> col_indices2{0,2,4,6,7};
    //     global_row = 2;
    //     for(const auto col : col_indices2){
    //         mat_a.InsertGlobalValues(global_row, 1, &values2[index], &col);
    //         ++index;
    //     }
    //     index = 0;


    //     std::vector<double> values3{3.3, 1.50, 5.111};
    //     std::vector<int> col_indices3{1,2,4};
    //     global_row = 4;
    //     for(const auto col : col_indices3){
    //         mat_a.InsertGlobalValues(global_row, 1, &values3[index], &col);
    //         ++index;
    //     }
    //     index = 0;
    // }
    // if(mpi_rank == 1){
    //     int global_row = 0;
    //     int index = 0;

    //     std::vector<double> values{2.0, 5.50, 2.111};
    //     std::vector<int> col_indices{5,7,8};
    //     global_row = 7;
    //     for(const auto col : col_indices){
    //         mat_a.InsertGlobalValues(global_row, 1, &values[index], &col);
    //         ++index;
    //     }
    //     index = 0;


    //     std::vector<double> values2{3.3, 1.50, 5.111, 2.11, 8.24, 9.45};
    //     std::vector<int> col_indices2{6,5,8, 0, 2, 3};
    //     global_row = 9;
    //     for(const auto col : col_indices2){
    //         mat_a.InsertGlobalValues(global_row, 1, &values[index], &col);
    //         ++index;
    //     }
    //     index = 0;
    // }
    // // mat_a.FillComplete();
    // {
    //     EpetraExt::RowMatrixToMatrixMarketFile("amat.mm", mat_a);
    // }
    // PRINT()<<"Setup of matrix A finished !"<<std::endl;

    // // // Making dummy epetra matrix 10x10
    // // std::vector<int> my_local_rows_b = my_local_rows;
    // // std::vector<int> my_local_cols_b{1,2,4,5,6,7,8};
    // // const Epetra_Map local_row_map_b(-1, my_local_rows_b.size(), my_local_rows_b.data(), 0, comm);
    // if(mpi_rank == 0)
    //     total_num_local_elems-=2;
    // if(mpi_rank == 1)
    //     total_num_local_elems+=2;

    // const Epetra_Map local_row_map_b(total_global_rows, total_num_local_elems, 0, comm);
    // // const Epetra_Map local_col_map_b(-1, my_local_cols_b.size(), my_local_cols_b.data(), 0, comm);

    // // Epetra_FECrsGraph b_graph(Copy, local_row_map_b, local_col_map_b, 10);
    // // ierr = b_graph.InsertGlobalIndices(my_local_rows_b.size(), my_local_rows_b.data(), my_local_cols_b.size(), my_local_cols_b.data());
    // // ierr = b_graph.GlobalAssemble();

    // // Epetra_FECrsMatrix mat_b(Copy,b_graph);
    // Epetra_FECrsMatrix mat_b(Copy, local_row_map_b, nnz_per_row);
    // if(mpi_rank == 1){
    //     int global_row = 0;
    //     int index = 0;
    //     std::vector<double> values{2.0, 5.50, 2.111};
    //     std::vector<int> col_indices{5,6,8};

    //     global_row = 6;
    //     for(const auto col : col_indices){
    //         mat_b.InsertGlobalValues(global_row, 1, &values[index], &col);
    //         ++index;
    //     }
    //     index = 0;

    //     std::vector<double> values2{3.3, 1.50};
    //     std::vector<int> col_indices2{5,7};
    //     global_row = 8;
    //     for(const auto col : col_indices2){
    //         mat_b.InsertGlobalValues(global_row, 1, &values2[index], &col);
    //         ++index;
    //     }
    //     index = 0;
    // }

    // if(mpi_rank == 0){
    //     int global_row = 0;
    //     int index = 0;
    //     std::vector<double> values{2.0, 5.50, 2.111};
    //     std::vector<int> col_indices{5,6,8};
    //     global_row = 5;
    //     for(const auto col : col_indices){
    //         mat_b.InsertGlobalValues(global_row, 1, &values[index], &col);
    //         ++index;
    //     }
    //     index = 0;

    //     std::vector<double> values2{3.3, 1.50, 5.111};
    //     std::vector<int> col_indices2{1,2,4};
    //     global_row = 2;
    //     for(const auto col : col_indices2){
    //         mat_b.InsertGlobalValues(global_row, 1, &values2[index], &col);
    //         ++index;
    //     }
    //     index = 0;


    //     std::vector<double> values3{3.3, 1.50, 5.111};
    //     std::vector<int> col_indices3{1,2,4};
    //     global_row = 4;
    //     for(const auto col : col_indices3){
    //         mat_b.InsertGlobalValues(global_row, 1, &values3[index], &col);
    //         ++index;
    //     }
    //     index = 0;
    // }

    // mat_b.FillComplete();
    // {
    //     EpetraExt::RowMatrixToMatrixMarketFile("bmat.mm", mat_b);
    // }
    // PRINT()<<"Setup of matrix B finished !"<<std::endl;

    // // Now multiplying c = b'*a
    // Epetra_FECrsMatrix mat_c(Copy,mat_b.RowMap(), 0);
    // // c <- B'*A
    // ierr = EpetraExt::MatrixMatrix::Multiply(mat_b, true, mat_a, false, mat_c, false, true);
    // PRINT()<<"Error when multiplying B'A is : "<<ierr<<std::endl;
    // mat_c.FillComplete();
    // {
    //     EpetraExt::RowMatrixToMatrixMarketFile("cmat.mm", mat_c);
    // }
    // PRINT()<<"Setup of matrix C finished !"<<std::endl;
    // // d <- C*B
    // PRINT()<<" ONE "<<std::endl;
    // Epetra_FECrsMatrix mat_d(Copy, mat_b.RowMap(), 0);
    // ierr = EpetraExt::MatrixMatrix::Multiply(mat_c, false, mat_b, false, mat_d, false, true);
    // PRINT()<<"Error when multiplying C*B is : "<<ierr<<std::endl;
    // // mat_d.FillComplete();
    // {
    //     EpetraExt::RowMatrixToMatrixMarketFile("dmat.mm", mat_a);
    // }
    // PRINT()<<"Setup of matrix D finished !"<<std::endl;

    // // ierr = time_matrix_matrix_multiply(comm, true);





    /// New TEST

    {
      int total_rows = 10;

      if(mpi_rank == 0){
          int ierr = 0;
          int local_rows = 3;
          std::vector<int> row_indices{0,1,2,3,4};
          std::vector<int> col_indices{0,1,2,3,4};

          const Epetra_Map row_map_n(-1, local_rows, 0, comm);
          // const Epetra_Map row_map_n(-1, row_indices.size(),row_indices.data(), 0, comm);
          const Epetra_Map col_map_n(-1, col_indices.size(),col_indices.data(), 0, comm);

          // Epetra_FECrsGraph n_graph(Copy, row_map_n, col_map_n, 10);

          // std::vector<int> test_row_indices{0,1,2,5,6};
          // std::vector<int> test_col_indices{0,1,2,5,6};

          // ierr = n_graph.InsertGlobalIndices(test_row_indices.size(), test_row_indices.data(),
          //                                         test_col_indices.size(), test_col_indices.data());
          // if(ierr != 0)
          //     std::cout << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
          //     << std::endl;

          // ierr = n_graph.GlobalAssemble();
          // if(ierr != 0)
          //     std::cout << ": Epetra failure in Graph.GlobalAssemble. Error code: " << ierr
          //     << std::endl;

          // Epetra_FECrsMatrix mat_b(Copy,n_graph);
          Epetra_FECrsMatrix mat_b(Copy, row_map_n, 10);
          // Epetra_FECrsMatrix mat_b(Copy, row_map_n, col_map_n, 10);

          std::vector<double> values{3.3, 1.50, 5.111};
          std::vector<int> input_col_indices{1,2,7};
          int global_row = 2;
          int index = 0;
          for(const auto col : input_col_indices){
              ierr = mat_b.InsertGlobalValues(global_row, 1, &values[index], &col);
              ++index;
              if(ierr != 0)
                  PRINT() << ": Epetra failure in matrix.InsertGlobalValues. Error code: " << ierr
                  << std::endl;
          }
          mat_b.FillComplete();
          {
              EpetraExt::RowMatrixToMatrixMarketFile("bmat.mm", mat_b);
          }
          PRINT()<<"Setup of matrix B finished !"<<std::endl;
      }

      if(mpi_rank == 1){
          int ierr = 0;
          int local_rows = 7;
          std::vector<int> row_indices{2,5,6,7,8,9};
          std::vector<int> col_indices{2,5,6,7,8,9};

          const Epetra_Map row_map_n(-1, local_rows, 0, comm);
          // const Epetra_Map row_map_n(-1, row_indices.size(),row_indices.data(), 0, comm);
          const Epetra_Map col_map_n(-1, col_indices.size(),col_indices.data(), 0, comm);

          // Epetra_FECrsGraph n_graph(Copy, row_map_n, col_map_n, 10);

          // std::vector<int> test_row_indices{7,8,2,5,6};
          // std::vector<int> test_col_indices{7,8,2,5,6};

          // ierr = n_graph.InsertGlobalIndices(test_row_indices.size(), test_row_indices.data(),
          //                                         test_col_indices.size(), test_col_indices.data());
          // if(ierr != 0)
          //     std::cout << ": Epetra failure in Graph.InsertGlobalIndices. Error code: " << ierr
          //     << std::endl;

          // ierr = n_graph.GlobalAssemble();
          // if(ierr != 0)
          //     std::cout << ": Epetra failure in Graph.GlobalAssemble. Error code: " << ierr
          //     << std::endl;

          // Epetra_FECrsMatrix mat_b(Copy,n_graph);
          // Epetra_FECrsMatrix mat_b(Copy, row_map_n, col_map_n, 10);
          Epetra_FECrsMatrix mat_b(Copy, row_map_n, 10);

          std::vector<double> values{-1.3, 12.50, 0.147};
          std::vector<int> input_col_indices{5,6,2};
          int global_row = 7;
          int index = 0;
          for(const auto col : input_col_indices){
              ierr = mat_b.InsertGlobalValues(global_row, 1, &values[index], &col);
              ++index;
              if(ierr != 0)
                  PRINT() << ": Epetra failure in matrix.InsertGlobalValues. Error code: " << ierr
                  << std::endl;
          }
          mat_b.FillComplete();
          {
              EpetraExt::RowMatrixToMatrixMarketFile("bmat.mm", mat_b);
          }
          PRINT()<<"Setup of matrix B finished !"<<std::endl;
      }

    }




















    // MPI_Finalize after you are done using MPI.
    (void) MPI_Finalize ();
    return 0;
}

Epetra_FECrsMatrix* create_epetra_crsmatrix(int numProcs,
                                          int localProc,
                                          int local_n,
                                          bool callFillComplete,
                                          bool symmetric,
                                          int factor)
{
  (void)localProc;
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  int global_num_rows = numProcs*local_n;

  Epetra_Map rowmap(global_num_rows, local_n, 0, comm);

  int nnz_per_row = 9;
  Epetra_FECrsMatrix* matrix =
    new Epetra_FECrsMatrix(Copy, rowmap, nnz_per_row);

  // Add  rows one-at-a-time
  double negOne = -1.0;
  double posTwo = 2.0;
  double val_L = symmetric ? negOne : -0.5;

  for (int i=0; i<local_n; i++) {
    int GlobalRow = matrix->GRID(i);
    int RowLess1 = GlobalRow - 1*factor;
    int RowPlus1 = GlobalRow + 1*factor;
    int RowLess5 = GlobalRow - 5*factor;
    int RowPlus5 = GlobalRow + 5*factor;
    int RowLess9 = GlobalRow - 9*factor;
    int RowPlus9 = GlobalRow + 9*factor;
    int RowLess24 = GlobalRow - 24*factor;
    int RowPlus24 = GlobalRow + 24*factor;
    int RowLess48 = GlobalRow - 48*factor;
    int RowPlus48 = GlobalRow + 48*factor;

//    if (!symmetric) RowLess5 -= 2;

    if (RowLess48>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess48);
    }
    if (RowLess24>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess24);
    }
    if (RowLess9>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess9);
    }
    if (RowLess5>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess5);
    }
    if (RowLess1>=0) {
      matrix->InsertGlobalValues(GlobalRow, 1, &val_L, &RowLess1);
    }
    if (RowPlus1<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus1);
    }
    if (RowPlus5<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus5);
    }
    if (RowPlus9<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus9);
    }
    if (RowPlus24<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus24);
    }
    if (RowPlus48<global_num_rows) {
      matrix->InsertGlobalValues(GlobalRow, 1, &negOne, &RowPlus48);
    }

    matrix->InsertGlobalValues(GlobalRow, 1, &posTwo, &GlobalRow);
  }

  if (callFillComplete) {
    int err = matrix->FillComplete();
    if (err != 0) {
      std::cout << "create_epetra_matrix: error in matrix.FillComplete()"
         <<std::endl;
    }
  }

  return(matrix);
}

int time_matrix_matrix_multiply(Epetra_Comm& Comm, bool verbose)
{

  const int magic_num = 3000;
  // 2009/02/23: rabartl: If you are going to do a timing test you need to
  // make this number adjustable form the command-line and you need to put in
  // a real test that compares against hard numbers for pass/fail.

  int localn = magic_num/Comm.NumProc();

  Epetra_FECrsMatrix* A = create_epetra_crsmatrix(Comm.NumProc(),
                                                Comm.MyPID(),
                                                localn);

  Epetra_FECrsMatrix* B = create_epetra_crsmatrix(Comm.NumProc(),
                                                Comm.MyPID(),
                                                localn,
                                                true, false, 2);

    {
        EpetraExt::RowMatrixToMatrixMarketFile("amat.mm", *A);
        EpetraExt::RowMatrixToMatrixMarketFile("bmat.mm", *B);
    }

  Epetra_FECrsMatrix* C = new Epetra_FECrsMatrix(Copy, A->RowMap(), 0);
  Epetra_FECrsMatrix* D = new Epetra_FECrsMatrix(Copy, C->RowMap(), 0);

  Epetra_Time timer(Comm);
  double start_time = timer.WallTime();

  int err = EpetraExt::MatrixMatrix::Multiply(*B, true, *A, false, *C);

  if (err != 0) {
    std::cout << "B'*A returned err=="<<err<<std::endl;
    return(err);
  }

  int globaln = localn*Comm.NumProc();
  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", C = B'*A, time: "
       << timer.WallTime()-start_time << std::endl;
  }

  C->FillComplete();

  start_time = timer.WallTime();

  err = EpetraExt::MatrixMatrix::Multiply(*C, false, *B, false, *D);

  if (err != 0) {
    std::cout << "C*A returned err=="<<err<<std::endl;
    return(err);
  }
  D->FillComplete();

  if (verbose) {
    std::cout << "size: " << globaln << "x"<<globaln<<", D = C*A, time: "
       << timer.WallTime()-start_time << std::endl;
  }

    {
        EpetraExt::RowMatrixToMatrixMarketFile("cmat.mm", *C);
        EpetraExt::RowMatrixToMatrixMarketFile("dmat.mm", *D);
    }

  delete C;
  delete D;

  return(err);
}