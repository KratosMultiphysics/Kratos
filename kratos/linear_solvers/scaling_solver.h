/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_SCALING_SOLVER_H_INCLUDED )
#define  KRATOS_SCALING_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class ScalingSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType,  TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ScalingSolver
      KRATOS_CLASS_POINTER_DEFINITION(ScalingSolver);

      typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ScalingSolver()
      {
      }

      ScalingSolver(typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer p_linear_solver,
                    bool symmetric_scaling )
      {
          msymmetric_scaling = true;
          mp_linear_solver = p_linear_solver;
      }

      /// Copy constructor.
      ScalingSolver(const ScalingSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~ScalingSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      ScalingSolver& operator=(const ScalingSolver& Other)
      {
        BaseType::operator=(Other);
	return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  if(IsNotConsistent(rA, rX, rB))
	    return false;

          VectorType scaling_vector(rX.size());

          //obtain the scaling matrix
	  GetScalingWeights(rA,scaling_vector);

	  //scale system
          if(msymmetric_scaling == false)
          {
              KRATOS_ERROR(std::logic_error,"not yet implemented","")
          }
          else
          {
              #pragma omp parallel for
              for(int i=0; i< static_cast<int>(scaling_vector.size()); i++)
                  scaling_vector[i] = sqrt(scaling_vector[i]);
              
              SymmetricScaling(rA,scaling_vector);
          }

          //solve the problem
          bool is_solved = mp_linear_solver->Solve(rA,rX,rB);

	  //backscale the solution
          if(msymmetric_scaling == true)
          {
               #pragma omp parallel for
              for(int i=0; i< static_cast<int>(scaling_vector.size()); i++)
                  rX[i] /= scaling_vector[i];
          }

	  return is_solved;
	}
      
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "Composite Linear Solver. Uses internally the following linear solver " << mp_linear_solver->Info();
	  return  buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  BaseType::PrintData(rOStream);
	}
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
      typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mp_linear_solver;
      bool msymmetric_scaling;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        static void SymmetricScaling( SparseMatrixType& A, const VectorType& aux)
        {
            typedef  unsigned int size_type;
            typedef  double value_type;

            //create partition
            OpenMPUtils::PartitionVector partition;
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(A.size1(),number_of_threads,  partition);
            //parallel loop

            #pragma omp parallel
            {
                int thread_id = OpenMPUtils::ThisThread();
                int number_of_rows = partition[thread_id+1] - partition[thread_id];
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
                 typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
                 typename compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::iterator value_begin = A.value_data().begin()+*row_iter_begin;

                perform_matrix_scaling(    number_of_rows,
                                    row_iter_begin,
                                    index_2_begin,
                                    value_begin,
                                    partition[thread_id], 
                                    aux
                                );
            }
        }

	/**
	 * calculates partial product resetting to Zero the output before
	 */
        static void perform_matrix_scaling(
                int number_of_rows,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator row_begin,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::iterator index2_begin,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::iterator value_begin,
		unsigned int output_begin_index,
                const VectorType& weights
                 )
        {
            int row_size;
            typename SparseMatrixType::index_array_type::const_iterator row_it = row_begin;
            int kkk = output_begin_index;
            for(int k = 0; k < number_of_rows; k++)
            {
                row_size= *(row_it+1)-*row_it;
                row_it++;
                const typename TDenseSpaceType::DataType row_weight = weights[kkk++];

                for(int i = 0; i<row_size; i++)
                {
                    const typename TDenseSpaceType::DataType col_weight = weights[*index2_begin];
                    typename TDenseSpaceType::DataType t = (*value_begin);
                    t /= (row_weight*col_weight);
                    (*value_begin) = t; //check if this is correcct!!
                    value_begin++;
		    index2_begin++;
                }

            }
        }


        static void GetScalingWeights( const SparseMatrixType& A, VectorType& aux)
        {
            typedef  unsigned int size_type;
            typedef  double value_type;

            //create partition
            OpenMPUtils::PartitionVector partition;
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::DivideInPartitions(A.size1(),number_of_threads,  partition);
            //parallel loop

            #pragma omp parallel
            {
                int thread_id = OpenMPUtils::ThisThread();
                int number_of_rows = partition[thread_id+1] - partition[thread_id];
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
                 typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
                 typename compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;

		 GS2weights(    number_of_rows,
                                    row_iter_begin,
                                    index_2_begin,
                                    value_begin,
                                    partition[thread_id],
                                    aux
                                );
            }
        }

	/**
	 * calculates partial product resetting to Zero the output before
	 */
        static void GS2weights(
                int number_of_rows,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator row_begin,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::index_array_type::const_iterator index2_begin,
                typename compressed_matrix<typename TDenseSpaceType::DataType>::value_array_type::const_iterator value_begin,
		unsigned int output_begin_index,
                VectorType& weights
                 )
        {
	    int row_size;
            typename SparseMatrixType::index_array_type::const_iterator row_it = row_begin;
            int kkk = output_begin_index;
            for(int k = 0; k < number_of_rows; k++)
            {
                row_size= *(row_it+1)-*row_it;
                row_it++;
                double t = typename TDenseSpaceType::DataType();

                for(int i = 0; i<row_size; i++)
		{
		  double tmp = *value_begin;
                    t += tmp*tmp;
		    value_begin++;
		}
                t = sqrt(t);
                weights[kkk++] = t;
            }
        }


      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
        
      ///@}    
        
    }; // Class ScalingSolver

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::istream& operator >> (std::istream& IStream, 
				      ScalingSolver<TSparseSpaceType, TDenseSpaceType,
				       TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const ScalingSolver<TSparseSpaceType, TDenseSpaceType,
				       TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_SCALING_SOLVER_H_INCLUDED  defined


