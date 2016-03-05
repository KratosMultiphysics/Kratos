/*
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: rrossi $  				 *
* Date:                $Date: 2008-12-05 13:40:39 $			 *
* Revision:            $Revision: 1.2 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_MULTILEVEL_SOLVER_H_INCLUDED )
#define  KRATOS_MULTILEVEL_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"


namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MultiLevelSolver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of MultiLevelSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(MultiLevelSolver);

    enum ScalingType {NoScaling, LeftScaling};

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::SparseMatrixPointerType SparseMatrixPointerType;
  
    typedef typename BaseType::VectorPointerType VectorPointerType;

    typedef typename boost::shared_ptr< ML_Epetra::MultiLevelPreconditioner > MLPreconditionerPointerType;

    /**
     * Default constructor
     */
    MultiLevelSolver(Teuchos::ParameterList& aztec_parameter_list, Teuchos::ParameterList& ml_parameter_list, double tol, int nit_max)
    {
        mAztecParameterList = aztec_parameter_list;
        mMLParameterList = ml_parameter_list;
        mtol = tol;
        mmax_iter = nit_max;
        mScalingType = LeftScaling;
        mMLPrecIsInitialized = false;
        mReformPrecAtEachStep = true;
    }

    /**
     * Destructor
     */
    virtual ~MultiLevelSolver()
    {
    }

    void SetScalingType(ScalingType val)
    {
      mScalingType = val;
    }

    ScalingType GetScalingType()
    {
      return mScalingType;
    }

    void SetReformPrecAtEachStep(bool val)
    {
      mReformPrecAtEachStep = val;
    }

    void ResetPreconditioner()
    {
      mpMLPrec.reset();
      mMLPrecIsInitialized = false;
    }

    void Clear()
    {
      ResetPreconditioner();
    }

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result in rX.
     * rX is also the initial guess for iterative methods.
     * @param rA. System matrix.
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        KRATOS_TRY
        Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);

        if (this->GetScalingType() == LeftScaling)
        {
          // don't use this with conjugate gradient
          // it destroys the symmetry
          Epetra_Vector scaling_vect(rA.RowMap());
          rA.InvColSums(scaling_vect);
          AztecProblem.LeftScale(scaling_vect);
        }
	
        mMLParameterList.set("PDE equations", mndof);

        // create the preconditioner now. this is expensive.
        // the preconditioner stores a pointer to the system
        // matrix. if the system matrix is freed from heap
        // before the preconditioner, a memory error can occur
        // when the preconditioner is freed. the strategy
        // should take care to Clear() the linear solver
        // before the system matrix.
        if (mReformPrecAtEachStep == true ||
            mMLPrecIsInitialized == false)
        {
          MLPreconditionerPointerType tmp(new ML_Epetra::MultiLevelPreconditioner(rA, mMLParameterList, true));
          mpMLPrec.swap(tmp);
          mMLPrecIsInitialized = true;
        }

        // create an AztecOO solver
        AztecOO aztec_solver(AztecProblem);
        aztec_solver.SetParameters(mAztecParameterList);

        // set preconditioner and solve
        aztec_solver.SetPrecOperator(&(*mpMLPrec));
 
        aztec_solver.Iterate(mmax_iter, mtol);

        return true;
        KRATOS_CATCH("");
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {

        return false;
    }
    
    	/** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {
      int old_ndof = -1;
      unsigned int old_node_id = rdof_set.begin()->Id();
      int ndof=0;
      for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
      {
        //			if(it->EquationId() < rdof_set.size() )
        //			{
        unsigned int id = it->Id();
        if(id != old_node_id)
        {
          old_node_id = id;
          if(old_ndof == -1) old_ndof = ndof;
          else if(old_ndof != ndof) //if it is different than the block size is 1
          {
            old_ndof = -1;
            break;
          }
          
          ndof=1;
        }
        else
        {
          ndof++;
        }
        //			}
      }
      
      r_model_part.GetCommunicator().MinAll(old_ndof);
		
      if(old_ndof == -1) 
        mndof = 1;
      else
        mndof = ndof;
      // 		KRATOS_WATCH(mndof);
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "trilinos ML solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:
    Teuchos::ParameterList mAztecParameterList;
    Teuchos::ParameterList mMLParameterList;
    SparseMatrixPointerType mpA;
    MLPreconditionerPointerType mpMLPrec;
    ScalingType mScalingType;
    bool mMLPrecIsInitialized;
    bool mReformPrecAtEachStep;
    double mtol;
    int mmax_iter;
    int mndof;

    /**
     * Assignment operator.
     */
    MultiLevelSolver& operator=(const MultiLevelSolver& Other);

    /**
     * Copy constructor.
     */
    MultiLevelSolver(const MultiLevelSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, MultiLevelSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MultiLevelSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_MULTILEVEL_SOLVER_H_INCLUDED  defined 


