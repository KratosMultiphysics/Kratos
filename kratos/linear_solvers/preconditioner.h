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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_PRECONDITIONER_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "includes/model_part.h"


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

/// Preconditioner class.
/** Base class for preconditioners for linesr system solvers which defining
    standard interface for all the preconditioners derived from it.

    Considering a linear solver FooSolver with a method FooSolver::Solve.
    A typical code using this type of Preconditioners would be:

    \begin{verbatim}
     FooSolver::Solve(A,b,x,preconditioner)
     {
        preconditioner.Initialize(A,x,b);
    ...
    ...
    while(...) // Start iteration.
    {
            preconditioner.ApplyLeft(x);
        mult(a,x)
        preconditioner.ApplyRight(x)
    } // End iteration

    preconditioner.Finalize(A,x,b);
     }
     \end{verbatim}
*/
template<class TSparseSpaceType, class TDenseSpaceType>
class Preconditioner
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Preconditioner
    KRATOS_CLASS_POINTER_DEFINITION(Preconditioner);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Preconditioner() {}

    /// Copy constructor.
    Preconditioner(const Preconditioner& Other) {}

    /// Destructor.
    virtual ~Preconditioner() {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Preconditioner& operator=(const Preconditioner& Other)
    {
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** Preconditioner Initialize
    Initialize preconditioner for linear system rA*rX=rB
    @param rA  system matrix.
    @param rX Unknows vector
    @param rB Right side linear system of equations.
    */
    virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) {}

    virtual void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
        VectorType x(TDenseSpaceType::Size1(rX));
        VectorType b(TDenseSpaceType::Size1(rB));

        TDenseSpaceType::GetColumn(0,rX, x);
        TDenseSpaceType::GetColumn(0,rB, b);

        Initialize(rA, x, b);
    }

    /** This function is designed to be called every time the coefficients change in the system
    * that is, normally at the beginning of each solve.
    * For example if we are implementing a direct solver, this is the place to do the factorization
    * so that then the backward substitution can be performed effectively more than once
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    virtual void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    virtual void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    virtual void Clear()
    {
    }

    /** Some preconditioners may require a minimum degree of knowledge of the structure of the matrix. To make an example
    * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
    * another example is the automatic prescription of rotation null-space for smoothed-aggregation preconditioners
    * which require knowledge on the spatial position of the nodes associated to a given dof.
    * This function tells if the solver requires such data
    */
    virtual bool AdditionalPhysicalDataIsNeeded()
    {
        return false;
    }

    /** Some preconditioners may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation preconditioners
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    virtual void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    )
    {}

    virtual void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        VectorType z = rX;
        ApplyRight(z);
        TSparseSpaceType::Mult(rA,z, rY);
        ApplyLeft(rY);
    }

    virtual void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        VectorType z = rX;
        ApplyTransposeLeft(z);
        TSparseSpaceType::TransposeMult(rA,z, rY);
        ApplyTransposeRight(rY);
    }

    virtual VectorType& ApplyLeft(VectorType& rX)
    {
        return rX;
    }

    virtual VectorType& ApplyRight(VectorType& rX)
    {
        return rX;
    }

    /** Preconditioner transpose solver.
    Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
    @param rX  Unknows of preconditioner suystem
    */
    virtual VectorType& ApplyTransposeLeft(VectorType& rX)
    {
// 	KRATOS_ERROR(std::logic_error,
// 		     " virtual TVectorType& ApplyTransposeLeft(TVectorType& rX)",
// 		     "This preconditioner dosn't have ApplyTransposeLeft defined.", "");
        return rX;
    }

    virtual VectorType& ApplyTransposeRight(VectorType& rX)
    {
// 	KRATOS_ERROR(std::logic_error,
// 		     " virtual TVectorType& ApplyTransposeRight(TVectorType& rX)",
// 		     "This preconditioner dosn't have ApplyTransposeRight defined.", "");
        return rX;
    }

    virtual VectorType& ApplyInverseRight(VectorType& rX)
    {
        return rX;
    }

    /* The method Finalize is used to recover the value of rX.
       In principle, it is enough to multiply by the right preconditioner.
    See the diagoinal preconditioner for a nontrivial example. */
    virtual VectorType& Finalize(VectorType& rX)
    {
        return ApplyRight(rX);
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
        return "Preconditioner";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Preconditioner";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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


    ///@}
    ///@name Private Operators
    ///@{


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

}; // Class Preconditioner

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::istream& operator >> (std::istream& IStream,
                                  Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PRECONDITIONER_H_INCLUDED  defined 


