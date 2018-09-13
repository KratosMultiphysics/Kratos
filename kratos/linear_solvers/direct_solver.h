//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


#if !defined(KRATOS_DIRECT_SOLVER_H_INCLUDED )
#define  KRATOS_DIRECT_SOLVER_H_INCLUDED



// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solver.h"


namespace Kratos
{



// Base class for all direct solvers in Kratos.
/* This class define the general interface for direct solvers in Kratos.
   direct solver is a template class with this parameter:
   - TSparseSpaceType which specify type
     of the unknowns, coefficients, sparse matrix, vector of
   unknowns, right hand side vector and their respective operators.
   - TDenseMatrixType which specify type of the
     matrices used as temporary matrices or multi solve unknowns and
   right hand sides and their operators.
   - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
*/
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DirectSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    /// Counted pointer of DirectSolver
    KRATOS_CLASS_POINTER_DEFINITION(DirectSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;


    /// Default constructor.
    DirectSolver() {}
    DirectSolver(Parameters settings) {}
    /// Destructor.
    ~DirectSolver() override {}

    /// Copy constructor.
    DirectSolver(const DirectSolver& Other) {}

    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Direct solver";
    }

    /// Print object's data.
    void  PrintData(std::ostream& rOStream) const override
    {
    }

private:

    /// Assignment operator.
    DirectSolver& operator=(const DirectSolver& Other);



}; // Class DirectSolver



/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream,
                                  DirectSolver<TSparseSpaceType, TDenseSpaceType,TReordererType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DirectSolver<TSparseSpaceType,TDenseSpaceType,TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_DIRECT_SOLVER_H_INCLUDED  defined 


