// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:
//
// ==============================================================================

#ifndef HEAT_METHOD_UTILITIES_H
#define HEAT_METHOD_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

// ==============================================================================

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

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) HeatMethodUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // For better reading
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef ModelPart::ElementType::GeometryType GeometryType;
    typedef std::size_t SizeType;
    typedef boost::numeric::ublas::vector<unsigned int> Vector_int;
    typedef zero_vector<unsigned int> ZeroVector_int;
    typedef matrix<unsigned int> Matrix_int;
    typedef zero_matrix<unsigned int> ZeroMatrix_int;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpace;
    typedef UblasSpace<double, SparseMatrix, Vector> SparseSpace;
    using NodeType = Node;

    /// Pointer definition of HeatMethodUtilities
    KRATOS_CLASS_POINTER_DEFINITION(HeatMethodUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HeatMethodUtilities( ModelPart& modelPart, LinearSolver<SparseSpaceType, LocalSpaceType>& Solver )
        : mrModelPart( modelPart ), mrSolver ( Solver )
    {
    }

    /// Destructor.
    virtual ~HeatMethodUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // TODO: structure algorithm

    void ComputeGeodesicDistance();

    // etc....

    // --------------------------------------------------------------------------

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
        return "HeatMethodUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "HeatMethodUtilities";
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

    ModelPart& mrModelPart;
    LinearSolver<SparseSpaceType, LocalSpaceType>& mrSolver;

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

    /// Assignment operator.
//      HeatMethodUtilities& operator=(HeatMethodUtilities const& rOther);

    /// Copy constructor.
//      HeatMethodUtilities(HeatMethodUtilities const& rOther);


    ///@}

}; // Class HeatMethodUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // HEAT_METHOD_UTILITIES_H