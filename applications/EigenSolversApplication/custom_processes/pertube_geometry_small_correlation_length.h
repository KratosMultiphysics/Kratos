// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//

#if !defined(KRATOS_PERTUBE_GEOMETRY_SMALL_CORRELATION_LENGTH_PROCESS)
#define KRATOS_PERTUBE_GEOMETRY_SMALL_CORRELATION_LENGTH_PROCESS

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "custom_solvers/eigensystem_solver.h"
#include "custom_solvers/eigen_direct_solver.h"
#include "custom_utilities/omp_node_search.h"


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

/**
 * @class PertubeGeometrySmallCorrelationLength
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the center of gravity
 * @details It takes into account all elements in the ModelPart
 *
 * @author Manuel Messmer
 */
class KRATOS_API(EIGEN_SOLVERS_APPLICATION) PertubeGeometrySmallCorrelationLength
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PertubeGeometrySmallCorrelationLength
    KRATOS_CLASS_POINTER_DEFINITION(PertubeGeometrySmallCorrelationLength);

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TSparseSpaceType::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixPointerType DenseMatrixPointerType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef Kratos::intrusive_ptr<Kratos::Node<3UL> >                     NodesType;
    



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PertubeGeometrySmallCorrelationLength(
        ModelPart& rThisModelPart
        ):mrThisModelPart(rThisModelPart)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~PertubeGeometrySmallCorrelationLength() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int CreateEigenvectors(double minDistance, double correlationLength);

    void AssembleEigenvectors( const std::vector<double>& variables, double maxDisplacement );
   
    double Kernel( double x, double sigma );

    double CorrelationFunction( ModelPart::NodeIterator itNode1, NodesType itNode2, double CorrelationLenth);
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
    std::string Info() const override
    {
        return "PertubeGeometrySmallCorrelationLength";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PertubeGeometrySmallCorrelationLength";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPart& mrThisModelPart;              // The main model part

    DenseMatrixType mDisplacement;

    EigensystemSolver<>* mpEigenSolver;

    OMP_NodeSearch* searcher;

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
    PertubeGeometrySmallCorrelationLength& operator=(PertubeGeometrySmallCorrelationLength const& rOther) = delete;

    /// Copy constructor.
    PertubeGeometrySmallCorrelationLength(PertubeGeometrySmallCorrelationLength const& rOther) = delete;


    ///@}

}; // Class PertubeGeometrySmallCorrelationLength

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   PertubeGeometrySmallCorrelationLength& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const PertubeGeometrySmallCorrelationLength& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_PERTUBE_GEOMETRY_PROCESS defined */


       