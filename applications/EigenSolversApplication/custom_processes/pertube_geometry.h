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

#if !defined(KRATOS_PERTUBE_GEOMETRY_PROCESS)
#define KRATOS_PERTUBE_GEOMETRY_PROCESS

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
#include "utilities/mortar_utilities.h"


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
 * @class PertubeGeometryProcess
 *
 * @ingroup StructuralMechanicsApplication
 *
 * @brief This method computes the center of gravity
 * @details It takes into account all elements in the ModelPart
 *
 * @author Manuel Messmer
 */
class KRATOS_API(EIGEN_SOLVERS_APPLICATION) PertubeGeometryProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PertubeGeometryProcess
    KRATOS_CLASS_POINTER_DEFINITION(PertubeGeometryProcess);

    typedef TUblasSparseSpace<double> TSparseSpaceType;
    typedef TUblasDenseSpace<double> TDenseSpaceType;

    typedef typename TSparseSpaceType::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;



    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PertubeGeometryProcess( ModelPart& rInitialModelPart, double MaximalDisplacement) : 
        mrInitialModelPart(rInitialModelPart),
        mMaximalDisplacement(MaximalDisplacement)
    {
        KRATOS_TRY
        MortarUtilities::ComputeNodesMeanNormalModelPart( mrInitialModelPart, false );
        // Remove this
        //###############################################################################
        CorrelationMatrix_check = Eigen::MatrixXd::Zero(mrInitialModelPart.NumberOfNodes(), mrInitialModelPart.NumberOfNodes());
        CorrelationMatrix_check_orig = Eigen::MatrixXd::Zero(mrInitialModelPart.NumberOfNodes(), mrInitialModelPart.NumberOfNodes());
        //#################################################################################
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~PertubeGeometryProcess() override
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

    int CreateEigenvectors(ModelPart& rThisModelPart, double minDistance, double correlationLength, double truncationTolerance);

    void AssembleEigenvectors(ModelPart& rThisModelPart, const std::vector<double>& variables, double correlationLength );
   
    double Kernel( double x, double sigma );

    double CorrelationFunction( ModelPart::NodeIterator itNode1, ModelPart::NodeIterator itNode2, double CorrelationLenth);
    //Remove this #################
    void Average(int number);
    //##############################

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
        return "PertubeGeometryProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PertubeGeometryProcess";
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

    ModelPart& mrInitialModelPart;

    Eigen::MatrixXd Displacement;

    OMP_NodeSearch* searcher;

    double mMaximalDisplacement;

    // Remove this again ################################
    Eigen::MatrixXd CorrelationMatrix_check;
    Eigen::MatrixXd CorrelationMatrix_check_orig;
    //###################################################

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
    PertubeGeometryProcess& operator=(PertubeGeometryProcess const& rOther) = delete;

    /// Copy constructor.
    PertubeGeometryProcess(PertubeGeometryProcess const& rOther) = delete;


    ///@}

}; // Class PertubeGeometryProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   PertubeGeometryProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const PertubeGeometryProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_PERTUBE_GEOMETRY_PROCESS defined */


       