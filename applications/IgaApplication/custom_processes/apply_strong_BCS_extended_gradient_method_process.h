//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti

#if !defined(KRATOS_APPLY_STRONG_BCS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_STRONG_BCS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED

// System includes

// External includes
#include <chrono>

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "containers/array_1d.h" 
#include "utilities/function_parser_utility.h"
#include "iga_application_variables.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/mls_shape_functions_utility.h"

// Linear solvers
#include "linear_solvers/reorderer.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/amgcl_solver.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/* @class ApplyStrongBCSExtendedGradientMethodProcess
 * @ingroup IgaApplication
 **/
class KRATOS_API(IGA_APPLICATION) ApplyStrongBCSExtendedGradientMethodProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef typename GeometryType::Pointer GeometryPointerType;
    typedef array_1d<double, 3> CoordinatesArrayType;
    typedef typename GeometryType::IntegrationPointsArrayType IntegrationPointsArrayType;
    typedef typename GeometryType::GeometriesArrayType GeometriesArrayType;
    typedef Element BaseType;
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;

    /// Tests
    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
    
    // The direct solver
    using ReordererType = Reorderer<SparseSpaceType, LocalSpaceType>;
    using DirectSolverType = DirectSolver<SparseSpaceType, LocalSpaceType, ReordererType>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    //using AMGCLSolverType = AMGCLSolver<SparseSpaceType, LocalSpaceType, ReordererType>;
    using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SparseSpaceType, LocalSpaceType, ReordererType>;
    // using PreconditionerType = Preconditioner<SparseSpaceType, LocalSpaceType>;
    // using MixedULMLinearSolverType = MixedULMLinearSolver<SparseSpaceType, LocalSpaceType, PreconditionerType, ReordererType>;


    /// Pointer definition of ApplyStrongBCSExtendedGradientMethodProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyStrongBCSExtendedGradientMethodProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ApplyStrongBCSExtendedGradientMethodProcess(
        Model& rModel,
        Parameters ThisParameters);

    /// Destructor.
    ~ApplyStrongBCSExtendedGradientMethodProcess() = default;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitializeSolutionStep() override {
        Execute();
    };

    void ComputeLHS(CompressedMatrix& LHS);

    void ComputeLHSBoundaryContribution(CompressedMatrix& LHSBoundaryContribution);

    void ComputeLHSTrimmedElementsGradientContribution(CompressedMatrix& LHSTrimmedElementsGradientContribution);

    void ComputeRHSBoundaryContribution(Vector& RHSBoundaryContribution);

    void ComputeRHSTrimmedElementsGradientContribution(Vector& RHSTrimmedElementsGradientContribution);

    void VerifyAndModifyDiagonalLHS(CompressedMatrix& LHS);

    void DefineInterpolationPointsAndGradients(Matrix& InterpolationPoints, Vector& Solution, Vector& SolutionGradientX, Vector& SolutionGradientY,  std::vector<GeometryPointerType>& InterpolationPointsPointers);

    void InterpolateSolutionGradients(Vector& SolutionGradients, double& Solution, CoordinatesArrayType position, IndexType NumberOfClosestPoints);

    void FindTheNClosestInterpolationPoints(IndexType NumberOfClosestPoints, CoordinatesArrayType position);

    double CalculateKernelParameterMLS(CoordinatesArrayType quadrature_point_position);

    void ApplyStrongBoundaryConditions();

    void ComputeExactGradient(double x, double y, Vector& exact_solution_gradient);


    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "background_mesh_model_part_name": "please_specify_model_part_name",
            "skin_model_part_name": "please_specify_model_part_name",
            "variable_name": "please_specify_a_variable",
            "value"           : "please_specify_a_value",
            "initial_gradient_inside_trimmed_elements"           : 0.0,
            "interpolation_scheme": "please_specify_an_interpolation_scheme",
            "number_of_interpolation_points": 200,
            "MLS_polinomial_order": 1
        })" );

        return default_parameters;
    }

    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyStrongBCSExtendedGradientMethodProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyStrongBCSExtendedGradientMethodProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Iga functionalities
    ///@{

    Model* mpModel = nullptr;
    Parameters mParameters;
    SizeType mEchoLevel;
    IndexType mIterations = 0;
    ModelPart* mpSkinModelPart = nullptr;
    ModelPart* mpBackgroundMeshModelPart = nullptr;
    CompressedMatrix mLHS;
    Vector mRHS;
    Vector mRHSBoundaryContribution;
    Vector mRHSTrimmedElementsContribution;
    Vector mPhiDir;
    IndexType mNumberOfNodesBackgroundMesh; 
    std::string mUnknownVariableName;
    const Kratos::Variable<double>* mpUnknownVariable; 
    std::string mInterpolationScheme;
    Kratos::unique_ptr<GenericFunctionUtility> mpEvalFunction;
    IndexType mMLSPolinomialOrder = 0;
    IndexType mInterpolationPointsNumber = 0;
    double mInitialGradientInsideTrimmedElements = 0.0;

    // Member variables for the interpolation
    Matrix mInterpolationPointsCoordinates;
    Vector mSolutionGradientX;
    Vector mSolutionGradientY;
    Vector mSolution;
    std::vector<GeometryPointerType> mInterpolationPointsPointers;

    Matrix mClosestInterpolationPointsCoordinates;
    Vector mClosestPointsSolutionGradientX;
    Vector mClosestPointsSolutionGradientY;
    Vector mClosestPointsSolution;
    std::vector<GeometryPointerType> mClosestInterpolationPointsPointers;

    ///@}
    ///@name Iga functionalities
    ///@{

    ///@}
    ///@}

    ///@}
    ///@name Utility
    ///@{


    ///@}
    ///@name Input and output
    ///@{


   

    ///@}

}; // Class ApplyStrongBCSExtendedGradientMethodProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyStrongBCSExtendedGradientMethodProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyStrongBCSExtendedGradientMethodProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_APPLY_STRONG_BCS_EXTENDED_GRADIENT_METHOD_PROCESS_H_INCLUDED 
