// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "optimization_utilities.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "shape_optimization_application.h"
#include "linear_solvers/linear_solver.h"

// ==============================================================================

namespace Kratos
{

void OptimizationUtilities::ComputeControlPointUpdate(ModelPart& rModelPart, const double StepSize, const bool Normalize)
{
    KRATOS_TRY;

    // Normalize if specified
    if(Normalize)
    {
        const double max_norm_search_dir = ComputeMaxNormOfNodalVariable(rModelPart, SEARCH_DIRECTION);
        if(max_norm_search_dir>1e-10)
            for (auto & node_i : rModelPart.Nodes())
            {
                array_3d& search_dir = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);
                search_dir/=max_norm_search_dir;
            }
        else
            KRATOS_WARNING("ShapeOpt::ComputeControlPointUpdate") << "Normalization of search direction by max norm activated but max norm is < 1e-10. Hence normalization is omitted!" << std::endl;
    }

    // Compute update
    for (auto & node_i : rModelPart.Nodes())
        noalias(node_i.FastGetSolutionStepValue(CONTROL_POINT_UPDATE)) = StepSize * node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);

    KRATOS_CATCH("");
}

void OptimizationUtilities::AddFirstVariableToSecondVariable( ModelPart& rModelPart, const Variable<array_3d> &rFirstVariable, const Variable<array_3d> &rSecondVariable )
{
    for (auto & node_i : rModelPart.Nodes())
        noalias(node_i.FastGetSolutionStepValue(rSecondVariable)) += node_i.FastGetSolutionStepValue(rFirstVariable);
}

double OptimizationUtilities::ComputeL2NormOfNodalVariable( ModelPart& rModelPart, const Variable<array_3d> &rVariable)
{
    double l2_norm = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
        l2_norm += inner_prod(variable_vector,variable_vector);
    }
    return std::sqrt(l2_norm);
}

double OptimizationUtilities::ComputeL2NormOfNodalVariable(ModelPart& rModelPart, const Variable<double> &rVariable)
{
    double l2_norm = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        double &value = node_i.FastGetSolutionStepValue(rVariable);
        l2_norm += value*value;
    }
    return std::sqrt(l2_norm);
}

double OptimizationUtilities::ComputeMaxNormOfNodalVariable(ModelPart& rModelPart, const Variable<array_3d> &rVariable)
{
    double max_norm = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
        double squared_value = inner_prod(variable_vector,variable_vector);

        max_norm = std::max(squared_value,max_norm);
    }
    return std::sqrt(max_norm);
}

double OptimizationUtilities::ComputeMaxNormOfNodalVariable(ModelPart& rModelPart, const Variable<double> &rVariable)
{
    double max_norm = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        double &value = node_i.FastGetSolutionStepValue(rVariable);
        double squared_value = value*value;

        max_norm = std::max(squared_value,max_norm);
    }
    return std::sqrt(max_norm);
}

void OptimizationUtilities::ComputeSearchDirectionSteepestDescent(ModelPart& rModelPart)
{
    KRATOS_TRY;

    // Some output for information
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "No constraints given or active. The negative objective gradient is chosen as search direction..." << std::endl;

    // search direction is negative of filtered gradient
    for (auto & node_i : rModelPart.Nodes())
    {
        node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) = -1.0 * node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
    }

    KRATOS_CATCH("");
}

void OptimizationUtilities::ComputeProjectedSearchDirection(ModelPart& rModelPart)
{
    KRATOS_TRY;

    // Some output for information
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Constraint is active. Modified search direction on the constraint hyperplane is computed..." << std::endl;

    // Compute norm of constraint gradient
    double norm_2_dCds_i = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
        norm_2_dCds_i += inner_prod(dCds_i,dCds_i);
    }
    norm_2_dCds_i = std::sqrt(norm_2_dCds_i);

    // Avoid division by zero
    if(std::abs(norm_2_dCds_i)<1e-12)
        norm_2_dCds_i = 1.0;

    // Compute dot product of objective gradient and normalized constraint gradient
    double dot_dFds_dCds = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d dFds_i = node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
        array_3d dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
        dot_dFds_dCds += inner_prod(dFds_i,(dCds_i / norm_2_dCds_i));
    }

    // Compute and assign projected search direction
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& dFds_i = node_i.FastGetSolutionStepValue(DF1DX_MAPPED);
        array_3d& dCds_i = node_i.FastGetSolutionStepValue(DC1DX_MAPPED);

        array_3d projection_term = dot_dFds_dCds * (dCds_i / norm_2_dCds_i);

        node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) = -1 * (dFds_i - projection_term);
    }

    KRATOS_CATCH("");
}

double OptimizationUtilities::CorrectProjectedSearchDirection(ModelPart& rModelPart, const double PrevConstraintValue, const double ConstraintValue, const double CorrectionScaling, const bool IsAdaptive )
{
    // Check correction necessary
    if(ConstraintValue==0)
        return CorrectionScaling;

    // Perform correction
    double correction_scaling = CorrectionScaling;
    double correction_factor = ComputeCorrectionFactor(rModelPart, PrevConstraintValue, ConstraintValue, correction_scaling, IsAdaptive);
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d correction_term = correction_factor * ConstraintValue * node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
        node_i.FastGetSolutionStepValue(SEARCH_DIRECTION) -= correction_term;
    }

    return correction_scaling;
}

double OptimizationUtilities::ComputeCorrectionFactor(ModelPart& rModelPart, const double PrevConstraintValue, const double ConstraintValue, double& CorrectionScaling, const bool IsAdaptive)
{
    double norm_correction_term = 0.0;
    double norm_search_direction = 0.0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d correction_term = ConstraintValue * node_i.FastGetSolutionStepValue(DC1DX_MAPPED);
        norm_correction_term += inner_prod(correction_term,correction_term);

        array_3d ds = node_i.FastGetSolutionStepValue(SEARCH_DIRECTION);
        norm_search_direction += inner_prod(ds,ds);
    }
    norm_correction_term = std::sqrt(norm_correction_term);
    norm_search_direction = std::sqrt(norm_search_direction);

    if(IsAdaptive)
    {
        // Adapt constraint scaling

        // Three cases need to be covered
        // 1) In case we have two subsequently decreasing constraint values --> correction is fine --> leave current correction scaling
        // 2) In case the correction jumps over the constraint (change of sign) --> correction was too big --> reduce
        if(ConstraintValue*PrevConstraintValue<0.0)
        {
            CorrectionScaling *= 0.5;
            KRATOS_INFO("ShapeOpt") << "Correction scaling needs to decrease...." << std::endl;
        }
        // 3) In case we have subsequently increasing constraint value --> correction was too low --> increase
        if(std::abs(ConstraintValue)>std::abs(PrevConstraintValue) && ConstraintValue*PrevConstraintValue>0)
        {
            KRATOS_INFO("ShapeOpt") << "Correction scaling needs to increase...." << std::endl;
            CorrectionScaling = std::min(CorrectionScaling*2,1.0);
        }
    }

    return CorrectionScaling * norm_search_direction / norm_correction_term;
}

void OptimizationUtilities::AssembleVector( ModelPart& rModelPart,
    Vector& rVector,
    const Variable<array_3d> &rVariable)
{
    if (rVector.size() != rModelPart.NumberOfNodes()*3){
        rVector.resize(rModelPart.NumberOfNodes()*3);
    }

    int i=0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
        rVector[i*3+0] = variable_vector[0];
        rVector[i*3+1] = variable_vector[1];
        rVector[i*3+2] = variable_vector[2];
        ++i;
    }
}

void OptimizationUtilities::AssignVectorToVariable(ModelPart& rModelPart,
    const Vector& rVector,
    const Variable<array_3d> &rVariable)
{
    KRATOS_ERROR_IF(rVector.size() != rModelPart.NumberOfNodes()*3)
        << "AssignVectorToVariable: Vector size does not mach number of Nodes!" << std::endl;

    int i=0;
    for (auto & node_i : rModelPart.Nodes())
    {
        array_3d& variable_vector = node_i.FastGetSolutionStepValue(rVariable);
        variable_vector[0] = rVector[i*3+0];
        variable_vector[1] = rVector[i*3+1];
        variable_vector[2] = rVector[i*3+2];
        ++i;
    }
}

void OptimizationUtilities::AssembleMatrix(ModelPart& rModelPart,
    Matrix& rMatrix,
    const std::vector<Variable<array_3d>*>& rVariables
)
{
    if ((rMatrix.size1() != rModelPart.NumberOfNodes()*3 || rMatrix.size2() !=  rVariables.size())){
        rMatrix.resize(rModelPart.NumberOfNodes()*3, rVariables.size());
    }

    int i=0;
    for (auto & node_i : rModelPart.Nodes())
    {
        int j=0;
        for (Variable<array_3d>* p_variable_j : rVariables)
        {
            const Variable<array_3d>& r_variable_j = *p_variable_j;
            array_3d& variable_vector = node_i.FastGetSolutionStepValue(r_variable_j);
            rMatrix(i*3+0, j) = variable_vector[0];
            rMatrix(i*3+1, j) = variable_vector[1];
            rMatrix(i*3+2, j) = variable_vector[2];
            ++j;
        }
        ++i;
    }
}

void OptimizationUtilities::CalculateProjectedSearchDirectionAndCorrection(
    Vector& rObjectiveGradient,
    Matrix& rConstraintGradients,
    Vector& rConstraintValues,
    LinearSolver<DenseSpace, DenseSpace>& rSolver,
    Vector& rProjectedSearchDirection,
    Vector& rRestoration
    )
{
    // local variable naming according to https://msulaiman.org/onewebmedia/GradProj_2.pdf
    Vector& nabla_f = rObjectiveGradient;
    Matrix& N = rConstraintGradients;
    Vector& g_a = rConstraintValues;
    Vector& s = rProjectedSearchDirection;
    Vector& c = rRestoration;

    Matrix NTN = prod(trans(N), N);
    Matrix I = IdentityMatrix(N.size2());
    Matrix NTN_inv(NTN.size1(), NTN.size2());

    rSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

    s = - (nabla_f - prod(N, Vector(prod(NTN_inv, Vector(prod(trans(N), nabla_f))))));

    c = - prod(N, Vector(prod(NTN_inv, g_a)));
}

void OptimizationUtilities::ConstructActiveSubspace(
    Vector& rObjectiveGradient,
    Vector& rObjectiveGradientSubSpace,
    LinearSolver<TDenseSpaceTypeIn, TDenseSpaceTypeIn>& rSolver,
    Matrix& rEigenVectors,
    Vector& rEigenValues,
    double& rShapeFraction)
{
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace starts" << std::endl;
    long int node_number = rObjectiveGradient.size() / 3;
    TDenseSpaceTypeIn::MatrixType C(node_number, node_number, 0.0);
    // C.resize(node_number, node_number);

    for (int i = 0; i < node_number; ++i) {
        for (int j = 0; j < node_number; ++j) {
            C(i, j) = 0.0;
            for (int dim = 0; dim < 3; ++dim) {
                C(i, j) += rObjectiveGradient(3*i+dim) * rObjectiveGradient(3*j+dim);
            }
        }
    }
    // TSparseSpaceTypeIn::MatrixType C = outer_prod(rObjectiveGradient, rObjectiveGradient);
    TDenseSpaceTypeIn::MatrixType I = IdentityMatrix(C.size2());

    KRATOS_INFO("ShapeOpt") << "Covariance matrix:\n" << std::endl;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            KRATOS_INFO("ShapeOpt") << "C(" << i << "," << j << "): " << C(i,j) << std::endl;
        }
    }
    Vector eigenvalues;
    // Matrix eigenvectors;
    rSolver.Solve(C, I, rEigenValues, rEigenVectors);
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rEigenValues size: " << rEigenValues.size() << std::endl;
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rObjectiveGradient size: " << rObjectiveGradient.size() << std::endl;

    // rObjectiveGradientSubSpace = prod(trans(rEigenVectors), rObjectiveGradient);




    // long int eigen_value_number = rEigenVectors.size1();
    // double shape_fraction = 0.15;
    long int eigen_value_number = rShapeFraction*node_number;
    // long int eigen_value_number = 3;
    Matrix eigen_vectors_3(3*eigen_value_number, 3*node_number, 0.0);
    for (int i = 0; i < eigen_value_number; ++i) {
        for (int j = 0; j < node_number; ++j) {
            for (int dim = 0; dim < 3; ++dim) {
                eigen_vectors_3(3*i+dim, 3*j+dim) = rEigenVectors(i, j);
            }
        }
    }
    // KRATOS_INFO("ShapeOpt") << "rEigenVectors:\n" << std::endl;
    // for (int i = 0; i < 3; ++i) {
    //     for (int j = 0; j < node_number; ++j) {
    //         KRATOS_INFO("ShapeOpt") << "rEigenVectors(" << i << "," << j << "): " << rEigenVectors(i,j) << std::endl;
    //     }
    // }

    rObjectiveGradientSubSpace = prod(eigen_vectors_3, rObjectiveGradient);
    rEigenVectors = eigen_vectors_3;
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rEigenVectors size: " << rEigenVectors.size1() << " by " << rEigenVectors.size2() << std::endl;
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rObjectiveGradientSubSpace size: " << rObjectiveGradientSubSpace.size() << std::endl;
    // KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rObjectiveGradientSubSpace: " << rObjectiveGradientSubSpace << std::endl;
    // KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace eigenvalues: " << eigenvalues << std::endl;
    // KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace eigenvalues: " << eigenvectors << std::endl;
}

void OptimizationUtilities::ActiveSubspaceInverseMap(
    Vector& rSubspaceInputVector,
    Vector& rFullspaceOutputVector,
    LinearSolver<DenseSpace, DenseSpace>& rSolver,
    Matrix& rEigenVectors)
{
    Matrix eigen_vectors_t = trans(rEigenVectors);
    // for (unsigned i = 0; i < rFullspaceOutputVector.size (); ++ i)
    //     rFullspaceOutputVector (i) = 0;

    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace eigen_vectors_t size: " << eigen_vectors_t.size1() << " by " << eigen_vectors_t.size2() << std::endl;
    KRATOS_INFO("ShapeOpt") << "ConstructActiveSubspace rSubspaceInputVector size: " << rSubspaceInputVector.size() << std::endl;
    Vector output_vector;
    // rSolver.Solve(eigen_vectors_t, rFullspaceOutputVector, rSubspaceInputVector);

    // Method 1: Using transpose
    rFullspaceOutputVector = prod(eigen_vectors_t, rSubspaceInputVector);

    // Method 2: Ordinary Least Squares
    // Matrix A = prod(eigen_vectors_t, rEigenVectors);
    // Vector b = prod(eigen_vectors_t, rSubspaceInputVector);

    // rSolver.Solve(A, rFullspaceOutputVector, b);


    // long int node_number = rFullspaceOutputVector.size() / 3;
    // for (int i = 0; i < node_number; ++i) {
    //     for (int dim = 0; dim < 3; ++dim) {
    //         rFullspaceOutputVector(3*i+dim) = output_vector(i);
    //     }
    // }
}

}  // namespace Kratos.
