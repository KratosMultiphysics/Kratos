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
#include "utilities/auxiliar_model_part_utilities.h"
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
    const Variable<double> &rVariable)
{
    AuxiliarModelPartUtilities(rModelPart).GetScalarData<Vector>(rVariable, Globals::DataLocation::NodeHistorical, rVector);
}

void OptimizationUtilities::AssembleVector( ModelPart& rModelPart,
    Vector& rVector,
    const Variable<array_3d> &rVariable)
{
    AuxiliarModelPartUtilities(rModelPart).GetVectorData<Vector>(rVariable, Globals::DataLocation::NodeHistorical, rVector);
}

void OptimizationUtilities::AssignVectorToVariable(ModelPart& rModelPart,
    const Vector& rVector,
    const Variable<double> &rVariable)
{
    AuxiliarModelPartUtilities(rModelPart).SetScalarData<Vector>(rVariable, Globals::DataLocation::NodeHistorical, rVector);
}

void OptimizationUtilities::AssignVectorToVariable(ModelPart& rModelPart,
    const Vector& rVector,
    const Variable<array_3d> &rVariable)
{
    AuxiliarModelPartUtilities(rModelPart).SetVectorData<Vector>(rVariable, Globals::DataLocation::NodeHistorical, rVector);
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

void OptimizationUtilities::AssembleMatrix(ModelPart& rModelPart,
    Matrix& rMatrix,
    const std::vector<Variable<double>*>& rVariables
)
{
    if ((rMatrix.size1() != rModelPart.NumberOfNodes() || rMatrix.size2() !=  rVariables.size())){
        rMatrix.resize(rModelPart.NumberOfNodes(), rVariables.size());
    }

    int i=0;
    for (auto & node_i : rModelPart.Nodes())
    {
        int j=0;
        for (Variable<double>* p_variable_j : rVariables)
        {
            const Variable<double>& r_variable_j = *p_variable_j;
            double& variable = node_i.FastGetSolutionStepValue(r_variable_j);
            rMatrix(i, j) = variable;
            ++j;
        }
        ++i;
    }
}

void OptimizationUtilities::AssembleMatrix(ModelPart& rModelPart,
    Matrix& rMatrix,
    const std::vector<Vector*>& rGradientVectors
)
{
    const int number_of_design_variables = rGradientVectors[0]->size();

    if ((rMatrix.size1() != number_of_design_variables || rMatrix.size2() !=  rGradientVectors.size())){
        rMatrix.resize(number_of_design_variables, rGradientVectors.size());
    }

    for (int i = 0; i < number_of_design_variables; ++i)
    {
        int j=0;
        for (Vector* p_variable_j : rGradientVectors)
        {
            rMatrix(i, j) = (*p_variable_j)[i];
            ++j;
        }
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

}  // namespace Kratos.
