//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Guillermo Casas (gcasas@cimne.upc.edu)
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/geometry_utilities.h"
// Application includes
#include "../../swimming_DEM_application.h"
#include "derivative_recovery_utility.h"
#include "standard_recovery_utility.h"

namespace Kratos
{

using ScalarVariableType = StandardRecoveryUtility::ScalarVariableType;
using ComponentVariableType = StandardRecoveryUtility::ComponentVariableType;
using VectorVariableType = StandardRecoveryUtility::VectorVariableType;

StandardRecoveryUtility::StandardRecoveryUtility(
    ModelPart& rModelPart,
    Parameters rParameters,
    RecoveryVariablesContainer& rVariablesContainer)
    : DerivativeRecoveryUtility(rModelPart, rParameters, rVariablesContainer)
{
    this->CheckDefaultsAndSettings(rParameters);
    mStoreFullGradient = rParameters["store_full_gradient_option"].GetBool();
}

StandardRecoveryUtility::StandardRecoveryUtility(
    Model& rModel,
    Parameters rParameters,
    RecoveryVariablesContainer& rVariablesContainer)
    : DerivativeRecoveryUtility(rModel, rParameters, rVariablesContainer)
{
    this->CheckDefaultsAndSettings(rParameters);
    mStoreFullGradient = rParameters["store_full_gradient_option"].GetBool();
}

void StandardRecoveryUtility::CheckDefaultsAndSettings(Parameters rParameters)
{
    Parameters default_parameters( R"(
    {
        "model_part_name" : "FluidModelPart",
        "recoverer_name" : "StandardRecoveryUtility",
        "store_full_gradient_option" : true
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);
}

void StandardRecoveryUtility::Initialize()
{

}

void StandardRecoveryUtility::AddPartialTimeDerivative(const VariableData& rVariable, const VariableData& rTimeDerivativeVariable)
{
    this->AddPartialTimeDerivative<VariableData, VariableData>(rVariable, rTimeDerivativeVariable);
}

void StandardRecoveryUtility::CalculateGradient(const VariableData& rVariable, const VariableData& rGradientVariable)
{
    this->CalculateScalarGradient<ScalarVariableType>(rScalarVariable, rGradientVariable);
}

void StandardRecoveryUtility::CalculateGradient(const ComponentVariableType& rScalarComponent, const VectorVariableType& rGradientVariable)
{
    this->CalculateScalarGradient<ComponentVariableType>(rScalarComponent, rGradientVariable);
}
// void StandardRecoveryUtility::CalculateGradient(const VectorVariableType& rVectorVariable,
//                                                 const VectorVariableType& rComponent0GradientVariable,
//                                                 const VectorVariableType& rComponent1GradientVariable,
//                                                 const VectorVariableType& rComponent2GradientVariable)
// {

// }

void StandardRecoveryUtility::CalculateDivergence(const VectorVariableType& rVectorVariable, const ScalarVariableType& rDivergenceVariable)
{

}
void StandardRecoveryUtility::CalculateLaplacian(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rLaplacianVariable)
{
    this->CalculateScalarLaplacian<ScalarVariableType>(rScalarVariable, rLaplacianVariable);
}
void StandardRecoveryUtility::CalculateLaplacian(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rLaplacianVariable)
{
    this->CalculateScalarLaplacian<ComponentVariableType>(rScalarComponent, rLaplacianVariable);
}
void StandardRecoveryUtility::CalculateLaplacian(const VectorVariableType& rVectorComponent, const VectorVariableType& rLaplacianVariable)
{

}
void StandardRecoveryUtility::CalculateMaterialDerivative(const ScalarVariableType& rScalarVariable, const ScalarVariableType& rMaterialDerivativeVariable)
{
    this->CalculateScalarMaterialDerivative<ScalarVariableType>(rScalarVariable, rMaterialDerivativeVariable);
}
void StandardRecoveryUtility::CalculateMaterialDerivative(const ComponentVariableType& rScalarComponent, const ScalarVariableType& rMaterialDerivativeVariable)
{
    this->CalculateScalarMaterialDerivative<ComponentVariableType>(rScalarComponent, rMaterialDerivativeVariable);
}
void StandardRecoveryUtility::CalculateMaterialDerivative(const VectorVariableType& rVectorVariable,
                                                          const VectorVariableType& rMaterialDerivativeVariable)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(rMaterialDerivativeVariable)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }

    std::vector<array_1d <double, 3> > convective_contributions_to_the_derivative;
    convective_contributions_to_the_derivative.resize(entry);

    array_1d <double, 3> grad = ZeroVector(3);
    const int TDim = 3;
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;

    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value

        // for each element, constructing the gradient contribution (to its nodes) of the component v_j
        // and storing it in rMaterialDerivativeVariable

        for (auto ielem = mrModelPart.ElementsBegin(); ielem != mrModelPart.ElementsEnd(); ++ielem){
            // computing the shape function derivatives
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(rVectorVariable)[j];
            }

            array_1d <double, 3> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }

            double nodal_area = Volume / static_cast<double>(TDim + 1);
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(rMaterialDerivativeVariable) += grad; // we use rMaterialDerivativeVariable to store the gradient of one component at a time
            }
        }

        // normalizing the contributions to the gradient and getting the j-component of the material derivative

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            array_1d <double, 3>& stored_gradient_of_component_j = inode->FastGetSolutionStepValue(rMaterialDerivativeVariable);
            stored_gradient_of_component_j /= inode->FastGetSolutionStepValue(NODAL_AREA);

            if (mStoreFullGradient){
                if (j == 0){
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }

                else if (j == 1){
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }

                else {
                    array_1d <double, 3>& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }
            }

            const array_1d <double, 3>& velocity = inode->FastGetSolutionStepValue(VELOCITY);
            convective_contributions_to_the_derivative[id_to_position[inode->Id()]][j] = SWIMMING_INNER_PRODUCT_3(velocity, stored_gradient_of_component_j);
            stored_gradient_of_component_j = ZeroVector(3);
        }

    }

    // Adding convective part

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        const array_1d <double, 3>& stored_convective_contribution = convective_contributions_to_the_derivative[id_to_position[inode->Id()]];
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(rMaterialDerivativeVariable);
        material_derivative = stored_convective_contribution;
    }

    // Adding Eulerian time derivative contribution

    AddPartialTimeDerivative(VELOCITY, MATERIAL_ACCELERATION);

    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by derivating nodal averages..." << std::endl;
}

void StandardRecoveryUtility::CalculateRotational(const VectorVariableType rVectorVariable, const VectorVariableType& rRotationalVariable)
{

}

/* Private functions ****************************************************/

template<class TVariable, class TDerivedVariable>
void StandardRecoveryUtility::AddPartialTimeDerivative(const TVariable& rVariable, const TDerivedVariable& rTimeDerivativeVariable)
{
    KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
}

template<class ScalarVariableType, class ScalarVariableType>
void StandardRecoveryUtility::AddPartialTimeDerivative(const ScalarVariableType& rVariable, const ScalarVariableType& rTimeDerivativeVariable)
{
    this->AddScalarPartialTimeDerivative<ScalarVariableType>(rVariable, rTimeDerivativeVariable);
}

template<class ComponentVariableType, class ScalarVariableType>
void StandardRecoveryUtility::AddPartialTimeDerivative(const ComponentVariableType& rVariable, const ScalarVariableType& rTimeDerivativeVariable)
{
    this->AddScalarPartialTimeDerivative<ComponentVariableType>(rVariable, rTimeDerivativeVariable);
}

template<class VectorVariableType, class VectorVariableType>
void StandardRecoveryUtility::AddPartialTimeDerivative(const VectorVariableType& rVariable, const VectorVariableType& rTimeDerivativeVariable)
{
    const double delta_time_inv = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(rTimeDerivativeVariable);
        const array_1d <double, 3> eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(rVariable)
                                                                               - inode->FastGetSolutionStepValue(rVariable, 1));
        noalias(material_derivative) += eulerian_rate_of_change;
    }
}

template<class TScalarVariable>
void StandardRecoveryUtility::AddScalarPartialTimeDerivative(const TScalarVariable& rVariable, const ScalarVariableType& rTimeDerivativeVariable)
{
    const double delta_time_inv = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        double& material_derivative = inode->FastGetSolutionStepValue(rTimeDerivativeVariable);
        const double eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(rVariable)
                                                                 - inode->FastGetSolutionStepValue(rVariable, 1));
        material_derivative += eulerian_rate_of_change;
    }
}

template<class TVariable>
void StandardRecoveryUtility::CalculateScalarGradient(const TVariable& rScalarVariable, const TVariable& rGradientVariable)
{
}

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarGradient(const TScalarVariable& rScalarVariable, const VectorVariableType& rGradientVariable)
{
}

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarLaplacian(const TScalarVariable& rScalarVariable, const ScalarVariableType& rLaplacianVariable)
{
}

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarMaterialDerivative(const TScalarVariable& rScalarVariable, const ScalarVariableType& rMaterialDerivativeVariable)
{
}

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<ScalarVariableType>(const ScalarVariableType&, const ScalarVariableType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<ComponentVariableType>(const ComponentVariableType&, const ScalarVariableType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<ScalarVariableType>(const ScalarVariableType&, const VectorVariableType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<ComponentVariableType>(const ComponentVariableType&, const VectorVariableType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<ScalarVariableType>(const ScalarVariableType&, const ScalarVariableType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<ComponentVariableType>(const ComponentVariableType&, const ScalarVariableType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<ScalarVariableType>(const ScalarVariableType&, const ScalarVariableType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<ComponentVariableType>(const ComponentVariableType&, const ScalarVariableType&);

}  // namespace Kratos.
