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

#include "swimming_dem_application_variables.h"
namespace Kratos
{

using Vector3 = StandardRecoveryUtility::Vector3;
using DoubleVarType = StandardRecoveryUtility::DoubleVarType;
using ComponentVarType = StandardRecoveryUtility::ComponentVarType;
using ArrayVarType = StandardRecoveryUtility::ArrayVarType;
using TensorComponentVarType = StandardRecoveryUtility::ArrayVarType;

StandardRecoveryUtility::StandardRecoveryUtility(
    ModelPart& rModelPart,
    Parameters rParameters): DerivativeRecoveryUtility(rModelPart, rParameters)
{
}

StandardRecoveryUtility::StandardRecoveryUtility(
    Model& rModel,
    Parameters rParameters): DerivativeRecoveryUtility(rModel, rParameters)
{
}

void StandardRecoveryUtility::CheckDefaultsAndSettings(Parameters rParameters)
{
}

void StandardRecoveryUtility::Initialize()
{

}

void StandardRecoveryUtility::AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable)
{
    //this->AddPartialTimeDerivative<rVariable.StaticObject(), rTimeDerivativeVariable.StaticObject()>(rVariable, rTimeDerivativeVariable);
}

void StandardRecoveryUtility::AddPartialTimeDerivative(const ArrayVarType& rVariable, const ArrayVarType& rTimeDerivativeVariable)
{
    const double delta_time_inv = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        Vector3& material_derivative = inode->FastGetSolutionStepValue(rTimeDerivativeVariable);
        const Vector3 eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(rVariable)
                                                                               - inode->FastGetSolutionStepValue(rVariable, 1));
        noalias(material_derivative) += eulerian_rate_of_change;
    }
}

// void StandardRecoveryUtility::CalculateGradient(const VariableData& rVariable, const VariableData& rGradientVariable)
// {
//     //this->CalculateScalarGradient<rVariable.StaticObject(), rGradientVariable.StaticObject()>(rScalarVariable, rGradientVariable);
// }

void StandardRecoveryUtility::CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable)
{
    this->CalculateScalarGradient<DoubleVarType>(rScalarVariable, rGradientVariable);
}

void StandardRecoveryUtility::CalculateGradient(const ArrayVarType& rVectorVariable, const Tensor3 rGradientVariable)
{
    for (auto comp_i : {"X", "Y", "Z"}){
        const auto& variable_i = dynamic_cast<ComponentVarType&>(KratosComponents<VariableData>::Get(rVectorVariable.Name() + "_" + comp_i));
        this->CalculateScalarGradient<ComponentVarType>(variable_i, SCALAR_GRADIENT);
        int j = 0;

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            const auto& component_gradient = inode->FastGetSolutionStepValue(SCALAR_GRADIENT);
            j = 0;
            for (auto comp_j : {"X", "Y", "Z"}){
                const auto& variable_ij = dynamic_cast<TensorComponentVarType&>(KratosComponents<VariableData>::Get(rVectorVariable.Name() + "_" + comp_i + comp_j));
                auto& vector_gradient_ij = inode->FastGetSolutionStepValue(variable_ij);
                vector_gradient_ij = component_gradient[j];
                ++j;
            }
        }
    }
}

// }

// void StandardRecoveryUtility::CalculateDivergence(const VariableData& rVectorVariable, const VariableData& rDivergenceVariable)
// {
// }
// void StandardRecoveryUtility::CalculateLaplacian(const VariableData& rVariable, const VariableData& rLaplacianVariable)
// {
//     // this->CalculateScalarLaplacian<rVariable.StaticObject(), rLaplacianVariable.StaticObject()>(rScalarVariable, rLaplacianVariable);
// }
// // void StandardRecoveryUtility::CalculateLaplacian(const ComponentVarType& rScalarComponent, const DoubleVarType& rLaplacianVariable)
// // {
// //     this->CalculateScalarLaplacian<rVariable.StaticObject(), rLaplacianVariable.StaticObject()>(rScalarComponent, rLaplacianVariable);
// // }
// // void StandardRecoveryUtility::CalculateLaplacian(const ArrayVarType& rVectorComponent, const ArrayVarType& rLaplacianVariable)
// // {

// // }
void StandardRecoveryUtility::CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable)
{
    this->CalculateScalarMaterialDerivative<double>(rVariable, rMaterialDerivativeVariable);
}
//TODO: generalize for any element
void StandardRecoveryUtility::CalculateMaterialDerivative(const ArrayVarType& rVectorVariable,
                                                          const ArrayVarType& rMaterialDerivativeVariable)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(rMaterialDerivativeVariable)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }

    std::vector<Vector3 > convective_contributions_to_the_derivative;
    convective_contributions_to_the_derivative.resize(entry);

    Vector3 grad = ZeroVector(3);
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

            Vector3 grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

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
            Vector3& stored_gradient_of_component_j = inode->FastGetSolutionStepValue(rMaterialDerivativeVariable);
            stored_gradient_of_component_j /= inode->FastGetSolutionStepValue(NODAL_AREA);

            if (mStoreFullGradient){
                if (j == 0){
                    Vector3& gradient = inode->FastGetSolutionStepValue(VELOCITY_X_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }

                else if (j == 1){
                    Vector3& gradient = inode->FastGetSolutionStepValue(VELOCITY_Y_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }

                else {
                    Vector3& gradient = inode->FastGetSolutionStepValue(VELOCITY_Z_GRADIENT);
                    noalias(gradient) = stored_gradient_of_component_j;
                }
            }

            const Vector3& velocity = inode->FastGetSolutionStepValue(VELOCITY);
            convective_contributions_to_the_derivative[id_to_position[inode->Id()]][j] = SWIMMING_INNER_PRODUCT_3(velocity, stored_gradient_of_component_j);
            stored_gradient_of_component_j = ZeroVector(3);
        }

    }

    // Adding convective part

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        const Vector3& stored_convective_contribution = convective_contributions_to_the_derivative[id_to_position[inode->Id()]];
        Vector3& material_derivative = inode->FastGetSolutionStepValue(rMaterialDerivativeVariable);
        material_derivative = stored_convective_contribution;
    }

    // Adding Eulerian time derivative contribution

    AddPartialTimeDerivative(VELOCITY, MATERIAL_ACCELERATION);

    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by derivating nodal averages..." << std::endl;
}

// void StandardRecoveryUtility::CalculateMaterialDerivative(const ComponentVarType& rScalarComponent, const DoubleVarType& rMaterialDerivativeVariable)
// {
//     this->CalculateScalarMaterialDerivative<ComponentVarType>(rScalarComponent, rMaterialDerivativeVariable);
// }

// void StandardRecoveryUtility::CalculateRotational(const VariableData rVectorVariable, const VariableData& rRotationalVariable)
// {

// }

/* Private functions ****************************************************/

// template<class TVariable, class TDerivedVariable>
// void StandardRecoveryUtility::AddPartialTimeDerivative(const TVariable& rVariable, const TDerivedVariable& rTimeDerivativeVariable)
// {
//     KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
// }

// template<>
// void StandardRecoveryUtility::AddPartialTimeDerivative(const DoubleVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable)
// {
//     this->AddScalarPartialTimeDerivative<DoubleVarType>(rVariable, rTimeDerivativeVariable);
// }

// template<class ComponentVarType, class DoubleVarType>
// void StandardRecoveryUtility::AddPartialTimeDerivative(const ComponentVarType& rVariable, const DoubleVarType& rTimeDerivativeVariable)
// {
//     this->AddScalarPartialTimeDerivative<ComponentVarType>(rVariable, rTimeDerivativeVariable);
// }



template<class TScalarVariable>
void StandardRecoveryUtility::AddScalarPartialTimeDerivative(const TScalarVariable& rVariable, const DoubleVarType& rTimeDerivativeVariable)
{
    const double delta_time_inv = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        double& material_derivative = inode->FastGetSolutionStepValue(rTimeDerivativeVariable);
        const double eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(rVariable)
                                                                 - inode->FastGetSolutionStepValue(rVariable, 1));
        material_derivative += eulerian_rate_of_change;
    }
}

// template<class TVariable>
// void StandardRecoveryUtility::CalculateScalarGradient(const TVariable& rScalarVariable, const TVariable& rGradientVariable)
// {
//     KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
// }

template<class TScalarVariable> //TODO: generalize for any element
void StandardRecoveryUtility::CalculateScalarGradient(const TScalarVariable& rScalarVariable, const ArrayVarType& rGradientVariable)
{
    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(rGradientVariable)) = ZeroVector(3);
    }

    const int TDim = 3;
    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;

    for (auto ielem = mrModelPart.ElementsBegin(); ielem != mrModelPart.ElementsEnd(); ++ielem){

        // computing the shape function derivatives

        Geometry<Node<3> >& geom = ielem->GetGeometry();
        double Volume;

        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

        // getting the gradients;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            elemental_values[i] = geom[i].FastGetSolutionStepValue(rScalarVariable);
        }

        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

        for (unsigned int i = 0; i < TDim; ++i){
            grad[i] = grad_aux[i];
        }

        double nodal_area = Volume / static_cast<double>(TDim + 1);
        grad *= nodal_area;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(rGradientVariable) += grad;
        }
    }

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        inode->FastGetSolutionStepValue(rGradientVariable) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}

// template<class TVariable, class TDerivedVariable>
// void StandardRecoveryUtility::CalculateLaplacian(const TVariable& rScalarVariable, const TDerivedVariable& rLaplacianVariable)
// {
//     KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
// }

// template<class TVariable, class TDerivedVariable>
// void StandardRecoveryUtility::CalculateLaplacian(const TVariable& rScalarVariable, const TDerivedVariable& rLaplacianVariable)
// {
//     KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
// }

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarMaterialDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rMaterialDerivativeVariable)
{
}

// template<class TVariable, class TDerivedVariable>
// void StandardRecoveryUtility::CalculateMaterialDerivative(const TVariable& rVariable, const TDerivedVariable& rMaterialDerivativeVariable)
// {
//     KRATOS_THROW_ERROR(std::invalid_argument, "Wrong combination.", "");
// }

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<DoubleVarType>(const DoubleVarType&, const ArrayVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<ComponentVarType>(const ComponentVarType&, const ArrayVarType&);

// template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
// template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

}  // namespace Kratos.
