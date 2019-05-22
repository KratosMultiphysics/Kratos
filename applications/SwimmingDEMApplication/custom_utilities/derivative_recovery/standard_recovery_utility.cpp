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


void StandardRecoveryUtility::InitializeRecovery()
{

}

void StandardRecoveryUtility::FinalizeRecovery()
{

}

void StandardRecoveryUtility::AddPartialTimeDerivative(const DoubleVarType& rScalarVariable, const DoubleVarType& rTimeDerivativeVariable)
{
    this->AddScalarPartialTimeDerivative<DoubleVarType>(rScalarVariable, rTimeDerivativeVariable);
}

void StandardRecoveryUtility::AddPartialTimeDerivative(const ComponentVarType& rScalarVariable, const DoubleVarType& rTimeDerivativeVariable)
{
    this->AddScalarPartialTimeDerivative<ComponentVarType>(rScalarVariable, rTimeDerivativeVariable);
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

void StandardRecoveryUtility::CalculateGradient(const DoubleVarType& rScalarVariable, const ArrayVarType& rGradientVariable)
{
    this->CalculateScalarGradient<DoubleVarType>(rScalarVariable, rGradientVariable);
}

void StandardRecoveryUtility::CalculateGradient(const ComponentVarType& rScalarVariable, const ArrayVarType& rGradientVariable)
{
    this->CalculateScalarGradient<ComponentVarType>(rScalarVariable, rGradientVariable);
}

void StandardRecoveryUtility::CalculateGradient(const ArrayVarType& rVectorVariable, const TensorVarType& rGradientVariable)
{
    // dvj/dxi = (dv0/dx | dv1/dx | dv2/dx )_ij
    for (auto comp_j : {"X", "Y", "Z"}){
        const auto& variable_j = dynamic_cast<ComponentVarType&>(KratosComponents<VariableData>::Get(rVectorVariable.Name() + "_" + comp_j));

        this->CalculateScalarGradient<ComponentVarType>(variable_j, SCALAR_GRADIENT);

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            const auto& component_gradient = inode->FastGetSolutionStepValue(SCALAR_GRADIENT);
            int i = 0;
            for (auto comp_i : {"X", "Y", "Z"}){
                const auto& variable_ij = dynamic_cast<TensorComponentVarType&>(KratosComponents<VariableData>::Get(rGradientVariable.Name() + "_" + comp_i + comp_j));
                auto& vector_gradient_ij = inode->FastGetSolutionStepValue(variable_ij);
                vector_gradient_ij = component_gradient[i];
                ++i;
            }
        }
    }
}

void StandardRecoveryUtility::CalculateDivergence(const ArrayVarType& rVectorVariable, const DoubleVarType& rDivergenceVariable)
{
    this->CalculateDivergenceAsScalar<DoubleVarType>(rVectorVariable, rDivergenceVariable);
}

void StandardRecoveryUtility::CalculateDivergence(const ArrayVarType& rVectorVariable, const ComponentVarType& rDivergenceVariable)
{
    this->CalculateDivergenceAsScalar<ComponentVarType>(rVectorVariable, rDivergenceVariable);
}

void StandardRecoveryUtility::CalculateRotational(const ArrayVarType rVectorVariable, const ArrayVarType& rRotationalVariable)
{
    this->CalculateGradient(rVectorVariable, VECTOR_GRADIENT);
    constexpr auto &GetIndex = StandardRecoveryUtility::GetVectorizedMatrixIndex;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        const auto& vector_gradient = inode->FastGetSolutionStepValue(VECTOR_GRADIENT);
        auto& rotational = inode->FastGetSolutionStepValue(rRotationalVariable);
        rotational[0] = vector_gradient[GetIndex(1, 2)] - vector_gradient[GetIndex(2, 1)];
        rotational[1] = vector_gradient[GetIndex(2, 0)] - vector_gradient[GetIndex(0, 2)];
        rotational[2] = vector_gradient[GetIndex(0, 1)] - vector_gradient[GetIndex(1, 0)];
    }
}

void StandardRecoveryUtility::CalculateMaterialDerivative(const DoubleVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable)
{
    this->CalculateScalarMaterialDerivative<DoubleVarType>(rVariable, rMaterialDerivativeVariable);
}

void StandardRecoveryUtility::CalculateMaterialDerivative(const ComponentVarType& rVariable, const DoubleVarType& rMaterialDerivativeVariable)
{
    this->CalculateScalarMaterialDerivative<ComponentVarType>(rVariable, rMaterialDerivativeVariable);
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

void StandardRecoveryUtility::CalculateLaplacian(const DoubleVarType& rVariable, const DoubleVarType& rLaplacianVariable)
{
    this->CalculateScalarLaplacian<DoubleVarType>(rVariable, rLaplacianVariable);
}

void StandardRecoveryUtility::CalculateLaplacian(const ComponentVarType& rVariable, const DoubleVarType& rLaplacianVariable)
{
    this->CalculateScalarLaplacian<ComponentVarType>(rVariable, rLaplacianVariable);
}

void StandardRecoveryUtility::CalculateLaplacian(const ArrayVarType& rVectorComponent, const ArrayVarType& rLaplacianVariable)
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the Laplacian by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(rLaplacianVariable)) = ZeroVector(3);
        id_to_position[inode->Id()] = entry;
        ++entry;
    }

    std::vector<array_1d <double, 3> > laplacians;
    laplacians.resize(entry);
    std::fill(laplacians.begin(), laplacians.end(), ZeroVector(3));

    const int TDim = 3;
    array_1d <double, 3> grad = ZeroVector(3);
    array_1d <double, TDim + 1 > elemental_values;
    array_1d <double, TDim + 1 > N; // shape functions vector
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    BoundedMatrix<double, TDim + 1, TDim> elemental_vectors; // They carry the nodal gradients of the corresponding component v_j
    const double nodal_area_share = 1.0 / static_cast<double>(TDim + 1);

    for (unsigned int j = 0; j < TDim; ++j){ // for each component of the original vector value

        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in rLaplacianVariable

        for (auto ielem = mrModelPart.ElementsBegin(); ielem != mrModelPart.ElementsEnd(); ++ielem){

            // computing the shape function derivatives

            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(rVectorComponent)[j];
            }

            array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }

            double nodal_area = Volume * nodal_area_share;
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(rLaplacianVariable) += grad; // we use rLaplacianVariable to store the gradient of one component at a time
            }
        }

        // normalizing the contributions to the gradient

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            inode->FastGetSolutionStepValue(rLaplacianVariable) /= inode->FastGetSolutionStepValue(NODAL_AREA);
        }

        // for each element, constructing the divergence contribution (to its nodes) of the gradient of component v_j

        for (auto ielem = mrModelPart.ElementsBegin(); ielem != mrModelPart.ElementsEnd(); ++ielem){
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                for (unsigned int k = 0; k < TDim; ++k){
                    elemental_vectors(i, k) = geom[i].FastGetSolutionStepValue(rLaplacianVariable)[k]; // it is actually the gradient of component v_j
                }
            }

            BoundedMatrix<double, TDim, TDim> grad_aux = prod(trans(DN_DX), elemental_vectors); // its dimension may be 2
            double divergence_of_vi = 0.0;

            for (unsigned int k = 0; k < TDim; ++k){ // the divergence is the trace of the gradient
                divergence_of_vi += grad_aux(k, k);
            }

            double nodal_area = Volume / static_cast<double>(TDim + 1);
            divergence_of_vi *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                laplacians[id_to_position[geom[i].Id()]][j] += divergence_of_vi; // adding the contribution of the elemental divergence to each of its nodes
            }
        }

        // clearing the values stored in rLaplacianVariable for the next component

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            array_1d <double, 3>& current_gradient_of_vi = inode->FastGetSolutionStepValue(rLaplacianVariable);
            noalias(current_gradient_of_vi) = ZeroVector(3);
        }

    } // for each component of the vector value

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        array_1d <double, 3>& stored_laplacian = laplacians[id_to_position[inode->Id()]];
        array_1d <double, 3>& laplacian = inode->FastGetSolutionStepValue(rLaplacianVariable);
        noalias(laplacian) = stored_laplacian / inode->FastGetSolutionStepValue(NODAL_AREA);
    }

    KRATOS_INFO("SwimmingDEM") << "Finished constructing the Laplacian by derivating nodal averages..." << std::endl;
}

/* Private functions ****************************************************/

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

template<class TScalarVariable> //TODO: generalize for any element
void StandardRecoveryUtility::CalculateDivergenceAsScalar(const ArrayVarType& rVectorVariable, const TScalarVariable& rDivergenceVariable)
{
    this->CalculateGradient(rVectorVariable, VECTOR_GRADIENT);
    constexpr auto &GetIndex = StandardRecoveryUtility::GetVectorizedMatrixIndex;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        const auto& vector_gradient = inode->FastGetSolutionStepValue(VECTOR_GRADIENT);
        auto& divergence = inode->FastGetSolutionStepValue(rDivergenceVariable);
        divergence = 0.0;
        for (int i = 0; i < 3; ++i){
            divergence += vector_gradient[GetIndex(i, i)];
        }
    }
}

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarMaterialDerivative(const TScalarVariable& rScalarVariable, const DoubleVarType& rMaterialDerivativeVariable)
{

}

template<class TScalarVariable>
void StandardRecoveryUtility::CalculateScalarLaplacian(const TScalarVariable& rScalarVariable, const DoubleVarType& rLaplacianVariable)
{
    this->CalculateGradient(rScalarVariable, SCALAR_GRADIENT_ERROR);
    this->CalculateDivergence(SCALAR_GRADIENT_ERROR, rLaplacianVariable);
}


template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::AddScalarPartialTimeDerivative<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<DoubleVarType>(const DoubleVarType&, const ArrayVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarGradient<ComponentVarType>(const ComponentVarType&, const ArrayVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateDivergenceAsScalar<DoubleVarType>(const ArrayVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateDivergenceAsScalar<ComponentVarType>(const ArrayVarType&, const ComponentVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarMaterialDerivative<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<DoubleVarType>(const DoubleVarType&, const DoubleVarType&);
template void KRATOS_API(SWIMMING_DEM_APPLICATION) StandardRecoveryUtility::CalculateScalarLaplacian<ComponentVarType>(const ComponentVarType&, const DoubleVarType&);

}  // namespace Kratos.
