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
#include "../swimming_DEM_application.h"
#include "derivative_recovery_process.h"
#include "standard_recovery_process.h"

namespace Kratos
{

StandardRecoveryProcess::StandardRecoveryProcess(
    ModelPart& rModelPart,
    Parameters Param)
    : DerivativeRecoveryProcess(rModelPart, Param)
{
    this->CheckDefaultsAndProcessSettings(Param);
    mStoreFullGradient = Param["store_full_gradient_option"].GetBool();
}

StandardRecoveryProcess::StandardRecoveryProcess(
    Model& rModel,
    Parameters Param)
    : DerivativeRecoveryProcess(rModel, Param)
{
    this->CheckDefaultsAndProcessSettings(Param);
    mStoreFullGradient = Param["store_full_gradient_option"].GetBool();
}

void StandardRecoveryProcess::CheckDefaultsAndProcessSettings(Parameters Param)
{
    Parameters default_parameters( R"(
    {
        "model_part_name" : "FluidModelPart",
        "recoverer_name" : "StandardRecoveryProcess",
        "store_full_gradient_option" : true
    }  )" );

    Param.ValidateAndAssignDefaults(default_parameters);
}

void StandardRecoveryProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void StandardRecoveryProcess::ExecuteInitialize() {

    KRATOS_TRY;


    KRATOS_CATCH("");
}

void StandardRecoveryProcess::ExecuteBeforeSolutionLoop() {
    this->ExecuteInitializeSolutionStep();
    this->ExecuteFinalizeSolutionStep();
}

void StandardRecoveryProcess::ExecuteInitializeSolutionStep() {

}

void StandardRecoveryProcess::ExecuteFinalizeSolutionStep() {

}

/* Protected functions ****************************************************/

void StandardRecoveryProcess::CalculateScalarGradients() {
}

/* Private functions ****************************************************/

void StandardRecoveryProcess::AddTimeDerivative()
{
    const double delta_time_inv = 1.0 / mrModelPart.GetProcessInfo()[DELTA_TIME];

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(mMaterialDerivativeContainer);
        const array_1d <double, 3> eulerian_rate_of_change = delta_time_inv * (inode->FastGetSolutionStepValue(VELOCITY)
                                                                               - inode->FastGetSolutionStepValue(VELOCITY, 1));
        noalias(material_derivative) += eulerian_rate_of_change;
    }
}

void StandardRecoveryProcess::CalculateVectorMaterialDerivative()
{
    KRATOS_INFO("SwimmingDEM") << "Constructing the material derivative by derivating nodal averages..." << std::endl;
    std::map <std::size_t, unsigned int> id_to_position;
    unsigned int entry = 0;

    for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
        noalias(inode->FastGetSolutionStepValue(mMaterialDerivativeContainer)) = ZeroVector(3);
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

        // for each element, constructing the gradient contribution (to its nodes) of the component v_j and storing it in mMaterialDerivativeContainer

        for (auto ielem = mrModelPart.ElementsBegin(); ielem != mrModelPart.ElementsEnd(); ++ielem){
            // computing the shape function derivatives
            Geometry<Node<3> >& geom = ielem->GetGeometry();
            double Volume;
            GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

            for (unsigned int i = 0; i < TDim + 1; ++i){
                elemental_values[i] = geom[i].FastGetSolutionStepValue(mMaterialDerivativeContainer)[j];
            }

            array_1d <double, 3> grad_aux = prod(trans(DN_DX), elemental_values); // its dimension may be 2

            for (unsigned int i = 0; i < TDim; ++i){
                grad[i] = grad_aux[i];
            }

            double nodal_area = Volume / static_cast<double>(TDim + 1);
            grad *= nodal_area;

            for (unsigned int i = 0; i < TDim + 1; ++i){
                geom[i].FastGetSolutionStepValue(mMaterialDerivativeContainer) += grad; // we use mMaterialDerivativeContainer to store the gradient of one component at a time
            }
        }

        // normalizing the contributions to the gradient and getting the j-component of the material derivative

        for (auto inode = mrModelPart.NodesBegin(); inode != mrModelPart.NodesEnd(); ++inode){
            array_1d <double, 3>& stored_gradient_of_component_j = inode->FastGetSolutionStepValue(mMaterialDerivativeContainer);
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
        array_1d <double, 3>& material_derivative = inode->FastGetSolutionStepValue(mMaterialDerivativeContainer);
        material_derivative = stored_convective_contribution;
    }

    // Adding Eulerian time derivative contribution

    AddTimeDerivative();

    KRATOS_INFO("SwimmingDEM") << "Finished constructing the material derivative by derivating nodal averages..." << std::endl;
}

};  // namespace Kratos.
