//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "includes/mat_variables.h"
#include "utilities/parallel_utilities.h"
#include "processes/set_initial_state_process.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

    template<std::size_t TDim>
    SetInitialStateProcess<TDim>::SetInitialStateProcess(
        ModelPart& rModelPart)
            : mrModelPart(rModelPart)
    {
        KRATOS_TRY

        const std::size_t voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        noalias(mInitialStrain) = ZeroVector(voigt_size);
        noalias(mInitialStress) = ZeroVector(voigt_size);
        noalias(mInitialF)      = ZeroMatrix(TDim, TDim);

        KRATOS_CATCH("")
    }

/***********************************************************************************/
/***********************************************************************************/

    template<std::size_t TDim>
    SetInitialStateProcess<TDim>::SetInitialStateProcess(
        ModelPart& rModelPart,
        const Vector& rInitialStrain,
        const Vector& rInitialStress,
        const Matrix& rInitialF) :
        mrModelPart(rModelPart), mInitialStrain(rInitialStrain),
        mInitialStress(rInitialStress), mInitialF(rInitialF)
    {
    }

/***********************************************************************************/
/***********************************************************************************/

    template<std::size_t TDim>
    SetInitialStateProcess<TDim>::SetInitialStateProcess(
        ModelPart& rModelPart,
        const Vector& rInitialStateVector,
        const int InitialStateType)
            : mrModelPart(rModelPart)
    {
        KRATOS_TRY

        const std::size_t voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        if (InitialStateType == 0) {
            noalias(mInitialStrain) = rInitialStateVector;
            noalias(mInitialStress) = ZeroVector(voigt_size);
        } else if (InitialStateType == 1) {
            noalias(mInitialStrain) = ZeroVector(voigt_size);
            noalias(mInitialStress) = rInitialStateVector;
        } else {
            noalias(mInitialStrain) = ZeroVector(voigt_size);
            noalias(mInitialStress) = ZeroVector(voigt_size);
        }
        noalias(mInitialF) = ZeroMatrix(TDim, TDim);

        KRATOS_CATCH("")
    }

/***********************************************************************************/
/***********************************************************************************/

    template<std::size_t TDim>
    SetInitialStateProcess<TDim>::SetInitialStateProcess(
        ModelPart& rModelPart,
        const Matrix& rInitialStateF)
            : mrModelPart(rModelPart)
    {
        KRATOS_TRY

        const std::size_t voigt_size = (TDim == 3) ? 6 : 3;

        mInitialStrain.resize(voigt_size, false);
        mInitialStress.resize(voigt_size, false);
        mInitialF.resize(TDim, TDim, false);

        noalias(mInitialStrain) = ZeroVector(voigt_size);
        noalias(mInitialStress) = ZeroVector(voigt_size);
        noalias(mInitialF)      = rInitialStateF;

        KRATOS_CATCH("")
    }

/***********************************************************************************/
/***********************************************************************************/

    template<std::size_t TDim>
    void SetInitialStateProcess<TDim>::ExecuteInitializeSolutionStep()
    {
        KRATOS_TRY

        Vector aux_initial_strain = mInitialStrain;
        Vector aux_initial_stress = mInitialStress;
        Matrix aux_initial_F      = mInitialF;
        InitialState::Pointer p_initial_state = Kratos::make_intrusive<InitialState>
            (aux_initial_strain, aux_initial_stress, aux_initial_F);

        block_for_each(mrModelPart.Elements(), [&](Element& r_element){
            // If the values are set element-wise have priority
            bool requires_unique_initial_state = false;
            if (r_element.GetGeometry().Has(INITIAL_STRAIN_VECTOR)) {
                noalias(aux_initial_strain) = (r_element.GetGeometry()).GetValue(INITIAL_STRAIN_VECTOR);
                requires_unique_initial_state = true;
            }
            if (r_element.GetGeometry().Has(INITIAL_STRESS_VECTOR)) {
                noalias(aux_initial_stress) = (r_element.GetGeometry()).GetValue(INITIAL_STRESS_VECTOR);
                requires_unique_initial_state = true;
            }
            if (r_element.GetGeometry().Has(INITIAL_DEFORMATION_GRADIENT_MATRIX)) {
                noalias(aux_initial_F) = (r_element.GetGeometry()).GetValue(INITIAL_DEFORMATION_GRADIENT_MATRIX);
                requires_unique_initial_state = true;
            }

            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector;
            r_element.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector, mrModelPart.GetProcessInfo());
            const auto& r_integration_points = r_element.GetGeometry().IntegrationPoints(r_element.GetIntegrationMethod());

            if (requires_unique_initial_state) {
                InitialState::Pointer p_initial_state_custom = Kratos::make_intrusive<InitialState>(aux_initial_strain, aux_initial_stress, aux_initial_F);
                for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                    constitutive_law_vector[point_number]->SetInitialState(p_initial_state_custom);
                }
            } else {
                for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
                    constitutive_law_vector[point_number]->SetInitialState(p_initial_state);
                }
            }
        });
        KRATOS_CATCH("")
    }

/***********************************************************************************/
/***********************************************************************************/

    template class SetInitialStateProcess<2>;
    template class SetInitialStateProcess<3>;

}  // namespace Kratos.
