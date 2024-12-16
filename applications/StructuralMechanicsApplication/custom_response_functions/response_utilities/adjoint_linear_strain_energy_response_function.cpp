// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "utilities/entities_utilities.h"
#include "adjoint_linear_strain_energy_response_function.h"

namespace Kratos
{
    AdjointLinearStrainEnergyResponseFunction::AdjointLinearStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
    }

    AdjointLinearStrainEnergyResponseFunction::~AdjointLinearStrainEnergyResponseFunction()
    {
    }

    void AdjointLinearStrainEnergyResponseFunction::Initialize()
    {
        KRATOS_TRY;

        BaseType::Initialize();

        // It is necessary to initialize the elements/conditions since no adjoint problem is solved for this response type.
        EntitiesUtilities::InitializeAllEntities(mrModelPart);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        this->CheckForBodyForces(rAdjointElement);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityMatrix.size1() == 0)
        {
            if (rSensitivityGradient.size() != 0)
                rSensitivityGradient.resize(0, false);
            return;
        }

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rSensitivityMatrix \frac{\partial F}{\partial s}

        Vector adjoint_variables;
        const auto& r_const_adjoint_condition_ref = rAdjointCondition;
        r_const_adjoint_condition_ref.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rSensitivityMatrix.size2())
             << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

        noalias(rSensitivityGradient) = prod(rSensitivityMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        this->CheckForBodyForces(rAdjointElement);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rSensitivityMatrix.size1() == 0)
        {
            if (rSensitivityGradient.size() != 0)
                rSensitivityGradient.resize(0, false);
            return;
        }

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rSensitivityMatrix \frac{\partial F}{\partial s}

        Vector adjoint_variables;
        const auto& r_const_adjoint_condition_ref = rAdjointCondition;
        r_const_adjoint_condition_ref.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rSensitivityMatrix.size2())
            << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

        noalias(rSensitivityGradient) = prod(rSensitivityMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    double AdjointLinearStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        const ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        double response_value = 0.0;

        // Check if there are at primal elements, because the primal state is required
        KRATOS_ERROR_IF( r_current_process_info.Has(IS_ADJOINT) && r_current_process_info[IS_ADJOINT] )
             << "Calculate value for strain energy response is only available when using primal elements" << std::endl;

        // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
        Matrix LHS;
        Vector RHS;
        Vector disp;

        for (auto& elem_i : rModelPart.Elements())
        {
            const auto& r_const_elem_ref = elem_i;
            r_const_elem_ref.GetValuesVector(disp,0);

            elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

            // Compute linear strain energy 0.5*u*K*u
            response_value += 0.5 * inner_prod(disp, prod(LHS,disp));
         }

        return response_value;

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CheckForBodyForces(Element& rAdjointElement)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        Vector acc = ZeroVector(3);
        if (rAdjointElement.GetProperties().Has( VOLUME_ACCELERATION ))
            acc+= rAdjointElement.GetProperties()[VOLUME_ACCELERATION];

        if( rAdjointElement.GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION))
            acc += rAdjointElement.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        KRATOS_ERROR_IF( norm_2(acc)>numerical_limit )
                << "linear strain energy response is not able to treat structures with self-weight correctly!" << std::endl;
    }

} // namespace Kratos.

