// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

// System includes

// External includes

// Project includes
#include "adjoint_strain_energy_response_function.h"

namespace Kratos
{
    AdjointStrainEnergyResponseFunction::AdjointStrainEnergyResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
    }

    AdjointStrainEnergyResponseFunction::~AdjointStrainEnergyResponseFunction()
    {
    }

    void AdjointStrainEnergyResponseFunction::Initialize()
    {
        KRATOS_TRY;

        BaseType::Initialize();

        // It is necessary to initialize the elements/conditions since no adjoint problem is solved for this response type.
        ModelPart& r_model_part = this->GetModelPart();

        #pragma omp parallel for
        for(int i=0; i< static_cast<int>(r_model_part.Elements().size()); ++i)
        {
            auto it = r_model_part.ElementsBegin() + i;
            it->Initialize();
        }
        #pragma omp parallel for
        for(int i=0; i< static_cast<int>(r_model_part.Conditions().size()); ++i)
        {
            auto it = r_model_part.ConditionsBegin() + i;
            it->Initialize();
        }

        KRATOS_CATCH("");
    }

    double AdjointStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        ProcessInfo &r_current_process_info = rModelPart.GetProcessInfo();
        double response_value = 0.0;

        // Check if there are at the time of calling adjoint or primal elements
        KRATOS_ERROR_IF( r_current_process_info[IS_ADJOINT] )
             << "Calculate value for strain energy response is not available when using adjoint elements" << std::endl;

        // Sum all elemental strain energy values calculated as: W_e = u_e^T K_e u_e
        Matrix LHS;
        Vector RHS;
        Vector disp;

        for (auto& elem_i : rModelPart.Elements())
        {
            // Get state solution relevant for energy calculation
            elem_i.GetValuesVector(disp,0);

            elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

            // Compute strain energy
            response_value += 0.5 * inner_prod(disp, prod(LHS,disp));
         }

        return response_value;

        KRATOS_CATCH("");
    }

    void AdjointStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        // There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

          KRATOS_CATCH("")
    }

    void AdjointStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        if (rResponseGradient.size() != rDerivativesMatrix.size1())
              rResponseGradient.resize(rDerivativesMatrix.size1(), false);
        rResponseGradient.clear();

        // There will be a mistake, if body forces are considered. Because the elements are responsible for the body forces!

        KRATOS_CATCH("")
    }

    void AdjointStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        Vector adjoint_variables;

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rDerivativesMatrix \frac{\partial F}{\partial s}

        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
            << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    void AdjointStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        Vector adjoint_variables;

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rDerivativesMatrix \frac{\partial F}{\partial s}

        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
             << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

} // namespace Kratos.


