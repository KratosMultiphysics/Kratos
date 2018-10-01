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
#include "adjoint_linear_strain_energy_response_function.h"
#include "linear_strain_energy_response_function.h"

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

    double AdjointLinearStrainEnergyResponseFunction::CalculateValue(ModelPart& rModelPart)
    {
        KRATOS_TRY;

        LinearStrainEnergyResponseFunction zero_order_response(rModelPart, mResponseSettings);

        return zero_order_response.CalculateValue();

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        if (rResponseGradient.size() != 0)
            rResponseGradient.resize(0, false);

        this->CheckForBodyForces(rAdjointElem);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Element& rAdjointElem,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the elements don't have F, they do not contribute here.

        if (rResponseGradient.size() != 0)
            rResponseGradient.resize(0, false);

        this->CheckForBodyForces(rAdjointElem);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<array_1d<double,3>>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rDerivativesMatrix.size2() == 0)
        {
            if (rResponseGradient.size() != 0)
                rResponseGradient.resize(0, false);
            return;
        }

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rDerivativesMatrix \frac{\partial F}{\partial s}

        Vector adjoint_variables;
        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
            << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CalculateSensitivityGradient(Condition& rAdjointCondition,
                                              const Variable<double>& rVariable,
                                              const Matrix& rDerivativesMatrix,
                                              Vector& rResponseGradient,
                                              ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        if (rDerivativesMatrix.size2() == 0)
        {
            if (rResponseGradient.size() != 0)
                rResponseGradient.resize(0, false);
            return;
        }

        // The partial derivative of the linear strain energy is 0.5*u*\frac{\partial F}{\partial s}
        // Assuming that the conditions don't have K, the remaining content of rDerivativesMatrix \frac{\partial F}{\partial s}

        Vector adjoint_variables;
        rAdjointCondition.GetValuesVector(adjoint_variables); // = 0.5*u

        KRATOS_ERROR_IF(adjoint_variables.size() != rDerivativesMatrix.size2())
             << "Size of adjoint vector does not fit to the size of the pseudo load!" << std::endl;

        if (rResponseGradient.size() != rDerivativesMatrix.size2())
            rResponseGradient.resize(adjoint_variables.size(), false);

        noalias(rResponseGradient) = prod(rDerivativesMatrix, adjoint_variables);

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyResponseFunction::CheckForBodyForces(Element& rAdjointElem)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        Vector acc = ZeroVector(3);
        if (rAdjointElem.GetProperties().Has( VOLUME_ACCELERATION ))
            acc+= rAdjointElem.GetProperties()[VOLUME_ACCELERATION];

        if( rAdjointElem.GetGeometry()[0].SolutionStepsDataHas(VOLUME_ACCELERATION))
            acc += rAdjointElem.GetGeometry()[0].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        KRATOS_ERROR_IF( norm_2(acc)>numerical_limit )
                << "linear strain energy response is not able to treat structures with self-weight correctly!" << std::endl;
    }

} // namespace Kratos.


