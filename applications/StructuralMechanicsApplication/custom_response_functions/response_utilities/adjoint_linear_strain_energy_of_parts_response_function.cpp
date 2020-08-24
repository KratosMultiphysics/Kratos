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
#include "adjoint_linear_strain_energy_of_parts_response_function.h"

namespace Kratos
{
    AdjointLinearStrainEnergyOfPartsResponseFunction::AdjointLinearStrainEnergyOfPartsResponseFunction(ModelPart& rModelPart, Parameters ResponseSettings)
    : AdjointStructuralResponseFunction(rModelPart, ResponseSettings)
    {
        mCriticalPartName = ResponseSettings["critical_part_name"].GetString();

        // Collect relevant elements
        ModelPart& adjoint_aggregation_part = rModelPart.GetSubModelPart(mCriticalPartName);
        for(auto& elem : adjoint_aggregation_part.Elements()) {
            mAggregatedElementIds.push_back(elem.Id());
        }
    }

    AdjointLinearStrainEnergyOfPartsResponseFunction::~AdjointLinearStrainEnergyOfPartsResponseFunction()
    {
    }


    void AdjointLinearStrainEnergyOfPartsResponseFunction::CalculateGradient(const Element& rAdjointElement,
                                                             const Matrix& rResidualGradient,
                                                             Vector& rResponseGradient,
                                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        if (rResponseGradient.size() != rResidualGradient.size1()) {
            rResponseGradient.resize(rResidualGradient.size1(), false);
        }
        rResponseGradient.clear();

        Vector primal_disp;
        Vector RHS;
        Matrix LHS;

        if(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end()) {
            Element::Pointer pElem = mrModelPart.pGetElement(rAdjointElement.Id());
            pElem->Calculate(PRIMAL_VALUES_VECTOR, primal_disp, rProcessInfo);
            pElem->CalculateLocalSystem(LHS, RHS, rProcessInfo);
            noalias(rResponseGradient) = -1.0 * prod(LHS, primal_disp);
        }

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        //KRATOS_ERROR_IF(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        //        << "AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity not implemented yet!" << std::endl;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<double>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity(Element& rAdjointElement,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        //KRATOS_ERROR_IF(std::find(mAggregatedElementIds.begin(), mAggregatedElementIds.end(), rAdjointElement.Id()) != mAggregatedElementIds.end())
        //        << "AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity not implemented yet!" << std::endl;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    void AdjointLinearStrainEnergyOfPartsResponseFunction::CalculatePartialSensitivity(Condition& rAdjointCondition,
                                             const Variable<array_1d<double, 3>>& rVariable,
                                             const Matrix& rSensitivityMatrix,
                                             Vector& rSensitivityGradient,
                                             const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());

        KRATOS_CATCH("");
    }

    double AdjointLinearStrainEnergyOfPartsResponseFunction::CalculateValue(ModelPart& rModelPart)
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

        ModelPart& primal_aggregation_part = rModelPart.GetSubModelPart(mCriticalPartName);

        for (auto& elem_i : primal_aggregation_part.Elements()) {
            elem_i.GetValuesVector(disp,0);

            elem_i.CalculateLocalSystem(LHS, RHS, r_current_process_info);

            // Compute linear strain energy 0.5*u*K*u
            response_value += 0.5 * inner_prod(disp, prod(LHS,disp));
         }

        return response_value;

        KRATOS_CATCH("");
    }
} // namespace Kratos.


