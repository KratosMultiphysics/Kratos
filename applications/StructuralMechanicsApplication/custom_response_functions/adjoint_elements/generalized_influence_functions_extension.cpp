//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Martin Fusseder, https://github.com/MFusseder
//
//


// Project includes
#include "generalized_influence_functions_extension.h"
#include "custom_response_functions/response_utilities/element_finite_difference_utility.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/compare_elements_and_conditions_utility.h"


namespace Kratos
{
    GeneralizedInfluenceFunctionsExtension::GeneralizedInfluenceFunctionsExtension(Parameters AnalysisSettings)
    {
        const std::string variable_type = AnalysisSettings["variable_type"].GetString();
        KRATOS_ERROR_IF_NOT(variable_type == "element_property")
            << "The method of generalized influence functions is currently only implemented for stiffness related element properties" << std::endl;

        mDesignVariableName = AnalysisSettings["design_variable_name"].GetString();
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mDesignVariableName))
            << "Chosen design variable " << mDesignVariableName << " is not available or no variable of type double!" << std::endl;

        mDelta = AnalysisSettings["delta"].GetDouble();

        mNormalize = AnalysisSettings["normalize"].GetBool();

        mAdaptStepSize = AnalysisSettings["adapt_step_size"].GetBool();
    }

    GeneralizedInfluenceFunctionsExtension::~GeneralizedInfluenceFunctionsExtension()
    {
    }

    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityOnIntegrationPoints(Element& rElement,
                        const Variable<array_1d<double, 3>>& rPseudoQuantityVariable,
                        std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        if (rPseudoQuantityVariable == PSEUDO_MOMENT)
        {
            const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);

            double perturbation_size = mDelta;
            if (mAdaptStepSize) { perturbation_size *= rElement.GetProperties()[r_design_variable];}

            ElementFiniteDifferenceUtility::CalculateIntegrationPointsResultsDerivative(rElement, MOMENT,
                r_design_variable, perturbation_size, rOutput, rCurrentProcessInfo);
        }
        else if (rPseudoQuantityVariable == PSEUDO_FORCE)
        {
            const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);

            double perturbation_size = mDelta;
            if (mAdaptStepSize) { perturbation_size *= rElement.GetProperties()[r_design_variable];}

            ElementFiniteDifferenceUtility::CalculateIntegrationPointsResultsDerivative(rElement, FORCE,
                r_design_variable, perturbation_size, rOutput, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "It is not possible to provide a pseudo quantity for: " << rPseudoQuantityVariable.Name() << "!" << std::endl;

        KRATOS_CATCH("");
    }

    void GeneralizedInfluenceFunctionsExtension::CalculateSensitivityOnIntegrationPoints(Element& rPrimalElement, Element& rAdjointElement,
    std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const SizeType write_points_number = rAdjointElement.GetGeometry().IntegrationPointsNumber(rAdjointElement.GetIntegrationMethod());
        if (rOutput.size() != write_points_number)
            rOutput.resize(write_points_number);

        std::string primal_element_name;
        CompareElementsAndConditionsUtility::GetRegisteredName(rPrimalElement, primal_element_name);

        double response_value = 1.0;
        double variable_value = 1.0;
        if(mNormalize)
        {
            if (rCurrentProcessInfo.Has(RESPONSE_VALUE)) {response_value = std::abs(rCurrentProcessInfo.GetValue(RESPONSE_VALUE));}
            else {KRATOS_ERROR << "Can't normalize variational sensitivity since no response value is provided!" << std::endl;}
            const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);
            variable_value = rAdjointElement.GetProperties()[r_design_variable];
        }

        if(primal_element_name == "CrLinearBeamElement3D2N")
        {
            std::vector< array_1d<double, 3> > pseudo_moment;
            pseudo_moment.resize(write_points_number);
            std::vector< array_1d<double, 3> > pseudo_force;
            pseudo_force.resize(write_points_number);
            std::vector< array_1d<double, 3> > adjoint_curvature;
            adjoint_curvature.resize(write_points_number);
            std::vector< array_1d<double, 3> > adjoint_strain;
            adjoint_strain.resize(write_points_number);
            this->CalculatePseudoQuantityOnIntegrationPoints(rPrimalElement, PSEUDO_MOMENT, pseudo_moment, rCurrentProcessInfo);
            this->CalculatePseudoQuantityOnIntegrationPoints(rPrimalElement, PSEUDO_FORCE, pseudo_force, rCurrentProcessInfo);
            rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_CURVATURE, adjoint_curvature, rCurrentProcessInfo);
            rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_STRAIN, adjoint_strain, rCurrentProcessInfo);

            // MFusseder TODO investigate signs!!
            for(IndexType i = 0; i < write_points_number; ++i)
            {
                rOutput[i] = (pseudo_moment[i][0] * adjoint_curvature[i][0] - pseudo_force[i][0] * adjoint_strain[i][0] +
                              pseudo_moment[i][1] * adjoint_curvature[i][1] + pseudo_force[i][1] * adjoint_strain[i][1] +
                              pseudo_moment[i][2] * adjoint_curvature[i][2] + pseudo_force[i][2] * adjoint_strain[i][2])*
                              variable_value / response_value;
            }
        }
        else if(primal_element_name == "TrussLinearElement3D2N")
        {
            std::vector< array_1d<double, 3> > pseudo_force;
            pseudo_force.resize(write_points_number);
            std::vector< array_1d<double, 3> > adjoint_strain;
            adjoint_strain.resize(write_points_number);
            this->CalculatePseudoQuantityOnIntegrationPoints(rPrimalElement, PSEUDO_FORCE, pseudo_force, rCurrentProcessInfo);
            rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_STRAIN, adjoint_strain, rCurrentProcessInfo);

            for(IndexType i = 0; i < write_points_number; ++i)
            {
                rOutput[i] = (pseudo_force[i][0] * adjoint_strain[i][0] +
                              pseudo_force[i][1] * adjoint_strain[i][1] +
                              pseudo_force[i][2] * adjoint_strain[i][2])*
                              variable_value / response_value;

            }
        }
        else
            KRATOS_ERROR << "CalculateSensitivityOnIntegrationPoints not available for " << primal_element_name << "!" << std::endl;




        KRATOS_CATCH("");
    }

};  // namespace Kratos.

