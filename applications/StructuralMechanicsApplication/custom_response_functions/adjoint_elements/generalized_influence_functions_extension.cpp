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

        mNormalize = AnalysisSettings["normalize"].GetBool();

        const std::string differentiation_method = AnalysisSettings["differentiation_method"].GetString();
        if(differentiation_method == "finite_differences")
        {
            mDifferentiationMethod = 1;
            mDelta = AnalysisSettings["delta"].GetDouble();
            mAdaptStepSize = AnalysisSettings["adapt_step_size"].GetBool();
        }
        else if (differentiation_method == "modify_material_matrix")
        {
            mDifferentiationMethod = 2;
        }
        else if (differentiation_method == "chain_rule")
        {
            mDifferentiationMethod = 3;
            mDelta = AnalysisSettings["delta"].GetDouble();
            mAdaptStepSize = AnalysisSettings["adapt_step_size"].GetBool();
        }
        else
            KRATOS_ERROR << "unknown differentiation method provided!" << std::endl;
    }

    GeneralizedInfluenceFunctionsExtension::~GeneralizedInfluenceFunctionsExtension()
    {
    }

    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityOnIntegrationPoints(Element& rElement,
                        const Variable<array_1d<double, 3>>& rPseudoQuantityVariable,
                        std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);

        if (rPseudoQuantityVariable == PSEUDO_MOMENT)
        {
            if (mDifferentiationMethod == 1)
            {
                this->CalculatePseudoQuantityWithFiniteDifferences(rElement, MOMENT, r_design_variable, rOutput, rCurrentProcessInfo);
            }
            else if(mDifferentiationMethod == 2)
            {
                this->CalculatePseudoQuantityByModificationOfMaterialMatrix(rElement, MOMENT, rOutput, rCurrentProcessInfo);
            }
            else if (mDifferentiationMethod == 3)
            {
                this->CalculatePseudoQuantityWithChainRule(rElement, MOMENT, rOutput, rCurrentProcessInfo);
            }

        }
        else if (rPseudoQuantityVariable == PSEUDO_FORCE)
        {
            if (mDifferentiationMethod == 1)
            {
                this->CalculatePseudoQuantityWithFiniteDifferences(rElement, FORCE, r_design_variable, rOutput, rCurrentProcessInfo);
            }
            else if(mDifferentiationMethod == 2)
            {
                this->CalculatePseudoQuantityByModificationOfMaterialMatrix(rElement, FORCE, rOutput, rCurrentProcessInfo);
            }
            else if (mDifferentiationMethod == 3)
            {
                this->CalculatePseudoQuantityWithChainRule(rElement, FORCE, rOutput, rCurrentProcessInfo);
            }
        }
        else
            KRATOS_ERROR << "It is not possible to provide a pseudo quantity for: " << rPseudoQuantityVariable.Name() << "!" << std::endl;

        if(mNormalize)
        {
            const SizeType write_points_number = rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod());
            const double variable_value = this->GetVariableValue(rElement, r_design_variable, rCurrentProcessInfo);
            for(IndexType i = 0; i < write_points_number; ++i)
                rOutput[i] *= variable_value;
        }

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
                rOutput[i] = pseudo_moment[i][0] * adjoint_curvature[i][0] - pseudo_force[i][0] * adjoint_strain[i][0] +
                             pseudo_moment[i][1] * adjoint_curvature[i][1] + pseudo_force[i][1] * adjoint_strain[i][1] +
                             pseudo_moment[i][2] * adjoint_curvature[i][2] + pseudo_force[i][2] * adjoint_strain[i][2];
            }

            // This is something special for normal stress response
            if (rAdjointElement.Has(RESPONSE_PREFACTOR_MOMENT_DERIVED) || rAdjointElement.Has(RESPONSE_PREFACTOR_FORCE_DERIVED))
            {
                const auto prefactor_moment = rAdjointElement.GetValue(RESPONSE_PREFACTOR_MOMENT_DERIVED);
                const auto prefactor_force = rAdjointElement.GetValue(RESPONSE_PREFACTOR_FORCE_DERIVED);

                std::vector< array_1d<double, 3> > moment;
                moment.resize(write_points_number);
                std::vector< array_1d<double, 3> > force;
                force.resize(write_points_number);
                std::vector< array_1d<double, 3> > adjoint_particular_curvature;
                adjoint_particular_curvature.resize(write_points_number);
                std::vector< array_1d<double, 3> > adjoint_particular_strain;
                adjoint_particular_strain.resize(write_points_number);

                rPrimalElement.CalculateOnIntegrationPoints(MOMENT, moment, rCurrentProcessInfo);
                rPrimalElement.CalculateOnIntegrationPoints(FORCE, force, rCurrentProcessInfo);
                rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_PARTICULAR_CURVATURE, adjoint_particular_curvature, rCurrentProcessInfo);
                rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_PARTICULAR_STRAIN, adjoint_particular_strain, rCurrentProcessInfo);

                double response_value = 1.0;
                double variable_value = 1.0;
                if(mNormalize)
                {
                    if (rCurrentProcessInfo.Has(RESPONSE_VALUE)) {response_value = std::abs(rCurrentProcessInfo.GetValue(RESPONSE_VALUE));}
                    else {KRATOS_ERROR << "Can't normalize variational sensitivity since no response value is provided!" << std::endl;}
                    const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);
                    variable_value = this->GetVariableValue(rAdjointElement, r_design_variable, rCurrentProcessInfo);
                }

                for(IndexType i = 0; i < write_points_number; ++i)
                { // TODO evaluate signs!
                    rOutput[i] += (moment[i][0] * adjoint_particular_curvature[i][0] * prefactor_moment[0]
                                + moment[i][1] * adjoint_particular_curvature[i][1] * prefactor_moment[1]
                                + moment[i][2] * adjoint_particular_curvature[i][2] * prefactor_moment[2]
                                - force[i][0] * adjoint_particular_strain[i][0] * prefactor_force[0]
                                + force[i][1] * adjoint_particular_strain[i][1] * prefactor_force[1]
                                + force[i][2] * adjoint_particular_strain[i][2] * prefactor_force[2])
                                * variable_value / response_value;
                }
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
                rOutput[i] = pseudo_force[i][0] * adjoint_strain[i][0] +
                             pseudo_force[i][1] * adjoint_strain[i][1] +
                             pseudo_force[i][2] * adjoint_strain[i][2];
            }
        }
        else
            KRATOS_ERROR << "CalculateSensitivityOnIntegrationPoints not available for " << primal_element_name << "!" << std::endl;

        KRATOS_CATCH("");
    }

    void GeneralizedInfluenceFunctionsExtension::CalculateAdjointWorkContributionOnIntegrationPoints(Element& rPrimalElement, Element& rAdjointElement,
                                                                                                const Variable<array_1d<double, 3>>& rAdjointWorkVariable,
                                                                                                std::vector< array_1d<double, 3> >& rOutput,
                                                                                                const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const SizeType write_points_number = rAdjointElement.GetGeometry().IntegrationPointsNumber(rAdjointElement.GetIntegrationMethod());
        if (rOutput.size() != write_points_number)
            rOutput.resize(write_points_number);
        std::vector< array_1d<double, 3> > pseudo_quantity;
        pseudo_quantity.resize(write_points_number);
        std::vector< array_1d<double, 3> > adjoint_quantity;
        adjoint_quantity.resize(write_points_number);

        if (rAdjointWorkVariable == ADJOINT_WORK_FORCE_CONTRIBUTION)
        {
            this->CalculatePseudoQuantityOnIntegrationPoints(rPrimalElement, PSEUDO_FORCE, pseudo_quantity, rCurrentProcessInfo);
            rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_STRAIN, adjoint_quantity, rCurrentProcessInfo);
        }
        else if (rAdjointWorkVariable == ADJOINT_WORK_MOMENT_CONTRIBUTION)
        {
            this->CalculatePseudoQuantityOnIntegrationPoints(rPrimalElement, PSEUDO_MOMENT, pseudo_quantity, rCurrentProcessInfo);
            rAdjointElement.CalculateOnIntegrationPoints(ADJOINT_CURVATURE, adjoint_quantity, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "CalculateAdjointWorkContributionOnIntegrationPoints not available for given variable!" << std::endl;

        for(IndexType i = 0; i < write_points_number; ++i)
        {
            rOutput[i][0] = pseudo_quantity[i][0] * adjoint_quantity[i][0];
            rOutput[i][1] = pseudo_quantity[i][1] * adjoint_quantity[i][1];
            rOutput[i][2] = pseudo_quantity[i][2] * adjoint_quantity[i][2];
        }

        KRATOS_CATCH("");
    }

    void GeneralizedInfluenceFunctionsExtension::NormalizeAdjointFieldIfRequested(Element& rElement, std::vector< array_1d<double, 3> >& rOutput,
                                                                                            const ProcessInfo& rCurrentProcessInfo) const
    {
        if(mNormalize)
        {
            if (rCurrentProcessInfo.Has(RESPONSE_VALUE))
            {
                const double response_value = std::abs(rCurrentProcessInfo.GetValue(RESPONSE_VALUE));
                for(IndexType i = 0; i < rOutput.size(); ++i)
                    rOutput[i] /= response_value;
            }
            else
                KRATOS_ERROR << "Can't normalize adjoint field since no response value is provided!" << std::endl;
        }
    }

    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityWithFiniteDifferences(Element& rElement,
                        const Variable<array_1d<double, 3>>& rQuantityVariable, const Variable<double>& rDesignVariable,
                        std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        double perturbation_size = mDelta;
        if (mAdaptStepSize)
        {
            if(rElement.GetProperties().Has(rDesignVariable))
                perturbation_size *= rElement.GetProperties()[rDesignVariable];
        }

        ElementFiniteDifferenceUtility::CalculateIntegrationPointsResultsDerivative(rElement, rQuantityVariable,
            rDesignVariable, perturbation_size, rOutput, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }


    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityWithChainRule(Element& rElement,
                        const Variable<array_1d<double, 3>>& rQuantityVariable,
                        std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        std::string primal_element_name;
        CompareElementsAndConditionsUtility::GetRegisteredName(rElement, primal_element_name);

        if(primal_element_name == "CrLinearBeamElement3D2N")
        {
            std::vector<Variable<double> > explicit_variables = {YOUNG_MODULUS, CROSS_AREA, I22, I33,
                                            AREA_EFFECTIVE_Y, AREA_EFFECTIVE_Z, TORSIONAL_INERTIA, POISSON_RATIO};
            std::vector<Variable<double> > derived_variables = {YOUNG_MODULUS_DERIVED, CROSS_AREA_DERIVED, I22_DERIVED, I33_DERIVED,
                                            AREA_EFFECTIVE_Y_DERIVED, AREA_EFFECTIVE_Z_DERIVED, TORSIONAL_INERTIA_DERIVED, POISSON_RATIO_DERIVED};

            const SizeType write_points_number = rElement.GetGeometry().IntegrationPointsNumber(rElement.GetIntegrationMethod());
            if (rOutput.size() != write_points_number)
                rOutput.resize(write_points_number);
            for(IndexType j = 0; j < write_points_number; ++j)
                rOutput[j].clear();

            for (IndexType i = 0; i < explicit_variables.size(); ++i)
            {
                if(rElement.GetProperties().Has(explicit_variables[i]) && rElement.GetProperties().Has(derived_variables[i]))
                {
                    std::vector< array_1d<double, 3> > pseudo_quantity;
                    pseudo_quantity.resize(write_points_number);
                    this->CalculatePseudoQuantityWithFiniteDifferences(rElement, rQuantityVariable, explicit_variables[i], pseudo_quantity, rCurrentProcessInfo);

                    const double variable_derivative = rElement.GetProperties()[derived_variables[i]];

                    for(IndexType gp = 0; gp < write_points_number; ++gp)
                    {
                        for(IndexType comp = 0; comp < 3; ++comp)
                            rOutput[gp][comp] += pseudo_quantity[gp][comp] * variable_derivative; // here the chain rule is applied (dM/dX * DX/ds)
                    }
                }
            }
        }
        else
            KRATOS_ERROR << "CalculatePseudoQuantityWithChainRule is currently only implemented for CrLinearBeamElement3D2N!" << std::endl;

        KRATOS_CATCH("");
    }


    void GeneralizedInfluenceFunctionsExtension::CalculatePseudoQuantityByModificationOfMaterialMatrix(Element& rElement,
                        const Variable<array_1d<double, 3>>& rQuantityVariable,
                        std::vector< array_1d<double, 3> >& rOutput, const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;

        const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);

        // Save property pointer
        Properties::Pointer p_global_properties = rElement.pGetProperties();

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        rElement.SetProperties(p_local_property);

        // Modify material matrix
        if (r_design_variable == SECTION_HEIGTH_SENSITIVITY || r_design_variable == SECTION_WIDTH_SENSITIVITY)
        {
            KRATOS_ERROR_IF(rElement.GetProperties().Has(AREA_EFFECTIVE_Y) || rElement.GetProperties().Has(AREA_EFFECTIVE_Z))
                << "It not possible to use CalculatePseudoQuantityByModificationOfMaterialMatrix for Timoshenko beam formulation!" << std::endl;

            std::string primal_element_name;
            CompareElementsAndConditionsUtility::GetRegisteredName(rElement, primal_element_name);

            if(rElement.GetProperties().Has(CROSS_AREA_DERIVED))
            {
                const double value_derived = rElement.GetProperties()[CROSS_AREA_DERIVED];
                p_local_property->SetValue(CROSS_AREA, value_derived);
            }
            else
                KRATOS_ERROR << "It is not possible to compute pseudo stress quantity since CROSS_AREA_DERIVED is not available" << std::endl;

            if(primal_element_name == "CrLinearBeamElement3D2N")
            {
                if(rElement.GetProperties().Has(I22_DERIVED))
                {
                    const double value_derived = rElement.GetProperties()[I22_DERIVED];
                    p_local_property->SetValue(I22, value_derived);
                }
                else
                    KRATOS_ERROR << "It is not possible to compute pseudo stress quantity since I22_DERIVED is not available" << std::endl;

                if(rElement.GetProperties().Has(I33_DERIVED))
                {
                    const double value_derived = rElement.GetProperties()[I33_DERIVED];
                    p_local_property->SetValue(I33, value_derived);
                }
                else
                    KRATOS_ERROR << "It is not possible to compute pseudo stress quantity since I33_DERIVED is not available" << std::endl;

                if(rElement.GetProperties().Has(TORSIONAL_INERTIA_DERIVED))
                {
                    const double value_derived = rElement.GetProperties()[TORSIONAL_INERTIA_DERIVED];
                    p_local_property->SetValue(TORSIONAL_INERTIA, value_derived);
                }
                else
                    KRATOS_ERROR << "It is not possible to compute pseudo stress quantity since TORSIONAL_INERTIA_DERIVED is not available" << std::endl;
            }
        }
        else
            KRATOS_ERROR << "CalculatePseudoQuantityByModificationOfMaterialMatrix not available for given variable!" << std::endl;

        rElement.CalculateOnIntegrationPoints(rQuantityVariable, rOutput, rCurrentProcessInfo);

        // Give element original properties back
        rElement.SetProperties(p_global_properties);

        KRATOS_CATCH("");
    }

    double GeneralizedInfluenceFunctionsExtension::GetVariableValue(Element& rElement, const Variable<double>& rDesignVariable,
                                                                    const ProcessInfo& rCurrentProcessInfo) const
    {
        double variable_value = 1.0;
        if (rDesignVariable == SECTION_HEIGTH_SENSITIVITY)
        {
            if (rElement.GetProperties().Has(SECTION_HEIGTH))
                variable_value = rElement.GetProperties()[SECTION_HEIGTH];
        }
        else if (rDesignVariable == SECTION_WIDTH_SENSITIVITY)
        {
            if (rElement.GetProperties().Has(SECTION_WIDTH))
                variable_value = rElement.GetProperties()[SECTION_WIDTH];
        }
        else
            variable_value = rElement.GetProperties()[rDesignVariable];

        return variable_value;

    }

};  // namespace Kratos.

