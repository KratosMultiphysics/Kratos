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

        const std::string differentiation_method = AnalysisSettings["differentiation_method"].GetString();
        if(differentiation_method == "finite_differences")
        {
            mDifferentiationMethod = 1;
            mDelta = AnalysisSettings["delta"].GetDouble();
            mAdaptStepSize = AnalysisSettings["adapt_step_size"].GetBool();
            mNormalize = AnalysisSettings["normalize"].GetBool(); //normalization is only possible for this differentiation method
        }
        else if (differentiation_method == "modify_material_matrix")
        {
            mDifferentiationMethod = 2;
            mNormalize = false;
        }
        else if (differentiation_method == "chain_rule")
        {
            mDifferentiationMethod = 3;
            mNormalize = false;
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

        if (rPseudoQuantityVariable == PSEUDO_MOMENT)
        {
            if (mDifferentiationMethod == 1)
            {
                const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);
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
                const Variable<double>& r_design_variable = KratosComponents<Variable<double>>::Get(mDesignVariableName);
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
            }
            // MFusseder TODO add also modification of the torsional inertia!
        }
        else
            KRATOS_ERROR << "CalculatePseudoQuantityByModificationOfMaterialMatrix not available for given variable!" << std::endl;

        rElement.CalculateOnIntegrationPoints(rQuantityVariable, rOutput, rCurrentProcessInfo);

        // Give element original properties back
        rElement.SetProperties(p_global_properties);

        KRATOS_CATCH("");
    }

};  // namespace Kratos.

