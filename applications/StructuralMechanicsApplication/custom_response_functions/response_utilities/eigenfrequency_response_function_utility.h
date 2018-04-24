// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Fusseder Martin, Armin Geiser, Daniel Baumgärtner
//

#ifndef EIGENFREQUENCY_RESPONSE_FUNCTION_UTILITY_H
#define EIGENFREQUENCY_RESPONSE_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.

 */

//template<class TDenseSpace>

class EigenfrequencyResponseFunctionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EigenfrequencyResponseFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(EigenfrequencyResponseFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EigenfrequencyResponseFunctionUtility(ModelPart& model_part, Parameters response_settings)
    : mrModelPart(model_part)
    {
        DetermineGradientMode(response_settings);
        DetermineTracedEigenfrequencies(response_settings);

        if(AreSeveralEigenfrequenciesTraced())
        {
            if(!response_settings.Has("weighting_method"))
                KRATOS_ERROR  << "Several eigenfrequencies are traced but no weighting method specified in the parameters!" << std::endl;
            if(response_settings["weighting_method"].GetString().compare("linear_scaling") == 0)
                CalculateLinearWeights(response_settings);
            else
                KRATOS_ERROR  << "Specified weighting method for eigenfrequencies is not implemented! Available weighting methods are: 'linear_scaling'." << std::endl;
        }
        else
            UseDefaultWeight();
    }

    /// Destructor.
    virtual ~EigenfrequencyResponseFunctionUtility()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void Initialize()
    {
        //not needed because only semi-analytical sensitivity analysis is implemented yet
    }

    // --------------------------------------------------------------------------
    double CalculateValue()
    {
        CheckIfAllNecessaryEigenvaluesAreComputed();

        double m_resp_function_value = 0.0;

        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
            m_resp_function_value += mWeightingFactors[i] * std::sqrt(GetEigenvalue(mTracedEigenfrequencyIds[i])) / (2*Globals::Pi);

        return m_resp_function_value;
    }

    // --------------------------------------------------------------------------
    void CalculateGradient()
    {
        CheckIfAllNecessaryEigenvaluesAreComputed();

        VariableUtils().SetToZero_VectorVar(SHAPE_SENSITIVITY, mrModelPart.Nodes());

        switch (mGradientMode)
        {
            // Semi-analytic sensitivities
            case 1:
            {
                CalculateCompleteResponseDerivativeByFiniteDifferencing();
                break;
            }
        }
    }

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "EigenfrequencyResponseFunctionUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "EigenfrequencyResponseFunctionUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // ==============================================================================
    void DetermineGradientMode(Parameters& rResponseSettings)
    {
        const std::string gradient_mode = rResponseSettings["gradient_mode"].GetString();

        if (gradient_mode.compare("semi_analytic") == 0)
        {
            mGradientMode = 1;
            mDelta = rResponseSettings["step_size"].GetDouble();
        }
        else
            KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;
    }

    // --------------------------------------------------------------------------
    void DetermineTracedEigenfrequencies(Parameters& rResponseSettings)
    {
        if(rResponseSettings["traced_eigenfrequency"].IsArray())
        {
            mTracedEigenfrequencyIds.resize(rResponseSettings["traced_eigenfrequency"].size(),false);
            for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
                mTracedEigenfrequencyIds[i] = rResponseSettings["traced_eigenfrequency"][i].GetInt();
        }
        else
        {
            mTracedEigenfrequencyIds.resize(1,false);
            mTracedEigenfrequencyIds[0] = rResponseSettings["traced_eigenfrequency"].GetInt();
        }
    }

    // --------------------------------------------------------------------------
    bool AreSeveralEigenfrequenciesTraced()
    {
        if(mTracedEigenfrequencyIds.size()>1)
            return true;
        else
            return false;
    }

    // --------------------------------------------------------------------------
    void CalculateLinearWeights(Parameters& rResponseSettings)
    {
        KRATOS_ERROR_IF_NOT(rResponseSettings.Has("weighting_factors")) << "No weighting factors defined for given eigenfrequency response!" << std::endl;
        KRATOS_ERROR_IF_NOT(rResponseSettings["weighting_factors"].size() == mTracedEigenfrequencyIds.size()) << "The number of chosen eigenvalues does not fit to the number of weighting factors!" << std::endl;

        mWeightingFactors.resize(rResponseSettings["weighting_factors"].size(),false);

        // Read weighting factors and sum them up
        double test_sum_weight_facs = 0.0;
        for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
        {
            mWeightingFactors[i] = rResponseSettings["weighting_factors"][i].GetDouble();
            test_sum_weight_facs += mWeightingFactors[i];
        }

        // Check the weighting factors and modify them for the case that their sum is unequal to one
        if(std::abs(test_sum_weight_facs - 1.0) > 1e-12)
        {
            for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
                mWeightingFactors[i] /= test_sum_weight_facs;

            std::cout << "> The sum of the chosen weighting factors is unequal to one. A scaling process was executed for them!" << std::endl;
        }
    }

    // --------------------------------------------------------------------------
    void UseDefaultWeight()
    {
        mWeightingFactors.resize(1,false);
        mWeightingFactors[0] = 1.0;
    }

    // --------------------------------------------------------------------------
    void CheckIfAllNecessaryEigenvaluesAreComputed()
    {
        const int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR]).size();
        const int max_required_eigenfrequency = *std::max_element(mTracedEigenfrequencyIds.begin(), mTracedEigenfrequencyIds.end());
        KRATOS_ERROR_IF(max_required_eigenfrequency > num_of_computed_eigenvalues) << "The following Eigenfrequency shall be traced but was not computed by the eigenvalue analysies: " << max_required_eigenfrequency << std::endl;
    }

    // --------------------------------------------------------------------------
    void CalculateCompleteResponseDerivativeByFiniteDifferencing()
    {
        ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

        // Predetermine all necessary eigenvalues and prefactors for gradient calculation
        const unsigned int num_of_traced_eigenfrequencies = mTracedEigenfrequencyIds.size();
        Vector traced_eigenvalues(num_of_traced_eigenfrequencies,0.0);
        Vector gradient_prefactors(num_of_traced_eigenfrequencies,0.0);
        for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
        {
            traced_eigenvalues[i] = GetEigenvalue(mTracedEigenfrequencyIds[i]);
            gradient_prefactors[i] = 1 / (2*Globals::Pi * std::sqrt(traced_eigenvalues[i]));
        }

        // Element-wise computation of gradients
        for (auto& elem_i : mrModelPart.Elements())
        {
            Matrix mass_matrix_org;
            Matrix LHS_org;
            Vector dummy;
            Matrix aux_matrix = Matrix(0,0);
            Vector aux_vector = Vector(0);
            elem_i.CalculateMassMatrix(mass_matrix_org, CurrentProcessInfo);
            elem_i.CalculateLocalSystem(LHS_org, dummy ,CurrentProcessInfo);

            // Predetermine the necessary eigenvectors
            int num_dofs_element = mass_matrix_org.size1();
            std::vector<Vector> eigenvectors_of_element(num_of_traced_eigenfrequencies,Vector(0));
            for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
            {
                eigenvectors_of_element[i].resize(num_dofs_element,false);
                eigenvectors_of_element[i] = GetEigenvectorOfElement(elem_i, mTracedEigenfrequencyIds[i], num_dofs_element);
            }

            // Semi-analytic computation of partial derivative of state equation w.r.t. node coordinates
            for (auto& node_i : elem_i.GetGeometry())
            {
                Vector gradient_contribution(3, 0.0);
                Matrix perturbed_LHS = Matrix(0,0);
                Matrix perturbed_mass_matrix = Matrix(0,0);

                // Derivative of response w.r.t. x-coord ------------------------
                node_i.X0() += mDelta;

                elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
                perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

                elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
                perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

                for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
                {
                    aux_matrix = perturbed_LHS;
                    aux_matrix -= (perturbed_mass_matrix * traced_eigenvalues[i]);
                    aux_vector = prod(aux_matrix , eigenvectors_of_element[i]);
                    gradient_contribution[0] += gradient_prefactors[i] * inner_prod(eigenvectors_of_element[i] , aux_vector) * mWeightingFactors[i];

                    aux_matrix = Matrix(0,0);
                    aux_vector = Vector(0);
                }

                node_i.X0() -= mDelta;

                // End derivative of response w.r.t. x-coord --------------------

                // Reset pertubed RHS and mass matrix
                perturbed_LHS = Matrix(0,0);
                perturbed_mass_matrix = Matrix(0,0);

                // Derivative of response w.r.t. y-coord ------------------------
                node_i.Y0() += mDelta;

                elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
                perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

                elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
                perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

                // Loop over eigenvalues
                for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
                {
                    aux_matrix = perturbed_LHS;
                    aux_matrix -= (perturbed_mass_matrix * traced_eigenvalues[i]);
                    aux_vector = prod(aux_matrix , eigenvectors_of_element[i]);
                    gradient_contribution[1] += gradient_prefactors[i] * inner_prod(eigenvectors_of_element[i] , aux_vector) * mWeightingFactors[i];

                    aux_matrix = Matrix(0,0);
                    aux_vector = Vector(0);;
                }

                node_i.Y0() -= mDelta;
                // End derivative of response w.r.t. y-coord --------------------

                // Reset pertubed RHS and mass matrix
                perturbed_LHS = Matrix(0,0);
                perturbed_mass_matrix= Matrix(0,0);

                // Derivative of response w.r.t. z-coord ------------------------
                node_i.Z0() += mDelta;

                elem_i.CalculateMassMatrix(perturbed_mass_matrix, CurrentProcessInfo);
                perturbed_mass_matrix = (perturbed_mass_matrix - mass_matrix_org) / mDelta;

                elem_i.CalculateLocalSystem(perturbed_LHS, dummy ,CurrentProcessInfo);
                perturbed_LHS = (perturbed_LHS - LHS_org) / mDelta;

                // Loop over eigenvalues
                for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
                {
                    aux_matrix = perturbed_LHS;
                    aux_matrix -= (perturbed_mass_matrix * traced_eigenvalues[i]);
                    aux_vector = prod(aux_matrix , eigenvectors_of_element[i]);
                    gradient_contribution[2] += gradient_prefactors[i] * inner_prod(eigenvectors_of_element[i] , aux_vector) * mWeightingFactors[i];

                    aux_matrix = Matrix(0,0);
                    aux_vector = Vector(0);
                }

                node_i.Z0() -= mDelta;
                // End derivative of response w.r.t. z-coord --------------------

                // Assemble sensitivity to node
                noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;

            }// End loop over nodes of element

        }// End loop over elements
    }

    // --------------------------------------------------------------------------
    double GetEigenvalue(int eigenfrequency_id)
    {
        const int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR]).size();
        return (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR])[eigenfrequency_id-1];
    }

    // --------------------------------------------------------------------------
    Vector GetEigenvectorOfElement(ModelPart::ElementType& traced_element, int eigenfrequency_id, int size_of_eigenvector)
    {
        Vector eigenvector_of_element(size_of_eigenvector,0.0);

        // Get eigenvector of element
        int k = 0;
        const int NumNodeDofs = size_of_eigenvector/traced_element.GetGeometry().size();
        for (auto& node_i : traced_element.GetGeometry())
        {
            Matrix& rNodeEigenvectors = node_i.GetValue(EIGENVECTOR_MATRIX);

            for (int i = 0; i < NumNodeDofs; i++)
                eigenvector_of_element(i+NumNodeDofs*k) = rNodeEigenvectors((eigenfrequency_id-1),i);

            k++;
        }

        return eigenvector_of_element;
    }

    // ==============================================================================

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    unsigned int mGradientMode;
    double mDelta;
    std::vector<int> mTracedEigenfrequencyIds;
    std::vector<double> mWeightingFactors;

    ///@}
///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //      EigenfrequencyResponseFunctionUtility& operator=(EigenfrequencyResponseFunctionUtility const& rOther);

    /// Copy constructor.
    //      EigenfrequencyResponseFunctionUtility(EigenfrequencyResponseFunctionUtility const& rOther);

    ///@}

}; // Class EigenfrequencyResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // EIGENFREQUENCY_RESPONSE_FUNCTION_UTILITY_H
