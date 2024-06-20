// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Talhah Ansari, Suneth Warnakalusuriya
//

#pragma once

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/finite_difference_utility.h"
#include "system_identification_application_variables.h"

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

class MacResponseFunctionUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MacResponseFunctionUtility
    KRATOS_CLASS_POINTER_DEFINITION(MacResponseFunctionUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MacResponseFunctionUtility(ModelPart& model_part, Parameters sensitivity_settings)
    : mrModelPart(model_part)
    {
        
        CheckSettingsForGradientAnalysis(sensitivity_settings);
        
        DetermineTracedEigenfrequencies(sensitivity_settings);
        
        if(AreSeveralEigenfrequenciesTraced())
        {
            KRATOS_ERROR_IF_NOT(sensitivity_settings.Has("weighting_method")) << "Several eigenfrequencies are traced but no weighting method specified in the parameters!" << std::endl;
            if(sensitivity_settings["weighting_method"].GetString() == "linear_scaling")
                CalculateLinearWeights(sensitivity_settings);
            else
                KRATOS_ERROR  << "Specified weighting method for eigenfrequencies is not implemented! Available weighting methods are: 'linear_scaling'." << std::endl;
        }
        else
            UseDefaultWeight();

        

        mDesignVariable = sensitivity_settings["element_data_value_sensitivity_variables"][0].GetString();
        
    }

    /// Destructor.
    virtual ~MacResponseFunctionUtility()
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

        double resp_function_value = 0.0;
        KRATOS_INFO("MacResponseFunctionUtility") << "CalculateValue:" << std::endl
        << "    #    ComputedEigenfrequency [Hz]    MeasuredEigenfrequency [Hz]    MAC    weighting factor" <<std::endl;
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        {
            const double eigencomputedfrequency = std::sqrt(GetComputedEigenvalue(mTracedEigenfrequencyIds[i])) / (2*Globals::Pi);
            const double eigenmeasuredfrequency = std::sqrt(GetMeasuredEigenvalue(mTracedEigenfrequencyIds[i])) / (2*Globals::Pi);
            
            const Vector eigencomputedvector = GetComputedEigenvector(mTracedEigenfrequencyIds[i]);
            const Vector eigenmeasuredvector = GetMeasuredEigenvector(mTracedEigenfrequencyIds[i]);
            double num = inner_prod(eigencomputedvector, eigenmeasuredvector);
            double den1 =  inner_prod(eigencomputedvector, eigencomputedvector);
            double den2 = inner_prod(eigenmeasuredvector, eigenmeasuredvector);
            double MAC = ( (num*num) / (den1 * den2 ));

            resp_function_value += mWeightingFactors[i] * std::pow((1 - MAC), 2) * 0.5;
            std::stringstream msg;
            msg << std::setw(5) << mTracedEigenfrequencyIds[i]
                << std::setw(23) << eigencomputedfrequency
                << std::setw(23) << eigenmeasuredfrequency
                << std::setw(25) << MAC
                << std::setw(16) << mWeightingFactors[i]<< std::endl;
            KRATOS_INFO("") << msg.str();
        }

        return resp_function_value;
    }

    // --------------------------------------------------------------------------
    void CalculateGradient()
    {
        CheckIfAllNecessaryEigenvaluesAreComputed();
        VariableUtils().SetNonHistoricalVariableToZero(YOUNG_MODULUS_SENSITIVITY, mrModelPart.Elements());
        PerformSemiAnalyticSensitivityAnalysis();
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
        return "MacResponseFunctionUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "MacResponseFunctionUtility";
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
    void CheckSettingsForGradientAnalysis(Parameters rResponseSettings)
    {
        const std::string gradient_mode = rResponseSettings["gradient_mode"].GetString();

        if (gradient_mode == "semi_analytic")
            mDelta = rResponseSettings["step_size"].GetDouble();
        else
            KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;
    }

    // --------------------------------------------------------------------------
    void DetermineTracedEigenfrequencies(Parameters rResponseSettings)
    {
        mTracedEigenfrequencyIds.resize(rResponseSettings["traced_eigenfrequencies"].size(),false);
        for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
            mTracedEigenfrequencyIds[i] = rResponseSettings["traced_eigenfrequencies"][i].GetInt();
    }

    // --------------------------------------------------------------------------
    bool AreSeveralEigenfrequenciesTraced()
    {
        return (mTracedEigenfrequencyIds.size()>1);
    }

    // --------------------------------------------------------------------------
    void CalculateLinearWeights(Parameters rResponseSettings)
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
            for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
                mWeightingFactors[i] /= test_sum_weight_facs;

            KRATOS_INFO("\n> MacResponseFunctionUtility") << "The sum of the chosen weighting factors is unequal one. A corresponding scaling process was exected!" << std::endl;
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
        const int num_of_computed_eigenvalues = (mrModelPart.GetValue(COMPUTED_EIGENVALUE_VECTOR)).size();
        const int max_required_eigenvalue = *(std::max_element(mTracedEigenfrequencyIds.begin(), mTracedEigenfrequencyIds.end()));
        KRATOS_ERROR_IF(max_required_eigenvalue > num_of_computed_eigenvalues) << "The following Eigenfrequency shall be traced but their corresponding eigenvalue was not computed by the eigenvalue analysies: " << max_required_eigenvalue << std::endl;
    }

    // --------------------------------------------------------------------------
    void PerformSemiAnalyticSensitivityAnalysis()
    {
        const ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

        const auto& rDesignVariable = KratosComponents<Variable<double>>::Get(mDesignVariable);

        // Predetermine all necessary eigenvalues and prefactors for gradient calculation
        const std::size_t num_of_traced_eigenfrequencies = mTracedEigenfrequencyIds.size();
        Vector traced_eigenvalues(num_of_traced_eigenfrequencies,0.0);
        Vector gradient_prefactors(num_of_traced_eigenfrequencies,0.0);
        const int dummy_size = GetComputedEigenvector(1).size();

        std::vector<Vector> traced_eigenvectors(num_of_traced_eigenfrequencies,Vector(dummy_size));
        for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
        {
            traced_eigenvalues[i] = GetComputedEigenvalue(mTracedEigenfrequencyIds[i]);
            gradient_prefactors[i] = 1.0 / (4.0 * Globals::Pi * std::sqrt(traced_eigenvalues[i]));
            traced_eigenvectors[i] = GetComputedEigenvector(mTracedEigenfrequencyIds[i]);            
        }
        
        std::vector<Vector> measured_eigenvectors(10,Vector(dummy_size));
        for(std::size_t i = 0; i < 10; i++)
        {  
            measured_eigenvectors[i] = GetMeasuredEigenvector(i+1);            
        }

        // std::cout << "computed measured values" << std::endl;
        // Vector eigen_vector_i = GetComputedEigenvector(1);
        // eigen_vector_i.clear();
        // Vector eigen_vector_j = Vector(eigen_vector_i.size());
        // eigen_vector_j.clear();
        Vector eigen_vector_derivative = Vector(dummy_size);
        eigen_vector_derivative.clear();

        // Element-wise computation of gradients
        for (auto& elem_i : mrModelPart.Elements())
        {
            Matrix LHS;
            Matrix mass_matrix;
            elem_i.CalculateLeftHandSide(LHS, CurrentProcessInfo);
            elem_i.CalculateMassMatrix(mass_matrix, CurrentProcessInfo);

            const std::size_t num_dofs_element = mass_matrix.size1();
            //const std::size_t domain_size = CurrentProcessInfo.GetValue(DOMAIN_SIZE);

            Matrix aux_matrix = Matrix(num_dofs_element,num_dofs_element);
            Vector aux_vector = Vector(num_dofs_element);

            // Predetermine the necessary eigenvectors
            std::vector<Vector> eigenvectors_of_element(num_of_traced_eigenfrequencies,Vector(0));
            for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++){
                DetermineEigenvectorOfElement(elem_i, mTracedEigenfrequencyIds[i], eigenvectors_of_element[i], CurrentProcessInfo);
            }
        
            // Computation of derivative of state equation w.r.t. design variable
            double gradient_contribution = 0.0;
            Matrix derived_LHS = Matrix(num_dofs_element,num_dofs_element);
            Matrix derived_mass_matrix = Matrix(num_dofs_element,num_dofs_element);
            
            FiniteDifferenceUtility::CalculateLeftHandSideDerivative(elem_i, LHS, rDesignVariable, mDelta, derived_LHS, CurrentProcessInfo);
            FiniteDifferenceUtility::CalculateMassMatrixDerivative(elem_i, mass_matrix, rDesignVariable, mDelta, derived_mass_matrix, CurrentProcessInfo);

            eigen_vector_derivative.clear();

            for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
            {
                aux_matrix.clear();
                aux_vector.clear();

                noalias(aux_matrix) = derived_LHS - derived_mass_matrix * traced_eigenvalues[i];
                noalias(aux_vector) = prod(aux_matrix , eigenvectors_of_element[i]);

                const double eigen_value_i = traced_eigenvalues[i];
                //eigen_vector_i.clear();
                //eigen_vector_i = GetComputedEigenvector(i+1);
                
                for (std::size_t j = 0; j < 10; j++){
                    if (j == i){
                        const double prefactor1 = (-1/(2*eigen_value_i)) * inner_prod(eigenvectors_of_element[i], prod(derived_LHS, eigenvectors_of_element[i]));
                        noalias(eigen_vector_derivative) += traced_eigenvectors[i] * prefactor1;                       
                    }
                    else{
                        double eigen_value_j = GetComputedEigenvalue(j+1);
                        //eigen_vector_j.clear();
                        //eigen_vector_j = GetComputedEigenvector(j+1);

                        double phi_j_phi_j = 0.0;
                        const int num_node = elem_i.GetGeometry().size();
                        const int num_node_dofs = num_dofs_element / num_node;
                        // for (const auto& r_node : mrModelPart.Nodes()) {
                        //     //auto& r_node = rModelPart.GetNode(i_node+1);
                        //     const Matrix& r_nodal_eigenvector_matrix = r_node.GetValue(COMPUTED_EIGENVECTOR_MATRIX);
                        //     //phi_j_phi_j += inner_prod(r_nodal_eigenvector_matrix(j,-1), r_nodal_eigenvector_matrix(j,-1));
                        //     for (int k = 0; k < num_node_dofs; k++) {
                        //         phi_j_phi_j += r_nodal_eigenvector_matrix(j, k) * r_nodal_eigenvector_matrix(j, k); 
                        //     }
                        // }
                        phi_j_phi_j += inner_prod(measured_eigenvectors[j], measured_eigenvectors[j]);
                        double prefactor2 = (1/(eigen_value_i - eigen_value_j)) * phi_j_phi_j;

                        int current_node = 0;
                        for (auto& r_node_i : elem_i.GetGeometry())
                        {                        
                            int node_index = r_node_i.Id() - 1;
                            for (int dof_index = 0; dof_index < num_node_dofs; ++dof_index)
                            {
                                eigen_vector_derivative[node_index * num_node_dofs + dof_index] += prefactor2 * aux_vector[current_node * num_node_dofs + dof_index]; //correct to add at right place of global system
                            }
                            current_node +=1;
                        }
                    }
                }

                const Vector& c_der = eigen_vector_derivative;
                const Vector& c = traced_eigenvectors[i];
                const Vector& m = measured_eigenvectors[i];

                const double cc = inner_prod(c,c); // scalar product of computed eigenvector
                const double mm = inner_prod(m,m); //scalar product of measured eigenvector
                const double cm = inner_prod(c,m); //scalar product of computed and measured eigenvector
                const double mc_der = inner_prod(m,c_der); //scalar product of measured eigenvector and computed eigenvector derivative
                const double cc_der = inner_prod(c,c_der); //scalar product of computed eigenvector and computed eigenvector derivative

                const double dMAC = 2* ((cm * cc * mc_der) - ( cm * cm * cc_der)) / (cc * cc * mm); //derivative of the Modal Assurance Criteria

                gradient_contribution += -1* dMAC * mWeightingFactors[i];
                // std::cout << "gradient is " << gradient_contribution << std::endl;
            }

            elem_i.GetValue(YOUNG_MODULUS_SENSITIVITY) += gradient_contribution;
        }

    }

    // --------------------------------------------------------------------------
    double GetComputedEigenvalue(const int eigenfrequency_id)
    {
        return (mrModelPart.GetValue(COMPUTED_EIGENVALUE_VECTOR))[eigenfrequency_id-1];
    }

    double GetMeasuredEigenvalue(const int eigenfrequency_id)
    {
        return (mrModelPart.GetValue(MEASURED_EIGENVALUE_VECTOR))[eigenfrequency_id-1];
    }

    Vector GetComputedEigenvector(const int eigenfrequency_id)
    {
        const int number_of_nodes = mrModelPart.NumberOfNodes();

        const Node& rnode = mrModelPart.GetNode(1);
        const Node::DofsContainerType& rnode_dofs = rnode.GetDofs();
        const int num_dofs_per_node = rnode_dofs.size();
        const int global_size = number_of_nodes * num_dofs_per_node;

        Vector eigenvector;
        SetVectorToZero(eigenvector, global_size);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {

            const auto& r_node = mrModelPart.GetNode(i_node+1);
            const Matrix nodal_eigenvector = r_node.GetValue(COMPUTED_EIGENVECTOR_MATRIX);

            for (int j = 0; j < num_dofs_per_node; j++) {
                eigenvector[i_node * num_dofs_per_node + j] = nodal_eigenvector((eigenfrequency_id-1), j);
            }
        }

        return eigenvector;
    }

    Vector GetMeasuredEigenvector(const int eigenfrequency_id)
    {
        const int number_of_nodes = mrModelPart.NumberOfNodes();

        const Node& rnode = mrModelPart.GetNode(1);
        const Node::DofsContainerType& rnode_dofs = rnode.GetDofs();
        const int num_dofs_per_node = rnode_dofs.size();
        const int global_size = number_of_nodes * num_dofs_per_node;

        Vector eigenvector;
        SetVectorToZero(eigenvector, global_size);

        for (int i_node = 0; i_node < number_of_nodes; ++i_node) {

            const auto& r_node = mrModelPart.GetNode(i_node+1);
            const Matrix nodal_eigenvector = r_node.GetValue(MEASURED_EIGENVECTOR_MATRIX);

            for (int j = 0; j < num_dofs_per_node; j++) {
                eigenvector[i_node * num_dofs_per_node + j] = nodal_eigenvector((eigenfrequency_id-1), j);
            }
        }

        return eigenvector;
    }

    // --------------------------------------------------------------------------
    void DetermineEigenvectorOfElement(ModelPart::ElementType& rElement, const int eigenfrequency_id, Vector& rEigenvectorOfElement, const ProcessInfo& CurrentProcessInfo)
    {
        std::vector<std::size_t> eq_ids;
        rElement.EquationIdVector(eq_ids, CurrentProcessInfo);

        if (rEigenvectorOfElement.size() != eq_ids.size())
        {
            rEigenvectorOfElement.resize(eq_ids.size(), false);
        }

        // sort the values of the eigenvector into the rEigenvectorOfElement according to the dof ordering at the element
        for (auto& r_node_i : rElement.GetGeometry())
        {
            const auto& r_node_dofs = r_node_i.GetDofs();

            const Matrix& rNodeEigenvectors = r_node_i.GetValue(COMPUTED_EIGENVECTOR_MATRIX);
            for (std::size_t dof_index = 0; dof_index < r_node_dofs.size(); dof_index++)
            {
                const auto& current_dof = *(std::begin(r_node_dofs) + dof_index);
                const std::size_t dof_index_at_element = std::distance(eq_ids.begin(), std::find(eq_ids.begin(), eq_ids.end(), current_dof->EquationId()));
                rEigenvectorOfElement(dof_index_at_element) = rNodeEigenvectors((eigenfrequency_id-1), dof_index);
            }
        }
    }


    void SetVectorToZero(
        Vector &rVector,
        const IndexType Size)
    {
        if (rVector.size() != Size)
        {
            rVector.resize(Size, false);
        }
        rVector.clear();
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
    double mDelta;
    std::vector<int> mTracedEigenfrequencyIds;
    std::vector<double> mWeightingFactors;
    std::string mDesignVariable;

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
    //      MacResponseFunctionUtility& operator=(MacResponseFunctionUtility const& rOther);

    /// Copy constructor.
    //      MacResponseFunctionUtility(MacResponseFunctionUtility const& rOther);

    ///@}

}; // Class MacResponseFunctionUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.
