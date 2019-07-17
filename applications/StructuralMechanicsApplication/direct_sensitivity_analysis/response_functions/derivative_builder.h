// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#ifndef DERIVATIVE_BUILDER_H
#define DERIVATIVE_BUILDER_H

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "structural_mechanics_application_variables.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/define.h"
#include "direct_sensitivity_analysis/math_functions/vector_math.h"
#include "direct_sensitivity_analysis/variable_utilities/direct_sensitivity_variable.h"
#include "utilities/compare_elements_and_conditions_utility.h"


namespace Kratos
{

/** \brief DerivativeBuilder
*
* This class computes the derivative of a certain response variable (e.g. MOMENT, FORCE etc.) derived by
* either the displacement or the design variable. It also assembles the result in a vector
* matching the data type of the response variable (e.g. array_1d<double, 3>, Matrix etc.) 
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DerivativeBuilder
{
public:

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    
    // Calculate the derivative of the response variable wrt. the displacement and assemble the results in a vector
    template <typename TDataType>
    static void ComputeStressDisplacementDerivative( Element& rDirectElement, 
                                    Variable<TDataType> const& rResponseVariable,
                                    std::vector<std::vector<TDataType>>& rOutput, 
                                    const ProcessInfo rCurrentProcessInfo)
    {         
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();
        std::vector<TDataType> dummy_vector;
        std::vector<std::string> traced_stresses_list;
        GetTracedStressesList( rResponseVariable, traced_stresses_list );

        // To get the number of Dofs
        DofsVectorType dofs_of_element;    
        ProcessInfo process_info = rCurrentProcessInfo;    
        rDirectElement.GetDofList(dofs_of_element, process_info);

        // To get the number of integration points
        rDirectElement.CalculateOnIntegrationPoints(rResponseVariable, dummy_vector, rCurrentProcessInfo);        
        const SizeType num_gp = dummy_vector.size();

        // Get the name of the current element
        std::string name_elem;
        CompareElementsAndConditionsUtility::GetRegisteredName(rDirectElement, name_elem); 

        // Size rOutput        
        rOutput.resize( dofs_of_element.size() );        
                
        for (IndexType i = 0; i < rOutput.size(); ++i)
        {
            if (num_gp > 0)
                rOutput[i].resize( num_gp );
            else
            {
                if (name_elem == "AdjointFiniteDifferenceCrBeamElementLinear3D2N" || name_elem == "AdjointFiniteDifferencingShellThinElement3D3N")
                    rOutput[i].resize(3);     
                else
                    rOutput[i].resize(1);
            }
        }

        VectorMath::SizeVectorComponents(rOutput);        

        // Set rOutput to zero
        VectorMath::SetToZero(rOutput);
        
        if (num_gp == 1 && rResponseVariable == FORCE)
            DeriveTrussForceDisplacementDerivative(rDirectElement, rResponseVariable ,rOutput, rCurrentProcessInfo);    
        else if (num_gp > 1)
        {                     
            for (IndexType dir_it = 0; dir_it < traced_stresses_list.size(); ++dir_it)
            { 
                TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType(traced_stresses_list[dir_it]);  
                rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
                            
                rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, DerivativeMatrix, rCurrentProcessInfo);                
                
                AssembleDerivationVector(rOutput, DerivativeMatrix, dir_it);
               
                DerivativeMatrix.clear();
            }            
        }                 
    }

    // Calculate the derivative of the response variable wrt. the design variable and assemble the results in a vector
    template <typename TDataType>
    static void ComputeStressDesignVariableDerivative( Element& rDirectElement, 
                                    Variable<TDataType> const& rResponseVariable,
                                    DirectSensitivityVariable& rDesignVariable,
                                    std::vector<std::vector<TDataType>>& rOutput, 
                                    const ProcessInfo rCurrentProcessInfo)
    { 
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();
        Matrix extracted_derivative_matrix;
        extracted_derivative_matrix.clear();
        std::vector<TDataType> dummy_vector;
        std::vector<std::string> traced_stresses_list;
        GetTracedStressesList( rResponseVariable, traced_stresses_list );
        
        // To get the number of integration points
        rDirectElement.CalculateOnIntegrationPoints(rResponseVariable, dummy_vector, rCurrentProcessInfo);        
        const SizeType num_gp = dummy_vector.size();

        // Get the name of the current element
        std::string name_elem;
        CompareElementsAndConditionsUtility::GetRegisteredName(rDirectElement, name_elem);
        
        // Size rOutput        
        rOutput.resize(1);        
                
        for (IndexType i = 0; i < rOutput.size(); ++i)
        {
            if (num_gp > 0)
                rOutput[i].resize( num_gp );
            else
            {
                if (name_elem == "AdjointFiniteDifferenceCrBeamElementLinear3D2N" || name_elem == "AdjointFiniteDifferencingShellThinElement3D3N")
                    rOutput[i].resize(3);     
                else
                    rOutput[i].resize(1);
            }
        }

        VectorMath::SizeVectorComponents(rOutput);

        // Set rOutput to zero
        VectorMath::SetToZero(rOutput);

        if (num_gp == 1 && rResponseVariable == FORCE)
            DeriveTrussForceDesignVariableDerivative(rDirectElement, rResponseVariable, rDesignVariable, rOutput, rCurrentProcessInfo);    
        else if (num_gp > 1)
        {
            for (IndexType dir_it = 0; dir_it < traced_stresses_list.size(); ++dir_it)
            { 
                TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType(traced_stresses_list[dir_it]);  
                rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
                            
                rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, DerivativeMatrix, rCurrentProcessInfo);
                
                if ( rDesignVariable.GetDesignVariableType() != "nodal_coordinates_type")
                    AssembleDerivationVector(rOutput, DerivativeMatrix, dir_it);    
                else    
                {                    
                    rDesignVariable.ExtractDataFromDerivativeMatrix(rDirectElement, extracted_derivative_matrix, DerivativeMatrix);
                    AssembleDerivationVector(rOutput, extracted_derivative_matrix, dir_it);                    
                }
                DerivativeMatrix.clear();
            }
        }         
    } 
    
private: 

    template <typename TDataType>
    static void GetTracedStressesList(Variable<TDataType> const& rResponseVariable,
                                    std::vector<std::string>& rTracedStresses)
    {  
        if (rResponseVariable == MOMENT)
            rTracedStresses = { std::string("MX"), std::string("MY"), std::string("MZ") };

        else if (rResponseVariable == FORCE)
            rTracedStresses = { std::string("FX"), std::string("FY"), std::string("FZ") };

        else if (rResponseVariable == SHELL_MOMENT_GLOBAL)
            rTracedStresses = { std::string("MXX"), std::string("MXY"), std::string("MXZ"), std::string("MYX"), 
                std::string("MYY"), std::string("MYZ"), std::string("MZX"), std::string("MZY"), std::string("MZZ") };

        else if (rResponseVariable == SHELL_FORCE_GLOBAL)
            rTracedStresses = { std::string("FXX"), std::string("FXY"), std::string("FXZ"), std::string("FYX"), 
                std::string("FYY"), std::string("FYZ"), std::string("FZX"), std::string("FZY"), std::string("FZZ") };

        else
            KRATOS_ERROR << "ResponseVariable " << rResponseVariable.Name() << " is not yet implemented in DerivativeBuilder" << std::endl;                 
    }    

    template <typename TDataType>
    static void DeriveTrussForceDisplacementDerivative(Element& rDirectElement,
                                Variable<TDataType> const& rResponseVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    { 
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();

        // Calculate Fx of the truss
        TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType("FX");
        rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
           
        rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, DerivativeMatrix, rCurrentProcessInfo);

        // Assemble results for the response function FORCE
        AssembleDerivationVector(rOutput, DerivativeMatrix, 0);               
        DerivativeMatrix.clear();
        AssembleDerivationVector(rOutput, DerivativeMatrix, 1);
        AssembleDerivationVector(rOutput, DerivativeMatrix, 2);                     
    }

    template <typename TDataType>
    static void DeriveTrussForceDesignVariableDerivative(Element& rDirectElement,
                                Variable<TDataType> const& rResponseVariable,
                                DirectSensitivityVariable& rDesignVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    { 
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();
        Matrix extracted_derivative_matrix;
        extracted_derivative_matrix.clear();

        // Calculate Fx of the truss
        TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType("FX");  
        rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
           
        rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, DerivativeMatrix, rCurrentProcessInfo);

        // Assemble results for the response function FORCE
        if ( rDesignVariable.GetDesignVariableType() != "nodal_coordinates_type")
        {
            AssembleDerivationVector(rOutput, DerivativeMatrix, 0);               
            DerivativeMatrix.clear();
            AssembleDerivationVector(rOutput, DerivativeMatrix, 1);
            AssembleDerivationVector(rOutput, DerivativeMatrix, 2); 
        }
        else    
        {                    
            rDesignVariable.ExtractDataFromDerivativeMatrix(rDirectElement, extracted_derivative_matrix, DerivativeMatrix);
            AssembleDerivationVector(rOutput, extracted_derivative_matrix, 0);
            extracted_derivative_matrix.clear();
            AssembleDerivationVector(rOutput, extracted_derivative_matrix, 1);
            AssembleDerivationVector(rOutput, extracted_derivative_matrix, 2);                     
        }                     
    }
    

    template <typename TDataType>
    static void AssembleDerivationVector(std::vector<std::vector<TDataType>>& rOutput,
                                    Matrix& DerivativeMatrix, unsigned int dir_it)
    {   
        // Define sizes        
        const SizeType num_derivatives = rOutput.size();
        const SizeType num_gp = rOutput[0].size();        
        
        // Assemble the results from DerivativeMatrix (results for 1 of the traced stresses) in the output vector
        for (IndexType deriv_it = 0; deriv_it < num_derivatives; ++deriv_it)
            for(IndexType gp_it = 0; gp_it < num_gp; ++gp_it)
                AssembleComponentOfDerivationVector(rOutput[deriv_it][gp_it], DerivativeMatrix(deriv_it, gp_it), dir_it);
    }


    static void AssembleComponentOfDerivationVector(array_1d<double,3>& rDerivationVectorComponent, const double& rValue, unsigned int dir_it)
    {   
        rDerivationVectorComponent[dir_it] = rValue;                    
    }    
    
    
    static void AssembleComponentOfDerivationVector(Matrix& rDerivationVectorComponent, const double& rValue, unsigned int dir_it)
    {          
        if (dir_it == 0)            
            rDerivationVectorComponent(0, 0) = rValue;
        if (dir_it == 1)            
            rDerivationVectorComponent(1, 0) = rValue;
        if (dir_it == 2)            
            rDerivationVectorComponent(2, 0) = rValue;
        if (dir_it == 3)            
            rDerivationVectorComponent(0, 1) = rValue;
        if (dir_it == 4)            
            rDerivationVectorComponent(1, 1) = rValue;
        if (dir_it == 5)            
            rDerivationVectorComponent(2, 1) = rValue;
        if (dir_it == 6)            
            rDerivationVectorComponent(0, 2) = rValue;
        if (dir_it == 7)            
            rDerivationVectorComponent(1, 2) = rValue;
        if (dir_it == 8)            
            rDerivationVectorComponent(2, 2) = rValue;                    
    }
    
    
}; // class DerivativeBuilder


} /* namespace Kratos.*/

#endif /* DERIVATIVE_BUILDER_H defined */