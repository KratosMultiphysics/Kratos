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
#include "state_derivative/output_utilities/output_utility.h"
#include "state_derivative/math_functions/vector_math.h"
#include "utilities/compare_elements_and_conditions_utility.h"


namespace Kratos
{

/** \brief DerivativeBuilder
*
* This class computes the derivative of a certain response variable (MOMENT etc.) derived by
* either the displacement or the design variable. It also assembles the result in a vector
* matching the data type of the response variable (array_1d<double, 3>, Matrix etc.) 
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DerivativeBuilder
{
public:

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>> VariableComponentType;

    
    // Calculate the derivative of the response variable and assemble the results in a vector        
    template <typename TDataType>
    static void ComputeDerivative(const std::string& DerivativeFlag,
                                Element& rDirectElement, 
                                Variable<TDataType> const& rResponseVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    {  
        if (rResponseVariable == MOMENT)
        {
            std::vector<std::string> moments = { std::string("MX"), std::string("MY"), std::string("MZ") };
            DeriveStressVariable(DerivativeFlag, rDirectElement, moments , rResponseVariable, rOutput, rCurrentProcessInfo);
        }
        if (rResponseVariable == FORCE)
        {
            std::vector<std::string> forces = { std::string("FX"), std::string("FY"), std::string("FZ") };
            DeriveStressVariable(DerivativeFlag, rDirectElement, forces , rResponseVariable, rOutput, rCurrentProcessInfo);
        }
        if (rResponseVariable == SHELL_MOMENT_GLOBAL)
        {
            std::vector<std::string> shell_moments = { std::string("MXX"), std::string("MXY"), std::string("MXZ"), std::string("MYX"), 
                                    std::string("MYY"), std::string("MYZ"), std::string("MZX"), std::string("MZY"), std::string("MZZ") };
            DeriveStressVariable(DerivativeFlag, rDirectElement, shell_moments , rResponseVariable, rOutput, rCurrentProcessInfo);
        }  
        if (rResponseVariable == SHELL_FORCE_GLOBAL)
        {
            std::vector<std::string> shell_forces = { std::string("FXX"), std::string("FXY"), std::string("FXZ"), std::string("FYX"), 
                                    std::string("FYY"), std::string("FYZ"), std::string("FZX"), std::string("FZY"), std::string("FZZ") };
            DeriveStressVariable(DerivativeFlag, rDirectElement, shell_forces , rResponseVariable, rOutput, rCurrentProcessInfo);
        }                       
    }
    
    
    
private:

    template <typename TDataType>
    static void DeriveStressVariable(const std::string& DerivativeFlag,
                                Element& rDirectElement,                                
                                const std::vector<std::string>& rTracedStresses,
                                Variable<TDataType> const& rResponseVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    {   
        
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();
        std::vector<TDataType> dummy_vector;

        // To get the number of Dofs
        DofsVectorType dofs_of_element;    
        ProcessInfo process_info = rCurrentProcessInfo;    
        rDirectElement.GetDofList(dofs_of_element, process_info);

        // To get the number of integration points
        rDirectElement.CalculateOnIntegrationPoints(rResponseVariable, dummy_vector, rCurrentProcessInfo);        

        // Size rOutput
        if(DerivativeFlag == "DISPLACEMENT_DERIVATIVE")
            rOutput.resize( dofs_of_element.size() );
        else if (DerivativeFlag == "DESIGN_VARIABLE_DERIVATIVE")
            rOutput.resize(1);
        else
            KRATOS_ERROR << "Response Variable can only get derived by the initial state results or by a design variable!" << std::endl;
                
        for (IndexType i = 0; i < rOutput.size(); ++i)
            rOutput[i].resize( dummy_vector.size() );

        VectorMath::SizeVectorComponents(rOutput);

        // Set rOutput to zero
        VectorMath::SetToZero(rOutput);

        if (dummy_vector.size() == 1 && rResponseVariable == FORCE)
            DeriveTrussForce(DerivativeFlag, rDirectElement, rResponseVariable ,rOutput, rCurrentProcessInfo);    
        else
        {
            for (IndexType dir_it = 0; dir_it < rTracedStresses.size(); ++dir_it)
            {   
                TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType(rTracedStresses[dir_it]);  
                rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
            
                if(DerivativeFlag == "DISPLACEMENT_DERIVATIVE")
                    rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, DerivativeMatrix, rCurrentProcessInfo);
                else if (DerivativeFlag == "DESIGN_VARIABLE_DERIVATIVE")
                    rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, DerivativeMatrix, rCurrentProcessInfo);

                AssembleDerivationVector(rOutput, DerivativeMatrix, dir_it);
               
                DerivativeMatrix.clear();
            }
        }        
    }

    template <typename TDataType>
    static void DeriveTrussForce(const std::string& DerivativeFlag,
                                Element& rDirectElement,
                                Variable<TDataType> const& rResponseVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    { 
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();

        TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType("FX");  
        rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));            
            
        if(DerivativeFlag == "DISPLACEMENT_DERIVATIVE")
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, DerivativeMatrix, rCurrentProcessInfo);
        else if (DerivativeFlag == "DESIGN_VARIABLE_DERIVATIVE")
            rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, DerivativeMatrix, rCurrentProcessInfo);

        AssembleDerivationVector(rOutput, DerivativeMatrix, 0);
               
        DerivativeMatrix.clear();

        AssembleDerivationVector(rOutput, DerivativeMatrix, 1);

        AssembleDerivationVector(rOutput, DerivativeMatrix, 2);                       
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