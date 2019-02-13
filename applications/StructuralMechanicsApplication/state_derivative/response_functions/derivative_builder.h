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

           
    template <typename TDataType>
    static void ComputeDerivative(Element& rDirectElement, 
                                Variable<TDataType> const& rResponseVariable,
                                const Variable<Matrix >& rDerivativeVariable,
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
        
        // Define seizes
        const SizeType num_gp = dummy_vector.size();
        SizeType num_derivatives;
        if(rDerivativeVariable == STRESS_DISP_DERIV_ON_GP)
            num_derivatives = dofs_of_element.size();
        else if (rDerivativeVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
            num_derivatives = 1;
        else
        {
            num_derivatives = 0;
            KRATOS_ERROR << "Response Variable can only get derived by the initial state results or by a design variable!" << std::endl;
        }
           
        // Size rOutput
        rOutput.resize(num_derivatives);
        for (IndexType i = 0; i < rOutput.size(); ++i)
            rOutput[i].resize(num_gp);

        // Calculate the derivative of the response variable and assemble the results in a vector 
        if (rResponseVariable == MOMENT)
        {
            std::array<std::string, 3> moments = { std::string("MX"), std::string("MY"), std::string("MZ") };
            DeriveStressVariable(moments, rDirectElement, rResponseVariable, rDerivativeVariable, rOutput, rCurrentProcessInfo);
        }
    }
    
    
    
private:

    template <typename TDataType>
    static void DeriveStressVariable(std::array<std::string, 3>& rTracedStresses,
                                Element& rDirectElement, 
                                Variable<TDataType> const& rResponseVariable,
                                const Variable<Matrix >& rDerivativeVariable,
                                std::vector<std::vector<TDataType>>& rOutput, 
                                const ProcessInfo rCurrentProcessInfo)
    {
        // Define working variables
        Matrix DerivativeMatrix;
        DerivativeMatrix.clear();

        // Define sizes
        const SizeType num_derivatives = rOutput.size();
        const SizeType num_gp = rOutput[0].size();

        for (IndexType dir_it = 0; dir_it < 3; ++dir_it)
        {   
            TracedStressType traced_stress = StressResponseDefinitions::ConvertStringToTracedStressType(rTracedStresses[dir_it]);  
            rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress));

            if(rDerivativeVariable == STRESS_DISP_DERIV_ON_GP)
                rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, DerivativeMatrix, rCurrentProcessInfo);
            else if (rDerivativeVariable == STRESS_DESIGN_DERIVATIVE_ON_GP)
                rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, DerivativeMatrix, rCurrentProcessInfo);

            for (IndexType deriv_it = 0; deriv_it < num_derivatives; ++deriv_it)
                for(IndexType gp_it = 0; gp_it < num_gp; ++gp_it)
                    rOutput[deriv_it][gp_it][dir_it] = DerivativeMatrix(deriv_it, gp_it);
                
            DerivativeMatrix.clear();
        }        
    }
    
    
}; // class DerivativeBuilder


} /* namespace Kratos.*/

#endif /* DERIVATIVE_BUILDER_H defined */