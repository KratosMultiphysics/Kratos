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
#include "utilities/openmp_utils.h"
#include "direct_sensitivity_local_stress_response_function.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{

    /// Constructor.
    DirectSensitivityLocalStressResponseFunction::DirectSensitivityLocalStressResponseFunction(ModelPart& rModelPart, 
                            Parameters ResponseSettings)
      :  DirectSensitivityResponseFunction(rModelPart, ResponseSettings)
    {
        KRATOS_TRY;   

        // Tell the traced forces 
        SizeType num_traced_forces = ResponseSettings["stress_types"]["forces"].size();
        mTracedForcesVector.resize(num_traced_forces);
        
        for ( IndexType i = 0; i < num_traced_forces; i++ )        
            mTracedForcesVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_types"]["forces"][i].GetString());
        
        // Tell the traced moments
        SizeType num_traced_moments = ResponseSettings["stress_types"]["moments"].size();
        mTracedMomentsVector.resize(num_traced_moments);
        
        for ( IndexType i = 0; i < num_traced_moments; i++ )        
            mTracedMomentsVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_types"]["moments"][i].GetString());
        
        SizeType num_traced_stresses = num_traced_forces + num_traced_moments;

        if (mTracedStressesVector.size() != num_traced_stresses)
            mTracedStressesVector.resize(num_traced_stresses);

        for ( IndexType i = 0; i < num_traced_forces; i++ )
        {
            mTracedStressesVector[i] = mTracedForcesVector[i];
            for ( IndexType j = 0; j < num_traced_moments; j++ )
            mTracedStressesVector[j+num_traced_forces] = mTracedMomentsVector[j];
        }
                      
        // Get info how and where to treat the stress
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );
        
        if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node)
        {
            mIdOfLocation = ResponseSettings["stress_location"].GetInt();
            KRATOS_ERROR_IF(mIdOfLocation < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
        }
        
        KRATOS_CATCH("");
    }


    DirectSensitivityLocalStressResponseFunction::~DirectSensitivityLocalStressResponseFunction(){}

    void DirectSensitivityLocalStressResponseFunction::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }    


    // Derivative of the response function "local stress" w.r.t the displacement
    void DirectSensitivityLocalStressResponseFunction::CalculateGradient(Element& rDirectElement,
                                   const Matrix& rLHS, 
                                   Matrix& rResponseGradientMatrix, 
                                   const ProcessInfo& rProcessInfo)
    {
        // Define sizes
        const SizeType num_traced_stresses = mTracedStressesVector.size();
        const SizeType num_dofs = rLHS.size1();
               
        ProcessInfo process_info = rProcessInfo;
        
        // Sizing of the Output Matrix
        if ( rResponseGradientMatrix.size1() != num_traced_stresses || rResponseGradientMatrix.size2() != num_dofs )
            rResponseGradientMatrix.resize( num_traced_stresses, num_dofs);                  
        
        for ( IndexType i = 0; i < num_traced_stresses; i++ )
        {  
            // Define working variable
            Vector rResponseGradientVector;
            TracedStressType traced_stress_type = mTracedStressesVector[i]; 

            // change the stress type for which the sensitivity will be calculated
            rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(traced_stress_type) );

            this->CalculateElementContributionToGradient(rDirectElement, rLHS, rResponseGradientVector, process_info);

            // Assembling of the gradient vectors for different traced stress types into a Matrix that can be given to postprocess
            for (IndexType j = 0; j < num_dofs; j++)
                rResponseGradientMatrix(i, j) = rResponseGradientVector[j];
        }
    }


     void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToGradient(Element& rDirectElement,
                                    const Matrix& rLHS,
                                    Vector& rResponseGradientVector,
                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        Matrix stress_displacement_derivative;

        if(mStressTreatment == StressTreatment::Mean)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradientVector);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_displacement_derivative, rResponseGradientVector);
        }
    
        KRATOS_ERROR_IF(rResponseGradientVector.size() != rLHS.size1())
            << "Size of stress displacement derivative does not fit!" << std::endl;

        rResponseGradientVector *= (-1);

        KRATOS_CATCH("");
    }
     
    void DirectSensitivityLocalStressResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable, 
                                    Vector& rSensitivityGradient, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY;

        // Define sizes
        const SizeType num_traced_stresses = mTracedStressesVector.size();
                
        // Sizing of rSensitivityGradient
        if ( rSensitivityGradient.size() != num_traced_stresses )
            rSensitivityGradient.resize(num_traced_stresses);

        if( rDirectElement.Id() == rDesignVariable.GetTracedElementId() )
        {
            for ( IndexType i = 0; i < num_traced_stresses; i++ )
            {   
                // Define working variable
                Vector sensitivity_gradient;

                // change the stress type for which the sensitivity will be calculated
                rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressesVector[i]) );
   
                const std::string variable_name = rDesignVariable.GetDesignVariableName();
                ProcessInfo process_info = rProcessInfo;
                
                this->CalculateElementContributionToPartialSensitivity(rDirectElement, variable_name, 
                                                                    sensitivity_gradient, process_info);

                rSensitivityGradient[i] = sensitivity_gradient[0];
            }
        }
        else
        {
            if (rSensitivityGradient.size() != 0)
                rSensitivityGradient.resize(0, false);
        }

        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                                    const std::string& rVariableName,
                                    Vector& rSensitivityGradient,
                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rDirectElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        Matrix stress_design_variable_derivative;

        if(mStressTreatment == StressTreatment::Mean)
        {
            rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_design_variable_derivative, rSensitivityGradient);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_design_variable_derivative, rSensitivityGradient);
        }
        
        KRATOS_ERROR_IF(rSensitivityGradient.size() != 1)
             << "Size of partial stress design variable does not fit!" << std::endl;

        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();
        double stress_derivative_value = 0.0;

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
        {
            for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
                stress_derivative_value += rStressDerivativesMatrix(deriv_it, stress_it);

            stress_derivative_value /= num_of_stress_positions;

            rResponseGradient[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::ExtractGaussPointStressDerivative(const Matrix& rStressDerivativesMatrix, Vector& rResponseGradient)
    {
        KRATOS_TRY;

        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        if(rResponseGradient.size() != num_of_derivatives_per_stress)
            rResponseGradient.resize(num_of_derivatives_per_stress, false);

        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= mIdOfLocation ) <<
                "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                                num_of_stress_positions  << "!"<< std::endl;

        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rResponseGradient[deriv_it] = rStressDerivativesMatrix(deriv_it, (mIdOfLocation-1));

        KRATOS_CATCH("");
    }

};


