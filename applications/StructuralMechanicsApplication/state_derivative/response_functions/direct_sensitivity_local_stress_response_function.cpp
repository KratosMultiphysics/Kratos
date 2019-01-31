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
        /*SizeType num_traced_forces = ResponseSettings["stress_type"]["forces"].size();
        mTracedForcesVector.resize(num_traced_forces);
        
        for (SizeType i = 0; i < num_traced_forces; i++)        
            mTracedForcesVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"]["forces"][i].GetString());
        
        // Tell the traced moments
        SizeType num_traced_moments = ResponseSettings["stress_type"]["moments"].size();
        mTracedMomentsVector.resize(num_traced_moments);
        
        for (SizeType i = 0; i < num_traced_moments; i++)        
            mTracedMomentsVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"]["moments"][i].GetString());*/
        

        // Tell  the stress type*/
        mTracedStressType = StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["stress_type"].GetString());
        
        // Get info how and where to treat the stress
        mStressTreatment = StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );
        
        if(mStressTreatment == StressTreatment::GaussPoint || mStressTreatment == StressTreatment::Node)
        {
            mIdOfLocation = ResponseSettings["stress_location"].GetInt();
            KRATOS_ERROR_IF(mIdOfLocation < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " << mIdOfLocation << std::endl;
        }
        
        // Get info how and where to treat the stress
       

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
                                   Vector& rResponseGradient, 
                                   const ProcessInfo& rProcessInfo)
    {
        Matrix stress_displacement_derivative;

        rDirectElement.SetValue(TRACED_STRESS_TYPE, static_cast<int>(mTracedStressType) );

        if(mStressTreatment == StressTreatment::Mean)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradient);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_displacement_derivative, rResponseGradient);
        }
        
        KRATOS_ERROR_IF(rResponseGradient.size() != rLHS.size1())
             << "Size of stress displacement derivative does not fit!" << std::endl;

        rResponseGradient *= (-1);
    }   
      

    void DirectSensitivityLocalStressResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable, 
                                    Vector& rSensitivityGradient, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY;

        if( rDirectElement.Id() == rDesignVariable.GetTracedElementId() )
        {   
            const std::string variable_name = rDesignVariable.GetDesignVariableName();
            ProcessInfo process_info = rProcessInfo;
            this->CalculateElementContributionToPartialSensitivity(rDirectElement, variable_name, 
                                                                    rSensitivityGradient, process_info);
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


