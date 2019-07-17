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


// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "direct_sensitivity_local_stress_response_function.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "derivative_builder.h"


namespace Kratos
{

    /// Constructor.
    DirectSensitivityLocalStressResponseFunction::DirectSensitivityLocalStressResponseFunction(ModelPart& rModelPart, 
                            Parameters ResponseSettings, std::string& ResponseVariableName)
      :  DirectSensitivityResponseFunction(rModelPart, ResponseSettings, ResponseVariableName)
    {
        KRATOS_TRY;   
                            
        // Get info how and where to treat the stress
        mStressTreatment = 
            StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["stress_treatment"].GetString() );
        
        if( mStressTreatment == StressTreatment::GaussPoint )
        {   
            const SizeType num_traced_gp = ResponseSettings["stress_location"].size();
            if (mIdOfLocationVector.size() != num_traced_gp)    
                mIdOfLocationVector.resize(num_traced_gp, false);
            for ( IndexType i = 0; i < num_traced_gp; i++ )
            {
                mIdOfLocationVector[i] = ResponseSettings["stress_location"][i].GetInt();
                KRATOS_ERROR_IF(mIdOfLocationVector[i] < 1) << "Chose a 'stress_location' > 0. Specified 'stress_location': " 
                    << mIdOfLocationVector[i] << std::endl;
            }
        }
        else if( mStressTreatment == StressTreatment::Mean )
        {
            if (mIdOfLocationVector.size() != 1)
                mIdOfLocationVector.resize(1, false);
            
            mIdOfLocationVector.clear();
        }
        
        KRATOS_CATCH("");
    }

    // Destructor
    DirectSensitivityLocalStressResponseFunction::~DirectSensitivityLocalStressResponseFunction(){}

    
    void DirectSensitivityLocalStressResponseFunction::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }    


    // Derivative of the response function "local stress" w.r.t the displacement
    void DirectSensitivityLocalStressResponseFunction::CalculateGradient(Element& rDirectElement,                            
                            Variable<array_1d<double, 3>> const& rStressVariable,
                            std::vector<std::vector<array_1d<double, 3>>>& rOutput, 
                            const ProcessInfo& rProcessInfo)
    {
        this->CalculateElementContributionToGradient(rDirectElement, rStressVariable, rOutput, rProcessInfo);
    }

    void DirectSensitivityLocalStressResponseFunction::CalculateGradient(Element& rDirectElement,                            
                            Variable<Matrix> const& rStressVariable,
                            std::vector<std::vector<Matrix>>& rOutput, 
                            const ProcessInfo& rProcessInfo)
    {
        this->CalculateElementContributionToGradient(rDirectElement, rStressVariable, rOutput, rProcessInfo);
    }

    template <typename TDataType>
    void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToGradient(Element& rDirectElement,
                                    Variable<TDataType> const& rStressVariable,
                                    std::vector<std::vector<TDataType>>& rOutput,
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        
        // To get the number of Dofs
        DofsVectorType dofs_of_element;    
        ProcessInfo process_info = rProcessInfo;    
        rDirectElement.GetDofList(dofs_of_element, process_info);

        // Define sizes
        const SizeType num_dofs = dofs_of_element.size();

        // Define working variables
        std::vector<std::vector<TDataType>> response_gradient;
        VectorMath::SetToZero(response_gradient);

        // Compute the derivative of the local stress wrt. the displacement.
        DerivativeBuilder::ComputeStressDisplacementDerivative(rDirectElement, rStressVariable, response_gradient, rProcessInfo);        
        
        // Size rOutput
        if (rOutput.size() != num_dofs)
            rOutput.resize(num_dofs);

        //  Extract the searched gauss point values        
        for (IndexType dof_it = 0; dof_it < num_dofs; ++dof_it)
        {   
            if( mStressTreatment == StressTreatment::Mean )
                this->ExtractMeanStressDerivative(response_gradient[dof_it], rOutput[dof_it]);
            else if( mStressTreatment == StressTreatment::GaussPoint )
                this->ExtractGaussPointStressDerivative(response_gradient[dof_it], rOutput[dof_it]);
        }
        
        KRATOS_CATCH("");
    }


    // Derivative of the response function "local stress" w.r.t the design variable 
    void DirectSensitivityLocalStressResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<array_1d<double, 3>> const& rStressVariable, 
                                    std::vector<array_1d<double, 3>>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY; 

        this->PartialSensitivityBuilder(rDirectElement, rDesignVariable, rStressVariable, rOutput, rProcessInfo);

        KRATOS_CATCH("");
    }

    void DirectSensitivityLocalStressResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<Matrix> const& rStressVariable, 
                                    std::vector<Matrix>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY; 

        this->PartialSensitivityBuilder(rDirectElement, rDesignVariable, rStressVariable, rOutput, rProcessInfo);

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void DirectSensitivityLocalStressResponseFunction::PartialSensitivityBuilder(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<TDataType> const& rStressVariable, 
                                    std::vector<TDataType>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY;
        
        // Define working variables
        std::vector<std::vector<TDataType>> sensitivity_gradient;
        std::vector<TDataType> extracted_sensitivity_gradient;

        rOutput.resize(0);        

        // Compute Partial Sensitivity
        std::vector<unsigned int> traced_elem_ids = rDesignVariable.GetTracedElementId();
        if (std::find(traced_elem_ids.begin(), traced_elem_ids.end(),rDirectElement.Id() ) != traced_elem_ids.end())
        {            
            CalculateElementContributionToPartialSensitivity(rDirectElement, rDesignVariable, rStressVariable, 
                                            sensitivity_gradient, rProcessInfo); 
            
            // Extract the searched gauss point value from sensitivity_gradient
            if(mStressTreatment == StressTreatment::Mean)
                this->ExtractMeanStressDerivative(sensitivity_gradient[0], rOutput);
            else if(mStressTreatment == StressTreatment::GaussPoint)
                this->ExtractGaussPointStressDerivative(sensitivity_gradient[0], rOutput);
        } 

        KRATOS_CATCH("");
    }

    template <typename TDataType>
    void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<TDataType> const& rStressVariable,
                                    std::vector<std::vector<TDataType>>& rOutput,
                                    const ProcessInfo& rProcessInfo)
    {
        // Get the name of the design variable 
        const std::string design_variable_name = rDesignVariable.GetDesignVariableName(); 
        rDirectElement.SetValue(DESIGN_VARIABLE_NAME, design_variable_name);
        
        
        // Get perturbation size  
        double delta = rDesignVariable.GetPerturbationSize(); 
        rDirectElement.SetValue(PERTURBATION_SIZE, delta);

        // Compute the derivative of the local stress wrt. the design variable. 
        DerivativeBuilder::ComputeStressDesignVariableDerivative(rDirectElement, rStressVariable, rDesignVariable, rOutput, rProcessInfo);        
        
    }

    template <typename TDataType>
    void DirectSensitivityLocalStressResponseFunction::ExtractMeanStressDerivative(
                                    const std::vector<TDataType>& rStressDerivativesMatrix,
                                    std::vector<TDataType>& rOutput)
    {
        KRATOS_TRY;
                
        // Define sizes
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size();
        
        // Define working variables
        TDataType stress_derivative_value;
        VectorMath::SizeVectorComponents(stress_derivative_value);
        VectorMath::SetToZero(stress_derivative_value);        
    
        // Sizing of rOutput
        if(rOutput.size() != 1)
            rOutput.resize(1);

        VectorMath::SizeVectorComponents(rOutput);        
        VectorMath::SetToZero(rOutput);
        
        // Compute mean value of all gauss points          
        for(IndexType gp_it = 0; gp_it < num_of_stress_positions; ++gp_it)            
            VectorMath::Addition(stress_derivative_value, rStressDerivativesMatrix[gp_it]);                
            
        stress_derivative_value /= num_of_stress_positions;
                        
        VectorMath::Addition(rOutput[0], stress_derivative_value);
        VectorMath::SetToZero(stress_derivative_value);   

        KRATOS_CATCH("");
    }
    
    template <typename TDataType>
    void DirectSensitivityLocalStressResponseFunction::ExtractGaussPointStressDerivative( 
                                    const std::vector<TDataType>& rStressDerivativeVector,
                                    std::vector<TDataType>& rOutput)
    {
        KRATOS_TRY;
        

        // Define sizes
        const SizeType num_gp = rStressDerivativeVector.size();
        const SizeType num_traced_gp = mIdOfLocationVector.size();
        
        // Sizing of rOutput
        if(rOutput.size() != num_traced_gp)
            rOutput.resize(num_traced_gp);
        
        VectorMath::SizeVectorComponents(rOutput);
        VectorMath::SetToZero(rOutput);

        for (IndexType i = 0; i < num_traced_gp; ++i)
        {    
            // Check if choosen gauss point is available
            if (num_gp > 1)
            {
                KRATOS_ERROR_IF_NOT(num_gp >= mIdOfLocationVector[i] ) <<
                    "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " << num_gp  << "!"<< std::endl;

                // Extract the values for the choosen gauss point from the derivative matrix
                VectorMath::Addition( rOutput[i], rStressDerivativeVector[mIdOfLocationVector[i]-1] );
            }
            else
                VectorMath::Addition( rOutput[i], rStressDerivativeVector[0] );
        }

        KRATOS_CATCH("");
    }

    
    std::string DirectSensitivityLocalStressResponseFunction::GetEvaluationFlag()
    {
        std::string flag = "on_gauss_point";
        return flag;
    } 
    
};


