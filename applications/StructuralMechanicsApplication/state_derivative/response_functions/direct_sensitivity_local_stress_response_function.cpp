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

        // Get list of traced forces 
        std::vector<TracedStressType> TracedForcesVector;
        const SizeType num_traced_forces = ResponseSettings["compute_on_gp"]["stresses"]["forces"].size();
        
        if (TracedForcesVector.size() != num_traced_forces)
            TracedForcesVector.resize(num_traced_forces);
        
        for ( IndexType i = 0; i < num_traced_forces; i++ )        
            TracedForcesVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["compute_on_gp"]["stresses"]["forces"][i].GetString());
        
        // Get list of traced moments
        std::vector<TracedStressType> TracedMomentsVector;
        const SizeType num_traced_moments = ResponseSettings["compute_on_gp"]["stresses"]["moments"].size();
        
        if (TracedMomentsVector.size() != num_traced_moments)
            TracedMomentsVector.resize(num_traced_moments);
        
        for ( IndexType i = 0; i < num_traced_moments; i++ )        
            TracedMomentsVector[i] = 
                StressResponseDefinitions::ConvertStringToTracedStressType(ResponseSettings["compute_on_gp"]["stresses"]["moments"][i].GetString());
        
        // Get list of all traced stresses
        const SizeType num_traced_stresses = num_traced_forces + num_traced_moments;

        if (mTracedStressesVector.size() != num_traced_stresses)
            mTracedStressesVector.resize(num_traced_stresses);

        for ( IndexType i = 0; i < num_traced_forces; i++ )
        {
            mTracedStressesVector[i] = TracedForcesVector[i];
            for ( IndexType j = 0; j < num_traced_moments; j++ )
                mTracedStressesVector[j+num_traced_forces] = TracedMomentsVector[j];
        }
                      
        // Get info how and where to treat the stress
        mStressTreatment = 
            StressResponseDefinitions::ConvertStringToStressTreatment( ResponseSettings["output_definition"]["stress_treatment"].GetString() );
        
        if( mStressTreatment == StressTreatment::GaussPoint )
        {   
            const SizeType num_traced_gp = ResponseSettings["output_definition"]["stress_location"].size();
            if (mIdOfLocationVector.size() != num_traced_gp)    
                mIdOfLocationVector.resize(num_traced_gp, false);
            for ( IndexType i = 0; i < num_traced_gp; i++ )
            {
                mIdOfLocationVector[i] = ResponseSettings["output_definition"]["stress_location"][i].GetInt();
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
                                   const Matrix& rLHS, 
                                   Matrix& rResponseGradientMatrix, 
                                   const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // Define sizes
        const SizeType num_traced_gp = mIdOfLocationVector.size();
        const SizeType num_traced_stresses = mTracedStressesVector.size();
        const SizeType num_dofs = rLHS.size1();
        const SizeType total_size = num_traced_stresses * num_traced_gp;

        // Define working variables
        Vector rResponseGradientVector;
        ProcessInfo process_info = rProcessInfo;
        
        // Sizing of the Output Matrix
        if ( rResponseGradientMatrix.size1() != total_size || rResponseGradientMatrix.size2() != num_dofs )
            rResponseGradientMatrix.resize(total_size, num_dofs, false);                        
        
        // Compute rResponseGradientMatrix
        for ( IndexType k = 0; k < num_traced_gp; k++ )
        {
            IndexType index = k * num_traced_stresses;
            for ( IndexType i = 0; i < num_traced_stresses; i++ )
            {  
                // Define working variables
                TracedStressType traced_stress_type = mTracedStressesVector[i];
                unsigned int id_of_gauss_point = mIdOfLocationVector[k];
                
                // change the stress type for which the sensitivity will be calculated
                rDirectElement.SetValue( TRACED_STRESS_TYPE, static_cast<int>(traced_stress_type) );

                this->CalculateElementContributionToGradient(rDirectElement, rLHS, rResponseGradientVector, id_of_gauss_point, process_info);
                
                // Assembling of the gradient vectors for different traced stress types into a Matrix
                for (IndexType j = 0; j < num_dofs; j++)
                    rResponseGradientMatrix(index + i, j) = rResponseGradientVector[j];
                
                rResponseGradientVector.clear();
            }
        }

        KRATOS_CATCH("");
    }


     void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToGradient(Element& rDirectElement,
                                    const Matrix& rLHS,
                                    Vector& rResponseGradientVector,
                                    unsigned int& rIdOfGaussPoint,
                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // Define working variables
        Matrix stress_displacement_derivative;

        // Compute derivative matrix and extract the traced gauss point values
        if(mStressTreatment == StressTreatment::Mean)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_displacement_derivative, rResponseGradientVector);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rDirectElement.Calculate(STRESS_DISP_DERIV_ON_GP, stress_displacement_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_displacement_derivative, rIdOfGaussPoint, rResponseGradientVector);
        }

        // Check if rResponseGradientVector has the right size for the post-process
        KRATOS_ERROR_IF(rResponseGradientVector.size() != rLHS.size1())
            << "Size of stress displacement derivative does not fit!" << std::endl;

        KRATOS_CATCH("");
    }


    // Derivative of the response function "local stress" w.r.t the design variable 
    void DirectSensitivityLocalStressResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable, 
                                    Matrix& rSensitivityGradient, 
                                    const ProcessInfo& rProcessInfo)
    {        
        KRATOS_TRY;

        // Define sizes
        const SizeType num_traced_stresses = mTracedStressesVector.size();
        const SizeType num_traced_gp = mIdOfLocationVector.size();
                        
        // Sizing of rSensitivityGradient
        if ( rSensitivityGradient.size1() != num_traced_gp  || rSensitivityGradient.size2() != num_traced_stresses )
            rSensitivityGradient.resize( num_traced_gp, num_traced_stresses, false );

        // Define working variables
        const std::string variable_name = rDesignVariable.GetDesignVariableName();
        ProcessInfo process_info = rProcessInfo;        

        // Compute rSensitivityGradient
        if( rDirectElement.Id() == rDesignVariable.GetTracedElementId() )
            for ( IndexType j = 0; j < num_traced_gp; j++ )
                for ( IndexType i = 0; i < num_traced_stresses; i++ )
                {   
                    // Define working variables
                    Vector sensitivity_gradient;
                    unsigned int id_of_gauss_point = mIdOfLocationVector[j];

                    // Change the stress type for which the sensitivity will be calculated
                    rDirectElement.SetValue( TRACED_STRESS_TYPE, static_cast<int>(mTracedStressesVector[i]) );
                                       
                    this->CalculateElementContributionToPartialSensitivity(rDirectElement, variable_name, sensitivity_gradient,
                                                id_of_gauss_point, process_info);

                    rSensitivityGradient(j, i) = -sensitivity_gradient[0];
                }  
        else
            if (rSensitivityGradient.size1() != 0  || rSensitivityGradient.size2() != 0)
                rSensitivityGradient.resize(0, 0, false);
        
        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::CalculateElementContributionToPartialSensitivity(Element& rDirectElement,
                                    const std::string& rVariableName,
                                    Vector& rSensitivityGradientVector,
                                    unsigned int& rIdOfGaussPoint,
                                    ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        rDirectElement.SetValue(DESIGN_VARIABLE_NAME, rVariableName);

        // Define working variables
        Matrix stress_design_variable_derivative;

        // Compute derivative matrix and extract the traced gauss point values
        if(mStressTreatment == StressTreatment::Mean)
        {
            rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractMeanStressDerivative(stress_design_variable_derivative, rSensitivityGradientVector);
        }
        else if(mStressTreatment == StressTreatment::GaussPoint)
        {
            rDirectElement.Calculate(STRESS_DESIGN_DERIVATIVE_ON_GP, stress_design_variable_derivative, rProcessInfo);
            this->ExtractGaussPointStressDerivative(stress_design_variable_derivative, rIdOfGaussPoint, rSensitivityGradientVector);
        }
        
        // Check if rSensitivityGradientVector is a scalar value
        KRATOS_ERROR_IF(rSensitivityGradientVector.size() != 1)
             << "Size of partial stress design variable does not fit!" << std::endl;

        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::ExtractMeanStressDerivative(const Matrix& rStressDerivativesMatrix, 
                                    Vector& rOutput)
    {
        KRATOS_TRY;
        
        // Define sizes
        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();
        double stress_derivative_value = 0.0;

        // Sizing of rOutput
        if(rOutput.size() != num_of_derivatives_per_stress)
            rOutput.resize(num_of_derivatives_per_stress, false);

        // Compute mean value of all gauss points
        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
        {
            for(IndexType stress_it = 0; stress_it < num_of_stress_positions; ++stress_it)
                stress_derivative_value += rStressDerivativesMatrix(deriv_it, stress_it);

            stress_derivative_value /= num_of_stress_positions;

            rOutput[deriv_it] = stress_derivative_value;
            stress_derivative_value = 0.0;
        }

        KRATOS_CATCH("");
    }


    void DirectSensitivityLocalStressResponseFunction::ExtractGaussPointStressDerivative(const Matrix& rStressDerivativesMatrix, 
                                    unsigned int& rIdOfGaussPoint,
                                    Vector& rOutput)
    {
        KRATOS_TRY;

        // Define sizes
        const SizeType num_of_derivatives_per_stress = rStressDerivativesMatrix.size1();
        const SizeType num_of_stress_positions = rStressDerivativesMatrix.size2();

        // Sizing of rOutput
        if(rOutput.size() != num_of_derivatives_per_stress)
            rOutput.resize(num_of_derivatives_per_stress, false);

        // Check if choosen gauss point is available
        KRATOS_ERROR_IF_NOT(num_of_stress_positions >= rIdOfGaussPoint ) <<
                "Chosen Gauss-Point is not available. Chose 'stress_location' between 1 and " <<
                            num_of_stress_positions  << "!"<< std::endl;

        // Extract the values for the choosen gauss point from the derivative matrix
        for (IndexType deriv_it = 0 ; deriv_it < num_of_derivatives_per_stress; ++deriv_it)
            rOutput[deriv_it] = rStressDerivativesMatrix(deriv_it, (rIdOfGaussPoint-1));

        KRATOS_CATCH("");
    }


    int DirectSensitivityLocalStressResponseFunction::GetNumberOfOutputPositions()
    {
        const SizeType num_traced_gp  = mIdOfLocationVector.size();
        return num_traced_gp;
    }
};


