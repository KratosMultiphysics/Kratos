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
#include "direct_sensitivity_nodal_displacement_response_function.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"




namespace Kratos
{
    /// Constructor.
    DirectSensitivityNodalDisplacementResponseFunction::DirectSensitivityNodalDisplacementResponseFunction(ModelPart& rModelPart, 
                            Parameters ResponseSettings, std::string& ResponseVariableName)
      :  DirectSensitivityResponseFunction(rModelPart, ResponseSettings, ResponseVariableName)
    {
        KRATOS_TRY;
        
        KRATOS_CATCH("");
    }

    
    // Destructor    
    DirectSensitivityNodalDisplacementResponseFunction::~DirectSensitivityNodalDisplacementResponseFunction(){}

    
    void DirectSensitivityNodalDisplacementResponseFunction::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    void DirectSensitivityNodalDisplacementResponseFunction::CalculateGradient(Node<3>& rNode,                            
                                    Variable<array_1d<double, 3>> const& rResponseVariable,
                                    std::vector<array_1d<double, 3>>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        
        // Size rOutput
        rOutput.resize(3);        

        VectorMath::SetToZero(rOutput);
        
        // Compute derivative vector
        for (IndexType dir_it = 0; dir_it < 3; ++dir_it)
            for(IndexType dof_it = 0; dof_it < 3; ++dof_it)
                if (dir_it == dof_it)
                    rOutput[dof_it][dir_it] = 1;

        KRATOS_CATCH("");         
    }
    

    // Derivative of the response function "nodal displacement" w.r.t the design parameter
    void DirectSensitivityNodalDisplacementResponseFunction::CalculatePartialSensitivity(Node<3>& rNode, 
                                    DirectSensitivityVariable& rDesignVariable,
                                    Variable<array_1d<double, 3>> const& rStressVariable, 
                                    array_1d<double, 3>& rOutput, 
                                    const ProcessInfo& rProcessInfo)
    {   
        KRATOS_TRY;
        
        if (rOutput.size() != 0) 
            rOutput.resize(0);
        
        VectorMath::SetToZero(rOutput);

        KRATOS_CATCH("");
    }


    std::string DirectSensitivityNodalDisplacementResponseFunction::GetEvaluationFlag()
    {
        std::string flag = "on_node";
        return flag;
    } 

    
};


