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
#include "direct_sensitivity_nodal_displacement_response_function.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"



namespace Kratos
{
    /// Constructor.
    DirectSensitivityNodalDisplacementResponseFunction::DirectSensitivityNodalDisplacementResponseFunction(ModelPart& rModelPart, 
                            Parameters ResponseSettings)
      :  DirectSensitivityResponseFunction(rModelPart, ResponseSettings)
    {
        KRATOS_TRY;

        // Get list of traced displacements 
        std::vector<std::string> TracedDisplacementsVector;
        const SizeType num_traced_displacements = ResponseSettings["compute_on_node"]["dofs"]["displacements"].size();
        
        if (TracedDisplacementsVector.size() != num_traced_displacements)
            TracedDisplacementsVector.resize(num_traced_displacements);
        
        for ( IndexType i = 0; i < num_traced_displacements; i++ )        
            TracedDisplacementsVector[i] = ResponseSettings["compute_on_node"]["dofs"]["displacements"][i].GetString();
        
        // Get list of traced rotations
        std::vector<std::string> TracedRotationsVector;
        const SizeType num_traced_rotations = ResponseSettings["compute_on_node"]["dofs"]["rotations"].size();
        
        if (TracedRotationsVector.size() != num_traced_rotations)
            TracedRotationsVector.resize(num_traced_rotations);
        
        for ( IndexType i = 0; i < num_traced_rotations; i++ )        
            TracedRotationsVector[i] = ResponseSettings["compute_on_node"]["dofs"]["rotations"][i].GetString();
        
        // Get list of all traced dofs
        const SizeType num_traced_dofs = num_traced_displacements + num_traced_rotations;

        if (mTracedDofsVector.size() != num_traced_dofs)
            mTracedDofsVector.resize(num_traced_dofs);

        for ( IndexType i = 0; i < num_traced_displacements; i++ )
        {
            mTracedDofsVector[i] = TracedDisplacementsVector[i];
            for ( IndexType j = 0; j < num_traced_rotations; j++ )
                mTracedDofsVector[j+num_traced_displacements] = TracedRotationsVector[j];
        }

        // Get the locations for which the sensitivity should be computed
        const SizeType num_traced_nodes = ResponseSettings["output_definition"]["node_location"].size();
        if (mIdOfLocationVector.size() != num_traced_nodes)    
            mIdOfLocationVector.resize(num_traced_nodes, false);
        for ( IndexType i = 0; i < num_traced_nodes; i++ )
        {
            mIdOfLocationVector[i] = ResponseSettings["output_definition"]["node_location"][i].GetInt();
            KRATOS_ERROR_IF(mIdOfLocationVector[i] < 1) << "Chose a 'node_location' > 0. Specified 'stress_location': " 
                << mIdOfLocationVector[i] << std::endl;
        }

        KRATOS_CATCH("");
    }

    
    // Destructor    
    DirectSensitivityNodalDisplacementResponseFunction::~DirectSensitivityNodalDisplacementResponseFunction(){}

    
    void DirectSensitivityNodalDisplacementResponseFunction::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }
    

    // Derivative of the response function "nodal displacement" w.r.t the displacement
    void DirectSensitivityNodalDisplacementResponseFunction::CalculateGradient(Element& rDirectElement,
                                   const Matrix& rLHS, 
                                   Matrix& rResponseGradient, 
                                   const ProcessInfo& rProcessInfo)
    {   
        KRATOS_TRY;

        // Define the sizes        
        DofsVectorType dofs_of_element;    
        ProcessInfo process_info = rProcessInfo;    
        rDirectElement.GetDofList(dofs_of_element, process_info);        
        const SizeType num_dofs = dofs_of_element.size();
        const SizeType num_traced_nodes = mIdOfLocationVector.size();
        const SizeType num_traced_dofs = mTracedDofsVector.size();
        const SizeType total_size = num_traced_nodes * num_traced_dofs;
        
        // Sizing and clearing of the response gradient
        if ( rResponseGradient.size1() != total_size || rResponseGradient.size2() != num_dofs ) 
            rResponseGradient.resize(total_size, num_dofs, false);            

        rResponseGradient.clear();

        // Compute rResponseGradient      
        for(IndexType k = 0; k < num_traced_nodes; ++k)
        {
            IndexType index = k * num_traced_dofs;            
            for(IndexType j = 0; j < num_traced_dofs; ++j)
                for(IndexType i = 0; i < num_dofs; ++i)
                {
                    const VariableComponentType& r_traced_dof = 
                        KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Get(std::string("ADJOINT_") + mTracedDofsVector[j]);
                    
                    if(rDirectElement.GetGeometry()[mIdOfLocationVector[k]-1].pGetDof(r_traced_dof) == dofs_of_element[i])
                        rResponseGradient(index + j, i) = 1; 
                }
        }

        KRATOS_CATCH("");                  
    }


    // Derivative of the response function "nodal displacement" w.r.t the design parameter
    void DirectSensitivityNodalDisplacementResponseFunction::CalculatePartialSensitivity(Element& rDirectElement, 
                                    DirectSensitivityVariable& rDesignVariable, 
                                    Matrix& rSensitivityGradient, 
                                    const ProcessInfo& rProcessInfo)
    {   
        if (rSensitivityGradient.size1() != 0 || rSensitivityGradient.size2() != 0 ) 
            rSensitivityGradient.resize(0, 0, false);
        
        rSensitivityGradient.clear(); 
    }


    int DirectSensitivityNodalDisplacementResponseFunction::GetNumberOfOutputPositions()
    {
        const SizeType num_traced_nodes  = mIdOfLocationVector.size();
        return num_traced_nodes;
    }
};


