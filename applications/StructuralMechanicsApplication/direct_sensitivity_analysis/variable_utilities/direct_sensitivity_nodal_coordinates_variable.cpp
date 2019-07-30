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
#include "direct_sensitivity_nodal_coordinates_variable.h"
#include "processes/find_nodal_neighbours_process.h"


namespace Kratos
{
    // Constructor
    DirectSensitivityNodalCoordinatesVariable::DirectSensitivityNodalCoordinatesVariable(ModelPart& rModelPart, Parameters VariableSettings)
    : DirectSensitivityVariable(rModelPart, VariableSettings)
    {        

        FindNodalNeighboursProcess N_FINDER = FindNodalNeighboursProcess(mrModelPart, 10, 10);
		N_FINDER.Execute();

        // Get the Id of the elements that adjoin the traced nodes
        this->GetNeighboringElementIdVector(mNeighboringElementIdVector);

        // Get the Id of the condition that adjoin the traced nodes
        this->GetNeighboringConditionIdVector(mNeighboringConditionIdVector);

        // Get the direction in which the node coordinates should be perturbed
        mCoordinateDirection = VariableSettings["coordinate_direction"].GetString();       
        
    }

    // Destructor
    DirectSensitivityNodalCoordinatesVariable::~DirectSensitivityNodalCoordinatesVariable(){}


    void DirectSensitivityNodalCoordinatesVariable::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

   
    void DirectSensitivityNodalCoordinatesVariable::CalculatePseudoLoadVector(Element& rDirectElement, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        
        KRATOS_TRY;        
        
        // Define working variable        
        ProcessInfo process_info = rProcessInfo;
        Vector init_position;        
        Vector coord_value;        
        unsigned int coord_dir = GetCoordinateDirection();
        
        // Sizing and clearing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size1())
            rPseudoLoadVector.resize(rLHS.size1(), false);

        rPseudoLoadVector.clear();                

        if (std::find(mNeighboringElementIdVector.begin(), mNeighboringElementIdVector.end(), rDirectElement.Id() ) != mNeighboringElementIdVector.end())    
        {
            // Calculate right hand side vector before perturbation
            Vector RHS_unperturbed;
            rDirectElement.CalculateRightHandSide(RHS_unperturbed, process_info);
            
            const SizeType number_of_nodes = rDirectElement.GetGeometry().PointsNumber();

            init_position.resize(number_of_nodes);
            init_position.clear();
            coord_value.resize(number_of_nodes);
            coord_value.clear();

            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id) 
            {
                // Save unperturbed nodal coordinates
                init_position[node_id] = rDirectElement.GetGeometry()[node_id].GetInitialPosition()[coord_dir]; 
                coord_value[node_id] = rDirectElement.GetGeometry()[node_id].Coordinates()[coord_dir];

                if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), rDirectElement.GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() ) 
                {
                    // perturb the nodal coordinates
                    rDirectElement.GetGeometry()[node_id].GetInitialPosition()[coord_dir] += mDelta;
                    rDirectElement.GetGeometry()[node_id].Coordinates()[coord_dir] += mDelta;
                }
            }

            // Calculate right hand side vector after perturbation
            Vector RHS_perturbed;
            rDirectElement.CalculateRightHandSide(RHS_perturbed, process_info);
               
            // Calculate the pseudo load
            noalias(rPseudoLoadVector) = (RHS_perturbed - RHS_unperturbed ) / mDelta;

            // Set the nodal coordinates back to the unperturbed values
            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id)                   
                if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), rDirectElement.GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() ) 
                {                    
                    rDirectElement.GetGeometry()[node_id].GetInitialPosition()[coord_dir] = init_position[node_id];
                    rDirectElement.GetGeometry()[node_id].Coordinates()[coord_dir] = coord_value[node_id];
                }           
        }
        else
            rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size());      

        KRATOS_CATCH("");           
    }

    void DirectSensitivityNodalCoordinatesVariable::CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        
        // Define working variable        
        ProcessInfo process_info = rProcessInfo;
        Vector init_position;        
        Vector coord_value;        
        unsigned int coord_dir = GetCoordinateDirection();
        
        // Sizing and clearing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size1())
            rPseudoLoadVector.resize(rLHS.size1(), false);

        rPseudoLoadVector.clear();                

        if (std::find(mNeighboringConditionIdVector.begin(), mNeighboringConditionIdVector.end(), rDirectCondition.Id() ) != mNeighboringConditionIdVector.end())    
        {
            // Calculate right hand side vector before perturbation
            Vector RHS_unperturbed;
            rDirectCondition.CalculateRightHandSide(RHS_unperturbed, process_info);
            
            const SizeType number_of_nodes = rDirectCondition.GetGeometry().PointsNumber();

            init_position.resize(number_of_nodes);
            init_position.clear();
            coord_value.resize(number_of_nodes);
            coord_value.clear();

            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id) 
            {
                // Save unperturbed nodal coordinates
                init_position[node_id] = rDirectCondition.GetGeometry()[node_id].GetInitialPosition()[coord_dir]; 
                coord_value[node_id] = rDirectCondition.GetGeometry()[node_id].Coordinates()[coord_dir];

                if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), rDirectCondition.GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() ) 
                {
                    // perturb the nodal coordinates
                    rDirectCondition.GetGeometry()[node_id].GetInitialPosition()[coord_dir] += mDelta;
                    rDirectCondition.GetGeometry()[node_id].Coordinates()[coord_dir] += mDelta;
                }
            }

            // Calculate right hand side vector after perturbation
            Vector RHS_perturbed;
            rDirectCondition.CalculateRightHandSide(RHS_perturbed, process_info);
               
            // Calculate the pseudo load
            noalias(rPseudoLoadVector) = (RHS_perturbed - RHS_unperturbed ) / mDelta;

            // Set the nodal coordinates back to the unperturbed values
            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id)                   
                if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), rDirectCondition.GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() ) 
                {                    
                    rDirectCondition.GetGeometry()[node_id].GetInitialPosition()[coord_dir] = init_position[node_id];
                    rDirectCondition.GetGeometry()[node_id].Coordinates()[coord_dir] = coord_value[node_id];
                }           
        }
        else
            rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size());
        
        KRATOS_CATCH("");
    }

    void DirectSensitivityNodalCoordinatesVariable::ExtractDataFromDerivativeMatrix(Element& rDirectElement,
                                                    Matrix& rExtractedDerivativeMatrix,
                                                    const Matrix& rDerivativeMatrix)
    { 
        if ( rExtractedDerivativeMatrix.size1() != 1 || rExtractedDerivativeMatrix.size2() != rDerivativeMatrix.size2() )
            rExtractedDerivativeMatrix.resize(1, rDerivativeMatrix.size2());

        rExtractedDerivativeMatrix.clear();
                     
        unsigned int coord_dir = this->GetCoordinateDirection();
        const SizeType number_of_nodes = rDirectElement.GetGeometry().PointsNumber();
        for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id)
        { 
            IndexType index = node_id * 3;           
            if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), rDirectElement.GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() )
                for(IndexType i = 0; i < rDerivativeMatrix.size2(); ++i)                    
                    rExtractedDerivativeMatrix(0, i) += rDerivativeMatrix( index + coord_dir, i );
        }  
    }    
    
    std::vector<unsigned int>  DirectSensitivityNodalCoordinatesVariable::GetTracedElementId() 
    {              
        return mNeighboringElementIdVector;
    }

    void DirectSensitivityNodalCoordinatesVariable::GetNeighboringElementIdVector(std::vector<unsigned int>& rNeighboringElementIdVector)
    {
        KRATOS_TRY;
        
        if (rNeighboringElementIdVector.size() != 0)
            rNeighboringElementIdVector.resize(0);
        
        for (auto& node_i : mpTracedModelPart->Nodes())
        {
            GlobalPointersVector<Element >& r_neighbouring_elements = node_i.GetValue(NEIGHBOUR_ELEMENTS);            
            for (IndexType j = 0; j < r_neighbouring_elements.size(); ++j)
                if( !( std::find( rNeighboringElementIdVector.begin(), rNeighboringElementIdVector.end(), r_neighbouring_elements[j].Id() ) != rNeighboringElementIdVector.end() ) )
                    rNeighboringElementIdVector.push_back( r_neighbouring_elements[j].Id() );
        }
        
        KRATOS_CATCH("");
    }

    void DirectSensitivityNodalCoordinatesVariable::GetNeighboringConditionIdVector(std::vector<unsigned int>& rNeighboringConditionIdVector)
    {
        KRATOS_TRY;
        
        if (rNeighboringConditionIdVector.size() != 0)
            rNeighboringConditionIdVector.resize(0);
                
        for (auto cond_it = mrModelPart.Conditions().ptr_begin(); cond_it != mrModelPart.Conditions().ptr_end(); ++cond_it)
        {
            const SizeType number_of_nodes = (*cond_it)->GetGeometry().PointsNumber();
            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id)
                if( std::find( mpTracedModelPart->Nodes().begin(), mpTracedModelPart->Nodes().end(), (*cond_it)->GetGeometry()[node_id] ) != mpTracedModelPart->Nodes().end() )
                    if( !( std::find( rNeighboringConditionIdVector.begin(), rNeighboringConditionIdVector.end(), (*cond_it)->Id() ) != rNeighboringConditionIdVector.end() ) )
                        rNeighboringConditionIdVector.push_back( (*cond_it)->Id() );
        }
        
        KRATOS_CATCH("");
    }        

    unsigned int DirectSensitivityNodalCoordinatesVariable::GetCoordinateDirection()
    {        
        if( mCoordinateDirection == "X" )
            return 0;
        else if( mCoordinateDirection == "Y" )
            return 1;
        else if( mCoordinateDirection == "Z" )
            return 2;
        else
            KRATOS_ERROR << "Invalid direction: " << mCoordinateDirection <<
                ". Available is only 'X','Y' or 'Z'-Direction " << std::endl;        
    } 
}    