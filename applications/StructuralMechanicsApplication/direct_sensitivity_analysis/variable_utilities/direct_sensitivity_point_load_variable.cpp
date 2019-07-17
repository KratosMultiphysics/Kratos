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
#include "direct_sensitivity_point_load_variable.h"


namespace Kratos
{ 
    // Constructor
    DirectSensitivityPointLoadVariable::DirectSensitivityPointLoadVariable(ModelPart& rModelPart, Parameters VariableSettings)
    : DirectSensitivityVariable(rModelPart, VariableSettings)
    {
        // Get traced conditions
        if (mTracedConditionIdVector.size() != 0)
            mTracedConditionIdVector.resize(0);

        for (auto& cond_i : mpTracedModelPart->Conditions())
            mTracedConditionIdVector.push_back(cond_i.Id());

        // Set mTracedElementIdVector to zero
        mTracedElementIdVector.resize(0);

        // Direction in which the intensity of the traced point loads should be examined
        mCoordinateDirection = VariableSettings["coordinate_direction"].GetString();             
        
    }

    // Destructor
    DirectSensitivityPointLoadVariable::~DirectSensitivityPointLoadVariable(){}


    void DirectSensitivityPointLoadVariable::Initialize()
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    void DirectSensitivityPointLoadVariable::CalculatePseudoLoadVector(Condition& rDirectCondition, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        
        KRATOS_TRY;

        // Define working variables
       unsigned int coord_dir = GetCoordinateDirection();

        // Define sizes
        const SizeType number_of_nodes = rDirectCondition.GetGeometry().PointsNumber();
        const SizeType dimension = rDirectCondition.GetGeometry().WorkingSpaceDimension();
        const SizeType vec_size = number_of_nodes * dimension;        

        // Sizing and clearing of the pseudo load vector
        if (rPseudoLoadVector.size() != vec_size)
            rPseudoLoadVector.resize(vec_size, false);

        rPseudoLoadVector.clear();        

        // Calculate pseudo load
        if (std::find(mTracedConditionIdVector.begin(), mTracedConditionIdVector.end(), rDirectCondition.Id() ) != mTracedConditionIdVector.end())    
            for(IndexType node_id = 0; node_id < number_of_nodes; ++node_id)
            { 
                IndexType index = node_id * dimension;                
                rPseudoLoadVector[index + coord_dir] = 1; 
            }
        else
            rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size()); 
        
        KRATOS_CATCH("");          
    }

    void DirectSensitivityPointLoadVariable::CalculatePseudoLoadVector(Element& rDirectElement, const Matrix& rLHS, Vector& rPseudoLoadVector, 
                                            const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;

        // Sizing of the pseudo load vector
        if (rPseudoLoadVector.size() != rLHS.size2())
            rPseudoLoadVector.resize(rLHS.size2(), false);        

        // Set pseudo load to zero    
        rPseudoLoadVector = ZeroVector(rPseudoLoadVector.size());
      

        KRATOS_CATCH("");
    }
    
    std::vector<unsigned int>  DirectSensitivityPointLoadVariable::GetTracedElementId() 
    {              
        return mTracedElementIdVector;
    } 

    unsigned int DirectSensitivityPointLoadVariable::GetCoordinateDirection()
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