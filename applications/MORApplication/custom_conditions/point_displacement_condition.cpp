// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: MOR_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_displacement_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "custom_elements/acoustic_element.h"
#include "mor_application_variables.h"


namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

PointDisplacementCondition::PointDisplacementCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseDisplacementCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

PointDisplacementCondition::PointDisplacementCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseDisplacementCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer PointDisplacementCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PointDisplacementCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointDisplacementCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PointDisplacementCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PointDisplacementCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PointDisplacementCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

PointDisplacementCondition::~PointDisplacementCondition()
{
}

//************************************************************************************
//************************************************************************************

void PointDisplacementCondition:: CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag,
        const bool CalculateMassMatrixFlag
        )
{
    KRATOS_TRY

    const unsigned int NumberOfNodes = GetGeometry().size();
    auto& r_geometry   = GetGeometry();
    const auto& r_prop = GetProperties();
    const double density   = r_prop[DENSITY];
    
    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //Check required to validate the need
    if ( CalculateMassMatrixFlag ==true ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); 
    }
    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
        {
            rRightHandSideVector.resize( MatSize, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    // Vector with a loading applied to the condition
    double AcousticDisplacement = 0.0;
    if( this->Has( ACOUSTIC_DISPLACEMENT ) )
        {
            AcousticDisplacement = this->GetValue( ACOUSTIC_DISPLACEMENT );
        }

    for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
    {
        

        if( GetGeometry()[ii].SolutionStepsDataHas( ACOUSTIC_DISPLACEMENT ) )
        {
            AcousticDisplacement += GetGeometry()[ii].FastGetSolutionStepValue( ACOUSTIC_DISPLACEMENT );
         
        }
       
        rRightHandSideVector[ii] += GetPointLoadIntegrationWeight() * 1.21 * rCurrentProcessInfo[FREQUENCY] * rCurrentProcessInfo[FREQUENCY] * AcousticDisplacement;
    }

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

double PointDisplacementCondition::GetPointLoadIntegrationWeight() const
{
    return 1.0 ;
}


} // Namespace Kratos


