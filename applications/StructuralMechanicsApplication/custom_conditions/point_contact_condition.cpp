// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_contact_condition.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

PointContactCondition::PointContactCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

PointContactCondition::PointContactCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer PointContactCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PointContactCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointContactCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PointContactCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PointContactCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PointContactCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

PointContactCondition::~PointContactCondition()
{
}

//************************************************************************************
//************************************************************************************

void PointContactCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int NumberOfNodes = GetGeometry().size();
    const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes * Dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
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

    //obtain distance, gradient of distance and nodal normal and displacement at which the distance was measured
    const double d_reference = GetGeometry()[0].GetValue(DISTANCE);
    const array_1d<double,3>& grad_d = GetGeometry()[0].GetValue(DISTANCE_GRADIENT);
    const array_1d<double,3>& n = GetGeometry()[0].FastGetSolutionStepValue(NORMAL); //this normal is proportional to the area
    const array_1d<double,3>& reference_disp = GetGeometry()[0].GetValue(DISPLACEMENT); //this is the displacement measured at the moment in which the d_reference was computed
//    KRATOS_WATCH(grad_d)
//    KRATOS_WATCH(n)
//    KRATOS_WATCH(reference_disp)
    const array_1d<double,3> ddisp = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - reference_disp;
    const double d_current = d_reference + inner_prod(grad_d,ddisp);
//    KRATOS_WATCH(ddisp)
//    KRATOS_WATCH(d_current)
//        if(d_reference > 0) //
        const double E = GetProperties()[YOUNG_MODULUS]/10000.0;
        const double characteristic_lenght = 1.0e-2; //should get this from a variable
        const double spring_stiffness = E/characteristic_lenght;
        const double contact_pressure = d_current*spring_stiffness;

    if(d_current > 0)
    {

// KRATOS_WATCH(Dimension)
        // Vector with a loading applied to the condition
        array_1d<double, 3 > PointLoad = contact_pressure*n;
// KRATOS_WATCH(PointLoad)

        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            const unsigned int base = ii*Dimension;

            for(unsigned int k = 0; k < Dimension; ++k)
                rRightHandSideVector[base + k] = -PointLoad[k];

            noalias(GetGeometry()[0].FastGetSolutionStepValue(FORCE)) = -PointLoad; //TODO: remove
            GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = d_current;
            GetGeometry()[0].FastGetSolutionStepValue(NODAL_PAUX) = d_reference;


//                 if ( CalculateStiffnessMatrixFlag == true )
//                 {
//                     for(unsigned int k = 0; k < Dimension; ++k)
//                     {
// //                 KRATOS_WATCH(rRightHandSideVector)
// //                 KRATOS_WATCH(rLeftHandSideMatrix)
//                         for(unsigned int l = 0; l < Dimension; ++l)
//                         {
// //                             rLeftHandSideMatrix(base+k,base+l) = spring_stiffness*(n[k]*n[l])/norm_2(n);
//                              rLeftHandSideMatrix(base+k,base+l) = spring_stiffness*(n[k]*grad_d[l]); //this would be the proper tangent
//                         }
//                     }
//                 }
        }
    }
    else
    {
            GetGeometry()[0].FastGetSolutionStepValue(FORCE).clear(); //TODO: remove
            GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = 0.0;
            GetGeometry()[0].FastGetSolutionStepValue(NODAL_PAUX) = d_reference;
    }
//             KRATOS_WATCH(n)
//             KRATOS_WATCH(grad_d)
//             KRATOS_WATCH(rLeftHandSideMatrix)

    if(d_current > 0 )//d_reference > 0)
    {
        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            const unsigned int base = ii*Dimension;

            if ( CalculateStiffnessMatrixFlag == true )
            {
                for(unsigned int k = 0; k < Dimension; ++k)
                {
//                 KRATOS_WATCH(rRightHandSideVector)
//                 KRATOS_WATCH(rLeftHandSideMatrix)
                    for(unsigned int l = 0; l < Dimension; ++l)
                    {
//                             rLeftHandSideMatrix(base+k,base+l) = spring_stiffness*(n[k]*n[l])/norm_2(n);
                            rLeftHandSideMatrix(base+k,base+l) = spring_stiffness*(n[k]*grad_d[l]); //this would be the proper tangent
                    }
                }
            }
        }
    }
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double PointContactCondition::GetPointLoadIntegrationWeight() const
{
    return 1.0;
}

} // Namespace Kratos


