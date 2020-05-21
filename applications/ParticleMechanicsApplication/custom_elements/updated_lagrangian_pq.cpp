//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_pq.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"
#include "includes/checks.h"


#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianPQ::UpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;
}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ(UpdatedLagrangianPQ const& rOther)
    :UpdatedLagrangian(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianPQ& UpdatedLagrangianPQ::operator=(UpdatedLagrangianPQ const& rOther)
{
    UpdatedLagrangian::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianPQ::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianPQ( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer UpdatedLagrangianPQ::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< UpdatedLagrangianPQ >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianPQ::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    UpdatedLagrangianPQ NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new UpdatedLagrangianPQ(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianPQ::~UpdatedLagrangianPQ()
{
}

//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianPQ::CalculateAndAddExternalForces( // TODO Merge into base element
    VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for (IndexType int_p = 0; int_p < GetGeometry().IntegrationPointsNumber(); ++int_p)
    {
        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            if (r_N(int_p, i) * GetGeometry().IntegrationPoints()[int_p].Weight() > std::numeric_limits<double>::epsilon()) // skip inactive nodes
            {
                int index = dimension * i;

                for (unsigned int j = 0; j < dimension; j++)
                {
                    rRightHandSideVector[index + j] += r_N(int_p, i) * rVolumeForce[j] * GetGeometry().IntegrationPoints()[int_p].Weight();
                }
            }
        }
    }

    KRATOS_CATCH( "" )
}

//*******************************************************************************************
//*******************************************************************************************
void UpdatedLagrangianPQ::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo ) // TODO merge this into normal element
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    mFinalizedStep = false;

    const bool is_explicit_central_difference = (rCurrentProcessInfo.Has(IS_EXPLICIT_CENTRAL_DIFFERENCE))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE)
        : false;

    // Calculating shape functions
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    array_1d<double, 3> nodal_momentum = ZeroVector(3);
    array_1d<double, 3>  nodal_inertia = ZeroVector(3);

    // Here MP contribution in terms of momentum, inertia and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (IndexType int_p = 0; int_p < GetGeometry().IntegrationPointsNumber(); ++int_p)
        {
            if (r_N(int_p, i)* GetGeometry().IntegrationPoints()[int_p].Weight() > std::numeric_limits<double>::epsilon())
            {
                for (unsigned int j = 0; j < dimension; j++)
                {
                    nodal_momentum[j] = r_N(int_p, i) * mMP.velocity[j] * mMP.mass* GetGeometry().IntegrationPoints()[int_p].Weight();
                    nodal_inertia[j] = r_N(int_p, i) * mMP.acceleration[j] * mMP.mass * GetGeometry().IntegrationPoints()[int_p].Weight();
                }

                // Add in the predictor velocity increment for central difference explicit
                // This is the 'previous grid acceleration', which is actually
                // be the initial particle acceleration mapped to the grid.
                if (is_explicit_central_difference) {
                    const double& delta_time = rCurrentProcessInfo[DELTA_TIME];
                    for (unsigned int j = 0; j < dimension; j++) {
                        nodal_momentum[j] += 0.5 * delta_time * (r_N(int_p, i) * mMP.acceleration[j]) * mMP.mass * GetGeometry().IntegrationPoints()[int_p].Weight();
                    }
                }

                r_geometry[i].SetLock();
                r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
                r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0) += nodal_inertia;
                r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) += r_N(int_p, i) * mMP.mass * GetGeometry().IntegrationPoints()[int_p].Weight();
                r_geometry[i].UnSetLock();
            }
        }
    }
}


void UpdatedLagrangianPQ::InitializeMaterial() // TODO keep for clarity and error catching
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        Vector N_dummy; // send uninitialized dummy shape functions to throw error if actually needed
        //N_dummy = row(GetGeometry().ShapeFunctionsValues(), 0);
        mConstitutiveLawVector->InitializeMaterial( 
            GetProperties(), GetGeometry(), N_dummy);

        mMP.almansi_strain_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());
        mMP.cauchy_stress_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());

        // Resize the deformation gradient if we are axisymmetric
        if (mConstitutiveLawVector->GetStrainSize() == 4) mDeformationGradientF0 = IdentityMatrix(3);
    }
    else
        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}



///@}
///@name Access Get Values
///@{

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1) rValues.resize(1);
    if (rVariable == MP_SUB_POINTS)   rValues[0] = GetGeometry().IntegrationPointsNumber();
    else UpdatedLagrangian::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}


///@}

void UpdatedLagrangianPQ::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void UpdatedLagrangianPQ::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}
} // Namespace Kratos

