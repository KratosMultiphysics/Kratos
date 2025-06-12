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
#include "custom_elements/mpm_updated_lagrangian_PQ.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"


#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/mpm_math_utilities.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"


namespace Kratos
{

MPMUpdatedLagrangianPQ::MPMUpdatedLagrangianPQ( )
    : MPMUpdatedLagrangian( )
{ }//DO NOT CALL IT: only needed for Register and Serialization!!!

MPMUpdatedLagrangianPQ::MPMUpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMUpdatedLagrangian( NewId, pGeometry )
{ }//DO NOT ADD DOFS HERE!!!

MPMUpdatedLagrangianPQ::MPMUpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : MPMUpdatedLagrangian( NewId, pGeometry, pProperties )
{ mFinalizedStep = true; }

MPMUpdatedLagrangianPQ::MPMUpdatedLagrangianPQ(MPMUpdatedLagrangianPQ const& rOther)
    :MPMUpdatedLagrangian(rOther)
{ }

MPMUpdatedLagrangianPQ& MPMUpdatedLagrangianPQ::operator=(MPMUpdatedLagrangianPQ const& rOther)
{
    MPMUpdatedLagrangian::operator=(rOther);
    return *this;
}

Element::Pointer MPMUpdatedLagrangianPQ::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{ return Element::Pointer( new MPMUpdatedLagrangianPQ( NewId, GetGeometry().Create( ThisNodes ), pProperties ) ); }

Element::Pointer MPMUpdatedLagrangianPQ::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{ return Kratos::make_intrusive< MPMUpdatedLagrangianPQ >(NewId, pGeom, pProperties); }

Element::Pointer MPMUpdatedLagrangianPQ::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    MPMUpdatedLagrangianPQ NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );
    return Element::Pointer( new MPMUpdatedLagrangianPQ(NewElement) );
}

MPMUpdatedLagrangianPQ::~MPMUpdatedLagrangianPQ()
{ }

void MPMUpdatedLagrangianPQ::CalculateAndAddExternalForces(
    VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    for (IndexType int_p = 0; int_p < GetGeometry().IntegrationPointsNumber(); ++int_p)
    {
        double weight = (GetGeometry().IntegrationPointsNumber() > 1) ? GetGeometry().IntegrationPoints()[int_p].Weight() : 1.0;
        for (IndexType i = 0; i < number_of_nodes; ++i)
        {
            if (GetGeometry().ShapeFunctionValue(int_p, i) >= 0.0) // skip inactive nodes
            {
                for (unsigned int j = 0; j < dimension; ++j)
                {
                    rRightHandSideVector[dimension * i + j] += GetGeometry().ShapeFunctionValue(int_p, i) *
                        rVolumeForce[j] * weight;
                }
            }
        }
    }

    KRATOS_CATCH( "" )
}

void MPMUpdatedLagrangianPQ::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    mFinalizedStep = false;

    // Calculating shape functions
    array_1d<double, 3> nodal_momentum = ZeroVector(3);
    array_1d<double, 3>  nodal_inertia = ZeroVector(3);
    // Here MP contribution in terms of momentum, inertia and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (IndexType int_p = 0; int_p < GetGeometry().IntegrationPointsNumber(); ++int_p)
        {
            double weight = (GetGeometry().IntegrationPointsNumber() > 1) ? GetGeometry().IntegrationPoints()[int_p].Weight() : 1.0;
            if (r_geometry.ShapeFunctionValue(int_p, i) >= 0.0) // skip inactive nodes
            {
                for (unsigned int j = 0; j < dimension; j++)
                {
                    nodal_momentum[j] = r_geometry.ShapeFunctionValue(int_p, i) * mMP.velocity[j] *
                        mMP.mass * weight;
                    nodal_inertia[j] = r_geometry.ShapeFunctionValue(int_p, i) * mMP.acceleration[j] *
                        mMP.mass * weight;
                }

                // Add in the predictor velocity increment for central difference explicit
                // This is the 'previous grid acceleration', which is actually
                // be the initial material point acceleration mapped to the grid.
                if (rCurrentProcessInfo.Has(IS_EXPLICIT_CENTRAL_DIFFERENCE)) {
                    if (rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE)) {
                        const double& delta_time = rCurrentProcessInfo[DELTA_TIME];
                        for (unsigned int j = 0; j < dimension; j++) {
                            nodal_momentum[j] += 0.5 * delta_time * (r_geometry.ShapeFunctionValue(int_p, i) *
                                mMP.acceleration[j]) * mMP.mass * weight;
                        }
                    }
                }

                r_geometry[i].SetLock();
                r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
                r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0) += nodal_inertia;
                r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) += r_geometry.ShapeFunctionValue(int_p, i)
                    * mMP.mass * weight;
                r_geometry[i].UnSetLock();
            }
        }
    }
}


void MPMUpdatedLagrangianPQ::InitializeMaterial(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        Vector N_dummy; // this is because the shape functions are not explicitly defined at the 'master' material point anymore
        GetGeometry().SetValue(MP_VOLUME, mMP.volume);
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


void MPMUpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1) rValues.resize(1);
    if (rVariable == MP_SUB_POINTS) rValues[0] = GetGeometry().IntegrationPointsNumber();
    else MPMUpdatedLagrangian::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}


void MPMUpdatedLagrangianPQ::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
}

void MPMUpdatedLagrangianPQ::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangian )
}
} // Namespace Kratos

