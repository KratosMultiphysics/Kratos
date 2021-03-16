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
#include "custom_elements/updated_lagrangian_CPDI.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"
#include "includes/checks.h"


#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"


namespace Kratos
{

UpdatedLagrangianCPDI::UpdatedLagrangianCPDI( )
    : UpdatedLagrangian( )
{ }//DO NOT CALL IT: only needed for Register and Serialization!!!

UpdatedLagrangianCPDI::UpdatedLagrangianCPDI( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
{ }//DO NOT ADD DOFS HERE!!!

UpdatedLagrangianCPDI::UpdatedLagrangianCPDI( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
{ mFinalizedStep = true; }

UpdatedLagrangianCPDI::UpdatedLagrangianCPDI(UpdatedLagrangianCPDI const& rOther)
    :UpdatedLagrangian(rOther)
{ }

UpdatedLagrangianCPDI& UpdatedLagrangianCPDI::operator=(UpdatedLagrangianCPDI const& rOther)
{
    UpdatedLagrangian::operator=(rOther);
    return *this;
}

Element::Pointer UpdatedLagrangianCPDI::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{ return Element::Pointer( new UpdatedLagrangianCPDI( NewId, GetGeometry().Create( ThisNodes ), pProperties ) ); }

Element::Pointer UpdatedLagrangianCPDI::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{ return Kratos::make_intrusive< UpdatedLagrangianCPDI >(NewId, pGeom, pProperties); }

Element::Pointer UpdatedLagrangianCPDI::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    UpdatedLagrangianCPDI NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );
    return Element::Pointer( new UpdatedLagrangianCPDI(NewElement) );
}

UpdatedLagrangianCPDI::~UpdatedLagrangianCPDI()
{ }

void UpdatedLagrangianCPDI::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    KRATOS_WATCH(number_of_nodes)

    mFinalizedStep = false;

    // Calculating shape functions
    array_1d<double, 3> nodal_momentum = ZeroVector(3);
    array_1d<double, 3>  nodal_inertia = ZeroVector(3);
    // Here MP contribution in terms of momentum, inertia and mass are added
    
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (IndexType int_p = 0; int_p < GetGeometry().IntegrationPointsNumber(); ++int_p)
        {
            double weight = 1.0;
            KRATOS_WATCH(r_geometry.ShapeFunctionValue(int_p, i))
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
                // be the initial particle acceleration mapped to the grid.
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
                KRATOS_WATCH(weight)
                
                r_geometry[i].UnSetLock();
            }
        }
    }
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        KRATOS_WATCH(r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0))
    }

}


void UpdatedLagrangianCPDI::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1) rValues.resize(1);
    if (rVariable == MP_SUB_POINTS) rValues[0] = GetGeometry().IntegrationPointsNumber();
    else UpdatedLagrangian::CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
}


void UpdatedLagrangianCPDI::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, UpdatedLagrangian )
}

void UpdatedLagrangianCPDI::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, UpdatedLagrangian )
}
} // Namespace Kratos

