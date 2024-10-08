//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_elements/Particle_Contact_Element.h"
#include "utilities/math_utils.h"
#include "DEM_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
ParticleContactElement::ParticleContactElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
ParticleContactElement::ParticleContactElement( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
}

//create contact elements instances.

Element::Pointer ParticleContactElement::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new ParticleContactElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

ParticleContactElement::~ParticleContactElement()
{
}

std::string ParticleContactElement::Info() const
{
    std::stringstream buffer;
    buffer << "Particle Contact Element" << std::endl;
    return buffer.str();
}

void ParticleContactElement::Initialize(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    mFailureCriterionState = 0.0;
    mLocalContactForce[0] = 0.0;
    mLocalContactForce[1] = 0.0;
    mLocalContactForce[2] = 0.0;
    mElasticLocalRotationalMoment[0] = 0.0;
    mElasticLocalRotationalMoment[1] = 0.0;
    mElasticLocalRotationalMoment[2] = 0.0;
    mUnidimendionalDamage = 0.0;
    mContactFailure = 0.0;
    mContactSigma = 0.0;
    mContactTau = 0.0;
    mContactRadius = 0.0;

    array_1d<double, 3> vector_of_zeros(3,0.0);
    this->SetValue(LOCAL_CONTACT_FORCE, vector_of_zeros);
    this->SetValue(GLOBAL_CONTACT_FORCE, vector_of_zeros);
    this->SetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT, vector_of_zeros);
    this->SetValue(CONTACT_SIGMA, 0.0);
    this->SetValue(CONTACT_TAU, 0.0);
    this->SetValue(CONTACT_FAILURE, 0.0);
    this->SetValue(FAILURE_CRITERION_STATE, 0.0);
    this->SetValue(UNIDIMENSIONAL_DAMAGE, 0.0);
    this->SetValue(CONTACT_RADIUS, 0.0);

    KRATOS_CATCH( "" )
}

void ParticleContactElement::PrepareForPrinting() {
    KRATOS_TRY

    this->GetValue(LOCAL_CONTACT_FORCE)[0] = mLocalContactForce[0];
    this->GetValue(LOCAL_CONTACT_FORCE)[1] = mLocalContactForce[1];
    this->GetValue(LOCAL_CONTACT_FORCE)[2] = mLocalContactForce[2];
    this->GetValue(GLOBAL_CONTACT_FORCE)[0] = mGlobalContactForce[0];
    this->GetValue(GLOBAL_CONTACT_FORCE)[1] = mGlobalContactForce[1];
    this->GetValue(GLOBAL_CONTACT_FORCE)[2] = mGlobalContactForce[2];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[0] = mElasticLocalRotationalMoment[0];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[1] = mElasticLocalRotationalMoment[1];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[2] = mElasticLocalRotationalMoment[2];
    this->GetValue(CONTACT_SIGMA)          = mContactSigma;
    this->GetValue(CONTACT_TAU)            = mContactTau;
    this->GetValue(CONTACT_FAILURE)        = mContactFailure;
    this->GetValue(FAILURE_CRITERION_STATE)= mFailureCriterionState;
    this->GetValue(UNIDIMENSIONAL_DAMAGE)  = mUnidimendionalDamage;
    this->GetValue(CONTACT_RADIUS)         = mContactRadius;

    KRATOS_CATCH( "" )
}

void ParticleContactElement::CalculateOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable, std::vector<array_1d<double,3> >& rOutput, const ProcessInfo& r_process_info)
{
    //if(rVariable == LOCAL_CONTACT_FORCE) {  //3D VARIABLE WITH COMPONENTS
    rOutput.resize(1);
    const ParticleContactElement* const_this = dynamic_cast< const ParticleContactElement* >(this); //To ensure we don't set the value here
    rOutput[0][0] = const_this->GetValue(rVariable)[0];
    rOutput[0][1] = const_this->GetValue(rVariable)[1];
    rOutput[0][2] = const_this->GetValue(rVariable)[2];
    //}
}

void ParticleContactElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& r_process_info) {
    Output.resize(1);
    const ParticleContactElement* const_this = dynamic_cast< const ParticleContactElement* >(this); //To ensure we don't set the value here
    Output[0] = double(const_this->GetValue(rVariable));
}

void ParticleContactElement::InitializeSolutionStep(const ProcessInfo& r_process_info )
{
    mContactTau           = 0.0;
    mContactSigma         = 0.0;
    mLocalContactForce[0] = 0.0;
    mLocalContactForce[1] = 0.0;
    mLocalContactForce[2] = 0.0;
    mElasticLocalRotationalMoment[0] = 0.0;
    mElasticLocalRotationalMoment[1] = 0.0;
    mElasticLocalRotationalMoment[2] = 0.0;
    if (mFailureCriterionState<1.0) {
        mFailureCriterionState = 0.0;
    } // else we keep it at 1.0.
}

////************************************************************************************
////************************************************************************************
void ParticleContactElement::FinalizeSolutionStep(const ProcessInfo& r_process_info) {}

void ParticleContactElement::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {}

} // Namespace Kratos


