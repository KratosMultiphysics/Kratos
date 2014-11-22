//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Nelson $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_SPHERIC_PARTICLE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"
#include "../custom_conditions/RigidFace.h"
#include "../custom_conditions/dem_wall.h"
#include "kratos_DEMApplication_export_dll.h"
#include "../custom_utilities/properties_proxies.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class DEMWall;
class KRATOS_DEMAPPLICATION_EXPORT_DLL SphericParticle : public DiscreteElement
{
public:

///@name Type Definitions
///@{

/// Pointer definition of SphericParticle
KRATOS_CLASS_POINTER_DEFINITION(SphericParticle);

typedef WeakPointerVector<Condition> ConditionWeakVectorType;
typedef WeakPointerVector<Condition >::iterator ConditionWeakIteratorType;

typedef WeakPointerVector<Element> ParticleWeakVectorType;
typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

///@}
///@name Life Cycle
///@{

/// Default constructor.

SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

/// Destructor.
virtual ~SphericParticle();

///@}
///@name Operations
///@{
void Initialize();
void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity);
void FirstCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo, double dt);
void CollectCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo);
void FinalCalculateRightHandSide(ProcessInfo& rCurrentProcessInfo, double dt, const array_1d<double,3>& gravity);
void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);
void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);
void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo);
void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);
void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);
void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo);
void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo);

void CalculateMaxBallToBallIndentation(double& rCurrentMaxIndentation);
void CalculateMaxBallToFaceIndentation(double& rCurrentMaxIndentation);

int   GetClusterId();
void  SetClusterId(const int Id);

double GetRadius();
void   SetRadius(double radius);
double GetSqrtOfRealMass();
void   SetSqrtOfRealMass(double sqrt_of_real_mass);
double GetYoung();
void   SetYoungFromProperties(double* young);
double GetRollingFriction();
void   SetRollingFrictionFromProperties(double* rolling_friction);
double GetPoisson();
void   SetPoissonFromProperties(double* poisson);
double GetTgOfFrictionAngle();
void   SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle);
double GetLnOfRestitCoeff();
void   SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff);
double GetDensity();
void   SetDensityFromProperties(double* density);
int    GetParticleMaterial();
void   SetParticleMaterialFromProperties(int* particle_material);

PropertiesProxy* GetFastProperties();
void   SetFastProperties(PropertiesProxy* pProps);
void   SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies);

double SlowGetYoung();
double SlowGetRollingFriction();
double SlowGetPoisson();
double SlowGetTgOfFrictionAngle();
double SlowGetLnOfRestitCoeff();
double SlowGetDensity();
double SlowGetCohesion();
int SlowGetParticleMaterial();


///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.
virtual std::string Info() const
{
std::stringstream buffer;
buffer << "SphericParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SphericParticle";}

/// Print object's data.
virtual void PrintData(std::ostream& rOStream) const {}

virtual void ComputeNewNeighboursHistoricalData( std::vector<int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                           std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces);

virtual void ComputeNewRigidFaceNeighboursHistoricalData();

std::vector<SphericParticle*> mNeighbourElements;
std::vector<DEMWall*>         mNeighbourRigidFaces;
std::vector<double>           mNeighbourRigidFacesPram;
std::vector<double>           mNeighbourRigidFacesTotalContactForce;
std::vector<double>           mNeighbourRigidFacesElasticContactForce;

///@}
///@name Friends
///@{
virtual void ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& rCurrentProcessInfo, const array_1d<double,3>& gravity);
virtual void MemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);

protected:

SphericParticle();

virtual void ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rElasticForce,
                                          array_1d<double, 3>& InitialRotaMoment,
                                          ProcessInfo& rCurrentProcessInfo,
                                          double mTimeStep);

virtual void ComputeRigidFaceToMeVelocity(DEMWall* rObj_2, std::size_t ino, double LocalCoordSystem[3][3],double & DistPToB, array_1d<double, 3 > &other_to_me_vel, int & ContactType);

virtual void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);
void CalculateKineticEnergy(double& rKineticEnergy);
void CalculateElasticEnergyOfContacts(double& rElasticEnergy);
void CalculateMomentum(array_1d<double, 3>& rMomentum);
void CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum);
virtual void ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                     array_1d<double, 3>& InitialRotaMoment,
                                     ProcessInfo& rCurrentProcessInfo,
                                     double dt,
                                     const bool multi_stage_RHS);

virtual void CalculateEquivalentConstitutiveParameters(array_1d<double, 3>& other_to_me_vect,
                                                 const double& other_radius,
                                                 const double& radius_sum,
                                                 double& kn,
                                                 double& kt,
                                                 double& equiv_visco_damp_coeff_normal,
                                                 double& equiv_visco_damp_coeff_tangential,
                                                 double& equiv_tg_of_fri_ang,
                                                 SphericParticle* neighbour_iterator);

virtual void EvaluateDeltaDisplacement(double DeltDisp[3],
                    double RelVel[3],
                    double LocalCoordSystem[3][3],
                    double OldLocalCoordSystem[3][3],
                    array_1d<double, 3> &other_to_me_vect,
                    const array_1d<double, 3> &vel,
                    const array_1d<double, 3> &delta_displ,
                    SphericParticle* neighbour_iterator,
                    double& distance);

virtual void NormalForceCalculation(double& normal_force, double kn, double indentation);

virtual void CohesionCalculation(double& cohesion_force, double cohesion);

virtual void TangentialForceCalculation(const double normal_force, double LocalElasticContactForce[3], double LocalDeltDisp[3], const double& kt, const double& equiv_tg_of_fri_ang, bool& sliding);

virtual void DisplacementDueToRotation(double DeltDesp[3],
                    double OldLocalCoordSystem[3][3],
                    const double &other_radius,
                    const double &dt,
                    const array_1d<double, 3> &angl_vel,
                    SphericParticle* neighbour_iterator);

virtual void ComputeMoments(double normalLocalElasticContactForce,
                      array_1d<double, 3>& GlobalElasticContactForces,
                      array_1d<double, 3>& rInitialRotaMoment,
                      double LocalCoordSystem_2[3],
                      SphericParticle* neighbour_iterator);

virtual void CustomInitialize();

virtual double GetInitialDelta(int index);

virtual void AddUpForcesAndProject(double OldCoordSystem[3][3],
                    double LocalCoordSystem[3][3],
                    double normal_force,
                    double cohesion_force,
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    array_1d<double, 3> &rElasticForce,
                    const double &i_neighbour_count);

virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    array_1d<double, 3> &rElasticForce,
                    const double &iRigidFaceNeighbour);

virtual void CreateDiscontinuumConstitutiveLaws(const ProcessInfo& rCurrentProcessInfo);

virtual void CalculateViscoDamping(double LocalRelVel[3],
                                      double ViscoDampingLocalContactForce[3],
                                      double indentation,
                                      double equiv_visco_damp_coeff_normal,
                                      double equiv_visco_damp_coeff_tangential,
                                      bool sliding);

virtual void AdditionalMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo);
virtual void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo);

      DEMDiscontinuumConstitutiveLaw::Pointer mDiscontinuumConstitutiveLaw;


//const int mParticleId; // (NOT YET ACTIVE!!) Identifies the particle biunivocally if it has been properly created (i.e., a non-repeated NewId is passed to the constructor)
int mDimension;
int mDampType;
int mElasticityType;

array_1d<double, 3> mContactForce;
array_1d<double, 3> mContactMoment;

double mRadius;
double mSqrtOfRealMass;

PropertiesProxy* mFastProperties;

std::vector<int> mOldNeighbourIds;
std::vector< array_1d<double, 3> > mOldNeighbourElasticContactForces;
std::vector< array_1d<double, 3> > mOldNeighbourTotalContactForces;

std::vector<int> mFemOldNeighbourIds;
std::vector< array_1d<double, 3> >  mFemOldNeighbourContactForces;

int mClusterId;

///@name Protected static Member Variables
///@{


///@}
///@name Protected member Variables
///@{


///@}
///@name Protected Operators
///@{


///@}
///@name Protected Operations
///@{


///@}
///@name Protected  Access
///@{


///@}
///@name Protected Inquiry
///@{


///@}
///@name Protected LifeCycle
///@{


///@}

private:


///@name Static Member Variables
///@{


///@}
///@name Member Variables
///@{


///@}
///@name Private Operators
///@{


///@}
///@name Private Operations
///@{


///@}
///@name Private  Access
///@{


///@}
///@name Private Inquiry
///@{


///@}
///@name Un accessible methods
///@{


///@}
///@name Serialization
///@{

friend class Serializer;

virtual void save(Serializer& rSerializer) const
{
KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
rSerializer.save("mRadius",mRadius);
rSerializer.save("mSqrtOfRealMass",mSqrtOfRealMass);
//rSerializer.save("mFastProperties",mFastProperties);

/*rSerializer.save("mRollingFriction",mRollingFriction);
rSerializer.save("mYoung",mYoung);
rSerializer.save("mPoisson",mPoisson);
rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);  */
}

virtual void load(Serializer& rSerializer)
{
KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
rSerializer.load("mRadius",mRadius);
rSerializer.load("mSqrtOfRealMass",mSqrtOfRealMass);
//rSerializer.load("mFastProperties",mFastProperties);

/*rSerializer.load("mRollingFriction",mRollingFriction);
rSerializer.load("mYoung",mYoung);
rSerializer.load("mPoisson",mPoisson);
rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);  */
}

/*
/// Assignment operator.
SphericParticle& operator=(SphericParticle const& rOther)
{
return *this;
}

/// Copy constructor.
SphericParticle(SphericParticle const& rOther)
{
*this = rOther;
}
*/

///@}

}; // Class SphericParticle

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
            SphericParticle& rThis){ return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
            const SphericParticle& rThis)
{
rThis.PrintInfo(rOStream);
rOStream << std::endl;
rThis.PrintData(rOStream);

return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_PARTICLE_H_INCLUDED  defined


