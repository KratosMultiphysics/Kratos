//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
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
#include "../custom_conditions/RigidFace.h"
#include "../custom_conditions/dem_wall.h"
#include "kratos_DEMApplication_export_dll.h"
#include "../custom_utilities/properties_proxies.h"


namespace Kratos
{

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
virtual void Initialize();
virtual void FullInitialize(const ProcessInfo& r_process_info);
virtual void CreateDiscontinuumConstitutiveLaws(const ProcessInfo& r_process_info);
virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity, int search_control);
virtual void FirstCalculateRightHandSide(ProcessInfo& r_process_info, double dt, int search_control);
virtual void CollectCalculateRightHandSide(ProcessInfo& r_process_info);
virtual void FinalCalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity);
virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info);
virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info);
virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info);
virtual void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo );
virtual void FinalizeSolutionStep(ProcessInfo& r_process_info);
virtual void PrepareForPrinting(ProcessInfo& r_process_info);
virtual void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info);

virtual void CalculateMaxBallToBallIndentation(double& rCurrentMaxIndentation);
virtual void CalculateMaxBallToFaceIndentation(double& rCurrentMaxIndentation);

int   GetClusterId();
void  SetClusterId(const int Id);

double GetRadius();
void   SetRadius(double radius);
double GetVolume();
virtual double GetSearchRadius();
double GetMass();
void   SetMass(double real_mass);
double GetYoung();
void   SetYoungFromProperties(double* young);
double GetRollingFriction();
void   SetRollingFrictionFromProperties(double* rolling_friction);
double GetPoisson();
void   SetPoissonFromProperties(double* poisson);
double GetTgOfFrictionAngle();
void   SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle);
double GetCoefficientOfRestitution();
void   SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution);
double GetLnOfRestitCoeff();
void   SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff);
double GetDensity();
void   SetDensityFromProperties(double* density);
int    GetParticleMaterial();
void   SetParticleMaterialFromProperties(int* particle_material);
double GetParticleCohesion();
void   SetParticleCohesionFromProperties(double* particle_cohesion);

PropertiesProxy* GetFastProperties();
void   SetFastProperties(PropertiesProxy* pProps);
void   SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies);

double SlowGetYoung();
double SlowGetRollingFriction();
double SlowGetPoisson();
double SlowGetTgOfFrictionAngle();
double SlowGetCoefficientOfRestitution();
double SlowGetLnOfRestitCoeff();
double SlowGetDensity();
double SlowGetParticleCohesion();
int    SlowGetParticleMaterial();

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

virtual void ComputeNewNeighboursHistoricalData( std::vector<unsigned int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
                                           std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces);

virtual void ComputeNewRigidFaceNeighboursHistoricalData();

virtual void CalculateKinematicEnergy(double& rKinematicEnergy);
virtual void CalculateGravitationalEnergy(const array_1d<double,3>& gravity, double& r_gravitational_energy);

std::vector<SphericParticle*> mNeighbourElements;
std::vector<DEMWall*>         mNeighbourRigidFaces;
std::vector<double>           mNeighbourRigidFacesPram;
std::vector<array_1d<double, 3> >           mNeighbourRigidFacesTotalContactForce;
std::vector<array_1d<double, 3> >           mNeighbourRigidFacesElasticContactForce;

///@}
///@name Friends
///@{
virtual void ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, ProcessInfo& r_process_info, const array_1d<double,3>& gravity);
virtual void MemberDeclarationFirstStep(const ProcessInfo& r_process_info);

array_1d<double, 3> mContactForce; //SLS
array_1d<double, 3> mContactMoment; //SLS

Matrix* mStressTensor;
Matrix* mSymmStressTensor;

protected:

SphericParticle();

virtual void ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rElasticForce,
                                          array_1d<double, 3>& rContactForce,
                                          array_1d<double, 3>& InitialRotaMoment,
                                          array_1d<double, 3>& rigid_element_force,
                                          ProcessInfo& r_process_info,
                                          double mTimeStep,
                                          int search_control);

virtual void ComputeRigidFaceToMeVelocity(DEMWall* rObj_2, 
                                        std::size_t ino, 
                                        double LocalCoordSystem[3][3],
                                        double & DistPToB, 
                                        double Weight[4], 
                                        array_1d<double,3>& wall_delta_disp_at_contact_point,
                                        array_1d<double, 3 > &wall_velocity_at_contact_point, 
                                        int & ContactType);

virtual void UpdateDistanceToWall(DEMWall* const wall, 
                                    const int neighbour_index, 
                                    double LocalCoordSystem[3][3], 
                                    double& DistPToB, 
                                    double Weight[4], 
                                    int& ContactType);

virtual void UpdateRF_Pram(DEMWall* rObj_2, 
                            const std::size_t ino,
                            const double LocalCoordSystem[3][3], 
                            const double DistPToB,
                            const double Weight[4], 
                            const int ContactType);

virtual void InitializeSolutionStep(ProcessInfo& r_process_info);

virtual void CalculateMomentum(array_1d<double, 3>& rMomentum);
virtual void CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum);
virtual void ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                     array_1d<double, 3>& rContactForce,
                                     array_1d<double, 3>& InitialRotaMoment,
                                     ProcessInfo& r_process_info,
                                     double dt,
                                     const bool multi_stage_RHS);

virtual void EvaluateDeltaDisplacement(double DeltDisp[3],
                    double RelVel[3],
                    double LocalCoordSystem[3][3],
                    double OldLocalCoordSystem[3][3],
                    array_1d<double, 3> &other_to_me_vect,
                    const array_1d<double, 3> &vel,
                    const array_1d<double, 3> &delta_displ,
                    SphericParticle* neighbour_iterator,
                    double& distance);

virtual void DisplacementDueToRotation(const double indentation,
                    double DeltDesp[3],
                    double RelVel[3],
                    double OldLocalCoordSystem[3][3],
                    const double &other_radius,
                    const double &dt,
                    const array_1d<double, 3> &angl_vel,
                    SphericParticle* neighbour_iterator);

virtual void DisplacementDueToRotationMatrix(double DeltDisp[3],
                                                double RelVel[3],
                                                double OldLocalCoordSystem[3][3],
                                                const double& other_radius,
                                                const double& dt,
                                                const array_1d<double, 3>& ang_vel,
                                                SphericParticle* p_neighbour);

virtual void ComputeMoments(double normalLocalElasticContactForce,
                      array_1d<double, 3>& GlobalElasticContactForces,
                      array_1d<double, 3>& rInitialRotaMoment,
                      double LocalCoordSystem_2[3],
                      SphericParticle* neighbour_iterator,
                      double indentation,
                      bool wall=false);

virtual void CustomInitialize();

virtual double GetInitialDeltaWithFEM(int index);

virtual void AddUpForcesAndProject(double OldCoordSystem[3][3],
                    double LocalCoordSystem[3][3],
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    const double cohesive_force,
                    array_1d<double, 3>& rElasticForce,
                    array_1d<double, 3>& rContactForce,
                    const unsigned int i_neighbour_count);

virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    const double cohesive_force,
                    array_1d<double, 3>& rElasticForce,
                    array_1d<double, 3>& rContactForce,
                    const unsigned int iRigidFaceNeighbour);

virtual void ComputeWear(double LocalCoordSystem[3][3], array_1d<double, 3>& vel, double tangential_vel[3],
                         double mTimeStep, double density, bool sliding, double inverse_of_volume,
                         double LocalElasticContactForce, DEMWall* cast_neighbour);

virtual void AdditionalMemberDeclarationFirstStep(const ProcessInfo& r_process_info);
virtual void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);

virtual void AddNeighbourContributionToStressTensor(const double GlobalElasticContactForce[3],
                                                    const double other_to_me_vect[3],
                                                    const double distance,
                                                    const double radius_sum);

virtual void AddContributionToRepresentativeVolume(const double distance,
                                                    const double radius_sum,
                                                    const double contact_area);

virtual void AddWallContributionToStressTensor(const double GlobalElasticContactForce[3],
                                                const double other_to_me_vect[3],
                                                const double distance,
                                                const double contact_area);

DEMDiscontinuumConstitutiveLaw::Pointer mDiscontinuumConstitutiveLaw;


//const int mParticleId; // (NOT YET ACTIVE!!) Identifies the particle biunivocally if it has been properly created (i.e., a non-repeated NewId is passed to the constructor)

//const double* mSearchControl;

//array_1d<double, 3> mContactForce; //SLS

double mRadius;
double mRealMass;

PropertiesProxy* mFastProperties;

std::vector<unsigned int> mOldNeighbourIds;
std::vector< array_1d<double, 3> > mNeighbourElasticContactForces;
std::vector< array_1d<double, 3> > mNeighbourTotalContactForces;

std::vector<int> mFemOldNeighbourIds;

int mClusterId;

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const
{
KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
rSerializer.save("mRadius",mRadius);
rSerializer.save("mRealMass",mRealMass);
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
rSerializer.load("mRealMass",mRealMass);
//rSerializer.load("mFastProperties",mFastProperties);

/*rSerializer.load("mRollingFriction",mRollingFriction);
rSerializer.load("mYoung",mYoung);
rSerializer.load("mPoisson",mPoisson);
rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);  */
}


}; // Class SphericParticle


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


