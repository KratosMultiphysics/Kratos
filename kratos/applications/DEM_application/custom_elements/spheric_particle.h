//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "discrete_element.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "../custom_conditions/RigidFace.h"
#include "../custom_conditions/dem_wall.h"
#include "includes/kratos_export_api.h"
//#include "../kratos_DEMApplication_export_dll.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/kratos_flags.h"


namespace Kratos
{

class DEMWall;

class KRATOS_API(DEM_APPLICATION) SphericParticle : public DiscreteElement
{
public:

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

using DiscreteElement::Initialize; //To avoid Clang Warning. We tell the compiler that we are aware of the existence of this function, but we overload it still.
virtual void Initialize(const ProcessInfo& r_process_info);
virtual void CreateDiscontinuumConstitutiveLaws(const ProcessInfo& r_process_info);
using DiscreteElement::CalculateRightHandSide; //To avoid Clang Warning. We tell the compiler that we are aware of the existence of this function, but we overload it still.
virtual void CalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity, int search_control);
virtual void FirstCalculateRightHandSide(ProcessInfo& r_process_info, double dt, int search_control);
virtual void CollectCalculateRightHandSide(ProcessInfo& r_process_info);
virtual void FinalCalculateRightHandSide(ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity);
virtual void InitializeForceComputation(ProcessInfo& r_process_info);
virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info);
virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info);
virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info);
virtual void GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& r_process_info );
virtual void FinalizeSolutionStep(ProcessInfo& r_process_info);
virtual void SymmetrizeStressTensor();
virtual void CorrectRepresentativeVolume(double& rRepresentative_Volume/*, bool& is_smaller_than_sphere*/);
virtual void ComputeReactions();
virtual void PrepareForPrinting(ProcessInfo& r_process_info);
virtual void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info);
virtual void CalculateMaxBallToBallIndentation(double& rCurrentMaxIndentation);
virtual void CalculateMaxBallToFaceIndentation(double& rCurrentMaxIndentation);
virtual double CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info);

virtual void ComputeConditionRelativeData(int rigid_neighbour_index,
                                          DEMWall* const wall,
                                            double LocalCoordSystem[3][3],
                                            double& DistPToB,
                                            array_1d<double, 4>& Weight,
                                            array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                            array_1d<double, 3>& wall_velocity_at_contact_point,
                                            int& ContactType
                                            );

virtual void RenewData();
virtual void SendForcesToFEM();
int   GetClusterId();
void  SetClusterId(const int Id);

virtual double GetRadius();
virtual void   SetRadius(double radius);
virtual void   SetRadius();
virtual double CalculateVolume();
virtual double GetInteractionRadius();
virtual void SetInteractionRadius(const double radius);
virtual double GetSearchRadius();
virtual double GetSearchRadiusWithFem();
virtual void SetSearchRadius(const double radius);
virtual void SetSearchRadiusWithFem(const double radius);
virtual double GetMass();
virtual void   SetMass(double real_mass);
virtual double   CalculateMomentOfInertia();
virtual double GetYoung();
void   SetYoungFromProperties(double* young);
virtual double GetRollingFriction();
void   SetRollingFrictionFromProperties(double* rolling_friction);
virtual double GetPoisson();
void   SetPoissonFromProperties(double* poisson);
virtual double GetTgOfFrictionAngle();
void   SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle);
virtual double GetCoefficientOfRestitution();
void   SetCoefficientOfRestitutionFromProperties(double* coefficient_of_restitution);
virtual double GetLnOfRestitCoeff();
void   SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff);
virtual double GetDensity();
void   SetDensityFromProperties(double* density);
virtual int    GetParticleMaterial();
void   SetParticleMaterialFromProperties(int* particle_material);
virtual double GetParticleCohesion();
void   SetParticleCohesionFromProperties(double* particle_cohesion);
virtual double GetParticleKNormal();
void   SetParticleKNormalFromProperties(double* particle_k_normal);
virtual double GetParticleKTangential();
void   SetParticleKTangentialFromProperties(double* particle_k_tangential);

//Conical damage
virtual double GetParticleContactRadius();
void   SetParticleContactRadiusFromProperties(double* particle_contact_radius);
virtual double GetParticleMaxStress();
void   SetParticleMaxStressFromProperties(double* particle_max_stress);
virtual double GetParticleAlpha();
void   SetParticleAlphaFromProperties(double* particle_alpha);
virtual double GetParticleGamma();
void   SetParticleGammaFromProperties(double* particle_gamma);

array_1d<double, 3>& GetForce();
double mElasticEnergy;
virtual double& GetElasticEnergy();
double mInelasticFrictionalEnergy;
virtual double& GetInelasticFrictionalEnergy();
double mInelasticViscodampingEnergy;
virtual double& GetInelasticViscodampingEnergy();

PropertiesProxy* GetFastProperties();
void   SetFastProperties(PropertiesProxy* pProps);
void   SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies);

double SlowGetYoung();
double SlowGetRollingFriction();
double SlowGetPoisson();
double SlowGetTgOfFrictionAngle();
double SlowGetCoefficientOfRestitution();
double SlowGetDensity();
double SlowGetParticleCohesion();
int    SlowGetParticleMaterial();

double GetBoundDeltaDispSq();

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

virtual void ComputeNewNeighboursHistoricalData(boost::numeric::ublas::vector<int>& mTempNeighboursIds, std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces);

virtual void ComputeNewRigidFaceNeighboursHistoricalData();
std::vector<SphericParticle*>     mNeighbourElements;
std::vector<DEMWall*>             mNeighbourRigidFaces;
std::vector<DEMWall*>             mNeighbourPotentialRigidFaces;

std::vector<array_1d<double, 4> > mContactConditionWeights;
std::vector<int>                  mContactConditionContactTypes;
std::vector< array_1d<double,3> > mConditionContactPoints;

std::vector<array_1d<double, 3> > mNeighbourRigidFacesTotalContactForce;
std::vector<array_1d<double, 3> > mNeighbourRigidFacesElasticContactForce;
std::vector<array_1d<double, 3> > mNeighbourElasticContactForces;
std::vector<array_1d<double, 3> > mNeighbourElasticExtraContactForces;

virtual void ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, const ProcessInfo& r_process_info, const array_1d<double,3>& gravity);
virtual array_1d<double,3> ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info);

virtual void MemberDeclarationFirstStep(const ProcessInfo& r_process_info);

array_1d<double, 3> mContactMoment; //SLS

Matrix* mStressTensor;
Matrix* mSymmStressTensor;
double mPartialRepresentativeVolume;

std::vector<int> mFemOldNeighbourIds;

protected:

SphericParticle();

virtual void ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rElasticForce,
                                                array_1d<double, 3>& rContactForce,
                                                array_1d<double, 3>& InitialRotaMoment,
                                                array_1d<double, 3>& rigid_element_force,
                                                ProcessInfo& r_process_info,
                                                double mTimeStep,
                                                int search_control) final;


virtual void InitializeSolutionStep(ProcessInfo& r_process_info);

virtual void CalculateMomentum(array_1d<double, 3>& rMomentum);
virtual void CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum);
virtual void ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                     array_1d<double, 3>& rContactForce,
                                     array_1d<double, 3>& InitialRotaMoment,
                                     ProcessInfo& r_process_info,
                                     const double dt,
                                     const bool multi_stage_RHS);

virtual void EvaluateDeltaDisplacement(double DeltDisp[3],
                    double RelVel[3],
                    double LocalCoordSystem[3][3],
                    double OldLocalCoordSystem[3][3],
                    const array_1d<double, 3> &other_to_me_vect,
                    const array_1d<double, 3> &vel,
                    const array_1d<double, 3> &delta_displ,
                    SphericParticle* neighbour_iterator,
                    double& distance);

virtual void RelativeDisplacementAndVelocityOfContactPointDueToRotation(const double indentation,
                    double DeltDesp[3],
                    double RelVel[3],
                    double OldLocalCoordSystem[3][3],
                    const double &other_radius,
                    const double &dt,
                    const array_1d<double, 3> &angl_vel,
                    SphericParticle* neighbour_iterator);

virtual void RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                            double DeltDisp[3], //IN GLOBAL AXES
                                                                            double RelVel[3], //IN GLOBAL AXES
                                                                            double OldLocalCoordSystem[3][3],
                                                                            double LocalCoordSystem[3][3],
                                                                            SphericParticle* neighbour_iterator);

virtual void RelativeDisplacementAndVelocityOfContactPointDueToRotationMatrix(double DeltDisp[3],
                                                                                double RelVel[3],
                                                                                double OldLocalCoordSystem[3][3],
                                                                                const double& other_radius,
                                                                                const double& dt,
                                                                                const array_1d<double, 3>& ang_vel,
                                                                                SphericParticle* p_neighbour);

virtual void ComputeMoments(double normalLocalElasticContactForce,
                      double GlobalElasticContactForces[3],
                      array_1d<double, 3>& rInitialRotaMoment,
                      double LocalCoordSystem_2[3],
                      SphericParticle* neighbour_iterator,
                      double indentation,
                      bool wall=false) final;

virtual double GetInitialDeltaWithFEM(int index);

virtual void ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces);

virtual void AddUpForcesAndProject(double OldCoordSystem[3][3],
                    double LocalCoordSystem[3][3],
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double LocalElasticExtraContactForce[3],                    
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double GlobalElasticExtraContactForce[3],
                    double TotalGlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    const double cohesive_force,
                    array_1d<double, 3>& other_ball_to_ball_forces,
                    array_1d<double, 3>& rElasticForce,
                    array_1d<double, 3>& rContactForce,
                    const unsigned int i_neighbour_count,
                    ProcessInfo& r_process_info) final;

virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                    double LocalContactForce[3],
                    double LocalElasticContactForce[3],
                    double GlobalContactForce[3],
                    double GlobalElasticContactForce[3],
                    double ViscoDampingLocalContactForce[3],
                    const double cohesive_force,
                    array_1d<double, 3>& rElasticForce,
                    array_1d<double, 3>& rContactForce,
                    const unsigned int iRigidFaceNeighbour) final;

virtual void AddUpMomentsAndProject(double LocalCoordSystem[3][3],
                                    double ElasticLocalRotationalMoment[3],
                                    double ViscoLocalRotationalMoment[3]) final;

virtual void ComputeWear(double LocalCoordSystem[3][3], array_1d<double, 3>& vel, double tangential_vel[3],
                         double mTimeStep, double density, bool sliding, double inverse_of_volume,
                         double LocalElasticContactForce, DEMWall* cast_neighbour);

virtual void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);

virtual void AddNeighbourContributionToStressTensor(const double GlobalElasticContactForce[3],
                                                    const double other_to_me_vect[3],
                                                    const double distance,
                                                    const double radius_sum,
                                                    SphericParticle* element);

virtual void AddWallContributionToStressTensor(const double GlobalElasticContactForce[3],
                                                const double other_to_me_vect[3],
                                                const double distance,
                                                const double contact_area);

virtual void RotateOldContactForces(const double LocalCoordSystem[3][3], const double OldLocalCoordSystem[3][3], array_1d<double, 3>& mNeighbourElasticContactForces) final;

DEMDiscontinuumConstitutiveLaw::Pointer mDiscontinuumConstitutiveLaw;

//const int mParticleId; // (NOT YET ACTIVE!!) Identifies the particle biunivocally if it has been properly created (i.e., a non-repeated NewId is passed to the constructor)

//const double* mSearchControl;

//array_1d<double, 3> mContactForce; //SLS

double mRadius;
double mSearchRadius;
double mSearchRadiusWithFem;
double mRealMass;
PropertiesProxy* mFastProperties;
int mClusterId;
double mBoundDeltaDispSq;

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );
    rSerializer.save("mRadius",mRadius);
    rSerializer.save("mSearchRadius", mSearchRadius);
    rSerializer.save("mSearchRadiusWithFem", mSearchRadiusWithFem);
    rSerializer.save("mRealMass",mRealMass);
    rSerializer.save("mClusterId",mClusterId);
    rSerializer.save("mBoundDeltaDispSq",mBoundDeltaDispSq);
    rSerializer.save("HasStressTensor", (int)this->Is(DEMFlags::HAS_STRESS_TENSOR));
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        rSerializer.save("mSymmStressTensor", mSymmStressTensor);
    }
}

virtual void load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
    rSerializer.load("mRadius",mRadius);
    rSerializer.load("mSearchRadius", mSearchRadius);
    rSerializer.load("mSearchRadiusWithFem", mSearchRadiusWithFem);
    rSerializer.load("mRealMass",mRealMass);
    rSerializer.load("mClusterId",mClusterId);
    rSerializer.load("mBoundDeltaDispSq",mBoundDeltaDispSq); 
    int aux_int=0;
    rSerializer.load("HasStressTensor", aux_int);
    if(aux_int) this->Set(DEMFlags::HAS_STRESS_TENSOR, true);
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        mStressTensor  = new Matrix(3,3);
        *mStressTensor = ZeroMatrix(3,3);
        mSymmStressTensor  = new Matrix(3,3);
        *mSymmStressTensor = ZeroMatrix(3,3);
        rSerializer.load("mSymmStressTensor", mSymmStressTensor);
    }
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

}  // namespace Kratos.

#endif // KRATOS_SPHERIC_PARTICLE_H_INCLUDED  defined
