//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_SPHERIC_PARTICLE_H_INCLUDED)
#define  KRATOS_SPHERIC_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <map>

// Project includes
#include "includes/define.h"
#include "discrete_element.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "custom_constitutive/DEM_rolling_friction_model.h"
#include "custom_constitutive/DEM_global_damping.h"
#include "custom_conditions/RigidFace.h"
#include "custom_conditions/dem_wall.h"
#include "custom_strategies/schemes/dem_integration_scheme.h"
#include "includes/kratos_export_api.h"
#include "custom_utilities/properties_proxies.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

class DEMWall;

class KRATOS_API(DEM_APPLICATION) SphericParticle : public DiscreteElement
{

public:

/// Pointer definition of SphericParticle
KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SphericParticle);

typedef GlobalPointersVector<Condition> ConditionWeakVectorType;
typedef GlobalPointersVector<Condition >::iterator ConditionWeakIteratorType;

typedef GlobalPointersVector<Element> ParticleWeakVectorType;
typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
/// Default constructor.
ModelPart* mpInlet;
SphericParticle();
SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry );
SphericParticle( IndexType NewId, NodesArrayType const& ThisNodes);
SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

/// Destructor.
virtual ~SphericParticle();

SphericParticle& operator=(const SphericParticle& rOther);

class ParticleDataBuffer
{
public:
    ParticleDataBuffer(SphericParticle* p_this_particle): mpThisParticle(p_this_particle)
    {}

    virtual ~ParticleDataBuffer(){}

virtual bool SetNextNeighbourOrExit(int& i)
{
    if (i < int(mpThisParticle->mNeighbourElements.size())){
        SetCurrentNeighbour(mpThisParticle->mNeighbourElements[i]);
        mpOtherParticleNode = &(mpOtherParticle->GetGeometry()[0]);
        return true;
    }

    else { // other_neighbour is nullified upon exiting loop
        mpOtherParticle = NULL;
        mpOtherParticleNode = NULL;
        return false;
    }
}

void SetCurrentNeighbour(SphericParticle* p_neighbour)
{
    mpOtherParticle = p_neighbour;
}

void SetBoundingBox(const bool periodicity, const array_1d<double, 3> domain_min, const array_1d<double, 3> domain_max)
{
    mDomainIsPeriodic = periodicity;
    mDomainMin = domain_min;
    mDomainMax = domain_max;
}

bool mMultiStageRHS;
bool mDomainIsPeriodic;
double mDistance;
double mRadiusSum;
double mDt;
double mTime;
double mOtherRadius;
double mIndentation;
double mMyCoors[3];
double mOtherCoors[3];
double mLocalRelVel[3];
array_1d<double, 3> mOtherToMeVector;
array_1d<double, 3> mDomainMin;
array_1d<double, 3> mDomainMax;
SphericParticle* mpThisParticle;
SphericParticle* mpOtherParticle;
Node* mpOtherParticleNode;
DEMWall* mpOtherRigidFace;
double mLocalCoordSystem[3][3];
double mOldLocalCoordSystem[3][3];

std::vector<DEMWall*> mNeighbourRigidFaces; // why repeated? it is in the sphere as well!

};

typedef std::unique_ptr<ParticleDataBuffer> BufferPointerType;

virtual std::unique_ptr<ParticleDataBuffer> CreateParticleDataBuffer(SphericParticle* p_this_particle)
{
    return std::unique_ptr<ParticleDataBuffer>(new ParticleDataBuffer(p_this_particle));
}

void TransformNeighbourCoorsToClosestInPeriodicDomain(ParticleDataBuffer & data_buffer);
void TransformNeighbourCoorsToClosestInPeriodicDomain(ParticleDataBuffer & data_buffer,
                                                    const array_1d<double, 3>& coors,
                                                    array_1d<double, 3>& neighbour_coors);
void TransformNeighbourCoorsToClosestInPeriodicDomain(const ProcessInfo& r_process_info,
                                                    const double coors[3],
                                                    double neighbour_coors[3]);

virtual bool CalculateRelativePositionsOrSkipContact(ParticleDataBuffer & data_buffer);

void Initialize(const ProcessInfo& r_process_info) override;
virtual void MemberDeclarationFirstStep(const ProcessInfo& r_process_info);
using DiscreteElement::CalculateRightHandSide; //To avoid Clang Warning. We tell the compiler that we are aware of the existence of this function, but we overload it still.
virtual void CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity);
virtual void FirstCalculateRightHandSide(const ProcessInfo& r_process_info, double dt);
virtual void CollectCalculateRightHandSide(const ProcessInfo& r_process_info);
virtual void FinalCalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double,3>& gravity);
virtual void InitializeForceComputation(const ProcessInfo& r_process_info);
virtual void FinalizeForceComputation(ParticleDataBuffer & data_buffer){}
virtual void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& r_process_info) const override;
virtual void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& r_process_info) override;
virtual void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& r_process_info) override;
virtual void GetDofList( DofsVectorType& ElementalDofList, const ProcessInfo& r_process_info ) const override;
virtual void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids, std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces);
virtual void ComputeNewRigidFaceNeighboursHistoricalData();
virtual void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;
virtual void InitializeSolutionStep(const ProcessInfo& r_process_info) override;
virtual void FinalizeStressTensor(const ProcessInfo& r_process_info, double& rRepresentative_Volume){};
virtual void SymmetrizeStressTensor();
virtual void ComputeStrainTensor(const ProcessInfo& r_process_info);
virtual void ComputeDifferentialStrainTensor(const ProcessInfo& r_process_info);
virtual void SymmetrizeDifferentialStrainTensor();
virtual void CorrectRepresentativeVolume(double& rRepresentative_Volume/*, bool& is_smaller_than_sphere*/);
virtual void ComputeReactions();
virtual void PrepareForPrinting(const ProcessInfo& r_process_info);
virtual void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;
virtual void Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& r_process_info) override;
virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info) override;
virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& r_process_info) override;
virtual void CalculateMaxBallToBallIndentation(double& rCurrentMaxIndentation, const ProcessInfo& r_process_info);
virtual void CalculateMaxBallToFaceIndentation(double& rCurrentMaxIndentation);

virtual void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag);
virtual void SetIntegrationScheme(DEMIntegrationScheme::Pointer& translational_integration_scheme, DEMIntegrationScheme::Pointer& rotational_integration_scheme);
virtual bool SwapIntegrationSchemeToGluedToWall(Condition* p_wall);
virtual DEMIntegrationScheme& GetTranslationalIntegrationScheme() { return *mpTranslationalIntegrationScheme; }
virtual DEMIntegrationScheme& GetRotationalIntegrationScheme() { return *mpRotationalIntegrationScheme; }

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
double GetInitializationTime() const;
double GetProgrammedDestructionTime() const;
void SetProgrammedDestructionTime(const double time);
virtual double GetRadius();
virtual void   SetRadius(double radius);
virtual void   SetRadius();
virtual void   SetRadius(bool is_radius_expansion, double radius_expansion_rate, double radius_multiplier_max, double radius_multiplier, double radius_multiplier_old);
virtual double CalculateVolume();
virtual double GetInteractionRadius(const int radius_index = 0);
virtual void SetInteractionRadius(const double radius, const int radius_index = 0);
virtual double GetSearchRadius();
virtual void SetDefaultRadiiHierarchy(const double radius);
virtual void SetSearchRadius(const double radius);
virtual double GetMass();
virtual void   SetMass(double real_mass);
virtual double   CalculateMomentOfInertia();
virtual double GetYoung();
void   SetYoungFromProperties(double* young);
virtual double GetPoisson();
void   SetPoissonFromProperties(double* poisson);
virtual double GetDensity();
void   SetDensityFromProperties(double* density);
virtual int    GetParticleMaterial();
void   SetParticleMaterialFromProperties(int* particle_material);


array_1d<double, 3>& GetForce();

virtual double& GetElasticEnergy();
virtual double& GetInelasticFrictionalEnergy();
virtual double& GetInelasticViscodampingEnergy();
virtual double& GetInelasticRollingResistanceEnergy();
virtual double& GetMaxNormalBallToBallForceTimesRadius();

PropertiesProxy* GetFastProperties();
void   SetFastProperties(PropertiesProxy* pProps);
void   SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies);

double SlowGetYoung() const;
double SlowGetPoisson() const;
double SlowGetDensity() const;
int    SlowGetParticleMaterial() const;

/// Turn back information as a string.
virtual std::string Info() const override
{
std::stringstream buffer;
buffer << "SphericParticle" ;
return buffer.str();
}

/// Print information about this object.
virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SphericParticle";}

/// Print object's data.
//virtual void PrintData(std::ostream& rOStream) const override {}

double mElasticEnergy;
double mInelasticFrictionalEnergy;
double mInelasticViscodampingEnergy;
double mInelasticRollingResistanceEnergy;
double mPartialRepresentativeVolume;
double mMaxNormalBallToBallForceTimesRadius;

std::vector<ParticleContactElement*> mBondElements;
std::vector<SphericParticle*>     mNeighbourElements;
std::vector<int>                  mContactingNeighbourIds;
std::vector<int>                  mContactingFaceNeighbourIds;
std::vector<DEMWall*>             mNeighbourRigidFaces;
std::vector<DEMWall*>             mNeighbourNonContactRigidFaces;
std::vector<DEMWall*>             mNeighbourPotentialRigidFaces;

std::vector<array_1d<double, 4> > mContactConditionWeights;
std::vector<int>                  mContactConditionContactTypes;
std::vector< array_1d<double,3> > mConditionContactPoints;

std::vector<array_1d<double, 3> > mNeighbourRigidFacesTotalContactForce;
std::vector<array_1d<double, 3> > mNeighbourRigidFacesElasticContactForce;
std::vector<array_1d<double, 3> > mNeighbourElasticContactForces;
std::vector<array_1d<double, 3> > mNeighbourElasticExtraContactForces;
std::vector<array_1d<double, 3> > mNeighbourRollingFrictionMoments;
std::vector<int> mFemOldNeighbourIds;
array_1d<double, 3> mContactMoment; //SLS

BoundedMatrix<double, 3, 3>* mStressTensor;
BoundedMatrix<double, 3, 3>* mSymmStressTensor;
BoundedMatrix<double, 3, 3>* mStrainTensor;
BoundedMatrix<double, 3, 3>* mDifferentialStrainTensor;

virtual void ComputeAdditionalForces(array_1d<double, 3>& externally_applied_force, array_1d<double, 3>& externally_applied_moment, const ProcessInfo& r_process_info, const array_1d<double,3>& gravity);
virtual array_1d<double,3> ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info);
virtual void CalculateOnContactElements(size_t i_neighbour_count, double LocalContactForce[3], double GlobalContactForce[3]);

std::unique_ptr<DEMDiscontinuumConstitutiveLaw> pCloneDiscontinuumConstitutiveLawWithNeighbour(SphericParticle* neighbour);

std::unique_ptr<DEMDiscontinuumConstitutiveLaw> pCloneDiscontinuumConstitutiveLawWithFEMNeighbour(Condition* neighbour);

std::unique_ptr<DEMRollingFrictionModel> pCloneRollingFrictionModel(SphericParticle* element);

std::unique_ptr<DEMRollingFrictionModel> pCloneRollingFrictionModelWithNeighbour(SphericParticle* neighbour);

std::unique_ptr<DEMRollingFrictionModel> pCloneRollingFrictionModelWithFEMNeighbour(Condition* neighbour);

std::unique_ptr<DEMGlobalDampingModel> pCloneGlobalDampingModel(SphericParticle* element);

protected:

virtual void ComputeBallToRigidFaceContactForceAndMoment(ParticleDataBuffer & data_buffer,
                                                array_1d<double, 3>& rElasticForce,
                                                array_1d<double, 3>& rContactForce,
                                                array_1d<double, 3>& rigid_element_force,
                                                const ProcessInfo& r_process_info) ;

virtual void CalculateMomentum(array_1d<double, 3>& rMomentum);

virtual void CalculateLocalAngularMomentum(array_1d<double, 3>& rAngularMomentum);

virtual void ComputeBallToBallContactForceAndMoment(ParticleDataBuffer & data_buffer,
                                        const ProcessInfo& r_process_info,
                                        array_1d<double, 3>& rElasticForce,
                                        array_1d<double, 3>& rContactForce);

virtual void EvaluateDeltaDisplacement(ParticleDataBuffer & data_buffer,
                                    double DeltDisp[3],
                                    double RelVel[3],
                                    double LocalCoordSystem[3][3],
                                    double OldLocalCoordSystem[3][3],
                                    const array_1d<double, 3> &vel,
                                    const array_1d<double, 3> &delta_displ);

virtual void RelativeDisplacementAndVelocityOfContactPointDueToRotation(const double indentation,
                                                                        double DeltDesp[3],
                                                                        double RelVel[3],
                                                                        const double OldLocalCoordSystem[3][3],
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
                                                                            const double OldLocalCoordSystem[3][3],
                                                                            const double& other_radius,
                                                                            const double& dt,
                                                                            const array_1d<double, 3>& ang_vel,
                                                                            SphericParticle* p_neighbour);

virtual void RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(double DeltDesp[3],
                                                                                double RelVel[3],
                                                                                const double OldLocalCoordSystem[3][3],
                                                                                const double &other_radius,
                                                                                const double &dt,
                                                                                const array_1d<double, 3> &angl_vel,
                                                                                SphericParticle* neighbour_iterator,
                                                                                ParticleDataBuffer & data_buffer);

virtual void ComputeMoments(double normalLocalContactForce,
                            double GlobalContactForce[3],
                            double LocalCoordSystem_2[3],
                            SphericParticle* neighbour_iterator,
                            double indentation,
                            unsigned int i);

virtual void ComputeMomentsWithWalls(double normalLocalContactForce,
                            double GlobalContactForce[3],
                            double LocalCoordSystem_2[3],
                            Condition* wall,
                            double indentation,
                            unsigned int i);

virtual double GetInitialDeltaWithFEM(int index);

virtual void ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces);

virtual void StoreBallToBallContactInfo(const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding);

virtual void StoreBallToRigidFaceContactInfo(const ProcessInfo& r_process_info, SphericParticle::ParticleDataBuffer& data_buffer, double GlobalContactForceTotal[3], double LocalContactForceTotal[3], double LocalContactForceDamping[3], bool sliding);

virtual void EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                            const ProcessInfo& r_process_info,
                                                            double LocalElasticContactForce[3],
                                                            double DeltDisp[3],
                                                            double LocalDeltDisp[3],
                                                            double RelVel[3],
                                                            double indentation,
                                                            double ViscoDampingLocalContactForce[3],
                                                            double& cohesive_force,
                                                            SphericParticle* element2,
                                                            bool& sliding,
                                                            double LocalCoordSystem[3][3],
                                                            double OldLocalCoordSystem[3][3],
                                                            array_1d<double, 3>& neighbour_elastic_contact_force);


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
                                const ProcessInfo& r_process_info) final;

virtual void AddUpFEMForcesAndProject(double LocalCoordSystem[3][3],
                                    double LocalContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double GlobalContactForce[3],
                                    double GlobalElasticContactForce[3],
                                    double ViscoDampingLocalContactForce[3],
                                    const double cohesive_force,
                                    array_1d<double, 3>& rElasticForce,
                                    array_1d<double, 3>& rContactForce,
                                    array_1d<double, 3>& elastic_force_backup,
                                    array_1d<double, 3>& total_force_backup) final;

virtual void AddUpMomentsAndProject(double LocalCoordSystem[3][3],
                                    double ElasticLocalRotationalMoment[3],
                                    double ViscoLocalRotationalMoment[3]) final;

virtual void ComputeWear(double LocalRelVel[3],
                        double mTimeStep,
                        bool sliding,
                        double inverse_of_volume,
                        double LocalElasticContactForce,
                        DEMWall* cast_neighbour);

virtual void AdditionalCalculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info);

virtual void AddNeighbourContributionToStressTensor(const ProcessInfo& r_process_info,
                                                    const double GlobalContactForce[3],
                                                    const double other_to_me_vect[3],
                                                    const double distance,
                                                    const double radius_sum,
                                                    SphericParticle* element);

virtual void AddWallContributionToStressTensor(const double GlobalContactForce[3],
                                               const double other_to_me_vect[3],
                                               const double distance,
                                               const double contact_area);

virtual void RotateOldContactForces(const double LocalCoordSystem[3][3], const double OldLocalCoordSystem[3][3], array_1d<double, 3>& mNeighbourElasticContactForces) final;
virtual void ApplyGlobalDampingToContactForcesAndMoments(array_1d<double,3>& total_forces, array_1d<double,3>& total_moment);

std::unique_ptr<DEMDiscontinuumConstitutiveLaw> mDiscontinuumConstitutiveLaw;
std::unique_ptr<DEMRollingFrictionModel> mRollingFrictionModel;
std::unique_ptr<DEMRollingFrictionModel> mRollingFrictionModelNeighbour;
std::unique_ptr<DEMGlobalDampingModel> mGlobalDampingModel;

double mInitializationTime;
double mIndentationInitialOption;
std::map<int, double> mIndentationInitial;
std::map<int, double> mIndentationInitialWall;
double mProgrammedDestructionTime=-1.0; // set to a negative value, so that when marked TO_ERASE, elimination is by default.
double mRadius;
double mInitialRadius;
double mSearchRadius;
double mRealMass;
PropertiesProxy* mFastProperties;
int mClusterId;
DEMIntegrationScheme* mpTranslationalIntegrationScheme;
DEMIntegrationScheme* mpRotationalIntegrationScheme;
double mBondedScalingFactor[3] = {0.0};

private:

friend class Serializer;

virtual void save(Serializer& rSerializer) const override
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DiscreteElement );

    // public members
    rSerializer.save("mpInlet", mpInlet);
    rSerializer.save("mElasticEnergy", mElasticEnergy);
    rSerializer.save("mInelasticFrictionalEnergy", mInelasticFrictionalEnergy);
    rSerializer.save("mInelasticViscodampingEnergy", mInelasticViscodampingEnergy);
    rSerializer.save("mInelasticRollingResistanceEnergy", mInelasticRollingResistanceEnergy);
    rSerializer.save("mPartialRepresentativeVolume", mPartialRepresentativeVolume);
    rSerializer.save("mMaxNormalBallToBallForceTimesRadius", mMaxNormalBallToBallForceTimesRadius);
    rSerializer.save("mBondElements", mBondElements);
    rSerializer.save("mNeighbourElements", mNeighbourElements);
    rSerializer.save("mContactingNeighbourIds", mContactingNeighbourIds);
    rSerializer.save("mContactingFaceNeighbourIds", mContactingFaceNeighbourIds);
    rSerializer.save("mNeighbourRigidFaces", mNeighbourRigidFaces);
    rSerializer.save("mNeighbourNonContactRigidFaces", mNeighbourNonContactRigidFaces);
    rSerializer.save("mNeighbourPotentialRigidFaces", mNeighbourPotentialRigidFaces);
    rSerializer.save("mContactConditionWeights", mContactConditionWeights);
    rSerializer.save("mContactConditionContactTypes", mContactConditionContactTypes);
    rSerializer.save("mConditionContactPoints", mConditionContactPoints);
    rSerializer.save("mNeighbourRigidFacesTotalContactForce", mNeighbourRigidFacesTotalContactForce);
    rSerializer.save("mNeighbourRigidFacesElasticContactForce", mNeighbourRigidFacesElasticContactForce);
    rSerializer.save("mNeighbourElasticContactForces", mNeighbourElasticContactForces);
    rSerializer.save("mNeighbourElasticExtraContactForces", mNeighbourElasticExtraContactForces);
    rSerializer.save("mNeighbourRollingFrictionMoments", mNeighbourRollingFrictionMoments);
    rSerializer.save("mFemOldNeighbourIds", mFemOldNeighbourIds);
    rSerializer.save("mContactMoment", mContactMoment);

    rSerializer.save("HasStressTensor", (int)this->Is(DEMFlags::HAS_STRESS_TENSOR));
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        rSerializer.save("mStressTensor", mStressTensor);
        rSerializer.save("mSymmStressTensor", mSymmStressTensor);
        rSerializer.save("mStrainTensor", mStrainTensor);
        rSerializer.save("mDifferentialStrainTensor", mDifferentialStrainTensor);        
    }

    // protected members
    rSerializer.save("mRadius", mRadius);
    rSerializer.save("mSearchRadius", mSearchRadius);
    rSerializer.save("mRealMass", mRealMass);
    rSerializer.save("mClusterId", mClusterId);
}

virtual void load(Serializer& rSerializer) override
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DiscreteElement );
    // public members
    rSerializer.load("mpInlet", mpInlet);
    rSerializer.load("mElasticEnergy", mElasticEnergy);
    rSerializer.load("mInelasticFrictionalEnergy", mInelasticFrictionalEnergy);
    rSerializer.load("mInelasticViscodampingEnergy", mInelasticViscodampingEnergy);
    rSerializer.load("mInelasticRollingResistanceEnergy", mInelasticRollingResistanceEnergy);
    rSerializer.load("mPartialRepresentativeVolume", mPartialRepresentativeVolume);
    rSerializer.load("mMaxNormalBallToBallForceTimesRadius", mMaxNormalBallToBallForceTimesRadius);
    rSerializer.load("mBondElements", mBondElements);
    rSerializer.load("mNeighbourElements", mNeighbourElements);
    rSerializer.load("mContactingNeighbourIds", mContactingNeighbourIds);
    rSerializer.load("mContactingFaceNeighbourIds", mContactingFaceNeighbourIds);
    rSerializer.load("mNeighbourRigidFaces", mNeighbourRigidFaces);
    rSerializer.load("mNeighbourNonContactRigidFaces", mNeighbourNonContactRigidFaces);
    rSerializer.load("mNeighbourPotentialRigidFaces", mNeighbourPotentialRigidFaces);
    rSerializer.load("mContactConditionWeights", mContactConditionWeights);
    rSerializer.load("mContactConditionContactTypes", mContactConditionContactTypes);
    rSerializer.load("mConditionContactPoints", mConditionContactPoints);
    rSerializer.load("mNeighbourRigidFacesTotalContactForce", mNeighbourRigidFacesTotalContactForce);
    rSerializer.load("mNeighbourRigidFacesElasticContactForce", mNeighbourRigidFacesElasticContactForce);
    rSerializer.load("mNeighbourElasticContactForces", mNeighbourElasticContactForces);
    rSerializer.load("mNeighbourElasticExtraContactForces", mNeighbourElasticExtraContactForces);
    rSerializer.load("mNeighbourRollingFrictionMoments", mNeighbourRollingFrictionMoments);
    rSerializer.load("mFemOldNeighbourIds", mFemOldNeighbourIds);
    rSerializer.load("mContactMoment", mContactMoment);

    int aux_int=0;
    rSerializer.load("HasStressTensor", aux_int);
    if(aux_int) this->Set(DEMFlags::HAS_STRESS_TENSOR, true);
    if (this->Is(DEMFlags::HAS_STRESS_TENSOR)){
        mStressTensor = new BoundedMatrix<double, 3, 3>(3,3);
        *mStressTensor = ZeroMatrix(3,3);
        mSymmStressTensor = new BoundedMatrix<double, 3, 3>(3,3);
        *mSymmStressTensor = ZeroMatrix(3,3);
        rSerializer.load("mStressTensor", mStressTensor);
        rSerializer.load("mSymmStressTensor", mSymmStressTensor);
        mStrainTensor = new BoundedMatrix<double, 3, 3>(3,3);
        *mStrainTensor = ZeroMatrix(3,3);
        rSerializer.load("mStrainTensor", mStrainTensor);
        mDifferentialStrainTensor = new BoundedMatrix<double, 3, 3>(3,3);
        *mDifferentialStrainTensor = ZeroMatrix(3,3);
        rSerializer.load("mDifferentialStrainTensor", mDifferentialStrainTensor);
    }

    // protected members
    rSerializer.load("mRadius", mRadius);
    rSerializer.load("mSearchRadius", mSearchRadius);
    rSerializer.load("mRealMass", mRealMass);
    rSerializer.load("mClusterId", mClusterId);
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
