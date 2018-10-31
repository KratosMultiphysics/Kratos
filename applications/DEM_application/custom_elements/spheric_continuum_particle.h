//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

#if !defined(KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"
#include "utilities/timer.h"
#include "Particle_Contact_Element.h"
#include "../custom_constitutive/DEM_continuum_constitutive_law.h"
#include "containers/vector_component_adaptor.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) SphericContinuumParticle : public SphericParticle
    {
    public:

        /// Pointer definition of SphericContinuumParticle
        KRATOS_CLASS_POINTER_DEFINITION(SphericContinuumParticle);

        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

        /// Default constructor
        SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry);
        SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes);
        SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~SphericContinuumParticle();


        class ParticleDataBuffer : public SphericParticle::ParticleDataBuffer
        {
        public:
            ParticleDataBuffer(SphericParticle* p_this_particle): SphericParticle::ParticleDataBuffer(p_this_particle){}

            virtual ~ParticleDataBuffer(){}

            bool SetNextNeighbourOrExit(int& i) override
            {
                while (i < int(mpThisParticle->mNeighbourElements.size()) && (mpThisParticle->mNeighbourElements[i]==NULL)){
                    i++;
                }

                if (i < int(mpThisParticle->mNeighbourElements.size())) {
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
        };

        std::unique_ptr<SphericParticle::ParticleDataBuffer> CreateParticleDataBuffer(SphericParticle* p_this_particle) override
        {
            return std::unique_ptr<SphericParticle::ParticleDataBuffer>(new ParticleDataBuffer(p_this_particle));
        }

        void SetInitialSphereContacts(ProcessInfo& r_process_info);
        void SetInitialFemContacts();
        void CreateContinuumConstitutiveLaws();
        //void InitializeSolutionStep(ProcessInfo& r_process_info) override;
        void FinalizeSolutionStep(ProcessInfo& r_process_info) override;
        void GetStressTensorFromNeighbourStep1();
        void GetStressTensorFromNeighbourStep2();
        void GetStressTensorFromNeighbourStep3();
        //void CorrectRepresentativeVolume(double& rRepresentative_Volume, bool& is_smaller_than_sphere) override;

        void Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) override;

        void ReorderAndRecoverInitialPositionsAndFilter(std::vector<SphericParticle*>& mTempNeighbourElements);
        void ReorderFEMneighbours();
        virtual void UpdateContinuumNeighboursVector(ProcessInfo& r_process_info);
        virtual void ComputeForceWithNeighbourFinalOperations();

        virtual double CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info) override;
        virtual double CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info);
        virtual bool OverlappedParticleRemoval();
        virtual void CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info);
        virtual void CalculateOnContactElements(size_t i_neighbour_count, double LocalElasticContactForce[3],
                                                double contact_sigma, double contact_tau, double failure_criterion_state, double acumulated_damage, int time_steps);


        virtual void FilterNonSignificantDisplacements(double DeltDisp[3], //IN GLOBAL AXES
                                                       double RelVel[3], //IN GLOBAL AXES
                                                       double& indentation);



        virtual void ContactAreaWeighting();
        virtual double EffectiveVolumeRadius();
        virtual double GetInitialDelta(int index);
        virtual bool IsSkin() { return (bool)*mSkinSphere; }
        void MarkNewSkinParticlesDueToBreakage();

        /// Turn back information as a string
        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "SphericCosntinuumParticle" ;
            return buffer.str();
        }

        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "SphericContinuumParticle";}

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const override {}

        //member variables DEM_CONTINUUM
        int mContinuumGroup;
        std::vector<int> mIniNeighbourIds;
        std::vector<int> mIniNeighbourFailureId;
        std::vector<double> mIniNeighbourDelta;

        unsigned int mContinuumInitialNeighborsSize;
        std::vector<Kratos::DEMContinuumConstitutiveLaw::Pointer> mContinuumConstitutiveLawArray;

    protected:

        SphericContinuumParticle();

        virtual void Initialize(const ProcessInfo& r_process_info) override;
        virtual double GetInitialDeltaWithFEM(int index) override;
        virtual void ComputeBallToBallContactForce(SphericParticle::ParticleDataBuffer &,
                                                   ProcessInfo& r_process_info,
                                                   array_1d<double, 3>& rElasticForce,
                                                   array_1d<double, 3>& rContactForce,
                                                   double& RollingResistance) override final;

        virtual void ComputeBrokenBondsRatio();
        virtual void AddContributionToRepresentativeVolume(const double distance,
                                                    const double radius_sum,
                                                    const double contact_area);

        double*                     mSkinSphere;
        std::vector<int>            mFemIniNeighbourIds;
        std::vector<double>         mFemIniNeighbourDelta;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
            //rSerializer.save("mContinuumGroup",mContinuumGroup);
            //rSerializer.save("mIniNeighbourIds",mIniNeighbourIds);
            //rSerializer.save("mSymmStressTensor",mSymmStressTensor);
            rSerializer.save("mContinuumInitialNeighborsSize",mContinuumInitialNeighborsSize);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
            //rSerializer.load("mContinuumGroup",mContinuumGroup);
            //rSerializer.load("mIniNeighbourIds",mIniNeighbourIds);
            //rSerializer.load("mSymmStressTensor",mSymmStressTensor);
            rSerializer.load("mContinuumInitialNeighborsSize",mContinuumInitialNeighborsSize);
            mContinuumGroup = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);
            mSkinSphere     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
        }

        /* Assignment operator
        SphericContinuumParticle& operator=(SphericContinuumParticle const& rOther) { return *this; }
        Copy constructor
        SphericContinuumParticle(SphericContinuumParticle const& rOther) { *this = rOther; }
        */

    }; // Class SphericContinuumParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, SphericContinuumParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const SphericContinuumParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }
} // namespace Kratos

#endif // KRATOS_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED defined
