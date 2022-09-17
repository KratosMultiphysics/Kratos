//
// Author: Joaquín Irazábal jirazabal@cimne.upc.edu
//

#if !defined(KRATOS_BEAM_PARTICLE_H_INCLUDED)
#define KRATOS_BEAM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "spheric_continuum_particle.h"
#include "custom_constitutive/DEM_continuum_constitutive_law.h"
#include "custom_constitutive/DEM_beam_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) BeamParticle : public SphericContinuumParticle {

        public:

        /// Pointer definition of BeamParticle
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamParticle);

        typedef SphericContinuumParticle BaseType;
        typedef BaseType::ParticleDataBuffer BaseBufferType;
        typedef std::unique_ptr<BaseType::ParticleDataBuffer> BaseBufferPointerType;

        /// Default constructor.
        BeamParticle();
        BeamParticle( IndexType NewId, GeometryType::Pointer pGeometry);
        BeamParticle( IndexType NewId, NodesArrayType const& ThisNodes);
        BeamParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        BeamParticle(Element::Pointer p_continuum_spheric_particle);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~BeamParticle(){}

        /// Turn back information as a string
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "BeamParticle" ;
            return buffer.str();
        }

        /// Print information about this object
        void PrintInfo(std::ostream& rOStream) const override {rOStream << "BeamParticle";}

        /// Print object's data
        void PrintData(std::ostream& rOStream) const override {}

        void CreateContinuumConstitutiveLaws() override;
        void ContactAreaWeighting() override;
        void CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info) override;
        double CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info) override;

        virtual void Initialize(const ProcessInfo& r_process_info) override;

        virtual void InitializeSolutionStep(const ProcessInfo& r_process_info) override;

        virtual void ComputeBallToBallContactForceAndMoment(SphericParticle::ParticleDataBuffer &,
                                                   const ProcessInfo& r_process_info,
                                                   array_1d<double, 3>& rElasticForce,
                                                   array_1d<double, 3>& rContactForce) override;

        void Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) override;

        void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;

        using SphericContinuumParticle::CalculateOnContinuumContactElements;

        std::vector<Kratos::DEMBeamConstitutiveLaw::Pointer> mBeamConstitutiveLawArray;

        private:

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

}; // Class BeamParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, BeamParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const BeamParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

}  // namespace Kratos

#endif // KRATOS_BEAM_PARTICLE_H_INCLUDED defined
