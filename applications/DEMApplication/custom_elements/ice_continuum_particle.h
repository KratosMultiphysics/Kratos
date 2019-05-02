//
// Author: Salva Latorre, latorre@cimne.upc.edu
//

#if !defined(KRATOS_ICECONTINUUMPARTICLE_H_INCLUDED)
#define KRATOS_ICECONTINUUMPARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// Project includes
#include "includes/define.h"
#include "spheric_continuum_particle.h"

namespace Kratos {
    
    class KRATOS_API(DEM_APPLICATION) IceContinuumParticle : public SphericContinuumParticle {
        
        public:

        /// Pointer definition of IceContinuumParticle
        KRATOS_CLASS_POINTER_DEFINITION(IceContinuumParticle);

        IceContinuumParticle() : SphericContinuumParticle() {}
        IceContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry) {}
        IceContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : SphericContinuumParticle(NewId, ThisNodes) {}
        IceContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SphericContinuumParticle(NewId, pGeometry, pProperties) {}

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return SphericContinuumParticle::Pointer(new IceContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Destructor
        virtual ~IceContinuumParticle() {};

        /// Turn back information as a string
        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "IceContinuumParticle" ;
            return buffer.str();
        }

        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "IceContinuumParticle";}

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const override {}    

        virtual array_1d<double,3> ComputeWeight(const array_1d<double,3>& gravity, const ProcessInfo& r_process_info) override;
        
        private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericContinuumParticle);
        }

        /*
        /// Assignment operator
        IceContinuumParticle& operator=(IceContinuumParticlev const& rOther)
        {
        return *this;
        }

        /// Copy constructor
        IceContinuumParticle(IceContinuumParticle const& rOther)
        {
        *this = rOther;
        }
        */

        ///@}

    }; // Class IceContinuumParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, IceContinuumParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const IceContinuumParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_ICECONTINUUMPARTICLE_H_INCLUDED defined
