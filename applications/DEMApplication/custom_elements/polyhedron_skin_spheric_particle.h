//
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu
//

#if !defined(KRATOS_POLYHEDRON_SKIN_SPHERIC_PARTICLE_H_INCLUDED )
#define  KRATOS_POLYHEDRON_SKIN_SPHERIC_PARTICLE_H_INCLUDED

// System includes


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"


namespace Kratos
{
    class KRATOS_API(DEM_APPLICATION) PolyhedronSkinSphericParticle : public SphericParticle
    {
    public:

        /// Pointer definition of PolyhedronSkinSphericParticle
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PolyhedronSkinSphericParticle);

        typedef GlobalPointersVector<Element> ParticleWeakVectorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
        typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

        /// Default constructor
        PolyhedronSkinSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry);
        PolyhedronSkinSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes);
        PolyhedronSkinSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~PolyhedronSkinSphericParticle();


        /// Turn back information as a string
        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "PolyhedronSkinSphericParticle" ;
            return buffer.str();
        }

        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "PolyhedronSkinSphericParticle";}

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const override {}


    protected:

        PolyhedronSkinSphericParticle();


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
        }

        /* Assignment operator
        PolyhedronSkinSphericParticle& operator=(PolyhedronSkinSphericParticle const& rOther) { return *this; }
        Copy constructor
        PolyhedronSkinSphericParticle(PolyhedronSkinSphericParticle const& rOther) { *this = rOther; }
        */

    }; // Class PolyhedronSkinSphericParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, PolyhedronSkinSphericParticle& rThis) {return rIStream;}

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const PolyhedronSkinSphericParticle& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }
} // namespace Kratos

#endif // KRATOS_POLYHEDRON_SKIN_SPHERIC_PARTICLE_H_INCLUDED defined
