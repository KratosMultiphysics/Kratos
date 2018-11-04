//
// Author: Miguel Angel Celigueta    maceli@cimne.upc.edu
//

#if !defined(KRATOS_NANOPARTICLE_H_INCLUDED )
#define  KRATOS_NANOPARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"

namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) NanoParticle : public SphericParticle
{
public:

    /// Pointer definition of NanoParticle
    KRATOS_CLASS_POINTER_DEFINITION(NanoParticle);

    using SphericParticle::GetGeometry;
    using SphericParticle::GetDensity;
    using SphericParticle::mRadius;

    NanoParticle():SphericParticle()
    {
        mThicknessOverRadius = 0.01; // Hard-coded but should go into node
    }

    NanoParticle( IndexType NewId, GeometryType::Pointer pGeometry ):SphericParticle(NewId, pGeometry){mThicknessOverRadius = 0.01;}
    NanoParticle( IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes){mThicknessOverRadius = 0.01;}
    NanoParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):SphericParticle(NewId, pGeometry, pProperties){mThicknessOverRadius = 0.01;}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return SphericParticle::Pointer(new NanoParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor.
    virtual ~NanoParticle();

    /// Assignment operator.
    NanoParticle& operator=(NanoParticle const& rOther); 

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "NanoParticle" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "NanoParticle";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}    
    void Initialize(const ProcessInfo& r_process_info) override;

    void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                 array_1d<double, 3>& additionally_applied_moment,
                                 const ProcessInfo& r_current_process_info,
                                 const array_1d<double,3>& gravity) override;

    void MemberDeclarationFirstStep(const ProcessInfo& r_process_info) override;
    
    double CalculateVolume() override;

    double GetInteractionRadius(const int radius_index = 0) override;
    void SetInteractionRadius(const double radius, const int radius_index = 0) override;
    void SetDefaultRadiiHierarchy(const double radius) override;

    double GetCationConcentration();

protected:

    double mThicknessOverRadius;
    double mInteractionRadius;

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

    /// Copy constructor.
    NanoParticle(NanoParticle const& rOther)
    {
    *this = rOther;
    }

    ///@}

}; // Class NanoParticle

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, NanoParticle& rThis) {return rIStream;}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const NanoParticle& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos

#endif // KRATOS_NANOPARTICLE_H_INCLUDED defined 


