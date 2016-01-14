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
class KRATOS_DEMAPPLICATION_EXPORT_DLL NanoParticle : public SphericParticle
{
public:

    /// Pointer definition of NanoParticle
    KRATOS_CLASS_POINTER_DEFINITION(NanoParticle);

    NanoParticle():SphericParticle(){}
    NanoParticle( IndexType NewId, GeometryType::Pointer pGeometry ):SphericParticle(NewId, pGeometry){}
    NanoParticle( IndexType NewId, NodesArrayType const& ThisNodes):SphericParticle(NewId, ThisNodes){}
    NanoParticle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ):SphericParticle(NewId, pGeometry, pProperties){}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return SphericParticle::Pointer(new NanoParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor.
    virtual ~NanoParticle(){};    


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "NanoParticle" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "NanoParticle";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    void ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
                                 array_1d<double, 3>& additionally_applied_moment,
                                 ProcessInfo& r_current_process_info,
                                 const array_1d<double,3>& gravity){
        KRATOS_TRY

        array_1d<double, 3> brownian_motion_force; brownian_motion_force.clear();
        array_1d<double, 3> van_der_waals_force; van_der_waals_force.clear();
        array_1d<double, 3> double_layer_force; double_layer_force.clear();

        this->ComputeBrownianMotionForce(brownian_motion_force, r_current_process_info);
        this->ComputeVanDerWaalsForce(van_der_waals_force, r_current_process_info);
        this->ComputeDoubleLayerForce(double_layer_force, r_current_process_info);

        additionally_applied_force += brownian_motion_force + van_der_waals_force + double_layer_force;


        //Now add the contribution of base class function (gravity or other forces added in upper levels):
        SphericParticle::ComputeAdditionalForces(additionally_applied_force, additionally_applied_moment, r_current_process_info, gravity);

        KRATOS_CATCH( "" )
    }

protected:

    void ComputeBrownianMotionForce(array_1d<double, 3>& brownian_motion_force, ProcessInfo& r_current_process_info){};
    void ComputeVanDerWaalsForce(array_1d<double, 3>& van_der_waals_force, ProcessInfo& r_current_process_info){};
    void ComputeDoubleLayerForce(array_1d<double, 3>& double_layer_force, ProcessInfo& r_current_process_info){};


private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
    }

    /*
    /// Assignment operator.
    NanoParticle& operator=(NanoParticle const& rOther)
    {
    return *this;
    }

    /// Copy constructor.
    NanoParticle(NanoParticle const& rOther)
    {
    *this = rOther;
    }
    */

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


