// Last Modified by: Salva, latorre@cimne.upc.edu

#if !defined KRATOS_SHIP_ELEMENT_3D_H_INCLUDED
#define KRATOS_SHIP_ELEMENT_3D_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"
#include "custom_elements/rigid_body_element.h"

namespace Kratos
{
    class Element;
    class ShipElement3D : public RigidBodyElement3D {
        
    public:
        /// Pointer definition of ShipElement3D
        KRATOS_CLASS_POINTER_DEFINITION(ShipElement3D);
       
        ShipElement3D();
        ShipElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
        ShipElement3D(IndexType NewId, NodesArrayType const& ThisNodes);
        ShipElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;      

        /// Destructor
        virtual ~ShipElement3D();

        void ComputeBuoyancyEffects();
        void ComputeEngineForce();
        void ComputeWaterDragForce();
        void ComputeAdditionalForces(const array_1d<double,3>& gravity);

        virtual std::string Info() const override
        {
	    std::stringstream buffer;
	    buffer << "Discrete Element #" << Id();
	    return buffer.str();
        }
      
        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
	    rOStream << "Discrete Element #" << Id();
        }
      
        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
	    //mpGeometry->PrintData(rOStream);
        }
 
    protected:
     
      
    private:
       
        friend class Serializer;
        virtual void save(Serializer& rSerializer) const { KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RigidBodyElement3D); }
        virtual void load(Serializer& rSerializer) { KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RigidBodyElement3D); }

    }; // Class ShipElement3D
   
    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, ShipElement3D& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const ShipElement3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }
 
}  // namespace Kratos

#endif // KRATOS_SHIP_ELEMENT_3D_H_INCLUDED defined
