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
    class KRATOS_API(DEM_APPLICATION) ShipElement3D : public RigidBodyElement3D {

    public:
        /// Pointer definition of ShipElement3D
        KRATOS_CLASS_POINTER_DEFINITION(ShipElement3D);

        ShipElement3D();
        ShipElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
        ShipElement3D(IndexType NewId, NodesArrayType const& ThisNodes);
        ShipElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~ShipElement3D();

        void CustomInitialize(ModelPart& rigid_body_element_sub_model_part) override;

        void ComputeBuoyancyEffects();
        void ComputeEngineForce();
        void ComputeWaterDragForce();

        void ComputeExternalForces(const array_1d<double,3>& gravity) override;

        // Engine characteristics
        double mEnginePower; //60000000; // 60MW for the Arktika-class icebreaker
        double mMaxEngineForce; //60000000; // with 20MN the ship almost couldn't make it through the ice
        double mThresholdVelocity; //1.0; // It was set to 3.0 m/s before, which corresponded to a maximum force of 20MN
        double mEnginePerformance;

        // Water drag. Applied to de center of mass. We are assuming the ship is moving in the X direction (in GiD)
        // Drag constant values in the 3 axes
        // drag_X = 500000; // Such that the X maximum velocity is 11 m/s, which corresponds to Arktika-class icebreakers
        // drag_Y = 240000000; // Such that the Y maximum velocity is 0.5 m/s
        // drag_Z = 240000000; // Such that the Z maximum velocity is 0.5 m/s
        array_1d<double,3> mDragConstantVector;

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
        virtual void save(Serializer& rSerializer) const override{ KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RigidBodyElement3D); }
        virtual void load(Serializer& rSerializer) override{ KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RigidBodyElement3D); }

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
