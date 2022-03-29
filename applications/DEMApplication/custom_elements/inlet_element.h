#if !defined KRATOS_INLET_ELEMENT_H_INCLUDED
#define KRATOS_INLET_ELEMENT_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "custom_elements/rigid_body_element.h"

namespace Kratos
{
    class Element;
    class KRATOS_API(DEM_APPLICATION) InletElement3D : public RigidBodyElement3D {

    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(InletElement3D);

        InletElement3D();
        InletElement3D(IndexType NewId, GeometryType::Pointer pGeometry);
        InletElement3D(IndexType NewId, NodesArrayType const& ThisNodes);
        InletElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        ~InletElement3D();

        void Initialize(const ProcessInfo& r_process_info) override;
        void UpdateLinearDisplacementAndVelocityOfNodes() override;
        void UpdateAngularDisplacementAndVelocityOfNodes() override;

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

    };

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, InletElement3D& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const InletElement3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

}  // namespace Kratos

#endif
