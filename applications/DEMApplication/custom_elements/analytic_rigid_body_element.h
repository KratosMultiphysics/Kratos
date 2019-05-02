//   $Author: Guillermo Casas $

#if !defined(KRATOS_ANALYTIC_RIGID_BODY_ELEMENT_H_INCLUDED)
#define KRATOS_ANALYTIC_RIGID_BODY_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "rigid_body_element.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) AnalyticRigidBodyElement : public RigidBodyElement3D {

    public:
        /// Pointer definition of AnalyticRigidBodyElement
        KRATOS_CLASS_POINTER_DEFINITION(AnalyticRigidBodyElement);

        AnalyticRigidBodyElement();
        AnalyticRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry);
        AnalyticRigidBodyElement(IndexType NewId, NodesArrayType const& ThisNodes);
        AnalyticRigidBodyElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

        //Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~AnalyticRigidBodyElement();

        using Element::Initialize;

        virtual std::string Info() const override
        {
	    std::stringstream buffer;
        buffer << "AnalyticRigidBodyElement #" << Id();
	    return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
        rOStream << "AnalyticRigidBodyElement #" << Id();
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override
        {
	    //mpGeometry->PrintData(rOStream);
        }

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RigidBodyElement3D);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RigidBodyElement3D);
        }

    }; // Class AnalyticRigidBodyElement

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, AnalyticRigidBodyElement& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const AnalyticRigidBodyElement& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

} // namespace Kratos

#endif // KRATOS_ANALYTIC_RIGID_BODY_ELEMENT_H_INCLUDED  defined
