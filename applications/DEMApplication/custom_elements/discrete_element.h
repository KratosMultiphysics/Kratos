//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_DISCRETE_ELEMENT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/indexed_object.h"
#include "containers/global_pointers_vector.h"
#include "includes/constitutive_law.h"

//Cfeng,RigidFace
#include "includes/condition.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DiscreteElement : public Element {
    public:
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(DiscreteElement);
        DiscreteElement(IndexType NewId = 0) : Element(NewId) {}
        DiscreteElement(IndexType NewId, const NodesArrayType& ThisNodes) : Element(NewId, ThisNodes) {}
        DiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {}
        DiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties) {}
        /// Copy constructor.
        DiscreteElement(DiscreteElement const& rOther) : Element(rOther) {
        }
        /// Destructor.
        virtual ~DiscreteElement() {
        }

        DiscreteElement& operator=(const DiscreteElement& rOther) {
            Element::operator=(rOther);
            return *this;
        }

        virtual void InitializeSolutionStep(const ProcessInfo& r_process_info) override {}

        virtual void FinalizeSolutionStep(const ProcessInfo& r_process_info) override {}

        virtual std::string Info() const override {
            std::stringstream buffer;
            buffer << "Discrete Element #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "Discrete Element #" << Id();}

        /// Print object's data.
        //virtual void PrintData(std::ostream& rOStream) const override { /*mpGeometry->PrintData(rOStream);*/ }

    protected:


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }


    }; // Class DiscreteElement


    /// input stream function
    inline std::istream& operator>>(std::istream& rIStream,
            DiscreteElement& rThis);

    /// output stream function

    inline std::ostream& operator<<(std::ostream& rOStream,
            const DiscreteElement& rThis) {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }


} // namespace Kratos.

#endif // KRATOS_DISCRETE_ELEMENT_H_INCLUDED  defined
