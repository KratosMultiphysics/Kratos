//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//


#if !defined(KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED

/* Project includes */
#include "includes/serializer.h"
#include "custom_elements/small_displacement_element.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

class SmallDisplacementThermoMechanicElement : public SmallDisplacementElement
{



public:

    ///Type for element variables
    typedef SmallDisplacementElement::ElementDataType ElementDataType;



    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallDisplacementThermoMechanicElement );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    SmallDisplacementThermoMechanicElement();

    // Constructor 1
    SmallDisplacementThermoMechanicElement(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor 2
    SmallDisplacementThermoMechanicElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Destructor
    virtual ~SmallDisplacementThermoMechanicElement();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallDisplacementElement )
    }


}; // Class SmallDisplacementThermoMechanicElement

} // namespace Kratos

#endif // KRATOS_SMALL_DISPLACEMENT_THERMO_MECHANIC_ELEMENT_H_INCLUDED  defined
