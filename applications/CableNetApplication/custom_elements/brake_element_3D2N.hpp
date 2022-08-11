// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//

#if !defined(KRATOS_Brake_ELEMENT_3D2N_H_INCLUDED )
#define  KRATOS_Brake_ELEMENT_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
/**
 * @class BrakeElement3D2N
 *
 * @brief This is a 3D-2node brake element with 3 translational dofs per node inheriting from the TrussElement3D2N
 *
 * @author Klaus B Sautter
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) BrakeElement3D2N : public TrussElement3D2N
{

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BrakeElement3D2N);



    BrakeElement3D2N() {};
    BrakeElement3D2N(IndexType NewId,
                     GeometryType::Pointer pGeometry);
    BrakeElement3D2N(IndexType NewId,
                     GeometryType::Pointer pGeometry,
                     PropertiesType::Pointer pProperties);


    ~BrakeElement3D2N() override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param ThisNodes The array containing nodes
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    virtual double ReturnReferenceLength() const override;
    virtual double ReturnCurrentLength() const override;

    
private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

}

#endif
