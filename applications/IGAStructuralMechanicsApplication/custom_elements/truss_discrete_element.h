#if !defined(KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED )
#define  KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

// External includes

// Project includes
#include "custom_elements/curve_base_discrete_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Truss element.
*/
class TrussDiscreteElement
    : public CurveBaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of TrussDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(TrussDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    TrussDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : CurveBaseDiscreteElement(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    TrussDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : CurveBaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    TrussDiscreteElement() : CurveBaseDiscreteElement() {};

    /// Destructor.
    virtual ~TrussDiscreteElement() override
    {};

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< TrussDiscreteElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * This functions calculates both the RHS and the LHS
    * @param rLeftHandSideMatrix: The LHS
    * @param rRightHandSideVector: The RHS
    * @param rCurrentProcessInfo: The current process info instance
    * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
    * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
    */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    ) override;

    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

    ///@}

};  // Class TrussDiscreteElement
///@}

}  // namespace Kratos.

#endif // KRATOS_TRUSS_DISCRETE_ELEMENT_H_INCLUDED  defined 