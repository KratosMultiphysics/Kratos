#if !defined(KRATOS_MESHLESS_MEMBRANE_ELEMENT_H_INCLUDED )
#define  KRATOS_MESHLESS_MEMBRANE_ELEMENT_H_INCLUDED


// System includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "custom_elements/surface_base_discrete_element.h"

namespace Kratos
{
///@name Kratos Classes
///@{
/// Short class definition.
/** Discrete Membrane element.
*/
class MembraneDiscreteElement
    : public SurfaceBaseDiscreteElement
{
public:
    ///@name Type Definitions
    ///@{
    /// Counted pointer of MembraneDiscreteElement
    KRATOS_CLASS_POINTER_DEFINITION(MembraneDiscreteElement);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    // Constructor using an array of nodes
    MembraneDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : SurfaceBaseDiscreteElement(NewId, pGeometry)
    {};
    // Constructor using an array of nodes with properties
    MembraneDiscreteElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : SurfaceBaseDiscreteElement(NewId, pGeometry, pProperties)
    {};

    // default constructor necessary for serialization
    MembraneDiscreteElement() : SurfaceBaseDiscreteElement() {};

    /// Destructor.
    virtual ~MembraneDiscreteElement() override
    {};

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_shared< MembraneDiscreteElement >(NewId, GetGeometry().Create(ThisNodes), pProperties);
    };

    ///@}
    ///@name Operations
    ///@{

    /**
    * Called to initialize the element.
    * Must be called before any calculation is done
    */
    void Initialize() override;


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


    /**
    * Calculate a double Variable on the Element Constitutive Law
    * @param rVariable: The variable we want to get
    * @param rOutput: The values obtained int the integration points
    * @param rCurrentProcessInfo: the current process info instance
    */
    void Calculate(
        const Variable<double>& rVariable,
        double& rOutput,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * Calculate a Vector Variable on the Element
    * @param rVariable: The variable we want to get
    * @param rOutput: The values obtained int the integration points
    * @param rCurrentProcessInfo: the current process info instance
    */
    void Calculate(
        const Variable<Vector>& rVariable,
        Vector& rValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    /**
    * This is called during the assembling process in order to calculate the elemental mass matrix
    * @param rMassMatrix: the elemental mass matrix
    * @param rCurrentProcessInfo: the current process info instance
    */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo
    ) override;


    void CalculateStresses(
        Vector& rStresses,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculatePresstressTensor(
        Vector& rPrestressTensor,
        MetricVariables& rMetric);
    ///@}

protected:

private:
    ///@name Static Member Variables
    ///@{
    ///@name Operations
    ///@{
    void CalculateMetric( MetricVariables& rMetric ) override;


    /**
    * This functions updates the constitutive variables
    * @param rActualMetric: The actual metric
    * @param rThisConstitutiveVariables: The constitutive variables to be calculated
    * @param rValues: The CL parameters
    * @param ThisStressMeasure: The stress measure considered
    */
    void CalculateConstitutiveVariables(
        MetricVariables& rActualMetric,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) override;
    ///@}

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

};     // Class MembraneDiscreteElement
///@}

}  // namespace Kratos.

#endif // KRATOS_MESHLESS_MESHLESS_MEMBRANE_ELEMENT_H_INCLUDED  defined 