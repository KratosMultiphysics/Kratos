// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#if !defined(KRATOS_SMALL_DISPLACEMENT_ANALYTIC_SENSITIVITIES_ELEMENT_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_ANALYTIC_SENSITIVITIES_ELEMENT_H_INCLUDED

// Project includes
#include "structural_mechanics_application.h"
#include "custom_elements/small_displacement.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{

class SmallDisplacementAnalyticSensitivityElement :
        public SmallDisplacement
{
public:

    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of SmallDisplacementAnalyticSensitivityElement
    KRATOS_CLASS_POINTER_DEFINITION( SmallDisplacementAnalyticSensitivityElement );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    SmallDisplacementAnalyticSensitivityElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementAnalyticSensitivityElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    SmallDisplacementAnalyticSensitivityElement(SmallDisplacementAnalyticSensitivityElement const& rOther);

    /// Destructor.
    virtual ~SmallDisplacementAnalyticSensitivityElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SmallDisplacementAnalyticSensitivityElement& operator=(SmallDisplacementAnalyticSensitivityElement const& rOther);

    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

    /// Overwritten function to calculate sensitivities
    using SmallDisplacement::Calculate;    
    void Calculate(const Variable<Vector> &rVariable, Vector &rOutput, const ProcessInfo &rCurrentProcessInfo);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    SmallDisplacementAnalyticSensitivityElement() : SmallDisplacement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    void CalculateDerivedDeformationMatrix(Matrix& rDB, const Matrix& rDN_DX, const int node_index, const int direction);


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization
    virtual void save(Serializer& rSerializer) const;
    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class SmallDisplacementAnalyticSensitivityElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_ANALYTIC_SENSITIVITIES_ELEMENT_H_INCLUDED  defined
