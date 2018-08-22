//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               MCaicedo $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SMALL_DISPLACEMENT_BBAR_ELEMENT_H_INCLUDED)
#define KRATOS_SMALL_DISPLACEMENT_BBAR_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/small_displacement_element.hpp"

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

/// Small Displacement Bbar Element for 3D and 2D geometries.

/**
 * Implements a Small Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) SmallDisplacementBbarElement
  : public SmallDisplacementElement
{
public:

    ///@name Type Definitions
    ///@{
    /// Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    /// Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    /// StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    /// Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;
    ///Type for size
    typedef GeometryData::SizeType SizeType;
    ///Type for element variables
    typedef SmallDisplacementElement::ElementDataType ElementDataType;

    /// Counted pointer of SmallDisplacementBbarElement
    KRATOS_CLASS_POINTER_DEFINITION(SmallDisplacementBbarElement);
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    SmallDisplacementBbarElement();

    /// Default constructors
    SmallDisplacementBbarElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SmallDisplacementBbarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Copy constructor
    SmallDisplacementBbarElement(SmallDisplacementBbarElement const& rOther);

    /// Destructor.
    ~SmallDisplacementBbarElement() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SmallDisplacementBbarElement& operator=(SmallDisplacementBbarElement const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

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
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber) override;

    /**
     * Initialize Element General Variables
     */
    void InitializeElementData(ElementDataType& rVariables,
                                    const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculation of the Volumetric Deformation Matrix  H
     */
    void CalculateVolumetricDeformationMatrix(ElementDataType& rVariables);

    /**
     * Calculation of the Deformation Matrix B_bar B
     */
    void CalculateDeformationMatrixBbar(Matrix& rB, const Matrix& rH, const Matrix& rDN_DX);

    /**
     * Calculation of the Infinitesimal B_bar Strain
     */
    void CalculateInfinitesimalStrainBbar(Vector& rStrainVector, const Matrix& rH, const Matrix& rDN_DX);


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
    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class SmallDisplacementBbarElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_SMALL_DISPLACEMENT_BBAR_ELEMENT_H_INCLUDED  defined
