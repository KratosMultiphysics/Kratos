//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SPATIAL_LAGRANGIAN_U_wP_STAB_LAG_ELEMENT_H_INCLUDED )
#define  KRATOS_SPATIAL_LAGRANGIAN_U_wP_STAB_LAG_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/spatial_lagrangian_U_wP_Stab_element.hpp"

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

/// Spatial Lagrangian Large Displacement Lagrangian U-Pw Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)
// ONLY THE Lagrangian Permeability tensor terms


class SpatialLagrangianUwPStabLagElement
    : public SpatialLagrangianUwPStabElement
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

    /// Counted pointer of LargeDisplacementUPElement
    KRATOS_CLASS_POINTER_DEFINITION( SpatialLagrangianUwPStabLagElement );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    SpatialLagrangianUwPStabLagElement();

    /// Default constructors
    SpatialLagrangianUwPStabLagElement(IndexType NewId, GeometryType::Pointer pGeometry);

    SpatialLagrangianUwPStabLagElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    SpatialLagrangianUwPStabLagElement(SpatialLagrangianUwPStabLagElement const& rOther);


    /// Destructor.
    virtual ~SpatialLagrangianUwPStabLagElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    SpatialLagrangianUwPStabLagElement& operator=(SpatialLagrangianUwPStabLagElement const& rOther);


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


    virtual void GetPermeabilityTensor( const double& rPermeability, const Matrix& rF, Matrix& rPermeabilityTensor);

    virtual double GetPermeabilityLDTerm( const Matrix& rPermeability, const Matrix& rF, const int i, const int j, const int k, const int l);

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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class SpatialLagrangianUwPStabLagElement



} // namespace Kratos
#endif // KRATOS_____

