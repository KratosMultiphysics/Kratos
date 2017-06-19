// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_AXISYM_KINEMATIC_LINEAR_H_INCLUDED )
#define  KRATOS_AXISYM_KINEMATIC_LINEAR_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/kinematic_linear.h"

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

/// Axisymmetric Kinemaic Linear element 

/**
 * Implements a Axisymmetric Kinemaic Linear definition for structural analysis.
 */

class AxisymKinematicLinear
    : public KinematicLinear
{
public:
    ///@name Type Definitions
    ///@{
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of AxisymKinematicLinear
    KRATOS_CLASS_POINTER_DEFINITION(AxisymKinematicLinear);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AxisymKinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry);
    AxisymKinematicLinear(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /// Destructor.
    virtual ~AxisymKinematicLinear();

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    //TODO: ADD THE OTHER CREATE FUNCTION
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    //std::string Info() const;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    //      virtual String Info() const;

    /// Print information about this object.
    //      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    //      virtual void PrintData(std::ostream& rOStream) const;
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
    
    AxisymKinematicLinear() : KinematicLinear()
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

    /**
     * Calculation of the Deformation Matrix B
     * @param B: The deformation matrix
     * @param DN_DX: The derivatives of the shape functions
     */
    void CalculateB(
        Matrix& rB,
        const Matrix& DN_DX,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber
        ) override;
    
    /**
     * Calculation of the equivalent deformation gradient
     * @param StrainVector: The strain tensor (Voigt notation)
     * @return The deformation gradient F
     */
    Matrix ComputeEquivalentF(const Vector& rStrainVector) override;
    
    /**
     * This functions computes the integration weight to consider
     * @param IntegrationPoints: The array containing the integration points
     * @param PointNumber: The id of the integration point considered
     * @param detJ: The determinant of the jacobian of the element
     */
    double GetIntegrationWeight(
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const unsigned int PointNumber,
        const double detJ
        ) override;

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    
    friend class Serializer;

    // A private default constructor necessary for serialization

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //AxisymKinematicLinear& operator=(const AxisymKinematicLinear& rOther);
    /// Copy constructor.
    //AxisymKinematicLinear(const AxisymKinematicLinear& rOther);
    ///@}

}; // Class AxisymKinematicLinear

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_AXISYM_KINEMATIC_LINEAR_H_INCLUDED  defined 
