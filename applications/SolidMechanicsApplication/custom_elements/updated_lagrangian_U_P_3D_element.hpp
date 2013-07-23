//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_P_3D_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_P_3D_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/spatial_lagrangian_U_P_3D_element.hpp"


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

/// Spatial Lagrangian element for 3D geometries.

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D
 */

class UpdatedLagrangianUP3DElement
    : public SpatialLagrangianUP3DElement
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
 
    /// Counted pointer of UpdatedLagrangianUP3DElement
    KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangianUP3DElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    UpdatedLagrangianUP3DElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianUP3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianUP3DElement(UpdatedLagrangianUP3DElement const& rOther);

    /// Destructor.
    virtual ~UpdatedLagrangianUP3DElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianUP3DElement& operator=(UpdatedLagrangianUP3DElement const& rOther);

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

    //************* STARTING - ENDING  METHODS

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);


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
    UpdatedLagrangianUP3DElement() : SpatialLagrangianUP3DElement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS 
     */

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
				    GeneralVariables& rVariables, 
				    double& rIntegrationWeight);
  
    /**
     * Calculation and addition of the vectors of the RHS 
     */

    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector, 
				    GeneralVariables& rVariables, 
				    Vector& rVolumeForce, 
				    double& rIntegrationWeight);

    /**
     * Initialize Element General Variables
     */ 
    virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    void SetGeneralVariables  (GeneralVariables& rVariables,
			       ConstitutiveLaw::Parameters& rValues,
			       const int & rPointNumber);

    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(GeneralVariables& rVariables,
			     const double& rPointNumber);

    
    /**
     * Calculation of the Deformation Matrix  BL
     */
    virtual void CalculateDeformationMatrix(Matrix& rB,
					    Matrix& rF,
					    Matrix& rDN_DX);


    /**
     * Calculation of the DN_DX in the updated configuration
     */
    void CalculatePushForwardDN_DX(GeneralVariables& rVariables);

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

}; // Class UpdatedLagrangianUP3DElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_U_P_3D_ELEMENT_H_INCLUDED  defined 
