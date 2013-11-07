//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/spatial_lagrangian_U_P_element.hpp"


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

/// Updated Lagrangian U-P Element for 3D and 2D geometries. Linear Triangles and Tetrahedra

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D
 */

class UpdatedLagrangianUPElement
    : public SpatialLagrangianUPElement
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

    /// Counted pointer of UpdatedLagrangianUPElement
    KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangianUPElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    UpdatedLagrangianUPElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianUPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianUPElement(UpdatedLagrangianUPElement const& rOther);

    /// Destructor.
    virtual ~UpdatedLagrangianUPElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianUPElement& operator=(UpdatedLagrangianUPElement const& rOther);

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
    UpdatedLagrangianUPElement() : SpatialLagrangianUPElement()
    {
    }

    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */
    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                            GeneralVariables& rVariables,
                            double& rIntegrationWeight);

    /**
     * Calculation and addition of the vectors of the RHS
     */
    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                            GeneralVariables& rVariables,
                            Vector& rVolumeForce,
                            double& rIntegrationWeight);

    /**
     * Initialize Element General Variables
     */
    void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);


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
    void CalculateDeformationMatrix(Matrix& rB,
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

}; // Class UpdatedLagrangianUPElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_U_P_ELEMENT_H_INCLUDED  defined 
