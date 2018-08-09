//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_J_W_WP_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_J_W_WP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_element.hpp"

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

/// Updated Lagrangian Large Displacement Lagrangian U-W Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)


class UpdatedLagrangianUJWwPElement
    : public UpdatedLagrangianUJElement
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
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUJWwPElement );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianUJWwPElement();

    /// Default constructors
    UpdatedLagrangianUJWwPElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianUJWwPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianUJWwPElement(UpdatedLagrangianUJWwPElement const& rOther);


    /// Destructor.
    virtual ~UpdatedLagrangianUJWwPElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianUJWwPElement& operator=(UpdatedLagrangianUJWwPElement const& rOther);


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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    //************* STARTING - ENDING  METHODS


    /**
    * Sets on rElementalDofList the degrees of freedom of the considered element geometry
    */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Sets on rValues the nodal displacements
     */
    void GetValuesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal velocities
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /**
     * Sets on rValues the nodal accelerations
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;


    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;


    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix, 
		    ProcessInfo& rCurrentProcessInfo) override;

    /**
      * this is called during the assembling process in order
      * to calculate the elemental damping matrix
      * @param rDampingMatrix: the elemental damping matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, 
		    ProcessInfo& rCurrentProcessInfo) override;

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

    
   /**** 
       the time step (requiered). It shall be somewhere else.
    ****/    
    double mTimeStep;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

 
    /**
     * Calculation and addition of the matrices of the LHS
     */

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ElementDataType& rVariables,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation of the water Material Stiffness Matrix. 
     */
    virtual void CalculateAndAddKWwP(MatrixType& rK,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight
                                    );

    /**
     * Calculation of the water pressure contrib to the internal forces 
     */
    virtual void CalculateAndAddKUwP(MatrixType& rK,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight
                                    );
    /**
     * Calculation of the stabilization at the matrix
     */
    virtual void CalculateAndAddKPPStab(MatrixType& rK,
                                     ElementDataType & rVariables,
                                     double& rIntegrationWeight
                                     );
    /**
     * Calculation and addition of the vectors of the RHS
     */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                                    ElementDataType& rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight) override;

    /**
     * Calculation of the Internal Forces due to sigma. Fi = B * (-pW * I)
     */
    virtual void CalculateAndAddInternalWaterForces(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Volumetric loads
     */
    void CalculateAndAddExternalForcesUJWwP(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          Vector & rVolumeForces,
          double& rIntegrationWeight
          );
    /**
     * Fluid Linear Momentum balance equation
     */
    void CalculateAndAddFluidLinearMomentum(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );
    /**
     * Mass balance for the mixture
     */
    void CalculateAndAddMassBalanceEquation(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Calculation of the Internal Forces due to stabilization
     */
    virtual void CalculateAndAddStabilizationRHS(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );
    /**
     * Part of the mass matrix due to the stabilization
     */
    virtual void CalculateAndAddMassStabilizationMatrix(MatrixType& rMassMatrix,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Part of the damping matrix due to the stabilization
     */
    virtual void CalculateAndAddDampingStabilizationMatrix(MatrixType& rDampingMatrix,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Part of the damping matrix due to the high order terms
     */
    virtual void CalculateAndAddHighOrderDampingMatrix(MatrixType& rDampingMatrix,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Calculation and addition of the KPP due to high order terms
     */
    virtual void CalculateAndAddHighOrderKPP(MatrixType& rK,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );

    /**
     * Calculation of the Internal Forces due to high order terms
     */
    virtual void CalculateAndAddHighOrderRHS(VectorType& rRightHandSideVector,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          );


    /**
     * Initialize Element General Variables
     */
    virtual void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  Flags& rCalculationFlags) override;


    /**
      * Calculate an estimate of the stabilization factor
      */
    double & CalculateStabilizationFactor( ElementDataType & rVariables, double & rStabilizationFactor);

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

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class UpdatedLagrangianUJWwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_W_wP_ELEMENT_H_INCLUDED

