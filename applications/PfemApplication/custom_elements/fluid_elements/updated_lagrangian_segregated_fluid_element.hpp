//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED)
#define  KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/updated_lagrangian_segregated_fluid_element.hpp"


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

/// Large Displacement Lagrangian VP Element for 3D and 2D geometries

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for linear Triangles and Tetrahedra (base class)
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) UpdatedLagrangianSegregatedFluidElement
    : public LargeDisplacementVElement
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

    /// Counted pointer of UpdatedLagrangianSegregatedFluidElement
    KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianSegregatedFluidElement );
    ///@}

    enum StepType {VELOCITY_STEP = 1, PRESSURE_STEP = 2}
    
    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianSegregatedFluidElement();

    /// Default constructors
    UpdatedLagrangianSegregatedFluidElement(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianSegregatedFluidElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianSegregatedFluidElement(UpdatedLagrangianSegregatedFluidElement const& rOther);


    /// Destructor.
    virtual ~UpdatedLagrangianSegregatedFluidElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianSegregatedFluidElement& operator=(UpdatedLagrangianSegregatedFluidElement const& rOther);

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

    //************* GETTING METHODS

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

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    

    /**
     * Called at the end of eahc solution step
     */
    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;
    
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Updated Lagrangian Fluid Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Updated Lagrangian Fluid Element #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      GetGeometry().PrintData(rOStream);
    }
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

    /**
     * Finalize and Initialize label
     */
    bool mFinalizedStep;
   
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculation and addition of the matrices of the LHS
     */

    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                            ElementVariables& rVariables,
                            double& rIntegrationWeight) override;
    
    /**
     * Get element size from the dofs
     */    
    unsigned int GetDofsSize() override;

    /**
     * Calculation of the Material Stiffness Matrix. Kvvm = BT * C * B
     */

    virtual void CalculateAndAddKvvm(MatrixType& rLeftHandSideMatrix,
                                     ElementVariables& rVariables,
                                     double& rIntegrationWeight);

    /**
     * Calculation of the Geometric Stiffness Matrix. Kvvg = BT * S
     */
    virtual void CalculateAndAddKvvg(MatrixType& rLeftHandSideMatrix,
                                     ElementVariables& rVariables,
                                     double& rIntegrationWeight);    

    /**
     * Calculation of the Pressure Stiffness Matrix. Kpp
     */
    virtual void CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                                    ElementVariables& rVariables,
                                    double& rIntegrationWeight);
    /**
     * Set Variables of the Element to the Parameters of the Constitutive Law
     */
    void SetElementVariables(ElementVariables& rVariables,
                             ConstitutiveLaw::Parameters& rValues,
                             const int & rPointNumber) override;
    
    /**
     * Calculation of the velocity gradient
     */
    void CalculateVelocityGradient(Matrix& rH,
                                   const Matrix& rDN_DX,
                                   unsigned int step = 0);

    /**
     * Calculation of the velocity gradient
     */
    void CalculateVelocityGradientVector(Vector& rH,
                                         const Matrix& rDN_DX,
                                         unsigned int step = 0);

    
    /**
     * Calculation of the symmetric velocity gradient Vector
     */
    void CalculateSymmetricVelocityGradient(const Matrix& rH,
                                            Vector& rStrainVector);

    /**
     * Calculation of the skew symmetric velocity gradient Vector
     */
    void CalculateSkewSymmetricVelocityGradient(const Matrix& rH,
                                                Vector& rStrainVector);

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

}; // Class UpdatedLagrangianSegregatedFluidElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_SEGREGATED_FLUID_ELEMENT_H_INCLUDED  defined 
