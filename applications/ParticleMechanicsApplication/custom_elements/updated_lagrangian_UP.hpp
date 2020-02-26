//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_UP_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_UP_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian.hpp"

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

/// Large Displacement Lagrangian Element for 3D and 2D geometries. (base class)

/**
 * Implements a Large Displacement Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */

class UpdatedLagrangianUP
    : public UpdatedLagrangian
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

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( UpdatedLagrangianUP );
    ///@}




public:


    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    UpdatedLagrangianUP();


    /// Default constructors
    UpdatedLagrangianUP(IndexType NewId, GeometryType::Pointer pGeometry);

    UpdatedLagrangianUP(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    UpdatedLagrangianUP(UpdatedLagrangianUP const& rOther);

    /// Destructor.
    ~UpdatedLagrangianUP() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    UpdatedLagrangianUP& operator=(UpdatedLagrangianUP const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
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

    //************* STARTING - ENDING  METHODS

    /**
      * Called to initialize the element.
      * Must be called before any calculation is done
      */
    void Initialize() override;

    /**
     * Called at the beginning of each solution step
     */
    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    //************* COMPUTING  METHODS

    /**
      * this is called during the assembling process in order
      * to calculate the elemental mass matrix
      * @param rMassMatrix: the elemental mass matrix
      * @param rCurrentProcessInfo: the current process info instance
      */
    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             ProcessInfo& rCurrentProcessInfo) override;


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
        buffer << "MPM Element #" << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MPM Element #" << Id();
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
    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    //virtual void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
    //ProcessInfo& rCurrentProcessInfo);
    ///@}
    ///@name Protected Operations
    ///@{
    void FinalizeStepVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculation and addition of the matrices of the LHS
     */

    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                            GeneralVariables& rVariables,
                            const double& rIntegrationWeight,
                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */

    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
                            GeneralVariables& rVariables,
                            Vector& rVolumeForce,
                            const double& rIntegrationWeight) override;


    /**
     * Calculation of the Material Stiffness Matrix. Kuum = BT * C * B
     */

    void CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
                             GeneralVariables& rVariables,
                             const double& rIntegrationWeight) override;

    /**
     * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
     */
    void CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
                             GeneralVariables& rVariables,
                             const double& rIntegrationWeight) override;

    /**
     * Calculation of the Kup matrix
     */
    virtual void CalculateAndAddKup (MatrixType& rK,
                                     GeneralVariables & rVariables,
                                     const double& rIntegrationWeight
                                    );

    /**
     * Calculation of the Kpu matrix
     */
    virtual void CalculateAndAddKpu(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    const double& rIntegrationWeight
                                   );

    /**
     * Calculation of the Kpp matrix
     */
    virtual void CalculateAndAddKpp(MatrixType& rK,
                                    GeneralVariables & rVariables,
                                    const double& rIntegrationWeight
                                   );

    /**
     * Calculation of the Kpp Stabilization Term matrix
     */
    virtual void CalculateAndAddKppStab(MatrixType& rK,
                                        GeneralVariables & rVariables,
                                        const double& rIntegrationWeight
                                       );
    /**
     * Calculation of the External Forces Vector. Fe = N * t + N * b
     */
    void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
                                       GeneralVariables& rVariables,
                                       Vector& rVolumeForce,
                                       const double& rIntegrationWeight) override;

    /**
      * Calculation of the Internal Forces Vector. Fi = B * sigma
      */
    void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                       GeneralVariables & rVariables,
                                       const double& rIntegrationWeight) override;

    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            const double& rIntegrationWeight
                                              );

    /**
     * Calculation of the Internal Forces due to Pressure-Balance
     */
    virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            const double& rIntegrationWeight
                                                  );

    /**
     * Initialize System Matrices
     */
    void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
                                  VectorType& rRightHandSideVector,
                                  Flags& rCalculationFlags) override;

    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Initialize Element General Variables
     */
    void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;
    /**
     * Update the position of the MP or Gauss point when Finalize Element Internal Variables is called
     */

    void UpdateGaussPoint(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * Calculation of the constitutive coefficient for pressure of the Element
     */
    virtual double& CalculatePUCoefficient(double& rCoefficient, GeneralVariables & rVariables);

    /**
     * Calculation of the constitutive coefficient derivative for pressure  of the Element
     */
    virtual double& CalculatePUDeltaCoefficient(double& rCoefficient, GeneralVariables & rVariables);

    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    void GetHistoricalVariables( GeneralVariables& rVariables) override;

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
                                    Matrix& rF,
                                    Matrix& rDN_DX) override;

    /**
     * Calculation of the Volume Change of the Element
     */
    double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables) override;

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

}; // Class UpdatedLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED  defined
