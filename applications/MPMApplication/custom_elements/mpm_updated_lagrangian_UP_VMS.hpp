//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martínez
//

#if !defined(KRATOS_mpm_updated_lagrangian_UP_VMS_H_INCLUDED )
#define  KRATOS_mpm_updated_lagrangian_UP_VMS_H_INCLUDED


#include "includes/element.h"
#include "custom_elements/mpm_updated_lagrangian_UP.hpp"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * Updated Lagrangian U-P MPM element with ASGS stabilization.
 */
class MPMUpdatedLagrangianUPVMS
    : public MPMUpdatedLagrangianUP
{
public:

    ///@name Type Definitions
    ///@{
    ///base type: Element
    typedef Element BaseType;
    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;
    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// Counted pointer of LargeDisplacementElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MPMUpdatedLagrangianUPVMS );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    MPMUpdatedLagrangianUPVMS();


    /// Default constructors
    MPMUpdatedLagrangianUPVMS(IndexType NewId, GeometryType::Pointer pGeometry);

    MPMUpdatedLagrangianUPVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    MPMUpdatedLagrangianUPVMS(MPMUpdatedLagrangianUPVMS const& rOther);

    /// Destructor.
    ~MPMUpdatedLagrangianUPVMS() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MPMUpdatedLagrangianUPVMS& operator=(MPMUpdatedLagrangianUPVMS const& rOther);

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

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

protected:
    ///@name Protected Static Member Variables
    ///@{

    static const unsigned int msIndexVoigt3D6C [6][2];
    static const unsigned int msIndexVoigt2D4C [4][2];
    static const unsigned int msIndexVoigt2D3C [3][2];

    ///@}
    ///@name Protected Operations
    ///@{

    double CalculateElementSize() const;

    /**
     * Calculates the elemental contributions
     * \f$ K^e = w\,B^T\,D\,B \f$ and
     * \f$ r^e \f$
     */
    void CalculateElementalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag) override;

    // Calculation of stabilization parameters
    void CalculateTaus(GeneralVariables& rVariables);

    double CalculateShearModulus() const;

    double CalculateBulkModulus() const;

    // To compute identity tensor
    void CalculateTensorIdentityMatrix(Matrix& rTensorIdentityMatrix);

    double CalculateTensorIdentityComponent(
        const unsigned int a,
        const unsigned int b,
        const unsigned int c,
        const unsigned int d) const;

    // To compute vector in voigt notation to multiply
    void ConvertPressureGradientInVoigt(const Vector& rPressureGradient, Vector& rPressureGradientVoigt);

    // Calculation of ASGS stabilization variables
    void CalculateStabilizationVariables(
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculatePressureVariables(GeneralVariables& rVariables);

    void CalculateDynamicStabilizationVariables(
        GeneralVariables& rVariables,
        const ProcessInfo& rCurrentProcessInfo);
    /*
     * Operations related with the osgs stabilization
     */
    void ComputeResidual(GeneralVariables& rVariables, Vector& rVolumeForce, Vector& rResidualU, double& rResidualP);

    virtual void CalculateProjections(const ProcessInfo &rCurrentProcessInfo);

    /**
     * Calculation and addition of the matrices of the LHS
     */
    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                            GeneralVariables& rVariables,
                            const double& rIntegrationWeight,
                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculation and addition of the vectors of the RHS
     */
    void CalculateAndAddRHS(VectorType& rRightHandSideVector,
                            GeneralVariables& rVariables,
                            Vector& rVolumeForce,
                            const double& rIntegrationWeight,
                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculation of the Kuu Stabilization Term matrix
     */
    void CalculateAndAddKuuStab(MatrixType& rK,
                                GeneralVariables& rVariables,
                                const double& rIntegrationWeight);

    /**
     * Calculation of the Kup Stabilization Term matrix
     */
    void CalculateAndAddKupStab(MatrixType& rK,
                                GeneralVariables& rVariables,
                                const double& rIntegrationWeight);

    /**
     * Calculation of the Kpu Stabilization Term matrix
     */
    void CalculateAndAddKpuStab(MatrixType& rK,
                                GeneralVariables& rVariables,
                                const double& rIntegrationWeight);

    /**
     * Calculation of the Kpp Stabilization Term matrix
     */
    void CalculateAndAddKppStab(MatrixType& rK,
                                GeneralVariables& rVariables,
                                const double& rIntegrationWeight) override;

    /**
     * Calculation of stabilization terms for the continuity equation
     */
    void CalculateAndAddStabilizedPressureVMS(VectorType& rRightHandSideVector,
                                              GeneralVariables& rVariables,
                                              Vector& rVolumeForce,
                                              const double& rIntegrationWeight);

    /**
     * Calculation stabilization terms for the momentum equation
     */
    void CalculateAndAddStabilizedDisplacementVMS(VectorType& rRightHandSideVector,
                                                  GeneralVariables& rVariables,
                                                  Vector& rVolumeForce,
                                                  const double& rIntegrationWeight);

    /**
     * Initialize Element General Variables
     */
    void InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) override;

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    // A private default constructor necessary for serialization

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class MPMUpdatedLagrangianUPVMS

///@}

} // namespace Kratos.
#endif // KRATOS_mpm_updated_lagrangian_UP_VMS_H_INCLUDED defined
