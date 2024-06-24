// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainElement : public UPwBaseElement
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwSmallStrainElement);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;
    /// The definition of the sizetype
    using SizeType = std::size_t;
    using UPwBaseElement::CalculateDerivativesOnInitialConfiguration;
    using UPwBaseElement::mConstitutiveLawVector;
    using UPwBaseElement::mIsInitialised;
    using UPwBaseElement::mRetentionLawVector;
    using UPwBaseElement::mStateVariablesFinalized;
    using UPwBaseElement::mStressVector;
    using UPwBaseElement::mThisIntegrationMethod;

    explicit UPwSmallStrainElement(IndexType NewId = 0) : UPwBaseElement(NewId) {}

    /// Constructor using an array of nodes
    UPwSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    UPwSmallStrainElement(IndexType                          NewId,
                          GeometryType::Pointer              pGeometry,
                          PropertiesType::Pointer            pProperties,
                          std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    ~UPwSmallStrainElement() override                              = default;
    UPwSmallStrainElement(const UPwSmallStrainElement&)            = delete;
    UPwSmallStrainElement& operator=(const UPwSmallStrainElement&) = delete;
    UPwSmallStrainElement(UPwSmallStrainElement&&)                 = delete;
    UPwSmallStrainElement& operator=(UPwSmallStrainElement&&)      = delete;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    using UPwBaseElement::SetValuesOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    using UPwBaseElement::CalculateOnIntegrationPoints;

    std::string Info() const override
    {
        return "U-Pw small strain Element #" + std::to_string(this->Id()) +
               "\nConstitutive law: " + mConstitutiveLawVector[0]->Info();
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:
    struct ElementVariables {
        /// Properties variables
        bool   IgnoreUndrained;
        bool   UseHenckyStrain;
        bool   ConsiderGeometricStiffness;
        double DynamicViscosityInverse;

        double                            BiotCoefficient;
        double                            BiotModulusInverse;
        BoundedMatrix<double, TDim, TDim> PermeabilityMatrix;

        /// ProcessInfo variables
        double VelocityCoefficient;
        double DtPressureCoefficient;

        /// Nodal variables
        array_1d<double, TNumNodes>        PressureVector;
        array_1d<double, TNumNodes>        DtPressureVector;
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        array_1d<double, TNumNodes * TDim> VelocityVector;
        array_1d<double, TNumNodes * TDim> VolumeAcceleration;

        /// Variables computed at each GP
        Matrix                                        B;
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu;
        array_1d<double, TDim>                        BodyAcceleration;
        array_1d<double, TDim>                        SoilGamma;

        /// Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix GradNpTInitialConfiguration;

        Matrix                                    F;
        Vector                                    detJContainer;
        Matrix                                    NContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

        /// Retention Law parameters
        double DegreeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;

        // needed for updated Lagrangian:
        double detJ;
        double detJInitialConfiguration;
        double IntegrationCoefficient;
        double IntegrationCoefficientInitialConfiguration;

        // Auxiliary Variables
        Matrix UVoigtMatrix;
    };

    void SaveGPStress(Matrix& rStressContainer, const Vector& rStressVector, unsigned int GPoint);

    void ExtrapolateGPValues(const Matrix& rStressContainer);

    void CalculateMaterialStiffnessMatrix(MatrixType&        rStiffnessMatrix,
                                          const ProcessInfo& CurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    virtual void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& CurrentProcessInfo);

    virtual void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    Matrix CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const;
    std::vector<Matrix> CalculateBMatrices(const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer,
                                           const Matrix& rNContainer) const;

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityMatrix(MatrixType&             rLeftHandSideMatrix,
                                                      const ElementVariables& rVariables);

    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint);

    void CalculateAndAddStiffnessForce(VectorType&             rRightHandSideVector,
                                       const ElementVariables& rVariables,
                                       unsigned int            GPoint);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, const ElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityFlow(VectorType&             rRightHandSideVector,
                                                    const ElementVariables& rVariables);

    virtual array_1d<double, TNumNodes> CalculateCompressibilityFlow(const ElementVariables& rVariables) const;

    [[nodiscard]] std::vector<double> CalculateRelativePermeabilityValues(const std::vector<double>& rFluidPressures) const;
    [[nodiscard]] std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const;
    virtual void CalculateAndAddPermeabilityFlow(VectorType&             rRightHandSideVector,
                                                 const ElementVariables& rVariables);
    virtual array_1d<double, TNumNodes> CalculatePermeabilityFlow(const ElementVariables& rVariables) const;

    virtual void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables);
    virtual array_1d<double, TNumNodes> CalculateFluidBodyFlow(const ElementVariables& rVariables) const;

    virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const;

    Matrix              CalculateDeformationGradient(unsigned int GPoint) const;
    std::vector<Matrix> CalculateDeformationGradients() const;

    void InitializeNodalDisplacementVariables(ElementVariables& rVariables);
    void InitializeNodalPorePressureVariables(ElementVariables& rVariables);
    void InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);
    std::vector<array_1d<double, TDim>> CalculateFluidFluxes(const std::vector<double>& rPermeabilityUpdateFactors,
                                                             const ProcessInfo& rCurrentProcessInfo);

    [[nodiscard]] std::vector<double> CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures) const;
    [[nodiscard]] std::vector<double> CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures) const;
    [[nodiscard]] virtual std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors) const;

    ///
    /// \brief This function calculates the constitutive matrices, stresses and strains depending on the
    ///        constitutive parameters. Note that depending on the settings in the rConstitutiveParameters
    ///        the function could calculate the stress, the constitutive matrix, the strains, or a combination.
    ///        In our elements we generally always calculate the constitutive matrix and sometimes the stress.
    ///
    void CalculateAnyOfMaterialResponse(const std::vector<Matrix>&   rDeformationGradients,
                                        ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                        const Matrix&                rNuContainer,
                                        const GeometryType::ShapeFunctionsGradientsType& rDNu_DXContainer,
                                        std::vector<Vector>& rStrainVectors,
                                        std::vector<Vector>& rStressVectors,
                                        std::vector<Matrix>& rConstitutiveMatrices);

    void CalculateExtrapolationMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rExtrapolationMatrix);

    void ResetHydraulicDischarge();
    void CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo);
    void CalculateSoilGamma(ElementVariables& rVariables);

    void CalculateAndAddGeometricStiffnessMatrix(MatrixType&   rLeftHandSideMatrix,
                                                 const Vector& rStressVector,
                                                 const Matrix& rGradNpt,
                                                 double        IntegrationCoefficient) const;

    VectorType GetPressureSolutionVector();

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

    template <class TValueType>
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }
};

// Class UPwSmallStrainElement

} // namespace Kratos