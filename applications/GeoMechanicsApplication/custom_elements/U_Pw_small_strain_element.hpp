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
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainElement : public UPwBaseElement<TDim, TNumNodes>
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
    using UPwBaseElement<TDim, TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mStressVector;
    using UPwBaseElement<TDim, TNumNodes>::mStateVariablesFinalized;
    using UPwBaseElement<TDim, TNumNodes>::mIsInitialised;
    using UPwBaseElement<TDim, TNumNodes>::CalculateDerivativesOnInitialConfiguration;
    using UPwBaseElement<TDim, TNumNodes>::mThisIntegrationMethod;

    explicit UPwSmallStrainElement(IndexType NewId = 0) : UPwBaseElement<TDim, TNumNodes>(NewId) {}

    /// Constructor using an array of nodes
    UPwSmallStrainElement(IndexType NewId, const NodesArrayType& ThisNodes, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    UPwSmallStrainElement(IndexType NewId, GeometryType::Pointer pGeometry, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    UPwSmallStrainElement(IndexType                          NewId,
                          GeometryType::Pointer              pGeometry,
                          PropertiesType::Pointer            pProperties,
                          std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
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

    using UPwBaseElement<TDim, TNumNodes>::SetValuesOnIntegrationPoints;

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

    std::string Info() const override
    {
        return "U-Pw small strain Element #" + std::to_string(this->Id()) +
               "\nConstitutive law: " + mConstitutiveLawVector[0]->Info();
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:
    static constexpr SizeType VoigtSize = (TDim == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN);
    static constexpr SizeType StressTensorSize =
        (TDim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D);

    struct ElementVariables {
        /// Properties variables
        bool   IgnoreUndrained;
        bool   UseHenckyStrain;
        bool   ConsiderGeometricStiffness;
        double DynamicViscosityInverse;
        double FluidDensity;
        double SolidDensity;
        double Density;
        double Porosity;
        double PermeabilityUpdateFactor;

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

        /// General elemental variables
        Vector VoigtVector;

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
        double                                    detF;
        Vector                                    detJContainer;
        Matrix                                    NContainer;
        GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

        /// Retention Law parameters
        double FluidPressure;
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;
        double EffectiveSaturation;

        // needed for updated Lagrangian:
        double detJ;
        double detJInitialConfiguration;
        double IntegrationCoefficient;
        double IntegrationCoefficientInitialConfiguration;

        // Auxiliary Variables
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes * TDim> UUMatrix;
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>        UPMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes * TDim>        PUMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes>               PPMatrix;
        Matrix                                                    UVoigtMatrix;
        BoundedMatrix<double, TNumNodes, TDim>                    PDimMatrix;
        array_1d<double, TNumNodes * TDim>                        UVector;
        array_1d<double, TNumNodes>                               PVector;
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

    void SetConstitutiveParameters(ElementVariables& rVariables, ConstitutiveLaw::Parameters& rConstitutiveParameters);

    void SetRetentionParameters(const ElementVariables& rVariables, RetentionLaw::Parameters& rRetentionParameters);

    virtual void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

    void InitializeBiotCoefficients(ElementVariables& rVariables, bool hasBiotCoefficient = false);

    void CalculatePermeabilityUpdateFactor(ElementVariables& rVariables);

    virtual void CalculateBMatrix(Matrix& rB, const Matrix& GradNpT, const Vector& Np);

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    virtual void CalculateCompressibilityMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                                const ElementVariables& rVariables) const;

    virtual void CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

    virtual void CalculatePermeabilityMatrix(BoundedMatrix<double, TNumNodes, TDim>& rPDimMatrix,
                                             BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                             const ElementVariables& rVariables) const;

    virtual void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint);

    void CalculateAndAddStiffnessForce(VectorType&       rRightHandSideVector,
                                       ElementVariables& rVariables,
                                       unsigned int      GPoint);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    virtual void CalculateCompressibilityFlow(BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                              array_1d<double, TNumNodes>&                 rPVector,
                                              const ElementVariables& rVariables) const;

    virtual void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    virtual void CalculatePermeabilityFlow(BoundedMatrix<double, TNumNodes, TDim>&      rPDimMatrix,
                                           BoundedMatrix<double, TNumNodes, TNumNodes>& rPMatrix,
                                           array_1d<double, TNumNodes>&                 rPVector,
                                           const ElementVariables& rVariables) const;

    virtual void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementVariables& rVariables);
    virtual void CalculateFluidBodyFlow(BoundedMatrix<double, TNumNodes, TDim>& rPDimMatrix,
                                        array_1d<double, TNumNodes>&            rPVector,
                                        const ElementVariables&                 rVariables) const;

    double CalculateBulkModulus(const Matrix& ConstitutiveMatrix) const;
    double CalculateBiotCoefficient(const ElementVariables& rVariables, bool hasBiotCoefficient) const;

    virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient);
    virtual void   CalculateCauchyStrain(ElementVariables& rVariables);
    virtual void   CalculateStrain(ElementVariables& rVariables, unsigned int GPoint);
    virtual void   CalculateDeformationGradient(ElementVariables& rVariables, unsigned int GPoint);

    void InitializeNodalDisplacementVariables(ElementVariables& rVariables);
    void InitializeNodalPorePressureVariables(ElementVariables& rVariables);
    void InitializeNodalVolumeAccelerationVariables(ElementVariables& rVariables);

    void   InitializeProperties(ElementVariables& rVariables);
    double CalculateFluidPressure(const ElementVariables& rVariables);

    void CalculateRetentionResponse(ElementVariables&         rVariables,
                                    RetentionLaw::Parameters& rRetentionParameters,
                                    unsigned int              GPoint);

    void CalculateExtrapolationMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rExtrapolationMatrix);

    void         ResetHydraulicDischarge();
    void         CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo);
    void         CalculateSoilGamma(ElementVariables& rVariables);
    virtual void CalculateSoilDensity(ElementVariables& rVariables);

    virtual void CalculateAndAddGeometricStiffnessMatrix(MatrixType&       rLeftHandSideMatrix,
                                                         ElementVariables& rVariables,
                                                         unsigned int      GPoint);

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

}; // Class UPwSmallStrainElement

} // namespace Kratos