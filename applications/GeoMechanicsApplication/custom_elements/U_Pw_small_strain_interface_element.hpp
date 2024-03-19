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
#include "custom_utilities/interface_element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainInterfaceElement
    : public UPwBaseElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwSmallStrainInterfaceElement);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;
    using UPwBaseElement<TDim, TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim, TNumNodes>::mStressVector;
    using UPwBaseElement<TDim, TNumNodes>::mStateVariablesFinalized;
    using UPwBaseElement<TDim, TNumNodes>::CalculateDerivativesOnInitialConfiguration;
    using UPwBaseElement<TDim, TNumNodes>::mThisIntegrationMethod;

    /// Default Constructor
    UPwSmallStrainInterfaceElement(IndexType NewId = 0) : UPwBaseElement<TDim, TNumNodes>(NewId) {}

    /// Constructor using an array of nodes
    UPwSmallStrainInterfaceElement(IndexType                          NewId,
                                   const NodesArrayType&              ThisNodes,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    UPwSmallStrainInterfaceElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    UPwSmallStrainInterfaceElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   PropertiesType::Pointer            pProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
        /// Lobatto integration method with the integration points located at the "mid plane nodes" of the interface
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    }

    ~UPwSmallStrainInterfaceElement() override                                       = default;
    UPwSmallStrainInterfaceElement(const UPwSmallStrainInterfaceElement&)            = delete;
    UPwSmallStrainInterfaceElement& operator=(const UPwSmallStrainInterfaceElement&) = delete;
    UPwSmallStrainInterfaceElement(UPwSmallStrainInterfaceElement&&)                 = delete;
    UPwSmallStrainInterfaceElement& operator=(UPwSmallStrainInterfaceElement&&)      = delete;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

protected:
    struct SFGradAuxVariables {
        array_1d<double, TDim> GlobalCoordinatesGradients;
        array_1d<double, TDim> LocalCoordinatesGradients;

        BoundedMatrix<double, TNumNodes, TDim - 1> ShapeFunctionsNaturalGradientsMatrix;
        BoundedMatrix<double, TDim - 1, TDim - 1>  LocalCoordinatesGradientsMatrix;
        BoundedMatrix<double, TDim - 1, TDim - 1>  LocalCoordinatesGradientsInvMatrix;
        BoundedMatrix<double, TNumNodes, TDim - 1> ShapeFunctionsGradientsMatrix;
    };

    struct InterfaceElementVariables {
        /// Properties variables
        bool   IgnoreUndrained;
        double DynamicViscosityInverse;
        double FluidDensity;
        double SolidDensity;
        double Density;
        double Porosity;
        double BiotCoefficient;
        double BiotModulusInverse;

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
        BoundedMatrix<double, TDim, TDim> RotationMatrix;
        array_1d<double, TDim>            VoigtVector;

        /// Variables computed at each GP

        /// Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Np;
        Matrix GradNpT;
        Matrix F;
        double detF;

        /// Auxiliary Variables
        BoundedMatrix<double, TDim, TNumNodes * TDim> Nu;
        BoundedMatrix<double, TDim, TDim>             LocalPermeabilityMatrix;
        array_1d<double, TDim>                        BodyAcceleration;
        array_1d<double, TDim>                        SoilGamma;

        double                                                    IntegrationCoefficient;
        double                                                    JointWidth;
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes * TDim> UUMatrix;
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>        UPMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes * TDim>        PUMatrix;
        BoundedMatrix<double, TNumNodes, TNumNodes>               PPMatrix;
        BoundedMatrix<double, TDim, TDim>                         DimMatrix;
        BoundedMatrix<double, TNumNodes * TDim, TDim>             UDimMatrix;
        BoundedMatrix<double, TNumNodes, TDim>                    PDimMatrix;
        array_1d<double, TNumNodes * TDim>                        UVector;
        array_1d<double, TNumNodes>                               PVector;

        /// Retention Law parameters
        double FluidPressure;
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;
    };

    /// Member Variables
    std::vector<double> mInitialGap;
    std::vector<bool>   mIsOpen;

    void ModifyInactiveElementStress(const double& JointWidth, Vector& StressVector);

    virtual void CalculateOnLobattoIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                     std::vector<array_1d<double, 3>>&    rOutput,
                                                     const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateOnLobattoIntegrationPoints(const Variable<Matrix>& rVariable,
                                                     std::vector<Matrix>&    rOutput,
                                                     const ProcessInfo&      rCurrentProcessInfo);

    virtual void CalculateOnLobattoIntegrationPoints(const Variable<Vector>& rVariable,
                                                     std::vector<Vector>&    rOutput,
                                                     const ProcessInfo&      rCurrentProcessInfo);

    void CalculateInitialGap(const GeometryType& Geom);

    void ExtrapolateGPValues(const std::vector<double>& JointWidthContainer);

    void CalculateMaterialStiffnessMatrix(MatrixType&        rStiffnessMatrix,
                                          const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      const bool         CalculateStiffnessMatrixFlag,
                      const bool         CalculateResidualVectorFlag) override;

    virtual void InitializeElementVariables(InterfaceElementVariables& rVariables,
                                            const GeometryType&        Geom,
                                            const PropertiesType&      Prop,
                                            const ProcessInfo&         CurrentProcessInfo);

    void CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix, const GeometryType& Geom);

    void CalculateJointWidth(double&             rJointWidth,
                             const double&       NormalRelDisp,
                             const double&       MinimumJointWidth,
                             const unsigned int& GPoint);

    void CheckAndCalculateJointWidth(double&                      rJointWidth,
                                     ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                     double&                      rNormalRelDisp,
                                     double                       MinimumJointWidth,
                                     unsigned int                 GPoint);

    template <class TMatrixType>
    void CalculateShapeFunctionsGradients(TMatrixType&                             rGradNpT,
                                          SFGradAuxVariables&                      rAuxVariables,
                                          const Matrix&                            Jacobian,
                                          const BoundedMatrix<double, TDim, TDim>& RotationMatrix,
                                          const Matrix&                            DN_De,
                                          const Matrix&                            Ncontainer,
                                          const double&                            JointWidth,
                                          const unsigned int&                      GPoint);

    virtual void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix,
                                                      InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddPermeabilityMatrix(MatrixType&                rLeftHandSideMatrix,
                                                   InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddRHS(VectorType&                rRightHandSideVector,
                                    InterfaceElementVariables& rVariables,
                                    unsigned int               GPoint);

    void CalculateAndAddStiffnessForce(VectorType&                rRightHandSideVector,
                                       InterfaceElementVariables& rVariables,
                                       unsigned int               GPoint);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddCompressibilityFlow(VectorType&                rRightHandSideVector,
                                                    InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddPermeabilityFlow(VectorType&                rRightHandSideVector,
                                                 InterfaceElementVariables& rVariables);

    virtual void CalculateAndAddFluidBodyFlow(VectorType&                rRightHandSideVector,
                                              InterfaceElementVariables& rVariables);

    void InterpolateOutputDoubles(std::vector<double>& rOutput, const std::vector<double>& GPValues);

    template <class TValueType>
    void InterpolateOutputValues(std::vector<TValueType>& rOutput, const std::vector<TValueType>& GPValues);

    void SetRetentionParameters(const InterfaceElementVariables& rVariables,
                                RetentionLaw::Parameters&        rRetentionParameters);

    double CalculateFluidPressure(const InterfaceElementVariables& rVariables);

    double CalculateBulkModulus(const Matrix& ConstitutiveMatrix);

    void InitializeBiotCoefficients(InterfaceElementVariables& rVariables, const bool& hasBiotCoefficient = false);

    void CalculateRetentionResponse(InterfaceElementVariables& rVariables,
                                    RetentionLaw::Parameters&  rRetentionParameters,
                                    unsigned int               GPoint);

    void CalculateSoilGamma(InterfaceElementVariables& rVariables);

    void CalculateSoilDensity(InterfaceElementVariables& rVariables);

    void SetConstitutiveParameters(InterfaceElementVariables&   rVariables,
                                   ConstitutiveLaw::Parameters& rConstitutiveParameters);

    Vector SetFullStressVector(const Vector& rStressVector);

private:
    /// Serialization
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }

}; // Class UPwSmallStrainInterfaceElement

} // namespace Kratos