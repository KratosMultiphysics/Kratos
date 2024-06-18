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

#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_retention/retention_law.h"
#include "geometries/geometry_data.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/kratos_export_api.h"
#include "includes/serializer.h"
#include "includes/smart_pointers.h"
#include "includes/ublas_interface.h"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUPwDiffOrderElement : public UPwBaseElement
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SmallStrainUPwDiffOrderElement);

    using UPwBaseElement::mConstitutiveLawVector;
    using UPwBaseElement::mIsInitialised;
    using UPwBaseElement::mRetentionLawVector;
    using UPwBaseElement::mStateVariablesFinalized;
    using UPwBaseElement::mStressVector;

    using UPwBaseElement::UPwBaseElement;

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   PropertiesType::Pointer            pProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    ~SmallStrainUPwDiffOrderElement() override = default;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    using Element::SetValuesOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                      std::vector<int>&    rValues,
                                      const ProcessInfo&   rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>&    rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    using Element::CalculateOnIntegrationPoints;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "U-Pw small strain different order Element #" << Id()
               << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "U-Pw small strain different order Element #" << Id()
                 << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

protected:
    struct ElementVariables {
        // Variables at all integration points
        Matrix                                    NuContainer;
        Matrix                                    NpContainer;
        GeometryType::ShapeFunctionsGradientsType DNu_DXContainer;
        GeometryType::ShapeFunctionsGradientsType DNp_DXContainer;
        Vector                                    detJuContainer;

        // Variables at each integration point
        Vector Nu;     // Contains the displacement shape functions at every node
        Vector Np;     // Contains the pressure shape functions at every node
        Matrix DNu_DX; // Contains the global derivatives of the displacement shape functions
        Matrix DNu_DXInitialConfiguration; // Contains the global derivatives of the displacement shape functions

        Matrix DNp_DX; // Contains the global derivatives of the pressure shape functions
        Matrix B;
        double IntegrationCoefficient;
        double IntegrationCoefficientInitialConfiguration;
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;

        // Variables needed for consistency with the general constitutive law
        Matrix F;

        // needed for updated Lagrangian:
        double detJ;                     // displacement
        double detJInitialConfiguration; // displacement

        // Nodal variables
        Vector BodyAcceleration;
        Vector DisplacementVector;
        Vector VelocityVector;
        Vector PressureVector;
        Vector DeltaPressureVector;
        Vector PressureDtVector;

        /// Retention Law parameters
        double DegreeOfSaturation;
        double DerivativeOfSaturation;
        double RelativePermeability;
        double BishopCoefficient;

        // Properties and processinfo variables
        bool IgnoreUndrained;
        bool UseHenckyStrain;
        bool ConsiderGeometricStiffness;

        // stress/flow variables
        double BiotCoefficient;
        double BiotModulusInverse;
        double DynamicViscosityInverse;
        Matrix IntrinsicPermeability;
        double VelocityCoefficient;
        double DtPressureCoefficient;
    };

    void CalculateMaterialStiffnessMatrix(MatrixType&        rStiffnessMatrix,
                                          const ProcessInfo& CurrentProcessInfo) override;

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void InitializeNodalVariables(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);

    virtual void CalculateKinematics(ElementVariables& rVariables, unsigned int GPoint);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables) const;

    void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const;

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const;

    void CalculateAndAddCompressibilityMatrix(MatrixType&             rLeftHandSideMatrix,
                                              const ElementVariables& rVariables) const;

    void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint);

    void CalculateAndAddStiffnessForce(VectorType&             rRightHandSideVector,
                                       const ElementVariables& rVariables,
                                       unsigned int            GPoint);

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables);

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const;

    void CalculateAndAddCompressibilityFlow(VectorType&             rRightHandSideVector,
                                            const ElementVariables& rVariables) const;

    [[nodiscard]] std::vector<double> CalculateRelativePermeabilityValues(const std::vector<double>& rFluidPressures) const;
    [[nodiscard]] std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const;
    void CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const;

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables);

    Matrix CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const;
    std::vector<Matrix> CalculateBMatrices(const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer,
                                           const Matrix& rNContainer) const;

    void AssignPressureToIntermediateNodes();

    virtual Vector CalculateGreenLagrangeStrain(const Matrix& rDeformationGradient) const;

    Matrix              CalculateDeformationGradient(unsigned int GPoint) const;
    std::vector<Matrix> CalculateDeformationGradients() const;

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

    Vector GetPressureSolutionVector();

    [[nodiscard]] std::vector<double> CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures);
    [[nodiscard]] std::vector<double> CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures);
    [[nodiscard]] virtual std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors) const;

    [[nodiscard]] SizeType GetNumberOfDOF() const override;

private:
    [[nodiscard]] DofsVectorType GetDofs() const;

    GeometryType::Pointer mpPressureGeometry;

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, UPwBaseElement)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, UPwBaseElement)
    }

    // Private Operations

    template <class TValueType>
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }

}; // Class SmallStrainUPwDiffOrderElement

} // namespace Kratos
