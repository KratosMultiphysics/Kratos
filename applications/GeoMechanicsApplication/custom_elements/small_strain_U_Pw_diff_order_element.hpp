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
    using UPwBaseElement::mRetentionLawVector;
    using UPwBaseElement::mStateVariablesFinalized;
    using UPwBaseElement::mStressVector;

    using UPwBaseElement::UPwBaseElement;

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : UPwBaseElement(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
        SetUpPressureGeometryPointer();
    }

    SmallStrainUPwDiffOrderElement(IndexType                          NewId,
                                   GeometryType::Pointer              pGeometry,
                                   PropertiesType::Pointer            pProperties,
                                   std::unique_ptr<StressStatePolicy> pStressStatePolicy,
                                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : UPwBaseElement(NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
    {
        SetUpPressureGeometryPointer();
    }

    ~SmallStrainUPwDiffOrderElement() override = default;

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

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

    using Element::CalculateOnIntegrationPoints;

    void Calculate(const Variable<Vector>& rVariable, Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) override;
    using Element::Calculate;

    // Turn back information as a string.
    std::string Info() const override
    {
        const std::string constitutive_info =
            !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
        return "U-Pw small strain different order Element #" + std::to_string(Id()) +
               "\nConstitutive law: " + constitutive_info;
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

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
        double detJ;

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

    void CalculateAndAddGeometricStiffnessMatrix(
        MatrixType& rLeftHandSideMatrix, const Vector& rStressVector, const Matrix& rDNuDx, double IntegrationCoefficient) const;

    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void InitializeNodalVariables(ElementVariables& rVariables);

    void InitializeProperties(ElementVariables& rVariables);

    virtual void ExtractShapeFunctionDataAtIntegrationPoint(ElementVariables& rVariables, unsigned int GPoint);

    void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const;

    static void CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables);

    void CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ElementVariables& rVariables) const;

    static void CalculateAndAddCompressibilityMatrix(MatrixType&             rLeftHandSideMatrix,
                                                     const ElementVariables& rVariables);

    void CalculateAndAddStiffnessForce(VectorType&             rRightHandSideVector,
                                       const ElementVariables& rVariables,
                                       unsigned int            GPoint) const;

    void CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementVariables& rVariables) const;

    void CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const;

    static void CalculateAndAddCompressibilityFlow(VectorType&             rRightHandSideVector,
                                                   const ElementVariables& rVariables);

    [[nodiscard]] std::vector<double> CalculateBishopCoefficients(const std::vector<double>& rFluidPressures) const;
    static void CalculateAndAddPermeabilityFlow(VectorType&             rRightHandSideVector,
                                                const ElementVariables& rVariables);

    void CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, const ElementVariables& rVariables) const;

    Matrix CalculateBMatrix(const Matrix& rDN_DX, const Vector& rN) const;
    std::vector<Matrix> CalculateBMatrices(const GeometryType::ShapeFunctionsGradientsType& rDN_DXContainer,
                                           const Matrix& rNContainer) const;

    Vector GetPressures(size_t n_nodes) const;
    void   AssignPressureToIntermediateNodes();

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

    [[nodiscard]] Vector GetPressureSolutionVector() const;

    [[nodiscard]] std::vector<double> CalculateDegreesOfSaturation(const std::vector<double>& rFluidPressures);
    [[nodiscard]] std::vector<double> CalculateDerivativesOfSaturation(const std::vector<double>& rFluidPressures);
    [[nodiscard]] virtual std::vector<double> GetOptionalPermeabilityUpdateFactors(const std::vector<Vector>& rStrainVectors, bool ConsiderGeometricStiffness) const;

    [[nodiscard]] SizeType GetNumberOfDOF() const override;

private:
    GeometryType::Pointer mpPressureGeometry;

    [[nodiscard]] DofsVectorType GetDofs() const override;

    /**
     * @brief Sets the up the pressure geometry pointer object
     * This function sets the pointer for the auxiliary geometry for the pressure problem
     * The pressure geometry pointer is set according to the element geometry number of nodes and dimension
     */
    void SetUpPressureGeometryPointer();

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    Vector CalculateInternalForces(ElementVariables&          rVariables,
                                   const std::vector<Matrix>& rBMatrices,
                                   const std::vector<double>& rIntegrationCoefficients,
                                   const std::vector<double>& rBiotCoefficients,
                                   const std::vector<double>& rDegreesOfSaturation,
                                   const std::vector<double>& rBiotModuliInverse,
                                   const std::vector<double>& rRelativePermeabilityValues,
                                   const std::vector<double>& rBishopCoefficients) const;

    Vector CalculateExternalForces(ElementVariables&          rVariables,
                                   const std::vector<double>& rIntegrationCoefficients,
                                   const std::vector<double>& rIntegrationCoefficientsOnInitialConfiguration,
                                   const std::vector<double>& rDegreesOfSaturation,
                                   const std::vector<double>& rRelativePermeabilityValues,
                                   const std::vector<double>& rBishopCoefficients) const;
}; // Class SmallStrainUPwDiffOrderElement

} // namespace Kratos
