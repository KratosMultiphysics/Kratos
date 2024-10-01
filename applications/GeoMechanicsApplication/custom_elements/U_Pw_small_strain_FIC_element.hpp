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

#if !defined(KRATOS_GEO_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwSmallStrainFICElement
    : public UPwSmallStrainElement<TDim, TNumNodes>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwSmallStrainFICElement);

    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;
    using UPwBaseElement::mConstitutiveLawVector;
    using UPwBaseElement::mStateVariablesFinalized;
    using UPwBaseElement::mStressVector;
    using UPwBaseElement::mThisIntegrationMethod;

    using ElementVariables = typename UPwSmallStrainElement<TDim, TNumNodes>::ElementVariables;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    explicit UPwSmallStrainFICElement(IndexType NewId = 0)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId)
    {
    }

    /// Constructor using an array of nodes
    UPwSmallStrainFICElement(IndexType NewId, const NodesArrayType& ThisNodes, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, ThisNodes, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Geometry
    UPwSmallStrainFICElement(IndexType NewId, GeometryType::Pointer pGeometry, std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, pGeometry, std::move(pStressStatePolicy))
    {
    }

    /// Constructor using Properties
    UPwSmallStrainFICElement(IndexType                          NewId,
                             GeometryType::Pointer              pGeometry,
                             PropertiesType::Pointer            pProperties,
                             std::unique_ptr<StressStatePolicy> pStressStatePolicy)
        : UPwSmallStrainElement<TDim, TNumNodes>(NewId, pGeometry, pProperties, std::move(pStressStatePolicy))
    {
    }

    ~UPwSmallStrainFICElement() = default;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    // Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        std::stringstream claw_buffer;
        if (mConstitutiveLawVector.size() != 0) {
            claw_buffer << mConstitutiveLawVector[0]->Info();
        } else {
            claw_buffer << "not defined";
        }
        buffer << "U-Pw smal strain FIC Element #" << this->Id()
               << "\nConstitutive law: " << claw_buffer.str();
        return buffer.str();
    }

    // Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    struct FICElementVariables {
        /// Properties variables
        double ShearModulus;

        /// Nodal variables
        array_1d<array_1d<double, TDim * TNumNodes>, TNumNodes> NodalShapeFunctionsGradients;

        /// General elemental variables
        double ElementLength;
        Matrix VoigtMatrix;

        /// Variables computed at each GP
        BoundedMatrix<double, TDim, TDim * TNumNodes>       StrainGradients;
        array_1d<Vector, TNumNodes>                         ShapeFunctionsSecondOrderGradients;
        array_1d<array_1d<double, TDim>, TDim>              DtStressGradients;
        array_1d<std::vector<array_1d<double, TDim>>, TDim> ConstitutiveTensorGradients;
        /// Auxiliary variables
        array_1d<double, TDim>                        DimVector;
        Matrix                                        DimVoigtMatrix;
        BoundedMatrix<double, TDim, TDim * TNumNodes> DimUMatrix;
    };

    /// Member Variables

    array_1d<std::vector<array_1d<double, TNumNodes>>, TDim> mNodalConstitutiveTensor;

    array_1d<array_1d<double, TNumNodes>, TDim> mNodalDtStress;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    double CalculateShearModulus(const Matrix& ConstitutiveMatrix) const;

    void SaveGPConstitutiveTensor(array_1d<Matrix, TDim>& rConstitutiveTensorContainer,
                                  const Matrix&           ConstitutiveMatrix,
                                  const unsigned int&     GPoint);

    void SaveGPDtStress(Matrix& rDtStressContainer, const Vector& StressVector, const unsigned int& GPoint);

    void ExtrapolateGPConstitutiveTensor(const array_1d<Matrix, TDim>& ConstitutiveTensorContainer);

    void ExtrapolateGPDtStress(const Matrix& DtStressContainer);

    void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                      VectorType&        rRightHandSideVector,
                      const ProcessInfo& CurrentProcessInfo,
                      bool               CalculateStiffnessMatrixFlag,
                      bool               CalculateResidualVectorFlag) override;

    void InitializeFICElementVariables(FICElementVariables& rFICVariables,
                                       const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer,
                                       const GeometryType&   Geom,
                                       const PropertiesType& Prop,
                                       const ProcessInfo&    CurrentProcessInfo);

    void ExtrapolateShapeFunctionsGradients(array_1d<array_1d<double, TDim * TNumNodes>, TNumNodes>& rNodalShapeFunctionsGradients,
                                            const GeometryType::ShapeFunctionsGradientsType& DN_DXContainer);

    void CalculateElementLength(double& rElementLength, const GeometryType& Geom);

    void InitializeSecondOrderTerms(FICElementVariables& rFICVariables);

    void CalculateShapeFunctionsSecondOrderGradients(FICElementVariables& rFICVariables,
                                                     ElementVariables&    rVariables);

    void CalculateAndAddLHSStabilization(MatrixType&          rLeftHandSideMatrix,
                                         ElementVariables&    rVariables,
                                         FICElementVariables& rFICVariables);

    void CalculateAndAddStrainGradientMatrix(MatrixType&                rLeftHandSideMatrix,
                                             const ElementVariables&    rVariables,
                                             const FICElementVariables& rFICVariables);

    void CalculateAndAddDtStressGradientMatrix(MatrixType&             rLeftHandSideMatrix,
                                               const ElementVariables& rVariables,
                                               FICElementVariables&    rFICVariables);

    void CalculateConstitutiveTensorGradients(FICElementVariables&    rFICVariables,
                                              const ElementVariables& Variables);

    void CalculateAndAddPressureGradientMatrix(MatrixType&                rLeftHandSideMatrix,
                                               const ElementVariables&    rVariables,
                                               const FICElementVariables& rFICVariables);

    void CalculateAndAddRHSStabilization(VectorType&          rRightHandSideVector,
                                         ElementVariables&    rVariables,
                                         FICElementVariables& rFICVariables);

    void CalculateAndAddStrainGradientFlow(VectorType&                rRightHandSideVector,
                                           const ElementVariables&    rVariables,
                                           const FICElementVariables& rFICVariables);

    void CalculateAndAddDtStressGradientFlow(VectorType&             rRightHandSideVector,
                                             const ElementVariables& rVariables,
                                             FICElementVariables&    rFICVariables);

    void CalculateDtStressGradients(FICElementVariables& rFICVariables, const ElementVariables& Variables);

    void CalculateAndAddPressureGradientFlow(VectorType&                rRightHandSideVector,
                                             const ElementVariables&    rVariables,
                                             const FICElementVariables& rFICVariables);

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    /// Member Variables

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Assignment operator.
    UPwSmallStrainFICElement& operator=(UPwSmallStrainFICElement const& rOther);

    /// Copy constructor.
    UPwSmallStrainFICElement(UPwSmallStrainFICElement const& rOther);

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

}; // Class UPwSmallStrainFICElement

} // namespace Kratos

#endif // KRATOS_GEO_U_PW_SMALL_STRAIN_FIC_ELEMENT_H_INCLUDED  defined
