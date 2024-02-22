// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_GEO_STRUCTURAL_BASE_ELEMENT_H_INCLUDED)
#define KRATOS_GEO_STRUCTURAL_BASE_ELEMENT_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/variables.h"

// Application includes
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoStructuralBaseElement : public Element
{
public:
    /// The definition of the sizetype
    typedef std::size_t SizeType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoStructuralBaseElement);

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    GeoStructuralBaseElement(IndexType NewId = 0) : Element(NewId) {}

    /// Constructor using an array of nodes
    GeoStructuralBaseElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    GeoStructuralBaseElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    GeoStructuralBaseElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    /// Destructor
    virtual ~GeoStructuralBaseElement() {}

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType               NewId,
                            NodesArrayType const&   ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void SetValuesOnIntegrationPoints(const Variable<double>&    rVariable,
                                      const std::vector<double>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    static constexpr SizeType N_DOF_NODE    = (TDim == 2 ? 3 : 6);
    static constexpr SizeType N_DOF_ELEMENT = N_DOF_NODE * TNumNodes;
    static constexpr SizeType VoigtSize = (TDim == 3 ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRESS);

    /// Member Variables
    struct ElementVariables {
        /// Properties variables

        /// ProcessInfo variables

        /// Nodal variables
        array_1d<double, TNumNodes * TDim> DisplacementVector;
        array_1d<double, TNumNodes * TDim> VelocityVector;
        array_1d<double, TNumNodes * TDim> NodalVolumeAcceleration;
        array_1d<double, TNumNodes * TDim> UVector;

        Vector DofValuesVector;

        /// General elemental variables
        Matrix NodalCrossDirection;

        /// Variables computed at each GP
        Matrix                                        B;
        BoundedMatrix<double, TDim, TNumNodes * TDim> NuTot;

        Matrix                 TransformationMatrix;
        array_1d<double, TDim> GaussVolumeAcceleration;
        double                 IntegrationCoefficient;
        /// Constitutive Law parameters
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        Vector Nu;
        Matrix GradNe;
        Matrix F;
        double detF;
        double HalfThickness;

        /// Auxiliary Variables
        Matrix UVoigtMatrix;
    };

    GeometryData::IntegrationMethod mThisIntegrationMethod;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<Vector>                   mStressVector;

    virtual SizeType GetTotalNumberIntegrationPoints() const;
    virtual SizeType GetCrossNumberIntegrationPoints() const;
    virtual SizeType GetAlongNumberIntegrationPoints() const;

    virtual void InitializeElementVariables(ElementVariables&            rVariables,
                                            ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                            const GeometryType&          Geom,
                                            const PropertiesType&        Prop,
                                            const ProcessInfo&           rCurrentProcessInfo) const;

    virtual void GetNodalDofValuesVector(Vector&             rNodalVariableVector,
                                         const GeometryType& Geom,
                                         IndexType           SolutionStepIndex = 0) const;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateStiffnessMatrix(MatrixType& rStiffnessMatrix, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo,
                              const bool         CalculateStiffnessMatrixFlag,
                              const bool         CalculateResidualVectorFlag);

    virtual void CalculateRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateNodalCrossDirection(Matrix& NodalCrossDirection) const;

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    /// Assignment operator.
    GeoStructuralBaseElement& operator=(GeoStructuralBaseElement const& rOther);

    /// Copy constructor.
    GeoStructuralBaseElement(GeoStructuralBaseElement const& rOther);

    [[nodiscard]] DofsVectorType GetDofs() const;

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

}; // Class GeoStructuralBaseElement

} // namespace Kratos

#endif // KRATOS_GEO_STRUCTURAL_BASE_ELEMENT_H_INCLUDED  defined
