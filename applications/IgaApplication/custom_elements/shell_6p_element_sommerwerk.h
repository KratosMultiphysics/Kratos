//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt


#if !defined(KRATOS_SHELL_6P_ELEMENT_SOMMERWERK_H_INCLUDED)
#define KRATOS_SHELL_6P_ELEMENT_SOMMERWERK_H_INCLUDED

#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"

#include "iga_application_variables.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) Shell6pElement_Sommerwerk final
    : public Element
{
public:

    using Matrix2d  = BoundedMatrix<double, 2, 2>;
    using Matrix23d = BoundedMatrix<double, 2, 3>;
    using Matrix3d  = BoundedMatrix<double, 3, 3>;


    struct PointGeometry
    {
        // Covariant base vectors a_alpha = phi_,alpha  (Diss. Gl. 2.34 / A.1)
        BoundedVector<double, 3> a1;
        BoundedVector<double, 3> a2;

        // Contravariant base vectors a^alpha = a^{alpha beta} a_beta  (Diss. Gl. A.8)
        BoundedVector<double, 3> aSup1;
        BoundedVector<double, 3> aSup2;

        // Second derivatives of the surface
        BoundedVector<double, 3> d11;
        BoundedVector<double, 3> d12;
        BoundedVector<double, 3> d22;

        // Rows 0,1 = local Cartesian frame e_1,e_2; row 2 = unit normal a_3 
        Matrix3d localCoordinateSystem;

        Matrix2d acov;   // first fundamental form
        Matrix2d acont;  // a^{alpha beta} = (a_{alpha beta})^{-1}

        Matrix2d bcov;     // second fundamental form
        Matrix2d bcovcont; // mixed curvature b^lambda_alpha (Diss. Gl. A.10)

        Matrix2d c;        // third fundamental form

        //Christoffel symbols
        Matrix2d chr_xi1;
        Matrix2d chr_xi2;

        Matrix2d T;     //strain basis transformation contravariant->local Cartesian
        Matrix3d T_VE;  //Voigt form of T for the in-plane strains
        double dA;      //surface differential

        Matrix N_loc_deriv;
    };

public:

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell6pElement_Sommerwerk);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef Node NodeType;
    typedef Geometry<NodeType> GeometryType;

    static constexpr SizeType DOFsPerNode = 6;

    Shell6pElement_Sommerwerk(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    Shell6pElement_Sommerwerk(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    Shell6pElement_Sommerwerk()
        : Element()
    {};

    virtual ~Shell6pElement_Sommerwerk() final {};

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const final
    {
        return Kratos::make_intrusive<Shell6pElement_Sommerwerk>(NewId, pGeom, pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const final
    {
        return Kratos::make_intrusive<Shell6pElement_Sommerwerk>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    void CalculateRightHandSide(VectorType& rRHS, const ProcessInfo& rInfo) final
    {
        const SizeType mat_size = GetGeometry().size() * DOFsPerNode;
        if (rRHS.size() != mat_size) rRHS.resize(mat_size);
        noalias(rRHS) = ZeroVector(mat_size);
        MatrixType lhs_dummy;
        CalculateAll(lhs_dummy, rRHS, rInfo, false, true);
    }

    void CalculateLeftHandSide(MatrixType& rLHS, const ProcessInfo& rInfo) final
    {
        const SizeType mat_size = GetGeometry().size() * DOFsPerNode;
        VectorType rhs_dummy;
        if (rLHS.size1() != mat_size) rLHS.resize(mat_size, mat_size);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
        CalculateAll(rLHS, rhs_dummy, rInfo, true, false);
    }

    void CalculateLocalSystem(MatrixType& rLHS, VectorType& rRHS, const ProcessInfo& rInfo) final
    {
        const SizeType mat_size = GetGeometry().size() * DOFsPerNode;
        if (rRHS.size() != mat_size) rRHS.resize(mat_size);
        noalias(rRHS) = ZeroVector(mat_size);
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size)
            rLHS.resize(mat_size, mat_size);
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
        CalculateAll(rLHS, rRHS, rInfo, true, true);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rInfo) const final;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rInfo) const final;

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    void Initialize(const ProcessInfo& rInfo) final;

    void GetValuesVector(Vector& rValues, int Step) const final;
    void GetFirstDerivativesVector(Vector& rValues, int Step) const final;
    void GetSecondDerivativesVector(Vector& rValues, int Step) const final;

    int Check(const ProcessInfo& rInfo) const final;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    double ParametricAngleToEulerDegrees(IndexType IntegrationPointIndex, double BetaRadians) const;

    array_1d<double, 3> GetLocalTangent1AtIP(IndexType IntegrationPointIndex) const;

    const std::vector<ConstitutiveLaw::Pointer>& GetConstitutiveLawsAtIPs() const
    {
        return mConstitutiveLawVector;
    }

    std::string Info() const final
    {
        std::stringstream s;
        s << "Shell6pElement_Sommerwerk #" << Id();
        return s.str();
    }

    void PrintInfo(std::ostream& rOStream) const final { rOStream << "Shell6pElement_Sommerwerk #" << Id(); }
    void PrintData(std::ostream& rOStream) const final { pGetGeometry()->PrintData(rOStream); }

private:

    std::vector<PointGeometry> mPointGeometry;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;

    void InitPoint(IndexType IntegrationPointIndex);

    void CalculateResultantMatrices(
        IndexType IntegrationPointIndex,
        const ProcessInfo& rInfo,
        Matrix3d& rA,
        Matrix3d& rB,
        Matrix3d& rD,
        Matrix2d& rSh) const;

    void CalculateContravariantMaterial(
        IndexType IntegrationPointIndex,
        const ProcessInfo& rInfo,
        Matrix3d& rA,
        Matrix3d& rB,
        Matrix3d& rD,
        Matrix2d& rSh) const;

    void CalculateGeneralizedStrain(
        IndexType IntegrationPointIndex,
        Vector& rMembrane,
        Vector& rCurvature,
        Vector& rShear) const;

    void CalculateAll(
        MatrixType& rLHS,
        VectorType& rRHS,
        const ProcessInfo& rInfo,
        bool CalculateStiffnessFlag,
        bool CalculateResidualFlag);

    friend class Serializer;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

};

}

#endif
