//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application


#if !defined(KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED)
#define KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/math_utils.h"
#include "geometries/nurbs_curve_geometry.h"


// Application includes
#include "iga_application_variables.h"

namespace Kratos
{

//class IsogeometricBeamElement
//    : public Element

class KRATOS_API(IGA_APPLICATION) IsogeometricBeamElement final
    : public Element

{
protected:
    using Vector3d = Kratos::array_1d<double, 3>;
    using Matrix3d = Kratos::BoundedMatrix<double, 3, 3>;
    using Matrix2d = Kratos::BoundedMatrix<double, 2, 2>;
    using Matrix32d = Kratos::BoundedMatrix<double, 3, 2>;
    using Matrix23d = Kratos::BoundedMatrix<double, 2, 3>;

    template <typename... Args>
    auto cross_prod(Args&&... args) {
        return MathUtils<double>::CrossProduct(std::forward<Args>(args)...);
    }

    struct KinematicVariables
    {
        // Reference configuration
        array_1d<double, 3> R1;     // Reference tangent vector
        array_1d<double, 3> R2;     // Reference curvature vector  
        array_1d<double, 3> R3;     // Reference third derivative
        double A;                    // Reference length measure
        double B;                    // Reference curvature measure
        
        // Current configuration
        array_1d<double, 3> r1;     // Current tangent vector
        array_1d<double, 3> r2;     // Current curvature vector
        array_1d<double, 3> r3;     // Current third derivative
        double a;                    // Current length measure
        double b;                    // Current curvature measure
        
        // Cross-section directors
        array_1d<double, 3> N0, V0;  // Reference directors
        array_1d<double, 3> n, v;    // Current directors
        
        // Curvature components
        double B_n, B_v;             // Reference curvatures
        double b_n, b_v;             // Current curvatures
        double C_12, C_13;           // Reference twist
        double c_12, c_13;           // Current twist
        
        // Rotations
        double Phi, Phi_der, Phi_der2;    // Reference rotation
        double phi, phi_der, phi_der2;    // Current rotation
        
        KinematicVariables(SizeType Dimension = 3)
        {
            R1 = ZeroVector(Dimension);
            R2 = ZeroVector(Dimension);
            R3 = ZeroVector(Dimension);
            r1 = ZeroVector(Dimension);
            r2 = ZeroVector(Dimension);
            r3 = ZeroVector(Dimension);
            N0 = ZeroVector(Dimension);
            V0 = ZeroVector(Dimension);
            n = ZeroVector(Dimension);
            v = ZeroVector(Dimension);
            A = B = a = b = 0.0;
            B_n = B_v = b_n = b_v = 0.0;
            C_12 = C_13 = c_12 = c_13 = 0.0;
            Phi = Phi_der = Phi_der2 = 0.0;
            phi = phi_der = phi_der2 = 0.0;
        }
    };

    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;
        
        ConstitutiveVariables(SizeType StrainSize = 3)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IsogeometricBeamElement);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Geometry<Node> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    IsogeometricBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    };

    /// Constructor using an array of nodes with properties
    IsogeometricBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    };

    /// Default constructor necessary for serialization
    IsogeometricBeamElement()
        : Element()
    {
    };

    /// Destructor.
    ~IsogeometricBeamElement() override = default;

    ///@}
    ///@name Element Creation
    ///@{

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<IsogeometricBeamElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<IsogeometricBeamElement>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    };


    ///@}
    ///@name Degrees of freedom
    ///@{

    /// Relates the degrees of freedom of the element geometry
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /// Sets the ID's of the element degrees of freedom
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    ///@}
    ///@name Element Calculations
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeMaterial();

    //Computes RHS
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 4;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        MatrixType left_hand_side_matrix;

        CalculateAll(left_hand_side_matrix, rRightHandSideVector,
            rCurrentProcessInfo, false, true);
    }

    //Computes LHS
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 4;

        VectorType right_hand_side_vector;
         
        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, right_hand_side_vector,
            rCurrentProcessInfo, true, false);
    }

    //Computes LHS and RHS
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const SizeType nb_nodes = GetGeometry().size();
        const SizeType nb_dofs = nb_nodes * 4;

        if (rRightHandSideVector.size() != nb_dofs)
            rRightHandSideVector.resize(nb_dofs);
        noalias(rRightHandSideVector) = ZeroVector(nb_dofs);

        if (rLeftHandSideMatrix.size1() != nb_dofs && rLeftHandSideMatrix.size2() != nb_dofs)
            rLeftHandSideMatrix.resize(nb_dofs, nb_dofs);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(nb_dofs, nb_dofs);

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
            rCurrentProcessInfo, true, true);
    }

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag);

    void CalculateKinematics(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematicVariables);

    void CalculateConstitutiveVariables(
        const IndexType IntegrationPointIndex,
        KinematicVariables& rKinematics,
        ConstitutiveVariables& rConstitutiveVars,
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure);

    void ComputeBMatrices(
        IndexType point_number,
        KinematicVariables& rKinematicVariables,
        Matrix& rBAxial,
        Matrix& rBBending1,
        Matrix& rBBending2,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2);

    void ComputeGMatrices(
        IndexType point_number,
        KinematicVariables& rKinematicVariables,
        Matrix& rGAxial,
        Matrix& rGBending1,
        Matrix& rGBending2,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2);

    /// Get Displacemnts
    void GetValuesVector(Vector& rValues, int Step) const override;

    ///@}
    ///@name Input and Output
    ///@{

    /// Check provided parameters
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "\"IsogeometricBeamElement\" #" << Id()
            << " with geometry #" << this->GetGeometry().Id() << " with center in: "
            << this->GetGeometry().Center() << std::endl;
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Isogeometric Beam Element #" << Id();
        return buffer.str();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }





    ///@}
private:

    double mTolerance = 1.0e-8;

    void CompMatLambda(Matrix3d& _mat_lambda, Vector3d _vec1, Vector3d _vec2);
    void CompMatLambdaDeriv(Matrix3d& _mat_lambda_der, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv);
    void CompMatLambdaDeriv2(Matrix3d& _mat_lambda_derder, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv, Vector3d _vec1_deriv2, Vector3d _vec2_deriv2);
    void CompMatLambdaVar(Matrix& _mat_lam_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var);
    void CompMatLambdaVarVar(Matrix& _mat_lam_var_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var, Matrix _vec2_var_var);
    void CompMatLambdaDerivVar(Matrix& _mat_lam_der_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var);
    void CompMatLambdaDeriv2Var(Matrix& _mat_lam_derder_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector3d _vec1_derder, Vector _vec2_var, Vector3d _vec2_der, Vector3d _vec2_derder, Vector _vec2_der_var, Vector _vec2_derder_var);
    void CompMatLambdaDerivVarVar(Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var);
    void CompMatLambdaAll(Matrix& _mat_lambda_var, Matrix& _mat_lam_der_var, Matrix& _mat_lam_var_var, Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var);
    
    void CompTVar(Vector& _t_var, Vector& _deriv, Vector3d& _r1);
    void CompTDerivVar(Vector& _t_deriv_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2);
    void CompTVarVar(Matrix& _t_var_var, Vector& _deriv, Vector3d& _r1);
    void CompTDerivVarVar(Matrix& _t_deriv_var_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2);
    void CompTDeriv2Var(Vector& _t_deriv2_var, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3);
    void CompGeometryReferenceCrossSection( KinematicVariables &kinematic_variables, Vector3d _R1, Vector3d _R2, Vector3d _T0_vec, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _n0, Vector3d& _v0, double& _B_n, double& _B_v, double& _C_12, double& _C_13, double& _Phi, double& _Phi_0_der);
    void CompGeometryActualCrossSection(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _N0, Vector3d& _V0, double& _b_n, double& _b_v, double& _c_12, double& _c_13, double _phi, double _phi_der, double _Phi, double _Phi_der);
    Vector CompPhiDof(Vector& _func);
    double GetDeltaPhi(KinematicVariables &kinematic_variables, Vector3d &n);
    void CompPhiRefProp(KinematicVariables &kinematic_variables, double& _Phi, double& _Phi_0_der);
    Vector CompEpsilonDof(Vector3d& _r1, Vector& _shape_func_deriv);
    Matrix CompEpsilonDof2(Vector3d& _r1, Vector& _shape_func_deriv);

    void CompMatRodrigues(Matrix3d& _mat_rod, Vector3d _vec, double _phi);
    void CompMatRodriguesDeriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, double _phi, double _phi_deriv);
    void CompMatRodriguesVar(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, double _phi);
    void comp_mat_rodrigues_var_var(Matrix& _mat_rod_var_var, Vector3d _vec, Vector _vec_var, Matrix _vec_var_var, Vector _func, double _phi);
    void comp_mat_rodrigues_deriv_var(Matrix& _mat_rod_der_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector _func, Vector _deriv, double _phi, double _phi_der);
    void comp_mat_rodrigues_deriv_var_var(Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der);
    void comp_mat_rodrigues_deriv2(Matrix3d& _mat_rod_derder, Vector3d _vec, Vector3d _vec_deriv, Vector3d _vec_deriv2, double _phi, double _phi_deriv, double _phi_deriv2);
    void comp_mat_rodrigues_deriv2_var(Matrix& _mat_rod_derder_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector3d _vec_derder, Vector _vec_derder_var, Vector _func, Vector _deriv, Vector _deriv2, double _phi, double _phi_der, double _phi_der2);
    void comp_mat_rodrigues_all(Matrix& _mat_rod_var, Matrix& _mat_rod_der_var, Matrix& _mat_rod_var_var, Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector3d _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der);

    Matrix3d cross_prod_vec_mat(const Vector3d& vec, const Matrix3d& mat)
    {
        Matrix3d mat_vec;
        mat_vec.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        mat_vec(i, j) += permutation[j][k][l] * mat(i, k) * vec(l);
                    }
                }
            }
        }
        return mat_vec;
    }

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    SizeType mNumberOfDofs;
    SizeType mDofsPerNode;
    Matrix mSMatRodVar;
    Matrix mSMatLamVar;
    Matrix mSMatLamVarRodLam;
    Matrix mSMatRodLamVarRodLam;
    Matrix mSMatRodVarLamRodLam;
    Matrix mSMatRodLamRodLamVar;
    Matrix mSMatRodDerVar;
    Matrix mSMatLamDerVar;
    Matrix mSMatLamDerVarRodLam;
    Matrix mSMatLamVarRodDerLam;
    Matrix mSMatLamVarRodLamDer;
    Matrix mSMatRodDerLamVarRodLam;
    Matrix mSMatRodLamDerVarRodLam;
    Matrix mSMatRodLamVarRodDerLam;
    Matrix mSMatRodLamVarRodLamDer;
    Matrix mSMatRodDerVarLamRodLam;
    Matrix mSMatRodVarLamDerRodLam;
    Matrix mSMatRodVarLamRodDerLam;
    Matrix mSMatRodVarLamRodLamDer;
    Matrix mSMatRodLamRodLamDerVar;
    Matrix mSMatRodVarVar;
    Matrix mSMatLamVarVar;
    Matrix mSMatLamVarVarRodLam;
    Matrix mSMatRodDerVarVar;
    Matrix mSMatLamDerVarVar;
    Matrix mSMatLamDerVarVarRodLam;
    Matrix mSMatLamVarVarRodDerLam;
    Matrix mSMatLamVarVarRodLamDer;

    Vector3d N;
    Vector3d n;
    Vector3d V;
    Vector3d v;
};

} // namespace Kratos

#endif // !defined(KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED)
