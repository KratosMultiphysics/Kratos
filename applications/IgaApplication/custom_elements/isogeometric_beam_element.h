//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application


//#if !defined(KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED)
//#define KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED

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
    //using Matrix3d = BoundedMatrix<double, 3, 3>;
    //using Matrix2d = BoundedMatrix<double, 2, 2>;
    //using Matrix32d = BoundedMatrix<double, 3, 2>;
    //using Matrix23d = BoundedMatrix<double, 2, 3>;
    //using Vector3d = BoundedVector<double, 3>;
    using Vector3d = Kratos::array_1d<double, 3>;
    using Matrix3d = Kratos::BoundedMatrix<double, 3, 3>;
    using Matrix2d = Kratos::BoundedMatrix<double, 2, 2>;
    using Matrix32d = Kratos::BoundedMatrix<double, 3, 2>;
    using Matrix23d = Kratos::BoundedMatrix<double, 2, 3>;
    //using Matrix23d = Kratos::matrix
    
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
        Vector StrainVector;         // [E¹¹, E¹², E¹³] - embedded beam formulation
        Vector StressVector;         // [S¹¹, S¹², S¹³] - stress resultants
        Matrix ConstitutiveMatrix;   // Material stiffness matrix
        
        ConstitutiveVariables(SizeType StrainSize = 3)
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    struct SecondVariations
    {
        Matrix Epsilon;    // Membrane strain
        Matrix Kappa_n;    // Bending about n
        Matrix Kappa_v;    // Bending about v
        Matrix Gamma_n;    // Torsion about n
        Matrix Gamma_v;    // Torsion about v
        
        SecondVariations(const SizeType mat_size)
        {
            Epsilon = ZeroMatrix(mat_size, mat_size);
            Kappa_n = ZeroMatrix(mat_size, mat_size);
            Kappa_v = ZeroMatrix(mat_size, mat_size);
            Gamma_n = ZeroMatrix(mat_size, mat_size);
            Gamma_v = ZeroMatrix(mat_size, mat_size);
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
    ///@name Life Cycle
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
    ///@name Analysis stages
    ///@{

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeMaterial();
    //new:
    void set_Memory();

    void comp_mat_lambda(Matrix3d& _mat_lambda, Vector3d _vec1, Vector3d _vec2);
    void comp_mat_lambda_deriv(Matrix3d& _mat_lambda_der, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv);
    void comp_mat_lambda_deriv2(Matrix3d& _mat_lambda_derder, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv, Vector3d _vec1_deriv2, Vector3d _vec2_deriv2);
    void comp_mat_lambda_var(Matrix& _mat_lam_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var);
    void comp_mat_lambda_var_var(Matrix& _mat_lam_var_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var, Matrix _vec2_var_var);
    void comp_mat_lambda_deriv_var(Matrix& _mat_lam_der_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var);
    void comp_mat_lambda_deriv2_var(Matrix& _mat_lam_derder_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector3d _vec1_derder, Vector _vec2_var, Vector3d _vec2_der, Vector3d _vec2_derder, Vector _vec2_der_var, Vector _vec2_derder_var);
    void comp_mat_lambda_deriv_var_var(Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var);
    void comp_mat_lambda_all(Matrix& _mat_lambda_var, Matrix& _mat_lam_der_var, Matrix& _mat_lam_var_var, Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var);
    
    
    void comp_T_var(Vector& _t_var, Vector& _deriv, Vector3d& _r1);
    void comp_T_deriv_var(Vector& _t_deriv_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2);
    void comp_T_var_var(Matrix& _t_var_var, Vector& _deriv, Vector3d& _r1);
    void comp_T_deriv_var_var(Matrix& _t_deriv_var_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2);
    void comp_T_deriv2_var(Vector& _t_deriv2_var, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3);

    // void comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, double& _a, double& _b);
    // void comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, double& _a, double& _b);
    // void comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector3d& _R1, Vector3d& _R2, double& _A_ref, double& _B_ref);
    // void comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector _deriv3, Vector3d& _R1, Vector3d& _R2, Vector3d& _R3, double& _A_ref, double& _B_ref);
    void comp_Geometry_reference_cross_section( Vector3d _R1, Vector3d _R2, Vector3d _T0_vec, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _n0, Vector3d& _v0, double& _B_n, double& _B_v, double& _C_12, double& _C_13, double& _Phi, double& _Phi_0_der);
    // void comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo,Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, double& _a, double& _b);
    // void comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo,Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, double& _a, double& _b);
    void comp_Geometry_actual_cross_section(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _N0, Vector3d& _V0, double& _b_n, double& _b_v, double& _c_12, double& _c_13, double _phi, double _phi_der, double _Phi, double _Phi_der);

    //void get_Dof_Types_Per_Node(std::vector<dof_type>& _act_dofs);
    Vector comp_phi_dof(Vector& _func);
    void comp_Phi_ref_prop(double& _Phi, double& _Phi_0_der);
    Vector comp_epsilon_dof(Vector3d& _r1, Vector& _shape_func_deriv);
    Matrix comp_epsilon_dof_2(Vector3d& _r1, Vector& _shape_func_deriv);

    void comp_mat_rodrigues(Matrix3d& _mat_rod, Vector3d _vec, double _phi);
    void comp_mat_rodrigues_deriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, double _phi, double _phi_deriv);
    void comp_mat_rodrigues_var(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, double _phi);
    void comp_mat_rodrigues_var_var(Matrix& _mat_rod_var_var, Vector3d _vec, Vector _vec_var, Matrix _vec_var_var, Vector _func, double _phi);
    void comp_mat_rodrigues_deriv_var(Matrix& _mat_rod_der_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector _func, Vector _deriv, double _phi, double _phi_der);
    void comp_mat_rodrigues_deriv_var_var(Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der);
    void comp_mat_rodrigues_deriv2(Matrix3d& _mat_rod_derder, Vector3d _vec, Vector3d _vec_deriv, Vector3d _vec_deriv2, double _phi, double _phi_deriv, double _phi_deriv2);
    void comp_mat_rodrigues_deriv2_var(Matrix& _mat_rod_derder_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector3d _vec_derder, Vector _vec_derder_var, Vector _func, Vector _deriv, Vector _deriv2, double _phi, double _phi_der, double _phi_der2);
    void comp_mat_rodrigues_all(Matrix& _mat_rod_var, Matrix& _mat_rod_der_var, Matrix& _mat_rod_var_var, Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector3d _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der);

    // void comp_dof_lin(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Vector& _shear_var_n, Vector& _shear_var_v, Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _r3, Vector3d& _R3, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, double _phi, double _phi_der, double _phi_der2, double _Phi, double _Phi_der, double _Phi_der2);
    // void comp_dof_lin(
    //     IndexType point_number,
    //     Vector& _cur_var_n,
    //     Vector& _cur_var_v,
    //     Vector& _tor_var_n,
    //     Vector& _tor_var_v,
    //     const KinematicVariables& rKinematicVariables);
    //void comp_dof_nln(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Matrix& _cur_var_n_2, Matrix& _cur_var_v_2, Matrix& _tor_var_n_2, Matrix& _tor_var_v_2,  Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, double _phi, double _phi_der, double _Phi, double _Phi_der);

    //void parametric_mapping(Vector& _knot_vector, double  u, double& u_mid);

    //void stiff_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie, double& _dL);
    //void stiff_mat_el_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie, double& _dL);
    //void mass_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, Matrix& _me);
    //void stiff_mat_el_geo(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke,  double& _dL);
    
    //void calc_Geo_Lin_Stiff(const ProcessInfo& rCurrentProcessInfo, Matrix& _stiffness_mtx, Vector& _f_Int);
    
    //void comp_transverse_shear_force_nln(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d _r3, Vector3d _R3,  Vector3d& _N0, Vector3d& _V0, double _phi, double _phi_der, double _phi_der2, double _Phi, double _Phi_der, double _Phi_der2, double& _shear_force_n, double& _shear_force_v);

    //void stress_res_lin(const ProcessInfo& rCurrentProcessInfo,  IndexType integration_point_index, Vector3d& _f, Vector3d& _m);
    //void stress_res_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Vector3d& _f, Vector3d& _m);


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
    /// Updates the constitutive law
    //void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Load functions
    ///@{

    //void CalculateBodyForces(Vector& rBodyForces);

    //bool HasSelfWeight() const;

    ///@}
    ///@name Explicit dynamic functions
    ///@{

    //void AddExplicitContribution(const VectorType& rRHSVector,const Variable<VectorType>& rRHSVariable,const Variable<double >& rDestinationVariable,const ProcessInfo& rCurrentProcessInfo) override;

    //void AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable,const Variable<array_1d<double, 3>>& rDestinationVariable,const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Dynamic functions
    ///@{

    //void CalculateDampingMatrix(MatrixType& rDampingMatrix,const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculates the mass matrix with use of the lumped mass vector
    //void CalculateMassMatrix(MatrixType& rMassMatrix,const ProcessInfo& rCurrentProcessInfo) override;

    /// Calculates lumped mass vector
    //void CalculateLumpedMassVector( VectorType& rLumpedMassVector,const ProcessInfo& rCurrentProcessInfo) const override;

    /// Get Displacemnts
    void GetValuesVector(Vector& rValues, int Step) const override;

    /// Get Velocities
    void GetFirstDerivativesVector(Vector& rValues,int Step = 0) const override;

    /// Get Accelerations
    void GetSecondDerivativesVector(Vector& rValues,int Step = 0) const override;

    ///@}
    ///@name Output functions
    ///@{

//void CalculateOnIntegrationPoints(const Variable<double>& rVariable,std::vector<double>& rOutput,const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Info
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

    double Tol = 1.0e-8;

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

    NurbsCurveGeometry<3, PointerVector<NodeType>>::Pointer pS_mat_rod_varCurve;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    SizeType N_Dof;
    SizeType Dof_Node;
    Matrix _gke;
    Vector _gfie;
    Matrix S_gke;
    Vector S_fie;
    Vector S_dL;
    Vector S_eps_var;
    Vector S_curv_n_var;
    Vector S_curv_v_var;
    Vector S_torsion_var;
    Vector S_torsion_n_var;
    Vector S_torsion_v_var;
    Vector S_shear_n_var;
    Vector S_shear_v_var;
    Matrix S_eps_var_var;
    Matrix S_curv_n_var_var;
    Matrix S_curv_v_var_var;
    Matrix S_torsion_var_var;
    Matrix S_torsion_n_var_var;
    Matrix S_torsion_v_var_var;
    Matrix S_kem;           //membrane stiffness
    Matrix S_keb_n;           //bending stiffness
    Matrix S_keb_v;           //bending stiffness
    Matrix S_ket_n;           //torsional stiffness
    Matrix S_ket_v;           //torsional stiffness
    Matrix S_me;            //mass matrix
    Vector S_fiem;          //internal force (membrane)
    Vector S_fieb_n;          //internal force (bending)
    Vector S_fieb_v;          //internal force (bending)
    Vector S_fiet_n;          //internal force (torsion)
    Vector S_fiet_v;          //internal force (torsion)
    Vector S_R;             //NURBS function
    Vector S_dR;            //1st derivatives
    Vector S_ddR;           //2nd derivatives
    //first variation of mapping matrices
    Matrix S_mat_rod_var;
    Matrix S_mat_lam_var;
    Matrix S_mat_lam_var_Rod_Lam;
    Matrix S_mat_rod_lam_var_Rod_Lam;
    Matrix S_mat_rod_var_lam_Rod_Lam;
    Matrix S_mat_rodlamRodLam_var;
    Matrix S_mat_rod_der_var;
    Matrix S_mat_lam_der_var;
    Matrix S_mat_lam_der_var_Rod_Lam;
    Matrix S_mat_lam_var_Rod_der_Lam;
    Matrix S_mat_lam_var_Rod_Lam_der;
    Matrix S_mat_rod_der_lam_var_Rod_Lam;
    Matrix S_mat_rod_lam_der_var_Rod_Lam;
    Matrix S_mat_rod_lam_var_Rod_der_Lam;
    Matrix S_mat_rod_lam_var_Rod_Lam_der;
    Matrix S_mat_rod_der_var_lam_Rod_Lam;
    Matrix S_mat_rod_var_lam_der_Rod_Lam;
    Matrix S_mat_rod_var_lam_Rod_der_Lam;
    Matrix S_mat_rod_var_lam_Rod_Lam_der;
    Matrix S_mat_rodlamRodLam_der_var;
    //second variation of mapping matrices
    Matrix S_mat_rod_var_var;
    Matrix S_mat_lam_var_var;
    Matrix S_mat_lam_var_var_Rod_Lam;
    Matrix S_mat_rod_lam_var_var_Rod_Lam;
    Matrix S_mat_rod_var_lam_var_Rod_Lam;
    Matrix S_mat_rod_var_var_lam_Rod_Lam;
    Matrix S_mat_rodlamRodLam_var_var;
    Matrix S_mat_rod_der_var_var;
    Matrix S_mat_lam_der_var_var;
    Matrix S_mat_lam_der_var_var_Rod_Lam;
    Matrix S_mat_lam_var_var_Rod_der_Lam;
    Matrix S_mat_lam_var_var_Rod_Lam_der;
    Vector S_U_Vec; 
    Vector S_weights;
    Vector3d N;
    Vector3d n;
    Vector3d V;
    Vector3d v;

    ///@name Private member
    ///@{


    std::vector<array_1d<double, 3>> mReferenceBaseVector;

    ///@}
    ///@name Private operations
    ///@{

    /// Initializes constitutive law vector and materials.
    //void InitializeMaterial();

    /// Computes the base vector at a integration point position.
    //array_1d<double, 3> CalculateActualBaseVector(IndexType IntegrationPointIndex) const;

    /// Computes Green Lagrange Strain for all integration points
    //void CalculateGreenLagrangeStrain(std::vector<double>& rGreenLagrangeVector) const;

    //void CalculateTangentModulus(std::vector<double>& rTangentModulusVector, const ProcessInfo& rCurrentProcessInfo);


    /// Computes prestress
    //double CalculatePrestressPK2(double reference_a, double actual_a) const;

    /// Computes PK2 stress
    //void CalculateStressPK2(std::vector<double>& rStressVector,const ProcessInfo& rCurrentProcessInfo) const;

    /// Computes cauchy stress
    //void CalculateStressCauchy(std::vector<double>& rStressVector,const ProcessInfo& rCurrentProcessInfo) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IsogeometricBeamElement);
        rSerializer.save("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IsogeometricBeamElement);
        rSerializer.load("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    }

    ///@}
};

} // namespace Kratos

//#endif // !defined(KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED)
