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

//class EmbeddedIsogeometricBeamElement
//    : public Element

class KRATOS_API(IGA_APPLICATION) EmbeddedIsogeometricBeamElement final
    : public Element

{
protected:
    using Matrix3d = BoundedMatrix<double, 3, 3>;
    using Matrix2d = BoundedMatrix<double, 2, 2>;
    using Matrix32d = BoundedMatrix<double, 3, 2>;
    using Matrix23d = BoundedMatrix<double, 2, 3>;
    using Vector3d = BoundedVector<double, 3>;
    
    
    enum class ConfigurationType {
        Current,
        Reference
    };

    template <typename... Args>
    auto cross_prod(Args&&... args) {
        return MathUtils<double>::CrossProduct(std::forward<Args>(args)...);
    }

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedIsogeometricBeamElement);

    /// Size types
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    // GometryType
    typedef Geometry<Node> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor using an array of nodes
    EmbeddedIsogeometricBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    };

    /// Constructor using an array of nodes with properties
    EmbeddedIsogeometricBeamElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    };

    /// Default constructor necessary for serialization
    EmbeddedIsogeometricBeamElement()
        : Element()
    {
    };

    /// Destructor.
    ~EmbeddedIsogeometricBeamElement() override = default;

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
        return Kratos::make_intrusive<EmbeddedIsogeometricBeamElement>(
            NewId, pGeom, pProperties);
    };

    /// Create with Id, pointer to geometry and pointer to property
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override
    {
        return Kratos::make_intrusive<EmbeddedIsogeometricBeamElement>(
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

    void comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, float& _a, float& _b);
    void comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, float& _a, float& _b);
    void comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector3d& _R1, Vector3d& _R2, float& _A_ref, float& _B_ref);
    void comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector _deriv3, Vector3d& _R1, Vector3d& _R2, Vector3d& _R3, float& _A_ref, float& _B_ref);
    void comp_Geometry_reference_cross_section( Vector3d _R1, Vector3d _R2, Vector3d _T0_vec, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _n0, Vector3d& _v0, float& _B_n, float& _B_v, float& _C_12, float& _C_13, float& _Phi, float& _Phi_0_der);
    void comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo,Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, float& _a, float& _b);
    void comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo,Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, float& _a, float& _b);
    void comp_Geometry_actual_cross_section(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _N0, Vector3d& _V0, float& _b_n, float& _b_v, float& _c_12, float& _c_13, float _phi, float _phi_der, float _Phi, float _Phi_der);

    //void get_Dof_Types_Per_Node(std::vector<dof_type>& _act_dofs);
    Vector comp_phi_dof(Vector& _func);
    void comp_Phi_ref_prop(float& _Phi, float& _Phi_0_der);
    Vector comp_epsilon_dof(Vector3d& _r1, Vector& _shape_func_deriv);
    Matrix comp_epsilon_dof_2(Vector3d& _r1, Vector& _shape_func_deriv);

    void comp_mat_rodrigues(Matrix3d& _mat_rod, Vector3d _vec, float _phi);
    void comp_mat_rodrigues_deriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, float _phi, float _phi_deriv);
    void comp_mat_rodrigues_var(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, float _phi);
    void comp_mat_rodrigues_var_var(Matrix& _mat_rod_var_var, Vector3d _vec, Vector _vec_var, Matrix _vec_var_var, Vector _func, float _phi);
    void comp_mat_rodrigues_deriv_var(Matrix& _mat_rod_der_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector _func, Vector _deriv, float _phi, float _phi_der);
    void comp_mat_rodrigues_deriv_var_var(Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, float _phi, float _phi_der);
    void comp_mat_rodrigues_deriv2(Matrix3d& _mat_rod_derder, Vector3d _vec, Vector3d _vec_deriv, Vector3d _vec_deriv2, float _phi, float _phi_deriv, float _phi_deriv2);
    void comp_mat_rodrigues_deriv2_var(Matrix& _mat_rod_derder_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector3d _vec_derder, Vector _vec_derder_var, Vector _func, Vector _deriv, Vector _deriv2, float _phi, float _phi_der, float _phi_der2);
    void comp_mat_rodrigues_all(Matrix& _mat_rod_var, Matrix& _mat_rod_der_var, Matrix& _mat_rod_var_var, Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector3d _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, float _phi, float _phi_der);

    void comp_dof_lin(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Vector& _shear_var_n, Vector& _shear_var_v, Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _r3, Vector3d& _R3, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, float _phi, float _phi_der, float _phi_der2, float _Phi, float _Phi_der, float _Phi_der2);
    void comp_dof_lin(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, float _phi, float _phi_der, float _Phi, float _Phi_der);
    void comp_dof_nln(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Matrix& _cur_var_n_2, Matrix& _cur_var_v_2, Matrix& _tor_var_n_2, Matrix& _tor_var_v_2,  Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, float _phi, float _phi_der, float _Phi, float _Phi_der);

    void parametric_mapping(Vector& _knot_vector, double  u, double& u_mid);

    void stiff_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie, float& _dL);
    void stiff_mat_el_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie, float& _dL);
    void mass_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, Matrix& _me);
    void stiff_mat_el_geo(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke,  float& _dL);
    
    void calc_Geo_Lin_Stiff(const ProcessInfo& rCurrentProcessInfo, Matrix& _stiffness_mtx, Vector& _f_Int);
    
    void comp_transverse_shear_force_nln(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d _r3, Vector3d _R3,  Vector3d& _N0, Vector3d& _V0, float _phi, float _phi_der, float _phi_der2, float _Phi, float _Phi_der, float _Phi_der2, float& _shear_force_n, float& _shear_force_v);

    void stress_res_lin(const ProcessInfo& rCurrentProcessInfo,  IndexType integration_point_index, Vector3d& _f, Vector3d& _m);
    void stress_res_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Vector3d& _f, Vector3d& _m);

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput,const ProcessInfo& rCurrentProcessInfo);
    void CalculateOnIntegrationPoints(const  Variable<Vector>& rVariable, std::vector <Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    void CalculateOnIntegrationPoints(const  Variable<double>& rVariable, std::vector <double>& rOutput, const ProcessInfo& rCurrentProcessInfo);
    
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
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide);

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
    void GetValuesVector(Vector& rValues,int Step = 0) const override;

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
        rOStream << "\"EmbeddedIsogeometricBeamElement\" #" << Id()
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

    float Tol = 1.0e-8;

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
   
    void GetElementOrientation(const Matrix& r_DN_De, const ConfigurationType& rConfiguration, Vector3d& B1, Vector3d& B2, Vector3d& B3, float& A, float& B);
    NurbsCurveGeometry<3, PointerVector<NodeType>>::Pointer pCurve;

    /// The vector containing the constitutive laws for all integration points.
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    SizeType _n_Dof;
    SizeType N_Dof;
    //list of matrices and vectors which size depends on the number of dofs
    SizeType Dof_Node;
    Matrix _gke;
    Vector _gfie;
    Matrix S_gke;           //Stiff matrix for Gausspoints
    Vector S_fie;           //internal force  (Gauss point)
    Vector S_dL;           //differential length
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
    Vector S_U_Vec;          //Knot vector in u-direction
    Vector S_weights;      //weights of the Control points

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, EmbeddedIsogeometricBeamElement);
        rSerializer.save("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, EmbeddedIsogeometricBeamElement);
        rSerializer.load("ReferenceBaseVector", mReferenceBaseVector);
        rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
    }

    void ExportMatrixToCSV(const MatrixType& matrix, const std::string& matrix_filename) {
        std::ofstream mat_file(matrix_filename);
        for (std::size_t i = 0; i < matrix.size1(); ++i) {
            for (std::size_t j = 0; j < matrix.size2(); ++j) {
                mat_file << matrix(i, j);
                if (j + 1 < matrix.size2()) mat_file << ",";
            }
            mat_file << "\n";
        }
        mat_file.close();
    }

    void ExportVectorToCSV(const VectorType& vector, const std::string& vector_filename) {
        std::ofstream vec_file(vector_filename);
        for (std::size_t i = 0; i < vector.size(); ++i) {
            vec_file << vector(i) << "\n";
        }
        vec_file.close();
    }

    ///@}
};

} // namespace Kratos

//#endif // !defined(KRATOS_ISOGEOMETRIC_BEAM_ELEMENT_H_INCLUDED)
