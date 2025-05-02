//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application

// System includes
// External includes
// Project includes
#include "custom_elements/isogeometric_beam_element.h"
#include <numeric>
#include "utilities/math_utils.h"
#include "geometries/nurbs_curve_geometry.h"

//debugging only
#include <sstream>
#include <iomanip>

namespace Kratos {
    ///@name Degrees of freedom
    ///@{

    void IsogeometricBeamElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        //KRATOS_WATCH("IsogeometricBeamElement::EquationIdVector");
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 4 * number_of_control_points)
            rResult.resize(4 * number_of_control_points, false);

        const IndexType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 4;
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(ROTATION_X, pos + 3).EquationId();
        }

        KRATOS_CATCH("")
    };

    void IsogeometricBeamElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        //KRATOS_WATCH("IsogeometricBeamElement::GetDofList");

        KRATOS_TRY;
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(4 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
        }
        KRATOS_CATCH("")
    };

    ///@}
    ///@name Analysis stages
    ///@{

    void IsogeometricBeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        KRATOS_CATCH("")
        //InitializeMaterial();
    }

    void IsogeometricBeamElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void IsogeometricBeamElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::CalculateAll");
        //KRATOS_WATCH(this->Id());
        set_Memory();
        const auto& r_geometry = GetGeometry();
        auto& r_integration_points = r_geometry.IntegrationPoints();
            
        int point_number = 0;
        float _dL;
        //// get integration data and weights, GP weight and dL and (u_vec[Knotspan_index+1]-u_vec[Knotspan_index])/2
        const double& integration_weight = r_integration_points[point_number].Weight();

        _gke.clear();
        _gfie.clear();
        //stiff_mat_el_lin(rCurrentProcessInfo, point_number, _gke, _gfie, _dL);
        stiff_mat_el_nln(rCurrentProcessInfo, point_number, _gke, _gfie, _dL);
        //stiff_mat_el_geo(rCurrentProcessInfo, point_number, _gke,  _dL);
        
        float mult = (integration_weight * _dL);
 

        if (ComputeLeftHandSide == true)
        {
            //KRATOS_WATCH(_gke);
            //KRATOS_WATCH(rLeftHandSideMatrix);
            //KRATOS_WATCH(r_integration_points[point_number].Coordinates()[0]);
            noalias(rLeftHandSideMatrix) += mult * _gke;
            //KRATOS_WATCH(rLeftHandSideMatrix);
            //std::ostringstream filename_stream;
            //filename_stream << "C:/Users/Max Friedrichs/Documents/Git/Kratos/owncode/BeamV1/BeamV1/stiffeval/" << std::fixed << std::setprecision(4) << this->Id() << ".csv";
            //std::string filename = filename_stream.str();
            //ExportMatrixToCSV(rLeftHandSideMatrix, filename);
        }

        if (ComputeRightHandSide == true)
        {
            //KRATOS_WATCH(_gfie);
            //KRATOS_WATCH(rRightHandSideVector);
            noalias(rRightHandSideVector) += mult * _gfie;
            //KRATOS_WATCH(this->Id());
            //KRATOS_WATCH(rRightHandSideVector);
        }
  
    }

    void IsogeometricBeamElement::GetValuesVector(Vector& rValues, int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const auto& rot = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, Step);
            rValues[index] = disp[0];
            rValues[index + IndexType(1)] = disp[1];
            rValues[index + IndexType(2)] = disp[2];
            rValues[index + IndexType(3)] = rot;
        }
    }

    void IsogeometricBeamElement::GetFirstDerivativesVector(Vector& rValues, int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& vel =r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
            const auto& ang_vel = r_geometry[i].FastGetSolutionStepValue(ANGULAR_VELOCITY_X, Step);
            
            rValues[index] = vel[0];
            rValues[index + IndexType(1)] = vel[1];
            rValues[index + IndexType(2)] = vel[2];
            rValues[index + IndexType(3)] = ang_vel;
        }
    }

    void IsogeometricBeamElement::GetSecondDerivativesVector(Vector& rValues, int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& acc = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const auto& ang_acc = r_geometry[i].FastGetSolutionStepValue(ANGULAR_ACCELERATION_X, Step); 

            rValues[index] = acc[0];
            rValues[index + IndexType(1)] = acc[1];
            rValues[index + IndexType(2)] = acc[2];
            rValues[index + IndexType(3)] = ang_acc;
        }
    }

    ///@}
    ///@name Output functions
    ///@{

    void IsogeometricBeamElement::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::CalculateOnIntegrationPoints");
        const auto& integration_points = GetGeometry().IntegrationPoints();

        Vector3d m_;
        Vector3d f_;

        if (rOutput.size() != integration_points.size()) {
            rOutput.resize(GetGeometry().IntegrationPointsNumber());
        }

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            if (rVariable == FORCE || rVariable == MOMENT)
            {
                stress_res_lin(rCurrentProcessInfo, point_number, f_, m_);
                //stress_res_nln(rCurrentProcessInfo, point_number, f_, m_ );

                if (rVariable == FORCE)
                {
                    rOutput[point_number][0] = f_[0];
                    rOutput[point_number][1] = f_[1];
                    rOutput[point_number][2] = f_[2];
                }
                if (rVariable == MOMENT)
                {
                    rOutput[point_number][0] = m_[0];
                    rOutput[point_number][1] = m_[1];
                    rOutput[point_number][2] = m_[2];
                }
            }
        }
    }

    void IsogeometricBeamElement::CalculateOnIntegrationPoints(const  Variable<Vector>& rVariable, std::vector <Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::CalculateOnIntegrationPoints");
        const auto& integration_points = GetGeometry().IntegrationPoints();

        if (rOutput.size() != integration_points.size()) {
            rOutput.resize(GetGeometry().IntegrationPointsNumber());
        }

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            if (rVariable == LOCAL_CS_N)
            {
                rOutput[point_number] = N;
            }
            if (rVariable == LOCAL_CS_n)
            {
                rOutput[point_number] = n;
            }
            if (rVariable == LOCAL_CS_V)
            {
                rOutput[point_number] = V;
            }
            if (rVariable == LOCAL_CS_v)
            {
                rOutput[point_number] = v;
            }
        }
    }

    void IsogeometricBeamElement::CalculateOnIntegrationPoints(const  Variable<double>& rVariable, std::vector <double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::CalculateOnIntegrationPoints");
        const auto& integration_points = GetGeometry().IntegrationPoints();

        if (rOutput.size() != integration_points.size()) {
            rOutput.resize(GetGeometry().IntegrationPointsNumber());
        }

        for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
            if (rVariable == INTEGRATION_COORDINATES_X)
            {
                rOutput[point_number] = integration_points[point_number].Coordinates()[0];
            }
        }
    }

    ///@}
    ///@name Info
    ///@{

    /// Check provided parameters
    int IsogeometricBeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_WATCH("IsogeometricBeamElement::Check");
        KRATOS_TRY
            const double numerical_limit = std::numeric_limits<double>::epsilon();

        KRATOS_ERROR_IF((GetGeometry().WorkingSpaceDimension() != 3) || (GetGeometry().size() != 2))
            << "The beam element works only in 3D and with 2 noded elements" << std::endl;

        // verify that the variables are correctly initialized
        KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0) << "DISPLACEMENT has Key zero! Check if the application is "
            "registered properly." << std::endl;
        KRATOS_ERROR_IF(CROSS_AREA.Key() == 0) << "CROSS_AREA has Key zero! Check if the application is "
            "registered properly." << std::endl;

        // verify that the dofs exist
        for (IndexType i = 0; i < GetGeometry().size(); ++i) {
            if (GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
            if (GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
                KRATOS_ERROR
                    << "missing one of the dofs for the variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
        }

        KRATOS_ERROR_IF(!GetProperties().Has(CROSS_AREA) ||
            GetProperties()[CROSS_AREA] <= numerical_limit)
            << "Please provide a reasonable value for \"CROSS_AREA\" for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(YOUNG_MODULUS) ||
            GetProperties()[YOUNG_MODULUS] <= numerical_limit)
            << "Please provide a reasonable value for \"YOUNG_MODULUS\" for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(DENSITY) ||
            GetProperties()[DENSITY] <= numerical_limit)
            << "Please provide a reasonable value for \"DENSITY\" for element #"
            << Id() << ". Provided density: " << GetProperties()[DENSITY] << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #" << Id() << std::endl;

        return 0;

        KRATOS_CATCH("Beam Element check.")
    }

    void IsogeometricBeamElement::set_Memory()
    {
        // definition of problem size
        Dof_Node = 4;
        N_Dof = this->GetGeometry().size() * 4;
        _n_Dof = this->GetGeometry().size() * 4;
        _gke.resize(_n_Dof, _n_Dof);
        _gfie.resize(_n_Dof);
        //resizing
        S_gke.resize(_n_Dof, _n_Dof);
        S_fie.resize(_n_Dof);
        S_eps_var.resize(_n_Dof);
        S_curv_n_var.resize(_n_Dof);
        S_curv_v_var.resize(_n_Dof);
        S_torsion_var.resize(_n_Dof);
        S_torsion_n_var.resize(_n_Dof);
        S_torsion_v_var.resize(_n_Dof);
        S_shear_n_var.resize(_n_Dof);
        S_shear_v_var.resize(_n_Dof);
        S_eps_var_var.resize(_n_Dof, _n_Dof);
        S_curv_n_var_var.resize(_n_Dof, _n_Dof);
        S_curv_v_var_var.resize(_n_Dof, _n_Dof);
        S_torsion_var_var.resize(_n_Dof, _n_Dof);
        S_torsion_n_var_var.resize(_n_Dof, _n_Dof);
        S_torsion_v_var_var.resize(_n_Dof, _n_Dof);
        S_kem.resize(_n_Dof, _n_Dof);
        S_keb_n.resize(_n_Dof, _n_Dof);
        S_keb_v.resize(_n_Dof, _n_Dof);
        S_ket_n.resize(_n_Dof, _n_Dof);
        S_ket_v.resize(_n_Dof, _n_Dof);
        S_me.resize(_n_Dof, _n_Dof);
        S_fiem.resize(_n_Dof);
        S_fieb_n.resize(_n_Dof);
        S_fieb_v.resize(_n_Dof);
        S_fiet_n.resize(_n_Dof);
        S_fiet_v.resize(_n_Dof);
        //first variation of mapping matrices
        S_mat_rod_var.resize(_n_Dof * 3, 3);
        S_mat_lam_var.resize(_n_Dof * 3, 3);
        S_mat_lam_var_Rod_Lam.resize(_n_Dof * 3, 3);//new
        S_mat_rod_lam_var_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rodlamRodLam_var.resize(_n_Dof * 3, 3);
        S_mat_rod_der_var.resize(_n_Dof * 3, 3);
        S_mat_lam_der_var.resize(_n_Dof * 3, 3);
        S_mat_lam_der_var_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_lam_var_Rod_der_Lam.resize(_n_Dof * 3, 3);
        S_mat_lam_var_Rod_Lam_der.resize(_n_Dof * 3, 3);
        S_mat_rod_der_lam_var_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_lam_der_var_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_lam_var_Rod_der_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_lam_var_Rod_Lam_der.resize(_n_Dof * 3, 3);
        S_mat_rod_der_var_lam_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_var_lam_der_Rod_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_der_Lam.resize(_n_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_Lam_der.resize(_n_Dof * 3, 3);
        S_mat_rodlamRodLam_der_var.resize(_n_Dof * 3, 3);
        //second variation of mapping matrices
        S_mat_rod_var_var.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_var_var.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_var_var_Rod_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_rod_lam_var_var_Rod_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_rod_var_lam_var_Rod_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_rod_var_var_lam_Rod_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_rodlamRodLam_var_var.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_rod_der_var_var.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_der_var_var.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_der_var_var_Rod_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_var_var_Rod_der_Lam.resize(_n_Dof * 3, _n_Dof * 3);
        S_mat_lam_var_var_Rod_Lam_der.resize(_n_Dof * 3, _n_Dof * 3);
    }


    void IsogeometricBeamElement::parametric_mapping(Vector& _knot_vector, double  u, double & u_mid  )
    {
        for (size_t i = 0; i < _knot_vector.size() - 1; ++i) {
            if (u >= _knot_vector(i) && u < _knot_vector(i + 1)) {
                u_mid = 0.5 * (_knot_vector(i) + _knot_vector(i + 1));
                return;
            }
        }
        u_mid = -1.0;
    }

    void IsogeometricBeamElement::comp_T_var(Vector& _t_var, Vector& _deriv, Vector3d& _r1)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);
        _t_var.resize(3 * N_Dof);
        _t_var.clear();

        float r1_dL = norm_2(_r1);
        float r1_dLpow3 = pow(r1_dL, 3);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 ->rot_tan
                int i = r / Dof_Node;     // index for the shape functions
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) += _deriv[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();

        for (int r = 0;r < N_Dof;r++) //in the case
        {
            for (int t = 0;t < 3;t++)
            {
                r1_r1_var(r) += r1_var[t * N_Dof + r] * _r1[t];
            }
        }

        for (int t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
                _t_var(t * N_Dof + r) += r1_var[t * N_Dof + r] / r1_dL - _r1[t] * r1_r1_var[r] / r1_dLpow3;
        }

    }

    void IsogeometricBeamElement::comp_T_var_var(Matrix& _t_var_var, Vector& _deriv, Vector3d& _r1)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        _t_var_var.resize(3 * N_Dof, N_Dof);
        _t_var_var.clear();

        float r1_dL = norm_2(_r1);
        float r1_pow3 = pow(r1_dL, 3);
        float r1_pow5 = pow(r1_dL, 5);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz && xyz < 3)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        Vector r1_r1var;
        r1_r1var.resize(N_Dof);
        r1_r1var.clear();


        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

            if (xyz > 2)
                r1_r1var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    //if(t==xyz)
                    r1_r1var(r) += r1_var(t * N_Dof + r) * _r1[t];
                }
            }
        }

        Matrix r1var_r1var;
        r1var_r1var.resize(N_Dof, N_Dof);
        r1var_r1var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                for (int s = 0;s < N_Dof;s++) //in the case
                {
                    r1var_r1var(r, s) += r1_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                }
            }
        }


        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                int i = r / Dof_Node;     // index for the shape functions
                for (int s = 0;s < N_Dof;s++)
                {
                    int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int j = s / Dof_Node;     // index for the shape functions
                    if (xyz_r > 2 || xyz_s > 2)
                        _t_var_var(r, s) += 0;
                    else
                    {
                        _t_var_var(t * N_Dof + r, s) += 3 * ((r1_r1var[r] * r1_r1var[s]) * _r1[t]) / r1_pow5;
                        _t_var_var(t * N_Dof + r, s) += (-(r1_var[t * N_Dof + r] * r1_r1var[s]) - (r1_var[t * N_Dof + s] * r1_r1var[r])) / r1_pow3;
                        _t_var_var(t * N_Dof + r, s) += -(r1var_r1var(r, s) * _r1[t]) / r1_pow3;
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_T_deriv_var(Vector& _t_deriv_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);
        _t_deriv_var.resize(3 * N_Dof);
        _t_deriv_var.clear();

        float r1r11 = inner_prod(_r1, _r2);
        //cfloat r1r111 = inner_prod(_r1,_r3);
        float r11r11 = inner_prod(_r2, _r2);
        float r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        Vector r11_var;
        r11_var.resize(3 * N_Dof);
        r11_var.clear();

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                    r11_var(t * N_Dof + r) = _deriv2[i];
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(N_Dof);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(N_Dof);
        r1_r11_var.clear();

        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

            if (xyz > 2)
            {
                r1_r1_var(r) = 0;
                r11_r1_var(r) = 0;
                r1_r11_var(r) = 0;
            }
            else
            {
                for (size_t t = 0;t < 3;t++) //in the case
                {
                    //if(t==xyz)
                    r1_r1_var(r) += r1_var[t * N_Dof + r] * _r1[t];
                    r11_r1_var(r) += r1_var[t * N_Dof + r] * _r2[t];
                    r1_r11_var(r) += r11_var[t * N_Dof + r] * _r1[t];
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                if (xyz > 2)
                    _t_deriv_var[t * N_Dof + r] = 0;
                else
                {
                    _t_deriv_var[t * N_Dof + r] += r11_var[t * N_Dof + r] / r1_dL - _r2(t) * r1_r1_var(r) / pow(r1_dL, 3);
                    _t_deriv_var[t * N_Dof + r] += +3 * r1_r1_var(r) * r1r11 * _r1[t] / pow(r1_dL, 5) - 1.0 / pow(r1_dL, 3) * (r1_var(t * N_Dof + r) * r1r11 + (r1_r11_var(r) + r11_r1_var(r)) * _r1[t]);
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_T_deriv_var_var(Matrix& _t_deriv_var_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        _t_deriv_var_var.resize(3 * N_Dof, N_Dof);
        _t_deriv_var_var.clear();

        float r1r11 = inner_prod(_r1, _r2);
        float r11r11 = inner_prod(_r2, _r2);
        float r1_dL = norm_2(_r1);
        float r1_pow3 = pow(r1_dL, 3);
        float r1_pow5 = pow(r1_dL, 5);
        float r1_pow7 = pow(r1_dL, 7);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        Vector r11_var;
        r11_var.resize(3 * N_Dof);
        r11_var.clear();
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) += _deriv[i];
                    r11_var(t * N_Dof + r) += _deriv2[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(N_Dof);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(N_Dof);
        r1_r11_var.clear();

        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z        
            if (xyz > 2)
                r1_r1_var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    r1_r1_var(r) += r1_var(t * N_Dof + r) * _r1[t];
                    r11_r1_var(r) += _r2[t] * r1_var(t * N_Dof + r);
                    r1_r11_var(r) += _r1[t] * r11_var(t * N_Dof + r);
                }
            }
        }

        Matrix r1_var_r1_var;
        r1_var_r1_var.resize(N_Dof, N_Dof);
        r1_var_r1_var.clear();
        Matrix r1_var_r11_var;
        r1_var_r11_var.resize(N_Dof, N_Dof);
        r1_var_r11_var.clear();
        Matrix r11_var_r1_var;
        r11_var_r1_var.resize(N_Dof, N_Dof);
        r11_var_r1_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                for (int s = 0;s < N_Dof;s++) //in the case
                {
                    r1_var_r1_var(r, s) += r1_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                    r1_var_r11_var(r, s) += r1_var[t * N_Dof + r] * r11_var[t * N_Dof + s];
                    r11_var_r1_var(r, s) += r11_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                }
            }
        }

        for (size_t t = 0;t < 3;t++)
        {
            for (int s = 0;s < N_Dof;s++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    if (xyzr > 2 || xyzs > 2)
                        _t_deriv_var_var(r, s) += 0;
                    else
                    {
                        _t_deriv_var_var(t * N_Dof + r, s) += (-(r11_var(t * N_Dof + r) * r1_r1_var(s)) - (r11_var(t * N_Dof + s) * r1_r1_var(r)) - _r2[t] * r1_var_r1_var(r, s)) / r1_pow3;         //BA auskommentiert schneller:-_r2[t]*r1_var_r1_var(r,s)
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(r) * _r2[t] * r1_r1_var(s))) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_var_r1_var(r, s) * r1r11 * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var[r] * (r11_r1_var[s] + r1_r11_var[s]) * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(r) * r1_var(t * N_Dof + s)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(s) * r1_var(t * N_Dof + r)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += -15 * (r1_r1_var(r) * r1_r1_var(s)) * r1r11 * _r1[t] / r1_pow7;
                        _t_deriv_var_var(t * N_Dof + r, s) += (-(r1_var_r11_var(r, s) + r11_var_r1_var(r, s)) * _r1[t] - ((r1_r11_var(r) + r11_r1_var(r)) * r1_var(t * N_Dof + s)) - ((r1_r11_var(s) + r11_r1_var(s)) * r1_var(t * N_Dof + r))) / r1_pow3;        // BA -(r1_var_r11_var(r,s)+r11_var_r1_var(r,s))*_r1[t] mit  + schneller
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(s) * (r11_r1_var(r) + r1_r11_var(r)) * _r1[t])) / r1_pow5;
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_T_deriv2_var(Vector& _t_deriv2_var, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);

        _t_deriv2_var.resize(3 * N_Dof);
        _t_deriv2_var.clear();

        float r1_r1der = inner_prod(_r1, _r2);
        float r1_r1der_der = inner_prod(_r2, _r2) + inner_prod(_r1, _r3);
        float r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        Vector r1der_var;
        r1der_var.resize(3 * N_Dof);
        r1der_var.clear();
        Vector r1derder_var;
        r1derder_var.resize(3 * N_Dof);
        r1derder_var.clear();
        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (int r = 0; r < N_Dof; r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) = _deriv[i];
                    r1der_var(t * N_Dof + r) = _deriv2[i];
                    r1derder_var(t * N_Dof + r) = _deriv3[i];
                }
            }
        }

        Vector r1_r1var;  //(t * t,r)
        r1_r1var.resize(N_Dof);
        r1_r1var.clear();
        Vector r1der_r1var;  //(t,1 * t,r)
        r1der_r1var.resize(N_Dof);
        r1der_r1var.clear();
        Vector r1_r1dervar;  //(t * t,1,r)
        r1_r1dervar.resize(N_Dof);
        r1_r1dervar.clear();
        Vector r1_r1derdervar;  //(t * t,1,1,r)
        r1_r1derdervar.resize(N_Dof);
        r1_r1derdervar.clear();
        Vector r1der_r1dervar;  //(t,1 * t,1,r)
        r1der_r1dervar.resize(N_Dof);
        r1der_r1dervar.clear();
        Vector r1derder_r1var;  //(t,1,1 * t,r)
        r1derder_r1var.resize(N_Dof);
        r1derder_r1var.clear();
        Vector r1_r1der_var;  //(t * t,1),r
        r1_r1der_var.resize(N_Dof);
        r1_r1der_var.clear();
        Vector r1_r1der_dervar;  //(t * t,1),1,r
        r1_r1der_dervar.resize(N_Dof);
        r1_r1der_dervar.clear();

        for (int r = 0; r < N_Dof; r++) //in the case
        {
            int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

            if (xyz > 2)
            {
                r1_r1var(r) = 0;
                r1der_r1var(r) = 0;
                r1_r1dervar(r) = 0;
                r1_r1derdervar[r] = 0;
                r1der_r1dervar[r] = 0;
                r1derder_r1var[r] = 0;
            }
            else
            {
                for (size_t t = 0; t < 3; t++) //in the case
                {
                    //if(t==xyz)
                    r1_r1var(r) += r1_var[t * N_Dof + r] * _r1[t];
                    r1der_r1var(r) += r1_var[t * N_Dof + r] * _r2[t];
                    r1_r1dervar(r) += r1der_var[t * N_Dof + r] * _r1[t];
                    r1_r1derdervar(r) += r1derder_var[t * N_Dof + r] * _r1[t];
                    r1der_r1dervar(r) += r1der_var[t * N_Dof + r] * _r2[t];
                    r1derder_r1var(r) += r1_var[t * N_Dof + r] * _r3[t];
                }
            }
        }
        r1_r1der_var = r1_r1dervar + r1der_r1var;
        r1_r1der_dervar = r1der_r1dervar * 2 + r1derder_r1var + r1_r1derdervar;

        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (int r = 0; r < N_Dof; r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                if (xyz > 2)
                    _t_deriv2_var[t * N_Dof + r] = 0;
                else
                {
                    _t_deriv2_var[t * N_Dof + r] += r1derder_var[t * N_Dof + r] / r1_dL - _r3(t) * r1_r1var(r) / pow(r1_dL, 3);
                    _t_deriv2_var[t * N_Dof + r] += -(r1der_var[t * N_Dof + r] * r1_r1der + _r2[t] * r1_r1der_var[r]) / pow(r1_dL, 3) + 3 * _r2(t) * r1_r1der * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * N_Dof + r] += -(r1der_var[t * N_Dof + r] * r1_r1der + _r2[t] * r1_r1der_var[r] + r1_var[t * N_Dof + r] * r1_r1der_der + _r1[t] * r1_r1der_dervar[r]) / pow(r1_dL, 3) + 3 * (_r2(t) * r1_r1der + _r1[t] * r1_r1der_der) * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * N_Dof + r] += 3 * (r1_var[t * N_Dof + r] * pow(r1_r1der, 2) + _r1[t] * 2 * r1_r1der * r1_r1der_var[r]) / pow(r1_dL, 5) - 15.0 / pow(r1_dL, 7) * (_r1[t] * pow(r1_r1der, 2) * r1_r1var[r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues(Matrix3d& _mat_rod, Vector3d _vec, float _phi)//here was an error
    {
        _mat_rod.clear();
        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  //initialization by 0 
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        for (int i = 0; i < 3; i++) { _mat_rod(i, i) = cos(_phi); }
        _mat_rod += cross_prod_vec_mat(_vec, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, float _phi, float _phi_deriv)
    {
        _mat_rod_der.clear();

        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  //initialization by 0 
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1; }

        for (int i = 0; i < 3; i++) { _mat_rod_der(i, i) = -_phi_deriv * sin(_phi); }
        _mat_rod_der += cross_prod_vec_mat(_vec, _mat_identity) * cos(_phi) * _phi_deriv;
        _mat_rod_der += cross_prod_vec_mat(_vec_deriv, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_var(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, float _phi)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_rod_var.resize(3*N_Dof,3);
        _mat_rod_var.clear();

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


        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int i = r / Dof_Node;     // index for the shape functions
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_var(t * N_Dof + r, u) += -sin(_phi) * _func[i];
                        else
                            _mat_rod_var(t * N_Dof + r, u) += 0;
                    }
                    else
                    {
                        if (xyz > 2)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_var(t * N_Dof + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            }
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_rod_var(t * N_Dof + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];//*_mat_identity(u,u)
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_var_var(Matrix& _mat_rod_var_var, Vector3d _vec, Vector _vec_var, Matrix _vec_var_var, Vector _func, float _phi)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_rod_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_rod_var_var.clear();

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

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();

        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            int i = r / Dof_Node;     // index for the shape functions

            if (xyz > 2)
                phi_var(r) = _func[i];
            else
                phi_var(r) = 0;
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int s = 0;s < N_Dof;s++)
                {
                    int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    for (int r = 0;r < N_Dof;r++)
                    {
                        int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                        if (t == u)
                        {
                            if (xyzr > 2 || xyzs > 2)
                                _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            else _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += 0;
                        }
                        //else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyzr > 2 || xyzs > 2)
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                                else
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s); //*_mat_identity(s,s)
                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv_var(Matrix& _mat_rod_der_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector _func, Vector _deriv, float _phi, float _phi_der)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_rod_der_var.resize(3*N_Dof,3);
        _mat_rod_der_var.clear();

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

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();

        for (int r = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func(r);
        }

        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (int r = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int i = r / Dof_Node;     // index for the shape functions
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_der_var(t * N_Dof + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        else
                            _mat_rod_der_var(t * N_Dof + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_der_var(t * N_Dof + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];//*_mat_identity(u,u)
                            else
                                _mat_rod_der_var(t * N_Dof + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]; //*_mat_identity(u,u)
                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv_var_var(Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, float _phi, float _phi_der)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_rod_der_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_rod_der_var_var.clear();

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


        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (int r = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func[r];
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        float cs;
        cs = cos(_phi);
        float sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                        int i = r / Dof_Node;     // index for the shape functions
                        if (t == u)
                        {
                            _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        //else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]
                                                    ;

                                                //else      
                                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * N_Dof + r, s); //*_mat_identity(u,u)

                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv2(Matrix3d& _mat_rod_derder, Vector3d _vec, Vector3d _vec_deriv, Vector3d _vec_deriv2, float _phi, float _phi_deriv, float _phi_deriv2)
    {
        _mat_rod_derder.clear();

        Matrix3d _mat_identity;
        _mat_identity.clear();  //initialization by 0 
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1; }

        for (int i = 0; i < 3; i++) { _mat_rod_derder(i, i) = -_phi_deriv2 * sin(_phi) - pow(_phi_deriv, 2) * cos(_phi); }
        _mat_rod_derder += cross_prod_vec_mat(_vec, _mat_identity) * (cos(_phi) * _phi_deriv2 - pow(_phi_deriv, 2) * sin(_phi));
        _mat_rod_derder += cross_prod_vec_mat(_vec_deriv, _mat_identity) * 2 * _phi_deriv * cos(_phi);
        _mat_rod_derder += cross_prod_vec_mat(_vec_deriv2, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv2_var(Matrix& _mat_rod_derder_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector3d _vec_derder, Vector _vec_derder_var, Vector _func, Vector _deriv, Vector _deriv2, float _phi, float _phi_der, float _phi_der2)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);
        _mat_rod_derder_var.resize(3 * N_Dof, 3);
        _mat_rod_derder_var.clear();

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

        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (size_t u = 0; u < 3; u++)
            {
                for (int r = 0; r < N_Dof; r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int i = r / Dof_Node;     // index for the shape functions
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_derder_var(t * N_Dof + r, u) += -_deriv2[i] * sin(_phi) - cos(_phi) * _phi_der2 * _func[i] - 2 * _deriv[i] * _phi_der * cos(_phi) + sin(_phi) * pow(_phi_der, 2) * _func[i];
                        else
                            _mat_rod_derder_var(t * N_Dof + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_derder_var(t * N_Dof + r, u) += (_deriv2[i] * cos(_phi) - _phi_der2 * _func[i] * sin(_phi) - 2 * _deriv[i] * _phi_der * sin(_phi) - pow(_phi_der, 2) * _func[i] * cos(_phi)) * permutation[t][k][u] * _vec[k] + 2 * (cos(_phi) * _deriv[i] - sin(_phi) * _func[i] * _phi_der) * permutation[t][k][u] * _vec_der[k] + cos(_phi) * _func[i] * permutation[t][k][u] * _vec_derder[k];//*_mat_identity(u,u)
                            else
                                _mat_rod_derder_var(t * N_Dof + r, u) += (cos(_phi) * _phi_der2 - pow(_phi_der, 2) * sin(_phi)) * permutation[t][k][u] * _vec_var[k * N_Dof + r] + 2 * _phi_der * cos(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_derder_var[k * N_Dof + r]; //*_mat_identity(u,u)
                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_all(Matrix& _mat_rod_var, Matrix& _mat_rod_der_var, Matrix& _mat_rod_var_var, Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector3d _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, float _phi, float _phi_der)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_rod_var.resize(3*N_Dof,3);
        _mat_rod_var.clear();
        //_mat_rod_der_var.resize(3*N_Dof,3);
        _mat_rod_der_var.clear();
        //_mat_rod_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_rod_var_var.clear();
        //_mat_rod_der_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_rod_der_var_var.clear();

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

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (int r = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func[r];
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        float cs;
        cs = cos(_phi);
        float sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int i = r / Dof_Node;     // index for the shape functions
                    if (t == u)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * N_Dof + r, u) += -sin(_phi) * _func[i];
                            _mat_rod_der_var(t * N_Dof + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * N_Dof + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            _mat_rod_der_var(t * N_Dof + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];//*_mat_identity(u,u)
                        }
                        else
                        {
                            _mat_rod_var(t * N_Dof + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];//*_mat_identity(u,u)
                            _mat_rod_der_var(t * N_Dof + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]; //*_mat_identity(u,u)
                        }
                    }
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                        int i = r / Dof_Node;     // index for the shape functions
                        if (t == u)
                        {
                            _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        //else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                                else
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s); //*_mat_identity(s,s)

                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]
                                                    ;

                                                //else      
                                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * N_Dof + r, s); //*_mat_identity(u,u)
                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda(Matrix3d& _mat_lambda, Vector3d  _vec1, Vector3d _vec2)
    {
        _mat_lambda.clear();  //initialization by 0
        Matrix3d _mat_lambda_tmp;
        _mat_lambda_tmp.clear();  //initialization by 0
        double tmp;

        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  //initialization by 0
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }
        Vector3d cross_vec1_vec2;
        cross_prod(cross_vec1_vec2, _vec1, _vec2);
        float l_cross_vec1_vec2 = norm_2(cross_vec1_vec2);
        Vector3d e_hat = cross_vec1_vec2;
        if (l_cross_vec1_vec2 > 0.000000000001) e_hat = e_hat / l_cross_vec1_vec2;

        if ((inner_prod(_vec1, _vec2) + 1) > Tol)
        {
            for (int i = 0; i < 3; i++) { _mat_lambda(i, i) = inner_prod(_vec1, _vec2); }
            _mat_lambda += cross_prod_vec_mat(cross_prod(_vec1, _vec2), _mat_identity);

            Vector3d v1_cr_v2 = cross_prod(_vec1, _vec2);
            _mat_lambda_tmp = outer_prod(v1_cr_v2, v1_cr_v2);
            tmp = 1.0 / (1.0 + inner_prod(_vec1, _vec2));
            _mat_lambda_tmp = _mat_lambda_tmp * tmp;
            _mat_lambda += _mat_lambda_tmp;
        }
        else
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    _mat_lambda(i, j) = (e_hat(i) * e_hat(j)) * (1 - inner_prod(_vec1, _vec2)) + inner_prod(_vec1, _vec2) * _mat_identity(i, j);
                }
            }
            _mat_lambda += cross_prod_vec_mat(cross_vec1_vec2, _mat_identity);
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv(Matrix3d& _mat_lambda_der, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv)
    {
        _mat_lambda_der.clear();  //initialization by 0

        float tmp;

        Matrix3d _mat_identity;
        _mat_identity.clear();  //initialization by 0
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        float T0_T = inner_prod(_vec1, _vec2);
        float T0_T1 = inner_prod(_vec1, _vec2_deriv);
        float T01_T = inner_prod(_vec1_deriv, _vec2);
        Vector3d T0xT = cross_prod(_vec1, _vec2);
        Vector3d T0xT1 = cross_prod(_vec1, _vec2_deriv);
        Vector3d T01xT = cross_prod(_vec1_deriv, _vec2);
        Vector3d T0xT_1 = T0xT1 + T01xT;

        for (int i = 0; i < 3; i++) { _mat_lambda_der(i, i) = T0_T1 + T01_T; }

        Matrix3d A = cross_prod_vec_mat(T0xT_1, _mat_identity);
        _mat_lambda_der += cross_prod_vec_mat(T0xT_1, _mat_identity);
        tmp = -(T0_T1 + T01_T) / pow((1.0 + T0_T), 2);
        _mat_lambda_der += outer_prod(T0xT, T0xT) * tmp;
        tmp = 1.0 / (1.0 + T0_T);
        _mat_lambda_der += outer_prod(T0xT_1, T0xT) * tmp;
        _mat_lambda_der += outer_prod(T0xT, T0xT_1) * tmp;

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv2(Matrix3d& _mat_lambda_derder, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv, Vector3d _vec1_deriv2, Vector3d _vec2_deriv2)
    {
        _mat_lambda_derder.clear();  //initialization by 0

        float tmp;

        Matrix3d _mat_identity;
        _mat_identity.clear();  //initialization by 0
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        float T0_T = inner_prod(_vec1, _vec2);
        float T0_Tder = inner_prod(_vec1, _vec2_deriv);
        float T0der_T = inner_prod(_vec1_deriv, _vec2);
        float T0der_Tder = inner_prod(_vec1_deriv, _vec2_deriv);
        float T0derder_T = inner_prod(_vec1_deriv2, _vec2);
        float T0_Tderder = inner_prod(_vec1, _vec2_deriv2);
        float T0_T_der = T0der_T + T0_Tder;
        float T0_T_derder = T0derder_T + 2 * T0der_Tder + T0_Tderder;
        Vector3d T0xT = cross_prod(_vec1, _vec2);
        Vector3d T0xTder = cross_prod(_vec1, _vec2_deriv);
        Vector3d T0derxT = cross_prod(_vec1_deriv, _vec2);
        Vector3d T0derxTder = cross_prod(_vec1_deriv, _vec2_deriv);
        Vector3d T0derderxT = cross_prod(_vec1_deriv2, _vec2);
        Vector3d T0xTderder = cross_prod(_vec1, _vec2_deriv2);
        Vector3d T0xT_der = T0xTder + T0derxT;
        Vector3d T0xT_derder = T0derderxT + 2 * T0derxTder + T0xTderder;

        for (int i = 0; i < 3; i++) { _mat_lambda_derder(i, i) = T0_T_derder; }

        _mat_lambda_derder += cross_prod_vec_mat(T0xT_derder, _mat_identity);
        tmp = 2 * pow(T0_T_der, 2) / pow((1.0 + T0_T), 3) - T0_T_derder / pow((1.0 + T0_T), 2);
        _mat_lambda_derder += outer_prod(T0xT, T0xT) * tmp;
        tmp = -T0_T_der / pow((1.0 + T0_T), 2) * 2;
        _mat_lambda_derder += outer_prod(T0xT_der, T0xT) * tmp;
        _mat_lambda_derder += outer_prod(T0xT, T0xT_der) * tmp;
        tmp = 1.0 / (1.0 + T0_T);
        _mat_lambda_derder += outer_prod(T0xT_derder, T0xT) * tmp;
        _mat_lambda_derder += outer_prod(T0xT_der, T0xT_der) * 2 * tmp;
        _mat_lambda_derder += outer_prod(T0xT, T0xT_derder) * tmp;

    }

    void IsogeometricBeamElement::comp_mat_lambda_var(Matrix& _mat_lam_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);

        float T0_T = inner_prod(_vec1, _vec2);
        //_mat_lam_var.resize(3*N_Dof,3);
        _mat_lam_var.clear();

        Matrix3d _mat_identity;
        _mat_identity.clear();  //initialization by 0 
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

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

        Vector cross_vec1_vec2_var;
        Vector cross_vec1_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz > 2)
                            cross_vec1_vec2_var[t * N_Dof + r] += 0;
                        else
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                //cint i = r/Dof_Node;     // index for the shape functions

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int i = r / Dof_Node;     // index for the shape functions
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_lam_var(t * N_Dof + r, u) += 0;
                        else
                            _mat_lam_var(t * N_Dof + r, u) += T0_T_var(r);
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_lam_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_var[r + k * N_Dof]; //*_mat_identity(u,u)
                    }
                    _mat_lam_var(t * N_Dof + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (cross_vec1_vec2[t] * cross_vec1_vec2[u]);
                    _mat_lam_var(t * N_Dof + r, u) += +1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_var_var(Matrix& _mat_lam_var_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var, Matrix _vec2_var_var)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
        //_mat_lam_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_lam_var_var.clear();  //initialization by 0

        float T0_T = inner_prod(_vec1, _vec2);

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

        Vector cross_vec1_vec2_var;
        Matrix cross_vec1_vec2_var_var;
        Vector cross_vec1_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2_var_var.resize(3 * N_Dof, N_Dof);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2_var_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                for (int s = 0;s < N_Dof;s++) //in the case
                {
                    int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                    if (xyzr > 2 || xyzs > 2)
                        T0_T_var_var(r, s) += 0;
                    else
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                }
            }
        }



        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz_r > 2)
                            cross_vec1_vec2_var[t * N_Dof + r] += 0;
                        else
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++) //in the case
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                        for (size_t k = 0;k < 3;k++)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                cross_vec1_vec2_var_var(t * N_Dof + r, s) += 0;
                            else
                                cross_vec1_vec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 ->rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 ->rot_tan
                        if (xyzr > 2 || xyzs > 2) _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) = 0;
                        else
                        {
                            if (t == u)
                                _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_var_var(r, s);
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * cross_vec1_vec2_var_var(k * N_Dof + r, s); //*_mat_identity(u,u)
                            }
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * N_Dof + s] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + r]);
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var_var(t * N_Dof + r, s) * cross_vec1_vec2[u] + cross_vec1_vec2_var(t * N_Dof + r) * cross_vec1_vec2_var(u * N_Dof + s) +
                                cross_vec1_vec2_var(t * N_Dof + s) * cross_vec1_vec2_var(u * N_Dof + r) + cross_vec1_vec2[t] * cross_vec1_vec2_var_var(u * N_Dof + s, r));
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv_var(Matrix& _mat_lam_der_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_lam_der_var.resize(3*N_Dof,3);
        _mat_lam_der_var.clear();  //initialization by 0


        float T0_T = inner_prod(_vec1, _vec2);
        float T0_T1 = inner_prod(_vec1, _vec2_der);
        float T01_T = inner_prod(_vec1_der, _vec2);

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

        Vector cross_vec1_vec2_var;
        Vector cross_vec1_vec2_der_var;
        Vector cross_vec1_der_vec2_var;
        Vector cross_vec1_vec2;
        Vector cross_vec1_vec2_der;
        Vector cross_vec1_der_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2_der_var.resize(N_Dof * 3);
        cross_vec1_der_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);
        cross_vec1_vec2_der.resize(3);
        cross_vec1_der_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2_der_var.clear();
        cross_vec1_der_vec2_var.clear();
        cross_vec1_vec2.clear();
        cross_vec1_vec2_der.clear();
        cross_vec1_der_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);
        cross_vec1_vec2_der = cross_prod(_vec1, _vec2_der);
        cross_vec1_der_vec2 = cross_prod(_vec1_der, _vec2);

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++)
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        int xyz = r % Dof_Node;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2_var[t * N_Dof + r] = 0;
                            cross_vec1_vec2_der_var[t * N_Dof + r] = 0;
                            cross_vec1_der_vec2_var[t * N_Dof + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1_vec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1_der_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                    T0_der_T_var(r) += 0;
                }
                else
                {
                    //if (t==xyz)
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                }
            }
        }


        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                    if (xyz > 2)
                        _mat_lam_der_var(t * N_Dof + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_der_var(t * N_Dof + r, u) += T0_T_der_var(r) + T0_der_T_var(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_der_var[r + k * N_Dof] + permutation[t][k][u] * cross_vec1_der_vec2_var[r + k * N_Dof]; //*_mat_identity(s,s)
                            }

                        }
                        _mat_lam_der_var(t * N_Dof + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_var[t * N_Dof + r]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_var[u * N_Dof + r]));
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]));
                        _mat_lam_der_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * ((cross_vec1_vec2_der_var[t * N_Dof + r] + cross_vec1_der_vec2_var[t * N_Dof + r]) * cross_vec1_vec2[u] + (cross_vec1_vec2_var[t * N_Dof + r]) * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]) + (cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * (cross_vec1_vec2_var[u * N_Dof + r]) + cross_vec1_vec2[t] * (cross_vec1_vec2_der_var[u * N_Dof + r] + cross_vec1_der_vec2_var[u * N_Dof + r]));
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv2_var(Matrix& _mat_lam_derder_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector3d _vec1_derder, Vector _vec2_var, Vector3d _vec2_der, Vector3d _vec2_derder, Vector _vec2_der_var, Vector _vec2_derder_var)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
        _mat_lam_derder_var.resize(3 * N_Dof, 3);
        _mat_lam_derder_var.clear();  //initialization by 0

        float T0_T = inner_prod(_vec1, _vec2);
        float T0_T1 = inner_prod(_vec1, _vec2_der);
        float T01_T = inner_prod(_vec1_der, _vec2);
        float T011_T = inner_prod(_vec1_derder, _vec2);
        float T01_T1 = inner_prod(_vec1_der, _vec2_der);
        float T0_T11 = inner_prod(_vec1, _vec2_derder);
        float T0_T_1 = T01_T + T0_T1;
        float T0_T_11 = T011_T + 2 * T01_T1 + T0_T11;

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

        Vector cross_vec1_vec2var;           // T x t,r = (T x t),r
        Vector cross_vec1_vec2dervar;        // T x t,1,r
        Vector cross_vec1der_vec2var;        // T,1 x t,r
        Vector cross_vec1derder_vec2var;     // T,1,1 x t,r
        Vector cross_vec1der_vec2dervar;     // T,1 x t,1,r
        Vector cross_vec1_vec2derdervar;     // T x t,1,1,r
        Vector cross_vec1_vec2;              // T x t
        Vector cross_vec1_vec2der;           // T x t,1
        Vector cross_vec1der_vec2;           // T,1 x t
        Vector cross_vec1_vec2derder;        // T x t,1,1
        Vector cross_vec1derder_vec2;        // T,1,1 x t
        Vector cross_vec1der_vec2der;        // T,1 x t,1
        Vector cross_vec1_vec2_der;          // (T x t),1
        Vector cross_vec1_vec2_derder;       // (T x t),1,1
        Vector cross_vec1_vec2_dervar;       // (T x t),1,r
        Vector cross_vec1_vec2_derdervar;    // (T x t),1,1,r

        cross_vec1_vec2var.resize(N_Dof * 3);
        cross_vec1_vec2dervar.resize(N_Dof * 3);
        cross_vec1der_vec2var.resize(N_Dof * 3);
        cross_vec1derder_vec2var.resize(N_Dof * 3);
        cross_vec1der_vec2dervar.resize(N_Dof * 3);
        cross_vec1_vec2derdervar.resize(N_Dof * 3);
        cross_vec1_vec2_dervar.resize(N_Dof * 3);
        cross_vec1_vec2_derdervar.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);
        cross_vec1_vec2der.resize(3);
        cross_vec1der_vec2.resize(3);
        cross_vec1derder_vec2.resize(3);
        cross_vec1der_vec2der.resize(3);
        cross_vec1_vec2derder.resize(3);
        cross_vec1_vec2_der.resize(3);
        cross_vec1_vec2_derder.resize(3);

        cross_vec1_vec2var.clear();
        cross_vec1_vec2dervar.clear();
        cross_vec1der_vec2var.clear();
        cross_vec1derder_vec2var.clear();
        cross_vec1der_vec2dervar.clear();
        cross_vec1_vec2derdervar.clear();
        cross_vec1_vec2.clear();
        cross_vec1_vec2der.clear();
        cross_vec1der_vec2.clear();
        cross_vec1derder_vec2.clear();
        cross_vec1der_vec2der.clear();
        cross_vec1_vec2derder.clear();
        cross_vec1_vec2_der.clear();
        cross_vec1_vec2_derder.clear();
        cross_vec1_vec2_dervar.clear();
        cross_vec1_vec2_derdervar.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);
        cross_vec1_vec2der = cross_prod(_vec1, _vec2_der);
        cross_vec1der_vec2 = cross_prod(_vec1_der, _vec2);
        cross_vec1derder_vec2 = cross_prod(_vec1_derder, _vec2);
        cross_vec1der_vec2der = cross_prod(_vec1_der, _vec2_der);
        cross_vec1_vec2derder = cross_prod(_vec1, _vec2_derder);
        cross_vec1_vec2_der = cross_vec1_vec2der + cross_vec1der_vec2;
        cross_vec1_vec2_derder = 2 * cross_vec1der_vec2der + cross_vec1derder_vec2 + cross_vec1_vec2derder;

        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (int r = 0; r < N_Dof; r++)
            {
                for (size_t u = 0; u < 3; u++)
                {
                    for (size_t k = 0; k < 3; k++)
                    {
                        int xyz = r % Dof_Node;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2var[t * N_Dof + r] = 0;
                            cross_vec1_vec2dervar[t * N_Dof + r] = 0;
                            cross_vec1der_vec2var[t * N_Dof + r] = 0;
                            cross_vec1derder_vec2var[t * N_Dof + r] = 0;
                            cross_vec1der_vec2dervar[t * N_Dof + r] = 0;
                            cross_vec1_vec2derdervar[t * N_Dof + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1_vec2dervar[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1der_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1derder_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1_derder[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1der_vec2dervar[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1_vec2derdervar[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_derder_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }
        cross_vec1_vec2_dervar = cross_vec1_vec2dervar + cross_vec1der_vec2var;
        cross_vec1_vec2_derdervar = 2 * cross_vec1der_vec2dervar + cross_vec1derder_vec2var + cross_vec1_vec2derdervar;

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_Tder_var;
        T0_Tder_var.resize(N_Dof);
        T0_Tder_var.clear();
        Vector T0der_T_var;
        T0der_T_var.resize(N_Dof);
        T0der_T_var.clear();
        Vector T0der_Tder_var;
        T0der_Tder_var.resize(N_Dof);
        T0der_Tder_var.clear();
        Vector T0derder_T_var;
        T0derder_T_var.resize(N_Dof);
        T0derder_T_var.clear();
        Vector T0_Tderder_var;
        T0_Tderder_var.resize(N_Dof);
        T0_Tderder_var.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int r = 0; r < N_Dof; r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_Tder_var(r) += 0;
                    T0der_T_var(r) += 0;
                    T0der_Tder_var(r) += 0;
                    T0derder_T_var(r) += 0;
                    T0_Tderder_var(r) += 0;
                }
                else
                {
                    //if (t==xyz)
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_Tder_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                    T0der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0derder_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_derder[t];
                    T0der_Tder_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1_der[t];
                    T0_Tderder_var(r) += _vec2_derder_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }
        Vector T0_T_dervar;
        T0_T_dervar.resize(N_Dof);
        T0_T_dervar.clear();
        Vector T0_T_derdervar;
        T0_T_derdervar.resize(N_Dof);
        T0_T_derdervar.clear();

        T0_T_dervar = T0der_T_var + T0_Tder_var;
        T0_T_derdervar = T0derder_T_var + 2 * T0der_Tder_var + T0_Tderder_var;

        float T0_T_plus_1_pow2_inv = 1.0 / pow(1.0 + T0_T, 2);
        float T0_T_plus_1_pow3_inv = 1.0 / pow(1.0 + T0_T, 3);
        float T0_T_plus_1_pow4_inv = 1.0 / pow(1.0 + T0_T, 4);

        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (size_t u = 0; u < 3; u++)
            {
                for (int r = 0; r < N_Dof; r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                    if (xyz > 2)
                        _mat_lam_derder_var(t * N_Dof + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_derder_var(t * N_Dof + r, u) += T0_T_derdervar(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_derder_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_derdervar[r + k * N_Dof]; //*_mat_identity(s,s)
                            }

                        }
                        _mat_lam_derder_var(t * N_Dof + r, u) += (-T0_T_derdervar(r) * T0_T_plus_1_pow2_inv + 2 * (T0_T_var(r) * T0_T_11 + 2 * T0_T_dervar[r] * T0_T_1) * T0_T_plus_1_pow3_inv - 6 * (pow(T0_T_1, 2) * T0_T_var[r]) * T0_T_plus_1_pow4_inv) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_derder_var(t * N_Dof + r, u) += (-T0_T_11 * T0_T_plus_1_pow2_inv + 2 * pow(T0_T_1, 2) * T0_T_plus_1_pow3_inv) * (cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2var[u * N_Dof + r]));
                        _mat_lam_derder_var(t * N_Dof + r, u) += (4 * T0_T_1 * T0_T_var[r] * T0_T_plus_1_pow3_inv - 2 * T0_T_dervar[r] * T0_T_plus_1_pow2_inv) * (cross_vec1_vec2_der[t] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_der[u]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += -2 * T0_T_1 * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_dervar[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2_der[t] * cross_vec1_vec2var[u * N_Dof + r] + cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_dervar[u * N_Dof + r]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += -T0_T_var(r) * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_derder[t] * cross_vec1_vec2[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derder[u]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_derdervar[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2_derder[t] * cross_vec1_vec2var[u * N_Dof + r] + 2 * cross_vec1_vec2_dervar[t * N_Dof + r] * cross_vec1_vec2_der[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_dervar[u * N_Dof + r] + cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2_derder[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derdervar[u * N_Dof + r]);
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv_var_var(Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
        _mat_lam_der_var_var.clear();  //initialization by 0

        float T0_T = inner_prod(_vec1, _vec2);
        float T0_T1 = inner_prod(_vec1, _vec2_der);
        float T01_T = inner_prod(_vec1_der, _vec2);
        float T0_T_1 = T0_T1 + T01_T;

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

        Vector T0xvec2_var;
        Vector T0xvec2_der_var;
        Vector T0xvec2;
        Vector T0xvec2_der;
        Vector T0_derxvec2;
        Vector T0_derxvec2_var;

        T0xvec2_var.resize(N_Dof * 3);
        T0xvec2_der_var.resize(N_Dof * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(N_Dof * 3);

        T0xvec2_var.clear();
        T0xvec2_der_var.clear();
        T0xvec2.clear();
        T0xvec2_der.clear();
        T0_derxvec2.clear();
        T0_derxvec2_var.clear();

        T0xvec2 = cross_prod(_vec1, _vec2);
        T0xvec2_der = cross_prod(_vec1, _vec2_der);
        T0_derxvec2 = cross_prod(_vec1_der, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(N_Dof, N_Dof);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(N_Dof, N_Dof);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                int i = r / Dof_Node;     // index for the shape functions
                for (int s = 0;s < N_Dof;s++) //in the case
                {
                    int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    int j = s / Dof_Node;     // index for the shape functions

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * N_Dof + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            T0xvec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            T0_derxvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++) //in the case
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                                T0xvec2_der_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * N_Dof + r, s);
                                T0_derxvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * N_Dof + r, s);
                            }
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++) //in the case
            {
                for (int s = 0;s < N_Dof;s++)
                {
                    int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan

                    for (int r = 0;r < N_Dof;r++)
                    {
                        int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                        if (t == u)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                _mat_lam_der_var_var(t * N_Dof + r, u * Dof_Node + s) += 0;
                            else
                                _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                        }
                        else
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * N_Dof + r, s) += 0;
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * N_Dof + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * N_Dof + r, s); //*_mat_identity(u,u) 
                            }
                        }
                    }
                }
            }
        }
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                        if (xyzr > 2 || xyzs > 2)
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 0;
                        else
                        {
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + r, s));   // change of r,s _var_var
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + s) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * N_Dof + r, s) + T0_derxvec2_var_var(t * N_Dof + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * N_Dof + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2_var(u * N_Dof + r) + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)) + (T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2_var(u * N_Dof + s) + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * N_Dof + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * N_Dof + r, s) + T0_derxvec2_var_var(u * N_Dof + r, s)));  // change of r,s _var_var
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_all(Matrix& _mat_lambda_var, Matrix& _mat_lam_der_var, Matrix& _mat_lam_var_var, Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);

        //_mat_lambda_var.resize(3*N_Dof,3);
        _mat_lambda_var.clear();  //initialization by 0
        //_mat_lam_der_var.resize(3*N_Dof,3);
        _mat_lam_der_var.clear();  //initialization by 0
        //_mat_lam_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_lam_var_var.clear();  //initialization by 0 
        //_mat_lam_der_var_var.resize(3*N_Dof,3*N_Dof);
        _mat_lam_der_var_var.clear();  //initialization by 0

        float T0_T = inner_prod(_vec1, _vec2);
        float T0_T1 = inner_prod(_vec1, _vec2_der);
        float T01_T = inner_prod(_vec1_der, _vec2);
        float T0_T_1 = T0_T1 + T01_T;

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

        Vector T0xvec2_var;
        Vector T0xvec2_der_var;
        Vector T0xvec2;
        Vector T0xvec2_der;
        Vector T0_derxvec2;
        Vector T0_derxvec2_var;

        T0xvec2_var.resize(N_Dof * 3);
        T0xvec2_der_var.resize(N_Dof * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(N_Dof * 3);

        T0xvec2_var.clear();
        T0xvec2_der_var.clear();
        T0xvec2.clear();
        T0xvec2_der.clear();
        T0_derxvec2.clear();
        T0_derxvec2_var.clear();

        T0xvec2 = cross_prod(_vec1, _vec2);
        T0xvec2_der = cross_prod(_vec1, _vec2_der);
        T0_derxvec2 = cross_prod(_vec1_der, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(N_Dof, N_Dof);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(N_Dof, N_Dof);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                int i = r / Dof_Node;     // index for the shape functions
                for (int s = 0;s < N_Dof;s++) //in the case
                {
                    int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    int j = s / Dof_Node;     // index for the shape functions

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++)
            {
                int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * N_Dof + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            T0xvec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            T0_derxvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++) //in the case
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                                T0xvec2_der_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * N_Dof + r, s);
                                T0_derxvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * N_Dof + r, s);
                            }
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyz = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z

                    if (xyz > 2)
                    {
                        //_mat_lam_der_var(t*N_Dof+r,u)+= 0;
                                      //_mat_lambda_var(t*N_Dof+r,u)+= 0;
                    }
                    else
                    {
                        if (t == u)
                        {
                            _mat_lambda_var(t * N_Dof + r, u) += T0_T_var(r);
                            _mat_lam_der_var(t * N_Dof + r, u) += T0_T_der_var(r) + T0_der_T_var(r);
                        }
                        else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_lambda_var(t * N_Dof + r, u) += permutation[t][k][u] * T0xvec2_var[r + k * N_Dof];
                                _mat_lam_der_var(t * N_Dof + r, u) += permutation[t][k][u] * T0xvec2_der_var[r + k * N_Dof] + permutation[t][k][u] * T0_derxvec2_var[r + k * N_Dof]; //*_mat_identity(u,u)
                            }
                        }
                        //_mat_lambda_var
                        _mat_lambda_var(t * N_Dof + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (T0xvec2[t] * T0xvec2[u]);
                        _mat_lambda_var(t * N_Dof + r, u) += +1.0 / (1.0 + T0_T) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                        //_mat_lam_der_var
                        _mat_lam_der_var(t * N_Dof + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((T0xvec2_var[t * N_Dof + r]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_var[u * N_Dof + r]));
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                        _mat_lam_der_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var[t * N_Dof + r] + T0_derxvec2_var[t * N_Dof + r]) * T0xvec2[u] + (T0xvec2_var[t * N_Dof + r]) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der[t] + T0_derxvec2[t]) * (T0xvec2_var[u * N_Dof + r]) + T0xvec2[t] * (T0xvec2_der_var[u * N_Dof + r] + T0_derxvec2_var[u * N_Dof + r]));
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    int xyzr = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
                    for (int s = 0;s < N_Dof;s++)
                    {
                        int xyzs = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 ->rot_tan
                        if (xyzr > 2 || xyzs > 2)
                        {
                            //_mat_lam_var_var(t*N_Dof+r,u*N_Dof+s)+=0;
                            //_mat_lam_der_var_var(t*N_Dof+r,u*N_Dof+s)+=0;
                        }
                        else
                        {
                            if (t == u)
                            {
                                _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_var_var(r, s);
                                _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                            }
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_var_var(k * N_Dof + r, s); //*_mat_identity(u,u)
                                    _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * N_Dof + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * N_Dof + r, s); //*_mat_identity(u,u) 
                                }
                            }
                            //_mat_lam_var_var
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + s, r));
                            //_mat_lam_der_var_var
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + r, s));   // change of r,s _var_var
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + s) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * N_Dof + r, s) + T0_derxvec2_var_var(t * N_Dof + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * N_Dof + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2_var(u * N_Dof + r) + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)) + (T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2_var(u * N_Dof + s) + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * N_Dof + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * N_Dof + r, s) + T0_derxvec2_var_var(u * N_Dof + r, s)));  // change of r,s _var_var

                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_Phi_ref_prop(float& _Phi, float& _Phi_0_der)
    {

        //_Phi = this->GetProperties()[PHI];
        //_Phi_0_der = this->GetProperties()[PHI_DER];
        //KRATOS_WATCH(_Phi);
        //KRATOS_WATCH(_Phi_0_der);
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        const float _u_act = r_geometry.IntegrationPoints()[0].Coordinates()[0];

        float phi_0;
        float phi_1;
        float diff_phi;
        float u_0;
        float u_1;
        int n_size;

        Matrix cross_section_orientation = this->GetProperties()[CENTER_LINE_ROTATION];
        n_size = cross_section_orientation.size1();

        // search cross section orientation n before and after _u_act

        u_0 = cross_section_orientation(0, 0);
        u_1 = cross_section_orientation(n_size - 1, 0);
        phi_0 = cross_section_orientation(0, 1);
        phi_1 = cross_section_orientation(n_size - 1, 1);

        for (int i = 1; i < n_size; i++)
        {
            if (cross_section_orientation(i, 0) > _u_act)
            {
                u_0 = cross_section_orientation(i - 1, 0);
                phi_0 = cross_section_orientation(i - 1, 1);
                break;
            }
        }

        for (int i = 1; i < n_size; i++)
        {
            if (cross_section_orientation(n_size - i - 1, 0) <= _u_act)
            {
                u_1 = cross_section_orientation(n_size - i, 0);
                phi_1 = cross_section_orientation(n_size - i, 1);
                break;
            }
        }
        float pi;
        pi = 4 * atan(1.0);

        diff_phi = (phi_1 - phi_0);
        if (fabs(phi_1 - phi_0) > pi)
        {
            diff_phi = diff_phi - (diff_phi) / fabs(diff_phi) * 2 * pi;
        }

        _Phi += phi_0 + (_u_act - u_0) / (u_1 - u_0) * diff_phi;

        _Phi_0_der += diff_phi / (u_1 - u_0);
    }

    void IsogeometricBeamElement::comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector3d& _R1, Vector3d& _R2, float& _A_ref, float& _B_ref)
    {
        const auto& r_geometry = GetGeometry();

        _R1.clear();                //Clear 1st derivative of the curve          
        _R2.clear();                //Clear 2nd derivative of the curve   

        Vector3d coords;  //coordinates of the Nodes
        //KRATOS_WATCH(this->Id());
        //Computation of the Basis functions
        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates(); //Question: Is this really giving the reference coordinates?
            _R1(0) += _deriv(i) * coords(0);
            _R1(1) += _deriv(i) * coords(1);
            _R1(2) += _deriv(i) * coords(2);
            //KRATOS_WATCH(_R1);
            _R2(0) += _deriv2(i) * coords(0);
            _R2(1) += _deriv2(i) * coords(1);
            _R2(2) += _deriv2(i) * coords(2);
            //KRATOS_WATCH(_R2);
        }

        _A_ref = norm_2(_R1);  //length of the base vector

        float tmp = inner_prod(_R2, _R2) - pow(inner_prod(_R1, _R2), 2) / pow(_A_ref, 2);

        if (fabs(tmp) > Tol)
        {
            _B_ref = sqrt(tmp);
        }
        else
            _B_ref = 0;

    }

    void IsogeometricBeamElement::comp_Geometry_reference(Vector _deriv, Vector _deriv2, Vector _deriv3, Vector3d& _R1, Vector3d& _R2, Vector3d& _R3, float& _A_ref, float& _B_ref)
    {
        const auto& r_geometry = GetGeometry();

        _R1.clear();                //Clear 1st derivative of the curve          
        _R2.clear();                //Clear 2nd derivative of the curve   
        _R3.clear();                //Clear 2nd derivative of the curve  

        Vector3d coords;  //coordinates of the Nodes
        //KRATOS_WATCH(this->Id());
        //Computation of the Basis functions
        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates(); //Question: Is this really giving the reference coordinates?
            _R1(0) += _deriv(i) * coords(0);
            _R1(1) += _deriv(i) * coords(1);
            _R1(2) += _deriv(i) * coords(2);
            //KRATOS_WATCH(_R1);
            _R2(0) += _deriv2(i) * coords(0);
            _R2(1) += _deriv2(i) * coords(1);
            _R2(2) += _deriv2(i) * coords(2);
            //KRATOS_WATCH(_R2);
            _R3(0) += _deriv3(i) * coords(0);
            _R3(1) += _deriv3(i) * coords(1);
            _R3(2) += _deriv3(i) * coords(2);
            //KRATOS_WATCH(_R2);
        }

        _A_ref = norm_2(_R1);  //length of the base vector

        float tmp = inner_prod(_R2, _R2) - pow(inner_prod(_R1, _R2), 2) / pow(_A_ref, 2);

        if (fabs(tmp) > Tol)
        {
            _B_ref = sqrt(tmp);
        }
        else
            _B_ref = 0;
    }


    void IsogeometricBeamElement::comp_Geometry_reference_cross_section( Vector3d _R1, Vector3d _R2, Vector3d _T0_vec, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _n0, Vector3d& _v0, float& _B_n, float& _B_v, float& _C_12, float& _C_13, float& _Phi, float& _Phi_0_der)
    {

        comp_Phi_ref_prop(_Phi, _Phi_0_der);

        Matrix3d mat_lamb;
        Matrix3d mat_lamb_deriv;
        Matrix3d mat_rod;
        Matrix3d mat_rod_deriv;
        Matrix3d mat_Ax1;
        mat_Ax1.clear();

        float R1_dL = norm_2(_R1);

        Vector3d T_deriv;
        Vector3d T0_deriv;
        T0_deriv.clear();

        T_deriv = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R1;

        Vector3d  _T_vec = _R1 / R1_dL;

        comp_mat_lambda(mat_lamb, _T0_vec, _T_vec);
        comp_mat_lambda_deriv(mat_lamb_deriv, _T0_vec, _T_vec, T0_deriv, T_deriv);
        float alpha;
        alpha = acos(inner_prod(_T0_vec, _T_vec) / norm_2(_T0_vec) / norm_2(_T_vec));
        float alpha_der;
        alpha_der = -1.0 / (sqrt(1 - pow(inner_prod(_T0_vec, _T_vec) / norm_2(_T0_vec) / norm_2(_T_vec), 2))) * (inner_prod(_T0_vec, T_deriv) / norm_2(_T0_vec) / norm_2(_T_vec) - inner_prod(_T0_vec, _T_vec) / norm_2(_T0_vec) / pow(norm_2(_T_vec), 3) * inner_prod(_T_vec, T_deriv));
        Matrix3d mat_test;
        comp_mat_rodrigues(mat_rod, _T_vec, _Phi);
        comp_mat_rodrigues_deriv(mat_rod_deriv, _T_vec, T_deriv, _Phi, _Phi_0_der);
        _n_act.clear();
        _n0 = this->GetProperties()[N_0];//initial normal vector of the beam's cross-section in the undeformed reference configuration

        float T0_L = norm_2(_T0_vec);
        Vector3d T0 = _T0_vec / T0_L;

        //projection in perpendicular area
        _n0 = _n0 - inner_prod(T0, _n0) * T0;

        _v0 = cross_prod(T0, _n0);

        _n0 = _n0 / norm_2(_n0);

        _v0 = _v0 / norm_2(_v0);

        Vector3d n_tmp;
        n_tmp.clear();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                n_tmp[i] += mat_lamb(i, j) * _n0[j];
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                _n_act(i) += mat_rod(i, j) * n_tmp[j];
            }
        }
        _n_act = _n_act / norm_2(_n_act);

        _v_act = cross_prod(_T_vec, _n_act);

        for (int i = 0; i < 3;i++)
        {
            for (int j = 0; j < 3;j++)
            {
                for (int k = 0; k < 3;k++)
                {
                    mat_Ax1(i, j) += mat_rod_deriv(i, k) * mat_lamb(k, j);
                    mat_Ax1(i, j) += mat_rod(i, k) * mat_lamb_deriv(k, j);
                }
            }
        }

        Vector3d A21;
        Vector3d A31;
        A21.clear();
        A31.clear();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A21[i] += mat_Ax1(i, j) * _n0[j];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A31[i] += mat_Ax1(i, j) * _v0[j];
            }
        }

        _B_n = inner_prod(A21, _R1);
        _B_v = inner_prod(A31, _R1);
        _C_12 = inner_prod(A31, _n_act);
        _C_13 = inner_prod(A21, _v_act);

    }

    void IsogeometricBeamElement::comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, float& _a, float& _b)
    {
        const auto& r_geometry = GetGeometry();
        Vector3d displacement;

        _r1.clear();                //Clear 1st derivative of the curve          
        _r2.clear();                //Clear 2nd derivative of the curve         
        _r3.clear();                //Clear 2nd derivative of the curve         

        Vector3d coords;  //coordinates of the Nodes

        // get previous results
        Vector3d tmp_dof;
        tmp_dof.resize(3, false);
        tmp_dof.clear();

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates();
            //displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, rCurrentProcessInfo.GetSolutionStepIndex()); //ERROR was r_geometry[0] before
            coords[0] += tmp_dof[0];
            coords[1] += tmp_dof[1];
            coords[2] += tmp_dof[2];
            _r1(0) += _deriv(i) * coords[0];
            _r1(1) += _deriv(i) * coords[1];
            _r1(2) += _deriv(i) * coords[2];
            _r2(0) += _deriv2(i) * coords[0];
            _r2(1) += _deriv2(i) * coords[1];
            _r2(2) += _deriv2(i) * coords[2];
            _r3(0) += _deriv3(i) * coords[0];
            _r3(1) += _deriv3(i) * coords[1];
            _r3(2) += _deriv3(i) * coords[2];
        }

        _a = norm_2(_r1);  //length of the base vector

        float tmp = inner_prod(_r2, _r2) - pow(inner_prod(_r1, _r2), 2) / pow(_a, 2);

        //bending
        if (fabs(tmp) > Tol)
            _b = sqrt(tmp);
        else
            _b = 0;
    }

    void IsogeometricBeamElement::comp_Geometry_actual(const ProcessInfo& rCurrentProcessInfo,Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, float& _a, float& _b)
    {
        const auto& r_geometry = GetGeometry();
        Vector3d displacement;

        _r1.clear();                //Clear 1st derivative of the curve          
        _r2.clear();                //Clear 2nd derivative of the curve         

        Vector3d coords;  //coordinates of the Nodes

        // get previous results
        Vector3d tmp_dof;
        tmp_dof.resize(3, false);
        tmp_dof.clear();

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates();
            displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, rCurrentProcessInfo.GetSolutionStepIndex()); //ERROR was r_geometry[0] before
            //tmp_dof = coords + displacement;
            coords[0] += displacement[0];
            coords[1] += displacement[1];
            coords[2] += displacement[2];
            _r1(0) += _deriv(i) * coords[0];
            _r1(1) += _deriv(i) * coords[1];
            _r1(2) += _deriv(i) * coords[2];
            _r2(0) += _deriv2(i) * coords[0];
            _r2(1) += _deriv2(i) * coords[1];
            _r2(2) += _deriv2(i) * coords[2];
        }

        _a = norm_2(_r1);  //length of the base vector

        float tmp = inner_prod(_r2, _r2) - pow(inner_prod(_r1, _r2), 2) / pow(_a, 2);

        //bending
        if (fabs(tmp) > Tol)
            _b = sqrt(tmp);
        else
            _b = 0;

    }

    void IsogeometricBeamElement::comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3, float& _a, float& _b)
    {
        const auto& r_geometry = GetGeometry();
        Vector3d displacement;

        _r1.clear();                //Clear 1st derivative of the curve          
        _r2.clear();                //Clear 2nd derivative of the curve         
        _r3.clear();                //Clear 3rd derivative of the curve   


        Vector3d coords;  //coordinates of the Nodes

        // get initial displacements
        Vector tmp_dof_ini;
        tmp_dof_ini.resize(3, false);
        tmp_dof_ini.clear();
        //
        // get previous results
        Vector3d tmp_dof;
        tmp_dof.resize(3, false);
        tmp_dof.clear();

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates();
            displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, 0); //maybe the 0 is wrong? //ERROR was r_geometry[0] before
            //tmp_dof = coords + displacement;
            coords[0] += displacement[0];
            coords[1] += displacement[1];
            coords[2] += displacement[2];
            _r1(0) += _deriv(i) * coords[0];
            _r1(1) += _deriv(i) * coords[1];
            _r1(2) += _deriv(i) * coords[2];
            _r2(0) += _deriv2(i) * coords[0];
            _r2(1) += _deriv2(i) * coords[1];
            _r2(2) += _deriv2(i) * coords[2];
            _r3(0) += _deriv3(i) * coords[0];
            _r3(1) += _deriv3(i) * coords[1];
            _r3(2) += _deriv3(i) * coords[2];
        }

        _a = norm_2(_r1);  //length of the base vector

        float tmp = inner_prod(_r2, _r2) - pow(inner_prod(_r1, _r2), 2) / pow(_a, 2);

        //bending
        if (fabs(tmp) > Tol)
            _b = sqrt(tmp);
        else
            _b = 0;

    }

    void IsogeometricBeamElement::comp_Geometry_initial(Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2, float& _a, float& _b)
    {
        const auto& r_geometry = GetGeometry();
        Vector3d displacement;

        _r1.clear();                //Clear 1st derivative of the curve          
        _r2.clear();                //Clear 2nd derivative of the curve         

        Vector3d coords;  //coordinates of the Nodes

        // get initial displacements
        Vector tmp_dof_ini;
        tmp_dof_ini.resize(3, false);
        tmp_dof_ini.clear();
        //
        // get previous results
        Vector3d tmp_dof;
        tmp_dof.resize(3, false);
        tmp_dof.clear();
        //KRATOS_WATCH(r_geometry.size());
        //KRATOS_WATCH(this->GetGeometry().size());
        for (size_t i = 0;i < r_geometry.size();i++)
        {
            coords = r_geometry[i].Coordinates();
            displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, 0); //ERROR was r_geometry[0] before
            //tmp_dof = coords + displacement;
            coords[0] += displacement[0];
            coords[1] += displacement[1];
            coords[2] += displacement[2];
            _r1(0) += _deriv(i) * coords[0];
            _r1(1) += _deriv(i) * coords[1];
            _r1(2) += _deriv(i) * coords[2];
            _r2(0) += _deriv2(i) * coords[0];
            _r2(1) += _deriv2(i) * coords[1];
            _r2(2) += _deriv2(i) * coords[2];
        }

        _a = norm_2(_r1);  //length of the base vector

        float tmp = inner_prod(_r2, _r2) - pow(inner_prod(_r1, _r2), 2) / pow(_a, 2);

        //bending
        if (fabs(tmp) > Tol)
            _b = sqrt(tmp);
        else
            _b = 0;

    }

    void IsogeometricBeamElement::comp_Geometry_actual_cross_section(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _N0, Vector3d& _V0, float& _b_n, float& _b_v, float& _c_12, float& _c_13, float _phi, float _phi_der, float _Phi, float _Phi_der)
    {

        Vector3d t0_0 = this->GetProperties()[T_0];

        _n_act.clear();
        _v_act.clear();

        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Ax1;
        mat_Ax1.clear();

        float r1_dL = norm_2(_r1);
        float R1_dL = norm_2(_R1);


        Vector3d t_deriv;
        Vector3d T_deriv;
        Vector3d T0_deriv;
        T0_deriv.clear();

        t_deriv = _r2 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r1;
        T_deriv = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R1;

        Vector3d  _t = _r1 / r1_dL;
        Vector3d  _T_vec = _R1 / R1_dL;   // HELMUT, 24. JULI

        comp_mat_lambda(mat_lam, _T_vec, _t);
        comp_mat_lambda_deriv(mat_lam_der, _T_vec, _t, T_deriv, t_deriv);
        comp_mat_lambda(mat_Lam, t0_0, _T_vec);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, _T_vec, T0_deriv, T_deriv);

        comp_mat_rodrigues(mat_rod, _t, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, _t, t_deriv, _phi, _phi_der);
        comp_mat_rodrigues(mat_Rod, _T_vec, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, _T_vec, T_deriv, _Phi, _Phi_der);

        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                }
            }
        }

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

        Vector3d A21;
        Vector3d A31;
        A21.clear();
        A31.clear();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A21[i] += mat_rodlamRodLam_der(i, j) * _N0[j];
                A31[i] += mat_rodlamRodLam_der(i, j) * _V0[j];
                _n_act[i] += mat_rod_lam_Rod_Lam(i, j) * _N0[j];
                _v_act[i] += mat_rod_lam_Rod_Lam(i, j) * _V0[j];
            }
        }

        _b_n = inner_prod(A21, _r1);
        _b_v = inner_prod(A31, _r1);
        _c_12 = inner_prod(A31, _n_act);
        _c_13 = inner_prod(A21, _v_act);

    }

    void IsogeometricBeamElement::comp_dof_lin(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, float _phi, float _phi_der, float _Phi, float _Phi_der)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::comp_dof_lin");

        KRATOS_TRY

        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1);
        Vector3d t0_0 = this->GetProperties()[T_0];

        //updates here have also to be applied in comp_dof_lin(...) below
          //_cur_var_n.resize(N_Dof);
        _cur_var_n.clear();
        //_cur_var_v.resize(N_Dof);
        _cur_var_v.clear();
        //_tor_var_n.resize(N_Dof);
        _tor_var_n.clear();
        //_tor_var_v.resize(N_Dof);
        _tor_var_v.clear();

        Vector3d t_;
        t_.clear();
        Vector3d t_der;
        t_der.clear();
        Vector3d T_;
        T_.clear();
        Vector3d T_der;
        T_der.clear();
        Vector3d T0_der;
        T0_der.clear();
        Vector t_var;
        Vector t_der_var;
        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;

        t_ = _r1 / norm_2(_r1);
        t_der = _r2 / norm_2(_r1) - inner_prod(_r1, _r2) / pow(norm_2(_r1), 3) * _r1;
        T_ = _R1 / norm_2(_R1);
        T_der = _R2 / norm_2(_R1) - inner_prod(_R1, _R2) / pow(norm_2(_R1), 3) * _R1;
        comp_T_var(t_var, _deriv, _r1);
        comp_T_deriv_var(t_der_var, _deriv, _deriv2, _r1, _r2);

        comp_mat_lambda(mat_lam, T_, t_);
        comp_mat_lambda_deriv(mat_lam_der, T_, t_, T_der, t_der);
        comp_mat_lambda_var(S_mat_lam_var, T_, t_, t_var);
        comp_mat_lambda_deriv_var(S_mat_lam_der_var, T_, t_, T_der, t_var, t_der, t_der_var);
        comp_mat_lambda(mat_Lam, t0_0, T_);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, T_, T0_der, T_der);

        comp_mat_rodrigues(mat_rod, t_, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, t_, t_der, _phi, _phi_der);
        comp_mat_rodrigues_var(S_mat_rod_var, t_, t_var, _func, _phi);
        comp_mat_rodrigues_deriv_var(S_mat_rod_der_var, t_, t_var, t_der, t_der_var, _func, _deriv, _phi, _phi_der);
        comp_mat_rodrigues(mat_Rod, T_, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, T_, T_der, _Phi, _Phi_der);

        // computation of A_i,1 ->der
        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();
        Matrix3d mat_RodLam_der;
        mat_RodLam_der.clear();

        //KRATOS_WATCH(mat_Rod_der);
        //KRATOS_WATCH(mat_Lam);

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                }
            }
        }
        mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();
        Matrix3d mat_lamRodLam_der;
        mat_lamRodLam_der.clear();

        //KRATOS_WATCH(mat_lam);
        //KRATOS_WATCH(mat_Rod_der_Lam);

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                }
            }
        }

        mat_lamRodLam_der = mat_lam_Rod_Lam_der + mat_lam_Rod_der_Lam + mat_lam_der_Rod_Lam;

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

        // variation of A_i,1 ->_der_var

       //matrix<cfloat> mat_lam_var_Rod_Lam_der;
       //mat_lam_var_Rod_Lam_der.resize(3*N_Dof,3);
        S_mat_lam_var_Rod_Lam_der.clear();
        //matrix<cfloat> mat_lam_var_Rod_der_Lam;
        //mat_lam_var_Rod_der_Lam.resize(3*N_Dof,3);
        S_mat_lam_var_Rod_der_Lam.clear();
        //matrix<cfloat> mat_lam_var_Rod_Lam;
        //mat_lam_var_Rod_Lam.resize(3*N_Dof,3);
        S_mat_lam_var_Rod_Lam.clear();
        //matrix<cfloat> mat_lam_der_var_Rod_Lam;
        //mat_lam_der_var_Rod_Lam.resize(3*N_Dof,3);
        S_mat_lam_der_var_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        S_mat_lam_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                        S_mat_lam_var_Rod_Lam_der(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam_der(k, u);
                        S_mat_lam_var_Rod_der_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_der_Lam(k, u);
                        S_mat_lam_der_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                    }
                }
            }
        }

        Matrix mat_lamRodLam_der_var;
        mat_lamRodLam_der_var.resize(3 * N_Dof, 3);
        //mat_lamRodLam_der_var.resize(N_Dof, 3);
        mat_lamRodLam_der_var.clear();

        mat_lamRodLam_der_var = S_mat_lam_var_Rod_Lam_der + S_mat_lam_var_Rod_der_Lam + S_mat_lam_der_var_Rod_Lam;

        S_mat_rod_var_lam_Rod_Lam_der.clear();
        S_mat_rod_var_lam_Rod_der_Lam.clear();
        S_mat_rod_var_lam_der_Rod_Lam.clear();
        S_mat_rod_der_var_lam_Rod_Lam.clear();
        S_mat_rod_der_lam_var_Rod_Lam.clear();
        S_mat_rod_lam_der_var_Rod_Lam.clear();
        S_mat_rod_lam_var_Rod_der_Lam.clear();
        S_mat_rod_lam_var_Rod_Lam_der.clear();
        S_mat_rod_lam_var_Rod_Lam.clear();
        S_mat_rod_var_lam_Rod_Lam.clear();

        //KRATOS_WATCH(S_mat_rod_var);
        //KRATOS_WATCH(mat_lam_Rod_der_Lam);


        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        S_mat_rod_var_lam_Rod_Lam_der(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam_der(k, u);
                        S_mat_rod_var_lam_Rod_der_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_der_Lam(k, u);
                        S_mat_rod_var_lam_der_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_der_Rod_Lam(k, u);
                        S_mat_rod_der_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        S_mat_rod_lam_var_Rod_Lam_der(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_der_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_der_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_der_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod_der(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }

        S_mat_rodlamRodLam_der_var.clear();
        S_mat_rodlamRodLam_var.clear();

        S_mat_rodlamRodLam_der_var = S_mat_rod_var_lam_Rod_Lam_der + S_mat_rod_var_lam_Rod_der_Lam + S_mat_rod_var_lam_der_Rod_Lam + S_mat_rod_der_var_lam_Rod_Lam
            + S_mat_rod_lam_var_Rod_Lam_der + S_mat_rod_lam_var_Rod_der_Lam + S_mat_rod_lam_der_var_Rod_Lam + S_mat_rod_der_lam_var_Rod_Lam;
        S_mat_rodlamRodLam_var = S_mat_rod_var_lam_Rod_Lam + S_mat_rod_lam_var_Rod_Lam;


        Vector3d vec_n;
        vec_n.clear();
        Vector3d vec_v;
        vec_v.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int k = 0;k < 3;k++)
            {
                vec_n(t) += mat_rod_lam_Rod_Lam(t, k) * _N0(k);
                vec_v(t) += mat_rod_lam_Rod_Lam(t, k) * _V0(k);
            }
        }

        Vector vec_n_var;
        vec_n_var.resize(3 * N_Dof);
        vec_n_var.clear();
        Vector vec_v_var;
        vec_v_var.resize(3 * N_Dof);
        vec_v_var.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                for (int k = 0;k < 3;k++)
                {
                    vec_n_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * _N0(k);
                    vec_v_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * _V0(k);
                }
            }
        }

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        for (int t = 0;t < 3;t++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    _cur_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _N0(k) * r1_var[t * N_Dof + r];
                    _cur_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _V0(k) * r1_var[t * N_Dof + r];
                    _tor_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * vec_n(t) + mat_rodlamRodLam_der(t, k) * _V0(k) * vec_n_var(t * N_Dof + r);
                    _tor_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * vec_v(t) + mat_rodlamRodLam_der(t, k) * _N0(k) * vec_v_var(t * N_Dof + r);
                }
            }
        }

        KRATOS_CATCH("")
    }

    void IsogeometricBeamElement::comp_dof_lin(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Vector& _shear_var_n, Vector& _shear_var_v, Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _r3, Vector3d& _R3, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, float _phi, float _phi_der, float _phi_der2, float _Phi, float _Phi_der, float _Phi_der2)
    {
        //KRATOS_WATCH("IsogeometricBeamElement::comp_dof_lin");
        KRATOS_TRY

        Vector3d t0_0 = this->GetProperties()[T_0];


        //updates here have also to be applied in comp_dof_lin(...) above
        //_cur_var_n.resize(N_Dof);
        _cur_var_n.clear();
        // _cur_var_v.resize(N_Dof);
        _cur_var_v.clear();
        //_tor_var_n.resize(N_Dof);
        _tor_var_n.clear();
        //_tor_var_v.resize(N_Dof);
        _tor_var_v.clear();
        //_shear_var_n.resize(N_Dof);
        _shear_var_n.clear();
        //_shear_var_v.resize(N_Dof);
        _shear_var_v.clear();

        Vector3d t_;
        t_.clear();
        Vector3d t_der;
        t_der.clear();
        Vector3d t_derder;
        t_derder.clear();
        Vector3d T_;
        T_.clear();
        Vector3d T_der;
        T_der.clear();
        Vector3d T_derder;
        T_derder.clear();
        Vector3d T0_der;
        T0_der.clear();
        Vector3d T0_derder;
        T0_derder.clear();
        Vector t_var;
        Vector t_der_var;
        Vector t_derder_var;

        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_lam_derder;
        //matrix<cfloat> mat_lam_var;
        //matrix<cfloat> mat_lam_der_var;
        Matrix mat_lam_derder_var;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_Lam_derder;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_rod_derder;
        //matrix<cfloat> mat_rod_var;
        //matrix<cfloat> mat_rod_der_var;
        Matrix mat_rod_derder_var;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Rod_derder;

        float r1_dL = norm_2(_r1);
        float R1_dL = norm_2(_R1);
        t_ = _r1 / r1_dL;
        t_der = _r2 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r1;
        T_ = _R1 / R1_dL;
        T_der = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R1;
        t_derder = _r3 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r2 - (inner_prod(_r2, _r2) + inner_prod(_r1, _r3)) / pow(r1_dL, 3) * _r1 + 3 * inner_prod(_r1, _r2) * inner_prod(_r1, _r2) / pow(r1_dL, 5) * _r1 - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r2;
        T_derder = _R3 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R2 - (inner_prod(_R2, _R2) + inner_prod(_R1, _R3)) / pow(R1_dL, 3) * _R1 + 3 * inner_prod(_R1, _R2) * inner_prod(_R1, _R2) / pow(R1_dL, 5) * _R1 - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R2;

        comp_T_var(t_var, _deriv, _r1);
        comp_T_deriv_var(t_der_var, _deriv, _deriv2, _r1, _r2);
        comp_T_deriv2_var(t_derder_var, _deriv, _deriv2, _deriv3, _r1, _r2, _r3);
        //comp_T_var_var(t_var_var, _deriv, _r1);
        //comp_T_deriv_var_var(t_der_var_var, _deriv, _deriv2, _r1, _r2);

        comp_mat_lambda(mat_lam, T_, t_);
        comp_mat_lambda_deriv(mat_lam_der, T_, t_, T_der, t_der);
        comp_mat_lambda_deriv2(mat_lam_derder, T_, t_, T_der, t_der, T_derder, t_derder);
        comp_mat_lambda_var(S_mat_lam_var, T_, t_, t_var);
        comp_mat_lambda_deriv_var(S_mat_lam_der_var, T_, t_, T_der, t_var, t_der, t_der_var);
        comp_mat_lambda_deriv2_var(mat_lam_derder_var, T_, t_, T_der, T_derder, t_var, t_der, t_derder, t_der_var, t_derder_var);
        comp_mat_lambda(mat_Lam, t0_0, T_);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, T_, T0_der, T_der);
        comp_mat_lambda_deriv2(mat_Lam_derder, t0_0, T_, T0_der, T_der, T0_derder, T_derder);

        comp_mat_rodrigues(mat_rod, t_, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, t_, t_der, _phi, _phi_der);
        comp_mat_rodrigues_deriv2(mat_rod_derder, t_, t_der, t_derder, _phi, _phi_der, _phi_der2);
        comp_mat_rodrigues_var(S_mat_rod_var, t_, t_var, _func, _phi);
        comp_mat_rodrigues_deriv_var(S_mat_rod_der_var, t_, t_var, t_der, t_der_var, _func, _deriv, _phi, _phi_der);
        comp_mat_rodrigues_deriv2_var(mat_rod_derder_var, t_, t_var, t_der, t_der_var, t_derder, t_derder_var, _func, _deriv, _deriv2, _phi, _phi_der, _phi_der2);
        comp_mat_rodrigues(mat_Rod, T_, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, T_, T_der, _Phi, _Phi_der);
        comp_mat_rodrigues_deriv2(mat_Rod_derder, T_, T_der, T_derder, _Phi, _Phi_der, _Phi_der2);

        // computation of A_i,1 ->der
        Matrix3d mat_Rod_Lam_derder;
        mat_Rod_Lam_derder.clear();
        Matrix3d mat_Rod_der_Lam_der;
        mat_Rod_der_Lam_der.clear();
        Matrix3d mat_Rod_derder_Lam;
        mat_Rod_derder_Lam.clear();
        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();
        Matrix3d mat_RodLam_der;
        mat_RodLam_der.clear();
        Matrix3d mat_RodLam_derder;
        mat_RodLam_derder.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                    mat_Rod_derder_Lam(t, u) += mat_Rod_derder(t, k) * mat_Lam(k, u);
                    mat_Rod_der_Lam_der(t, u) += mat_Rod_der(t, k) * mat_Lam_der(k, u);
                    mat_Rod_Lam_derder(t, u) += mat_Rod(t, k) * mat_Lam_derder(k, u);
                }
            }
        }
        mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;
        mat_RodLam_derder = mat_Rod_Lam_derder + 2 * mat_Rod_der_Lam_der + mat_Rod_derder_Lam;

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();
        Matrix3d mat_lamRodLam_der;
        mat_lamRodLam_der.clear();
        Matrix3d mat_lamRodLam_derder;
        mat_lamRodLam_derder.clear();
        Matrix3d mat_lam_RodLam_derder;
        mat_lam_RodLam_derder.clear();
        Matrix3d mat_lam_derder_RodLam;
        mat_lam_derder_RodLam.clear();
        Matrix3d mat_lam_der_RodLam_der;
        mat_lam_der_RodLam_der.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_RodLam_derder(t, u) += mat_lam(t, k) * mat_RodLam_derder(k, u);
                    mat_lam_derder_RodLam(t, u) += mat_lam_derder(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_der_RodLam_der(t, u) += mat_lam_der(t, k) * mat_RodLam_der(k, u);
                }
            }
        }

        mat_lamRodLam_der = mat_lam_Rod_Lam_der + mat_lam_Rod_der_Lam + mat_lam_der_Rod_Lam;
        mat_lamRodLam_derder = mat_lam_RodLam_derder + 2 * mat_lam_der_RodLam_der + mat_lam_derder_RodLam;

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();
        Matrix3d mat_rodlamRodLam_derder;
        mat_rodlamRodLam_derder.clear();
        Matrix3d mat_rod_lamRodLam_derder;
        mat_rod_lamRodLam_derder.clear();
        Matrix3d mat_rod_derder_lamRodLam;
        mat_rod_derder_lamRodLam.clear();
        Matrix3d mat_rod_der_lamRodLam_der;
        mat_rod_der_lamRodLam_der.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_derder_lamRodLam(t, u) += mat_rod_derder(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_der_lamRodLam_der(t, u) += mat_rod_der(t, k) * mat_lamRodLam_der(k, u);
                    mat_rod_lamRodLam_derder(t, u) += mat_rod(t, k) * mat_lamRodLam_derder(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;
        mat_rodlamRodLam_derder = mat_rod_lamRodLam_derder + 2 * mat_rod_der_lamRodLam_der + mat_rod_derder_lamRodLam;

        // variation of A_i,1 ->_der_var

        //matrix<cfloat> mat_lam_var_Rod_Lam_der;
        //mat_lam_var_Rod_Lam_der.resize(3 * N_Dof, 3);
        S_mat_lam_var_Rod_Lam_der.clear();
        //matrix<cfloat> mat_lam_var_Rod_der_Lam;
        //mat_lam_var_Rod_der_Lam.resize(3 * N_Dof, 3);
        S_mat_lam_var_Rod_der_Lam.clear();
        //matrix<cfloat> mat_lam_var_Rod_Lam;
        //mat_lam_var_Rod_Lam.resize(3 * N_Dof, 3);
        S_mat_lam_var_Rod_Lam.clear();
        //matrix<cfloat> mat_lam_der_var_Rod_Lam;
        //mat_lam_der_var_Rod_Lam.resize(3 * N_Dof, 3);
        S_mat_lam_der_var_Rod_Lam.clear();
        Matrix mat_lam_derder_var_RodLam;
        mat_lam_derder_var_RodLam.resize(3 * N_Dof, 3);
        mat_lam_derder_var_RodLam.clear();
        Matrix mat_lam_der_var_RodLam_der;
        mat_lam_der_var_RodLam_der.resize(3 * N_Dof, 3);
        mat_lam_der_var_RodLam_der.clear();
        Matrix mat_lam_var_RodLam_derder;
        mat_lam_var_RodLam_derder.resize(3 * N_Dof, 3);
        mat_lam_var_RodLam_derder.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int r = 0; r < N_Dof; r++)
                    {
                        S_mat_lam_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                        S_mat_lam_var_Rod_Lam_der(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam_der(k, u);
                        S_mat_lam_var_Rod_der_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_der_Lam(k, u);
                        S_mat_lam_der_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                        mat_lam_var_RodLam_derder(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_RodLam_derder(k, u);
                        mat_lam_der_var_RodLam_der(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_RodLam_der(k, u);
                        mat_lam_derder_var_RodLam(t * N_Dof + r, u) += mat_lam_derder_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                    }
                }
            }
        }

        Matrix mat_lamRodLam_der_var;
        mat_lamRodLam_der_var.resize(3 * N_Dof, 3);
        mat_lamRodLam_der_var.clear();
        Matrix mat_lamRodLam_derder_var;
        mat_lamRodLam_derder_var.resize(3 * N_Dof, 3);
        mat_lamRodLam_derder_var.clear();

        mat_lamRodLam_der_var = S_mat_lam_var_Rod_Lam_der + S_mat_lam_var_Rod_der_Lam + S_mat_lam_der_var_Rod_Lam;
        mat_lamRodLam_derder_var = mat_lam_var_RodLam_derder + 2 * mat_lam_der_var_RodLam_der + mat_lam_derder_var_RodLam;

        S_mat_rod_var_lam_Rod_Lam_der.clear();
        S_mat_rod_var_lam_Rod_der_Lam.clear();
        S_mat_rod_var_lam_der_Rod_Lam.clear();
        S_mat_rod_der_var_lam_Rod_Lam.clear();
        S_mat_rod_der_lam_var_Rod_Lam.clear();
        S_mat_rod_lam_der_var_Rod_Lam.clear();
        S_mat_rod_lam_var_Rod_der_Lam.clear();
        S_mat_rod_lam_var_Rod_Lam_der.clear();
        S_mat_rod_lam_var_Rod_Lam.clear();
        S_mat_rod_var_lam_Rod_Lam.clear();
        Matrix mat_rod_derder_var_lamRodLam;
        mat_rod_derder_var_lamRodLam.resize(3 * N_Dof, 3);
        mat_rod_derder_var_lamRodLam.clear();
        Matrix mat_rod_derder_lamRodLam_var;
        mat_rod_derder_lamRodLam_var.resize(3 * N_Dof, 3);
        mat_rod_derder_lamRodLam_var.clear();
        Matrix mat_rod_der_var_lamRodLam_der;
        mat_rod_der_var_lamRodLam_der.resize(3 * N_Dof, 3);
        mat_rod_der_var_lamRodLam_der.clear();
        Matrix mat_rod_der_lamRodLam_der_var;
        mat_rod_der_lamRodLam_der_var.resize(3 * N_Dof, 3);
        mat_rod_der_lamRodLam_der_var.clear();
        Matrix mat_rod_var_lamRodLam_derder;
        mat_rod_var_lamRodLam_derder.resize(3 * N_Dof, 3);
        mat_rod_var_lamRodLam_derder.clear();
        Matrix mat_rod_lamRodLam_derder_var;
        mat_rod_lamRodLam_derder_var.resize(3 * N_Dof, 3);
        mat_rod_lamRodLam_derder_var.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (int r = 0; r < N_Dof; r++)
                    {
                        S_mat_rod_var_lam_Rod_Lam_der(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam_der(k, u);
                        S_mat_rod_var_lam_Rod_der_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_der_Lam(k, u);
                        S_mat_rod_var_lam_der_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_der_Rod_Lam(k, u);
                        S_mat_rod_der_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        S_mat_rod_lam_var_Rod_Lam_der(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_der_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_der_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_der_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod_der(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        mat_rod_derder_var_lamRodLam(t * N_Dof + r, u) += mat_rod_derder_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        mat_rod_derder_lamRodLam_var(t * N_Dof + r, u) += mat_rod_derder(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        mat_rod_der_var_lamRodLam_der(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lamRodLam_der(k, u);
                        mat_rod_der_lamRodLam_der_var(t * N_Dof + r, u) += mat_rod_der(t, k) * mat_lamRodLam_der_var(k * N_Dof + r, u);
                        mat_rod_var_lamRodLam_derder(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lamRodLam_derder(k, u);
                        mat_rod_lamRodLam_derder_var(t * N_Dof + r, u) += mat_rod(t, k) * mat_lamRodLam_derder_var(k * N_Dof + r, u);
                    }
                }
            }
        }

        Matrix mat_rodlamRodLam_derder_var;
        mat_rodlamRodLam_derder_var.resize(3 * N_Dof, 3);
        mat_rodlamRodLam_derder_var.clear();
        Matrix mat_rodlamRodLam_der_var;
        mat_rodlamRodLam_der_var.resize(3 * N_Dof, 3);
        mat_rodlamRodLam_der_var.clear();
        Matrix mat_rodlamRodLam_var;
        mat_rodlamRodLam_var.resize(3 * N_Dof, 3);
        mat_rodlamRodLam_var.clear();

        mat_rodlamRodLam_derder_var = mat_rod_var_lamRodLam_derder + mat_rod_lamRodLam_derder_var + 2 * mat_rod_der_var_lamRodLam_der + 2 * mat_rod_der_lamRodLam_der_var + mat_rod_derder_var_lamRodLam + mat_rod_derder_lamRodLam_var;
        mat_rodlamRodLam_der_var = S_mat_rod_var_lam_Rod_Lam_der + S_mat_rod_var_lam_Rod_der_Lam + S_mat_rod_var_lam_der_Rod_Lam + S_mat_rod_der_var_lam_Rod_Lam
            + S_mat_rod_lam_var_Rod_Lam_der + S_mat_rod_lam_var_Rod_der_Lam + S_mat_rod_lam_der_var_Rod_Lam + S_mat_rod_der_lam_var_Rod_Lam;
        mat_rodlamRodLam_var = S_mat_rod_var_lam_Rod_Lam + S_mat_rod_lam_var_Rod_Lam;

        Vector3d vec_n;
        vec_n.clear();
        Vector3d vec_v;
        vec_v.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int k = 0; k < 3; k++)
            {
                vec_n(t) += mat_rod_lam_Rod_Lam(t, k) * _N0(k);
                vec_v(t) += mat_rod_lam_Rod_Lam(t, k) * _V0(k);
            }
        }

        Vector vec_n_var;
        vec_n_var.resize(3 * N_Dof);
        vec_n_var.clear();
        Vector vec_v_var;
        vec_v_var.resize(3 * N_Dof);
        vec_v_var.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int r = 0; r < N_Dof; r++)
            {
                for (int k = 0; k < 3; k++)
                {
                    vec_n_var(t * N_Dof + r) += mat_rodlamRodLam_var(t * N_Dof + r, k) * _N0(k);
                    vec_v_var(t * N_Dof + r) += mat_rodlamRodLam_var(t * N_Dof + r, k) * _V0(k);
                }
            }
        }

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        Vector r2_var;
        r2_var.resize(3 * N_Dof);
        r2_var.clear();
        for (size_t t = 0; t < 3; t++) //in the case
        {
            for (int r = 0; r < N_Dof; r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) = _deriv[i];
                    r2_var(t * N_Dof + r) = _deriv2[i];
                }
            }
        }

        for (int t = 0; t < 3; t++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int r = 0; r < N_Dof; r++)
                {
                    _cur_var_n(r) += mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _N0(k) * r1_var[t * N_Dof + r];
                    _cur_var_v(r) += mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _V0(k) * r1_var[t * N_Dof + r];
                    _tor_var_n(r) += mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * vec_n(t) + mat_rodlamRodLam_der(t, k) * _V0(k) * vec_n_var(t * N_Dof + r);
                    _tor_var_v(r) += mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * vec_v(t) + mat_rodlamRodLam_der(t, k) * _N0(k) * vec_v_var(t * N_Dof + r);
                    _shear_var_n(r) += mat_rodlamRodLam_derder_var(t * N_Dof + r, k) * _N0(k) * _r1[t] + mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * _r2[t] + mat_rodlamRodLam_derder(t, k) * _N0(k) * r1_var[t * N_Dof + r] + mat_rodlamRodLam_der(t, k) * _N0(k) * r2_var[t * N_Dof + r];
                    _shear_var_v(r) += mat_rodlamRodLam_derder_var(t * N_Dof + r, k) * _V0(k) * _r1[t] + mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * _r2[t] + mat_rodlamRodLam_derder(t, k) * _V0(k) * r1_var[t * N_Dof + r] + mat_rodlamRodLam_der(t, k) * _V0(k) * r2_var[t * N_Dof + r];
                }
            }
        }
        KRATOS_CATCH("")
    }

    void IsogeometricBeamElement::comp_dof_nln(Vector& _cur_var_n, Vector& _cur_var_v, Vector& _tor_var_n, Vector& _tor_var_v, Matrix& _cur_var_n_2, Matrix& _cur_var_v_2, Matrix& _tor_var_n_2, Matrix& _tor_var_v_2,  Vector3d& _r1, Vector3d& _R1, Vector3d& _r2, Vector3d& _R2, Vector3d& _N0, Vector3d& _V0, Vector& _func, Vector& _deriv, Vector& _deriv2, float _phi, float _phi_der, float _Phi, float _Phi_der)
    {
        KRATOS_TRY
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
        Vector3d t0_0 = this->GetProperties()[T_0];

        _cur_var_n.clear();
        _cur_var_v.clear();
        _tor_var_n.clear();
        _tor_var_v.clear();
        _cur_var_n_2.clear();
        _cur_var_v_2.clear();
        _tor_var_n_2.clear();
        _tor_var_v_2.clear();

        Vector3d t_;
        t_.clear();
        Vector3d t_der;
        t_der.clear();
        Vector3d T_;
        T_.clear();
        Vector3d T_der;
        T_der.clear();
        Vector3d T0_der;
        T0_der.clear();
        Vector t_var;
        Vector t_der_var;
        Matrix t_var_var;
        Matrix t_der_var_var;
        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;

        t_ = _r1 / norm_2(_r1);
        t_der = _r2 / norm_2(_r1) - inner_prod(_r1, _r2) / pow(norm_2(_r1), 3) * _r1;
        T_ = _R1 / norm_2(_R1);
        T_der = _R2 / norm_2(_R1) - inner_prod(_R1, _R2) / pow(norm_2(_R1), 3) * _R1;
        comp_T_var(t_var, _deriv, _r1);
        comp_T_deriv_var(t_der_var, _deriv, _deriv2, _r1, _r2);
        comp_T_var_var(t_var_var, _deriv, _r1);
        comp_T_deriv_var_var(t_der_var_var, _deriv, _deriv2, _r1, _r2);

        comp_mat_lambda(mat_lam, T_, t_);
        comp_mat_lambda_deriv(mat_lam_der, T_, t_, T_der, t_der);

        comp_mat_lambda(mat_Lam, t0_0, T_);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, T_, T0_der, T_der);
        comp_mat_lambda_all(S_mat_lam_var, S_mat_lam_der_var, S_mat_lam_var_var, S_mat_lam_der_var_var, T_, t_, T_der, t_var, t_der, t_der_var, t_var_var, t_der_var_var);

        comp_mat_rodrigues(mat_rod, t_, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, t_, t_der, _phi, _phi_der);
        comp_mat_rodrigues_var(S_mat_rod_var, t_, t_var, _func, _phi);
        comp_mat_rodrigues_deriv_var(S_mat_rod_der_var, t_, t_var, t_der, t_der_var, _func, _deriv, _phi, _phi_der);
        comp_mat_rodrigues_var_var(S_mat_rod_var_var, t_, t_var, t_var_var, _func, _phi);
        comp_mat_rodrigues_deriv_var_var(S_mat_rod_der_var_var, t_, t_var, t_der, t_der_var, t_var_var, t_der_var_var, _func, _deriv, _phi, _phi_der);
        //comp_mat_rodrigues_all(mat_rod_var,mat_rod_der_var,mat_rod_var_var,mat_rod_der_var_var, t_,t_var,t_der,t_der_var,t_var_var,t_der_var_var,_func,_deriv, _phi, _phi_der);
        comp_mat_rodrigues(mat_Rod, T_, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, T_, T_der, _Phi, _Phi_der);

        // computation of A_i,1 ->der
        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();
        Matrix3d mat_RodLam_der;
        mat_RodLam_der.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                }
            }
        }
        mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();
        Matrix3d mat_lamRodLam_der;
        mat_lamRodLam_der.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                }
            }
        }

        mat_lamRodLam_der = mat_lam_Rod_Lam_der + mat_lam_Rod_der_Lam + mat_lam_der_Rod_Lam;

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

        // variation of A_i,1 ->_der_var
        S_mat_lam_var_Rod_Lam_der.clear();
        S_mat_lam_var_Rod_der_Lam.clear();
        S_mat_lam_var_Rod_Lam.clear();
        S_mat_lam_der_var_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        S_mat_lam_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                        S_mat_lam_var_Rod_Lam_der(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam_der(k, u);
                        S_mat_lam_var_Rod_der_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_der_Lam(k, u);
                        S_mat_lam_der_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                    }
                }
            }
        }

        Matrix mat_lamRodLam_der_var;
        mat_lamRodLam_der_var.resize(3 * N_Dof, 3);
        mat_lamRodLam_der_var.clear();

        mat_lamRodLam_der_var = S_mat_lam_var_Rod_Lam_der + S_mat_lam_var_Rod_der_Lam + S_mat_lam_der_var_Rod_Lam;

        S_mat_rod_var_lam_Rod_Lam_der.clear();
        S_mat_rod_var_lam_Rod_der_Lam.clear();
        S_mat_rod_var_lam_der_Rod_Lam.clear();
        S_mat_rod_der_var_lam_Rod_Lam.clear();
        S_mat_rod_der_lam_var_Rod_Lam.clear();
        S_mat_rod_lam_der_var_Rod_Lam.clear();
        S_mat_rod_lam_var_Rod_der_Lam.clear();
        S_mat_rod_lam_var_Rod_Lam_der.clear();
        S_mat_rod_lam_var_Rod_Lam.clear();
        S_mat_rod_var_lam_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        S_mat_rod_var_lam_Rod_Lam_der(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam_der(k, u);
                        S_mat_rod_var_lam_Rod_der_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_der_Lam(k, u);
                        S_mat_rod_var_lam_der_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_der_Rod_Lam(k, u);
                        S_mat_rod_der_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        S_mat_rod_lam_var_Rod_Lam_der(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_der_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_der_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_der_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod_der(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }

        //matrix<cfloat> mat_rodlamRodLam_der_var;
        //mat_rodlamRodLam_der_var.resize(3*N_Dof,3);
        S_mat_rodlamRodLam_der_var.clear();
        //matrix<cfloat> mat_rodlamRodLam_var;
        //mat_rodlamRodLam_var.resize(3*N_Dof,3);
        S_mat_rodlamRodLam_var.clear();

        S_mat_rodlamRodLam_der_var = S_mat_rod_var_lam_Rod_Lam_der + S_mat_rod_var_lam_Rod_der_Lam + S_mat_rod_var_lam_der_Rod_Lam + S_mat_rod_der_var_lam_Rod_Lam
            + S_mat_rod_lam_var_Rod_Lam_der + S_mat_rod_lam_var_Rod_der_Lam + S_mat_rod_lam_der_var_Rod_Lam + S_mat_rod_der_lam_var_Rod_Lam;
        S_mat_rodlamRodLam_var = S_mat_rod_var_lam_Rod_Lam + S_mat_rod_lam_var_Rod_Lam;

        // 2nd variation of A_i,1 ->_der_var_var

        S_mat_lam_var_var_Rod_Lam.clear();
        S_mat_lam_der_var_var_Rod_Lam.clear();
        S_mat_lam_var_var_Rod_der_Lam.clear();
        S_mat_lam_var_var_Rod_Lam_der.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        for (int s = 0;s < N_Dof;s++)
                        {
                            S_mat_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam(k, u);
                            S_mat_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_der_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam(k, u);
                            S_mat_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_der_Lam(k, u);
                            S_mat_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam_der(k, u);
                        }
                    }
                }
            }
        }

        Matrix mat_rod_der_lam_var_var_Rod_Lam;
        mat_rod_der_lam_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_der_lam_var_var_Rod_Lam.clear();
        Matrix mat_rod_lam_der_var_var_Rod_Lam;
        mat_rod_lam_der_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_lam_der_var_var_Rod_Lam.clear();
        Matrix mat_rod_lam_var_var_Rod_der_Lam;
        mat_rod_lam_var_var_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_lam_var_var_Rod_der_Lam.clear();
        Matrix mat_rod_lam_var_var_Rod_Lam_der;
        mat_rod_lam_var_var_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_lam_var_var_Rod_Lam_der.clear();

        Matrix mat_rod_der_var_lam_var_Rod_Lam;
        mat_rod_der_var_lam_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_der_var_lam_var_Rod_Lam.clear();
        Matrix mat_rod_var_lam_der_var_Rod_Lam;
        mat_rod_var_lam_der_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_lam_der_var_Rod_Lam.clear();
        Matrix mat_rod_var_lam_var_Rod_der_Lam;
        mat_rod_var_lam_var_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_lam_var_Rod_der_Lam.clear();
        Matrix mat_rod_var_lam_var_Rod_Lam_der;
        mat_rod_var_lam_var_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_lam_var_Rod_Lam_der.clear();

        Matrix mat_rod_der_var_var_lam_Rod_Lam;
        mat_rod_der_var_var_lam_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_der_var_var_lam_Rod_Lam.clear();
        Matrix mat_rod_var_var_lam_der_Rod_Lam;
        mat_rod_var_var_lam_der_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_var_lam_der_Rod_Lam.clear();
        Matrix mat_rod_var_var_lam_Rod_der_Lam;
        mat_rod_var_var_lam_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_var_lam_Rod_der_Lam.clear();
        Matrix mat_rod_var_var_lam_Rod_Lam_der;
        mat_rod_var_var_lam_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_var_lam_Rod_Lam_der.clear();
        Matrix mat_rod_lam_var_var_Rod_Lam;
        mat_rod_lam_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_lam_var_var_Rod_Lam.clear();
        Matrix mat_rod_var_lam_var_Rod_Lam;
        mat_rod_var_lam_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_lam_var_Rod_Lam.clear();
        Matrix mat_rod_var_var_lam_Rod_Lam;
        mat_rod_var_var_lam_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
        mat_rod_var_var_lam_Rod_Lam.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    for (int r = 0;r < N_Dof;r++)
                    {
                        for (int s = 0;s < N_Dof;s++)
                        {
                            mat_rod_der_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod_der(t, k) * S_mat_lam_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                            mat_rod_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_der_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                            mat_rod_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_der_Lam(k * N_Dof + r, u * N_Dof + s);
                            mat_rod_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_Lam_der(k * N_Dof + r, u * N_Dof + s);
                            mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_der_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + s, u);
                            mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + s, u);
                            mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + s, u);
                            mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + s, u);
                            mat_rod_der_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_der_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam(k, u);
                            mat_rod_var_var_lam_der_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_der_Rod_Lam(k, u);
                            mat_rod_var_var_lam_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_der_Lam(k, u);
                            mat_rod_var_var_lam_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam_der(k, u);
                            mat_rod_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                            mat_rod_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + s, u);
                            mat_rod_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam(k, u);
                        }
                    }
                }
            }
        }

        Matrix mat_rodlamRodLam_der_var_var;
        mat_rodlamRodLam_der_var_var.resize(3 * N_Dof, 3 * N_Dof);
        mat_rodlamRodLam_der_var_var.clear();

        Vector3d vec_n;
        vec_n.clear();
        Vector3d vec_v;
        vec_v.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int k = 0;k < 3;k++)
            {
                vec_n(t) += mat_rod_lam_Rod_Lam(t, k) * _N0(k);
                vec_v(t) += mat_rod_lam_Rod_Lam(t, k) * _V0(k);
            }
        }

        Vector vec_n_var;
        vec_n_var.resize(3 * N_Dof);
        vec_n_var.clear();
        Vector vec_v_var;
        vec_v_var.resize(3 * N_Dof);
        vec_v_var.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                for (int k = 0;k < 3;k++)
                {
                    vec_n_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * _N0(k);
                    vec_v_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * _V0(k);
                }
            }
        }

        Matrix mat_rodlamRodLam_var_var;
        mat_rodlamRodLam_var_var.resize(3 * N_Dof, 3 * N_Dof);
        mat_rodlamRodLam_var_var.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    for (int s = 0;s < N_Dof;s++)
                    {
                        mat_rodlamRodLam_var_var(t * N_Dof + r, u * N_Dof + s) += mat_rod_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s);
                        mat_rodlamRodLam_der_var_var(t * N_Dof + r, u * N_Dof + s) += mat_rod_der_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s)
                            + mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s)
                            + mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + s, u * N_Dof + r)
                            + mat_rod_der_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_der_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s);
                    }
                }
            }
        }

        Matrix vec_n_var_var;
        vec_n_var_var.resize(3 * N_Dof, N_Dof);
        vec_n_var_var.clear();
        Matrix vec_v_var_var;
        vec_v_var_var.resize(3 * N_Dof, N_Dof);
        vec_v_var_var.clear();

        for (int t = 0;t < 3;t++)
        {
            for (int r = 0;r < N_Dof;r++)
            {
                for (int s = 0;s < N_Dof;s++)
                {
                    for (int k = 0;k < 3;k++)
                    {
                        vec_n_var_var(t * N_Dof + r, s) += mat_rodlamRodLam_var_var(t * N_Dof + r, k * N_Dof + s) * _N0(k);
                        vec_v_var_var(t * N_Dof + r, s) += mat_rodlamRodLam_var_var(t * N_Dof + r, k * N_Dof + s) * _V0(k);
                    }
                }
            }
        }

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) //in the case
        {
            for (int r = 0;r < N_Dof;r++) //in the case
            {
                int xyz = r % Dof_Node;
                int i = r / Dof_Node;
                if (t == xyz)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        for (int t = 0;t < 3;t++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (int r = 0;r < N_Dof;r++)
                {
                    _cur_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _N0(k) * r1_var[t * N_Dof + r];
                    _cur_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * _r1[t] + mat_rodlamRodLam_der(t, k) * _V0(k) * r1_var[t * N_Dof + r];
                    _tor_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * vec_n(t) + mat_rodlamRodLam_der(t, k) * _V0(k) * vec_n_var(t * N_Dof + r);
                    _tor_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * vec_v(t) + mat_rodlamRodLam_der(t, k) * _N0(k) * vec_v_var(t * N_Dof + r);
                    for (int s = 0;s < N_Dof;s++)
                    {
                        _cur_var_n_2(r, s) += 0;//mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s)* _N0(k)* _r1(t) + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * r1_var(t * N_Dof + s) + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * _N0(k) * r1_var(t * N_Dof + r);
                        _cur_var_v_2(r, s) += 0;//mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s)* _V0(k)* _r1(t) + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * r1_var(t * N_Dof + s) + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * _V0(k) * r1_var(t * N_Dof + r);

                        _tor_var_n_2(r, s) += mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s) * _V0(k) * vec_n(t)
                            + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _V0(k) * vec_n_var(t * N_Dof + s)
                            + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * _V0(k) * vec_n_var(t * N_Dof + r)
                            + mat_rodlamRodLam_der(t, k) * _V0(k) * vec_n_var_var(t * N_Dof + r, s);

                        _tor_var_v_2(r, s) += mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s) * _N0(k) * vec_v(t)
                            + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * _N0(k) * vec_v_var(t * N_Dof + s)
                            + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * _N0(k) * vec_v_var(t * N_Dof + r)
                            + mat_rodlamRodLam_der(t, k) * _N0(k) * vec_v_var_var(t * N_Dof + r, s);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    Vector IsogeometricBeamElement::comp_epsilon_dof(Vector3d& _r1, Vector& _shape_func_deriv)
    {
        //const unsigned int N_Dof = 4;//this->GetGeometry().PointsNumber()* (this->GetGeometry().WorkingSpaceDimension() + 1); changed could lead to an error
        Vector epsilon_var; //variations of the axial strain
        epsilon_var.resize(N_Dof);
        float r1_L2 = norm_2(_r1);
        Vector3d r1 = _r1;
        for (int r = 0;r < N_Dof;r++)
        {
            int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            int i = r / Dof_Node;     // index for the shape functions
            if (xyz_r > 2)
                epsilon_var[r] = 0;
            else
                epsilon_var[r] = r1[xyz_r] * _shape_func_deriv[i];
        }
        //KRATOS_WATCH(epsilon_var);
        return epsilon_var;

    }

    Matrix IsogeometricBeamElement::comp_epsilon_dof_2(Vector3d& _r1, Vector& _shape_func_deriv)
    {
        //const unsigned int N_Dof = this->GetGeometry().PointsNumber() * (this->GetGeometry().WorkingSpaceDimension() + 1);
        Matrix epsilon_var_2; //variations of the axial strain
        epsilon_var_2.resize(N_Dof, N_Dof);
        epsilon_var_2.clear();

        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz_r = r % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 -> rot_tan
            int i = r / Dof_Node;     // index for the shape functions
            if (xyz_r > 2)
                for (int s = 0;s < N_Dof;s++)
                    epsilon_var_2(r, s) = 0.0;
            else
            {
                for (int s = 0;s < N_Dof;s++)
                {
                    int xyz_s = s % Dof_Node; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                    int j = s / Dof_Node;     // index for the shape functions
                    if (xyz_s > 2)
                        epsilon_var_2(r, s) = 0;
                    else
                        if (xyz_r == xyz_s)
                            epsilon_var_2(r, s) = _shape_func_deriv[i] * _shape_func_deriv[j];
                        else
                            epsilon_var_2(r, s) = 0;
                }
            }
        }
        return epsilon_var_2;
    }




    void IsogeometricBeamElement::stiff_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie, float& _dL)
    {
        const float _emod = this->GetProperties()[YOUNG_MODULUS];
        const float poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const float _gmod = _emod / (2.0 * (1.0 + poisson_ratio));
        const float _area = this->GetProperties()[CROSS_AREA];
        const float _m_inert_z = this->GetProperties()[I_Z];
        const float _m_inert_y = this->GetProperties()[I_Y];
        const float _mt_iniert = this->GetProperties()[I_T];
        Vector3d t0_0 = this->GetProperties()[T_0];

        const auto& r_geometry = GetGeometry();

        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, integration_point_index, this->GetIntegrationMethod()), 0);

        //KRATOS_WATCH(R_vec);
        //KRATOS_WATCH(dR_vec);
        //KRATOS_WATCH(ddR_vec);

        //declarations
        Vector3d R_1;  //1st derivative of the curve undeformed config
        Vector3d R_2;  //2nd derivative of the curve undeformed config
        float A;                //length of the base vector
        float B;                //curvature of the curve undeformed config

        float B_n;
        float B_v;
        float C_12;
        float C_13;
        float Phi;
        float Phi_der;

        Vector3d r_1;  //1st derivative of the curve deformed config
        Vector3d r_2;  //2nd derivative of the curve deformed config
        float a;                //length of the base vector
        float b;                //curvature of the curve undeformed config
        float b_n;
        float b_v;
        float c_12;
        float c_13;
        float phi;
        float phi_der;

        //Vector3d n;        //principal axis 1 of cross section
        //Vector3d v;        //principal axis 2 of cross section,
        //Vector3d N;        //principal axis 1 of cross section
        //Vector3d V;        //principal axis 2 of cross section,
        Vector3d N0;       //principal axis 1 of cross section at u=0
        Vector3d V0;       //principal axis 2 of cross section at u=0

        // Material and cross section
        float emod_A = _emod * _area;
        float emod_I_n = _emod * _m_inert_y;
        float emod_I_v = _emod * _m_inert_z;
        float gmod_It = _gmod * _mt_iniert;

        //prestresses
        float prestress = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend1 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend2 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_tor = 0; //NOT YET IMPLEMENTED IN KRATOS

        // get previous results
        double tmp_ini_dof;

        Phi = 0;
        Phi_der = 0;
        phi = 0;
        phi_der = 0;

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            tmp_ini_dof = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, rCurrentProcessInfo.GetSolutionStepIndex()); //maybe this is wrong ERROR was r_geometry[integration_point_index] before
            phi += R_vec(i) * tmp_ini_dof;
            phi_der += dR_vec[i] * tmp_ini_dof;
        }

        //compute configurations
        comp_Geometry_reference(dR_vec, ddR_vec, R_1, R_2, A, B);
        comp_Geometry_initial(dR_vec, ddR_vec, r_1, r_2, a, b); 
        
        //phi += Phi;
        //phi_der += Phi_der;

        comp_Geometry_reference_cross_section(R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
        comp_Geometry_actual_cross_section(r_1, R_1, r_2, R_2, n, v, N0, V0, b_n, b_v, c_12, c_13, phi, phi_der, Phi, Phi_der); 

        _dL = A;
        float Apow2 = pow(A, 2);
        float apow2 = pow(a, 2);

        //stresses
        float E11_m = 0.5 * (apow2 - Apow2);   //Green Lagrange formulation (strain)
        float E11_cur_n = (b_n - B_n);
        float E11_cur_v = (b_v - B_v);
        float E12 = (c_12 - C_12);
        float E13 = (c_13 - C_13);

        //float S11_m = prestress * _area + E11_m * emod_A / Apow2;                      //normal force                  //initial displacement missing
        //float S11_n = prestress_bend1 + E11_cur_n * emod_I_v / Apow2;                        //bending moment n
        //float S11_v = prestress_bend2 + E11_cur_v * emod_I_n / Apow2;                        //bending moment v
        //float S12 = 0.5 * (-prestress_tor + E12 * gmod_It / A);                      //0.5 torsional moment
        //float S13 = 0.5 * (prestress_tor + E13 * gmod_It / A);                      //0.5 torsional moment


        float S11_m = prestress * _area /*+ E11_m * emod_A / Apow2*/;                      //normal force                  //initial displacement missing
        float S11_n = prestress_bend1 /*+ E11_cur_n * emod_I_v / Apow2*/;                        //bending moment n
        float S11_v = prestress_bend2 /*+ E11_cur_v * emod_I_n / Apow2*/;                        //bending moment v
        float S12 = 0.5 * (-prestress_tor /*+ E12 * gmod_It / A*/);                      //0.5 torsional moment
        float S13 = 0.5 * (prestress_tor /*+ E13 * gmod_It / A*/);                      //0.5 torsional moment



        // variation of the axial strain 
        S_eps_var = comp_epsilon_dof(r_1, dR_vec);
        comp_dof_lin(S_curv_n_var, S_curv_v_var, S_torsion_n_var, S_torsion_v_var, r_1, R_1, r_2, R_2, N0, V0, R_vec, dR_vec, ddR_vec, phi, phi_der, Phi, Phi_der);

        //get physical quantities
        //axial strain
        S_eps_var = S_eps_var / Apow2;

        //curvature
        S_curv_n_var = S_curv_n_var / Apow2;
        S_curv_v_var = S_curv_v_var / Apow2;

        //torsion
        S_torsion_n_var = S_torsion_n_var / A;
        S_torsion_v_var = S_torsion_v_var / A;
        

        S_kem = outer_prod(S_eps_var, S_eps_var);
        S_keb_n = outer_prod(S_curv_n_var, S_curv_n_var);
        S_keb_v = outer_prod(S_curv_v_var, S_curv_v_var);
        S_ket_n = outer_prod(S_torsion_n_var, S_torsion_n_var);
        S_ket_v = outer_prod(S_torsion_v_var, S_torsion_v_var);

        _gke.clear();
        _gke += emod_A * S_kem;
        _gke += emod_I_v * S_keb_n;
        _gke += emod_I_n * S_keb_v;
        _gke += 0.5 * gmod_It * S_ket_n;
        _gke += 0.5 * gmod_It * S_ket_v;

        _gfie.clear();
        _gfie += -S11_m * S_eps_var;
        _gfie += -S11_n * S_curv_n_var;
        _gfie += -S11_v * S_curv_v_var;
        _gfie += -S12 * S_torsion_n_var;
        _gfie += -S13 * S_torsion_v_var;


}

    void IsogeometricBeamElement::mass_mat_el_lin(const ProcessInfo& rCurrentProcessInfo, Matrix& _me)   ///!!!! nodal thickness and relative denstiy is missing
    {
        const auto& r_geometry = GetGeometry();
        int P_Deg = r_geometry.PolynomialDegree(0); //this will throw an error
        const int ne = (P_Deg + 1);
        Vector func = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        const float area = this->GetProperties()[CROSS_AREA];
        const float dens = this->GetProperties()[DENSITY];

        for (int r = 0;r < ne;r++)
        {
            for (int s = 0;s < ne;s++)
            {
                _me(4 * s, 4 * r) = func(s) * func(r);
                _me(4 * s + 1, 4 * r + 1) = _me(4 * s, 4 * r);
                _me(4 * s + 2, 4 * r + 2) = _me(4 * s, 4 * r);
                _me(4 * s + 3, 4 * r + 3) = _me(4 * s, 4 * r);
            }
        }
        _me = _me * area * dens;  //* dL is considered in the outer loop!!!!! relative density is not considered
    }

    void IsogeometricBeamElement::stiff_mat_el_geo(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, float& _dL)
    {
        const float _emod = this->GetProperties()[YOUNG_MODULUS];
        const float poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const float _gmod = _emod / (2.0 * (1.0 + poisson_ratio));
        const float _area = this->GetProperties()[CROSS_AREA];
        const float _m_inert_z = this->GetProperties()[I_Z];
        const float _m_inert_y = this->GetProperties()[I_Y];
        const float _mt_iniert = this->GetProperties()[I_T];
        Vector3d t0_0 = this->GetProperties()[T_0];

        const auto& r_geometry = GetGeometry();

        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, integration_point_index, this->GetIntegrationMethod()), 0);

        //declarations
        Vector3d R_1;  //1st derivative of the curve undeformed config
        Vector3d R_2;  //2nd derivative of the curve undeformed config
        float A;                //length of the base vector
        float B;                //curvature including metric
        Vector3d r_1;  //1st derivative of the curve deformed config
        Vector3d r_2;  //2nd derivative of the curve deformed config
        float a;                //length of the base vector
        float b;                //curvature including metric

        float B_n;
        float B_v;
        float C_12;
        float C_13;
        float Phi;
        float Phi_der;
        float b_n;
        float b_v;
        float c_12;
        float c_13;
        float phi;
        float phi_der;

        //Vector3d N;        //principal axis 1 of cross section
        //Vector3d V;        //principal axis 2 of cross section
        //Vector3d n;        //principal axis 1 of cross section
        //Vector3d v;        //principal axis 2 of cross section
        Vector3d N0;
        Vector3d V0;

        // Material and cross section
        float emod_A = _emod * _area;
        float emod_I_n = _emod * _m_inert_y;
        float emod_I_v = _emod * _m_inert_z;
        float gmod_It = _gmod * _mt_iniert;

        float prestress = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend1 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend2 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_tor = 0; //NOT YET IMPLEMENTED IN KRATOS

        //compute configurations
        comp_Geometry_reference(dR_vec, ddR_vec, R_1, R_2, A, B);
        comp_Geometry_actual(rCurrentProcessInfo, dR_vec, ddR_vec, r_1, r_2, a, b);
        // get previous results
        Vector tmp_dof;
        tmp_dof.resize(4, false);
        tmp_dof.clear();
        Vector tmp_ini_dof;
        tmp_ini_dof.resize(4, false);
        tmp_ini_dof.clear();

        phi = 0;
        phi_der = 0;
        Phi = 0;
        Phi_der = 0;

        for (size_t i = 0; i < r_geometry.size(); i++)
        {
            tmp_dof = r_geometry[integration_point_index].FastGetSolutionStepValue(DISPLACEMENT, rCurrentProcessInfo.GetSolutionStepIndex());
            phi += R_vec(i) * tmp_dof[3]; //this is wrong
            phi_der += dR_vec[i] * tmp_dof[3];
        }

        comp_Geometry_reference_cross_section(R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
        comp_Geometry_actual_cross_section(r_1, R_1, r_2, R_2, n, v, N0, V0, b_n, b_v, c_12, c_13, phi, phi_der, Phi, Phi_der);

        _dL = A;
        float Apow2 = pow(A, 2);
        float Apow4 = pow(A, 4);
        float apow2 = pow(a, 2);

        //stresses
        float E11_m = 0.5 * (apow2 - Apow2);   //Green Lagrange formulation (strain)
        float E11_cur_n = (b_n - B_n);
        float E11_cur_v = (b_v - B_v);
        float E12 = (c_12 - C_12);
        float E13 = (c_13 - C_13);

        float S11_m = prestress * _area + E11_m * emod_A / Apow2;                      //normal force
        float S11_n = prestress_bend1 + E11_cur_n * emod_I_v / Apow2;                        //bending moment n
        float S11_v = prestress_bend2 + E11_cur_v * emod_I_n / Apow2;                        //bending moment v
        float S12 = 0.5 * (-prestress_tor + E12 * gmod_It / A);                      //0.5 torsional moment
        float S13 = 0.5 * (prestress_tor + E13 * gmod_It / A);                      //0.5 torsional moment

        // 1st variation
        // variation of the axial strain 
        S_eps_var = comp_epsilon_dof(r_1, dR_vec);
        S_eps_var = S_eps_var / Apow2;

        // 2nd variation
        // variation of the axial strain 
        S_eps_var_var = comp_epsilon_dof_2(r_1, dR_vec);
        S_eps_var_var = S_eps_var_var / Apow2;

        comp_dof_nln(S_curv_n_var, S_curv_v_var, S_torsion_n_var, S_torsion_v_var, S_curv_n_var_var, S_curv_v_var_var, S_torsion_n_var_var, S_torsion_v_var_var, r_1, R_1, r_2, R_2, N0, V0, R_vec, dR_vec, ddR_vec, phi, phi_der, Phi, Phi_der);

        S_curv_n_var = S_curv_n_var / Apow2;
        S_curv_v_var = S_curv_v_var / Apow2;
        S_torsion_n_var = S_torsion_n_var / A;
        S_torsion_v_var = S_torsion_v_var / A;
        S_curv_n_var_var = S_curv_n_var_var / Apow2;
        S_curv_v_var_var = S_curv_v_var_var / Apow2;
        S_torsion_n_var_var = S_torsion_n_var_var / A;
        S_torsion_v_var_var = S_torsion_v_var_var / A;

        //stiffness matrix of the membran part
        for (int r = 0; r < N_Dof; r++)
            for (int s = 0; s < N_Dof; s++)
                S_kem(r, s) = S11_m * S_eps_var_var(r, s);

        //stiffness matrix of the bending part
        for (int r = 0; r < N_Dof; r++)
            for (int s = 0; s < N_Dof; s++)
                S_keb_n(r, s) = S_curv_n_var_var(r, s) * S11_n;

        //stiffness matrix of the bending part
        for (int r = 0; r < N_Dof; r++)
            for (int s = 0; s < N_Dof; s++)
                S_keb_v(r, s) = S_curv_v_var_var(r, s) * S11_v;

        //stiffness matrix of the torsion part
        for (int r = 0; r < N_Dof; r++)
            for (int s = 0; s < N_Dof; s++)
                S_ket_n(r, s) = S12 * S_torsion_n_var_var(r, s);

        //stiffness matrix of the torsion part
        for (int r = 0; r < N_Dof; r++)
            for (int s = 0; s < N_Dof; s++)
                S_ket_v(r, s) = S13 * S_torsion_v_var_var(r, s);

        //compute final element stiffness matrix
        _gke.clear();
        _gke += S_kem;
        _gke += S_keb_n;
        _gke += S_keb_v;
        _gke += S_ket_n;
        _gke += S_ket_v;

        //_gfie = -(S11_m*S_eps_var + S11_n * S_curv_n_var + S11_v * S_curv_v_var + S12 * S_torsion_n_var + S13 * S_torsion_v_var);
    }

    void IsogeometricBeamElement::stiff_mat_el_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Matrix& _gke, Vector& _gfie,  float& _dL)
    {
        const float _emod = this->GetProperties()[YOUNG_MODULUS];
        const float poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const float _gmod = _emod / (2.0 * (1.0 + poisson_ratio));
        const float _area = this->GetProperties()[CROSS_AREA];
        const float _m_inert_z = this->GetProperties()[I_Z];
        const float _m_inert_y = this->GetProperties()[I_Y];
        const float _mt_iniert = this->GetProperties()[I_T];
        Vector3d t0_0 = this->GetProperties()[T_0];

        const auto& r_geometry = GetGeometry();

        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, integration_point_index, this->GetIntegrationMethod()), 0);

        //declarations
        Vector3d R_1;  //1st derivative of the curve undeformed config
        Vector3d R_2;  //2nd derivative of the curve undeformed config
        float A;                //length of the base vector
        float B;                //curvature including metric
        Vector3d r_1;  //1st derivative of the curve deformed config
        Vector3d r_2;  //2nd derivative of the curve deformed config
        float a;                //length of the base vector
        float b;                //curvature including metric

        float B_n;
        float B_v;
        float C_12;
        float C_13;
        float Phi;
        float Phi_der;
        float b_n;
        float b_v;
        float c_12;
        float c_13;
        float phi;
        float phi_der;

        //Vector3d N;        //principal axis 1 of cross section
        //Vector3d V;        //principal axis 2 of cross section
        //Vector3d n;        //principal axis 1 of cross section
        //Vector3d v;        //principal axis 2 of cross section
        Vector3d N0;
        Vector3d V0;

        // Material and cross section
        float emod_A = _emod * _area;
        float emod_I_n = _emod * _m_inert_z;
        float emod_I_v = _emod * _m_inert_y;
        float gmod_It = _gmod * _mt_iniert;

        //prestresses
        float prestress = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend1 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_bend2 = 0; //NOT YET IMPLEMENTED IN KRATOS
        float prestress_tor = 0; //NOT YET IMPLEMENTED IN KRATOS
        bool prestress_bend1_auto = 0; //NOT YET IMPLEMENTED IN KRATOS
        bool prestress_bend2_auto = 0; //NOT YET IMPLEMENTED IN KRATOS
        bool prestress_tor_auto = 0; //NOT YET IMPLEMENTED IN KRATOS

        //compute configurations
        comp_Geometry_reference(dR_vec, ddR_vec, R_1, R_2, A, B);
        comp_Geometry_actual(rCurrentProcessInfo,dR_vec, ddR_vec, r_1, r_2, a, b);
        // get previous results
        double tmp_ini_dof;

        phi = 0;
        phi_der = 0;
        Phi = 0;
        Phi_der = 0;

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            tmp_ini_dof = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, rCurrentProcessInfo.GetSolutionStepIndex()); //maybe this is wrong ERROR was r_geometry[integration_point_index] before
            phi += R_vec(i) * tmp_ini_dof;
            phi_der += dR_vec[i] * tmp_ini_dof;
        }

        comp_Geometry_reference_cross_section( R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
        comp_Geometry_actual_cross_section(r_1, R_1, r_2, R_2, n, v, N0, V0, b_n, b_v, c_12, c_13, phi, phi_der, Phi, Phi_der);

        _dL = A;
        float Apow2 = pow(A, 2);
        float Apow4 = pow(A, 4);
        float apow2 = pow(a, 2);

        ////modified reference configuration
        //float lf = Prop_Ptr->get_Act_LoadFactor();
        //
        //if (prestress_bend1_auto)
        //    B_n = B_n * (1.0 - lf);
        //if (prestress_bend2_auto)
        //    B_v = B_v * (1.0 - lf);
        //if (prestress_tor_auto)
        //{
        //    C_12 = C_12 * (1.0 - lf);
        //    C_13 = C_13 * (1.0 - lf);
        //}
        // 
        //stresses
        float E11_m = 0.5 * (apow2 - Apow2);   //Green Lagrange formulation (strain)
        float E11_cur_n = (b_n - B_n);
        float E11_cur_v = (b_v - B_v);
        float E12 = (c_12 - C_12);
        float E13 = (c_13 - C_13);

        float S11_m = prestress * _area + E11_m * emod_A / Apow2;                      //normal force
        float S11_n = prestress_bend1 + E11_cur_n * emod_I_v / Apow2;                        //bending moment n
        float S11_v = prestress_bend2 + E11_cur_v * emod_I_n / Apow2;                        //bending moment v
        float S12 = 0.5 * (-prestress_tor + E12 * gmod_It / A);                      //0.5 torsional moment
        float S13 = 0.5 * (prestress_tor + E13 * gmod_It / A);                      //0.5 torsional moment

        // 1st variation
        // variation of the axial strain 
        S_eps_var = comp_epsilon_dof(r_1, dR_vec);
        S_eps_var = S_eps_var / Apow2;

        // 2nd variation
        // variation of the axial strain 
        S_eps_var_var = comp_epsilon_dof_2(r_1, dR_vec);
        S_eps_var_var = S_eps_var_var / Apow2;

        comp_dof_nln(S_curv_n_var, S_curv_v_var, S_torsion_n_var, S_torsion_v_var, S_curv_n_var_var, S_curv_v_var_var, S_torsion_n_var_var, S_torsion_v_var_var, r_1, R_1, r_2, R_2, N0, V0, R_vec, dR_vec, ddR_vec, phi, phi_der, Phi, Phi_der);

        S_curv_n_var = S_curv_n_var / Apow2;
        S_curv_v_var = S_curv_v_var / Apow2;
        S_torsion_n_var = S_torsion_n_var / A;
        S_torsion_v_var = S_torsion_v_var / A;
        S_curv_n_var_var = S_curv_n_var_var / Apow2;
        S_curv_v_var_var = S_curv_v_var_var / Apow2;
        S_torsion_n_var_var = S_torsion_n_var_var / A;
        S_torsion_v_var_var = S_torsion_v_var_var / A;


        //stiffness matrix of the membran part
        for (int r = 0;r < N_Dof;r++)
            for (int s = 0;s < N_Dof;s++)
                S_kem(r, s) = emod_A * S_eps_var[r] * S_eps_var[s] + S11_m * S_eps_var_var(r, s);

        //stiffness matrix of the bending part
        for (int r = 0;r < N_Dof;r++)
            for (int s = 0;s < N_Dof;s++)
                S_keb_n(r, s) = emod_I_v * S_curv_n_var[r] * S_curv_n_var[s] + S_curv_n_var_var(r, s) * S11_n;

        //stiffness matrix of the bending part
        for (int r = 0;r < N_Dof;r++)
            for (int s = 0;s < N_Dof;s++)
                S_keb_v(r, s) = emod_I_n * S_curv_v_var[r] * S_curv_v_var[s] + S_curv_v_var_var(r, s) * S11_v;

        //stiffness matrix of the torsion part
        for (int r = 0;r < N_Dof;r++)
            for (int s = 0;s < N_Dof;s++)
                S_ket_n(r, s) = 0.5 * gmod_It * S_torsion_n_var[r] * S_torsion_n_var[s] + S12 * S_torsion_n_var_var(r, s);

        //stiffness matrix of the torsion part
        for (int r = 0;r < N_Dof;r++)
            for (int s = 0;s < N_Dof;s++)
                S_ket_v(r, s) = 0.5 * gmod_It * S_torsion_v_var[r] * S_torsion_v_var[s] + S13 * S_torsion_v_var_var(r, s);
        //compute final element stiffness matrix
        _gke.clear();
        _gke += S_kem;
        _gke += S_keb_n;
        _gke += S_keb_v;
        _gke += S_ket_n;
        _gke += S_ket_v;

        _gfie = -(S11_m * S_eps_var + S11_n * S_curv_n_var + S11_v * S_curv_v_var + S12 * S_torsion_n_var + S13 * S_torsion_v_var);
    }
    

    void IsogeometricBeamElement::comp_transverse_shear_force_nln(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d _r3, Vector3d _R3,
        Vector3d& _N0, Vector3d& _V0, float _phi, float _phi_der, float _phi_der2, float _Phi, float _Phi_der, float _Phi_der2,
        float& _shear_force_n, float& _shear_force_v)
    {
        Vector3d t0_0 = this->GetProperties()[T_0];

        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_lam_derder;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_rod_derder;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_Lam_derder;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Rod_derder;

        float r1_dL = norm_2(_r1);
        float R1_dL = norm_2(_R1);

        Vector3d t_deriv;
        Vector3d T_deriv;
        Vector3d T0_deriv;   //dummy
        T0_deriv.clear();
        Vector3d T0_deriv2;   //dummy
        T0_deriv2.clear();
        Vector3d t_deriv2;
        Vector3d T_deriv2;

        t_deriv = _r2 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r1;
        T_deriv = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(r1_dL, 3) * _R1;
        t_deriv2 = _r3 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r2 - (inner_prod(_r2, _r2) + inner_prod(_r1, _r3)) / pow(r1_dL, 3) * _r1 + 3 * inner_prod(_r1, _r2) * inner_prod(_r1, _r2) / pow(r1_dL, 5) * _r1 - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r2;
        T_deriv2 = _R3 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R2 - (inner_prod(_R2, _R2) + inner_prod(_R1, _R3)) / pow(R1_dL, 3) * _R1 + 3 * inner_prod(_R1, _R2) * inner_prod(_R1, _R2) / pow(R1_dL, 5) * _R1 - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R2;

        Vector3d  _t = _r1 / r1_dL;
        Vector3d  _T_vec = _R1 / R1_dL;   // HELMUT, 24. JULI

        comp_mat_lambda(mat_lam, _T_vec, _t);
        comp_mat_lambda_deriv(mat_lam_der, _T_vec, _t, T_deriv, t_deriv);
        comp_mat_lambda_deriv2(mat_lam_derder, _T_vec, _t, T_deriv, t_deriv, T_deriv2, t_deriv2);
        comp_mat_lambda(mat_Lam, t0_0, _T_vec);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, _T_vec, T0_deriv, T_deriv);
        comp_mat_lambda_deriv2(mat_Lam_derder, t0_0, _T_vec, T0_deriv, T_deriv, T0_deriv2, T_deriv2);

        comp_mat_rodrigues(mat_rod, _t, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, _t, t_deriv, _phi, _phi_der);
        comp_mat_rodrigues_deriv2(mat_rod_derder, _t, t_deriv, t_deriv2, _phi, _phi_der, _phi_der2);
        comp_mat_rodrigues(mat_Rod, _T_vec, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, _T_vec, T_deriv, _Phi, _Phi_der);
        comp_mat_rodrigues_deriv2(mat_Rod_derder, _T_vec, T_deriv, T_deriv2, _Phi, _Phi_der, _Phi_der2);

        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();
        Matrix3d mat_Rod_der_Lam_der;
        mat_Rod_der_Lam_der.clear();
        Matrix3d mat_Rod_Lam_derder;
        mat_Rod_Lam_derder.clear();
        Matrix3d mat_Rod_derder_Lam;
        mat_Rod_derder_Lam.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                    mat_Rod_der_Lam_der(t, u) += mat_Rod_der(t, k) * mat_Lam_der(k, u);
                    mat_Rod_derder_Lam(t, u) += mat_Rod_derder(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_derder(t, u) += mat_Rod(t, k) * mat_Lam_derder(k, u);
                }
            }
        }

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam_der;
        mat_lam_der_Rod_Lam_der.clear();
        Matrix3d mat_lam_der_Rod_der_Lam;
        mat_lam_der_Rod_der_Lam.clear();
        Matrix3d mat_lam_derder_Rod_Lam;
        mat_lam_derder_Rod_Lam.clear();
        Matrix3d mat_lam_Rod_der_Lam_der;
        mat_lam_Rod_der_Lam_der.clear();
        Matrix3d mat_lam_Rod_derder_Lam;
        mat_lam_Rod_derder_Lam.clear();
        Matrix3d mat_lam_Rod_Lam_derder;
        mat_lam_Rod_Lam_derder.clear();

        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_der_Rod_Lam_der(t, u) += mat_lam_der(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_der_Rod_der_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_derder_Rod_Lam(t, u) += mat_lam_derder(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_der_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_der_Lam_der(k, u);
                    mat_lam_Rod_derder_Lam(t, u) += mat_lam(t, k) * mat_Rod_derder_Lam(k, u);
                    mat_lam_Rod_Lam_derder(t, u) += mat_lam(t, k) * mat_Rod_Lam_derder(k, u);
                }
            }
        }

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam_der;
        mat_rod_der_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_der_lam_Rod_der_Lam;
        mat_rod_der_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_der_lam_der_Rod_Lam;
        mat_rod_der_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_derder_lam_Rod_Lam;
        mat_rod_derder_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam_der;
        mat_rod_lam_der_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_der_Rod_der_Lam;
        mat_rod_lam_der_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_derder_Rod_Lam;
        mat_rod_lam_derder_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam_der;
        mat_rod_lam_Rod_der_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_derder_Lam;
        mat_rod_lam_Rod_derder_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam_derder;
        mat_rod_lam_Rod_Lam_derder.clear();
        for (int t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);

                    mat_rod_der_lam_Rod_Lam_der(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_der_lam_Rod_der_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_der_lam_der_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_derder_lam_Rod_Lam(t, u) += mat_rod_derder(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam_der(k, u);
                    mat_rod_lam_der_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_der_Lam(k, u);
                    mat_rod_lam_derder_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_derder_Rod_Lam(k, u);
                    mat_rod_lam_Rod_der_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam_der(k, u);
                    mat_rod_lam_Rod_derder_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_derder_Lam(k, u);
                    mat_rod_lam_Rod_Lam_derder(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_derder(k, u);
                }
            }
        }
        Matrix3d mat_RodLam_der;
        mat_RodLam_der.clear();
        mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;
        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;
        Matrix3d mat_RodLam_derder;
        mat_RodLam_derder.clear();
        mat_RodLam_derder = mat_Rod_Lam_derder + 2 * mat_Rod_der_Lam_der + mat_Rod_derder_Lam;
        Matrix3d mat_rodlamRodLam_derder;
        mat_rodlamRodLam_derder.clear();
        mat_rodlamRodLam_derder =
            2 * mat_rod_der_lam_Rod_Lam_der + 2 * mat_rod_der_lam_Rod_der_Lam + 2 * mat_rod_der_lam_der_Rod_Lam + mat_rod_derder_lam_Rod_Lam
            + 2 * mat_rod_lam_der_Rod_Lam_der + 2 * mat_rod_lam_der_Rod_der_Lam + mat_rod_lam_derder_Rod_Lam
            + 2 * mat_rod_lam_Rod_der_Lam_der + mat_rod_lam_Rod_derder_Lam
            + mat_rod_lam_Rod_Lam_derder;
        Vector3d A21;
        Vector3d A31;
        A21.clear();
        A31.clear();
        Vector3d A211;
        Vector3d A311;
        A211.clear();
        A311.clear();
        Vector3d a21;
        Vector3d a31;
        a21.clear();
        a31.clear();
        Vector3d a211;
        Vector3d a311;
        a211.clear();
        a311.clear();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A21[i] += mat_RodLam_der(i, j) * _N0[j];
                A31[i] += mat_RodLam_der(i, j) * _V0[j];
                a21[i] += mat_rodlamRodLam_der(i, j) * _N0[j];
                a31[i] += mat_rodlamRodLam_der(i, j) * _V0[j];
                A211[i] += mat_RodLam_derder(i, j) * _N0[j];
                A311[i] += mat_RodLam_derder(i, j) * _V0[j];
                a211[i] += mat_rodlamRodLam_derder(i, j) * _N0[j];
                a311[i] += mat_rodlamRodLam_derder(i, j) * _V0[j];
            }
        }

        float Bn_1 = inner_prod(A211, _R1) + inner_prod(A21, _R2);
        float bn_1 = inner_prod(a211, _r1) + inner_prod(a21, _r2);
        float Bv_1 = inner_prod(A311, _R1) + inner_prod(A31, _R2);
        float bv_1 = inner_prod(a311, _r1) + inner_prod(a31, _r2);

        _shear_force_n = bn_1 - Bn_1;
        _shear_force_v = bv_1 - Bv_1;
    }

    void IsogeometricBeamElement::stress_res_lin(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Vector3d&_f, Vector3d&_m)
    {
        const auto& r_geometry = GetGeometry();
        const float emod = this->GetProperties()[YOUNG_MODULUS];
        const float poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const float gmod = emod / (2.0 * (1.0 + poisson_ratio));
        const float area = this->GetProperties()[CROSS_AREA];
        const float height = this->GetProperties()[HEIGHT];
        const float width = this->GetProperties()[WIDTH];
        const float m_inert_z = this->GetProperties()[I_Z];
        const float m_inert_y = this->GetProperties()[I_Y];
        const float mt_iniert = this->GetProperties()[I_T];
        Vector3d t0_0 = this->GetProperties()[T_0];



        Vector func = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        Vector deriv = column(r_geometry.ShapeFunctionDerivatives(1, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector deriv2 = column(r_geometry.ShapeFunctionDerivatives(2, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector deriv3 = column(r_geometry.ShapeFunctionDerivatives(3, integration_point_index, this->GetIntegrationMethod()), 0);

        //declarations
        Vector3d R_1;  //1st derivative of the curve undeformed config
        Vector3d R_2;  //2nd derivative of the curve undeformed config
        Vector3d R_3;  //3rd derivative of the curve undeformed config
        float A;                //length of the base vector
        float B;                //length of the curvature of the undeformed config

        float B_n;
        float B_v;
        float C_12;
        float C_13;
        float Phi;
        float Phi_der;
        float Phi_der2;

        Vector3d r_1;  //1st derivative of the curve undeformed config
        Vector3d r_2;  //2nd derivative of the curve undeformed config
        Vector3d r_3;  //3rd derivative of the curve undeformed config
        float a;                //length of the base vector
        float b;                //length of the curvature of the undeformed config

        float b_n;
        float b_v;
        float c_12;
        float c_13;
        float phi;
        float phi_der;

        // Material and cross section
        float emod_A = emod * area;
        float emod_I_v = emod * m_inert_y;
        float emod_I_n = emod * m_inert_z;
        float gmod_It = gmod * mt_iniert;

        //prestresses
        float prestress = 0;//not yet implemented
        float prestress_bend1 = 0;//not yet implemented
        float prestress_bend2 = 0;//not yet implemented
        float prestress_tor = 0;//not yet implemented

        //Vector3d n;
        //Vector3d v;
        //Vector3d N;
        //Vector3d V;
        Vector3d N0;
        N0.clear(); //new
        Vector3d V0;
        V0.clear(); //new

        // get previous results
        double tmp_ini_dof;

        Phi = 0;
        Phi_der = 0;
        Phi_der2 = 0;
        phi = 0;
        phi_der = 0;

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            tmp_ini_dof = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, rCurrentProcessInfo.GetSolutionStepIndex()); //maybe this is wrong ERROR was r_geomtry[integration_point_index] before
            phi += func(i) * tmp_ini_dof;
            phi_der += deriv[i] * tmp_ini_dof;
        }

        //compute configurations
        comp_Geometry_reference(deriv, deriv2, deriv3, R_1, R_2, R_3, A, B);
        comp_Geometry_reference_cross_section(R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
        //phi=phi+Phi;
        //phi_der+=Phi_der;
        comp_Geometry_initial(deriv, deriv2, deriv3, r_1, r_2, r_3, a, b);
        comp_Geometry_actual_cross_section(r_1, R_1, r_2, R_2, n, v, N0, V0, b_n, b_v, c_12, c_13, phi, phi_der, Phi, Phi_der);

        float Apow2 = pow(A, 2);

        // variation of the axial strain 
        S_eps_var = comp_epsilon_dof(R_1, deriv);

        comp_dof_lin(S_curv_n_var, S_curv_v_var, S_torsion_n_var, S_torsion_v_var, S_shear_n_var, S_shear_v_var, r_1, R_1, r_2, R_2, r_3, R_3, N0, V0, func, deriv, deriv2, deriv3, 0.0, 0.0, 0.0, Phi, Phi_der, Phi_der2);

        //get physical quantities
        //axial strain
        S_eps_var = S_eps_var / Apow2;
        //curvature
        S_curv_n_var = S_curv_n_var / Apow2;
        S_curv_v_var = S_curv_v_var / Apow2;

        //derivative of curvature for shear force
        S_shear_n_var = S_shear_n_var / Apow2;
        S_shear_v_var = S_shear_v_var / Apow2;

        //torsion
        S_torsion_n_var = S_torsion_n_var / A;
        S_torsion_v_var = S_torsion_v_var / A;

        //strains
        float E11_n = 0;
        float E11_cur_n = 0;
        float E11_cur_v = 0;
        float E12E13 = 0;
        float E11_cur1_n = 0;
        float E11_cur1_v = 0;


        Vector3d displacement;
        float rot;

        Vector tmp_dof;
        tmp_dof.resize(Dof_Node, false);
        tmp_dof.clear();

        Vector3d coords;
        for (int r = 0;r < N_Dof;r++) //in the case
        {
            int xyz = r % Dof_Node;   //0 ->disp_x; 1 ->disp_y; 2 ->disp_z; 3 ->rot_tan
            int i = r / Dof_Node;     // index for the nodes
            
            coords = r_geometry[i].Coordinates(); //check if this really gives the coodrinates of nodes 1-6 and not the integration point coordinates!!
            displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, rCurrentProcessInfo.GetSolutionStepIndex());
            rot = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, rCurrentProcessInfo.GetSolutionStepIndex());

            tmp_dof[0] = displacement[0];
            tmp_dof[1] = displacement[1];
            tmp_dof[2] = displacement[2];
            tmp_dof[3] = rot;

            E11_n += S_eps_var[r] * tmp_dof[xyz];
            E11_cur_n += S_curv_n_var[r] * tmp_dof[xyz];
            E11_cur_v += S_curv_v_var[r] * tmp_dof[xyz];
            E12E13 += -0.5 * S_torsion_n_var[r] * tmp_dof[xyz];
            E12E13 += 0.5 * S_torsion_v_var[r] * tmp_dof[xyz];
            E11_cur1_n += S_shear_n_var[r] * tmp_dof[xyz];
            E11_cur1_v += S_shear_v_var[r] * tmp_dof[xyz];
        }

        float n_tmp = E11_n * area * emod;
        float mn_tmp = E11_cur_n * m_inert_y * emod;
        float mv_tmp = E11_cur_v * m_inert_z * emod;

        //stresses
        //_m.resize(3);
        //_q.resize(2);
        _f[0] = ((area + m_inert_y * pow(B_n / Apow2, 2) - m_inert_z * pow(B_v / Apow2, 2)) * E11_n + E11_cur_n * m_inert_y * B_n / Apow2 + E11_cur_v * m_inert_z * B_v / Apow2) * emod + prestress * area;
        _m[1] = (E11_cur_n * m_inert_y/*+m_inert_y*B_n/Apow2*E11_n*/) * emod + prestress_bend1;
        _m[2] = (E11_cur_v * m_inert_z/*+m_inert_z*B_v/Apow2*E11_n*/) * emod + prestress_bend2;
        _m[0] = E12E13 * (gmod_It/*+mn_tmp*(m_inert_z+m_inert_y)/area-mv_tmp*(m_inert_z+m_inert_y)/area*/)/*+emod/2.0*1.26394444444445e-09*E12E13pow3*/+prestress_tor;
        _f[1] = emod * m_inert_y * E11_cur1_n / A;
        _f[2] = emod * m_inert_z * E11_cur1_v / A;
    }

    void IsogeometricBeamElement::stress_res_nln(const ProcessInfo& rCurrentProcessInfo, IndexType integration_point_index, Vector3d& _f, Vector3d& _m)
    {
        //TODO: Also return values between the integration points through changing the shape functions
        const auto& r_geometry = GetGeometry();
        const float emod = this->GetProperties()[YOUNG_MODULUS];
        const float poisson_ratio = this->GetProperties()[POISSON_RATIO];
        const float gmod = emod / (2.0 * (1.0 + poisson_ratio));
        const float area = this->GetProperties()[CROSS_AREA];
        const float height = this->GetProperties()[HEIGHT];
        const float width = this->GetProperties()[WIDTH];
        const float m_inert_z = this->GetProperties()[I_Z];
        const float m_inert_y = this->GetProperties()[I_Y];
        const float mt_iniert = this->GetProperties()[I_T];
        Vector3d t0_0 = this->GetProperties()[T_0];
        
        Vector func = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), 0);
        Vector deriv = column(r_geometry.ShapeFunctionDerivatives(1, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector deriv2 = column(r_geometry.ShapeFunctionDerivatives(2, integration_point_index, this->GetIntegrationMethod()), 0);
        Vector deriv3 = column(r_geometry.ShapeFunctionDerivatives(3, integration_point_index, this->GetIntegrationMethod()), 0);

        //declarations
        Vector3d R_1;  //1st derivative of the curve undeformed config
        Vector3d R_2;  //2nd derivative of the curve undeformed config
        Vector3d R_3;  //3rd derivative of the curve undeformed config
        float A;                //length of the base vector
        float B;                //length of the curvature of the undeformed config
        float B_n;
        float B_v;
        float C_12;
        float C_13;
        float Phi;
        float Phi_der;
        float Phi_der2;

        Vector3d r_1;  //1st derivative of the curve undeformed config
        Vector3d r_2;  //2nd derivative of the curve undeformed config
        Vector3d r_3;  //3rd derivative of the curve undeformed config
        float a;                //length of the base vector
        float b;                //length of the curvature of the undeformed config
        float b_n;
        float b_v;
        float c_12;
        float c_13;
        float phi;
        float phi_der;
        float phi_der2;
        //cfloat phi_0;


        // Material and cross section
        float emod_A = emod * area;
        float emod_I_v = emod * m_inert_y;
        float emod_I_n = emod * m_inert_z;
        float gmod_It = gmod * mt_iniert;

        //prestresses
        float prestress = 0;//not yet implemented
        float prestress_bend1 = 0;//not yet implemented
        float prestress_bend2 = 0;//not yet implemented
        float prestress_tor = 0;//not yet implemented
        //bool prestress_bend1_auto = Prop_Ptr->get_Prestress_Bend1_Auto();
        //bool prestress_bend2_auto = Prop_Ptr->get_Prestress_Bend2_Auto();
        //bool prestress_tor_auto = Prop_Ptr->get_Prestress_Tor_Auto();

        //Vector3d n;
        //Vector3d v;
        //Vector3d N;
        //Vector3d V;
        //cfloat phi;
        //cfloat phi_der;

        //compute configurations
        comp_Geometry_reference(deriv, deriv2, deriv3, R_1, R_2, R_3, A, B);
        comp_Geometry_actual(rCurrentProcessInfo, deriv, deriv2, deriv3, r_1, r_2, r_3, a, b);

        Vector3d N0;
        N0.clear();
        Vector3d V0;
        V0.clear();

        // get previous results
        double tmp_ini_dof;

        phi = 0;
        phi_der = 0;
        phi_der2 = 0;
        Phi = 0;
        Phi_der = 0;
        Phi_der2 = 0;

        for (size_t i = 0;i < r_geometry.size();i++)
        {
            tmp_ini_dof = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, rCurrentProcessInfo.GetSolutionStepIndex()); //maybe this is wrong ERROR was r_geomtry[integration_point_index] before
            phi += func(i) * tmp_ini_dof;
            phi_der += deriv[i] * tmp_ini_dof;
            phi_der2 += deriv2[i] * tmp_ini_dof;
        }

        comp_Geometry_reference_cross_section(R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der);
        //phi=phi+Phi;
        //phi_der+=Phi_der;
        comp_Geometry_actual_cross_section(r_1, R_1, r_2, R_2, n, v, N0, V0, b_n, b_v, c_12, c_13, phi, phi_der, Phi, Phi_der);
        float shear_force_n = 0;
        float shear_force_v = 0;
        comp_transverse_shear_force_nln(r_1, R_1, r_2, R_2, r_3, R_3, N0, V0, phi, phi_der, phi_der2, Phi, Phi_der, Phi_der2, shear_force_n, shear_force_v);
        float Apow2 = pow(A, 2);
        float Apow4 = pow(Apow2, 2);

        //stresses
        float E11_m = 0.5 * (pow(a, 2) - pow(A, 2));   //Green Lagrange formulation (strain)
        float E11_cur_n = (b_n - B_n);
        float E11_cur_v = (b_v - B_v);
        float E12 = (c_12 - C_12);
        float E13 = (c_13 - C_13);

        float S11_m = prestress * area + E11_m * emod_A / Apow2;                      //normal force
        float S11_n = prestress_bend1 + E11_cur_n * emod_I_v / Apow2;                        //bending moment n
        float S11_v = prestress_bend2 + E11_cur_v * emod_I_n / Apow2;                        //bending moment v
        float S12 = 0.5 * (prestress_tor + E12 * gmod_It / A);                      //0.5 torsional moment
        float S13 = 0.5 * (prestress_tor + E13 * gmod_It / A);                      //0.5 torsional moment

        float S1213 = -0.5 * (E12 - E13) * gmod_It / A + prestress_tor;
        shear_force_n *= emod_I_v / Apow2;
        shear_force_v *= emod_I_n / Apow2;

        //// variation of the axial strain 



        _f[0] = S11_m;
        _m[1] = -S11_v;
        _m[2] = -S11_n;
        _m.resize(2);
        _f[0] = ((area + m_inert_y * pow(B_n / Apow2, 2) - m_inert_z * pow(B_v / Apow2, 2)) * E11_m) * emod / Apow2 + prestress * area;    
        _f[0] += E11_cur_n * m_inert_y * b_n / Apow4 * emod;
        _f[0] += E11_cur_v * m_inert_z * b_v / Apow4 * emod;
        _f[0] *= a / A;
        _m[1] = ((E11_cur_n * m_inert_y/*+m_inert_y*B_n/Apow2*E11_m*/) * emod / Apow2 + prestress_bend1) * a / A;
        _m[2] = ((E11_cur_v * m_inert_z/*+m_inert_z*B_v/Apow2*E11_m*/) * emod / Apow2 + prestress_bend2) * a / A;
        //_t[0] = 0.5*S12 + 0.5*S13;
        _m[0] = S1213 * a / A;
        _f[1] = shear_force_n / Apow2 * a;    // /A * a / A transformation to local cartesian and Cauchy
        _f[2] = shear_force_v / Apow2 * a;
    }



    ///@}
} // namespace Kratos
