//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/variables.h"
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_element_utilities.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


/**
    * @brief Auxiliary data container
    * This class serves as a temporary data container for the wall-related data
    */
template<std::size_t TDim, std::size_t TNumNodes>
struct NitscheNavierSlipWallDataContainer
{
    double DynamicViscosity;  //TODO remove?
    double ElementSize;
    array_1d<double,3> Normal;
    Vector GaussPtsWeights;
    Matrix ShapeFunctionsContainer;
    std::vector<Matrix> ShapeFunctionDerivativesContainer;
    array_1d<double, TNumNodes> NodalSlipLength;
    array_1d<double, TNumNodes> NodalPenalty;
    array_1d<double, TNumNodes> NodalDensity;
    array_1d<double, TNumNodes> NodalPressure;
    array_1d<array_1d<double,3>,TNumNodes> NodalVelocities;
    array_1d<array_1d<double,3>,TNumNodes> NodalMeshVelocities;

    // Store variables and parameters in the container to avoid recomputing them multiple times
    void Initialize(const Condition& rCondition)
    {
        static_assert(TDim == 2 || TDim == 3, "Unsupported problem dimension and condition node numbers.");
        static_assert(TDim == TNumNodes, "Unsupported condition node numbers for problem dimension.");

        const auto& r_geom = rCondition.GetGeometry();

        // Get dynamic viscosity from parent element properties
        // Note that this assumes constant viscosity within the element
        const auto& r_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)[0];
        DynamicViscosity = r_parent_element.GetProperties().GetValue(DYNAMIC_VISCOSITY);
        // NOTE that this only works for simplex elements (TNumNode+1)!
        //ElementSize = ElementSizeCalculator<TDim,TNumNodes+1>::MinimumElementSize(r_parent_element.GetGeometry());
        ElementSize = r_geom.Length();

        // Compute condition unit normal vector
        // Note that in here we are assuming a constant normal over the entire condition
        Normal = r_geom.UnitNormal(0, GeometryData::IntegrationMethod::GI_GAUSS_1);

        // Calculate the required Gauss points integration data and resize the containers
        const GeometryData::IntegrationMethod gp_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_pts = r_geom.IntegrationPoints(gp_method);
        const std::size_t n_gauss = r_integration_pts.size();
        ShapeFunctionDerivativesContainer.resize(n_gauss);
        GaussPtsWeights.resize(n_gauss);
        ShapeFunctionsContainer.resize(n_gauss, TNumNodes);
        Vector detJs(n_gauss);
        r_geom.DeterminantOfJacobian(detJs, gp_method);

        // Loop over the Gauss points and compute integration weights and Gauss point shape function values and shape function derivatives
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get weight of the current Gauss point
            GaussPtsWeights[i_gauss] = detJs[i_gauss] * r_integration_pts[i_gauss].Weight();

            // Get DN_DX of the element at the Gauss point
            Matrix aux_DN_DXi, aux_J, aux_J_pinv;
            r_geom.ShapeFunctionsLocalGradients(aux_DN_DXi, r_integration_pts[i_gauss]);  // 2D: [-0.5, 0.5]
            r_geom.Jacobian(aux_J, r_integration_pts[i_gauss]);  // 2D: 0.5*[[dx], [dy]]
            // NOTE that for line condition of 2D problem, the Jacobian is a 2x1 matrix and for triangle condition in 3D, the Jacobian is a 3x2 matrix, so the inverse does not exist.
            // Instead, we can compute the Moore-Penrose pseudo-inverse of the Jacobian using the formula: J^+ = (J^T * J)^-1 * J^T
            BoundedMatrix<double, TDim, TDim> aux_JT = prod(trans(aux_J), aux_J);
            double aux_detJT;
            BoundedMatrix<double, TDim, TDim> aux_JT_inv;
            MathUtils<double>::InvertMatrix(aux_JT, aux_JT_inv, aux_detJT);
            aux_J_pinv = prod(aux_JT_inv, trans(aux_J));
            ShapeFunctionDerivativesContainer[i_gauss] = prod(aux_DN_DXi, aux_J_pinv);
        }
        ShapeFunctionsContainer = r_geom.ShapeFunctionsValues(gp_method);  // matrix of SF values F_{ij}, where i is integration pt index and j is SF index

        // Save the nodal velocities and interpolate the slip length at current integration point
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
            const auto& r_node = r_geom[i_node];
            NodalSlipLength[i_node] = r_node.GetValue(SLIP_LENGTH);
            NodalPenalty[i_node] = r_node.GetValue(PENALTY_COEFFICIENT);
            NodalDensity[i_node] = r_node.FastGetSolutionStepValue(DENSITY);
            NodalPressure[i_node] = r_node.FastGetSolutionStepValue(PRESSURE);
            NodalVelocities[i_node] = r_node.FastGetSolutionStepValue(VELOCITY);
            NodalMeshVelocities[i_node] = r_node.FastGetSolutionStepValue(MESH_VELOCITY);
        }
    }
};

/**
 * @brief Navier-slip law LHS and RHS contribution implementation
 * This class implements the LHS and RHS contributions of the Navier-slip wall model using Nitsche's method.
 * This class should be used in combination with an incompressible Navier-Stokes (or Stokes) condition)
 * implementing the remaining terms (see @NavierStokesWallCondition).
 * More information about the Nitsche imposition of Navier-slip conditions can be found in https://doi.org/10.1016/j.cma.2017.10.023
 * The same implementation is used for CutFEM (@EmbeddedFluidElement and @EmbeddedFluidElementDiscontinuous)
 * and a similar one for the shifted boundary method (see @ShiftedBoundaryWallCondition).
 * @tparam TDim Number of dimensions
 * @tparam TNumNodes Number of nodes
 */
template<std::size_t TDim, std::size_t TNumNodes>
class NitscheNavierSlipWallLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NitscheNavierSlipWallLaw
    KRATOS_CLASS_POINTER_DEFINITION(NitscheNavierSlipWallLaw);

    static constexpr std::size_t BlockSize = TDim+1;

    static constexpr std::size_t VoigtSize = 3 * (TDim-1);  // 3 in 2D, 6 in 3D - StrainSize of FluidElementData

    static constexpr std::size_t LocalSize = TNumNodes*BlockSize;

    using SizeType = Condition::SizeType;

    using IndexType = Condition::IndexType;

    using VectorType = Condition::VectorType;

    using MatrixType = Condition::MatrixType;

    using WallDataContainerType = NitscheNavierSlipWallDataContainer<TDim,TNumNodes>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NitscheNavierSlipWallLaw() = delete;

    /// Copy constructor.
    NitscheNavierSlipWallLaw(NitscheNavierSlipWallLaw const& rOther) = delete;

    /// Destructor.
    ~NitscheNavierSlipWallLaw() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Add the LHS and RHS Navier-slip contributions using Nitsche's method
     * This method adds the Navier-slip LHS and RHS contributions to the provided matrices using Nitsche's method
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param rRightHandSideVector Reference to the output RHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo con
     */
    static void AddWallModelLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Create and fill the wall law auxiliary data container
        WallDataContainerType wall_law_data;
        wall_law_data.Initialize(*pCondition);

        // KRATOS_WATCH(wall_law_data.ElementSize);
        // KRATOS_WATCH(wall_law_data.Normal);
        // KRATOS_WATCH(wall_law_data.GaussPtsWeights);
        // KRATOS_WATCH(wall_law_data.ShapeFunctionsContainer);
        // KRATOS_WATCH(wall_law_data.ShapeFunctionDerivativesContainer.size());
        // KRATOS_WATCH(wall_law_data.ShapeFunctionDerivativesContainer[0]);
        // KRATOS_WATCH(wall_law_data.NodalSlipLength);
        // KRATOS_WATCH(wall_law_data.NodalPenalty);
        // KRATOS_WATCH(wall_law_data.NodalDensity);
        // KRATOS_WATCH(wall_law_data.NodalPressure);
        // KRATOS_WATCH(wall_law_data.NodalVelocities);
        // KRATOS_WATCH(wall_law_data.NodalMeshVelocities);

        // Set whether the shear stress term is adjoint consistent (1.0) or adjoint inconsistent (-1.0) for Nitsche imposition
        // NOTE that adjoint inconsistent formulation enjoys improved inf-sup stability for any value of the penalty constant (gamma) according to Winter et al. (2018)
        constexpr double adjoint_consistency = -1.0;


        // ---------------------------------------------------------------
        //   Obtain solution values and precompute projection matrices
        // ---------------------------------------------------------------

        // Obtain the previous iteration velocity and pressure solution for all nodes
        // and the embedded nodal velocity EMBEDDED_VELOCITY to account for a moving boundary (FM-ALE)
        Vector unknown_values(LocalSize);
        Vector embedded_velocity(LocalSize);
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {;
            const IndexType base = i_node * BlockSize;
            for (IndexType d = 0; d < TDim; ++d) {
                unknown_values[base + d] = wall_law_data.NodalVelocities[i_node][d];
                embedded_velocity[base + d] = wall_law_data.NodalMeshVelocities[i_node][d];
            }
            unknown_values[base + TDim] = wall_law_data.NodalPressure[i_node];
            embedded_velocity[base + TDim] = 0.0;
        }

        // Precompute projection matrices
        // Normal projection in Voigt notation (A matrix)
        BoundedMatrix<double, TDim, VoigtSize> voigt_normal_proj_matrix;
        FluidElementUtilities<3>::VoigtTransformForProduct(wall_law_data.Normal, voigt_normal_proj_matrix);
        // Normal projection matrix: n x n
        BoundedMatrix<double, TDim, TDim> Pnorm;
        FluidElementUtilities<3>::SetNormalProjectionMatrix(wall_law_data.Normal, Pnorm);
        //SetNormalProjectionMatrix(wall_law_data.Normal, Pnorm);
        // Tangential projection matrix: I - n x n
        BoundedMatrix<double, TDim, TDim> Ptang;
        FluidElementUtilities<3>::SetTangentialProjectionMatrix(wall_law_data.Normal, Ptang);
        //SetTangentialProjectionMatrix(wall_law_data.Normal, Ptang);

        // ---------------------------------------------------------------
        //   Iterate over the Gauss points of the Condition to add all contributions to the local system
        // ---------------------------------------------------------------
        const SizeType n_gauss = wall_law_data.GaussPtsWeights.size();
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Get current Gauss point data
            const double w_gauss = wall_law_data.GaussPtsWeights[i_gauss];
            const auto& N_gauss = row(wall_law_data.ShapeFunctionsContainer, i_gauss);
            const auto& dN_DX_gauss = wall_law_data.ShapeFunctionDerivativesContainer[i_gauss];

            // Set the shape functions auxiliary matrix
            Matrix N_matrix = ZeroMatrix(TDim, LocalSize);
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node){
                for (IndexType d = 0; d < TDim; ++d){
                    N_matrix(d, i_node*BlockSize + d) = N_gauss(i_node);
                }
            }


            // ---------------------------------------------------------------
            //   Get strain matrix, constitutive matrix and precompute frequently used products
            // ---------------------------------------------------------------

            // Set the current integration point's strain matrix (SF derivatives at all nodes)
            Matrix B_matrix = ZeroMatrix(VoigtSize, LocalSize);
            CalculateStrainMatrix<TDim>(dN_DX_gauss, B_matrix);  // see FluidElementUtilities

            // Get C Matrix from constitutive law, see FluidElementData and FluidElement::CalculateMaterialResponse
            const auto& r_parent_element = pCondition->GetValue(NEIGHBOUR_ELEMENTS)[0];
            const auto p_constitutive_law =  r_parent_element.GetProperties()[CONSTITUTIVE_LAW];
            ConstitutiveLaw::Parameters constitutive_law_values(r_parent_element.GetGeometry(), r_parent_element.GetProperties(), rCurrentProcessInfo);
            Vector strain_rate = prod(B_matrix, unknown_values);
            Vector shear_stress(VoigtSize);
            Matrix C_matrix = ZeroMatrix(VoigtSize, VoigtSize);

            // Set constitutive law flags
            Flags& r_cl_options = constitutive_law_values.GetOptions();
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);  // not being used (?)
            r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            // Set constitutive law values
            //constitutive_law_values.SetShapeFunctionsValues(N_gauss);
            //constitutive_law_values.SetShapeFunctionsDerivatives(dN_DX_gauss);
            constitutive_law_values.SetStrainVector(strain_rate);
            constitutive_law_values.SetStressVector(shear_stress);    //this is an output parameter
            constitutive_law_values.SetConstitutiveMatrix(C_matrix);  //this is an output parameter

            // Calculate material response and effective viscosity
            //NOTE that here we assume that only one constitutive law is employed for all of the integration points.
            //this is ok under the hypothesis that no history dependent behaviour is employed
            double effective_viscosity = 0.0;
            p_constitutive_law->CalculateMaterialResponseCauchy(constitutive_law_values);  // computes the constitutive matrix C (and the shear stress vector)
            p_constitutive_law->CalculateValue(constitutive_law_values, EFFECTIVE_VISCOSITY, effective_viscosity);

            // Precompute frequently used products
            // C * B  (VoigtSize x LocalSize)
            const Matrix CB = prod(C_matrix, B_matrix);
            // A * C * B  (TDim x LocalSize)
            const Matrix ACB  = prod(voigt_normal_proj_matrix, CB);
            // B^T * A^T  (LocalSize x TDim)
            const Matrix BtAt = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));


            // ---------------------------------------------------------------
            //   Build aux_LHS_traction by accumulating boundary traction contributions of the Gauss point
            // ---------------------------------------------------------------
            //TODO Is this somehow already included in the other contributions? (see AddBoundaryTraction of fluid element)
            //   OR use AddBoundaryTraction of the element?!
            // r_parent_element.UpdateIntegrationPointData(data, 0, w_gauss, N_parent, DN_DX_parent);
            // r_parent_element.AddBoundaryTraction(data, wall_law_data.Normal, rLeftHandSideMatrix, rRightHandSideVector);

            MatrixType aux_LHS_traction = ZeroMatrix(LocalSize, LocalSize);

            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                const double w_Ni = w_gauss * N_gauss[i_node];
                const IndexType base_i = i_node * BlockSize;

                // --- Shear stress contribution ---
                // aux_LHS -= w * N_u^T * ACB
                for (IndexType di = 0; di < TDim; ++di) {
                    const IndexType row = base_i + di;
                    for (IndexType col = 0; col < LocalSize; ++col) {
                        aux_LHS_traction(row, col) -= w_Ni * ACB(di, col);
                    }
                }

                // --- Pressure contribution ---
                // aux_LHS += w * N_u^T * n_d * N_p
                for (IndexType di = 0; di < TDim; ++di) {
                    const IndexType row = base_i + di;
                    const double w_Ni_nd = w_Ni * wall_law_data.Normal[di];
                    for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                        const IndexType pres_col = j_node * BlockSize + TDim;
                        aux_LHS_traction(row, pres_col) += w_Ni_nd * N_gauss[j_node];
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////
            // Add boundary traction contributions of the integration point to the local system (residual-based formulation with RHS=f_gamma-LHS*previous_solution)
            noalias(rLeftHandSideMatrix)  += aux_LHS_traction;
            noalias(rRightHandSideVector) -= prod(aux_LHS_traction, unknown_values);


            // ---------------------------------------------------------------
            //   Build aux_LHS_vel and aux_LHS_tau by accumulating Navier-slip contributions of the Gauss point regarding the velocity and shear stress
            // ---------------------------------------------------------------
            MatrixType aux_LHS_vel = ZeroMatrix(LocalSize, LocalSize);
            MatrixType aux_LHS_tau = ZeroMatrix(LocalSize, LocalSize);

            /////////////////////////////////////////////////////////////////////////////////
            // 1. Normal penalty: pen * w * (N_i n_di) (N_j n_dj)

            // Compute the Nitsche normal imposition penalty coefficient
            const double pen_coeff_normal = NitscheNavierSlipWallLaw<TDim, TNumNodes>::ComputeNormalPenaltyCoefficient(
                wall_law_data, N_gauss, effective_viscosity, rCurrentProcessInfo.GetValue(DELTA_TIME));
            const double pen_normal_w = pen_coeff_normal * w_gauss;

            // Add normal velocity penalty contribution of the integration point
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    for (IndexType di = 0; di < TDim; ++di) {
                        const IndexType row = i_node * BlockSize + di;
                        for (IndexType dj = 0; dj < TDim; ++dj){
                            const IndexType col = j_node * BlockSize + dj;
                            aux_LHS_vel(row, col) += pen_normal_w * N_gauss(i_node) * wall_law_data.Normal(di) * wall_law_data.Normal(dj) * N_gauss(j_node);
                        }
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////////////////
            // 2. Normal Nitsche symmetric counterpart stabilization (shear stress + pressure)

            // Fill the pressure to Voigt notation operator normal projected matrix
            Matrix trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, TDim);
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node){
                for (IndexType d = 0; d < TDim; ++d){
                    trans_pres_to_voigt_matrix_normal_op(i_node*BlockSize + TDim, d) = N_gauss[i_node]*wall_law_data.Normal(d);
                }
            }

            // Compute some integration point auxiliary matrices
            const Matrix AtPnorm = prod(trans(voigt_normal_proj_matrix), Pnorm);
            const Matrix BtCtAtPnorm = prod(trans(CB), AtPnorm);

            // Contribution coming from the shear stress operator (traction vector normal component)
            noalias(aux_LHS_vel) -= adjoint_consistency * w_gauss * prod(BtCtAtPnorm, N_matrix);

            // Contribution coming from the pressure terms
            const Matrix aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, Pnorm);
            noalias(aux_LHS_vel) -= w_gauss * prod(aux_matrix_VPnorm, N_matrix);

            /////////////////////////////////////////////////////////////////////////////////
            // 3. Tangential penalty

            // Compute the Nitsche slip tangential penalty coefficients
            auto [pen_coeff_tang_1, pen_coeff_tang_2] = NitscheNavierSlipWallLaw<TDim, TNumNodes>::ComputeTangentialPenaltyCoefficients(
                wall_law_data, N_gauss, effective_viscosity);
            const double pen_tang_w_1 = pen_coeff_tang_1 * w_gauss;
            const double pen_tang_w_2 = pen_coeff_tang_2 * w_gauss;

            const Matrix Nt_Ptang = prod(trans(N_matrix), Ptang);

            // Contribution coming from the tangential velocity (towards zero for slip)
            // pen * w * N^T * P_tang * N
            noalias(aux_LHS_vel) += pen_tang_w_2 * prod(Nt_Ptang, N_matrix);

            // Contribution coming from the traction vector tangential component (no contribution for no slip)
            // pen * w * N^T * P_tang * A*C*B
            noalias(aux_LHS_tau) += pen_tang_w_1 * prod(Nt_Ptang, ACB);

            /////////////////////////////////////////////////////////////////////////////////
            // 4. Tangential Nitsche symmetric counterpart stabilization

            // Compute the stabilization coefficients
            auto [nitsche_coeff_tang_1, nitsche_coeff_tang_2] = NitscheNavierSlipWallLaw<TDim, TNumNodes>::ComputeTangentialNitscheCoefficients(
                wall_law_data, N_gauss, effective_viscosity);
            const double nitsche_tang_w_1 = nitsche_coeff_tang_1 * w_gauss * adjoint_consistency * 2.0;
            const double nitsche_tang_w_2 = nitsche_coeff_tang_2 * w_gauss * adjoint_consistency * 2.0;

            // BtAtPtang = BtAt * Ptang  (LocalSize x TDim)
            const Matrix BtAtPtang = prod(BtAt, Ptang);

            // Contribution coming from the tangential velocity (towards zero for slip)
            // (-1) * pen * w * 2 * B^T*A^T*P_tang * N (exploit sparsity)
            for (IndexType i = 0; i < LocalSize; ++i) {
                for (IndexType j_node = 0; j_node < TNumNodes; ++j_node) {
                    const double Nj = N_gauss(j_node);
                    const IndexType base_j = j_node * BlockSize;
                    for (IndexType d = 0; d < TDim; ++d) {
                        aux_LHS_vel(i, base_j + d) -= nitsche_tang_w_2 * BtAtPtang(i, d) * Nj;
                    }
                }
            }

            // Contribution coming from the traction vector tangential component (no contribution for no slip)
            // (-1) * pen * w * 2 * B^T*A^T*P_tang * A*C*B
            noalias(aux_LHS_tau) -= nitsche_tang_w_1 * prod(BtAtPtang, ACB);

            /////////////////////////////////////////////////////////////////////////////////
            // Add Nitsche Navier-slip contributions of the integration point to the local system (residual-based formulation with RHS=f_gamma-LHS*previous_solution)
            // NOTE that we subtract the embedded velocity here to account for a moving boundary (FM-ALE)
            noalias(rLeftHandSideMatrix)  += aux_LHS_vel + aux_LHS_tau;
            noalias(rRightHandSideVector) -= prod(aux_LHS_vel, unknown_values - embedded_velocity);
            noalias(rRightHandSideVector) -= prod(aux_LHS_tau, unknown_values);
        }
    }

    /**
     * @brief Add the LHS Navier-slip contribution
     * This method adds the Navier-slip LHS contribution to the provided matrix
     * @param rLeftHandSideMatrix Reference to the output LHS matrix
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     */
    static void AddWallModelLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        VectorType rRightHandSideVector;
        AddWallModelLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, pCondition, rCurrentProcessInfo);
    }

    /**
     * @brief Add the RHS Navier-slip contribution
     * This method adds the Navier-slip RHS contribution to the provided vector
     * @param rRightHandSideVector Reference to the output RHS vector
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to current ProcessInfo container
     */
    static void AddWallModelRightHandSide(
        VectorType& rRightHandSideVector,
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType rLeftHandSideMatrix;
        AddWallModelLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, pCondition, rCurrentProcessInfo);
    }

    /**
     * @brief Check function
     * This function checks the current wall law input parameters
     * @param pCondition Pointer to the current condition
     * @param rCurrentProcessInfo Reference to the ProcessInfo data container
     * @return int Check status
     */
    static int Check(
        const Condition* pCondition,
        const ProcessInfo& rCurrentProcessInfo)
    {
        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "NitscheNavierSlipWallLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NitscheNavierSlipWallLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Private operations
    ///@{

    static double ComputeNormalPenaltyCoefficient(
        const WallDataContainerType& rWallLawData,
        const array_1d<double,TNumNodes>& rN,
        const double EffectiveViscosity,
        const double DeltaTime)
    {
        // Get the nodal magnitudes at the current Gauss point
        double gp_rho = 0.0;
        array_1d<double,TDim> gp_v = ZeroVector(TDim);
        double gp_penalty = 0.0;
        for (IndexType i_node = 0;  i_node < TNumNodes; ++i_node) {
            gp_rho += rN[i_node] * rWallLawData.NodalDensity[i_node];
            gp_v += rN[i_node] * rWallLawData.NodalVelocities[i_node];
            gp_penalty += rN[i_node] * rWallLawData.NodalPenalty[i_node];
        }
        const double gp_v_norm = norm_2(gp_v);
        const double gamma = 1.0 / gp_penalty;

        // Compute the Nitsche coefficient (including the Winter stabilization term)
        const double h = rWallLawData.ElementSize;
        const double eff_mu = EffectiveViscosity;  // rWallLawData.DynamicViscosity;

        // Compute the velocity stabilization constant of Schott et al. (doi: 10.1002/fld.4218)
        const double stab_constant_u = eff_mu + gp_rho*gp_v_norm*h / 6.0 + gp_rho*h*h/DeltaTime / 12.0;

        // Compute the Nitsche coefficient (including the stabilization term from Schott et al. (doi: 10.1002/fld.4218))
        const double coeff = (eff_mu + stab_constant_u) / (gamma * h);

        return coeff;
    }

    static std::pair<const double, const double> ComputeTangentialPenaltyCoefficients(
        const WallDataContainerType& rWallLawData,
        const array_1d<double,TNumNodes>& rN,
        const double EffectiveViscosity)
    {
        double gp_penalty = 0.0;
        double gp_slip_length = 0.0;
        for (IndexType i_node = 0;  i_node < TNumNodes; ++i_node) {
            gp_penalty += rN[i_node] * rWallLawData.NodalPenalty[i_node];
            gp_slip_length += rN[i_node] * rWallLawData.NodalSlipLength[i_node];
        }
        const double gamma = 1.0 / gp_penalty;

        const double h = rWallLawData.ElementSize;
        const double eff_mu = EffectiveViscosity;  // rWallLawData.DynamicViscosity;

        const double penalty_coeff = 1.0 / (gp_slip_length + gamma*h);
        const double coeff_1 = penalty_coeff * gp_slip_length;
        const double coeff_2 = penalty_coeff * eff_mu;

        return std::make_pair(coeff_1, coeff_2);
    }

    static std::pair<const double, const double> ComputeTangentialNitscheCoefficients(
        const WallDataContainerType& rWallLawData,
        const array_1d<double,TNumNodes>& rN,
        const double EffectiveViscosity)
    {
        double gp_penalty = 0.0;
        double gp_slip_length = 0.0;
        for (IndexType i_node = 0;  i_node < TNumNodes; ++i_node) {
            gp_penalty += rN[i_node] * rWallLawData.NodalPenalty[i_node];
            gp_slip_length += rN[i_node] * rWallLawData.NodalSlipLength[i_node];
        }
        const double gamma = 1.0 / gp_penalty;

        const double h = rWallLawData.ElementSize;
        const double eff_mu = EffectiveViscosity;  // rWallLawData.DynamicViscosity;

        const double penalty_coeff = gamma*h / (gp_slip_length + gamma*h);
        const double coeff_1 = penalty_coeff * gp_slip_length;
        const double coeff_2 = penalty_coeff * eff_mu;

        return std::make_pair(coeff_1, coeff_2);
    }

    // //TODO: Remove this by a common implementation (FluidElementUtilities cannot be used yet)
    // static void SetNormalProjectionMatrix(
    //     const array_1d<double,3>& rUnitNormal,
    //     BoundedMatrix<double,TDim,TDim>& rNormProjMat)
    // {
    //     noalias(rNormProjMat) = ZeroMatrix(TDim,TDim);
    //     for (IndexType d1 = 0; d1 < TDim; ++d1) {
    //         for (IndexType d2 = 0; d2 < TDim; ++d2) {
    //             rNormProjMat(d1,d2) += rUnitNormal[d1]*rUnitNormal[d2];
    //         }
    //     }
    // }

    // //TODO: Remove this by a common implementation (FluidElementUtilities cannot be used yet)
    // static void SetTangentialProjectionMatrix(
    //     const array_1d<double,3>& rUnitNormal,
    //     BoundedMatrix<double,TDim,TDim>& rTangProjMat)
    // {
    //     noalias(rTangProjMat) = IdentityMatrix(TDim,TDim);
    //     for (IndexType d1 = 0; d1 < TDim; ++d1) {
    //         for (IndexType d2 = 0; d2 < TDim; ++d2) {
    //             rTangProjMat(d1,d2) -= rUnitNormal[d1]*rUnitNormal[d2];
    //         }
    //     }
    // }

    /**
     * @brief This function builds the strain matrix from the shape function derivatives utilizing Voigt notation.
     * @param rDN_DX matrix of shape function derivatives at all nodes of the condition
     * @param rB B_matrix/ strain matrix
     */
    template<std::size_t TStrainDim>
    static void CalculateStrainMatrix(
        const Matrix& rDN_DX,
        Matrix& rB)
    {
        static_assert(TStrainDim == 2 || TStrainDim == 3, "Unsupported strain dimension");
        rB.clear();

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node ) {
            const std::size_t col = BlockSize * i_node;

            if constexpr (TStrainDim == 3) {
                rB(0, col  ) = rDN_DX(i_node, 0);
                rB(1, col+1) = rDN_DX(i_node, 1);
                rB(2, col+2) = rDN_DX(i_node, 2);
                rB(3, col  ) = rDN_DX(i_node, 1);
                rB(3, col+1) = rDN_DX(i_node, 0);
                rB(4, col+1) = rDN_DX(i_node, 2);
                rB(4, col+2) = rDN_DX(i_node, 1);
                rB(5, col  ) = rDN_DX(i_node, 2);
                rB(5, col+2) = rDN_DX(i_node, 0);
            } else {
                rB(0, col  ) = rDN_DX(i_node, 0);
                rB(1, col+1) = rDN_DX(i_node, 1);
                rB(2, col  ) = rDN_DX(i_node, 1);
                rB(2, col+1) = rDN_DX(i_node, 0);
            }
        }
    }

    ///@}
}; // Class NitscheNavierSlipWallLaw

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
