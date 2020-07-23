//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4)
#define KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4

// System includes
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/find_nodal_h_process.h"
#include "processes/compute_nodal_gradient_process.h"
#include "solving_strategies/strategies/explicit_solving_strategy_runge_kutta_4.h"
#include "utilities/element_size_calculator.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_processes/shock_detection_process.h"

namespace Kratos
{

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

/** @brief Explicit solving strategy base class
 * @details This is the base class from which we will derive all the explicit strategies (FE, RK4, ...)
 */
template <class TSparseSpace, class TDenseSpace>
class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 : public ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class definition
    typedef ExplicitSolvingStrategyRungeKutta4<TSparseSpace, TDenseSpace> BaseType;

    /// The explicit builder and solver definition
    typedef typename BaseType::ExplicitBuilderType ExplicitBuilderType;

    /// The local vector definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// Pointer definition of CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4
    KRATOS_CLASS_POINTER_DEFINITION(CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(SHOCK_CAPTURING);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
        , mShockCapturing(ThisParameters["shock_capturing"].GetBool())
    {
        // TODO: DO THE PARAMETERS CHECK
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, pExplicitBuilder, MoveMeshFlag, RebuildLevel)
    {
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : BaseType(rModelPart, MoveMeshFlag, RebuildLevel)
    {
    }

    /** Copy constructor.
     */
    CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4(const CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 &Other) = delete;

    /** Destructor.
     */
    virtual ~CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialization of member variables and prior operations
     * In this method we call the base strategy initialize and initialize the time derivatives
     * This is required to prevent OpenMP errors as the time derivatives are stored in the non-historical database
     */
    void Initialize() override
    {
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const unsigned int dim = r_process_info[DOMAIN_SIZE];

        // Call the base RK4 finalize substep method
        BaseType::Initialize();

        // Initialize the non-historical database values
        for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
            // Initialize the unknowns time derivatives to zero
            r_node.SetValue(DENSITY_TIME_DERIVATIVE, 0.0);
            r_node.SetValue(MOMENTUM_TIME_DERIVATIVE, MOMENTUM_TIME_DERIVATIVE.Zero());
            r_node.SetValue(TOTAL_ENERGY_TIME_DERIVATIVE, 0.0);
            // Initialize the shock capturing magnitudes
            // r_node.SetValue(SHOCK_SENSOR, 0.0);
            // r_node.SetValue(SHOCK_CAPTURING_VISCOSITY, 0.0);
            // r_node.SetValue(SHOCK_CAPTURING_CONDUCTIVITY, 0.0);
            // r_node.SetValue(VELOCITY, VELOCITY.Zero());
            // r_node.SetValue(CHARACTERISTIC_VELOCITY, CHARACTERISTIC_VELOCITY.Zero());
        }

        // If required, initialize the OSS projection variables
        if (r_process_info[OSS_SWITCH]) {
            for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
                r_node.SetValue(DENSITY_PROJECTION, 0.0);
                r_node.SetValue(TOTAL_ENERGY_PROJECTION, 0.0);
                r_node.SetValue(MOMENTUM_PROJECTION, ZeroVector(3));
            }
        }

        // If required, initialize the orthogonal projection shock capturing variables
        if (mShockCapturing) {
            // Initialize nodal values
            for (auto& r_node : r_model_part.GetCommunicator().LocalMesh().Nodes()) {
                r_node.SetValue(NODAL_AREA, 0.0);
                r_node.SetValue(MOMENTUM_GRADIENT, ZeroMatrix(dim,dim));
                r_node.SetValue(TOTAL_ENERGY_GRADIENT, ZeroVector(dim));
            }

            // Initialize elemental values
            for (auto& r_elem : r_model_part.GetCommunicator().LocalMesh().Elements()) {
                r_elem.SetValue(SHOCK_CAPTURING_VISCOSITY, 0.0);
                r_elem.SetValue(SHOCK_CAPTURING_CONDUCTIVITY, 0.0);
                r_elem.SetValue(MOMENTUM_GRADIENT, ZeroMatrix(dim,dim));
                r_elem.SetValue(TOTAL_ENERGY_GRADIENT, ZeroVector(dim));
            }
        }

        // // Requirements in case the shock capturing is needed
        // if (mShockCapturing) {
        //     InitializeShockCapturing();
        // }
    }

    /**
     * @brief Initialize the Runge-Kutta step
     * In case the mesh has been updated in the previous step we need to reinitialize the shock capturing
     * This includes the calculation of the nodal element size, nodal area and nodal neighbours
     */
    void InitializeSolutionStep() override
    {
        // Call the base RK4 initialize substep method
        BaseType::InitializeSolutionStep();

//         // Initialize acceleration variables
//         auto &r_model_part = BaseType::GetModelPart();
//         const int n_nodes = r_model_part.NumberOfNodes();
// #pragma omp parallel for
//         for (int i_node = 0; i_node < n_nodes; ++i_node) {
//             auto it_node = r_model_part.NodesBegin() + i_node;
//             it_node->GetValue(DENSITY_TIME_DERIVATIVE) = 0.0;
//             it_node->GetValue(MOMENTUM_TIME_DERIVATIVE) = ZeroVector(3);
//             it_node->GetValue(TOTAL_ENERGY_TIME_DERIVATIVE) = 0.0;
//         }

        // // If the mesh has changed, reinitialize the shock capturing method
        // if (BaseType::MoveMeshFlag() && mShockCapturing) {
        //     InitializeShockCapturing();
        // }

        // Calculate the magnitudes time derivatives
        UpdateUnknownsTimeDerivatives(1.0);

        // const auto& r_model_part = BaseType::GetModelPart();
        // const auto& r_process_info = r_model_part.GetProcessInfo();
        // if (r_process_info[OSS_SWITCH]) {
        //     CalculateOrthogonalSubScalesProjection();
        // }

        // // Perform orthogonal projection shock capturing
        // if (mShockCapturing) {
        //     CalculateOrthogonalProjectionShockCapturing();
        // }
    }

    /**
     * @brief Finalize the Runge-Kutta step
     * In this method we calculate the final linearised time derivatives after the final update
     * These will be the time derivatives employed in the first RK4 sub step of the next time step
     */
    void FinalizeSolutionStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeSolutionStep();

        // // Calculate the magnitudes time derivatives with the obtained solution
        // UpdateUnknownsTimeDerivatives(1.0);
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void SolveWithLumpedMassMatrix() override
    {
        // Call the base RK4 strategy to do the explicit update
        BaseType::SolveWithLumpedMassMatrix();

        // // If proceeds, do the intermediate PAD
        // if (mPAD) {
        //     // CLIPPING TEST
        //     // TODO: ADD OPENMP PRAGMAS
        //     const double aux_eps = 1.0e-6;
        //     for (auto& r_node : BaseType::GetModelPart().Nodes()) {
        //         double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
        //         auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
        //         double& r_enr = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);

        //         // Check that the density is above the minimum limit
        //         if (!r_node.IsFixed(DENSITY)) {
        //             if (r_rho < mMinDensity) {
        //                 // Update the residual vector accordingly
        //                 r_node.SetValue(DENSITY_PAD, true);
        //                 r_rho = mMinDensity;
        //             }
        //         }

        //         // Check that the energy is above the minimum limit
        //         if (!r_node.IsFixed(TOTAL_ENERGY)) {
        //             const double enr_lower_bound = inner_prod(r_mom, r_mom) / r_rho + aux_eps;
        //             if (r_enr < enr_lower_bound) {
        //                 r_enr = enr_lower_bound;
        //             }
        //         }
        //     }
        // }
    }

    void InitializeRungeKuttaIntermediateSubStep() override
    {
        // Call the base RK4 to perform the initialize intermediate RK sub step
        BaseType::InitializeRungeKuttaIntermediateSubStep();

        // Approximate the unknowns time derivatives with a FE scheme
        // These will be used in the next RK substep residual calculation to compute the subscales
        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        // const unsigned int rk_sub_step = r_process_info.GetValue(RUNGE_KUTTA_STEP);
        // if (rk_sub_step > 1) {
        //     const double sub_step_acc_coeff = 0.5;
        //     UpdateUnknownsTimeDerivatives(sub_step_acc_coeff);
        // }

        // Calculate the Orthogonal SubsScales projections
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }

        // Perform orthogonal projection shock capturing
        if (mShockCapturing) {
            CalculateOrthogonalProjectionShockCapturing();
        }

        // // Perform shock capturing
        // if (mShockCapturing) {
        //     // Call the shock detection process
        //     mpShockDetectionProcess->ExecuteInitializeSolutionStep();
        //     // Add the corresponding artificial magnitudes
        //     CalculateShockCapturingMagnitudes();
        // }
    }

    void InitializeRungeKuttaLastSubStep() override
    {
        // Call the base RK4 to perform the initialize intermediate RK sub step
        BaseType::InitializeRungeKuttaLastSubStep();

        // // Approximate the unknowns time derivatives with a FE scheme
        // // These will be used in the next RK substep residual calculation to compute the subscales
        // const double sub_step_acc_coeff = 1.0;
        // UpdateUnknownsTimeDerivatives(sub_step_acc_coeff);

        // Calculate the Orthogonal SubsScales projections
        auto& r_model_part = BaseType::GetModelPart();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        if (r_process_info[OSS_SWITCH]) {
            CalculateOrthogonalSubScalesProjection();
        }

        // Perform orthogonal projection shock capturing
        if (mShockCapturing) {
            CalculateOrthogonalProjectionShockCapturing();
        }

        // // Perform shock capturing
        // if (mShockCapturing) {
        //     // Call the shock detection process
        //     mpShockDetectionProcess->ExecuteInitializeSolutionStep();
        //     // Add the corresponding artificial magnitudes
        //     CalculateShockCapturingMagnitudes();
        // }
    }

    /**
     * @brief Finalize the Runge-Kutta intermediate substep
     * In this method we calculate the linearised time derivatives after the intemediate substep
     */
    void FinalizeRungeKuttaIntermediateSubStep() override
    {
        // Call the base RK4 finalize substep method
        BaseType::FinalizeRungeKuttaIntermediateSubStep();
    }

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    bool mShockCapturing;
    const bool mPAD = false;
    const double mMinDensity = 1.0e-2;
    const double mMinTotalEnergy = 1.0e-6;

    ShockDetectionProcess::UniquePointer mpShockDetectionProcess = nullptr;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Update the compressible Navier-Stokes unknowns time derivatives
     * This method approximates the compressible Navier-Stokes unknowns time derivatives
     * These are required to calculate the inertial stabilization terms in the compressible NS element
     * To that purpose a linear Forward-Euler interpolation is used
     */
    void UpdateUnknownsTimeDerivatives(const double SubStepAccCoefficient)
    {
        const double dt = BaseType::GetDeltaTime();
        KRATOS_ERROR_IF(dt < 1.0e-12) << "ProcessInfo DELTA_TIME is close to zero." << std::endl;
        auto &r_model_part = BaseType::GetModelPart();

#pragma omp parallel for
        for (int i_node = 0; i_node < static_cast<int>(r_model_part.NumberOfNodes()); ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;

            // Density DOF time derivative
            const double& r_rho = it_node->FastGetSolutionStepValue(DENSITY);
            const double& r_rho_old = it_node->FastGetSolutionStepValue(DENSITY, 1);
            it_node->GetValue(DENSITY_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_rho - r_rho_old) / dt;

            // Momentum DOF time derivative
            const auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
            const auto& r_mom_old = it_node->FastGetSolutionStepValue(MOMENTUM, 1);
            it_node->GetValue(MOMENTUM_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_mom - r_mom_old) / dt;

            // Total energy DOF time derivative
            const double& r_tot_enr = it_node->FastGetSolutionStepValue(TOTAL_ENERGY);
            const double& r_tot_enr_old = it_node->FastGetSolutionStepValue(TOTAL_ENERGY, 1);
            it_node->GetValue(TOTAL_ENERGY_TIME_DERIVATIVE) = SubStepAccCoefficient * (r_tot_enr - r_tot_enr_old) / dt;
        }
    }

    void CalculateOrthogonalSubScalesProjection()
    {
        // Get model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const int n_elem = r_model_part.NumberOfElements();
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const unsigned int block_size = r_process_info[DOMAIN_SIZE] + 2;

        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the projection values
#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->GetValue(NODAL_AREA) = 0.0;
            it_node->GetValue(DENSITY_PROJECTION) = 0.0;
            it_node->GetValue(MOMENTUM_PROJECTION) = ZeroVector(3);
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) = 0.0;
        }

        // Calculate the residuals projection
        double dens_proj;
        double tot_ener_proj;
        array_1d<double,3> mom_proj;
#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            // Calculate the projection values
            it_elem->Calculate(DENSITY_PROJECTION, dens_proj, r_process_info);
            it_elem->Calculate(MOMENTUM_PROJECTION, mom_proj, r_process_info);
            it_elem->Calculate(TOTAL_ENERGY_PROJECTION, tot_ener_proj, r_process_info);
            // Calculate the NODAL_AREA
            // TODO: This is not probably required each time
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();
            const double geom_domain_size = r_geom.DomainSize();
            const double aux_weight = geom_domain_size / static_cast<double>(n_nodes);
            for (auto& r_node : r_geom) {
#pragma omp atomic
                r_node.GetValue(NODAL_AREA) += aux_weight;
            }
        }

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            const double nodal_area = it_node->GetValue(NODAL_AREA);
            it_node->GetValue(DENSITY_PROJECTION) /= nodal_area;
            it_node->GetValue(MOMENTUM_PROJECTION) /= nodal_area;
            it_node->GetValue(TOTAL_ENERGY_PROJECTION) /= nodal_area;
        }
    }

    void CalculateOrthogonalProjectionShockCapturing()
    {
        // Calculate the model part data
        auto& r_model_part = BaseType::GetModelPart();
        const int n_nodes = r_model_part.NumberOfNodes();
        const int n_elems = r_model_part.NumberOfElements();
        const int dim = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
        // const unsigned int block_size = dim + 2;

        // // Get the required data from the explicit builder and solver
        // const auto p_explicit_bs = BaseType::pGetExplicitBuilder();
        // const auto& r_lumped_mass_vector = p_explicit_bs->GetLumpedMassMatrixVector();

        // Initialize the values to zero
#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            it_node->GetValue(NODAL_AREA) = 0.0;
            it_node->GetValue(MOMENTUM_GRADIENT) = ZeroMatrix(dim,dim);
            it_node->GetValue(TOTAL_ENERGY_GRADIENT) = ZeroVector(dim);
        }

        // Set the functor to calculate the element size
        // Note that this assumes a unique geometry in the computational mesh
        std::function<double(Geometry<Node<3>>& rGeometry)> avg_h_function;
        const GeometryData::KratosGeometryType geometry_type = (r_model_part.ElementsBegin()->GetGeometry()).GetGeometryType();
        switch (geometry_type) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                avg_h_function = [&](Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<2,3>::AverageElementSize(rGeometry);};
                break;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                avg_h_function = [&](Geometry<Node<3>>& rGeometry){return ElementSizeCalculator<3,4>::AverageElementSize(rGeometry);};
                break;
            default:
                KRATOS_ERROR << "Asking for a non-implemented geometry.";
        }

        // Loop the elements to project the gradients
        // Note that it is assumed that the gradient is constant within the element
        // Hence, only one Gauss point is used
        Geometry<Node<3>>::ShapeFunctionsGradientsType dNdX_container;
#pragma omp parallel for private(dNdX_container)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();
            const double geom_domain_size = r_geom.DomainSize();
            const double aux_weight = geom_domain_size / static_cast<double>(n_nodes);

            // Calculate the gradients in the center of the element
            auto& r_elem_mom_grad = it_elem->GetValue(MOMENTUM_GRADIENT);
            auto& r_elem_tot_ener_grad = it_elem->GetValue(TOTAL_ENERGY_GRADIENT);
            r_elem_mom_grad = ZeroMatrix(dim, dim);
            r_elem_tot_ener_grad = ZeroVector(dim);
            dNdX_container = r_geom.ShapeFunctionsIntegrationPointsGradients(dNdX_container, GeometryData::GI_GAUSS_1);
            const auto& r_dNdX = dNdX_container[0];

            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                auto& r_node = r_geom[i_node];
                const auto node_dNdX = row(r_dNdX, i_node);
                const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
                const double& r_tot_ener = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
                for (int d1 = 0; d1 < dim; ++d1) {
                    const double aux_val_ener = node_dNdX(d1) * r_tot_ener;
                    r_elem_tot_ener_grad[d1] += aux_val_ener;
                    for (int d2 = 0; d2 < dim; ++d2) {
                        const double aux_val_mom = node_dNdX(d1) * r_mom[d2];
                        r_elem_mom_grad(d1,d2) += aux_val_mom;
                    }
                }
            }

            // Project the computed gradients to the nodes
            const double midpoint_N = 1.0 / static_cast<double>(n_nodes);
            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                // Get nodal values
                auto& r_node = r_geom[i_node];
                auto& r_node_mom_grad = r_node.GetValue(MOMENTUM_GRADIENT);
                auto& r_node_tot_ener_grad = r_node.GetValue(TOTAL_ENERGY_GRADIENT);
                for (int d1 = 0; d1 < dim; ++d1) {
#pragma omp atomic
                    r_node_tot_ener_grad[d1] += midpoint_N * geom_domain_size * r_elem_tot_ener_grad[d1];
                    for (int d2 = 0; d2 < dim; ++d2) {
#pragma omp atomic
                        r_node_mom_grad(d1,d2) += midpoint_N * geom_domain_size * r_elem_mom_grad(d1,d2);
                    }
                }
#pragma omp atomic
                r_node.GetValue(NODAL_AREA) += aux_weight;
            }
        }

        #pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_model_part.NodesBegin() + i_node;
            // const double mass = r_lumped_mass_vector(i_node * block_size);
            // it_node->GetValue(MOMENTUM_GRADIENT) /= mass;
            // it_node->GetValue(TOTAL_ENERGY_GRADIENT) /= mass;
            const double weight = it_node->GetValue(NODAL_AREA);
            it_node->GetValue(MOMENTUM_GRADIENT) /= weight;
            it_node->GetValue(TOTAL_ENERGY_GRADIENT) /= weight;
        }

        // Calculate shock capturing values
        const double zero_tol = 1.0e-12;
        const double sc_visc_max_ratio = 1.0e3;
        const double sc_cond_max_ratio = 1.0e3;

        array_1d<double,3> midpoint_v;
        Matrix midpoint_mom_grad_proj;
        Vector midpoint_tot_ener_grad_proj;
        // Vector N_container, midpoint_tot_ener_grad_proj;
// #pragma omp parallel for private(N_container, midpoint_v, midpoint_mom_grad_proj, midpoint_tot_ener_grad_proj)
#pragma omp parallel for private(midpoint_v, midpoint_mom_grad_proj, midpoint_tot_ener_grad_proj)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            auto it_elem = r_model_part.ElementsBegin() + i_elem;
            auto& r_geom = it_elem->GetGeometry();
            const unsigned int n_nodes = r_geom.PointsNumber();

            // Interpolate the nodal projection values in the midpoint and calculate the average velocity norm
            midpoint_v = ZeroVector(3);
            midpoint_mom_grad_proj = ZeroMatrix(dim, dim);
            midpoint_tot_ener_grad_proj = ZeroVector(dim);
            const double midpoint_N = 1.0 / static_cast<double>(n_nodes);
            // r_geom.ShapeFunctionsValues(N_container, r_geom.Center());
            for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                const auto& r_node = r_geom[i_node];
                // Interpolate the nodal projection values in the midpoint
                const auto& r_node_mom_grad = r_node.GetValue(MOMENTUM_GRADIENT);
                const auto& r_node_tot_ener_grad = r_node.GetValue(TOTAL_ENERGY_GRADIENT);
                midpoint_mom_grad_proj += midpoint_N * r_node_mom_grad;
                midpoint_tot_ener_grad_proj += midpoint_N * r_node_tot_ener_grad;
                // midpoint_mom_grad_proj += N_container[i_node] * r_node_mom_grad;
                // midpoint_tot_ener_grad_proj += N_container[i_node] * r_node_tot_ener_grad;
                // Calculate the midpoint velocity
                const auto& r_mom = r_node.FastGetSolutionStepValue(MOMENTUM);
                const double& r_rho = r_node.FastGetSolutionStepValue(DENSITY);
                midpoint_v += midpoint_N * r_mom / r_rho;
                // midpoint_v += N_container[i_node] * r_mom / r_rho;
            }

            // Calculate the norms of the gradients
            // Total energy gradients
            const auto& r_tot_ener_elem_grad = it_elem->GetValue(TOTAL_ENERGY_GRADIENT);
            const double tot_ener_grad_norm = norm_2(r_tot_ener_elem_grad);
            const double tot_ener_grad_proj_norm = norm_2(r_tot_ener_elem_grad - midpoint_tot_ener_grad_proj);
            // Momentum gradients
            const auto& r_elem_mom_grad = it_elem->GetValue(MOMENTUM_GRADIENT);
            double mom_grad_norm = 0.0;
            double mom_grad_proj_norm = 0.0;
            for (unsigned int d1 = 0; d1 < dim; ++d1) {
                for (unsigned int d2 = 0; d2 < dim; ++d2) {
                    const double& r_aux = r_elem_mom_grad(d1,d2);
                    mom_grad_norm += std::pow(r_aux,2);
                    mom_grad_proj_norm += std::pow(r_aux - midpoint_mom_grad_proj(d1,d2), 2);
                }
            }
            mom_grad_norm = std::sqrt(mom_grad_norm);
            mom_grad_proj_norm = std::sqrt(mom_grad_proj_norm);

            // Calculate the shock capturing magnitudes
            const double c_a = 0.05;
            const double v_norm = norm_2(midpoint_v);
            const double aux = 0.5 * c_a * v_norm * avg_h_function(r_geom);
            it_elem->GetValue(SHOCK_CAPTURING_VISCOSITY) = mom_grad_norm > zero_tol && mom_grad_proj_norm / mom_grad_norm < sc_visc_max_ratio? aux * mom_grad_proj_norm / mom_grad_norm : 0.0;
            it_elem->GetValue(SHOCK_CAPTURING_CONDUCTIVITY) = tot_ener_grad_norm > zero_tol && tot_ener_grad_proj_norm / tot_ener_grad_norm < sc_cond_max_ratio ? aux * tot_ener_grad_proj_norm / tot_ener_grad_norm : 0.0;
        }
    }

    // /**
    //  * @brief Initializes the shock capturing method
    //  * This function performs all the operations required to initialize the shock capturing method
    //  * Note that if the mesh deforms (or changes the topology) it is required to repeat them again
    //  */
    // void InitializeShockCapturing()
    // {
    //     auto& r_model_part = BaseType::GetModelPart();

    //     // Calculate the nodal element size
    //     // This will be used in the calculation of the artificial shock capturing magnitudes
    //     auto nodal_h_process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(r_model_part);
    //     nodal_h_process.Execute();

    //     // Initialize the shock detection process
    //     mpShockDetectionProcess = Kratos::make_unique<ShockDetectionProcess>(r_model_part, DENSITY, DENSITY_GRADIENT);
    //     mpShockDetectionProcess->ExecuteInitialize();
    // }

//     /**
//      * @brief Calculates the artificial shock capturing values
//      * This method computes the values of the artificial shock capturing viscosity
//      * and conductivity (note that in here we assume dynamic artificial viscosity)
//      */
//     void CalculateShockCapturingMagnitudes()
//     {
//         // Get model part data
//         auto& r_model_part = BaseType::GetModelPart();
//         auto& r_comm = r_model_part.GetCommunicator();

//         // Calculate the corresponding artificial magnitudes
//         array_1d<double,3> aux_vel;
// #pragma omp parallel for private(aux_vel)
//         for (int i_node = 0; i_node < r_comm.LocalMesh().NumberOfNodes(); ++i_node) {
//             auto it_node = r_comm.LocalMesh().NodesBegin() + i_node;
//             const double& r_nodal_h = it_node->GetValue(NODAL_H);
//             const double& r_shock_sensor = it_node->GetValue(SHOCK_SENSOR);
//             const double& r_rho = it_node->FastGetSolutionStepValue(DENSITY);
//             const auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
//             aux_vel[0] = r_mom[0] / r_rho;
//             aux_vel[1] = r_mom[1] / r_rho;
//             aux_vel[2] = r_mom[2] / r_rho;
//             double& r_mu_sc = it_node->GetValue(SHOCK_CAPTURING_VISCOSITY);
//             r_mu_sc = r_shock_sensor * r_nodal_h * norm_2(aux_vel); // Kinematic shock capturing viscosity (this is the one implemented in the element right now)
//             // r_mu_sc = r_shock_sensor * r_nodal_h * norm_2(aux_vel) * r_rho; // Dynamic shock capturing viscosity
//         }

// LAPIDUS VISCOSITY
//         // If required recompute the NODAL_AREA
//         // This is required for the nodal gradients calculation
//         CalculateNodalAreaProcess<CalculateNodalAreaSettings::SaveAsNonHistoricalVariable>(
//         r_model_part,
//         r_model_part.GetProcessInfo().GetValue(DOMAIN_SIZE)).Execute();

//         array_1d<double,3> aux_vel;
// #pragma omp parallel for private(aux_vel)
//         for (int i_node = 0; i_node < r_comm.LocalMesh().NumberOfNodes(); ++i_node) {
//             auto it_node = r_comm.LocalMesh().NodesBegin() + i_node;
//             const double& r_rho = it_node->FastGetSolutionStepValue(DENSITY);
//             const auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
//             aux_vel[0] = r_mom[0] / r_rho;
//             aux_vel[1] = r_mom[1] / r_rho;
//             aux_vel[2] = r_mom[2] / r_rho;
//             it_node->GetValue(CHARACTERISTIC_VELOCITY) = norm_2(aux_vel);
//         }

//         // Calculate the shock variable nodal gradients
//         ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(
//             r_model_part,
//             CHARACTERISTIC_VELOCITY,
//             DENSITY_GRADIENT,
//             NODAL_AREA,
//             true).Execute();

// #pragma omp parallel for private(aux_vel)
//         for (int i_node = 0; i_node < r_comm.LocalMesh().NumberOfNodes(); ++i_node) {
//             auto it_node = r_comm.LocalMesh().NodesBegin() + i_node;
//             auto& r_l = it_node->GetValue(DENSITY_GRADIENT);
//             const double norm_l = norm_2(r_l);
//             if (norm_l > 1.0e-12) {
//                 r_l = r_l / norm_2(r_l);
//             } else {
//                 r_l = ZeroVector(3);
//             }
//             const double& r_rho = it_node->FastGetSolutionStepValue(DENSITY);
//             const auto& r_mom = it_node->FastGetSolutionStepValue(MOMENTUM);
//             aux_vel[0] = r_mom[0] / r_rho;
//             aux_vel[1] = r_mom[1] / r_rho;
//             aux_vel[2] = r_mom[2] / r_rho;
//             it_node->GetValue(CHARACTERISTIC_VELOCITY) = inner_prod(aux_vel, r_l);
//         }

//         // Calculate the shock variable nodal gradients
//         ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(
//             r_model_part,
//             CHARACTERISTIC_VELOCITY,
//             VELOCITY,
//             NODAL_AREA,
//             true).Execute();

//         const double c_lap = 1.0;
// #pragma omp parallel for firstprivate(c_lap)
//         for (int i_node = 0; i_node < r_comm.LocalMesh().NumberOfNodes(); ++i_node) {
//             auto it_node = r_comm.LocalMesh().NodesBegin() + i_node;
//             const double& r_nodal_h = it_node->GetValue(NODAL_H);
//             const auto& r_v_l = it_node->GetValue(VELOCITY);
//             const auto& r_l = it_node->GetValue(DENSITY_GRADIENT);
//             double& r_nu_sc = it_node->GetValue(SHOCK_CAPTURING_VISCOSITY);
//             r_nu_sc =  c_lap * std::pow(r_nodal_h,2) * std::abs(inner_prod(r_l, r_v_l));
//         }
//
//    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; /* Class CompressibleNavierStokesExplicitSolvingStrategyRungeKutta4 */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_COMPRESSIBLE_NAVIER_STOKES_EXPLICIT_SOLVING_STRATEGY_RUNGE_KUTTA_4  defined */
