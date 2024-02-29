//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//
//

#if !defined(KRATOS_PORO_EXPLICIT_VV_SCHEME_HPP_INCLUDED)
#define KRATOS_PORO_EXPLICIT_VV_SCHEME_HPP_INCLUDED

/* External includes */

/* Project includes */
#include "custom_strategies/schemes/one-phase_flow/poro_explicit_cd_scheme.hpp"
#include "utilities/variable_utils.h"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos {

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
 * @class PoroExplicitVVScheme
 * @ingroup StructuralMechanicsApplciation
 * @brief An explicit forward euler scheme with a split of the inertial term
 * @author Ignasi de Pouplana
 */
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class PoroExplicitVVScheme
    : public PoroExplicitCDScheme<TSparseSpace, TDenseSpace> {

public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef Scheme<TSparseSpace, TDenseSpace> BaseofBaseType;
    typedef PoroExplicitCDScheme<TSparseSpace, TDenseSpace> BaseType;

    /// Some definitions related with the base class
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition fo the node iterator
    typedef typename ModelPart::NodeIterator NodeIterator;

    /// The definition of the numerical limit
    static constexpr double numerical_limit = std::numeric_limits<double>::epsilon();

    using BaseType::mDeltaTime;
    using BaseType::mAlpha;
    using BaseType::mBeta;
    using BaseType::mTheta;
    using BaseType::mGCoefficient;

    /// Counted pointer of PoroExplicitVVScheme
    KRATOS_CLASS_POINTER_DEFINITION(PoroExplicitVVScheme);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details The PoroExplicitVVScheme method
     */
    PoroExplicitVVScheme()
        : PoroExplicitCDScheme<TSparseSpace, TDenseSpace>()
    {

    }

    /** Destructor.
    */
    virtual ~PoroExplicitVVScheme() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This method initializes some rutines related with the explicit scheme
     * @param rModelPart The model of the problem to solve
     * @param DomainSize The current dimention of the problem
     */
    void InitializeExplicitScheme(
        ModelPart& rModelPart,
        const SizeType DomainSize = 3
        ) override
    {
        KRATOS_TRY

        /// The array of ndoes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // The first iterator of the array of nodes
        const auto it_node_begin = rModelPart.NodesBegin();

        /// Initialise the database of the nodes
        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = (it_node_begin + i);
            it_node->SetValue(NODAL_MASS, 0.0);
            // TODO: Set Nodal AntiCompressibility to zero for mass-balance equation (C=1/Q, with Q being the compressibility coeff.)
            array_1d<double, 3>& r_force_residual = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            double& r_flux_residual = it_node->FastGetSolutionStepValue(FLUX_RESIDUAL);
            array_1d<double, 3>& r_external_force = it_node->FastGetSolutionStepValue(EXTERNAL_FORCE);
            array_1d<double, 3>& r_internal_force = it_node->FastGetSolutionStepValue(INTERNAL_FORCE);
            array_1d<double, 3>& r_damping_force = it_node->FastGetSolutionStepValue(DAMPING_FORCE);
            noalias(r_force_residual) = ZeroVector(3);
            r_flux_residual = 0.0;
            noalias(r_external_force) = ZeroVector(3);
            noalias(r_internal_force) = ZeroVector(3);
            noalias(r_damping_force) = ZeroVector(3);
            Matrix& rInitialStress = it_node->FastGetSolutionStepValue(INITIAL_STRESS_TENSOR);
            if(rInitialStress.size1() != 3)
                rInitialStress.resize(3,3,false);
            noalias(rInitialStress) = ZeroMatrix(3,3);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method initializes the residual in the nodes of the model part
     * @param rModelPart The model of the problem to solve
     */
    void InitializeResidual(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        // Auxiliar values
        const array_1d<double, 3> zero_array = ZeroVector(3);
        // Initializing the variables
        VariableUtils().SetVariable(FORCE_RESIDUAL, zero_array,r_nodes);
        VariableUtils().SetVariable(FLUX_RESIDUAL, 0.0,r_nodes);
        VariableUtils().SetVariable(EXTERNAL_FORCE, zero_array,r_nodes);
        VariableUtils().SetVariable(INTERNAL_FORCE, zero_array,r_nodes);
        VariableUtils().SetVariable(DAMPING_FORCE, zero_array,r_nodes);

        KRATOS_CATCH("")
    }

     void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY;

        this->CalculateAndAddRHS(rModelPart);

        // The current process info
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // The array of nodes
        NodesArrayType& r_nodes = rModelPart.Nodes();

        /// Working in 2D/3D (the definition of DOMAIN_SIZE is check in the Check method)
        const SizeType dim = r_current_process_info[DOMAIN_SIZE];

        // Step Update
        mDeltaTime = r_current_process_info[DELTA_TIME];

        // The iterator of the first node
        const auto it_node_begin = rModelPart.NodesBegin();

        // Getting dof position
        const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);

        #pragma omp parallel for schedule(guided,512)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            // Current step information "N+1" (before step update).
            this->PredictTranslationalDegreesOfFreedom(it_node_begin + i, disppos, dim);
        } // for Node parallel

        InitializeResidual(rModelPart);

        this->CalculateAndAddRHS(rModelPart);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    virtual void PredictTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        )
    {
        array_1d<double, 3>& r_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3>& r_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);

        const array_1d<double, 3>& r_external_force = itCurrentNode->FastGetSolutionStepValue(EXTERNAL_FORCE);
        const array_1d<double, 3>& r_internal_force = itCurrentNode->FastGetSolutionStepValue(INTERNAL_FORCE);
        const array_1d<double, 3>& r_damping_force = itCurrentNode->FastGetSolutionStepValue(DAMPING_FORCE);

        std::array<bool, 3> fix_displacements = {false, false, false};
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        // Solution of the explicit equation:
        if (nodal_mass > numerical_limit){
            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j] == false) {
                    r_displacement[j] += r_velocity[j]*mDeltaTime + 0.5 * (r_external_force[j]
                                                                            - r_internal_force[j]
                                                                            - r_damping_force[j])/nodal_mass * mDeltaTime * mDeltaTime;
                    r_velocity[j] += 0.5 * mDeltaTime * (r_external_force[j]
                                                        - r_internal_force[j]
                                                        - r_damping_force[j])/nodal_mass;
                }
            }
        }
        else {
            noalias(r_displacement) = ZeroVector(3);
            noalias(r_velocity) = ZeroVector(3);
        }
    }

    /**
     * @brief This method updates the translation DoF
     * @param itCurrentNode The iterator of the current node
     * @param DisplacementPosition The position of the displacement dof on the database
     * @param DomainSize The current dimention of the problem
     */
    void UpdateTranslationalDegreesOfFreedom(
        NodeIterator itCurrentNode,
        const IndexType DisplacementPosition,
        const SizeType DomainSize = 3
        ) override
    {
        array_1d<double, 3>& r_velocity = itCurrentNode->FastGetSolutionStepValue(VELOCITY);
        const double nodal_mass = itCurrentNode->GetValue(NODAL_MASS);

        double& r_current_liquid_pressure = itCurrentNode->FastGetSolutionStepValue(LIQUID_PRESSURE);
        double& r_current_dt_liquid_pressure = itCurrentNode->FastGetSolutionStepValue(DT_LIQUID_PRESSURE);

        const array_1d<double, 3>& r_external_force = itCurrentNode->FastGetSolutionStepValue(EXTERNAL_FORCE);
        const array_1d<double, 3>& r_internal_force = itCurrentNode->FastGetSolutionStepValue(INTERNAL_FORCE);
        const array_1d<double, 3>& r_damping_force = itCurrentNode->FastGetSolutionStepValue(DAMPING_FORCE);

        std::array<bool, 3> fix_displacements = {false, false, false};
        fix_displacements[0] = (itCurrentNode->GetDof(DISPLACEMENT_X, DisplacementPosition).IsFixed());
        fix_displacements[1] = (itCurrentNode->GetDof(DISPLACEMENT_Y, DisplacementPosition + 1).IsFixed());
        if (DomainSize == 3)
            fix_displacements[2] = (itCurrentNode->GetDof(DISPLACEMENT_Z, DisplacementPosition + 2).IsFixed());

        // Solution of the explicit equation:
        if (nodal_mass > numerical_limit){
            for (IndexType j = 0; j < DomainSize; j++) {
                if (fix_displacements[j] == false) {
                    r_velocity[j] += 0.5 * mDeltaTime * (r_external_force[j]
                                                            - r_internal_force[j]
                                                            - r_damping_force[j])/nodal_mass;
                }
            }
        }
        else {
            noalias(r_velocity) = ZeroVector(3);
        }
        // Solution of the darcy_equation
        if( itCurrentNode->IsFixed(LIQUID_PRESSURE) == false ) {
            // TODO: this is on standby
            r_current_liquid_pressure = 0.0;
            r_current_dt_liquid_pressure = 0.0;
        }

        const array_1d<double, 3>& r_velocity_old = itCurrentNode->FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double, 3>& r_acceleration = itCurrentNode->FastGetSolutionStepValue(ACCELERATION);

        noalias(r_acceleration) = (1.0/mDeltaTime) * (r_velocity - r_velocity_old);
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        // this->TCalculateRHSContribution(rCurrentElement, RHS_Contribution, rCurrentProcessInfo);
        // rCurrentElement.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
        rCurrentElement.AddExplicitContribution(RHS_Contribution, RESIDUAL_VECTOR, DAMPING_FORCE, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@}
    ///@name Protected Structs
    ///@{


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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

}; /* Class PoroExplicitVVScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_PORO_EXPLICIT_VV_SCHEME_HPP_INCLUDED  defined */
