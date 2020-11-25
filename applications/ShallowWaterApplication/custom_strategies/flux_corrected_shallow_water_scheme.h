//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_FLUX_CORRECTED_SHALLOW_WATER_SCHEME_H_INCLUDED
#define KRATOS_FLUX_CORRECTED_SHALLOW_WATER_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "shallow_water_residual_based_bdf_scheme.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"

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

/**
 * @class FluxCorrectedShallowWaterScheme
 * @ingroup KratosShallowWaterApplication
 * @brief BDF integration scheme (for dynamic problems) with flux correction for extra diffusion to ensure monotonic solutions
 * @details The \f$n\f$ order Backward Differentiation Formula (BDF) method is a two step \f$n\f$ order accurate method.
 * This scheme is designed to solve a system of the type:
 * \f[
 *   \mathbf{M} \frac{du_{n0}}{dt} + \mathbf{K} u_{n0} = \mathbf{f}_{ext}
 * \f]
 * @author Miguel Maso Sotomayor
 */
template<class TSparseSpace,  class TDenseSpace>
class FluxCorrectedShallowWaterScheme
    : public ShallowWaterResidualBasedBDFScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( FluxCorrectedShallowWaterScheme );

    typedef ShallowWaterResidualBasedBDFScheme<TSparseSpace,TDenseSpace> SWBaseType;

    typedef typename SWBaseType::BDFBaseType                            BDFBaseType;

    typedef typename BDFBaseType::ImplicitBaseType                 ImplicitBaseType;

    typedef typename ImplicitBaseType::BaseType                            BaseType;
  
    typedef typename BaseType::Pointer                              BaseTypePointer;

    typedef typename BaseType::DofsArrayType                          DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                  TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                  TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType          LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType          LocalSystemMatrixType;

    typedef typename ModelPart::NodeType                                   NodeType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor
    explicit FluxCorrectedShallowWaterScheme(const std::size_t Order = 2)
        : SWBaseType(Order)
    {}

    // Copy Constructor
    explicit FluxCorrectedShallowWaterScheme(FluxCorrectedShallowWaterScheme& rOther)
        : SWBaseType(rOther)
    {}

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new FluxCorrectedShallowWaterScheme(*this) );
    }

    // Destructor
    ~FluxCorrectedShallowWaterScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This is the place to initialize the Scheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void Initialize(ModelPart& rModelPart) override
    {
        // Memory allocation
        const IndexType num_threads = OpenMPUtils::GetNumThreads();
        mMl.resize(num_threads);
        mUn0.resize(num_threads);

        // Initialization of non-historical variables
        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            r_node.SetValue(POSITIVE_FLUX, 0.0);
            r_node.SetValue(NEGATIVE_FLUX, 0.0);
            r_node.SetValue(MAX_INCREMENT, 0.0);
            r_node.SetValue(MIN_DECREMENT, 0.0);
            r_node.SetValue(POSITIVE_RATIO, 0.0);
            r_node.SetValue(NEGATIVE_RATIO, 0.0);
        });

        // Initialization and execution of nodal neighbours search
        const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
        FindGlobalNodalNeighboursProcess(r_data_communicator, rModelPart).Execute();

        // Finalization of initialize
        SWBaseType::Initialize(rModelPart);
    }

    /**
     * @brief It initializes a non-linear iteration (for the element)
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        SWBaseType::InitializeNonLinIteration(rModelPart, rA, rDx, rb);

        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            r_node.GetValue(POSITIVE_FLUX) = 0.0;
            r_node.GetValue(NEGATIVE_FLUX) = 0.0;
            r_node.GetValue(MAX_INCREMENT) = 0.0;
            r_node.GetValue(MIN_DECREMENT) = 0.0;
            r_node.GetValue(POSITIVE_RATIO) = 0.0;
            r_node.GetValue(NEGATIVE_RATIO) = 0.0;
        });

        const ProcessInfo& r_const_process_info = rModelPart.GetProcessInfo();
        block_for_each(rModelPart.Elements(), [&](Element& r_elem){
            ComputeAntiFluxes(r_elem, r_const_process_info);
        });

        block_for_each(rModelPart.Nodes(), [&](NodeType& r_node){
            const auto& neigh_nodes = r_node.GetValue(NEIGHBOUR_NODES);
            double u_i = r_node.FastGetSolutionStepValue(HEIGHT);
            double& max_incr = r_node.GetValue(MAX_INCREMENT);
            double& min_decr = r_node.GetValue(MIN_DECREMENT);
            for (size_t j = 0; j < neigh_nodes.size(); ++j){
                double delta_ij = neigh_nodes[j].FastGetSolutionStepValue(HEIGHT) - u_i;
                max_incr = std::max(max_incr, delta_ij);
                min_decr = std::min(min_decr, delta_ij);
            }
            double pos_flux = r_node.GetValue(POSITIVE_FLUX);
            double neg_flux = r_node.GetValue(NEGATIVE_FLUX);
            if (pos_flux > 0.0){
                r_node.GetValue(POSITIVE_RATIO) = std::min(1.0, max_incr / pos_flux);
            } else {
                r_node.GetValue(POSITIVE_RATIO) = 1;
            }
            if (neg_flux < 0.0){
                r_node.GetValue(NEGATIVE_RATIO) = std::min(1.0, min_decr / neg_flux);
            } else {
                r_node.GetValue(NEGATIVE_RATIO) = 1;
            }
        });
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentElement The element to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(ImplicitBaseType::mMatrix.D[this_thread], rCurrentProcessInfo);

        ComputeLumpedMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], mMl[this_thread]);

        AddDynamicsToLHS(rLHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], mMl[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(rCurrentElement, rRHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], mMl[this_thread], rCurrentProcessInfo);
    
        AddFluxCorrection<Element>(
            rCurrentElement,
            rLHS_Contribution,
            rRHS_Contribution,
            ImplicitBaseType::mMatrix.M[this_thread],
            mMl[this_thread],
            ImplicitBaseType::mMatrix.D[this_thread],
            mUn0[this_thread],
            BDFBaseType::mVector.dotun0[this_thread]);

        SWBaseType::mRotationTool.Rotate(rLHS_Contribution, rRHS_Contribution, rCurrentElement.GetGeometry());
        SWBaseType::mRotationTool.ApplySlipCondition(rLHS_Contribution, rRHS_Contribution, rCurrentElement.GetGeometry());

        KRATOS_CATCH("FluxCorrectedShallowWaterScheme.CalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rCurrentElement.CalculateRightHandSide(rRHS_Contribution,rCurrentProcessInfo);

        rCurrentElement.CalculateMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        rCurrentElement.CalculateDampingMatrix(ImplicitBaseType::mMatrix.D[this_thread],rCurrentProcessInfo);

        rCurrentElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

        ComputeLumpedMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], mMl[this_thread]);

        AddDynamicsToRHS(rCurrentElement, rRHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        AddFluxCorrection<Element>(
            rCurrentElement,
            rRHS_Contribution,
            ImplicitBaseType::mMatrix.M[this_thread],
            mMl[this_thread],
            ImplicitBaseType::mMatrix.D[this_thread],
            mUn0[this_thread],
            BDFBaseType::mVector.dotun0[this_thread]);

        SWBaseType::mRotationTool.Rotate(rRHS_Contribution, rCurrentElement.GetGeometry());
        SWBaseType::mRotationTool.ApplySlipCondition(rRHS_Contribution, rCurrentElement.GetGeometry());

        KRATOS_CATCH("FluxCorrectedShallowWaterScheme.Calculate_RHS_Contribution");
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rCurrentCondition The condition to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateLocalSystem(rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(ImplicitBaseType::mMatrix.D[this_thread], rCurrentProcessInfo);

        ComputeLumpedMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], mMl[this_thread]);

        AddDynamicsToLHS(rLHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        AddFluxCorrection<Condition>(
            rCurrentCondition,
            rLHS_Contribution,
            rRHS_Contribution,
            ImplicitBaseType::mMatrix.M[this_thread],
            mMl[this_thread],
            ImplicitBaseType::mMatrix.D[this_thread],
            mUn0[this_thread],
            BDFBaseType::mVector.dotun0[this_thread]);

        SWBaseType::mRotationTool.Rotate(rLHS_Contribution, rRHS_Contribution, rCurrentCondition.GetGeometry());
        SWBaseType::mRotationTool.ApplySlipCondition(rLHS_Contribution, rRHS_Contribution, rCurrentCondition.GetGeometry());

        KRATOS_CATCH("FluxCorrectedShallowWaterScheme.CalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentCondition The condition to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rCurrentCondition.CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);

        rCurrentCondition.CalculateMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        rCurrentCondition.CalculateDampingMatrix(ImplicitBaseType::mMatrix.D[this_thread], rCurrentProcessInfo);

        ComputeLumpedMassMatrix(ImplicitBaseType::mMatrix.M[this_thread], mMl[this_thread]);

        AddDynamicsToRHS(rCurrentCondition, rRHS_Contribution, ImplicitBaseType::mMatrix.D[this_thread], ImplicitBaseType::mMatrix.M[this_thread], rCurrentProcessInfo);

        AddFluxCorrection<Condition>(
            rCurrentCondition,
            rRHS_Contribution,
            ImplicitBaseType::mMatrix.M[this_thread],
            mMl[this_thread],
            ImplicitBaseType::mMatrix.D[this_thread],
            mUn0[this_thread],
            BDFBaseType::mVector.dotun0[this_thread]);

        SWBaseType::mRotationTool.Rotate(rRHS_Contribution, rCurrentCondition.GetGeometry());
        SWBaseType::mRotationTool.ApplySlipCondition(rRHS_Contribution, rCurrentCondition.GetGeometry());

        KRATOS_CATCH("FluxCorrectedShallowWaterScheme.Calculate_RHS_Contribution");
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
    std::string Info() const override
    {
        return "FluxCorrectedShallowWaterScheme";
    }

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    std::vector<Vector> mUn0;
    std::vector<Matrix> mMl;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief It adds the dynamic LHS contribution of the elements
     * \f[ LHS = \frac{d(-RHS)}{d(u_{n0})} = c_0^2\mathbf{M} + c_0 \mathbf{D} + \mathbf{K} \f]
     * @param rLHS_Contribution The dynamic contribution for the LHS
     * @param rD The diffusion matrix
     * @param rM The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // Adding mass contribution to the dynamic stiffness
        if (rM.size1() != 0) { // if M matrix declared
            noalias(rLHS_Contribution) += rM * BDFBaseType::mBDF[0];
        }

        // Adding monotonic diffusion
        if (rD.size1() != 0) { // if D matrix declared
            noalias(rLHS_Contribution) += rD;
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements
     * @param rElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The diffusion matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const auto& r_const_element = rElement;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            r_const_element.GetFirstDerivativesVector(BDFBaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BDFBaseType::mVector.dotun0[this_thread]);
        }

        // Adding monotonic diffusion
        if (rD.size1() != 0) {
            r_const_element.GetValuesVector(mUn0[this_thread]);
            noalias(rRHS_Contribution) -= prod(rD, mUn0[this_thread]);
        }
    }

    /**
     * @brief It adds the dynamic RHS contribution of the elements
     * @param rElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The diffusion matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    void AddDynamicsToRHS(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        LocalSystemMatrixType& rD,
        LocalSystemMatrixType& rM,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        const auto& r_const_condition = rCondition;
        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (rM.size1() != 0) {
            r_const_condition.GetFirstDerivativesVector(BDFBaseType::mVector.dotun0[this_thread], 0);
            noalias(rRHS_Contribution) -= prod(rM, BDFBaseType::mVector.dotun0[this_thread]);
        }

        // Adding monotonic diffusion
        if (rD.size1() != 0) {
            r_const_condition.GetValuesVector(mUn0[this_thread]);
            noalias(rRHS_Contribution) -= prod(rD, mUn0[this_thread]);
        }
    }

    template<class EntityType>
    void ComputeAntiFluxes(EntityType& rEntity, const ProcessInfo& rProcessInfo)
    {
        const auto& r_const_entity = rEntity; // TODO: remove that statement as soon as deprecation warnings are removed
        const IndexType t = OpenMPUtils::ThisThread();

        rEntity.CalculateMassMatrix(ImplicitBaseType::mMatrix.M[t], rProcessInfo);

        rEntity.CalculateDampingMatrix(ImplicitBaseType::mMatrix.D[t], rProcessInfo);

        r_const_entity.GetValuesVector(mUn0[t]);

        r_const_entity.GetFirstDerivativesVector(BDFBaseType::mVector.dotun0[t]);

        ComputeLumpedMassMatrix(ImplicitBaseType::mMatrix.M[t], mMl[t]);

        auto aec = prod(ImplicitBaseType::mMatrix.D[t], mUn0[t]) + prod(mMl[t] - ImplicitBaseType::mMatrix.M[t], BDFBaseType::mVector.dotun0[t]);

        auto r_geom = rEntity.GetGeometry();
        const IndexType block_size = aec.size() / r_geom.size();
        for (IndexType i = 0; i < r_geom.size(); ++i)
        {
            const double flux = aec((i+1)*block_size -1);
            if (flux > 0) {
                r_geom[i].SetLock();
                r_geom[i].GetValue(POSITIVE_FLUX) += flux;
                r_geom[i].UnSetLock();
            } else {
                r_geom[i].SetLock();
                r_geom[i].GetValue(NEGATIVE_FLUX) += flux;
                r_geom[i].UnSetLock();
            }
        }
    }

    template<class EntityType>
    void AddFluxCorrection(
        EntityType& rEntity,
        LocalSystemVectorType& rRHS,
        const LocalSystemMatrixType& rMc,
        const LocalSystemMatrixType& rMl,
        const LocalSystemMatrixType& rD,
        const LocalSystemVectorType& rU,
        const LocalSystemVectorType& rDotU)
    {
        // Construction of the element contribution of anti-fluxes
        auto aec = prod(rD, rU) + prod(rMl - rMc, rDotU);

        // Checking the sign of the contribution
        IndexType block_size = 3;
        IndexType nodes = rU.size() / block_size;
        double element_contribution = 0.0;
        for (IndexType i = 0; i < nodes; ++i)
        {
            element_contribution += aec((i+1)*block_size -1);
        }

        // Getting the limiter
        double c = 1.0;
        if (element_contribution > 0.0) {
            for (auto& r_node : rEntity.GetGeometry())
            {
                c = std::min(c, r_node.GetValue(POSITIVE_RATIO));
            }
        } else {
            for (auto& r_node : rEntity.GetGeometry())
            {
                c = std::min(c, r_node.GetValue(NEGATIVE_RATIO));
            }   
        }

        // Adding the limited anti-diffusion
        rRHS += c * aec;
    }

    template<class EntityType>
    void AddFluxCorrection(
        EntityType& rEntity,
        LocalSystemMatrixType& rLHS,
        LocalSystemVectorType& rRHS,
        const LocalSystemMatrixType& rMc,
        const LocalSystemMatrixType& rMl,
        const LocalSystemMatrixType& rD,
        const LocalSystemVectorType& rU,
        const LocalSystemVectorType& rDotU)
    {
        // Construction of the element contribution of anti-fluxes
        auto aec = prod(rD, rU) + prod(rMl - rMc, rDotU);

        // Checking the sign of the contribution
        IndexType block_size = 3;
        IndexType nodes = rU.size() / block_size;
        double element_contribution = 0.0;
        for (IndexType i = 0; i < nodes; ++i)
        {
            element_contribution += aec((i+1)*block_size -1);
        }

        // Getting the limiter
        double c = 1.0;
        if (element_contribution > 0.0) {
            for (auto& r_node : rEntity.GetGeometry())
            {
                c = std::min(c, r_node.GetValue(POSITIVE_RATIO));
            }
        } else {
            for (auto& r_node : rEntity.GetGeometry())
            {
                c = std::min(c, r_node.GetValue(NEGATIVE_RATIO));
            }   
        }

        // Adding the limited anti-diffusion
        // rLHS += c * (BDFBaseType::mBDF[0]*(rMc - rMl) - rD);
        rRHS += c * aec;
    }

    void ComputeLumpedMassMatrix(
        const LocalSystemMatrixType& rConsistentMassMatrix,
        LocalSystemMatrixType& rLumpedMassMatrix)
    {
        const IndexType size = rConsistentMassMatrix.size1();

        if (rLumpedMassMatrix.size1() != size) {
            rLumpedMassMatrix.resize(size,size,false);
        }

        for (IndexType i = 0; i < size; ++i)
        {
            double l = 0.0;
            for (IndexType j = 0; j < size; ++j)
            {
                l += rConsistentMassMatrix(i,j);
                rLumpedMassMatrix(i,j) = 0.0;
            }
            rLumpedMassMatrix(i,i) = l;
        }
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
    ///@{

}; // Class FluxCorrectedShallowWaterScheme

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // Namespace Kratos

#endif // KRATOS_FLUX_CORRECTED_SHALLOW_WATER_SCHEME_H_INCLUDED defined
