//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $StructuralMechanicsApplication $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         September 2016   $
//   Revision:            $Revision:                0.0   $

#if !defined(KRATOS_EIGENSOLVER_DYNAMIC_SCHEME )
#define  KRATOS_EIGENSOLVER_DYNAMIC_SCHEME


// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h" //TODO: For the OSS_SWITCH. Do this better...
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "structural_mechanics_application_variables.h"

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

/// An adapter scheme for obtaining mass and stiffness matrices for dynamic eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace
         >
class EigensolverDynamicScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( EigensolverDynamicScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EigensolverDynamicScheme() : Scheme<TSparseSpace,TDenseSpace>() {}

    /// Destructor.
    ~EigensolverDynamicScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // TODO: We should place the OSS projections somewhere else
    void Initialize(ModelPart& rModelPart) override
    {
        // Allocate the OSS projection variables
        const auto &r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;
        if (oss_switch) {
            const array_1d<double,3> aux_zero = ZeroVector(3);
            block_for_each(rModelPart.Nodes(), aux_zero, [](Node<3>& rNode, array_1d<double,3>& rAuxZero){
                rNode.SetValue(DISPLACEMENT_PROJECTION, rAuxZero);
                rNode.SetValue(VOLUMETRIC_STRAIN_PROJECTION, 0.0);
            });
        }

        // Call the base Initialize method
        BaseType::Initialize(rModelPart);
    }

    //TODO: We should place the OSS projections somewhere else
    //TODO: To make this fully flexible, I'd promote the parameters based constructor with a list of variables to which their projection is to be computed
    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY

        // Check if the Orthogonal SubScales (OSS) are active
        const auto& r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;

        // Calculate the OSS projections
        if (oss_switch) {
            // Initialize the projection values
            block_for_each(rModelPart.Nodes(), [](Node<3>& rNode){
                rNode.GetValue(DISPLACEMENT_PROJECTION) = ZeroVector(3);
                rNode.GetValue(VOLUMETRIC_STRAIN_PROJECTION) = 0.0;
            });

            // Calculate the element residuals projection
            std::tuple<double, array_1d<double,3>> oss_proj_tls;
            block_for_each(rModelPart.Elements(), oss_proj_tls, [&](Element& rElement, std::tuple<double, array_1d<double,3>>& rOssProjTLS){
                double& r_eps_proj = std::get<0>(rOssProjTLS);
                array_1d<double,3>& r_disp_proj = std::get<1>(rOssProjTLS);
                rElement.Calculate(DISPLACEMENT_PROJECTION, r_disp_proj, r_process_info);
                rElement.Calculate(VOLUMETRIC_STRAIN_PROJECTION, r_eps_proj, r_process_info);
            });

            // Do the nodal weighting
            //TODO: We need to do the weighted L2 projection with the density for the multimaterial case
            block_for_each(rModelPart.Nodes(), [](Node<3>& rNode){
                const double nodal_area = rNode.GetValue(NODAL_AREA);
                rNode.GetValue(DISPLACEMENT_PROJECTION) /= nodal_area;
                rNode.GetValue(VOLUMETRIC_STRAIN_PROJECTION) /= nodal_area;
            });
        }

        // Call base class FinalizeNonLinIteration
        BaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo
    ) override
    {
        KRATOS_TRY

        if (CurrentProcessInfo[BUILD_LEVEL] == 1)
        { // mass matrix
            rCurrentElement.CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
                RHS_Contribution.resize(LocalSize,false);
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if (CurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
        {
            rCurrentElement.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        CalculateSystemContributions(
                rCurrentElement,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentElement.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        rCurrentElement.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        if (CurrentProcessInfo[BUILD_LEVEL] == 1)
        { // mass matrix
            rCurrentCondition.CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
            std::size_t LocalSize = LHS_Contribution.size1();
            if (RHS_Contribution.size() != LocalSize)
            {
                RHS_Contribution.resize(LocalSize,false);
            }
            noalias(RHS_Contribution) = ZeroVector(LocalSize);
        }
        else if (CurrentProcessInfo[BUILD_LEVEL] == 2) // stiffness matrix
        {
            rCurrentCondition.CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL" << std::endl;
        }

        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;
        RHS_Contribution.resize(LHS_Contribution.size1(), false);
        CalculateSystemContributions(
                rCurrentCondition,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        rCurrentCondition.CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        rCurrentCondition.EquationIdVector(EquationId,CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    ///@}

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

}  /* namespace Kratos.*/

#endif /* KRATOS_EIGENSOLVER_DYNAMIC_SCHEME  defined */

