//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_RANS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_WALL_CONDITION_H_INCLUDED)
#define KRATOS_RANS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_WALL_CONDITION_H_INCLUDED

// System includes
#include <algorithm>
#include <iterator>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/condition.h"
#include "includes/properties.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

#include "custom_conditions/evm_k_epsilon/rans_evm_epsilon_adjoint_wall_condition.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_k_adjoint_wall_condition.h"
#include "custom_conditions/evm_k_epsilon/rans_evm_vms_monolithic_adjoint_wall_condition.h"

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

template <unsigned int TDim>
class RansEvmMonolithicKEpsilonVMSAdjointWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    // defining the base type
    using BaseType = Condition;
    // defining the base adjoint base fluid condition type
    using AdjointFluidCondition = RansEvmVmsMonolithicAdjointWallCondition<TDim>;
    // defining the k condition type
    using AdjointKCondition = RansEvmKAdjointWallCondition<TDim>;
    // defining the epsilon condition type
    using AdjointEpsilonCondition = RansEvmEpsilonAdjointWallCondition<TDim>;

    constexpr static unsigned int TNumNodes = TDim;

    constexpr static unsigned int TFluidBlockSize = (TDim + 1);

    constexpr static unsigned int TFluidLocalSize = TFluidBlockSize * TNumNodes;

    constexpr static unsigned int TKBlockSize = 1;

    constexpr static unsigned int TKLocalSize = TKBlockSize * TNumNodes;

    constexpr static unsigned int TEpsilonBlockSize = 1;

    constexpr static unsigned int TEpsilonLocalSize = TEpsilonBlockSize * TNumNodes;

    constexpr static unsigned int TConditionBlockSize =
        (TFluidBlockSize + TKBlockSize + TEpsilonBlockSize);

    constexpr static unsigned int TConditionLocalSize = TConditionBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    // variable definitions
    typedef std::size_t IndexType;

    typedef Condition::NodeType NodeType;

    typedef Condition::NodesArrayType NodesArrayType;

    typedef Condition::GeometryType GeometryType;

    typedef Condition::PropertiesType PropertiesType;

    typedef Condition::VectorType VectorType;

    typedef Condition::MatrixType MatrixType;

    ///@name Pointer Definitions
    /// Pointer definition of RansEvmMonolithicKEpsilonVMSAdjointWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(RansEvmMonolithicKEpsilonVMSAdjointWallCondition);

    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit RansEvmMonolithicKEpsilonVMSAdjointWallCondition(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    RansEvmMonolithicKEpsilonVMSAdjointWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    RansEvmMonolithicKEpsilonVMSAdjointWallCondition(IndexType NewId,
                                                     GeometryType::Pointer pGeometry,
                                                     PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~RansEvmMonolithicKEpsilonVMSAdjointWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjointWallCondition>(
            NewId, Condition::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjointWallCondition>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new condition pointer and clones the previous condition data
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<RansEvmMonolithicKEpsilonVMSAdjointWallCondition>(
            NewId, Condition::GetGeometry().Create(ThisNodes), Condition::pGetProperties());
        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        BaseType::Check(rCurrentProcessInfo);

        AdjointFluidCondition adjoint_fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition adjoint_k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition adjoint_epsilon_condition(this->Id(), this->pGetGeometry());

        adjoint_fluid_condition.Check(rCurrentProcessInfo);
        adjoint_k_condition.Check(rCurrentProcessInfo);
        adjoint_epsilon_condition.Check(rCurrentProcessInfo);

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
        {
            NodeType& r_node = this->GetGeometry()[i_node];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_2, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_VECTOR_3, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_FLUID_SCALAR_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(AUX_ADJOINT_FLUID_VECTOR_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_2, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_1_ADJOINT_3, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUX_ADJOINT_SCALAR_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_1, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_2, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_SCALAR_2_ADJOINT_3, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUX_ADJOINT_SCALAR_2, r_node);
        }

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * this determines the conditional equation ID vector for all conditional
     * DOFs
     * @param rResult the conditional equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rResult.size() != TConditionLocalSize)
            rResult.resize(TConditionLocalSize);

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        k_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());

        std::array<std::size_t, TFluidLocalSize> fluid_ids;
        fluid_condition.EquationIdArray(fluid_ids, rCurrentProcessInfo);
        AssignSubArray(fluid_ids, rResult, VelPresBlock());
        std::array<std::size_t, TKLocalSize> k_ids;
        k_condition.EquationIdArray(k_ids, rCurrentProcessInfo);
        AssignSubArray(k_ids, rResult, KBlock());
        std::array<std::size_t, TEpsilonLocalSize> epsilon_ids;
        epsilon_condition.EquationIdArray(epsilon_ids, rCurrentProcessInfo);
        AssignSubArray(epsilon_ids, rResult, EpsilonBlock());
    }

    /**
     * determines the conditional list of DOFs
     * @param ConditionalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rConditionalDofList.size() != TConditionLocalSize)
            rConditionalDofList.resize(TConditionLocalSize);

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        k_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());

        std::array<Dof<double>::Pointer, TFluidLocalSize> fluid_dofs;
        fluid_condition.GetDofArray(fluid_dofs, rCurrentProcessInfo);
        AssignSubArray(fluid_dofs, rConditionalDofList, VelPresBlock());
        std::array<Dof<double>::Pointer, TKLocalSize> k_dofs;
        k_condition.GetDofArray(k_dofs, rCurrentProcessInfo);
        AssignSubArray(k_dofs, rConditionalDofList, KBlock());
        std::array<Dof<double>::Pointer, TEpsilonLocalSize> epsilon_dofs;
        epsilon_condition.GetDofArray(epsilon_dofs, rCurrentProcessInfo);
        AssignSubArray(epsilon_dofs, rConditionalDofList, EpsilonBlock());
    }

    /// Returns the adjoint values stored in this condition's nodes.
    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TConditionLocalSize)
            rValues.resize(TConditionLocalSize, false);

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        k_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());

        std::array<double, TFluidLocalSize> fluid_values;
        fluid_condition.GetValuesArray(fluid_values, Step);
        AssignSubArray(fluid_values, rValues, VelPresBlock());
        std::array<double, TKLocalSize> k_values;
        k_condition.GetValuesArray(k_values, Step);
        AssignSubArray(k_values, rValues, KBlock());
        std::array<double, TEpsilonLocalSize> epsilon_values;
        epsilon_condition.GetValuesArray(epsilon_values, Step);
        AssignSubArray(epsilon_values, rValues, EpsilonBlock());
    }

    /// Returns the adjoint velocity values stored in this condition's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TConditionLocalSize)
            rValues.resize(TConditionLocalSize, false);

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        k_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());

        std::array<double, TFluidLocalSize> fluid_values;
        fluid_condition.GetFirstDerivativesArray(fluid_values, Step);
        AssignSubArray(fluid_values, rValues, VelPresBlock());
        std::array<double, TKLocalSize> k_values;
        k_condition.GetFirstDerivativesArray(k_values, Step);
        AssignSubArray(k_values, rValues, KBlock());
        std::array<double, TEpsilonLocalSize> epsilon_values;
        epsilon_condition.GetFirstDerivativesArray(epsilon_values, Step);
        AssignSubArray(epsilon_values, rValues, EpsilonBlock());
    }

    void GetSecondDerivativesVector(VectorType& rValues, int Step) override
    {
        if (rValues.size() != TConditionLocalSize)
            rValues.resize(TConditionLocalSize, false);

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointKCondition k_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        k_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());

        std::array<double, TFluidLocalSize> fluid_values;
        fluid_condition.GetSecondDerivativesArray(fluid_values, Step);
        AssignSubArray(fluid_values, rValues, VelPresBlock());
        std::array<double, TKLocalSize> k_values;
        k_condition.GetSecondDerivativesArray(k_values, Step);
        AssignSubArray(k_values, rValues, KBlock());
        std::array<double, TEpsilonLocalSize> epsilon_values;
        epsilon_condition.GetSecondDerivativesArray(epsilon_values, Step);
        AssignSubArray(epsilon_values, rValues, EpsilonBlock());
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition::"
                        "CalculateLocalSystem method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
            rLeftHandSideMatrix.size2() != TConditionLocalSize)
            rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);

        rLeftHandSideMatrix.clear();
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition::"
                        "CalculateRightHandSide method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
        AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

        fluid_condition.SetData(this->Data());
        epsilon_condition.SetData(this->Data());
        fluid_condition.SetFlags(this->GetFlags());
        epsilon_condition.SetFlags(this->GetFlags());

        if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
            rLeftHandSideMatrix.size2() != TConditionLocalSize)
            rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);

        rLeftHandSideMatrix.clear();

        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> vms_vms;
        fluid_condition.CalculateFirstDerivativesLHS(vms_vms, rCurrentProcessInfo);
        AssignSubMatrix(vms_vms, rLeftHandSideMatrix, VelPresBlock(), VelPresBlock());

        BoundedMatrix<double, TNumNodes, TFluidLocalSize> vms_k;
        fluid_condition.CalculateConditionResidualTurbulentKineticEnergyDerivatives(
            vms_k, rCurrentProcessInfo);
        AssignSubMatrix(vms_k, rLeftHandSideMatrix, KBlock(), VelPresBlock());

        BoundedMatrix<double, TNumNodes, TNumNodes> epsilon_k;
        epsilon_condition.CalculateConditionResidualTurbulentKineticEnergyDerivatives(
            epsilon_k, rCurrentProcessInfo);
        AssignSubMatrix(epsilon_k, rLeftHandSideMatrix, KBlock(), EpsilonBlock());

        BoundedMatrix<double, TNumNodes, TNumNodes> epsilon_epsilon;
        epsilon_condition.CalculateFirstDerivativesLHS(epsilon_epsilon, rCurrentProcessInfo);
        AssignSubMatrix(epsilon_epsilon, rLeftHandSideMatrix, EpsilonBlock(),
                        EpsilonBlock());

        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TConditionLocalSize ||
            rLeftHandSideMatrix.size2() != TConditionLocalSize)
            rLeftHandSideMatrix.resize(TConditionLocalSize, TConditionLocalSize, false);

        rLeftHandSideMatrix.clear();

        KRATOS_CATCH("");
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition::"
                        "CalculateMassMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition::"
                        "CalculateDampingMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            if (rOutput.size1() != TCoordLocalSize || rOutput.size2() != TConditionLocalSize)
                rOutput.resize(TCoordLocalSize, TConditionLocalSize, false);

            rOutput.clear();

            AdjointFluidCondition fluid_condition(this->Id(), this->pGetGeometry());
            AdjointEpsilonCondition epsilon_condition(this->Id(), this->pGetGeometry());

            fluid_condition.SetData(this->Data());
            epsilon_condition.SetData(this->Data());
            fluid_condition.SetFlags(this->GetFlags());
            epsilon_condition.SetFlags(this->GetFlags());

            BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize> vms_residuals;
            fluid_condition.CalculateSensitivityMatrix(
                rSensitivityVariable, vms_residuals, rCurrentProcessInfo);
            AssignSubMatrix(vms_residuals, rOutput, CoordBlock(), VelPresBlock());

            BoundedMatrix<double, TCoordLocalSize, TNumNodes> epsilon_residuals;
            epsilon_condition.CalculateSensitivityMatrix(
                rSensitivityVariable, epsilon_residuals, rCurrentProcessInfo);
            AssignSubMatrix(epsilon_residuals, rOutput, CoordBlock(), EpsilonBlock());
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
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
        std::stringstream buffer;
        buffer << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition #"
               << Condition::Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansEvmMonolithicKEpsilonVMSAdjointWallCondition #"
                 << Condition::Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Condition::pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    struct SubBlockLayout
    {
        std::size_t SubBlockOffset;
        std::size_t SubBlockSize;
        std::size_t BlockSize;
        std::size_t NumBlocks;
    };

    constexpr SubBlockLayout CoordBlock()
    {
        return {0, TDim, TDim, TNumNodes};
    }

    constexpr SubBlockLayout VelBlock()
    {
        return {0, TDim, TDim + 3, TNumNodes};
    }

    constexpr SubBlockLayout VelPresBlock()
    {
        return {0, TDim + 1, TDim + 3, TNumNodes};
    }

    constexpr SubBlockLayout KBlock()
    {
        return {TDim + 1, 1, TDim + 3, TNumNodes};
    }

    constexpr SubBlockLayout EpsilonBlock()
    {
        return {TDim + 2, 1, TDim + 3, TNumNodes};
    }

    constexpr std::size_t MonolithicIndex(SubBlockLayout L, std::size_t SubIndex)
    {
        return SubIndex + (SubIndex / L.SubBlockSize) * (L.BlockSize - L.SubBlockSize) +
               L.SubBlockOffset;
    }

    template <class TMatrix1, class TMatrix2>
    void AssignSubMatrix(const TMatrix1& rSrc, TMatrix2& rDest, SubBlockLayout RowLayout, SubBlockLayout ColLayout)
    {
        for (std::size_t i = 0; i < rSrc.size1(); ++i)
            for (std::size_t j = 0; j < rSrc.size2(); ++j)
                rDest(MonolithicIndex(RowLayout, i), MonolithicIndex(ColLayout, j)) =
                    rSrc(i, j);
    }

    template <class TArray1, class TArray2>
    void AssignSubArray(const TArray1& rSrc, TArray2& rDest, SubBlockLayout Layout)
    {
        for (std::size_t i = 0; i < rSrc.size(); ++i)
            rDest[MonolithicIndex(Layout, i)] = rSrc[i];
    }

    ///@}
};

///@}

} // namespace Kratos

#endif