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

#if !defined(KRATOS_RANS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_RANS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes
#include <algorithm>
#include <iterator>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"

// Application includes
#include "custom_elements/evm_k_epsilon/rans_evm_epsilon_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_adjoint.h"
#include "custom_elements/evm_k_epsilon/rans_evm_k_epsilon_vms_adjoint.h"
#include "rans_application_variables.h"

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

template <unsigned int TDim, unsigned int TNumNodes = TDim + 1>
class KRATOS_API(RANS_APPLICATION) RansEvmMonolithicKEpsilonVMSAdjoint : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement) : mpElement{pElement}
        {
        }

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 3);
            std::size_t index = 0;
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
            if (TDim == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Z, Step);
            }
            rVector[index++] = IndirectScalar<double>{}; // pressure
            rVector[index++] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_2, Step);
            rVector[index] = MakeIndirectScalar(r_node, RANS_SCALAR_2_ADJOINT_2, Step);
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 3);
            std::size_t index = 0;
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
            rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
            if (TDim == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Z, Step);
            }
            rVector[index++] = IndirectScalar<double>{}; // pressure
            rVector[index++] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_3, Step);
            rVector[index] = MakeIndirectScalar(r_node, RANS_SCALAR_2_ADJOINT_3, Step);
        }

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 3);
            std::size_t index = 0;
            rVector[index++] =
                MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
            rVector[index++] =
                MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
            if (TDim == 3)
            {
                rVector[index++] =
                    MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Z, Step);
            }
            rVector[index++] = IndirectScalar<double>{}; // pressure
            rVector[index++] = MakeIndirectScalar(r_node, RANS_AUX_ADJOINT_SCALAR_1, Step);
            rVector[index] = MakeIndirectScalar(r_node, RANS_AUX_ADJOINT_SCALAR_2, Step);
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(3);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
            rVariables[1] = &RANS_SCALAR_1_ADJOINT_2;
            rVariables[2] = &RANS_SCALAR_2_ADJOINT_2;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(3);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
            rVariables[1] = &RANS_SCALAR_1_ADJOINT_3;
            rVariables[2] = &RANS_SCALAR_2_ADJOINT_3;
        }

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(3);
            rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
            rVariables[1] = &RANS_AUX_ADJOINT_SCALAR_1;
            rVariables[2] = &RANS_AUX_ADJOINT_SCALAR_2;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    // defining the base type
    typedef Element BaseType;
    // defining the base adjoint base fluid element type
    typedef RansEvmKEpsilonVMSAdjoint<TDim, TNumNodes> AdjointFluidElement;
    // defining the k element type
    typedef RansEvmKAdjoint<TDim, TNumNodes> AdjointKElement;
    // defining the epsilon element type
    typedef RansEvmEpsilonAdjoint<TDim, TNumNodes> AdjointEpsilonElement;

    constexpr static unsigned int TFluidBlockSize = (TDim + 1);

    constexpr static unsigned int TFluidLocalSize = TFluidBlockSize * TNumNodes;

    constexpr static unsigned int TKBlockSize = 1;

    constexpr static unsigned int TKLocalSize = TKBlockSize * TNumNodes;

    constexpr static unsigned int TEpsilonBlockSize = 1;

    constexpr static unsigned int TEpsilonLocalSize = TEpsilonBlockSize * TNumNodes;

    constexpr static unsigned int TElementBlockSize =
        (TFluidBlockSize + TKBlockSize + TEpsilonBlockSize);

    constexpr static unsigned int TElementLocalSize = TElementBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    // variable definitions
    typedef std::size_t IndexType;

    typedef Element::NodeType NodeType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::VectorType VectorType;

    typedef Element::MatrixType MatrixType;

    ///@name Pointer Definitions
    /// Pointer definition of RansEvmMonolithicKEpsilonVMSAdjoint
    KRATOS_CLASS_POINTER_DEFINITION(RansEvmMonolithicKEpsilonVMSAdjoint);

    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit RansEvmMonolithicKEpsilonVMSAdjoint(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    RansEvmMonolithicKEpsilonVMSAdjoint(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    RansEvmMonolithicKEpsilonVMSAdjoint(IndexType NewId,
                                        GeometryType::Pointer pGeometry,
                                        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~RansEvmMonolithicKEpsilonVMSAdjoint() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override;

    void GetSecondDerivativesVector(VectorType& rValues, int Step) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& /*rCurrentProcessInfo*/) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "RansEvmMonolithicKEpsilonVMSAdjoint #" << Element::Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansEvmMonolithicKEpsilonVMSAdjoint #" << Element::Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Element::pGetGeometry()->PrintData(rOStream);
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