//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_ADJOINT_FLUID_CONDITION_H_INCLUDED)
#define KRATOS_ADJOINT_FLUID_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "conditions/mesh_condition.h"
#include "includes/process_info.h"
#include "includes/properties.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes>
class AdjointFluidCondition : public MeshCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AdjointFluidCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFluidCondition);

    using BaseType = MeshCondition;

    using NodeType = Node<3>;

    using PropertiesType = Properties;

    using GeometryType = Geometry<NodeType>;

    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    static constexpr IndexType BlockSize = TDim + 1;

    static constexpr IndexType LocalSize = BlockSize * TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit AdjointFluidCondition(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    AdjointFluidCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    AdjointFluidCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    AdjointFluidCondition(IndexType NewId,
                          GeometryType::Pointer pGeometry,
                          PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    AdjointFluidCondition(AdjointFluidCondition const& rOther)
        : BaseType(rOther)
    {
    }

    /// Destructor.
    ~AdjointFluidCondition() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    AdjointFluidCondition& operator=(AdjointFluidCondition const& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFluidCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFluidCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition =
            Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void GetValuesVector(Vector& values, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& values, int Step = 0) const override;

    void CalculateFirstDerivativesLHS(Matrix& rLeftHandSideMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(
        const Variable<array_1d<double, 3>>& rDesignVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "AdjointFluidCondition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AdjointFluidCondition #" << this->Id();
    }

    ///@}

private:
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MeshCondition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MeshCondition);
    }

    ///@}

}; // Class AdjointFluidCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                AdjointFluidCondition<TDim, TNumNodes>& rThis);

/// output stream function

template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const AdjointFluidCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.
#endif // KRATOS_ADJOINT_FLUID_CONDITION_H_INCLUDED  defined
