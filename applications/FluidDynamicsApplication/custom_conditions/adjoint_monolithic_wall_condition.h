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

#if !defined(KRATOS_ADJOINT_MONOLITHIC_WALL_CONDITION_H_INCLUDED)
#define KRATOS_ADJOINT_MONOLITHIC_WALL_CONDITION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/serializer.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/indirect_scalar_fwd.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes = TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) AdjointMonolithicWallCondition : public Condition
{
    class ThisExtensions : public AdjointExtensions
    {
        Condition* mpCondition;

    public:
        explicit ThisExtensions(Condition* pCondition)
            : mpCondition{pCondition}
        {
        }

        void GetFirstDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetSecondDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetAuxiliaryVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override;

        void GetFirstDerivativesVariables(
            std::vector<VariableData const*>& rVariables) const override;

        void GetSecondDerivativesVariables(
            std::vector<VariableData const*>& rVariables) const override;

        void GetAuxiliaryVariables(
            std::vector<VariableData const*>& rVariables) const override;
    };

public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointMonolithicWallCondition);

    using NodeType = Node;

    using PropertiesType = Properties;

    using GeometryType = Geometry<NodeType>;

    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using VectorType = Vector;

    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using SizeType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using DofsArrayType = PointerVectorSet<Dof<double>, IndexedObject>;

    static constexpr IndexType TBlockSize = (TDim + 1);

    static constexpr IndexType TFluidLocalSize = TBlockSize * TNumNodes;

    static constexpr IndexType TCoordsLocalSize = TDim * TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointMonolithicWallCondition(
        IndexType NewId = 0) : Condition(NewId)
    {
    }

    AdjointMonolithicWallCondition(
        IndexType NewId,
        const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    AdjointMonolithicWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    AdjointMonolithicWallCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    AdjointMonolithicWallCondition(AdjointMonolithicWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    ~AdjointMonolithicWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    AdjointMonolithicWallCondition& operator=(AdjointMonolithicWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointMonolithicWallCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointMonolithicWallCondition>(NewId, pGeom, pProperties);
    }

    Condition::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override
    {
        Condition::Pointer pNewCondition =
            Create(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

        pNewCondition->SetData(this->GetData());
        pNewCondition->SetFlags(this->GetFlags());

        return pNewCondition;
    }

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(
        DofsVectorType& ConditionDofList,
        const ProcessInfo& CurrentProcessInfo) const override;

    void GetValuesVector(
        Vector& Values,
        int Step = 0) const override;

    void GetFirstDerivativesVector(
        Vector& Values,
        int Step = 0) const override;

    void GetSecondDerivativesVector(
        Vector& Values,
        int Step = 0) const override;

    /**
     * @brief This method is required to calculate residual
     *
     * This method is used in adjoint schemes to calculate residual
     * which is required to properly calculate adjoint contributions
     * when slip conditions are applied.
     *
     * @param rLeftHandSideMatrix       This is the left hand side, which is not used in adjoint schemes
     * @param rRightHandSideVector      This is the residual, which is used in adjoint schemes
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method is required to calculate residual
     *
     * This method is used in adjoint schemes to calculate residual
     * which is required to properly calculate adjoint contributions
     * when slip conditions are applied.
     *
     * @param rLeftHandSideMatrix       This is the left hand side, which is not used in adjoint schemes
     * @param rRightHandSideVector      This is the residual, which is used in adjoint schemes
     * @param rCurrentProcessInfo       Current process info
     */
    void CalculateLocalVelocityContribution(
        MatrixType& rDampingMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method gives the gradient matrix
     *
     * @param rLeftHandSideMatrix
     * @param rCurrentProcessInfo
     */
    void CalculateLeftHandSide(
        Matrix& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method gives the residual derivatives w.r.t. first derivatives
     *
     * @param rLeftHandSideMatrix
     * @param rCurrentProcessInfo
     */
    void CalculateFirstDerivativesLHS(
        Matrix& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method gives the residual derivatives w.r.t. second derivatives
     *
     * @param rLeftHandSideMatrix
     * @param rCurrentProcessInfo
     */
    void CalculateSecondDerivativesLHS(
        Matrix& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This method gives the residual derivatives w.r.t. sensitivity variable
     *
     * @param rDesignVariable
     * @param rOutput
     * @param rCurrentProcessInfo
     */
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
        buffer << "AdjointMonolithicWallCondition" << TDim << "D";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AdjointMonolithicWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    virtual void ApplyNeumannCondition(
        MatrixType& rLocalMatrix,
        VectorType& rLocalVector,
        const ProcessInfo& rCurrentProcessInfo) const;

    virtual void ApplyWallLaw(
        MatrixType& rLocalMatrix,
        VectorType& rLocalVector,
        const ProcessInfo& rCurrentProcessInfo) const;

    virtual void ApplyNeumannConditionShapeDerivatives(
        MatrixType& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo) const;

    virtual void ApplyWallLawStateDerivatives(
        MatrixType& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo) const;

    virtual void ApplyWallLawShapeDerivatives(
        MatrixType& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo) const;

    void CalculateData(
        double& rArea,
        double& rDetJ,
        array_1d<double, TDim>& rUnitNormal) const;

    void CalculateDataShapeDerivatives(
        BoundedVector<double, TCoordsLocalSize>& rAreaDerivatives,
        BoundedVector<double, TCoordsLocalSize>& rDetJDerivatives,
        BoundedMatrix<double, TCoordsLocalSize, TDim>& rUnitNormalDerivatives,
        const double Area,
        const double DetJ,
        const array_1d<double, TDim>& rUnitNormal) const;

    ///@}

private:
    ///@name Serialization
    ///@{

    double CalculateDetJ(const double Area) const;

    void CalculateDetJShapeDerivatives(
        BoundedVector<double, TCoordsLocalSize>& rDetJDerivatives,
        const BoundedVector<double, TCoordsLocalSize>& rAreaDerivatives) const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}

}; // Class AdjointMonolithicWallCondition

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(
    std::istream& rIStream,
    AdjointMonolithicWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const AdjointMonolithicWallCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_ADJOINT_MONOLITHIC_WALL_CONDITION_H_INCLUDED
