//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on A. Geiser, M. Fusseder, I. Lopez and R. Rossi work
//

#ifndef KRATOS_ADJOINT_POTENTIAL_WALL_CONDITION_H
#define KRATOS_ADJOINT_POTENTIAL_WALL_CONDITION_H

namespace Kratos
{


template <class TPrimalCondition>
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) AdjointPotentialWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    static constexpr int TNumNodes = TPrimalCondition::NumNodes;
    static constexpr int TDim = TPrimalCondition::Dim;

    /// Pointer definition of AdjointPotentialWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointPotentialWallCondition);

    typedef Element::WeakPointer ElementWeakPointerType;

    typedef Element::Pointer ElementPointerType;

    AdjointPotentialWallCondition(IndexType NewId = 0)
    : Condition(NewId),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGetGeometry()))
    {
    }

    AdjointPotentialWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGeometry))
    {
    }

    AdjointPotentialWallCondition(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties),
      mpPrimalCondition(Kratos::make_intrusive<TPrimalCondition>(NewId, pGeometry, pProperties))
    {
    }

    /// Destructor.
    ~AdjointPotentialWallCondition() override {}

    AdjointPotentialWallCondition & operator=(AdjointPotentialWallCondition const& rOther)
    {
        Condition::operator=(rOther);
        return *this;
    }


    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId,
                              Condition::GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

    // Find the condition's parent element.
    void Initialize() override;

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;


    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step=0)  override;

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo)  override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;


protected:

    Condition::Pointer mpPrimalCondition;


private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;



}; // Class AdjointPotentialWallCondition

///@}

/// input stream function
template <class TPrimalCondition>
inline std::istream& operator>>(std::istream& rIStream,
                                AdjointPotentialWallCondition<TPrimalCondition>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TPrimalCondition>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const AdjointPotentialWallCondition<TPrimalCondition>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_POTENTIAL_WALL_CONDITION_H
