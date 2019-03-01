//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#ifndef KRATOS_POTENTIAL_WALL_CONDITION_H
#define KRATOS_POTENTIAL_WALL_CONDITION_H


#include "includes/kratos_flags.h"

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/condition.h"
#include "includes/process_info.h"

// Application includes
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{


template <unsigned int TDim, unsigned int TNumNodes = TDim>
class IncompressibleAdjointPotentialWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IncompressibleAdjointPotentialWallCondition
    KRATOS_CLASS_POINTER_DEFINITION(IncompressibleAdjointPotentialWallCondition);

    typedef Node < 3 > NodeType;

    typedef Properties PropertiesType;

    typedef Geometry<NodeType> GeometryType;

    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector< Dof<double>::Pointer > DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef Element::WeakPointer ElementWeakPointerType;
    
    typedef Element::Pointer ElementPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId = 0):
        Condition(NewId)
    {
    }

    IncompressibleAdjointPotentialWallCondition(Condition::Pointer pPrimalCondition)
                    : Condition(pPrimalCondition->Id(), pPrimalCondition->pGetGeometry(), pPrimalCondition->pGetProperties())
                    , mpPrimalCondition(pPrimalCondition)
    {
    };


    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
                           const NodesArrayType& ThisNodes):
        Condition(NewId,ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry):
        Condition(NewId,pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    IncompressibleAdjointPotentialWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties):
        Condition(NewId,pGeometry,pProperties)
    {
    }

    /// Copy constructor.
    IncompressibleAdjointPotentialWallCondition(IncompressibleAdjointPotentialWallCondition const& rOther):
        Condition(rOther)
    {
    }

    /// Destructor.
    ~IncompressibleAdjointPotentialWallCondition() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor
    IncompressibleAdjointPotentialWallCondition & operator=(IncompressibleAdjointPotentialWallCondition const& rOther)
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

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo)  override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    Condition::Pointer mpPrimalCondition;


private:


    ///@}
    ///@name Member Variables
    ///@{

    bool mInitializeWasPerformed = false;
    ElementWeakPointerType mpElement;

    void CalculateNormal2D(array_1d<double, 3>& An) const;

    void CalculateNormal3D(array_1d<double, 3>& An) const;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Operators
    ///@{

    inline ElementPointerType pGetElement() const;

    void GetElementCandidates(WeakPointerVector<Element>& ElementCandidates,
                              const GeometryType& rGeom) const;

    void GetSortedIds(std::vector<IndexType>& Ids, const GeometryType& rGeom) const;

    void FindParentElement(std::vector<IndexType>& NodeIds,
                           std::vector<IndexType>& ElementNodeIds,
                           WeakPointerVector<Element> ElementCandidates);


}; // Class IncompressibleAdjointPotentialWallCondition

///@}

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>& rThis)
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
