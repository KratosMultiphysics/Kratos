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

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H

// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"

#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <int Dim, int NumNodes>
class IncompressiblePotentialFlowElement : public Element
{
public:
    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double, TNumNodes> potentials, distances;
        double vol;

        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        array_1d<double, TNumNodes> N;
    };
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;
    static constexpr int TNumNodes = NumNodes;
    static constexpr int TDim = Dim;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IncompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit IncompressiblePotentialFlowElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowElement(IncompressiblePotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    IncompressiblePotentialFlowElement(IncompressiblePotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowElement& operator=(IncompressiblePotentialFlowElement const& rOther) = delete;

    /// Move operator.
    IncompressiblePotentialFlowElement& operator=(IncompressiblePotentialFlowElement&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
                                     std::vector<int>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                     std::vector<array_1d<double, 3>>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

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

private:
    ///@name Private Operators
    ///@{

    void GetWakeDistances(array_1d<double, NumNodes>& distances) const;

    void GetEquationIdVectorNormalElement(EquationIdVectorType& rResult) const;

    void GetEquationIdVectorKuttaElement(EquationIdVectorType& rResult) const;

    void GetEquationIdVectorWakeElement(EquationIdVectorType& rResult) const;

    void GetDofListNormalElement(DofsVectorType& rElementalDofList) const;

    void GetDofListKuttaElement(DofsVectorType& rElementalDofList) const;

    void GetDofListWakeElement(DofsVectorType& rElementalDofList) const;

    void CalculateLocalSystemNormalElement(MatrixType& rLeftHandSideMatrix,
                                           VectorType& rRightHandSideVector,
                                           const ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystemWakeElement(MatrixType& rLeftHandSideMatrix,
                                         VectorType& rRightHandSideVector,
                                         const ProcessInfo& rCurrentProcessInfo);

    void CalculateLocalSystemSubdividedElement(BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
                                               BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
                                               const ProcessInfo& rCurrentProcessInfo);

    void ComputeLHSGaussPointContribution(const double weight,
                                          BoundedMatrix<double, NumNodes, NumNodes>& lhs,
                                          const ElementalData<NumNodes, Dim>& data) const;

    void CalculateBlockLeftHandSideWakeElement(BoundedMatrix<double, NumNodes, NumNodes>& rLhs_total,
                                            BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                            const ElementalData<NumNodes, Dim>& rData,
                                            const ProcessInfo& rCurrentProcessInfo);

    void AssignLocalSystemSubdividedElement(MatrixType& rLeftHandSideMatrix,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                            BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                            const ElementalData<NumNodes, Dim>& data) const;

    void AssignLocalSystemWakeElement(MatrixType& rLeftHandSideMatrix,
                                      BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                      BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                      const ElementalData<NumNodes, Dim>& data) const;

    void AssignLocalSystemWakeNode(MatrixType& rLeftHandSideMatrix,
                                   BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                    BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                   const ElementalData<NumNodes, Dim>& data,
                                   unsigned int& row) const;

    void ComputeElementInternalEnergy();

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class IncompressiblePotentialFlowElement

///@}
} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H  defined
