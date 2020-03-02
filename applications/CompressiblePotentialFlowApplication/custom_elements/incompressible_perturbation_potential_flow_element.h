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

#if !defined(KRATOS_INCOMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_INCOMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H

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
class IncompressiblePerturbationPotentialFlowElement : public Element
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
    /// Pointer definition of IncompressiblePerturbationPotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(IncompressiblePerturbationPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit IncompressiblePerturbationPotentialFlowElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePerturbationPotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    IncompressiblePerturbationPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    IncompressiblePerturbationPotentialFlowElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    IncompressiblePerturbationPotentialFlowElement(IncompressiblePerturbationPotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    IncompressiblePerturbationPotentialFlowElement(IncompressiblePerturbationPotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~IncompressiblePerturbationPotentialFlowElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePerturbationPotentialFlowElement& operator=(IncompressiblePerturbationPotentialFlowElement const& rOther) = delete;

    /// Move operator.
    IncompressiblePerturbationPotentialFlowElement& operator=(IncompressiblePerturbationPotentialFlowElement&& rOther) = delete;

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
                              ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                     std::vector<double>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<int>& rVariable,
                                     std::vector<int>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                     std::vector<array_1d<double, 3>>& rValues,
                                     const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Inquiry
    ///@{

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

    void CalculateLeftHandSideNormalElement(MatrixType& rLeftHandSideMatrix,
                                           const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSideNormalElement(VectorType& rRightHandSideVector,
                                           const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSideWakeElement(MatrixType& rLeftHandSideMatrix,
                                         const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSideWakeElement(VectorType& rRightHandSideVector,
                                         const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSideSubdividedElement(BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
                                               BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
                                               const ProcessInfo& rCurrentProcessInfo);

    void CalculateVolumesSubdividedElement(double& rUpper_vol,
                                           double& rLower_vol,
                                           const ProcessInfo& rCurrentProcessInfo);

    void ComputeLHSGaussPointContribution(const double weight,
                                          BoundedMatrix<double, NumNodes, NumNodes>& lhs,
                                          const ElementalData<NumNodes, Dim>& data) const;

    void AssignLeftHandSideSubdividedElement(MatrixType& rLeftHandSideMatrix,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_positive,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_negative,
                                            BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                            const ElementalData<NumNodes, Dim>& data) const;

    void AssignLeftHandSideWakeElement(MatrixType& rLeftHandSideMatrix,
                                      BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                      const ElementalData<NumNodes, Dim>& data) const;

    void AssignLeftHandSideWakeNode(MatrixType& rLeftHandSideMatrix,
                                   BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
                                   const ElementalData<NumNodes, Dim>& data,
                                   unsigned int& row) const;

    void AssignRightHandSideWakeNode(VectorType& rRightHandSideVector,
                                   const BoundedVector<double, NumNodes>& rUpper_rhs,
                                   const BoundedVector<double, NumNodes>& rLower_rhs,
                                   const BoundedVector<double, NumNodes>& rWake_rhs,
                                   const ElementalData<NumNodes, Dim>& rData,
                                   unsigned int& rRow) const;

    void ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo);

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class IncompressiblePerturbationPotentialFlowElement

///@}
} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H  defined
