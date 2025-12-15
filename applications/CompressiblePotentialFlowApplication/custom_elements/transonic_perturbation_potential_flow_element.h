//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:     Eloisa Baez Jones, Inigo Lopez, Marc Nuñez and Riccardo Rossi
//

#if !defined(KRATOS_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "custom_utilities/potential_flow_utilities.h"

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

template <int Dim, int NumNodes>
class TransonicPerturbationPotentialFlowElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    typedef PointerVector<GeometryType> GeometriesArrayType;
    typedef PotentialFlowUtilities::ElementalData<NumNodes, Dim> ElementalData;
    static constexpr int TNumNodes = NumNodes;
    static constexpr int TDim = Dim;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of TransonicPerturbationPotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransonicPerturbationPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constructor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit TransonicPerturbationPotentialFlowElement(IndexType NewId = 0)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    TransonicPerturbationPotentialFlowElement(IndexType NewId,
                                        const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    TransonicPerturbationPotentialFlowElement(IndexType NewId,
                                        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    TransonicPerturbationPotentialFlowElement(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    TransonicPerturbationPotentialFlowElement(TransonicPerturbationPotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    TransonicPerturbationPotentialFlowElement(TransonicPerturbationPotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~TransonicPerturbationPotentialFlowElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TransonicPerturbationPotentialFlowElement& operator=(TransonicPerturbationPotentialFlowElement const& rOther) = delete;

    /// Move operator.
    TransonicPerturbationPotentialFlowElement& operator=(TransonicPerturbationPotentialFlowElement&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId,
                           NodesArrayType const& ThisNodes) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& CurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rCurrentProcessInfo) const override;

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
protected:

    BoundedVector<double, NumNodes + 1> AssembleDensityDerivativeAndShapeFunctions(const double densityDerivativeWRTVelocitySquared, const double densityDerivativeWRTUpwindVelocitySquared, const array_1d<double, Dim>& velocity, const array_1d<double, Dim>& upwindVelocity,const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSideNormalElement(VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSideNormalElement(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSideWakeElement(MatrixType& rLeftHandSideMatrix,
                                          const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateRightHandSideWakeElement(VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo);


    virtual void AssembleSupersonicLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const double densityDerivativeWRTVelocity,
        const double densityDerivativeWRTUpwindVelocity,
        const array_1d<double, Dim> velocity,
        const array_1d<double, Dim> upwindVelocity,
        const ProcessInfo& rCurrentProcessInfo);

    virtual void CalculateLeftHandSideContribution(BoundedMatrix<double, NumNodes, NumNodes>& rLhs_total,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const array_1d<double, Dim>& rVelocity,
                                         const ElementalData& rData);

    virtual void CalculateRightHandSideContribution(BoundedVector<double, NumNodes>& rRhs_total,
                                        const double rDensity,
                                        const array_1d<double, Dim>& rVelocity);

    inline GlobalPointer<Element> pGetUpwindElement() const;

    bool CheckUpwindElement();

    void pSetUpwindElement(GlobalPointer<Element> pUpwindElement);

    array_1d<double, 3> GetEdgeNormal(const GeometryType& rEdge);

    virtual void FindUpwindElement(const ProcessInfo& rCurrentProcessInfo);

private:
    ///@}
    ///@name Member Variables
    ///@{

    GlobalPointer<Element> mpUpwindElement;

    ///@name Private Operators
    ///@{

    void GetWakeDistances(array_1d<double,
                         NumNodes>& distances) const;

    void GetEquationIdVectorExtendedElement(EquationIdVectorType& rResult) const;

    void AddUpwindEquationId(EquationIdVectorType& rResult) const;

    void GetEquationIdVectorNormalElement(EquationIdVectorType& rResult) const;

    void GetEquationIdVectorKuttaElement(EquationIdVectorType& rResult) const;

    void GetEquationIdVectorWakeElement(EquationIdVectorType& rResult) const;

    void GetDofListNormalElement(DofsVectorType& rElementalDofList) const;

    void GetDofListKuttaElement(DofsVectorType& rElementalDofList) const;

    void GetDofListWakeElement(DofsVectorType& rElementalDofList) const;

    void CalculateLeftHandSideSubsonicElement(MatrixType& rLeftHandSideMatrix,
                                            const ProcessInfo& rCurrentProcessInfo);

    // void CalculateRightHandSideSupersonicElement(VectorType& rRightHandSideVector,
    //                                         const ProcessInfo& rCurrentProcessInfo);

    BoundedMatrix<double, NumNodes, NumNodes> CalculateLeftHandSideWakeConditions(
                                            const ElementalData& rData,
                                            const ProcessInfo& rCurrentProcessInfo);

    BoundedVector<double, NumNodes> CalculateRightHandSideWakeConditions(
                                            const ElementalData& rData,
                                            const ProcessInfo& rCurrentProcessInfo,
                                            const array_1d<double, Dim>& rDiff_velocity);

    void CalculateLeftHandSideSubdividedElement(Matrix& lhs_positive,
                                               Matrix& lhs_negative,
                                               const ProcessInfo& rCurrentProcessInfo);
    void CalculateVolumesSubdividedElement(double& rUpper_vol,
                                           double& rLower_vol,
                                           const ProcessInfo& rCurrentProcessInfo);

    void ComputeLHSGaussPointContribution(const double weight,
                                          Matrix& lhs,
                                          const ElementalData& data) const;

    void AssignLeftHandSideSubdividedElement(
        Matrix& rLeftHandSideMatrix,
        Matrix& lhs_positive,
        Matrix& lhs_negative,
        const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
        const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
        const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
        const ElementalData& data) const;

    void AssignLeftHandSideWakeElement(MatrixType& rLeftHandSideMatrix,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                    const ElementalData& rData) const;

    void AssignLeftHandSideWakeNode(MatrixType& rLeftHandSideMatrix,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                    const ElementalData& rData,
                                    unsigned int row) const;

    void AssignRightHandSideWakeNode(VectorType& rRightHandSideVector,
                                    const BoundedVector<double, NumNodes>& rUpper_rhs,
                                    const BoundedVector<double, NumNodes>& rLower_rhs,
                                    const BoundedVector<double, NumNodes>& rWake_rhs,
                                    const ElementalData& rData,
                                    unsigned int& rRow) const;

    array_1d<size_t, NumNodes> GetAssemblyKey(const GeometryType& rGeom, const GeometryType& rUpwindGeom, const ProcessInfo& rCurrentProcessInfo);

    void FindUpwindEdge(GeometryType& rUpwindEdge,
                        const ProcessInfo& rCurrentProcessInfo);

    void GetElementGeometryBoundary(GeometriesArrayType& rElementGeometryBoundary);


    void SelectUpwindElement(std::vector<IndexType>& rUpwindElementNodesIds,
                             GlobalPointersVector<Element>& rUpwindElementCandidates);

    int GetAdditionalUpwindNodeIndex() const;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class TransonicPerturbationPotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H  defined
