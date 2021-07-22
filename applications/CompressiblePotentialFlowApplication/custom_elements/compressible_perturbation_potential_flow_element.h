//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez, Marc Nuñez and Riccardo Rossi
//

#if !defined(KRATOS_COMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_COMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"
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
class CompressiblePerturbationPotentialFlowElement : public Element
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
    /// Pointer definition of CompressiblePerturbationPotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(CompressiblePerturbationPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit CompressiblePerturbationPotentialFlowElement(IndexType NewId = 0)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    CompressiblePerturbationPotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    CompressiblePerturbationPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    CompressiblePerturbationPotentialFlowElement(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    CompressiblePerturbationPotentialFlowElement(CompressiblePerturbationPotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    CompressiblePerturbationPotentialFlowElement(CompressiblePerturbationPotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~CompressiblePerturbationPotentialFlowElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CompressiblePerturbationPotentialFlowElement& operator=(CompressiblePerturbationPotentialFlowElement const& rOther) = delete;

    /// Move operator.
    CompressiblePerturbationPotentialFlowElement& operator=(CompressiblePerturbationPotentialFlowElement&& rOther) = delete;

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

    BoundedMatrix<double, NumNodes, NumNodes> CalculateLeftHandSideWakeConditions(
                                            const ElementalData<NumNodes, Dim>& rData,
                                            const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSideWakeElement(VectorType& rRightHandSideVector,
                                         const ProcessInfo& rCurrentProcessInfo);

    BoundedVector<double, NumNodes> CalculateRightHandSideWakeConditions(
                                            const ElementalData<NumNodes, Dim>& rData,
                                            const ProcessInfo& rCurrentProcessInfo,
                                            const array_1d<double, Dim>& rDiff_velocity);

    void CalculateLeftHandSideContribution(BoundedMatrix<double, NumNodes, NumNodes>& rLhs_total,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const array_1d<double, Dim>& rVelocity,
                                         const ElementalData<NumNodes, Dim>& rData);

    void CalculateRightHandSideContribution(BoundedVector<double, NumNodes>& rRhs_total,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const array_1d<double, Dim>& rVelocity,
                                         const ElementalData<NumNodes, Dim>& rData);

    void CalculateLeftHandSideSubdividedElement(Matrix& lhs_positive,
                                               Matrix& lhs_negative,
                                               const ProcessInfo& rCurrentProcessInfo);
    void CalculateVolumesSubdividedElement(double& rUpper_vol,
                                           double& rLower_vol,
                                           const ProcessInfo& rCurrentProcessInfo);

    void ComputeLHSGaussPointContribution(const double weight,
                                          Matrix& lhs,
                                          const ElementalData<NumNodes, Dim>& data) const;

    void AssignLeftHandSideSubdividedElement(Matrix& rLeftHandSideMatrix,
                                             Matrix& lhs_positive,
                                             Matrix& lhs_negative,
                                             const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
                                             const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
                                             const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                             const ElementalData<NumNodes, Dim>& data) const;

    void AssignLeftHandSideWakeElement(MatrixType& rLeftHandSideMatrix,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                    const ElementalData<NumNodes, Dim>& rData) const;

    void AssignLeftHandSideWakeNode(MatrixType& rLeftHandSideMatrix,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rUpper_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLower_lhs_total,
                                    const BoundedMatrix<double, NumNodes, NumNodes>& rLhs_wake_condition,
                                    const ElementalData<NumNodes, Dim>& rData,
                                    unsigned int row) const;

    void AssignRightHandSideWakeNode(VectorType& rRightHandSideVector,
                                   const BoundedVector<double, NumNodes>& rUpper_rhs,
                                   const BoundedVector<double, NumNodes>& rLower_rhs,
                                   const BoundedVector<double, NumNodes>& rWake_rhs,
                                   const ElementalData<NumNodes, Dim>& rData,
                                   unsigned int& rRow) const;

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

}; // Class CompressiblePerturbationPotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H  defined
