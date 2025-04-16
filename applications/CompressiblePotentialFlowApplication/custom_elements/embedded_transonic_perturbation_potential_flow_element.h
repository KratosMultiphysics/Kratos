//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:     Marc Nuñez, Eloisa Baez Jones, Inigo Lopez and Riccardo Rossi
//

#if !defined(KRATOS_EMBEDDED_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H)
#define KRATOS_EMBEDDED_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "transonic_perturbation_potential_flow_element.h"
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

template <int TDim, int TNumNodes>
class EmbeddedTransonicPerturbationPotentialFlowElement : public TransonicPerturbationPotentialFlowElement<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{

    typedef TransonicPerturbationPotentialFlowElement<TDim,TNumNodes> BaseType;

    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;
    typedef PotentialFlowUtilities::ElementalData<TNumNodes, TDim> ElementalData;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of EmbeddedTransonicPerturbationPotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EmbeddedTransonicPerturbationPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constructor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit EmbeddedTransonicPerturbationPotentialFlowElement(IndexType NewId = 0)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    EmbeddedTransonicPerturbationPotentialFlowElement(IndexType NewId,
                                        const NodesArrayType& ThisNodes)
        : TransonicPerturbationPotentialFlowElement<TDim,TNumNodes>(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    EmbeddedTransonicPerturbationPotentialFlowElement(IndexType NewId,
                                        typename GeometryType::Pointer pGeometry)
        : TransonicPerturbationPotentialFlowElement<TDim,TNumNodes>(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    EmbeddedTransonicPerturbationPotentialFlowElement(IndexType NewId,
                                     typename GeometryType::Pointer pGeometry,
                                     typename PropertiesType::Pointer pProperties)
        : TransonicPerturbationPotentialFlowElement<TDim,TNumNodes>(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    EmbeddedTransonicPerturbationPotentialFlowElement(EmbeddedTransonicPerturbationPotentialFlowElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    EmbeddedTransonicPerturbationPotentialFlowElement(EmbeddedTransonicPerturbationPotentialFlowElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~EmbeddedTransonicPerturbationPotentialFlowElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    EmbeddedTransonicPerturbationPotentialFlowElement& operator=(EmbeddedTransonicPerturbationPotentialFlowElement const& rOther) = delete;

    /// Move operator.
    EmbeddedTransonicPerturbationPotentialFlowElement& operator=(EmbeddedTransonicPerturbationPotentialFlowElement&& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            typename GeometryType::Pointer pGeom,
                            typename PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId,
                           NodesArrayType const& ThisNodes) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rCurrentProcessInfo) override;

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


    void CalculateLeftHandSideKuttaWakeElement(MatrixType& rLeftHandSideMatrix,
                                          const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSideKuttaWakeElement(VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSideContribution(BoundedMatrix<double, TNumNodes, TNumNodes>& rLhs_total,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const array_1d<double, TDim>& rVelocity,
                                         const ElementalData& rData) override;

    void CalculateRightHandSideContribution(BoundedVector<double, TNumNodes>& rRhs_total,
                                        const double rDensity,
                                        const array_1d<double, TDim>& rVelocity) override;


    void AssembleSupersonicLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                        const double densityDerivativeWRTVelocity,
                                        const double densityDerivativeWRTUpwindVelocity,
                                        const array_1d<double, TDim> velocity,
                                        const array_1d<double, TDim> upwindVelocity,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    ModifiedShapeFunctions::Pointer pGetModifiedShapeFunctions(Vector& rDistances);

    void FindUpwindElement(const ProcessInfo& rCurrentProcessInfo) override;

private:
    ///@}
    ///@name Member Variables
    ///@{

    void SelectUpwindElement(std::vector<IndexType>& rUpwindElementNodesIds,
                            GlobalPointersVector<Element>& rUpwindElementCandidates,
                            const ProcessInfo& rCurrentProcessInfo);

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

}; // Class EmbeddedTransonicPerturbationPotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_TRANSONIC_PERTURBATION_POTENTIAL_FLOW_ELEMENT_H  defined
