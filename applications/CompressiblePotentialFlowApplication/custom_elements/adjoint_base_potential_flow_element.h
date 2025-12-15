//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nunez, based on A. Geiser, M. Fusseder, I. Lopez and R. Rossi work
//

#if !defined(KRATOS_ADJOINT_BASE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_ADJOINT_BASE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED


// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{

template <class TPrimalElement>
class AdjointBasePotentialFlowElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;
    typedef PointerVector<GeometryType> GeometriesArrayType;
    static constexpr int NumNodes = TPrimalElement::TNumNodes;
    static constexpr int Dim = TPrimalElement::TDim;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of AdjointBasePotentialFlowElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointBasePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit AdjointBasePotentialFlowElement(IndexType NewId = 0)
     : Element(NewId),
     mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId))
    {};

    AdjointBasePotentialFlowElement(IndexType NewId,
                        GeometryType::Pointer pGeometry)
     : Element(NewId, pGeometry),
      mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId, pGeometry))
    {
    }

    AdjointBasePotentialFlowElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
     : Element(NewId, pGeometry, pProperties),
      mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId, pGeometry, pProperties))
    {
    }
    /**
     * Copy Constructor
     */
    AdjointBasePotentialFlowElement(AdjointBasePotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~AdjointBasePotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AdjointBasePotentialFlowElement & operator=(AdjointBasePotentialFlowElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    IntegrationMethod GetIntegrationMethod() const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<int>& rVariable,
            std::vector<int>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step=0) const override;

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    Element::Pointer pGetPrimalElement();

    inline GlobalPointer<Element> pGetUpwindElement() const;

    bool CheckUpwindElement();

    void pSetUpwindElement(GlobalPointer<Element> pUpwindElement);

    array_1d<double, 3> GetEdgeNormal(const GeometryType& rEdge);

    virtual void FindUpwindElement(const ProcessInfo& rCurrentProcessInfo);

protected:

    Element::Pointer mpPrimalElement;

    GlobalPointer<Element> mpUpwindElement;

    void GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances) const;

private:

    void AddUpwindEquationId(EquationIdVectorType& rResult) const;

    void FindUpwindEdge(GeometryType& rUpwindEdge,
                        const ProcessInfo& rCurrentProcessInfo);

    void GetElementGeometryBoundary(GeometriesArrayType& rElementGeometryBoundary);


    void SelectUpwindElement(std::vector<IndexType>& rUpwindElementNodesIds,
                             GlobalPointersVector<Element>& rUpwindElementCandidates);

    int GetAdditionalUpwindNodeIndex() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class AdjointBasePotentialFlowElement


} // namespace Kratos.

#endif
