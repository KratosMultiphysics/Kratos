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

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_WAKE_ELEMENT_H)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_WAKE_ELEMENT_H

// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"
namespace Kratos
{
///@name Kratos Classes
///@{

template <int Dim, int NumNodes>
class IncompressiblePotentialFlowWakeElement : public Element
{
public:
    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double, TNumNodes> phis, distances;
        double vol;

        BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        array_1d<double, TNumNodes> N;
    };

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowWakeElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowWakeElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    explicit IncompressiblePotentialFlowWakeElement(IndexType NewId = 0){}

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowWakeElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes){}

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowWakeElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry){}

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowWakeElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties){}

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowWakeElement(IncompressiblePotentialFlowWakeElement const& rOther) = delete;

    /**
     * Move Constructor
     */
    IncompressiblePotentialFlowWakeElement(IncompressiblePotentialFlowWakeElement&& rOther) = delete;

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowWakeElement() override{}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowWakeElement& operator=(IncompressiblePotentialFlowWakeElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

    /// Move operator.
    IncompressiblePotentialFlowWakeElement& operator=(IncompressiblePotentialFlowWakeElement&& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

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

    void ComputeLHSGaussPointContribution(const double weight,
                                          Matrix& lhs,
                                          const ElementalData<NumNodes, Dim>& data) const;

    void AssignLocalSystemWakeElement(MatrixType& rLeftHandSideMatrix,
                                      Matrix& lhs_total,
                                      const ElementalData<NumNodes, Dim>& data) const;

    void AssignLocalSystemWakeNode(MatrixType& rLeftHandSideMatrix,
                                   Matrix& lhs_total,
                                   const ElementalData<NumNodes, Dim>& data,
                                   const unsigned int& row) const;

    void CheckWakeCondition() const;

    void ComputePotentialJump(ProcessInfo& rCurrentProcessInfo);

    void ComputeElementInternalEnergy();

    void GetPotential(Vector& split_element_values,
                      const array_1d<double, NumNodes>& distances) const;

    void GetPotentialOnUpperWakeElement(array_1d<double, NumNodes>& upper_phis,
                                        const array_1d<double, NumNodes>& distances) const;

    void GetPotentialOnLowerWakeElement(array_1d<double, NumNodes>& lower_phis,
                                        const array_1d<double, NumNodes>& distances) const;

    void ComputeUpperVelocity(array_1d<double, Dim>& velocity) const;

    void ComputeLowerVelocity(array_1d<double, Dim>& velocity) const;

    double ComputeUpperPressure(const ProcessInfo& rCurrentProcessInfo) const;

    double ComputeLowerPressure(const ProcessInfo& rCurrentProcessInfo) const;

    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

}; // Class IncompressiblePotentialFlowWakeElement

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_WAKE_ELEMENT_H  defined
