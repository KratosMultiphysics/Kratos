//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_INCOMPRESSIBLE_ADJOINT_POTENTIAL_FLOW_ELEMENT_H_INCLUDED )
#define KRATOS_INCOMPRESSIBLE_ADJOINT_POTENTIAL_FLOW_ELEMENT_H_INCLUDED


// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "utilities/enrichment_utilities.h"
namespace Kratos
{


template< int Dim, int NumNodes >
class IncompressibleAdjointPotentialFlowElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressibleAdjointPotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressibleAdjointPotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    IncompressibleAdjointPotentialFlowElement(IndexType NewId = 0) {};

    IncompressibleAdjointPotentialFlowElement(Element::Pointer pPrimalElement)
                    : Element(pPrimalElement->Id(), pPrimalElement->pGetGeometry(), pPrimalElement->pGetProperties())
                    , mpPrimalElement(pPrimalElement)
    {
    };

    /**
     * Constructor using an array of nodes
     */
    IncompressibleAdjointPotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes):Element(NewId, ThisNodes) {};

    /**
     * Constructor using Geometry
     */
    IncompressibleAdjointPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry):Element(NewId, pGeometry) {};

    /**
     * Constructor using Properties
     */
    IncompressibleAdjointPotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):Element(NewId, pGeometry, pProperties) {};

    /**
     * Copy Constructor
     */
    IncompressibleAdjointPotentialFlowElement(IncompressibleAdjointPotentialFlowElement const& rOther) {};

    /**
     * Destructor
     */
    ~IncompressibleAdjointPotentialFlowElement() override {};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressibleAdjointPotentialFlowElement & operator=(IncompressibleAdjointPotentialFlowElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    IntegrationMethod GetIntegrationMethod() const override;

    void Initialize() override;

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;
    
    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step=0) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

   Element::Pointer pGetPrimalElement();


protected:


    Element::Pointer mpPrimalElement;

    void GetWakeDistances(array_1d<double,NumNodes>& distances);

    void GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances);


private:

    double GetPerturbationSize();

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // Class IncompressibleAdjointPotentialFlowElement


} // namespace Kratos.

#endif 
