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

#if !defined(KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED)
#define KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED

// #define SYMMETRIC_CONSTRAINT_APPLICATION

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/element.h"
#include "includes/kratos_flags.h"
#include "compressible_potential_flow_application_variables.h"
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
class IncompressiblePotentialFlowElement : public Element
{
public:
    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double, TNumNodes> phis, distances;
        double rho;
        double vol;

        bounded_matrix<double, TNumNodes, TDim> DN_DX;
        array_1d<double, TNumNodes> N;
    };

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of IncompressiblePotentialFlowElement
    KRATOS_CLASS_POINTER_DEFINITION(IncompressiblePotentialFlowElement);

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    IncompressiblePotentialFlowElement(IndexType NewId = 0){};

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes){};

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry){};

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowElement(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties){};

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowElement(IncompressiblePotentialFlowElement const& rOther){};

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowElement() override{};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowElement& operator=(IncompressiblePotentialFlowElement const& rOther)
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
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    void GetWakeDistances(array_1d<double, NumNodes>& distances);

    void GetEquationIdVectorNormalElement(EquationIdVectorType& rResult);

    void GetEquationIdVectorKuttaElement(EquationIdVectorType& rResult);

    void GetEquationIdVectorWakeElement(EquationIdVectorType& rResult);

    void GetDofListNormalElement(DofsVectorType& rElementalDofList);

    void GetDofListKuttaElement(DofsVectorType& rElementalDofList);

    void GetDofListWakeElement(DofsVectorType& rElementalDofList);

    void CalculateLocalSystemNormalElement(MatrixType& rLeftHandSideMatrix,
                                           VectorType& rRightHandSideVector);

    void CalculateLocalSystemWakeElement(MatrixType& rLeftHandSideMatrix,
                                         VectorType& rRightHandSideVector);

    void CalculateLocalSystemSubdividedElement(Matrix& lhs_positive, Matrix& lhs_negative);

    void ComputeLHSGaussPointContribution(const double weight,
                                          Matrix& lhs,
                                          const ElementalData<NumNodes, Dim>& data);

    void AssignLocalSystemSubdividedElement(MatrixType& rLeftHandSideMatrix,
                                            Matrix& lhs_positive,
                                            Matrix& lhs_negative,
                                            Matrix& lhs_total,
                                            const ElementalData<NumNodes, Dim>& data);

    void AssignLocalSystemWakeElement(MatrixType& rLeftHandSideMatrix,
                                      Matrix& lhs_total,
                                      const ElementalData<NumNodes, Dim>& data);

    void AssignLocalSystemWakeNode(MatrixType& rLeftHandSideMatrix,
                                   Matrix& lhs_total,
                                   const ElementalData<NumNodes, Dim>& data,
                                   unsigned int& row);

    void CheckWakeCondition();

    void ComputePotentialJump(ProcessInfo& rCurrentProcessInfo);

    void ComputeElementInternalEnergy();

    void GetPotentialOnNormalElement(array_1d<double, NumNodes>& phis);

    void GetPotentialOnWakeElement(Vector& split_element_values,
                                   const array_1d<double, NumNodes>& distances);

    void GetPotentialOnUpperWakeElement(array_1d<double, NumNodes>& upper_phis,
                                        const array_1d<double, NumNodes>& distances);

    void GetPotentialOnLowerWakeElement(array_1d<double, NumNodes>& lower_phis,
                                        const array_1d<double, NumNodes>& distances);

    void ComputeVelocityUpper(array_1d<double, Dim>& velocity);

    void ComputeVelocityLower(array_1d<double, Dim>& velocity);

    void ComputeVelocityNormalElement(array_1d<double, Dim>& velocity);

    void ComputeVelocityUpperWakeElement(array_1d<double, Dim>& velocity);

    void ComputeVelocityLowerWakeElement(array_1d<double, Dim>& velocity);

    double ComputePressureUpper(const ProcessInfo& rCurrentProcessInfo);

    double ComputePressureLower(const ProcessInfo& rCurrentProcessInfo);

    double ComputePressureNormalElement(const ProcessInfo& rCurrentProcessInfo);

    double ComputePressureUpperWakeElement(const ProcessInfo& rCurrentProcessInfo);

    double ComputePressureLowerWakeElement(const ProcessInfo& rCurrentProcessInfo);

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

}; // Class IncompressiblePotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_INCOMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
