//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Iñigo Lopez and Riccardo Rossi
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

    //Constructors.

    /// Default constuctor.
    /**
     * @param NewId Index number of the new element (optional)
     */
    IncompressiblePotentialFlowElement(
        IndexType NewId = 0){};

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowElement(
        IndexType NewId,
        const NodesArrayType &ThisNodes) : Element(NewId, ThisNodes){};

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry) : Element(NewId, pGeometry){};

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties){};

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowElement(
        IncompressiblePotentialFlowElement const &rOther){};

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowElement() override{};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowElement &operator=(
        IncompressiblePotentialFlowElement const &rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressiblePotentialFlowElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const &ThisNodes) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix (output)
     * @param rRightHandSideVector: the elemental right hand side (output)
     * @param rCurrentProcessInfo: the current process info instance (input)
     */
    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief EquationIdVector Returns the global system rows corresponding to each local row.
     * @param rResult rResult[i] is the global index of local row i (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void EquationIdVector(
        EquationIdVectorType &rResult,
        ProcessInfo &CurrentProcessInfo) override;

    /**
     * @brief GetDofList Returns a list of the element's Dofs.
     * @param rElementalDofList List of DOFs. (output)
     * @param rCurrentProcessInfo Current ProcessInfo instance. (input)
     */
    void GetDofList(
        DofsVectorType &rElementalDofList,
        ProcessInfo &CurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override
    {
        //TODO: improve speed
        Matrix tmp;
        CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(
        ProcessInfo &rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo &rCurrentProcessInfo) override
    {

        KRATOS_TRY

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "IncompressiblePotentialFlowElement found with Id 0 or negative", "")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on IncompressiblePotentialFlowElement -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0", "")
        }

        for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(POSITIVE_FACE_PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing variable POSITIVE_FACE_PRESSURE on node ", this->GetGeometry()[i].Id())
        }

        return 0;

        KRATOS_CATCH("");
    }

    void GetValueOnIntegrationPoints(const Variable<double> &rVariable,
                                     std::vector<double> &rValues,
                                     const ProcessInfo &rCurrentProcessInfo) override
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == PRESSURE)
        {
            double p = 0.0;
            p = ComputePressureUpper(rCurrentProcessInfo);
            rValues[0] = p;
        }
        else if (rVariable == PRESSURE_LOWER)
        {
            double p = 0.0;
            p = ComputePressureLower(rCurrentProcessInfo);
            rValues[0] = p;
        }
        else if (rVariable == THICKNESS)
        {
            if (this->Is(THERMAL))
                rValues[0] = 30.0;
            else if (this->Is(MODIFIED))
                rValues[0] = 20.0;
            else if (this->Is(MARKER))
                rValues[0] = 10.0;
            else
                rValues[0] = 0.0;
        }
    }

    void GetValueOnIntegrationPoints(const Variable<int> &rVariable,
                                     std::vector<int> &rValues,
                                     const ProcessInfo &rCurrentProcessInfo) override
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if(rVariable == TRAILING_EDGE)
            rValues[0] = this->GetValue(TRAILING_EDGE);
        else if(rVariable == KUTTA)
            rValues[0] = this->GetValue(KUTTA);
        else if(rVariable == ALL_TRAILING_EDGE)
            rValues[0] = this->GetValue(ALL_TRAILING_EDGE);
        else if(rVariable == ZERO_VELOCITY_CONDITION)
            rValues[0] = this->GetValue(ZERO_VELOCITY_CONDITION);
        else if(rVariable == TRAILING_EDGE_ELEMENT)
            rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
        else if(rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
            rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
    }

    void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3>> &rVariable,
                                     std::vector<array_1d<double, 3>> &rValues,
                                     const ProcessInfo &rCurrentProcessInfo) override
    {
        if (rValues.size() != 1)
            rValues.resize(1);

        if (rVariable == VELOCITY)
        {
            array_1d<double, 3> v(3, 0.0);
            array_1d<double, Dim> vaux;
            ComputeVelocityUpper(vaux);
            for (unsigned int k = 0; k < Dim; k++)
                v[k] = vaux[k];
            rValues[0] = v;
        }
        else if (rVariable == VELOCITY_LOWER)
        {
            array_1d<double, 3> v(3, 0.0);
            array_1d<double, Dim> vaux;
            ComputeVelocityLower(vaux);
            for (unsigned int k = 0; k < Dim; k++)
                v[k] = vaux[k];
            rValues[0] = v;
        }
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IncompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "IncompressiblePotentialFlowElement #" << Id();
    }

    /// Print object's data.

    void PrintData(std::ostream &rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

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
    void GetWakeDistances(array_1d<double, NumNodes> &distances);

    void EquationIdVectorNormalElement(EquationIdVectorType &rResult);

    void EquationIdVectorKuttaElement(EquationIdVectorType &rResult);

    void EquationIdVectorWakeElement(EquationIdVectorType &rResult);

    void GetDofListNormalElement(DofsVectorType &rElementalDofList);

    void GetDofListKuttaElement(DofsVectorType &rElementalDofList);

    void GetDofListWakeElement(DofsVectorType &rElementalDofList);

    void CalculateLocalSystemNormalElement(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector);

    void CalculateLocalSystemWakeElement(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector);

    void CalculateLocalSystemSubdividedElement(
        Matrix &lhs_positive,
        Matrix &lhs_negative);

    void ComputeLHSGaussPointContribution(
        const double weight,
        Matrix &lhs,
        const ElementalData<NumNodes, Dim> &data);

    void AssignLocalSystemSubdividedElement(
        MatrixType &rLeftHandSideMatrix,
        Matrix &lhs_positive,
        Matrix &lhs_negative,
        Matrix &lhs_total,
        const ElementalData<NumNodes, Dim> &data);

    void AssignLocalSystemWakeElement(
        MatrixType &rLeftHandSideMatrix,
        Matrix &lhs_total,
        const ElementalData<NumNodes, Dim> &data);

    void AssignLocalSystemWakeNode(
        MatrixType &rLeftHandSideMatrix,
        Matrix &lhs_total,
        const ElementalData<NumNodes, Dim> &data,
        unsigned int &row);

    void ComputeRHSGaussPointContribution(
        const double weight,
        Vector &rhs,
        const ElementalData<NumNodes, Dim> &data)
    {
        array_1d<double, Dim> grad = prod(trans(data.DN_DX), data.phis);
        noalias(rhs) -= weight * prod(data.DN_DX, grad);
    }

    void GetPotentialOnNormalElement(array_1d<double, NumNodes> &phis);
    
    void GetValuesOnSplitElement(Vector &split_element_values, const array_1d<double, NumNodes> &distances)
    {
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] < 0)
                split_element_values[NumNodes + i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                split_element_values[NumNodes + i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
        }
    }

    void ComputeVelocityUpper(array_1d<double, Dim> &velocity)
    {
        velocity.clear();

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->IsNot(MARKER) && active == true)
            ComputeVelocityNormalElement(velocity);
        else if (active == true && this->Is(MARKER))
            ComputeVelocityUpperWakeElement(velocity);
    }

    void ComputeVelocityLower(array_1d<double, Dim> &velocity)
    {
        velocity.clear();

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->IsNot(MARKER) && active == true)
            ComputeVelocityNormalElement(velocity);
        else if (active == true && this->Is(MARKER))
            ComputeVelocityLowerWakeElement(velocity);
    }

    void ComputeVelocityNormalElement(array_1d<double, Dim> &velocity)
    {
        ElementalData<NumNodes, Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        const IncompressiblePotentialFlowElement &r_this = *this;
        const int &kutta = r_this.GetValue(KUTTA);

        if(kutta == 0)
        {
            //gather nodal data
            for (unsigned int i = 0; i < NumNodes; i++)
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        }
        else
        {
            // array_1d<double, NumNodes> distances;
            // GetWakeDistances(distances);
            //taking only negative part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(!GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                else
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
            }
        }
        noalias(velocity) = prod(trans(data.DN_DX), data.phis);
    }

    void ComputeVelocityUpperWakeElement(array_1d<double, Dim> &velocity)
    {
        ElementalData<NumNodes, Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        array_1d<double, NumNodes> distances;
        GetWakeDistances(distances);

        //taking only positive part
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if (distances[i] > 0)
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
            else
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
        }

        noalias(velocity) = prod(trans(data.DN_DX), data.phis);
    }

    void ComputeVelocityLowerWakeElement(array_1d<double, Dim> &velocity)
    {
        ElementalData<NumNodes, Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        array_1d<double, NumNodes> distances;
        GetWakeDistances(distances);

        if (this->Is(ISOLATED))
        {
            //taking only negative part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if (distances[i] < 0)
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                else
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NODAL_H);
            }
        }
        else
        {
            //taking only negative part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if (distances[i] < 0)
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
                else
                    data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
            }
        }

        noalias(velocity) = prod(trans(data.DN_DX), data.phis);
    }

    void CheckWakeCondition();
    void ComputePotentialJump(ProcessInfo &rCurrentProcessInfo);
    void ComputeElementInternalEnergy();

    double ComputePressureUpper(const ProcessInfo &rCurrentProcessInfo)
    {
        double pressure = 0.0;

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (active && !this->Is(MARKER))
            pressure = ComputePressureNormalElement(rCurrentProcessInfo);
        else if (active == true && this->Is(MARKER))
            pressure = ComputePressureUpperWakeElement(rCurrentProcessInfo);

        return pressure;
    }

    double ComputePressureLower(const ProcessInfo &rCurrentProcessInfo)
    {
        double pressure = 0.0;

        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (active && !this->Is(MARKER))
            pressure = ComputePressureNormalElement(rCurrentProcessInfo);
        else if (active == true && this->Is(MARKER))
            pressure = ComputePressureLowerWakeElement(rCurrentProcessInfo);

        return pressure;
    }

    double ComputePressureNormalElement(const ProcessInfo &rCurrentProcessInfo)
    {
        double pressure = 0.0;
        const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

        KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
            << "Error on element -> " << this->Id() << "\n"
            << "vinfinity_norm2 must be larger than zero." << std::endl;

        array_1d<double, Dim> v;
        ComputeVelocityNormalElement(v);

        pressure = (vinfinity_norm2 - inner_prod(v, v)) / vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

        return pressure;
    }

    double ComputePressureUpperWakeElement(const ProcessInfo &rCurrentProcessInfo)
    {
        double pressure = 0.0;
        const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

        KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
            << "Error on element -> " << this->Id() << "\n"
            << "vinfinity_norm2 must be larger than zero." << std::endl;

        array_1d<double, Dim> v;
        ComputeVelocityUpperWakeElement(v);

        pressure = (vinfinity_norm2 - inner_prod(v, v)) / vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

        return pressure;
    }

    double ComputePressureLowerWakeElement(const ProcessInfo &rCurrentProcessInfo)
    {
        double pressure = 0.0;
        const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
        const double vinfinity_norm2 = inner_prod(vinfinity, vinfinity);

        KRATOS_ERROR_IF(vinfinity_norm2 < std::numeric_limits<double>::epsilon())
            << "Error on element -> " << this->Id() << "\n"
            << "vinfinity_norm2 must be larger than zero." << std::endl;

        array_1d<double, Dim> v;
        ComputeVelocityLowerWakeElement(v);

        pressure = (vinfinity_norm2 - inner_prod(v, v)) / vinfinity_norm2; //0.5*(norm_2(vinfinity) - norm_2(v));

        return pressure;
    }

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

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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
