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
    IncompressiblePotentialFlowElement(IndexType NewId = 0){};

    /**
     * Constructor using an array of nodes
     */
    IncompressiblePotentialFlowElement(IndexType NewId, const NodesArrayType &ThisNodes) : Element(NewId, ThisNodes){};

    /**
     * Constructor using Geometry
     */
    IncompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry){};

    /**
     * Constructor using Properties
     */
    IncompressiblePotentialFlowElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties){};

    /**
     * Copy Constructor
     */
    IncompressiblePotentialFlowElement(IncompressiblePotentialFlowElement const &rOther){};

    /**
     * Destructor
     */
    ~IncompressiblePotentialFlowElement() override{};

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IncompressiblePotentialFlowElement &operator=(IncompressiblePotentialFlowElement const &rOther)
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
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes, PropertiesType::Pointer pProperties) const override
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
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
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
    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    /**
     * @brief EquationIdVector Returns the global system rows corresponding to each local row.
     * @param rResult rResult[i] is the global index of local row i (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void EquationIdVector(EquationIdVectorType &rResult, ProcessInfo &CurrentProcessInfo) override;

    /**
     * @brief GetDofList Returns a list of the element's Dofs.
     * @param rElementalDofList List of DOFs. (output)
     * @param rCurrentProcessInfo Current ProcessInfo instance. (input)
     */
    void GetDofList(DofsVectorType &rElementalDofList, ProcessInfo &CurrentProcessInfo) override;

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        ProcessInfo &rCurrentProcessInfo) override
    {
        ElementalData<NumNodes, Dim> data;

        //calculate shape functions
        GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

        //gather nodal data
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
        }

        if (this->IsNot(MARKER)) //normal element (non-wake) - eventually an embedded
        {
            if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
                rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
            if (rRightHandSideVector.size() != NumNodes)
                rRightHandSideVector.resize(NumNodes, false);
            rLeftHandSideMatrix.clear();

            ComputeLHSGaussPointContribution(data.vol, rLeftHandSideMatrix, data);

            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, data.phis);
        }
        else //it is a wake element
        {
            GetWakeDistances(data.distances);

            //note that the lhs and rhs have double the size!!
            if (rLeftHandSideMatrix.size1() != 2 * NumNodes || rLeftHandSideMatrix.size2() != 2 * NumNodes)
                rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
            if (rRightHandSideVector.size() != 2 * NumNodes)
                rRightHandSideVector.resize(2 * NumNodes, false);
            rLeftHandSideMatrix.clear();

            if (this->Is(STRUCTURE))
            {
                //subdivide the element
                constexpr unsigned int nvolumes = 3 * (Dim - 1);
                bounded_matrix<double, NumNodes, Dim> Points;
                array_1d<double, nvolumes> Volumes;
                bounded_matrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
                array_1d<double, nvolumes> PartitionsSign;
                std::vector<Matrix> GradientsValue(nvolumes);
                bounded_matrix<double, nvolumes, 2> NEnriched;

                for (unsigned int i = 0; i < GradientsValue.size(); ++i)
                    GradientsValue[i].resize(2, Dim, false);
                for (unsigned int i = 0; i < NumNodes; ++i)
                {
                    const array_1d<double, 3> &coords = GetGeometry()[i].Coordinates();
                    for (unsigned int k = 0; k < Dim; ++k)
                    {
                        Points(i, k) = coords[k];
                    }
                }

                const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(Points,
                                                                                                       data.DN_DX,
                                                                                                       data.distances,
                                                                                                       Volumes,
                                                                                                       GPShapeFunctionValues,
                                                                                                       PartitionsSign,
                                                                                                       GradientsValue,
                                                                                                       NEnriched);

                //compute the lhs and rhs that would correspond to it being divided
                Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
                Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

                for (unsigned int i = 0; i < nsubdivisions; ++i)
                {
                    if (PartitionsSign[i] > 0)
                        ComputeLHSGaussPointContribution(Volumes[i], lhs_positive, data);
                    else
                        ComputeLHSGaussPointContribution(Volumes[i], lhs_negative, data);
                }

                Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);
                ComputeLHSGaussPointContribution(data.vol, lhs_total, data);

                for (unsigned int i = 0; i < NumNodes; ++i)
                {
                    //No contribution to the TE node //and to extra dofs
                    if(GetGeometry()[i].FastGetSolutionStepValue(TRAILING_EDGE))
                    {
                        for (unsigned int j = 0; j < NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                            rLeftHandSideMatrix(i, j + NumNodes) = 0.0;

                            rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);
                            rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                        }
                    }
                    else
                    {
                        for (unsigned int j = 0; j < NumNodes; ++j)
                        {
                            rLeftHandSideMatrix(i, j) = lhs_total(i, j);
                            rLeftHandSideMatrix(i, j + NumNodes) = 0.0;

                            rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total(i, j);
                            rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                        }

                        //Applying wake condition
                        if (data.distances[i] < 0.0)
                        {
                            for (unsigned int j = 0; j < NumNodes; ++j)
                                rLeftHandSideMatrix(i, j + NumNodes) = -lhs_total(i, j);
                        }
                        else if (data.distances[i] > 0.0)
                        {
                            for (unsigned int j = 0; j < NumNodes; ++j)
                                rLeftHandSideMatrix(i + NumNodes, j) = -lhs_total(i, j);
                        }
                    }

                }
            }
            else
            {
                Matrix lhs_total = ZeroMatrix(NumNodes, NumNodes);
                ComputeLHSGaussPointContribution(data.vol, lhs_total, data);

                //Looping over rows
                for (unsigned int i = 0; i < NumNodes; ++i)
                { //Looping over columngs
                    for (unsigned int j = 0; j < NumNodes; ++j)
                    { //Filling the diagonal blocks (i.e. decoupling upper and lower fields)
                        rLeftHandSideMatrix(i, j) = lhs_total(i, j);
                        rLeftHandSideMatrix(i, j + NumNodes) = 0.0;

                        rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total(i, j);
                        rLeftHandSideMatrix(i + NumNodes, j) = 0.0;
                    }

                    if (data.distances[i] < 0.0 && !GetGeometry()[i].FastGetSolutionStepValue(DEACTIVATED_WAKE))
                    {                            //side1  -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs and not on the airfoil nodes
                        //Marking nodes where the wake constraint is applied
                        GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) = 30.0;
                        for (unsigned int j = 0; j < NumNodes; ++j)
                            rLeftHandSideMatrix(i, j + NumNodes) = -lhs_total(i, j);
                    }
                    else if (data.distances[i] > 0.0)// && !GetGeometry()[i].FastGetSolutionStepValue(DEACTIVATED_WAKE))
                    { //side2 -assign constraint only on the NEGATIVE_FACE_PRESSURE dofs and not on the airfoil nodes
                        GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) = 30.0;
                        for (unsigned int j = 0; j < NumNodes; ++j)
                            rLeftHandSideMatrix(i + NumNodes, j) = -lhs_total(i, j);
                    }
                }
            }
            
            Vector split_element_values(NumNodes * 2);
            GetValuesOnSplitElement(split_element_values, data.distances);
            noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
        }
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo) override
    {
        //TODO: improve speed
        Matrix tmp;
        CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo) override
    {
        bool active = true;
        if ((this)->IsDefined(ACTIVE))
            active = (this)->Is(ACTIVE);

        if (this->Is(MARKER) && active == true)
        {
            CheckWakeCondition();

            const array_1d<double, 3> vinfinity = rCurrentProcessInfo[VELOCITY_INFINITY];
            const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

            array_1d<double, NumNodes> distances;
            GetWakeDistances(distances);

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if (distances[i] > 0)
                {
                    GetGeometry()[i].GetSolutionStepValue(POTENTIAL_JUMP) = 2.0 / vinfinity_norm * (GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) - GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE));
                }
                else
                {
                    GetGeometry()[i].GetSolutionStepValue(POTENTIAL_JUMP) = 2.0 / vinfinity_norm * (GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) - GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE));
                }
            }
        }

        //Compute element internal energy
        VectorType rRightHandSideVector;
        MatrixType rLeftHandSideMatrix;
        this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

        ElementalData<NumNodes, Dim> data;

        double internal_energy = 0.0;

        if (this->IsNot(MARKER)) //normal element (non-wake) - eventually an embedded
        {
            VectorType tmp;
            tmp.resize(NumNodes, false);

            //gather nodal data
            for (unsigned int i = 0; i < NumNodes; i++)
                data.phis[i] = GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);

            noalias(tmp) = prod(rLeftHandSideMatrix, data.phis);
            internal_energy = 0.5 * inner_prod(tmp, data.phis);
        }
        else
        {
            VectorType tmp;
            tmp.resize(NumNodes * 2, false);

            GetWakeDistances(data.distances);
            Vector split_element_values(NumNodes * 2);
            GetValuesOnSplitElement(split_element_values, data.distances);

            noalias(tmp) = prod(rLeftHandSideMatrix, split_element_values);
            internal_energy = 0.5 * inner_prod(rRightHandSideVector, split_element_values);
        }
        this->SetValue(INTERNAL_ENERGY, internal_energy);
    }

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

    void ComputeLHSGaussPointContribution(
        const double weight,
        Matrix &lhs,
        const ElementalData<NumNodes, Dim> &data)
    {
        noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
    }

    void ComputeRHSGaussPointContribution(
        const double weight,
        Vector &rhs,
        const ElementalData<NumNodes, Dim> &data)
    {
        array_1d<double, Dim> grad = prod(trans(data.DN_DX), data.phis);
        noalias(rhs) -= weight * prod(data.DN_DX, grad);
    }

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

    void CheckWakeCondition()
    {
        array_1d<double, Dim> upper_wake_velocity;
        ComputeVelocityUpperWakeElement(upper_wake_velocity);
        const double vupnorm = inner_prod(upper_wake_velocity, upper_wake_velocity);

        array_1d<double, Dim> lower_wake_velocity;
        ComputeVelocityLowerWakeElement(lower_wake_velocity);
        const double vlownorm = inner_prod(lower_wake_velocity, lower_wake_velocity);

        if (std::abs(vupnorm - vlownorm) > 0.1)
            std::cout << "WAKE CONDITION NOT FULFILLED IN ELEMENT # " << this->Id() << std::endl;
    }

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
