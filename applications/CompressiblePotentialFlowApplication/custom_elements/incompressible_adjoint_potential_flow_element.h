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
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
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

template< int Dim, int NumNodes >
class IncompressibleAdjointPotentialFlowElement : public Element
{
public:

    template <unsigned int TNumNodes, unsigned int TDim>
    struct ElementalData
    {
        array_1d<double,TNumNodes> phis, distances;
        double rho;
        double vol;

        bounded_matrix<double, TNumNodes, TDim > DN_DX;
        array_1d<double, TNumNodes > N;
    };



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

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
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
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, pGeom, pProperties));
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Element::Pointer(new IncompressibleAdjointPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

     IntegrationMethod GetIntegrationMethod() const override
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

  

    void Initialize() override
    {   
        mpPrimalElement->Initialize();
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalElement->ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mpPrimalElement->CleanMemory();
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->Data() = this->Data();
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        //avoiding check wake condition
        //  mpPrimalElement->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrix,
                                               rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(std::vector< MatrixType >& rLeftHandSideMatrices,
					const std::vector< Variable< MatrixType > >& rLHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrices,
                                               rLHSVariables,
                                               rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
					const std::vector< Variable< VectorType > >& rRHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVectors,
                                                rRHSVariables,
                                                rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					        ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							    rRightHandSideVector,
							    rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateMassMatrix(rMassMatrix,rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                Variable<double >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                    rRHSVariable,
                                    rDestinationVariable,
                                    rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                Variable<array_1d<double,3> >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                rRHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void AddExplicitContribution(const MatrixType& rLHSMatrix,
                                const Variable<MatrixType>& rLHSVariable,
                                Variable<Matrix>& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rLHSMatrix,
                                rLHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
    
    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    void GetValuesVector(Vector& rValues, int Step=0) override
    {
        KRATOS_TRY

        if (this->Is(MARKER)){
            if(rValues.size() != 2*NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);
            GetValuesOnSplitElement(rValues,distances);
            
        }else{
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if (this->IsNot(STRUCTURE)){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] =GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (GetGeometry()[i].IsNot(STRUCTURE))
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        // KRATOS_WATCH("entering Sensitivity calculation")
        // Get perturbation size
        const double delta = this->GetPerturbationSize(rDesignVariable);
        ProcessInfo process_info = rCurrentProcessInfo;

        Vector RHS;
        Vector RHS_perturbed;
        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes || rOutput.size2() != RHS.size())
            rOutput.resize(NumNodes, RHS.size(), false);

        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            double distance_value = GetGeometry()[i_node].GetSolutionStepValue(LEVEL_SET_DISTANCE);
            //perturbate distance
            pGetPrimalElement()->GetGeometry()[i_node].GetSolutionStepValue(LEVEL_SET_DISTANCE) = distance_value+delta;
            //compute perturbated RHS
            pGetPrimalElement()->CalculateRightHandSide(RHS_perturbed, process_info);
            //recover distance value
            pGetPrimalElement()->GetGeometry()[i_node].GetSolutionStepValue(LEVEL_SET_DISTANCE) = distance_value;
            for (unsigned int i_dof =0;i_dof<RHS.size();i_dof++)
                rOutput(i_node,i_dof) = (RHS_perturbed(i_dof)-RHS(i_dof))/delta;                
        }
        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;
        const double delta = 1e-6;//this->GetPerturbationSize(rDesignVariable);
        ProcessInfo process_info = rCurrentProcessInfo;

        Vector RHS;
        Vector RHS_perturbed;

        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);
        bool kutta_node_spotted=false;
        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            if(GetGeometry()[i_node].Is(STRUCTURE)){
                kutta_node_spotted=true;
                break;
            }
        }

        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                if (GetGeometry()[i_node].Is(SOLID)){
                    pGetPrimalElement()->GetGeometry()[i_node].GetInitialPosition()[i_dim] += delta;
                    pGetPrimalElement()->GetGeometry()[i_node].Coordinates()[i_dim] += delta;

                    // compute LHS after perturbation
                    pGetPrimalElement()->CalculateRightHandSide(RHS_perturbed, process_info);

                    //compute derivative of RHS w.r.t. design variable with finite differences
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = (RHS_perturbed[i] - RHS[i]) / delta;

                    // unperturb the design variable
                    pGetPrimalElement()->GetGeometry()[i_node].GetInitialPosition()[i_dim] -= delta;
                    pGetPrimalElement()->GetGeometry()[i_node].Coordinates()[i_dim] -= delta;
                }else{
                    for(unsigned int i = 0; i < RHS.size(); ++i)
                        rOutput( (i_dim + i_node*Dim), i) = 0.0;
                }
            }
        }
        if (kutta_node_spotted){ //remove kutta node lines
            for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
                if(GetGeometry()[i_node].Is(STRUCTURE)){
                    for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                        for(unsigned int i = 0; i < RHS.size(); ++i)
                                rOutput( (i_dim + i_node*Dim), i) = 0.0;
                    }
                }
            }  
        }
        KRATOS_CATCH("")
    }
    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            if(this->IsNot(STRUCTURE)){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (GetGeometry()[i].IsNot(STRUCTURE))
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                    else
                        rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                }
            }
        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }
        }


    }


    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
    {
        if(this->IsNot(MARKER)) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(this->IsNot(STRUCTURE)){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (GetGeometry()[i].IsNot(STRUCTURE))
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }


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
    
    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {

        KRATOS_TRY

        if (this->Id() < 1)
        {
            KRATOS_THROW_ERROR(std::logic_error, "IncompressibleAdjointPotentialFlowElement found with Id 0 or negative", "")
        }

        if (this->GetGeometry().Area() <= 0)
        {
            std::cout << "error on IncompressibleAdjointPotentialFlowElement -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than or equal to 0", "")
        }

        for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
        {
            if (this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL) == false)
                KRATOS_THROW_ERROR(std::invalid_argument, "missing variable ADJOINT_VELOCITY_POTENTIAL on node ", this->GetGeometry()[i].Id())
        }

        return 0;

        KRATOS_CATCH("");
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
        buffer << "IncompressibleAdjointPotentialFlowElement #" << Id();
        return buffer.str();
    }

/// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IncompressibleAdjointPotentialFlowElement #" << Id();
    }

/// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    Element::Pointer pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    Element::Pointer mpPrimalElement;
    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{
    void GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

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

    void GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
    {

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] > 0)
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }

        //negative part - sign is opposite to the previous case
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            if(distances[i] < 0)
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                split_element_values[NumNodes+i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
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

    double GetPerturbationSize(const Variable<double>& rDesignVariable)
    {
        // const double correction_factor = this->GetPerturbationSizeModificationFactor(rDesignVariable);
        const double delta = 1e-9;//this->GetValue(PERTURBATION_SIZE);
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
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

}; // Class IncompressibleAdjointPotentialFlowElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_COMPRESSIBLE_POTENTIAL_FLOW_ELEMENT_H_INCLUDED  defined
