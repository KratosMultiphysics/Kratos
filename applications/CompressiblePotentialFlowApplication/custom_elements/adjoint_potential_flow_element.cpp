//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on A. Geiser, M. Fusseder, I. Lopez and R. Rossi work
//
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "compressible_potential_flow_element.h"
#include "adjoint_potential_flow_element.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::Pointer AdjointPotentialFlowElement<TPrimalElement>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
          return Kratos::make_intrusive<AdjointPotentialFlowElement<TPrimalElement>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointPotentialFlowElement<TPrimalElement>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const 
    {
        KRATOS_TRY
            return Kratos::make_intrusive<AdjointPotentialFlowElement<TPrimalElement>>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointPotentialFlowElement<TPrimalElement>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const 
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointPotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::IntegrationMethod AdjointPotentialFlowElement<TPrimalElement>::GetIntegrationMethod() const 
    {
        return mpPrimalElement->GetIntegrationMethod();
    }
  
    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::Initialize() 
    {   
        mpPrimalElement->Initialize();
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->Data() = this->Data();
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
    {

    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) 
    {
        MatrixType tmp;
        mpPrimalElement->CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);                                    
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
    
    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo) 
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step) 
    {
        KRATOS_TRY
        const AdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if (wake == 1) // wake element
        {
            if(rValues.size() != 2*NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances;
            GetWakeDistances(distances);
            GetValuesOnSplitElement(rValues,distances);
            
        }else{ // normal element
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if(kutta == 0){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] =GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) 
    {
        KRATOS_TRY;
        const double delta = this->GetPerturbationSize();
        ProcessInfo process_info = rCurrentProcessInfo;

        Vector RHS;
        Vector RHS_perturbed;

        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);

        for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
            for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                if ((GetGeometry()[i_node].Is(SOLID)) && (!GetGeometry()[i_node].GetValue(TRAILING_EDGE))){
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

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) 
    {
        const AdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0)//normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
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

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) 
    {
        const AdjointPotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);

        if(wake == 0) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
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

    template <class TPrimalElement>
    int AdjointPotentialFlowElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) 
    {

        KRATOS_TRY

        int Check = mpPrimalElement -> Check(rCurrentProcessInfo);

        if (Check != 0)
        {
            return Check;
        }
        else
        {
            for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
            {
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,
                                                    this->GetGeometry()[i]);

                return Check;
            }
        }

        return 0;

        KRATOS_CATCH("");
    }


    /// Turn back information as a string.
    template <class TPrimalElement>
    std::string AdjointPotentialFlowElement<TPrimalElement>::Info() const 
    {
        std::stringstream buffer;
        buffer << "AdjointPotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const 
    {
        rOStream << "AdjointPotentialFlowElement #" << Id();
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::PrintData(std::ostream& rOStream) const 
    {
        pGetGeometry()->PrintData(rOStream);
    }

    template <class TPrimalElement>
    Element::Pointer AdjointPotentialFlowElement<TPrimalElement>::pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    /*PROTECTED*/

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
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

    /*PRIVATE*/

    template <class TPrimalElement>
    double AdjointPotentialFlowElement<TPrimalElement>::GetPerturbationSize()
    {
        const double delta = this->GetValue(SCALE_FACTOR);
        KRATOS_DEBUG_ERROR_IF_NOT(delta > 0) << "The perturbation size is not > 0!";
        return delta;
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const 
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointPotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer) 
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointPotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    template class AdjointPotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.

