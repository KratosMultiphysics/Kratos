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
//  Main authors:    Iñigo Lopez based on M. Nuñez, A. Geiser, M. Fusseder and R. Rossi work
//
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "adjoint_analytical_incompressible_potential_flow_element.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
          return Kratos::make_intrusive<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
            return Kratos::make_intrusive<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointAnalyticalIncompressiblePotentialFlowElement(NewId, GetGeometry().Create(ThisNodes), pGetProperties()));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::IntegrationMethod AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetIntegrationMethod() const
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Initialize()
    {
        mpPrimalElement->Initialize();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->Data() = this->Data();
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {

    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType tmp;
        mpPrimalElement->CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step)
    {
        KRATOS_TRY
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
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
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        ProcessInfo process_info = rCurrentProcessInfo;
        Vector RHS;
        pGetPrimalElement()->CalculateRightHandSide(RHS, process_info);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);
        rOutput.clear();

        const int wake = pGetPrimalElement()->GetValue(WAKE);

        if (wake == 0) // Normal element (non-wake) - eventually an embedded
        {
            BoundedMatrix<double, NumNodes, Dim> x;
            for(unsigned int i_node = 0; i_node < NumNodes; i_node++){
                    x( i_node , 0 ) = GetGeometry()[i_node].X();
                    x( i_node , 1 ) = GetGeometry()[i_node].Y();
            }
            auto p = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*pGetPrimalElement());
            BoundedMatrix<double, Dim*NumNodes, NumNodes> test = ZeroMatrix(Dim*NumNodes, NumNodes);

            const double crOutput0 =             x(0,0) - x(1,0);
            const double crOutput1 =             -x(2,1);
            const double crOutput2 =             crOutput1 + x(0,1);
            const double crOutput3 =             -x(2,0);
            const double crOutput4 =             crOutput3 + x(0,0);
            const double crOutput5 =             x(0,1) - x(1,1);
            const double crOutput6 =             crOutput0*crOutput2 - crOutput4*crOutput5;
            const double crOutput7 =             pow(crOutput6, -2);
            const double crOutput8 =             0.5*crOutput7;
            const double crOutput9 =             crOutput3 + x(1,0);
            const double crOutput10 =             -p[2];
            const double crOutput11 =             crOutput6*(crOutput10 + p[1]);
            const double crOutput12 =             crOutput1 + x(1,1);
            const double crOutput13 =             crOutput0*crOutput9 + crOutput12*crOutput5;
            const double crOutput14 =             crOutput12*crOutput2 + crOutput4*crOutput9;
            const double crOutput15 =             crOutput13*p[2] - crOutput14*p[1] + p[0]*(pow(crOutput12, 2) + pow(crOutput9, 2));
            const double crOutput16 =             crOutput9*p[0];
            const double crOutput17 =             0.5*crOutput16;
            const double crOutput18 =             crOutput4*p[1];
            const double crOutput19 =             0.5*p[2];
            const double crOutput20 =             -2*x(0,0) + x(1,0) + x(2,0);
            const double crOutput21 =             -0.5*x(2,1);
            const double crOutput22 =             crOutput21 + 0.5*x(1,1);
            const double crOutput23 =             crOutput0*crOutput4 + crOutput2*crOutput5;
            const double crOutput24 =             crOutput14*p[0] + crOutput23*p[2] - p[1]*(pow(crOutput2, 2) + pow(crOutput4, 2));
            const double crOutput25 =             crOutput0*p[2];
            const double crOutput26 =             0.5*p[1];
            const double crOutput27 =             crOutput13*p[0] - crOutput23*p[1] + p[2]*(pow(crOutput0, 2) + pow(crOutput5, 2));
            const double crOutput28 =             crOutput12*p[0];
            const double crOutput29 =             0.5*crOutput28;
            const double crOutput30 =             crOutput2*p[1];
            const double crOutput31 =             -2*x(0,1) + x(1,1) + x(2,1);
            const double crOutput32 =             -0.5*x(2,0);
            const double crOutput33 =             crOutput32 + 0.5*x(1,0);
            const double crOutput34 =             crOutput5*p[2];
            const double crOutput35 =             0.5*crOutput18;
            const double crOutput36 =             x(0,0) - 2*x(1,0) + x(2,0);
            const double crOutput37 =             0.5*x(0,1);
            const double crOutput38 =             crOutput21 + crOutput37;
            const double crOutput39 =             crOutput6*(crOutput10 + p[0]);
            const double crOutput40 =             0.5*p[0];
            const double crOutput41 =             0.5*crOutput30;
            const double crOutput42 =             x(0,1) - 2*x(1,1) + x(2,1);
            const double crOutput43 =             0.5*x(0,0);
            const double crOutput44 =             crOutput32 + crOutput43;
            const double crOutput45 =             0.5*crOutput25;
            const double crOutput46 =             x(0,0) + x(1,0) - 2*x(2,0);
            const double crOutput47 =             crOutput37 - 0.5*x(1,1);
            const double crOutput48 =             crOutput6*(p[0] - p[1]);
            const double crOutput49 =             0.5*crOutput34;
            const double crOutput50 =             x(0,1) + x(1,1) - 2*x(2,1);
            const double crOutput51 =             crOutput43 - 0.5*x(1,0);
            rOutput(0,0)=crOutput8*(crOutput11*crOutput9 + crOutput12*crOutput15);
            rOutput(0,1)=-crOutput7*(crOutput22*crOutput24 + crOutput6*(-crOutput17 + 1.0*crOutput18 + crOutput19*crOutput20));
            rOutput(0,2)=crOutput7*(crOutput22*crOutput27 - crOutput6*(crOutput17 + crOutput20*crOutput26 + 1.0*crOutput25));
            rOutput(1,0)=crOutput8*(crOutput11*crOutput12 - crOutput15*crOutput9);
            rOutput(1,1)=crOutput7*(crOutput24*crOutput33 - crOutput6*(crOutput19*crOutput31 - crOutput29 + 1.0*crOutput30));
            rOutput(1,2)=-crOutput7*(crOutput27*crOutput33 + crOutput6*(crOutput26*crOutput31 + crOutput29 + 1.0*crOutput34));
            rOutput(2,0)=-crOutput7*(crOutput15*crOutput38 + crOutput6*(1.0*crOutput16 + crOutput19*crOutput36 - crOutput35));
            rOutput(2,1)=crOutput8*(crOutput2*crOutput24 + crOutput39*crOutput4);
            rOutput(2,2)=-crOutput7*(crOutput27*crOutput38 + crOutput6*(-crOutput25 + crOutput35 + crOutput36*crOutput40));
            rOutput(3,0)=crOutput7*(crOutput15*crOutput44 - crOutput6*(crOutput19*crOutput42 + 1.0*crOutput28 - crOutput41));
            rOutput(3,1)=crOutput8*(crOutput2*crOutput39 - crOutput24*crOutput4);
            rOutput(3,2)=crOutput7*(crOutput27*crOutput44 - crOutput6*(-crOutput34 + crOutput40*crOutput42 + crOutput41));
            rOutput(4,0)=crOutput7*(crOutput15*crOutput47 + crOutput6*(crOutput16 - crOutput26*crOutput46 + crOutput45));
            rOutput(4,1)=-crOutput7*(crOutput24*crOutput47 + crOutput6*(-crOutput18 + crOutput40*crOutput46 + crOutput45));
            rOutput(4,2)=crOutput8*(crOutput0*crOutput48 + crOutput27*crOutput5);
            rOutput(5,0)=-crOutput7*(crOutput15*crOutput51 - crOutput6*(-crOutput26*crOutput50 + crOutput28 + crOutput49));
            rOutput(5,1)=crOutput7*(crOutput24*crOutput51 - crOutput6*(-crOutput30 + crOutput40*crOutput50 + crOutput49));
            rOutput(5,2)=crOutput8*(-crOutput0*crOutput27 + crOutput48*crOutput5);

            for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
                for(unsigned int i_dim = 0; i_dim<Dim; i_dim++){
                    if ((!GetGeometry()[i_node].Is(SOLID)) || (GetGeometry()[i_node].GetValue(TRAILING_EDGE))){
                        for(unsigned int i = 0; i < RHS.size(); ++i)
                            rOutput( (i_dim + i_node*Dim), i) = 0.0;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
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
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        const AdjointAnalyticalIncompressiblePotentialFlowElement& r_this = *this;
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
    int AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo)
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
    std::string AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Info() const
    {
        std::stringstream buffer;
        buffer << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << Id();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::PrintData(std::ostream& rOStream) const
    {
        pGetGeometry()->PrintData(rOStream);
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    /*PROTECTED*/

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetWakeDistances(array_1d<double,NumNodes>& distances)
    {
        noalias(distances) = GetValue(ELEMENTAL_DISTANCES);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
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
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointAnalyticalIncompressiblePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    //template class AdjointAnalyticalIncompressiblePotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.

