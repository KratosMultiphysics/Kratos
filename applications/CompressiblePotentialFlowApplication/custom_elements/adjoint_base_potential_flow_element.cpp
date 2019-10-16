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
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "compressible_potential_flow_element.h"
#include "embedded_incompressible_potential_flow_element.h"
#include "embedded_compressible_potential_flow_element.h"
#include "adjoint_base_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::IntegrationMethod AdjointBasePotentialFlowElement<TPrimalElement>::GetIntegrationMethod() const
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::Initialize()
    {
        mpPrimalElement->Initialize();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->Data() = this->Data();
        mpPrimalElement->Set(Flags(*this));
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {

    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType tmp;
        mpPrimalElement->CalculateLeftHandSide(tmp, rCurrentProcessInfo);
        rLeftHandSideMatrix = trans(tmp);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
            std::vector< array_1d<double,3> >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
    {
        mpPrimalElement->GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValuesVector(Vector& rValues, int Step)
    {
        KRATOS_TRY
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if (wake == 1) // wake element
        {
            if(rValues.size() != 2*NumNodes)
                rValues.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);
            GetValuesOnSplitElement(rValues,distances);

        }else{ // normal element
            if(rValues.size() != NumNodes)
                rValues.resize(NumNodes, false);

            if(kutta == 0){
                for(unsigned int i=0; i<NumNodes; i++)
                    rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            }else{
                for(unsigned int i=0; i<NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling CalculateSensitivityMatrix from adjoint potential flow base element." << std::endl;

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR << "Calling CalculateSensitivityMatrix from adjoint potential flow base element." << std::endl;

        KRATOS_CATCH("")
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if(wake == 0)//normal element
        {
            if (rResult.size() != NumNodes)
                rResult.resize(NumNodes, false);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                    else
                        rResult[i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
                }
            }
        }
        else//wake element
        {
            if (rResult.size() != 2*NumNodes)
                rResult.resize(2*NumNodes, false);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rResult[NumNodes+i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
                else
                    rResult[NumNodes+i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,0).EquationId();
            }
        }


    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        const AdjointBasePotentialFlowElement& r_this = *this;
        const int wake = r_this.GetValue(WAKE);
        const int kutta = r_this.GetValue(KUTTA);
        const auto& r_geometry = GetGeometry();

        if(wake == 0) //normal element
        {
            if (rElementalDofList.size() != NumNodes)
                rElementalDofList.resize(NumNodes);

            if(kutta == 0){
                for (unsigned int i = 0; i < NumNodes; i++)
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            }
            else{
                for (unsigned int i = 0; i < NumNodes; i++){
                    if (!r_geometry[i].GetValue(TRAILING_EDGE))
                        rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                    else
                        rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
                }
            }
        }
        else//wake element
        {
            if (rElementalDofList.size() != 2*NumNodes)
                rElementalDofList.resize(2*NumNodes);

            array_1d<double,NumNodes> distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(r_this);

            //positive part
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] > 0)
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }

            //negative part - sign is opposite to the previous case
            for (unsigned int i = 0; i < NumNodes; i++)
            {
                if(distances[i] < 0)
                    rElementalDofList[NumNodes+i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
                else
                    rElementalDofList[NumNodes+i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
            }
        }
    }

    template <class TPrimalElement>
    int AdjointBasePotentialFlowElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo)
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
    std::string AdjointBasePotentialFlowElement<TPrimalElement>::Info() const
    {
        std::stringstream buffer;
        buffer << "AdjointBasePotentialFlowElement #" << Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AdjointBasePotentialFlowElement #" << Id();
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::PrintData(std::ostream& rOStream) const
    {
        pGetGeometry()->PrintData(rOStream);
    }

    template <class TPrimalElement>
    Element::Pointer AdjointBasePotentialFlowElement<TPrimalElement>::pGetPrimalElement()
    {
        return mpPrimalElement;
    }

    /*PROTECTED*/

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::GetValuesOnSplitElement(Vector& split_element_values, const array_1d<double,NumNodes>& distances )
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
    void AdjointBasePotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointBasePotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointBasePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<EmbeddedIncompressiblePotentialFlowElement<2,3>>;
    template class AdjointBasePotentialFlowElement<EmbeddedCompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.

