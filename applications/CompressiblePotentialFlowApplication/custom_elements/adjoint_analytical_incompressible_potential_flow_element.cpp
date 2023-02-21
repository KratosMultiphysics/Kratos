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
//  Main authors:    Inigo Lopez and Marc Nunez, based on A. Geiser, M. Fusseder and R. Rossi work
//
#include "compressible_potential_flow_application_variables.h"
#include "incompressible_potential_flow_element.h"
#include "adjoint_analytical_incompressible_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
          return Kratos::make_intrusive<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Create(IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
            return Kratos::make_intrusive<AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>>(NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    Element::Pointer AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
    {
        KRATOS_TRY
        return Element::Pointer(new AdjointAnalyticalIncompressiblePotentialFlowElement(NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties()));
        KRATOS_CATCH("");
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                            Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        Vector RHS;
        auto pPrimalElement = this->pGetPrimalElement();
        pPrimalElement->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        if (rOutput.size1() != NumNodes)
            rOutput.resize(Dim*NumNodes, RHS.size(), false);
        rOutput.clear();

        const int wake = pPrimalElement->GetValue(WAKE);

        if (wake == 0) // Normal element (non-wake) - eventually an embedded
        {
            BoundedMatrix<double, NumNodes, Dim> x;
            for(unsigned int i_node = 0; i_node < NumNodes; i_node++){
                    x( i_node , 0 ) = pPrimalElement->GetGeometry()[i_node].X();
                    x( i_node , 1 ) = pPrimalElement->GetGeometry()[i_node].Y();
            }
            auto p = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*pPrimalElement);
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
                    if ((!pPrimalElement->GetGeometry()[i_node].Is(SOLID)) || (pPrimalElement->GetGeometry()[i_node].GetValue(TRAILING_EDGE))){
                        for(unsigned int i = 0; i < RHS.size(); ++i)
                            rOutput( (i_dim + i_node*Dim), i) = 0.0;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    /// Turn back information as a string.
    template <class TPrimalElement>
    std::string AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::Info() const
    {
        std::stringstream buffer;
        buffer << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << this->Id();
        return buffer.str();
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "AdjointAnalyticalIncompressiblePotentialFlowElement #" << this->Id();
    }

    /*PRIVATE*/

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        rSerializer.save("mpPrimalElement", this->mpPrimalElement);
    }

    template <class TPrimalElement>
    void AdjointAnalyticalIncompressiblePotentialFlowElement<TPrimalElement>::load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );
        rSerializer.load("mpPrimalElement", this->mpPrimalElement);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    // Template class instantiation

    template class AdjointAnalyticalIncompressiblePotentialFlowElement<IncompressiblePotentialFlowElement<2,3>>;
    //template class AdjointAnalyticalIncompressiblePotentialFlowElement<CompressiblePotentialFlowElement<2,3>>;
} // namespace Kratos.

