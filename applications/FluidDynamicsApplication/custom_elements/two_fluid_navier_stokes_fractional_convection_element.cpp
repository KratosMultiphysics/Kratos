//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Uxue Chasco
//

#include "two_fluid_navier_stokes_fractional_convection_element.h"

namespace Kratos
{

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    TwoFluidNavierStokesFractionalConvectionElement<TElementData>::TwoFluidNavierStokesFractionalConvectionElement(IndexType NewId) : Element(NewId)
    {
    }

    template <class TElementData>
    TwoFluidNavierStokesFractionalConvectionElement<TElementData>::TwoFluidNavierStokesFractionalConvectionElement(IndexType NewId, const NodesArrayType &ThisNodes) : Element(NewId, ThisNodes)
    {
    }

    template <class TElementData>
    TwoFluidNavierStokesFractionalConvectionElement<TElementData>::TwoFluidNavierStokesFractionalConvectionElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry)
    {
    }

    template <class TElementData>
    TwoFluidNavierStokesFractionalConvectionElement<TElementData>::TwoFluidNavierStokesFractionalConvectionElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Element(NewId, pGeometry, pProperties)
    {
    }

    template <class TElementData>
    TwoFluidNavierStokesFractionalConvectionElement<TElementData>::~TwoFluidNavierStokesFractionalConvectionElement()
    {
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    Element::Pointer TwoFluidNavierStokesFractionalConvectionElement<TElementData>::Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        Properties::Pointer pProperties) const
    {
        return Kratos::make_intrusive<TwoFluidNavierStokesFractionalConvectionElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    /***********************************************************************************/

    template <class TElementData>
    Element::Pointer TwoFluidNavierStokesFractionalConvectionElement<TElementData>::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        Properties::Pointer pProperties) const
    {
        return Kratos::make_intrusive<TwoFluidNavierStokesFractionalConvectionElement>(NewId, pGeom, pProperties);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::CalculateRightHandSide(VectorType &rRightHandSideVector, const ProcessInfo &rCurrentProcessInfo)
    {
        KRATOS_THROW_ERROR(std::runtime_error, "CalculateRightHandSide not implemented", "");
    }

    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::CalculateLocalSystem(
        MatrixType &rLeftHandSideMatrix,
        VectorType &rRightHandSideVector,
        const ProcessInfo &rCurrentProcessInfo)
    {

        if (rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();
        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; ++g)
        {
            UpdateIntegrationPointData(data, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);
            this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::UpdateIntegrationPointData(
        TElementData &rData,
        unsigned int IntegrationPointIndex,
        double Weight,
        const typename TElementData::MatrixRowType &rN,
        const typename TElementData::ShapeDerivativesType &rDN_DX) const
    {
        rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::AddTimeIntegratedSystem(
        TElementData &rData,
        MatrixType &rLHS,
        VectorType &rRHS)
    {
        this->ComputeGaussPointLHSContribution(rData, rLHS);
        this->ComputeGaussPointRHSContribution(rData, rRHS);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                                                                   Matrix &rNContainer,
                                                                                   ShapeFunctionDerivativesArrayType &rDN_DX) const
    {
        const GeometryData::IntegrationMethod integration_method = this->GetIntegrationMethod();
        const GeometryType &r_geometry = this->GetGeometry();
        const unsigned int number_of_gauss_points = r_geometry.IntegrationPointsNumber(integration_method);

        Vector DetJ;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, integration_method);

        if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != NumNodes)
        {
            rNContainer.resize(number_of_gauss_points, NumNodes, false);
        }
        rNContainer = r_geometry.ShapeFunctionsValues(integration_method);

        const GeometryType::IntegrationPointsArrayType &IntegrationPoints = r_geometry.IntegrationPoints(integration_method);

        if (rGaussWeights.size() != number_of_gauss_points)
        {
            rGaussWeights.resize(number_of_gauss_points, false);
        }

        for (unsigned int g = 0; g < number_of_gauss_points; g++)
            rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<2, 3>>::ComputeGaussPointLHSContribution(
        TwoFluidNavierStokesFractionalConvectionElementData<2, 3> &rData,
        MatrixType &rLHS)
    {

        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double dyn_tau = rData.DynamicTau;
        const auto &vfrac = rData.FractionalVelocity;
        const auto &vmesh = rData.MeshVelocity;
        const auto &vconv = vfrac-vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        // Add LHS Gauss point contribution
        const double w_gauss =rData.Weight;

        const double crLHS0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS2 = DN(0,0)*crLHS0 + DN(0,1)*crLHS1;
const double crLHS3 = N[0]*bdf0;
const double crLHS4 = crLHS2 + crLHS3;
const double crLHS5 = 1.0*1.0/(dyn_tau*1.0/dt + stab_c2*1.0/h*sqrt(crLHS0*crLHS0 + crLHS1*crLHS1));
const double crLHS6 = crLHS5*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crLHS7 = N[0]*crLHS6;
const double crLHS8 = crLHS2*crLHS5;
const double crLHS9 = w_gauss*(N[0]*crLHS2 + bdf0*(N[0]*N[0]) + crLHS4*crLHS7 + crLHS4*crLHS8);
const double crLHS10 = N[1]*crLHS3;
const double crLHS11 = DN(1,0)*crLHS0 + DN(1,1)*crLHS1;
const double crLHS12 = N[1]*bdf0;
const double crLHS13 = crLHS11 + crLHS12;
const double crLHS14 = w_gauss*(N[0]*crLHS11 + crLHS10 + crLHS13*crLHS7 + crLHS13*crLHS8);
const double crLHS15 = N[2]*crLHS3;
const double crLHS16 = DN(2,0)*crLHS0 + DN(2,1)*crLHS1;
const double crLHS17 = N[2]*bdf0 + crLHS16;
const double crLHS18 = w_gauss*(N[0]*crLHS16 + crLHS15 + crLHS17*crLHS7 + crLHS17*crLHS8);
const double crLHS19 = N[1]*crLHS6;
const double crLHS20 = crLHS11*crLHS5;
const double crLHS21 = w_gauss*(N[1]*crLHS2 + crLHS10 + crLHS19*crLHS4 + crLHS20*crLHS4);
const double crLHS22 = w_gauss*(N[1]*crLHS11 + bdf0*(N[1]*N[1]) + crLHS13*crLHS19 + crLHS13*crLHS20);
const double crLHS23 = N[2]*crLHS12;
const double crLHS24 = w_gauss*(N[1]*crLHS16 + crLHS17*crLHS19 + crLHS17*crLHS20 + crLHS23);
const double crLHS25 = N[2]*crLHS6;
const double crLHS26 = crLHS16*crLHS5;
const double crLHS27 = w_gauss*(N[2]*crLHS2 + crLHS15 + crLHS25*crLHS4 + crLHS26*crLHS4);
const double crLHS28 = w_gauss*(N[2]*crLHS11 + crLHS13*crLHS25 + crLHS13*crLHS26 + crLHS23);
const double crLHS29 = w_gauss*(N[2]*crLHS16 + bdf0*(N[2]*N[2]) + crLHS17*crLHS25 + crLHS17*crLHS26);
rLHS(0,0)+=crLHS9;
rLHS(0,1)+=0;
rLHS(0,2)+=crLHS14;
rLHS(0,3)+=0;
rLHS(0,4)+=crLHS18;
rLHS(0,5)+=0;
rLHS(1,0)+=0;
rLHS(1,1)+=crLHS9;
rLHS(1,2)+=0;
rLHS(1,3)+=crLHS14;
rLHS(1,4)+=0;
rLHS(1,5)+=crLHS18;
rLHS(2,0)+=crLHS21;
rLHS(2,1)+=0;
rLHS(2,2)+=crLHS22;
rLHS(2,3)+=0;
rLHS(2,4)+=crLHS24;
rLHS(2,5)+=0;
rLHS(3,0)+=0;
rLHS(3,1)+=crLHS21;
rLHS(3,2)+=0;
rLHS(3,3)+=crLHS22;
rLHS(3,4)+=0;
rLHS(3,5)+=crLHS24;
rLHS(4,0)+=crLHS27;
rLHS(4,1)+=0;
rLHS(4,2)+=crLHS28;
rLHS(4,3)+=0;
rLHS(4,4)+=crLHS29;
rLHS(4,5)+=0;
rLHS(5,0)+=0;
rLHS(5,1)+=crLHS27;
rLHS(5,2)+=0;
rLHS(5,3)+=crLHS28;
rLHS(5,4)+=0;
rLHS(5,5)+=crLHS29;


    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<3, 4>>::ComputeGaussPointLHSContribution(
        TwoFluidNavierStokesFractionalConvectionElementData<3, 4> &rData,
        MatrixType &rLHS)
    {


        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;

        const double dyn_tau = rData.DynamicTau;
        const auto &vfrac = rData.FractionalVelocity;
        const auto &vmesh = rData.MeshVelocity;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Stabilization parameters
        constexpr double stab_c2 = 2.0;

        // Add LHS Gauss point contribution
        const double w_gauss =rData.Weight;

        const double crLHS0 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS1 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS2 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS3 = DN(0,0)*crLHS0 + DN(0,1)*crLHS1 + DN(0,2)*crLHS2;
const double crLHS4 = N[0]*bdf0;
const double crLHS5 = crLHS3 + crLHS4;
const double crLHS6 = 1.0*1.0/(dyn_tau*1.0/dt + stab_c2*1.0/h*sqrt(crLHS0*crLHS0 + crLHS1*crLHS1 + crLHS2*crLHS2));
const double crLHS7 = crLHS6*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crLHS8 = N[0]*crLHS7;
const double crLHS9 = crLHS3*crLHS6;
const double crLHS10 = w_gauss*(N[0]*crLHS3 + bdf0*(N[0]*N[0]) + crLHS5*crLHS8 + crLHS5*crLHS9);
const double crLHS11 = N[1]*crLHS4;
const double crLHS12 = DN(1,0)*crLHS0 + DN(1,1)*crLHS1 + DN(1,2)*crLHS2;
const double crLHS13 = N[1]*bdf0;
const double crLHS14 = crLHS12 + crLHS13;
const double crLHS15 = w_gauss*(N[0]*crLHS12 + crLHS11 + crLHS14*crLHS8 + crLHS14*crLHS9);
const double crLHS16 = N[2]*crLHS4;
const double crLHS17 = DN(2,0)*crLHS0 + DN(2,1)*crLHS1 + DN(2,2)*crLHS2;
const double crLHS18 = N[2]*bdf0;
const double crLHS19 = crLHS17 + crLHS18;
const double crLHS20 = w_gauss*(N[0]*crLHS17 + crLHS16 + crLHS19*crLHS8 + crLHS19*crLHS9);
const double crLHS21 = N[3]*crLHS4;
const double crLHS22 = DN(3,0)*crLHS0 + DN(3,1)*crLHS1 + DN(3,2)*crLHS2;
const double crLHS23 = N[3]*bdf0 + crLHS22;
const double crLHS24 = w_gauss*(N[0]*crLHS22 + crLHS21 + crLHS23*crLHS8 + crLHS23*crLHS9);
const double crLHS25 = N[1]*crLHS7;
const double crLHS26 = crLHS12*crLHS6;
const double crLHS27 = w_gauss*(N[1]*crLHS3 + crLHS11 + crLHS25*crLHS5 + crLHS26*crLHS5);
const double crLHS28 = w_gauss*(N[1]*crLHS12 + bdf0*(N[1]*N[1]) + crLHS14*crLHS25 + crLHS14*crLHS26);
const double crLHS29 = N[2]*crLHS13;
const double crLHS30 = w_gauss*(N[1]*crLHS17 + crLHS19*crLHS25 + crLHS19*crLHS26 + crLHS29);
const double crLHS31 = N[3]*crLHS13;
const double crLHS32 = w_gauss*(N[1]*crLHS22 + crLHS23*crLHS25 + crLHS23*crLHS26 + crLHS31);
const double crLHS33 = N[2]*crLHS7;
const double crLHS34 = crLHS17*crLHS6;
const double crLHS35 = w_gauss*(N[2]*crLHS3 + crLHS16 + crLHS33*crLHS5 + crLHS34*crLHS5);
const double crLHS36 = w_gauss*(N[2]*crLHS12 + crLHS14*crLHS33 + crLHS14*crLHS34 + crLHS29);
const double crLHS37 = w_gauss*(N[2]*crLHS17 + bdf0*(N[2]*N[2]) + crLHS19*crLHS33 + crLHS19*crLHS34);
const double crLHS38 = N[3]*crLHS18;
const double crLHS39 = w_gauss*(N[2]*crLHS22 + crLHS23*crLHS33 + crLHS23*crLHS34 + crLHS38);
const double crLHS40 = N[3]*crLHS7;
const double crLHS41 = crLHS22*crLHS6;
const double crLHS42 = w_gauss*(N[3]*crLHS3 + crLHS21 + crLHS40*crLHS5 + crLHS41*crLHS5);
const double crLHS43 = w_gauss*(N[3]*crLHS12 + crLHS14*crLHS40 + crLHS14*crLHS41 + crLHS31);
const double crLHS44 = w_gauss*(N[3]*crLHS17 + crLHS19*crLHS40 + crLHS19*crLHS41 + crLHS38);
const double crLHS45 = w_gauss*(N[3]*crLHS22 + bdf0*(N[3]*N[3]) + crLHS23*crLHS40 + crLHS23*crLHS41);
rLHS(0,0)+=crLHS10;
rLHS(0,1)+=0;
rLHS(0,2)+=0;
rLHS(0,3)+=crLHS15;
rLHS(0,4)+=0;
rLHS(0,5)+=0;
rLHS(0,6)+=crLHS20;
rLHS(0,7)+=0;
rLHS(0,8)+=0;
rLHS(0,9)+=crLHS24;
rLHS(0,10)+=0;
rLHS(0,11)+=0;
rLHS(1,0)+=0;
rLHS(1,1)+=crLHS10;
rLHS(1,2)+=0;
rLHS(1,3)+=0;
rLHS(1,4)+=crLHS15;
rLHS(1,5)+=0;
rLHS(1,6)+=0;
rLHS(1,7)+=crLHS20;
rLHS(1,8)+=0;
rLHS(1,9)+=0;
rLHS(1,10)+=crLHS24;
rLHS(1,11)+=0;
rLHS(2,0)+=0;
rLHS(2,1)+=0;
rLHS(2,2)+=crLHS10;
rLHS(2,3)+=0;
rLHS(2,4)+=0;
rLHS(2,5)+=crLHS15;
rLHS(2,6)+=0;
rLHS(2,7)+=0;
rLHS(2,8)+=crLHS20;
rLHS(2,9)+=0;
rLHS(2,10)+=0;
rLHS(2,11)+=crLHS24;
rLHS(3,0)+=crLHS27;
rLHS(3,1)+=0;
rLHS(3,2)+=0;
rLHS(3,3)+=crLHS28;
rLHS(3,4)+=0;
rLHS(3,5)+=0;
rLHS(3,6)+=crLHS30;
rLHS(3,7)+=0;
rLHS(3,8)+=0;
rLHS(3,9)+=crLHS32;
rLHS(3,10)+=0;
rLHS(3,11)+=0;
rLHS(4,0)+=0;
rLHS(4,1)+=crLHS27;
rLHS(4,2)+=0;
rLHS(4,3)+=0;
rLHS(4,4)+=crLHS28;
rLHS(4,5)+=0;
rLHS(4,6)+=0;
rLHS(4,7)+=crLHS30;
rLHS(4,8)+=0;
rLHS(4,9)+=0;
rLHS(4,10)+=crLHS32;
rLHS(4,11)+=0;
rLHS(5,0)+=0;
rLHS(5,1)+=0;
rLHS(5,2)+=crLHS27;
rLHS(5,3)+=0;
rLHS(5,4)+=0;
rLHS(5,5)+=crLHS28;
rLHS(5,6)+=0;
rLHS(5,7)+=0;
rLHS(5,8)+=crLHS30;
rLHS(5,9)+=0;
rLHS(5,10)+=0;
rLHS(5,11)+=crLHS32;
rLHS(6,0)+=crLHS35;
rLHS(6,1)+=0;
rLHS(6,2)+=0;
rLHS(6,3)+=crLHS36;
rLHS(6,4)+=0;
rLHS(6,5)+=0;
rLHS(6,6)+=crLHS37;
rLHS(6,7)+=0;
rLHS(6,8)+=0;
rLHS(6,9)+=crLHS39;
rLHS(6,10)+=0;
rLHS(6,11)+=0;
rLHS(7,0)+=0;
rLHS(7,1)+=crLHS35;
rLHS(7,2)+=0;
rLHS(7,3)+=0;
rLHS(7,4)+=crLHS36;
rLHS(7,5)+=0;
rLHS(7,6)+=0;
rLHS(7,7)+=crLHS37;
rLHS(7,8)+=0;
rLHS(7,9)+=0;
rLHS(7,10)+=crLHS39;
rLHS(7,11)+=0;
rLHS(8,0)+=0;
rLHS(8,1)+=0;
rLHS(8,2)+=crLHS35;
rLHS(8,3)+=0;
rLHS(8,4)+=0;
rLHS(8,5)+=crLHS36;
rLHS(8,6)+=0;
rLHS(8,7)+=0;
rLHS(8,8)+=crLHS37;
rLHS(8,9)+=0;
rLHS(8,10)+=0;
rLHS(8,11)+=crLHS39;
rLHS(9,0)+=crLHS42;
rLHS(9,1)+=0;
rLHS(9,2)+=0;
rLHS(9,3)+=crLHS43;
rLHS(9,4)+=0;
rLHS(9,5)+=0;
rLHS(9,6)+=crLHS44;
rLHS(9,7)+=0;
rLHS(9,8)+=0;
rLHS(9,9)+=crLHS45;
rLHS(9,10)+=0;
rLHS(9,11)+=0;
rLHS(10,0)+=0;
rLHS(10,1)+=crLHS42;
rLHS(10,2)+=0;
rLHS(10,3)+=0;
rLHS(10,4)+=crLHS43;
rLHS(10,5)+=0;
rLHS(10,6)+=0;
rLHS(10,7)+=crLHS44;
rLHS(10,8)+=0;
rLHS(10,9)+=0;
rLHS(10,10)+=crLHS45;
rLHS(10,11)+=0;
rLHS(11,0)+=0;
rLHS(11,1)+=0;
rLHS(11,2)+=crLHS42;
rLHS(11,3)+=0;
rLHS(11,4)+=0;
rLHS(11,5)+=crLHS43;
rLHS(11,6)+=0;
rLHS(11,7)+=0;
rLHS(11,8)+=crLHS44;
rLHS(11,9)+=0;
rLHS(11,10)+=0;
rLHS(11,11)+=crLHS45;


    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <>
    void TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<2, 3>>::ComputeGaussPointRHSContribution(
        TwoFluidNavierStokesFractionalConvectionElementData<2, 3> &rData,
        VectorType &rRHS)
    {

        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;
        const double dyn_tau = rData.DynamicTau;
        const auto &vn = rData.VelocityOldStep1;
        const auto &vnn = rData.VelocityOldStep2;
        // const auto &vnnn = rData.Velocity_OldStep3; #an bdf2
        const auto &vmesh = rData.MeshVelocity;
        const auto &vfrac = rData.FractionalVelocity;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Add RHS Gauss point contribution
        const double w_gauss = rData.Weight;
        
        const double crRHS0 = N[0]*(vn(0,0) - vnn(0,0));
const double crRHS1 = N[1]*(vn(1,0) - vnn(1,0));
const double crRHS2 = N[2]*(vn(2,0) - vnn(2,0));
const double crRHS3 = crRHS0 + crRHS1 + crRHS2;
const double crRHS4 = 1.0/dt;
const double crRHS5 = N[0]*crRHS4;
const double crRHS6 = N[0]*(bdf0*vfrac(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*vfrac(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*vfrac(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
const double crRHS7 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS9 = crRHS7*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0)) + crRHS8*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0));
const double crRHS10 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0);
const double crRHS11 = crRHS10*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0));
const double crRHS12 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1);
const double crRHS13 = crRHS12*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0));
const double crRHS14 = crRHS11 + crRHS13;
const double crRHS15 = -crRHS0*crRHS4 - crRHS1*crRHS4 - crRHS11 - crRHS13 - crRHS2*crRHS4 + crRHS6 + crRHS9;
const double crRHS16 = 1.0*1.0/(crRHS4*dyn_tau + stab_c2*1.0/h*sqrt(crRHS7*crRHS7 + crRHS8*crRHS8));
const double crRHS17 = crRHS16*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crRHS18 = N[0]*crRHS17;
const double crRHS19 = crRHS16*(DN(0,0)*crRHS7 + DN(0,1)*crRHS8);
const double crRHS20 = N[0]*(vn(0,1) - vnn(0,1));
const double crRHS21 = N[1]*(vn(1,1) - vnn(1,1));
const double crRHS22 = N[2]*(vn(2,1) - vnn(2,1));
const double crRHS23 = crRHS20 + crRHS21 + crRHS22;
const double crRHS24 = N[0]*(bdf0*vfrac(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*vfrac(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*vfrac(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
const double crRHS25 = crRHS7*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1)) + crRHS8*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1));
const double crRHS26 = crRHS10*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1));
const double crRHS27 = crRHS12*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1));
const double crRHS28 = crRHS26 + crRHS27;
const double crRHS29 = -crRHS20*crRHS4 - crRHS21*crRHS4 - crRHS22*crRHS4 + crRHS24 + crRHS25 - crRHS26 - crRHS27;
const double crRHS30 = N[1]*crRHS4;
const double crRHS31 = N[1]*crRHS17;
const double crRHS32 = crRHS16*(DN(1,0)*crRHS7 + DN(1,1)*crRHS8);
const double crRHS33 = N[2]*crRHS4;
const double crRHS34 = N[2]*crRHS17;
const double crRHS35 = crRHS16*(DN(2,0)*crRHS7 + DN(2,1)*crRHS8);
rRHS[0]+=-w_gauss*(-N[0]*crRHS14 + N[0]*crRHS6 + N[0]*crRHS9 + crRHS15*crRHS18 + crRHS15*crRHS19 - crRHS3*crRHS5);
rRHS[1]+=-w_gauss*(N[0]*crRHS24 + N[0]*crRHS25 - N[0]*crRHS28 + crRHS18*crRHS29 + crRHS19*crRHS29 - crRHS23*crRHS5);
rRHS[2]+=-w_gauss*(-N[1]*crRHS14 + N[1]*crRHS6 + N[1]*crRHS9 + crRHS15*crRHS31 + crRHS15*crRHS32 - crRHS3*crRHS30);
rRHS[3]+=-w_gauss*(N[1]*crRHS24 + N[1]*crRHS25 - N[1]*crRHS28 - crRHS23*crRHS30 + crRHS29*crRHS31 + crRHS29*crRHS32);
rRHS[4]+=-w_gauss*(-N[2]*crRHS14 + N[2]*crRHS6 + N[2]*crRHS9 + crRHS15*crRHS34 + crRHS15*crRHS35 - crRHS3*crRHS33);
rRHS[5]+=-w_gauss*(N[2]*crRHS24 + N[2]*crRHS25 - N[2]*crRHS28 - crRHS23*crRHS33 + crRHS29*crRHS34 + crRHS29*crRHS35);


    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <>
    void TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<3, 4>>::ComputeGaussPointRHSContribution(
        TwoFluidNavierStokesFractionalConvectionElementData<3, 4> &rData,
        VectorType &rRHS)
    {


        const double h = rData.ElementSize;
        const double dt = rData.DeltaTime;
        const double bdf0 = rData.bdf0;
        const double bdf1 = rData.bdf1;
        const double bdf2 = rData.bdf2;
        const double dyn_tau = rData.DynamicTau;
        const auto &vn = rData.VelocityOldStep1;
        const auto &vnn = rData.VelocityOldStep2;
        const auto &vmesh = rData.MeshVelocity;
        const auto &vfrac = rData.FractionalVelocity;
        const auto &vconv = vfrac - vmesh;

        // Get shape function values
        const auto &N = rData.N;
        const auto &DN = rData.DN_DX;

        // Add RHS Gauss point contribution
        const double w_gauss =rData.Weight;

        const double crRHS0 = N[0]*(vn(0,0) - vnn(0,0));
const double crRHS1 = N[1]*(vn(1,0) - vnn(1,0));
const double crRHS2 = N[2]*(vn(2,0) - vnn(2,0));
const double crRHS3 = N[3]*(vn(3,0) - vnn(3,0));
const double crRHS4 = crRHS0 + crRHS1 + crRHS2 + crRHS3;
const double crRHS5 = 1.0/dt;
const double crRHS6 = N[0]*crRHS5;
const double crRHS7 = N[0]*(bdf0*vfrac(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*vfrac(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*vfrac(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*vfrac(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
const double crRHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS11 = crRHS10*(DN(0,2)*vfrac(0,0) + DN(1,2)*vfrac(1,0) + DN(2,2)*vfrac(2,0) + DN(3,2)*vfrac(3,0)) + crRHS8*(DN(0,0)*vfrac(0,0) + DN(1,0)*vfrac(1,0) + DN(2,0)*vfrac(2,0) + DN(3,0)*vfrac(3,0)) + crRHS9*(DN(0,1)*vfrac(0,0) + DN(1,1)*vfrac(1,0) + DN(2,1)*vfrac(2,0) + DN(3,1)*vfrac(3,0));
const double crRHS12 = N[0]*vn(0,0) + N[1]*vn(1,0) + N[2]*vn(2,0) + N[3]*vn(3,0);
const double crRHS13 = crRHS12*(DN(0,0)*vn(0,0) + DN(1,0)*vn(1,0) + DN(2,0)*vn(2,0) + DN(3,0)*vn(3,0));
const double crRHS14 = N[0]*vn(0,1) + N[1]*vn(1,1) + N[2]*vn(2,1) + N[3]*vn(3,1);
const double crRHS15 = crRHS14*(DN(0,1)*vn(0,0) + DN(1,1)*vn(1,0) + DN(2,1)*vn(2,0) + DN(3,1)*vn(3,0));
const double crRHS16 = N[0]*vn(0,2) + N[1]*vn(1,2) + N[2]*vn(2,2) + N[3]*vn(3,2);
const double crRHS17 = crRHS16*(DN(0,2)*vn(0,0) + DN(1,2)*vn(1,0) + DN(2,2)*vn(2,0) + DN(3,2)*vn(3,0));
const double crRHS18 = crRHS13 + crRHS15 + crRHS17;
const double crRHS19 = -crRHS0*crRHS5 - crRHS1*crRHS5 + crRHS11 - crRHS13 - crRHS15 - crRHS17 - crRHS2*crRHS5 - crRHS3*crRHS5 + crRHS7;
const double crRHS20 = 1.0*1.0/(crRHS5*dyn_tau + stab_c2*1.0/h*sqrt(crRHS10*crRHS10 + crRHS8*crRHS8 + crRHS9*crRHS9));
const double crRHS21 = crRHS20*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crRHS22 = N[0]*crRHS21;
const double crRHS23 = crRHS20*(DN(0,0)*crRHS8 + DN(0,1)*crRHS9 + DN(0,2)*crRHS10);
const double crRHS24 = N[0]*(vn(0,1) - vnn(0,1));
const double crRHS25 = N[1]*(vn(1,1) - vnn(1,1));
const double crRHS26 = N[2]*(vn(2,1) - vnn(2,1));
const double crRHS27 = N[3]*(vn(3,1) - vnn(3,1));
const double crRHS28 = crRHS24 + crRHS25 + crRHS26 + crRHS27;
const double crRHS29 = N[0]*(bdf0*vfrac(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*vfrac(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*vfrac(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*vfrac(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
const double crRHS30 = crRHS10*(DN(0,2)*vfrac(0,1) + DN(1,2)*vfrac(1,1) + DN(2,2)*vfrac(2,1) + DN(3,2)*vfrac(3,1)) + crRHS8*(DN(0,0)*vfrac(0,1) + DN(1,0)*vfrac(1,1) + DN(2,0)*vfrac(2,1) + DN(3,0)*vfrac(3,1)) + crRHS9*(DN(0,1)*vfrac(0,1) + DN(1,1)*vfrac(1,1) + DN(2,1)*vfrac(2,1) + DN(3,1)*vfrac(3,1));
const double crRHS31 = crRHS12*(DN(0,0)*vn(0,1) + DN(1,0)*vn(1,1) + DN(2,0)*vn(2,1) + DN(3,0)*vn(3,1));
const double crRHS32 = crRHS14*(DN(0,1)*vn(0,1) + DN(1,1)*vn(1,1) + DN(2,1)*vn(2,1) + DN(3,1)*vn(3,1));
const double crRHS33 = crRHS16*(DN(0,2)*vn(0,1) + DN(1,2)*vn(1,1) + DN(2,2)*vn(2,1) + DN(3,2)*vn(3,1));
const double crRHS34 = crRHS31 + crRHS32 + crRHS33;
const double crRHS35 = -crRHS24*crRHS5 - crRHS25*crRHS5 - crRHS26*crRHS5 - crRHS27*crRHS5 + crRHS29 + crRHS30 - crRHS31 - crRHS32 - crRHS33;
const double crRHS36 = N[0]*(vn(0,2) - vnn(0,2));
const double crRHS37 = N[1]*(vn(1,2) - vnn(1,2));
const double crRHS38 = N[2]*(vn(2,2) - vnn(2,2));
const double crRHS39 = N[3]*(vn(3,2) - vnn(3,2));
const double crRHS40 = crRHS36 + crRHS37 + crRHS38 + crRHS39;
const double crRHS41 = N[0]*(bdf0*vfrac(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*vfrac(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*vfrac(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*vfrac(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
const double crRHS42 = crRHS10*(DN(0,2)*vfrac(0,2) + DN(1,2)*vfrac(1,2) + DN(2,2)*vfrac(2,2) + DN(3,2)*vfrac(3,2)) + crRHS8*(DN(0,0)*vfrac(0,2) + DN(1,0)*vfrac(1,2) + DN(2,0)*vfrac(2,2) + DN(3,0)*vfrac(3,2)) + crRHS9*(DN(0,1)*vfrac(0,2) + DN(1,1)*vfrac(1,2) + DN(2,1)*vfrac(2,2) + DN(3,1)*vfrac(3,2));
const double crRHS43 = crRHS12*(DN(0,0)*vn(0,2) + DN(1,0)*vn(1,2) + DN(2,0)*vn(2,2) + DN(3,0)*vn(3,2));
const double crRHS44 = crRHS14*(DN(0,1)*vn(0,2) + DN(1,1)*vn(1,2) + DN(2,1)*vn(2,2) + DN(3,1)*vn(3,2));
const double crRHS45 = crRHS16*(DN(0,2)*vn(0,2) + DN(1,2)*vn(1,2) + DN(2,2)*vn(2,2) + DN(3,2)*vn(3,2));
const double crRHS46 = crRHS43 + crRHS44 + crRHS45;
const double crRHS47 = -crRHS36*crRHS5 - crRHS37*crRHS5 - crRHS38*crRHS5 - crRHS39*crRHS5 + crRHS41 + crRHS42 - crRHS43 - crRHS44 - crRHS45;
const double crRHS48 = N[1]*crRHS5;
const double crRHS49 = N[1]*crRHS21;
const double crRHS50 = crRHS20*(DN(1,0)*crRHS8 + DN(1,1)*crRHS9 + DN(1,2)*crRHS10);
const double crRHS51 = N[2]*crRHS5;
const double crRHS52 = N[2]*crRHS21;
const double crRHS53 = crRHS20*(DN(2,0)*crRHS8 + DN(2,1)*crRHS9 + DN(2,2)*crRHS10);
const double crRHS54 = N[3]*crRHS5;
const double crRHS55 = N[3]*crRHS21;
const double crRHS56 = crRHS20*(DN(3,0)*crRHS8 + DN(3,1)*crRHS9 + DN(3,2)*crRHS10);
rRHS[0]+=-w_gauss*(N[0]*crRHS11 - N[0]*crRHS18 + N[0]*crRHS7 + crRHS19*crRHS22 + crRHS19*crRHS23 - crRHS4*crRHS6);
rRHS[1]+=-w_gauss*(N[0]*crRHS29 + N[0]*crRHS30 - N[0]*crRHS34 + crRHS22*crRHS35 + crRHS23*crRHS35 - crRHS28*crRHS6);
rRHS[2]+=-w_gauss*(N[0]*crRHS41 + N[0]*crRHS42 - N[0]*crRHS46 + crRHS22*crRHS47 + crRHS23*crRHS47 - crRHS40*crRHS6);
rRHS[3]+=-w_gauss*(N[1]*crRHS11 - N[1]*crRHS18 + N[1]*crRHS7 + crRHS19*crRHS49 + crRHS19*crRHS50 - crRHS4*crRHS48);
rRHS[4]+=-w_gauss*(N[1]*crRHS29 + N[1]*crRHS30 - N[1]*crRHS34 - crRHS28*crRHS48 + crRHS35*crRHS49 + crRHS35*crRHS50);
rRHS[5]+=-w_gauss*(N[1]*crRHS41 + N[1]*crRHS42 - N[1]*crRHS46 - crRHS40*crRHS48 + crRHS47*crRHS49 + crRHS47*crRHS50);
rRHS[6]+=-w_gauss*(N[2]*crRHS11 - N[2]*crRHS18 + N[2]*crRHS7 + crRHS19*crRHS52 + crRHS19*crRHS53 - crRHS4*crRHS51);
rRHS[7]+=-w_gauss*(N[2]*crRHS29 + N[2]*crRHS30 - N[2]*crRHS34 - crRHS28*crRHS51 + crRHS35*crRHS52 + crRHS35*crRHS53);
rRHS[8]+=-w_gauss*(N[2]*crRHS41 + N[2]*crRHS42 - N[2]*crRHS46 - crRHS40*crRHS51 + crRHS47*crRHS52 + crRHS47*crRHS53);
rRHS[9]+=-w_gauss*(N[3]*crRHS11 - N[3]*crRHS18 + N[3]*crRHS7 + crRHS19*crRHS55 + crRHS19*crRHS56 - crRHS4*crRHS54);
rRHS[10]+=-w_gauss*(N[3]*crRHS29 + N[3]*crRHS30 - N[3]*crRHS34 - crRHS28*crRHS54 + crRHS35*crRHS55 + crRHS35*crRHS56);
rRHS[11]+=-w_gauss*(N[3]*crRHS41 + N[3]*crRHS42 - N[3]*crRHS46 - crRHS40*crRHS54 + crRHS47*crRHS55 + crRHS47*crRHS56);


    }

    /***********************************************************************************/
    /***********************************************************************************/
    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::EquationIdVector(
        EquationIdVectorType &rResult,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const GeometryType &r_geometry = this->GetGeometry();

        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACTIONAL_VELOCITY_X);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_X, xpos).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_Y, xpos + 1).EquationId();
            if constexpr (Dim == 3)
                rResult[LocalIndex++] = r_geometry[i].GetDof(FRACTIONAL_VELOCITY_Z, xpos + 2).EquationId();
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    template <class TElementData>
    void TwoFluidNavierStokesFractionalConvectionElement<TElementData>::GetDofList(
        DofsVectorType &rElementalDofList,
        const ProcessInfo &rCurrentProcessInfo) const
    {
        const GeometryType &r_geometry = this->GetGeometry();

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACTIONAL_VELOCITY_X);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_X, xpos);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_Y, xpos + 1);
            if constexpr (Dim == 3)
                rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(FRACTIONAL_VELOCITY_Z, xpos + 2);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /***********************************************************************************/
    /***********************************************************************************/
    template class TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<2, 3>>;
    template class TwoFluidNavierStokesFractionalConvectionElement<TwoFluidNavierStokesFractionalConvectionElementData<3, 4>>;

    /***********************************************************************************/
    /***********************************************************************************/
}