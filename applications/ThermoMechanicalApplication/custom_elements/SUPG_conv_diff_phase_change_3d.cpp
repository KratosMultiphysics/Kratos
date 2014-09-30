/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/SUPG_conv_diff_phase_change_3d.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

SUPGConvDiffPhaseChange3D::SUPGConvDiffPhaseChange3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

SUPGConvDiffPhaseChange3D::SUPGConvDiffPhaseChange3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{

}

Element::Pointer SUPGConvDiffPhaseChange3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new SUPGConvDiffPhaseChange3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

SUPGConvDiffPhaseChange3D::~SUPGConvDiffPhaseChange3D()
{
}

//************************************************************************************
//************************************************************************************

void SUPGConvDiffPhaseChange3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    const unsigned int nodes_number = GetGeometry().size();
    const unsigned int matsize = nodes_number;

    //****************************************************
    //Get Vector of BDF coefficients
    const Vector& BDFVector = rCurrentProcessInfo[BDF_COEFFICIENTS];

    //get variables to be used in the computations
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
    const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
    const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
    const Variable<array_1d<double, 3 > >& rConvVar = my_settings->GetConvectionVariable();


    if (rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!

    //compute geometrical data of the element
    boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
    array_1d<double, 4 > Ncenter;
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, Ncenter, Volume);
    const double gauss_weight = 0.25*Volume;

    //gather variables on all nodes
    array_1d<double,4> temperatures, H, rho, c, k;
    for(unsigned int i=0; i<4; i++)
    {
        temperatures[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
        rho[i] = GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
        k[i] = GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
        H[i] = GetGeometry()[i].FastGetSolutionStepValue(ENTHALPY);;
    }

    const double liquidus_temperature = rCurrentProcessInfo[FLUID_TEMPERATURE];
    const double solidus_temperature = rCurrentProcessInfo[SOLID_TEMPERATURE];
    const double latent_heat = rCurrentProcessInfo[LATENT_HEAT];
//     KRATOS_WATCH(liquidus_temperature);
    ComputeEffectiveSpeficifHeat(c,H,temperatures,liquidus_temperature,solidus_temperature,latent_heat);

    //compute temperature gradient
    array_1d<double,3> gradT = prod(trans(DN_DX),temperatures);
    const double norm_gradT = norm_2(gradT);

    //estimate the element lenght h
    const double h = pow(6.0*Volume,0.3333333333333333333333333333);

    //compute DH_Dt (time variation of enthalpy )
    array_1d<double,4> DH_Dt;
    for(unsigned int i=0; i<4; i++) DH_Dt[i] = BDFVector[0]*H[i];

    for(unsigned int step=1; step<BDFVector.size(); step++)
        for(unsigned int i=0; i<4; i++)
            DH_Dt[i] += BDFVector[step]*GetGeometry()[i].FastGetSolutionStepValue(ENTHALPY,step);

    //add the conductivity. does not make sense to add it in the loop on gauss points since
    //the laplacian is all constant through the element
    double kavg = 0.0;
    for(unsigned int i=0; i<4; i++)
        kavg += 0.25*k[i];
    noalias(rLeftHandSideMatrix) = Volume*kavg*prod(DN_DX,trans(DN_DX));
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, temperatures);

    //we will do nodal integration and will consider the problem fully non linear thus
    //computing separately the RHS and LHS
    //NOTE THAT NODAL INTEGRATION IS ASSUMED THROUGH THE LOOP!!!!
    for(unsigned int igauss=0; igauss<4; igauss++)
    {
        array_1d<double, 4 > N = ZeroVector(4);
        N[igauss] = 1.0;  //NOTE THAT THE FACT THAT ONLY ONE NODE IS 1 and all the others zero is used in the following!!
        const double Qext = GetGeometry()[igauss].FastGetSolutionStepValue(rSourceVar);
        const double rho_c = rho[igauss]*c[igauss];

        //add external heat source to the rhs
        rRightHandSideVector[igauss] += Qext*N[igauss];

        //get convection velocity
        const array_1d<double,3>& a = GetGeometry()[igauss].FastGetSolutionStepValue(rConvVar)
                                      - GetGeometry()[igauss].FastGetSolutionStepValue(rMeshVelocityVar);

        //**************** LHS Galerkin ******************
        //add time derivatives to the lhs
        rLeftHandSideMatrix(igauss, igauss) += gauss_weight*rho_c*BDFVector[0];

        //add convective contribution to the lhs
        const array_1d<double,4> a_DN = prod(DN_DX,a);
        noalias(rLeftHandSideMatrix) +=  gauss_weight*rho_c * outer_prod(N,a_DN); //this could be optimized taking into account that N is only 1 on the igauss node

        //**************** RHS Galerkin ******************
        //add time derivatives to the rhs --> later on i should substitute this by the enthalpy!!
        rRightHandSideVector[igauss] -= gauss_weight*rho[igauss]*DH_Dt[igauss];

        //add convective contribution to the rhs
        const double a_DH = inner_prod(a_DN, H);
        noalias(rRightHandSideVector)  -= ( gauss_weight*rho[igauss]*a_DH ) * N;  //this could be optimized taking into account that N is only 1 on the igauss node

        //**************** stabilization *****************
        const double tau = ComputeTau(h,k[igauss],a);

        //contribution of stabilization to the lhs
        noalias(rLeftHandSideMatrix) += tau * gauss_weight*rho_c * outer_prod(a_DN,a_DN);
        noalias(rLeftHandSideMatrix) += tau * gauss_weight*rho_c * BDFVector[0] * outer_prod(a_DN,N);

        //contribution to the rhs
        const double residual = Qext - rho[igauss]*(DH_Dt[igauss] + a_DH); //residual = Qext - rho*DH_Dt - rho*a_DH  (other terms disappear for low order elements)
        noalias(rRightHandSideVector) += tau * gauss_weight * residual * a_DN;

        //**************** YZb discontinuity capturing *****************
        const double norm_a2 = inner_prod(a,a);
        const double norm_a = sqrt(norm_a2);        
        if(norm_gradT > 1e-3 && norm_a > 1e-9)
        {
//             const double reference_temperature  = rCurrentProcessInfo[AMBIENT_TEMPERATURE]; //TODO: i would set this to the "solidus" temperature, or to the ambient
//             const double kdc =1.0*ComputeDiscontinuityCapturingDiffusion(DN_DX, gradT, norm_gradT,residual,reference_temperature);
//             
//             //ISOTROPIC VERSION
//             noalias(rLeftHandSideMatrix) += gauss_weight*kdc*prod(DN_DX,trans(DN_DX));
//             noalias(rRightHandSideVector)-= gauss_weight*kdc*prod(DN_DX,gradT);
//             
            //ANISOTROPIC VERSION
//            bounded_matrix<double,3,3> D = 1.0/(norm_gradT*norm_gradT)*outer_prod(gradT,gradT);           
//            bounded_matrix<double,3,4> tmp = prod(D,trans(DN_DX));
//            noalias(rLeftHandSideMatrix) += gauss_weight*kdc*prod(DN_DX,tmp);
//            bounded_matrix<double,4,3> tmp2 = prod(DN_DX,D);
//            noalias(rRightHandSideVector)-= gauss_weight*kdc*prod(tmp2,gradT);
            
            //codina's version of shock capturing
            const double gamma = rho_c*norm_a*h/(2.0*k[igauss]+1e-20);
            const double C = 0.7;
            
            double alpha_dc = 0.0;
            if(gamma > 1e-3)
                alpha_dc = std::max( 0.0 ,  C - 1.0/gamma);
// KRATOS_WATCH(h);
// KRATOS_WATCH(norm_a);
//                 KRATOS_WATCH(gamma);
// KRATOS_WATCH(residual);
// KRATOS_WATCH(alpha_dc);
            const double kdc = /*0.5**/alpha_dc*h*fabs(residual)/fabs(norm_gradT);
            bounded_matrix<double,3,3> D = kdc*( IdentityMatrix(3,3));
            D += (std::max( kdc - tau*norm_a2 , 0.0) - kdc)/(norm_a2) * outer_prod(a,a);
//              if(kdc > 0)   KRATOS_WATCH(kdc);
//                 KRATOS_WATCH(std::max( kdc - tau*norm_a2 , 0.0));
            bounded_matrix<double,4,3> tmp = prod(DN_DX,D);
            noalias(rLeftHandSideMatrix) += gauss_weight*prod(tmp,trans(DN_DX));
            noalias(rRightHandSideVector)-= gauss_weight*prod(tmp,gradT);
            
        }
    }




    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void SUPGConvDiffPhaseChange3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();


    if (rResult.size() != number_of_nodes )
        rResult.resize(number_of_nodes , false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();

    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************

void SUPGConvDiffPhaseChange3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();



    if (ElementalDofList.size() != number_of_nodes )
        ElementalDofList.resize(number_of_nodes );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);
    }
    KRATOS_CATCH("");
}


//***********************************************************************************
//**************************************************************************************
double SUPGConvDiffPhaseChange3D::ComputeTau(const double h, const double k, const array_1d<double,3>& a)
{
    //follow the article DOI: 10.1002/fld.1484 by bazilevs and others
    const double norm_a = norm_2(a);
    const double Pe = norm_a * h / (2.0 * (k + 1e-10) );

    if( Pe < 3.0) return h*h/(12.0*(k + 1e-10));
    else return h/(2.0*norm_a);

}



//***********************************************************************************
//**************************************************************************************
//this function computes a "secant" specific heat in the attempt to improve the convergence of the phase change process
void SUPGConvDiffPhaseChange3D::ComputeEffectiveSpeficifHeat(array_1d<double,4>& c,
        const array_1d<double,4>& H,
        const array_1d<double,4>& temperatures,
        const double fluid_T,
        const double solid_T,
        const double latent_heat
         )
{
    
//      for(unsigned int i=0; i<GetGeometry().size(); i++)
//      {
//         //const double dT = temperatures[i] - GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1);
//         const double dT = GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1) - GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,2);
//         
//         c[i] = GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT,1);
//         if( fabs(dT) > 1e-5)
//         {
//             //const double dH = H[i] - GetGeometry()[i].FastGetSolutionStepValue(ENTHALPY,1);
//             const double dH = GetGeometry()[i].FastGetSolutionStepValue(ENTHALPY,1) - GetGeometry()[i].FastGetSolutionStepValue(ENTHALPY,2);
// //             std::cout << c[i] << " " << std::max(c[i],dH/dT) << std::endl;
//             c[i] = dH/dT;
//         }
//     }
         
         
     for(unsigned int i=0; i<GetGeometry().size(); i++)
     {         
       const double specific_heat = GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
       const double DF_DT =  GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION_RATE);


       if(  temperatures[i] <= fluid_T)
       {
           
           double tangent_DF_DT = DF_DT;

           if( !(solid_T <= temperatures[i] && temperatures[i] <= fluid_T) ) //if not in the range of phase change
           {
               double nd_FF =  GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION);
               double nd_FF_old =  GetGeometry()[i].FastGetSolutionStepValue(SOLID_FRACTION,1);
               double mid_T = 0.5*(solid_T + fluid_T);
               double aux_denom = temperatures[i] - mid_T;
               if(fabs(aux_denom) < 1e-6)
               {
                   if(aux_denom >= 0) aux_denom = 1e-6;
                   else aux_denom = -1e-6;
               }
               tangent_DF_DT = (nd_FF - nd_FF_old)/aux_denom;

           }
           c[i] = specific_heat + fabs(tangent_DF_DT*latent_heat);
       }
       else
       {
           c[i] = specific_heat;
       }
   }

}





//***********************************************************************************
//**************************************************************************************
    const double SUPGConvDiffPhaseChange3D::ComputeDiscontinuityCapturingDiffusion(
        const boost::numeric::ublas::bounded_matrix<double, 4, 3 >& DN_DX,
        const array_1d<double,3>& gradT,
        const double& norm_gradT,
        const double& residual,
        const double& reference_temperature
    )
    {
        //follow the article DOI: 10.1002/fld.1484 by bazilevs and others
        const double& Z = residual;
        const double Yinv = fabs(1.0/reference_temperature);

        //compute hdc_half
        array_1d<double,3> j = gradT / norm_gradT;
        double sum = 0.0;
        for(unsigned int i=0; i<4; i++)
        {
            double tmp = 0.0;
            for(unsigned int k=0; k<3; k++)
                tmp += DN_DX(i,k)*j[k];

            sum += fabs(tmp);
        }
        const double hdc_half = 1.0/sum;

        //compute kdc - case of beta = 1
        //note that we simplified Yinv in the formula!!
        const double kdc = fabs(Z) / (norm_gradT) * hdc_half;
//        const double kdc = fabs(Yinv*Z) / (Yinv * norm_gradT) * hdc_half;
        return kdc;
    }


//***********************************************************************************
//**************************************************************************************

    void SUPGConvDiffPhaseChange3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        MatrixType temp = Matrix();
        CalculateLocalSystem( temp,  rRightHandSideVector,  rCurrentProcessInfo);

        KRATOS_CATCH("")
    }


} // Namespace Kratos


