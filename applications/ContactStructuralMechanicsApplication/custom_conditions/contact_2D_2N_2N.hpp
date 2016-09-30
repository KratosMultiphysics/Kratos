// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_CONTACT2D2N2N)
#define KRATOS_CONTACT2D2N2N

// System includes

// External includes

// Project includes
#include "custom_conditions/mortar_contact_condition.h"
#include "custom_utilities/contact_utilities.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos 
{
    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{
        
    typedef Point<3>                                  PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
        
    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{
        
class Contact2D2N2N
{
public:
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHSInternalContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;

    const Matrix lm            = rContactData.LagrangeMultipliers;
    const Matrix NormalsMaster = rContactData.NormalsMaster;
    const Vector masternormal  = prod(NormalsMaster, N2);
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    
    const double clhs0 =             0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0);
    const double clhs1 =             Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double clhs2 =             0.5*X1(0,0);
    const double clhs3 =             0.5*X1(0,1);
    const double clhs4 =             0.5*X1(1,0);
    const double clhs5 =             0.5*X1(1,1);
    const double clhs6 =             std::sqrt(0.25*std::pow(X1(0,0), 2) + 0.25*std::pow(X1(0,1), 2) + 0.25*std::pow(X1(1,0), 2) - X1(1,0)*clhs2 + 0.25*std::pow(X1(1,1), 2) - X1(1,1)*clhs3 + clhs2*u1(0,0) - clhs2*u1(1,0) + clhs3*u1(0,1) - clhs3*u1(1,1) - clhs4*u1(0,0) + clhs4*u1(1,0) - clhs5*u1(0,1) + clhs5*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double clhs7 =             1.0/clhs6;
    const double clhs8 =             N2[0]*clhs1*clhs7;
    const double clhs9 =             clhs0*clhs8;
    const double clhs10 =             0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1);
    const double clhs11 =             clhs10*clhs8;
    const double clhs12 =             N2[0]*clhs6;
    const double clhs13 =             -Phi[0]*clhs12;
    const double clhs14 =             -Phi[1]*clhs12;
    const double clhs15 =             Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double clhs16 =             N2[0]*clhs15*clhs7;
    const double clhs17 =             clhs0*clhs16;
    const double clhs18 =             clhs10*clhs16;
    const double clhs19 =             N2[1]*clhs1*clhs7;
    const double clhs20 =             clhs0*clhs19;
    const double clhs21 =             clhs10*clhs19;
    const double clhs22 =             N2[1]*clhs6;
    const double clhs23 =             -Phi[0]*clhs22;
    const double clhs24 =             -Phi[1]*clhs22;
    const double clhs25 =             N2[1]*clhs15*clhs7;
    const double clhs26 =             clhs0*clhs25;
    const double clhs27 =             clhs10*clhs25;
    const double clhs28 =             N1[0]*clhs1*clhs7;
    const double clhs29 =             clhs0*clhs28;
    const double clhs30 =             clhs10*clhs28;
    const double clhs31 =             N1[0]*clhs6;
    const double clhs32 =             Phi[0]*clhs31;
    const double clhs33 =             Phi[1]*clhs31;
    const double clhs34 =             N1[0]*clhs15*clhs7;
    const double clhs35 =             clhs0*clhs34;
    const double clhs36 =             clhs10*clhs34;
    const double clhs37 =             N1[1]*clhs1*clhs7;
    const double clhs38 =             clhs0*clhs37;
    const double clhs39 =             clhs10*clhs37;
    const double clhs40 =             N1[1]*clhs6;
    const double clhs41 =             Phi[0]*clhs40;
    const double clhs42 =             Phi[1]*clhs40;
    const double clhs43 =             N1[1]*clhs15*clhs7;
    const double clhs44 =             clhs0*clhs43;
    const double clhs45 =             clhs10*clhs43;
    
    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=-clhs9;
    lhs(0,5)=-clhs11;
    lhs(0,6)=clhs9;
    lhs(0,7)=clhs11;
    lhs(0,8)=clhs13;
    lhs(0,9)=0;
    lhs(0,10)=clhs14;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=-clhs17;
    lhs(1,5)=-clhs18;
    lhs(1,6)=clhs17;
    lhs(1,7)=clhs18;
    lhs(1,8)=0;
    lhs(1,9)=clhs13;
    lhs(1,10)=0;
    lhs(1,11)=clhs14;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=-clhs20;
    lhs(2,5)=-clhs21;
    lhs(2,6)=clhs20;
    lhs(2,7)=clhs21;
    lhs(2,8)=clhs23;
    lhs(2,9)=0;
    lhs(2,10)=clhs24;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=-clhs26;
    lhs(3,5)=-clhs27;
    lhs(3,6)=clhs26;
    lhs(3,7)=clhs27;
    lhs(3,8)=0;
    lhs(3,9)=clhs23;
    lhs(3,10)=0;
    lhs(3,11)=clhs24;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=clhs29;
    lhs(4,5)=clhs30;
    lhs(4,6)=-clhs29;
    lhs(4,7)=-clhs30;
    lhs(4,8)=clhs32;
    lhs(4,9)=0;
    lhs(4,10)=clhs33;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=clhs35;
    lhs(5,5)=clhs36;
    lhs(5,6)=-clhs35;
    lhs(5,7)=-clhs36;
    lhs(5,8)=0;
    lhs(5,9)=clhs32;
    lhs(5,10)=0;
    lhs(5,11)=clhs33;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=clhs38;
    lhs(6,5)=clhs39;
    lhs(6,6)=-clhs38;
    lhs(6,7)=-clhs39;
    lhs(6,8)=clhs41;
    lhs(6,9)=0;
    lhs(6,10)=clhs42;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=clhs44;
    lhs(7,5)=clhs45;
    lhs(7,6)=-clhs44;
    lhs(7,7)=-clhs45;
    lhs(7,8)=0;
    lhs(7,9)=clhs41;
    lhs(7,10)=0;
    lhs(7,11)=clhs42;
    lhs(8,0)=0;
    lhs(8,1)=0;
    lhs(8,2)=0;
    lhs(8,3)=0;
    lhs(8,4)=0;
    lhs(8,5)=0;
    lhs(8,6)=0;
    lhs(8,7)=0;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=0;
    lhs(9,1)=0;
    lhs(9,2)=0;
    lhs(9,3)=0;
    lhs(9,4)=0;
    lhs(9,5)=0;
    lhs(9,6)=0;
    lhs(9,7)=0;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=0;
    lhs(10,1)=0;
    lhs(10,2)=0;
    lhs(10,3)=0;
    lhs(10,4)=0;
    lhs(10,5)=0;
    lhs(10,6)=0;
    lhs(10,7)=0;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=0;
    lhs(11,1)=0;
    lhs(11,2)=0;
    lhs(11,3)=0;
    lhs(11,4)=0;
    lhs(11,5)=0;
    lhs(11,6)=0;
    lhs(11,7)=0;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHSNormalContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix NormalsMaster = rContactData.NormalsMaster;
    const Vector masternormal  = prod(NormalsMaster, N2);
    const Matrix normalslave   = rContactData.NormalsSlave;
    const Matrix tan1slave     = rContactData.Tangent1Slave;
    const Matrix lm            = rContactData.LagrangeMultipliers;
    const double Dt            = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double clhs0 =             0.5*X1(0,1);
    const double clhs1 =             0.5*X1(1,1);
    const double clhs2 =             0.5*u1(0,1);
    const double clhs3 =             clhs0 - clhs1 + clhs2 - 0.5*u1(1,1);
    const double clhs4 =             std::pow(clhs3, 2);
    const double clhs5 =             N1[0] + N1[1];
    const double clhs6 =             0.5*X1(0,0);
    const double clhs7 =             0.5*X1(1,0);
    const double clhs8 =             0.5*u1(0,0);
    const double clhs9 =             clhs6 - clhs7 + clhs8 - 0.5*u1(1,0);
    const double clhs10 =             std::pow(clhs9, 2);
    const double clhs11 =             1.0/(clhs10 + clhs4);
    const double clhs12 =             std::sqrt(0.25*std::pow(X1(0,0), 2) + 0.25*std::pow(X1(0,1), 2) + 0.25*std::pow(X1(1,0), 2) - X1(1,0)*clhs6 + 0.25*std::pow(X1(1,1), 2) - X1(1,1)*clhs0 + clhs0*u1(0,1) - clhs0*u1(1,1) - clhs1*u1(0,1) + clhs1*u1(1,1) - clhs2*u1(1,1) + clhs6*u1(0,0) - clhs6*u1(1,0) - clhs7*u1(0,0) + clhs7*u1(1,0) - clhs8*u1(1,0) + 0.25*std::pow(u1(0,0), 2) + 0.25*std::pow(u1(0,1), 2) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double clhs13 =             N2[0]*Phi[0]*clhs11*clhs12*clhs5;
    const double clhs14 =             clhs3*clhs9;
    const double clhs15 =             -clhs13*clhs14;
    const double clhs16 =             N2[1]*Phi[0]*clhs11*clhs12*clhs5;
    const double clhs17 =             -clhs14*clhs16;
    const double clhs18 =             Phi[0]*clhs11;
    const double clhs19 =             0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0);
    const double clhs20 =             N1[0]*(X1(0,0) + u1(0,0)) + N1[1]*(X1(1,0) + u1(1,0)) - N2[0]*(X2(0,0) + u2(0,0)) - N2[1]*(X2(1,0) + u2(1,0));
    const double clhs21 =             clhs20*clhs3;
    const double clhs22 =             N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,1) + u1(1,1)) - N2[0]*(X2(0,1) + u2(0,1)) - N2[1]*(X2(1,1) + u2(1,1));
    const double clhs23 =             clhs22*clhs9;
    const double clhs24 =             clhs21 - clhs23;
    const double clhs25 =             clhs24*clhs5/clhs12;
    const double clhs26 =             clhs19*clhs25;
    const double clhs27 =             clhs11*clhs19*clhs5;
    const double clhs28 =             clhs12*clhs24;
    const double clhs29 =             clhs27*clhs28;
    const double clhs30 =             N1[0]*clhs5;
    const double clhs31 =             clhs21*clhs27;
    const double clhs32 =             0.5*N1[0];
    const double clhs33 =             0.5*N1[1];
    const double clhs34 =             clhs32 + clhs33;
    const double clhs35 =             clhs19*clhs9;
    const double clhs36 =             N1[0]*clhs11;
    const double clhs37 =             clhs35*clhs36;
    const double clhs38 =             N1[1]*clhs11;
    const double clhs39 =             clhs35*clhs38;
    const double clhs40 =             clhs12*(clhs22*(clhs34 - clhs37 - clhs39) - clhs3*clhs30 + clhs31);
    const double clhs41 =             clhs3*(-clhs26 + clhs29 + clhs40);
    const double clhs42 =             0.5*clhs12*clhs24*clhs5;
    const double clhs43 =             -clhs42;
    const double clhs44 =             0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1);
    const double clhs45 =             clhs3*clhs44;
    const double clhs46 =             clhs25*clhs45;
    const double clhs47 =             clhs11*clhs12*clhs24*clhs5;
    const double clhs48 =             clhs45*clhs47;
    const double clhs49 =             clhs30*clhs9;
    const double clhs50 =             clhs11*clhs44*clhs5;
    const double clhs51 =             clhs23*clhs50;
    const double clhs52 =             -clhs32 - clhs33;
    const double clhs53 =             clhs36*clhs45;
    const double clhs54 =             clhs38*clhs45;
    const double clhs55 =             clhs20*(clhs52 + clhs53 + clhs54);
    const double clhs56 =             -clhs12*clhs3*(-clhs49 + clhs51 - clhs55) + clhs43 - clhs46 + clhs48;
    const double clhs57 =             N1[1]*clhs5;
    const double clhs58 =             clhs3*clhs57;
    const double clhs59 =             clhs22*(clhs37 + clhs39 + clhs52);
    const double clhs60 =             clhs3*(clhs12*(-clhs31 - clhs58 + clhs59) + clhs26 - clhs29);
    const double clhs61 =             clhs12*(clhs20*(clhs34 - clhs53 - clhs54) + clhs51 + clhs57*clhs9);
    const double clhs62 =             clhs3*clhs61 + clhs42 + clhs46 - clhs48;
    const double clhs63 =             clhs25*clhs35;
    const double clhs64 =             clhs35*clhs47;
    const double clhs65 =             -clhs40*clhs9 + clhs42 + clhs63 - clhs64;
    const double clhs66 =             clhs25*clhs44;
    const double clhs67 =             clhs28*clhs50;
    const double clhs68 =             clhs9*(-clhs12*(clhs49 - clhs51 + clhs55) + clhs66 - clhs67);
    const double clhs69 =             clhs12*clhs9*(clhs31 + clhs58 - clhs59) + clhs43 - clhs63 + clhs64;
    const double clhs70 =             clhs9*(-clhs61 - clhs66 + clhs67);
    const double clhs71 =             N2[0]*Phi[1]*clhs11*clhs12*clhs5;
    const double clhs72 =             -clhs14*clhs71;
    const double clhs73 =             N2[1]*Phi[1]*clhs11*clhs12*clhs5;
    const double clhs74 =             -clhs14*clhs73;
    const double clhs75 =             Phi[1]*clhs11;
    
    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=0;
    lhs(0,9)=0;
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=0;
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=0;
    lhs(2,9)=0;
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=0;
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=0;
    lhs(4,9)=0;
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=0;
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=0;
    lhs(6,9)=0;
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=0;
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(8,0)=clhs13*clhs4;
    lhs(8,1)=clhs15;
    lhs(8,2)=clhs16*clhs4;
    lhs(8,3)=clhs17;
    lhs(8,4)=clhs18*clhs41;
    lhs(8,5)=clhs18*clhs56;
    lhs(8,6)=clhs18*clhs60;
    lhs(8,7)=clhs18*clhs62;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs15;
    lhs(9,1)=clhs10*clhs13;
    lhs(9,2)=clhs17;
    lhs(9,3)=clhs10*clhs16;
    lhs(9,4)=clhs18*clhs65;
    lhs(9,5)=clhs18*clhs68;
    lhs(9,6)=clhs18*clhs69;
    lhs(9,7)=clhs18*clhs70;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs4*clhs71;
    lhs(10,1)=clhs72;
    lhs(10,2)=clhs4*clhs73;
    lhs(10,3)=clhs74;
    lhs(10,4)=clhs41*clhs75;
    lhs(10,5)=clhs56*clhs75;
    lhs(10,6)=clhs60*clhs75;
    lhs(10,7)=clhs62*clhs75;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs72;
    lhs(11,1)=clhs10*clhs71;
    lhs(11,2)=clhs74;
    lhs(11,3)=clhs10*clhs73;
    lhs(11,4)=clhs65*clhs75;
    lhs(11,5)=clhs68*clhs75;
    lhs(11,6)=clhs69*clhs75;
    lhs(11,7)=clhs70*clhs75;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHSTangentContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix NormalsMaster = rContactData.NormalsMaster;
    const Vector masternormal  = prod(NormalsMaster, N2);
    const Matrix normalslave   = rContactData.NormalsSlave;
    const Matrix tan1slave     = rContactData.Tangent1Slave;
    const Matrix lm            = rContactData.LagrangeMultipliers;
    const double Dt            = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double clhs0 =             0.5*X1(0,0);
    const double clhs1 =             0.5*X1(1,0);
    const double clhs2 =             0.5*u1(0,0);
    const double clhs3 =             clhs0 - clhs1 + clhs2 - 0.5*u1(1,0);
    const double clhs4 =             std::pow(clhs3, 2);
    const double clhs5 =             1.0/Dt;
    const double clhs6 =             N1[0] + N1[1];
    const double clhs7 =             0.5*X1(0,1);
    const double clhs8 =             0.5*X1(1,1);
    const double clhs9 =             0.5*u1(0,1);
    const double clhs10 =             clhs7 - clhs8 + clhs9 - 0.5*u1(1,1);
    const double clhs11 =             std::pow(clhs10, 2);
    const double clhs12 =             1.0/(clhs11 + clhs4);
    const double clhs13 =             std::sqrt(0.25*std::pow(X1(0,0), 2) + 0.25*std::pow(X1(0,1), 2) + 0.25*std::pow(X1(1,0), 2) - X1(1,0)*clhs0 + 0.25*std::pow(X1(1,1), 2) - X1(1,1)*clhs7 + clhs0*u1(0,0) - clhs0*u1(1,0) - clhs1*u1(0,0) + clhs1*u1(1,0) - clhs2*u1(1,0) + clhs7*u1(0,1) - clhs7*u1(1,1) - clhs8*u1(0,1) + clhs8*u1(1,1) - clhs9*u1(1,1) + 0.25*std::pow(u1(0,0), 2) + 0.25*std::pow(u1(0,1), 2) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double clhs14 =             N2[0]*Phi[0]*clhs12*clhs13*clhs5*clhs6;
    const double clhs15 =             clhs10*clhs3;
    const double clhs16 =             clhs14*clhs15;
    const double clhs17 =             N2[1]*Phi[0]*clhs12*clhs13*clhs5*clhs6;
    const double clhs18 =             clhs15*clhs17;
    const double clhs19 =             Phi[0]*clhs12*clhs5;
    const double clhs20 =             N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - N2[0]*(Dt*v2(0,0)) - N2[1]*(Dt*v2(1,0));
    const double clhs21 =             clhs20*clhs3;
    const double clhs22 =             N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - N2[0]*(Dt*v2(0,1)) - N2[1]*(Dt*v2(1,1));
    const double clhs23 =             clhs10*clhs22;
    const double clhs24 =             clhs21 + clhs23;
    const double clhs25 =             0.5*clhs13*clhs24*clhs6;
    const double clhs26 =             -clhs25;
    const double clhs27 =             0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0);
    const double clhs28 =             clhs27*clhs3;
    const double clhs29 =             clhs24*clhs6/clhs13;
    const double clhs30 =             clhs28*clhs29;
    const double clhs31 =             clhs12*clhs13*clhs24*clhs6;
    const double clhs32 =             clhs28*clhs31;
    const double clhs33 =             N1[0]*clhs6;
    const double clhs34 =             clhs12*clhs27*clhs6;
    const double clhs35 =             clhs23*clhs34;
    const double clhs36 =             0.5*N1[0] + 0.5*N1[1];
    const double clhs37 =             N1[0]*clhs12;
    const double clhs38 =             N1[1]*clhs12;
    const double clhs39 =             clhs20*(-clhs28*clhs37 - clhs28*clhs38 + clhs36);
    const double clhs40 =             clhs13*(-clhs3*clhs33 + clhs35 - clhs39);
    const double clhs41 =             clhs26 + clhs3*clhs40 - clhs30 + clhs32;
    const double clhs42 =             Phi[0]*clhs12*clhs3*clhs5;
    const double clhs43 =             0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1);
    const double clhs44 =             clhs29*clhs43;
    const double clhs45 =             clhs12*clhs43*clhs6;
    const double clhs46 =             clhs13*clhs24;
    const double clhs47 =             clhs45*clhs46;
    const double clhs48 =             clhs21*clhs45;
    const double clhs49 =             clhs10*clhs43;
    const double clhs50 =             clhs22*(clhs36 - clhs37*clhs49 - clhs38*clhs49);
    const double clhs51 =             clhs13*(-clhs10*clhs33 + clhs48 - clhs50);
    const double clhs52 =             -clhs44 + clhs47 + clhs51;
    const double clhs53 =             N1[1]*clhs6;
    const double clhs54 =             clhs13*(-clhs3*clhs53 - clhs35 + clhs39);
    const double clhs55 =             clhs25 + clhs3*clhs54 + clhs30 - clhs32;
    const double clhs56 =             clhs13*(-clhs10*clhs53 - clhs48 + clhs50);
    const double clhs57 =             clhs44 - clhs47 + clhs56;
    const double clhs58 =             Phi[0]*clhs10*clhs12*clhs5;
    const double clhs59 =             clhs27*clhs29;
    const double clhs60 =             clhs34*clhs46;
    const double clhs61 =             clhs40 - clhs59 + clhs60;
    const double clhs62 =             clhs29*clhs49;
    const double clhs63 =             clhs31*clhs49;
    const double clhs64 =             clhs10*clhs51 + clhs26 - clhs62 + clhs63;
    const double clhs65 =             clhs54 + clhs59 - clhs60;
    const double clhs66 =             clhs10*clhs56 + clhs25 + clhs62 - clhs63;
    const double clhs67 =             N2[0]*Phi[1]*clhs12*clhs13*clhs5*clhs6;
    const double clhs68 =             clhs15*clhs67;
    const double clhs69 =             N2[1]*Phi[1]*clhs12*clhs13*clhs5*clhs6;
    const double clhs70 =             clhs15*clhs69;
    const double clhs71 =             Phi[1]*clhs12*clhs5;
    const double clhs72 =             Phi[1]*clhs12*clhs3*clhs5;
    const double clhs73 =             Phi[1]*clhs10*clhs12*clhs5;

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=0;
    lhs(0,5)=0;
    lhs(0,6)=0;
    lhs(0,7)=0;
    lhs(0,8)=0;
    lhs(0,9)=0;
    lhs(0,10)=0;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=0;
    lhs(1,5)=0;
    lhs(1,6)=0;
    lhs(1,7)=0;
    lhs(1,8)=0;
    lhs(1,9)=0;
    lhs(1,10)=0;
    lhs(1,11)=0;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=0;
    lhs(2,5)=0;
    lhs(2,6)=0;
    lhs(2,7)=0;
    lhs(2,8)=0;
    lhs(2,9)=0;
    lhs(2,10)=0;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=0;
    lhs(3,5)=0;
    lhs(3,6)=0;
    lhs(3,7)=0;
    lhs(3,8)=0;
    lhs(3,9)=0;
    lhs(3,10)=0;
    lhs(3,11)=0;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=0;
    lhs(4,5)=0;
    lhs(4,6)=0;
    lhs(4,7)=0;
    lhs(4,8)=0;
    lhs(4,9)=0;
    lhs(4,10)=0;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=0;
    lhs(5,5)=0;
    lhs(5,6)=0;
    lhs(5,7)=0;
    lhs(5,8)=0;
    lhs(5,9)=0;
    lhs(5,10)=0;
    lhs(5,11)=0;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=0;
    lhs(6,5)=0;
    lhs(6,6)=0;
    lhs(6,7)=0;
    lhs(6,8)=0;
    lhs(6,9)=0;
    lhs(6,10)=0;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=0;
    lhs(7,5)=0;
    lhs(7,6)=0;
    lhs(7,7)=0;
    lhs(7,8)=0;
    lhs(7,9)=0;
    lhs(7,10)=0;
    lhs(7,11)=0;
    lhs(8,0)=clhs14*clhs4;
    lhs(8,1)=clhs16;
    lhs(8,2)=clhs17*clhs4;
    lhs(8,3)=clhs18;
    lhs(8,4)=clhs19*clhs41;
    lhs(8,5)=clhs42*clhs52;
    lhs(8,6)=clhs19*clhs55;
    lhs(8,7)=clhs42*clhs57;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs16;
    lhs(9,1)=clhs11*clhs14;
    lhs(9,2)=clhs18;
    lhs(9,3)=clhs11*clhs17;
    lhs(9,4)=clhs58*clhs61;
    lhs(9,5)=clhs19*clhs64;
    lhs(9,6)=clhs58*clhs65;
    lhs(9,7)=clhs19*clhs66;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs4*clhs67;
    lhs(10,1)=clhs68;
    lhs(10,2)=clhs4*clhs69;
    lhs(10,3)=clhs70;
    lhs(10,4)=clhs41*clhs71;
    lhs(10,5)=clhs52*clhs72;
    lhs(10,6)=clhs55*clhs71;
    lhs(10,7)=clhs57*clhs72;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs68;
    lhs(11,1)=clhs11*clhs67;
    lhs(11,2)=clhs70;
    lhs(11,3)=clhs11*clhs69;
    lhs(11,4)=clhs61*clhs73;
    lhs(11,5)=clhs64*clhs71;
    lhs(11,6)=clhs65*clhs73;
    lhs(11,7)=clhs66*clhs71;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
        
    static inline array_1d<double,12> ComputeGaussPointRHSInternalContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;

    const Matrix lm           = rContactData.LagrangeMultipliers;

    const double crhs0 =             N2[0]*detJ;
    const double crhs1 =             Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double crhs2 =             Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double crhs3 =             N2[1]*detJ;
    const double crhs4 =             N1[0]*detJ;
    const double crhs5 =             N1[1]*detJ;
    
    rhs[0]=crhs0*crhs1;
    rhs[1]=crhs0*crhs2;
    rhs[2]=crhs1*crhs3;
    rhs[3]=crhs2*crhs3;
    rhs[4]=-crhs1*crhs4;
    rhs[5]=-crhs2*crhs4;
    rhs[6]=-crhs1*crhs5;
    rhs[7]=-crhs2*crhs5;
    rhs[8]=0;
    rhs[9]=0;
    rhs[10]=0;
    rhs[11]=0;

    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointRHSNormalContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Vector Gaps         = rContactData.Gaps;
    const double Dt           = rContactData.Dt;
    
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double crhs0 =             - inner_prod(Gaps, N1);
    const double crhs1 =             Phi[0]*crhs0*detJ;
    const double crhs2 =             Phi[1]*crhs0*detJ;
    
    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs1*normalslave(0,0);
    rhs[9]=crhs1*normalslave(0,1);
    rhs[10]=crhs2*normalslave(1,0);
    rhs[11]=crhs2*normalslave(1,1);
    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointRHSTangentContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Vector Gaps         = rContactData.Gaps;
    const double Dt           = rContactData.Dt;
    
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double crhs0 =             1.0/Dt;
    const double crhs1 =             (N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0))*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - N2[0]*(Dt*v2(0,0)) - N2[1]*(Dt*v2(1,0))) + (N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1))*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - N2[0]*(Dt*v2(0,1)) - N2[1]*(Dt*v2(1,1)));
    const double crhs2 =             Phi[0]*crhs0*crhs1*detJ;
    const double crhs3 =             Phi[1]*crhs0*crhs1*detJ;
    
    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs2*tan1slave(0,0);
    rhs[9]=crhs2*tan1slave(0,1);
    rhs[10]=crhs3*tan1slave(1,0);
    rhs[11]=crhs3*tan1slave(1,1);
    
    return rhs;
}

private:
    
    static inline Matrix GetCoordinates(
        const GeometryType& nodes,
        const bool current
        )
{
    /* DEFINITIONS */    
    const unsigned int dimension = nodes.WorkingSpaceDimension();
    const unsigned int number_of_nodes  =  nodes.PointsNumber( );
    
    Matrix Coordinates(number_of_nodes, dimension);
    
    for (unsigned int iNode = 0; iNode < number_of_nodes; iNode++)
    {
        array_1d<double, 3> coord = nodes[iNode].Coordinates();
        
        if (current == false)
        {
            coord -= nodes[iNode].FastGetSolutionStepValue(DISPLACEMENT, 0);
        }

        for (unsigned int iDof = 0; iDof < dimension; iDof++)
        {
            Coordinates(iNode, iDof) = coord[iDof];
        }
    }
    
    return Coordinates;
}

    /***********************************************************************************/
    /***********************************************************************************/

    static inline Matrix GetVariableMatrix(
        const GeometryType& nodes,
        const Variable<array_1d<double,3> >& rVarName,
        unsigned int step
        )
{
    /* DEFINITIONS */
    const unsigned int dimension = nodes.WorkingSpaceDimension();
    const unsigned int number_of_nodes  =  nodes.PointsNumber( );
    
    Matrix VarMatrix(number_of_nodes, dimension);
    
    for (unsigned int iNode = 0; iNode < number_of_nodes; iNode++)
    {
        array_1d<double, 3> Value = nodes[iNode].FastGetSolutionStepValue(rVarName, step);
        for (unsigned int iDof = 0; iDof < dimension; iDof++)
        {
            VarMatrix(iNode, iDof) = Value[iDof];
        }
    }
    
    return VarMatrix;
}

};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */