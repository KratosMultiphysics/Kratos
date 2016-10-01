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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
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
    const double clhs4 =             GPnormal[0]*clhs3;
    const double clhs5 =             0.5*X1(0,0);
    const double clhs6 =             0.5*X1(1,0);
    const double clhs7 =             0.5*u1(0,0);
    const double clhs8 =             clhs5 - clhs6 + clhs7 - 0.5*u1(1,0);
    const double clhs9 =             GPnormal[1]*clhs8;
    const double clhs10 =             clhs4 - clhs9;
    const double clhs11 =             N1[0] + N1[1];
    const double clhs12 =             1.0/(std::pow(clhs3, 2) + std::pow(clhs8, 2));
    const double clhs13 =             std::sqrt(0.25*std::pow(X1(0,0), 2) + 0.25*std::pow(X1(0,1), 2) + 0.25*std::pow(X1(1,0), 2) - X1(1,0)*clhs5 + 0.25*std::pow(X1(1,1), 2) - X1(1,1)*clhs0 + clhs0*u1(0,1) - clhs0*u1(1,1) - clhs1*u1(0,1) + clhs1*u1(1,1) - clhs2*u1(1,1) + clhs5*u1(0,0) - clhs5*u1(1,0) - clhs6*u1(0,0) + clhs6*u1(1,0) - clhs7*u1(1,0) + 0.25*std::pow(u1(0,0), 2) + 0.25*std::pow(u1(0,1), 2) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double clhs14 =             N2[0]*Phi[0]*clhs11*clhs12*clhs13*clhs3;
    const double clhs15 =             -clhs4 + clhs9;
    const double clhs16 =             N2[0]*Phi[0]*clhs11*clhs12*clhs13*clhs8;
    const double clhs17 =             N2[1]*Phi[0]*clhs11*clhs12*clhs13*clhs3;
    const double clhs18 =             N2[1]*Phi[0]*clhs11*clhs12*clhs13*clhs8;
    const double clhs19 =             Phi[0]*clhs12;
    const double clhs20 =             0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0);
    const double clhs21 =             N1[0]*(X1(0,0) + u1(0,0)) + N1[1]*(X1(1,0) + u1(1,0)) - N2[0]*(X2(0,0) + u2(0,0)) - N2[1]*(X2(1,0) + u2(1,0));
    const double clhs22 =             clhs21*clhs3;
    const double clhs23 =             N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,1) + u1(1,1)) - N2[0]*(X2(0,1) + u2(0,1)) - N2[1]*(X2(1,1) + u2(1,1));
    const double clhs24 =             clhs23*clhs8;
    const double clhs25 =             clhs22 - clhs24;
    const double clhs26 =             clhs11*clhs25/clhs13;
    const double clhs27 =             clhs20*clhs26;
    const double clhs28 =             clhs11*clhs12*clhs20;
    const double clhs29 =             clhs13*clhs25;
    const double clhs30 =             clhs28*clhs29;
    const double clhs31 =             N1[0]*clhs11;
    const double clhs32 =             clhs22*clhs28;
    const double clhs33 =             0.5*N1[0];
    const double clhs34 =             0.5*N1[1];
    const double clhs35 =             clhs33 + clhs34;
    const double clhs36 =             clhs20*clhs8;
    const double clhs37 =             N1[0]*clhs12;
    const double clhs38 =             clhs36*clhs37;
    const double clhs39 =             N1[1]*clhs12;
    const double clhs40 =             clhs36*clhs39;
    const double clhs41 =             clhs13*(clhs23*(clhs35 - clhs38 - clhs40) - clhs3*clhs31 + clhs32);
    const double clhs42 =             -clhs27 + clhs30 + clhs41;
    const double clhs43 =             0.5*clhs11*clhs13*clhs25;
    const double clhs44 =             clhs11*clhs12*clhs13*clhs25;
    const double clhs45 =             clhs26*clhs36 - clhs36*clhs44 + clhs43;
    const double clhs46 =             -clhs41*clhs8 + clhs45;
    const double clhs47 =             GPnormal[1]*clhs46 + clhs4*clhs42;
    const double clhs48 =             clhs31*clhs8;
    const double clhs49 =             0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1);
    const double clhs50 =             clhs11*clhs12*clhs49;
    const double clhs51 =             clhs24*clhs50;
    const double clhs52 =             -clhs33 - clhs34;
    const double clhs53 =             clhs3*clhs49;
    const double clhs54 =             clhs37*clhs53;
    const double clhs55 =             clhs39*clhs53;
    const double clhs56 =             clhs21*(clhs52 + clhs54 + clhs55);
    const double clhs57 =             clhs26*clhs49 - clhs29*clhs50;
    const double clhs58 =             -clhs13*(clhs48 - clhs51 + clhs56) + clhs57;
    const double clhs59 =             clhs26*clhs53 + clhs43 - clhs44*clhs53;
    const double clhs60 =             clhs13*clhs3*(-clhs48 + clhs51 - clhs56) + clhs59;
    const double clhs61 =             -GPnormal[0]*clhs60 + clhs58*clhs9;
    const double clhs62 =             N1[1]*clhs11;
    const double clhs63 =             clhs3*clhs62;
    const double clhs64 =             clhs23*(clhs38 + clhs40 + clhs52);
    const double clhs65 =             clhs13*(-clhs32 - clhs63 + clhs64) + clhs27 - clhs30;
    const double clhs66 =             -clhs13*clhs8*(clhs32 + clhs63 - clhs64) + clhs45;
    const double clhs67 =             -GPnormal[1]*clhs66 + clhs4*clhs65;
    const double clhs68 =             clhs13*(clhs21*(clhs35 - clhs54 - clhs55) + clhs51 + clhs62*clhs8);
    const double clhs69 =             clhs57 + clhs68;
    const double clhs70 =             clhs3*clhs68 + clhs59;
    const double clhs71 =             GPnormal[0]*clhs70 - clhs69*clhs9;
    const double clhs72 =             GPtangent1[0]*clhs3;
    const double clhs73 =             GPtangent1[1]*clhs8;
    const double clhs74 =             clhs72 - clhs73;
    const double clhs75 =             -clhs72 + clhs73;
    const double clhs76 =             GPtangent1[1]*clhs46 + clhs42*clhs72;
    const double clhs77 =             -GPtangent1[0]*clhs60 + clhs58*clhs73;
    const double clhs78 =             -GPtangent1[1]*clhs66 + clhs65*clhs72;
    const double clhs79 =             GPtangent1[0]*clhs70 - clhs69*clhs73;
    const double clhs80 =             N2[0]*Phi[1]*clhs11*clhs12*clhs13*clhs3;
    const double clhs81 =             N2[0]*Phi[1]*clhs11*clhs12*clhs13*clhs8;
    const double clhs82 =             N2[1]*Phi[1]*clhs11*clhs12*clhs13*clhs3;
    const double clhs83 =             N2[1]*Phi[1]*clhs11*clhs12*clhs13*clhs8;
    const double clhs84 =             Phi[1]*clhs12;
    
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
    lhs(8,0)=clhs10*clhs14;
    lhs(8,1)=clhs15*clhs16;
    lhs(8,2)=clhs10*clhs17;
    lhs(8,3)=clhs15*clhs18;
    lhs(8,4)=clhs19*clhs47;
    lhs(8,5)=clhs19*clhs61;
    lhs(8,6)=clhs19*clhs67;
    lhs(8,7)=clhs19*clhs71;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs14*clhs74;
    lhs(9,1)=clhs16*clhs75;
    lhs(9,2)=clhs17*clhs74;
    lhs(9,3)=clhs18*clhs75;
    lhs(9,4)=clhs19*clhs76;
    lhs(9,5)=clhs19*clhs77;
    lhs(9,6)=clhs19*clhs78;
    lhs(9,7)=clhs19*clhs79;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs10*clhs80;
    lhs(10,1)=clhs15*clhs81;
    lhs(10,2)=clhs10*clhs82;
    lhs(10,3)=clhs15*clhs83;
    lhs(10,4)=clhs47*clhs84;
    lhs(10,5)=clhs61*clhs84;
    lhs(10,6)=clhs67*clhs84;
    lhs(10,7)=clhs71*clhs84;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs74*clhs80;
    lhs(11,1)=clhs75*clhs81;
    lhs(11,2)=clhs74*clhs82;
    lhs(11,3)=clhs75*clhs83;
    lhs(11,4)=clhs76*clhs84;
    lhs(11,5)=clhs77*clhs84;
    lhs(11,6)=clhs78*clhs84;
    lhs(11,7)=clhs79*clhs84;
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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
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
    const double clhs4 =             1.0/Dt;
    const double clhs5 =             N1[0] + N1[1];
    const double clhs6 =             GPnormal[0]*clhs3;
    const double clhs7 =             0.5*X1(0,1);
    const double clhs8 =             0.5*X1(1,1);
    const double clhs9 =             0.5*u1(0,1);
    const double clhs10 =             clhs7 - clhs8 + clhs9 - 0.5*u1(1,1);
    const double clhs11 =             GPnormal[1]*clhs10;
    const double clhs12 =             clhs11 + clhs6;
    const double clhs13 =             1.0/(std::pow(clhs10, 2) + std::pow(clhs3, 2));
    const double clhs14 =             std::sqrt(0.25*std::pow(X1(0,0), 2) + 0.25*std::pow(X1(0,1), 2) + 0.25*std::pow(X1(1,0), 2) - X1(1,0)*clhs0 + 0.25*std::pow(X1(1,1), 2) - X1(1,1)*clhs7 + clhs0*u1(0,0) - clhs0*u1(1,0) - clhs1*u1(0,0) + clhs1*u1(1,0) - clhs2*u1(1,0) + clhs7*u1(0,1) - clhs7*u1(1,1) - clhs8*u1(0,1) + clhs8*u1(1,1) - clhs9*u1(1,1) + 0.25*std::pow(u1(0,0), 2) + 0.25*std::pow(u1(0,1), 2) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double clhs15 =             N2[0]*Phi[0]*clhs12*clhs13*clhs14*clhs4*clhs5;
    const double clhs16 =             N2[1]*Phi[0]*clhs12*clhs13*clhs14*clhs4*clhs5;
    const double clhs17 =             Phi[0]*clhs13*clhs4;
    const double clhs18 =             0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0);
    const double clhs19 =             N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - N2[0]*(Dt*v2(0,0)) - N2[1]*(Dt*v2(1,0));
    const double clhs20 =             clhs19*clhs3;
    const double clhs21 =             N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - N2[0]*(Dt*v2(0,1)) - N2[1]*(Dt*v2(1,1));
    const double clhs22 =             clhs10*clhs21;
    const double clhs23 =             clhs20 + clhs22;
    const double clhs24 =             clhs23*clhs5/clhs14;
    const double clhs25 =             clhs18*clhs24;
    const double clhs26 =             clhs13*clhs18*clhs5;
    const double clhs27 =             clhs14*clhs23;
    const double clhs28 =             clhs26*clhs27;
    const double clhs29 =             N1[0]*clhs5;
    const double clhs30 =             clhs22*clhs26;
    const double clhs31 =             -clhs29*clhs3 + clhs30;
    const double clhs32 =             0.5*N1[0];
    const double clhs33 =             0.5*N1[1];
    const double clhs34 =             -clhs32 - clhs33;
    const double clhs35 =             clhs18*clhs3;
    const double clhs36 =             N1[0]*clhs13;
    const double clhs37 =             clhs35*clhs36;
    const double clhs38 =             N1[1]*clhs13;
    const double clhs39 =             clhs35*clhs38;
    const double clhs40 =             clhs14*(clhs19*(clhs34 + clhs37 + clhs39) + clhs31) - clhs25 + clhs28;
    const double clhs41 =             0.5*clhs14*clhs23*clhs5;
    const double clhs42 =             -clhs41;
    const double clhs43 =             clhs24*clhs35;
    const double clhs44 =             clhs13*clhs14*clhs23*clhs5;
    const double clhs45 =             clhs35*clhs44;
    const double clhs46 =             clhs32 + clhs33;
    const double clhs47 =             clhs19*(-clhs37 - clhs39 + clhs46);
    const double clhs48 =             clhs14*clhs3*(clhs31 - clhs47) + clhs42 - clhs43 + clhs45;
    const double clhs49 =             GPnormal[0]*clhs48 + clhs11*clhs40;
    const double clhs50 =             0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1);
    const double clhs51 =             clhs24*clhs50;
    const double clhs52 =             clhs13*clhs5*clhs50;
    const double clhs53 =             clhs27*clhs52;
    const double clhs54 =             clhs20*clhs52;
    const double clhs55 =             -clhs10*clhs29 + clhs54;
    const double clhs56 =             clhs10*clhs50;
    const double clhs57 =             clhs36*clhs56;
    const double clhs58 =             clhs38*clhs56;
    const double clhs59 =             clhs14*(clhs21*(clhs34 + clhs57 + clhs58) + clhs55) - clhs51 + clhs53;
    const double clhs60 =             clhs24*clhs56;
    const double clhs61 =             clhs44*clhs56;
    const double clhs62 =             clhs21*(clhs46 - clhs57 - clhs58);
    const double clhs63 =             clhs10*clhs14*(clhs55 - clhs62) + clhs42 - clhs60 + clhs61;
    const double clhs64 =             GPnormal[1]*clhs63 + clhs59*clhs6;
    const double clhs65 =             N1[1]*clhs5;
    const double clhs66 =             clhs14*(-clhs3*clhs65 - clhs30 + clhs47);
    const double clhs67 =             clhs25 - clhs28 + clhs66;
    const double clhs68 =             clhs3*clhs66 + clhs41 + clhs43 - clhs45;
    const double clhs69 =             GPnormal[0]*clhs68 + clhs11*clhs67;
    const double clhs70 =             clhs14*(-clhs10*clhs65 - clhs54 + clhs62);
    const double clhs71 =             clhs51 - clhs53 + clhs70;
    const double clhs72 =             clhs10*clhs70 + clhs41 + clhs60 - clhs61;
    const double clhs73 =             GPnormal[1]*clhs72 + clhs6*clhs71;
    const double clhs74 =             GPtangent1[0]*clhs3;
    const double clhs75 =             GPtangent1[1]*clhs10;
    const double clhs76 =             clhs74 + clhs75;
    const double clhs77 =             N2[0]*Phi[0]*clhs13*clhs14*clhs4*clhs5*clhs76;
    const double clhs78 =             N2[1]*Phi[0]*clhs13*clhs14*clhs4*clhs5*clhs76;
    const double clhs79 =             GPtangent1[0]*clhs48 + clhs40*clhs75;
    const double clhs80 =             GPtangent1[1]*clhs63 + clhs59*clhs74;
    const double clhs81 =             GPtangent1[0]*clhs68 + clhs67*clhs75;
    const double clhs82 =             GPtangent1[1]*clhs72 + clhs71*clhs74;
    const double clhs83 =             N2[0]*Phi[1]*clhs12*clhs13*clhs14*clhs4*clhs5;
    const double clhs84 =             N2[1]*Phi[1]*clhs12*clhs13*clhs14*clhs4*clhs5;
    const double clhs85 =             Phi[1]*clhs13*clhs4;
    const double clhs86 =             N2[0]*Phi[1]*clhs13*clhs14*clhs4*clhs5*clhs76;
    const double clhs87 =             N2[1]*Phi[1]*clhs13*clhs14*clhs4*clhs5*clhs76;
           
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
    lhs(8,0)=clhs15*clhs3;
    lhs(8,1)=clhs10*clhs15;
    lhs(8,2)=clhs16*clhs3;
    lhs(8,3)=clhs10*clhs16;
    lhs(8,4)=clhs17*clhs49;
    lhs(8,5)=clhs17*clhs64;
    lhs(8,6)=clhs17*clhs69;
    lhs(8,7)=clhs17*clhs73;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs3*clhs77;
    lhs(9,1)=clhs10*clhs77;
    lhs(9,2)=clhs3*clhs78;
    lhs(9,3)=clhs10*clhs78;
    lhs(9,4)=clhs17*clhs79;
    lhs(9,5)=clhs17*clhs80;
    lhs(9,6)=clhs17*clhs81;
    lhs(9,7)=clhs17*clhs82;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs3*clhs83;
    lhs(10,1)=clhs10*clhs83;
    lhs(10,2)=clhs3*clhs84;
    lhs(10,3)=clhs10*clhs84;
    lhs(10,4)=clhs49*clhs85;
    lhs(10,5)=clhs64*clhs85;
    lhs(10,6)=clhs69*clhs85;
    lhs(10,7)=clhs73*clhs85;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs3*clhs86;
    lhs(11,1)=clhs10*clhs86;
    lhs(11,2)=clhs3*clhs87;
    lhs(11,3)=clhs10*clhs87;
    lhs(11,4)=clhs79*clhs85;
    lhs(11,5)=clhs80*clhs85;
    lhs(11,6)=clhs81*clhs85;
    lhs(11,7)=clhs82*clhs85;
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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
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
    rhs[8]=crhs1*(GPnormal[0]*normalslave(0,0) + GPnormal[1]*normalslave(0,1));
    rhs[9]=crhs1*(GPtangent1[0]*normalslave(0,0) + GPtangent1[1]*normalslave(0,1));
    rhs[10]=crhs2*(GPnormal[0]*normalslave(1,0) + GPnormal[1]*normalslave(1,1));
    rhs[11]=crhs2*(GPtangent1[0]*normalslave(1,0) + GPtangent1[1]*normalslave(1,1));
    
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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
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
    rhs[8]=crhs2*(GPnormal[0]*tan1slave(0,0) + GPnormal[1]*tan1slave(0,1));
    rhs[9]=crhs2*(GPtangent1[0]*tan1slave(0,0) + GPtangent1[1]*tan1slave(0,1));
    rhs[10]=crhs3*(GPnormal[0]*tan1slave(1,0) + GPnormal[1]*tan1slave(1,1));
    rhs[11]=crhs3*(GPtangent1[0]*tan1slave(1,0) + GPtangent1[1]*tan1slave(1,1));

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