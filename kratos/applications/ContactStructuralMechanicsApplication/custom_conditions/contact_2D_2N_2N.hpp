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
        const Matrix DN1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;

    const Matrix lm           = rContactData.LagrangeMultipliers;
    const Matrix masternormal  = rContactData.NormalsMaster;
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    
    const double clhs0 =             N2[0]*detJ;
const double clhs1 =             Phi[0]*clhs0;
const double clhs2 =             Phi[1]*clhs0;
const double clhs3 =             N2[1]*detJ;
const double clhs4 =             Phi[0]*clhs3;
const double clhs5 =             Phi[1]*clhs3;
const double clhs6 =             N1[0]*detJ;
const double clhs7 =             -Phi[0]*clhs6;
const double clhs8 =             -Phi[1]*clhs6;
const double clhs9 =             N1[1]*detJ;
const double clhs10 =             -Phi[0]*clhs9;
const double clhs11 =             -Phi[1]*clhs9;
            lhs(0,0)=0;
            lhs(0,1)=0;
            lhs(0,2)=0;
            lhs(0,3)=0;
            lhs(0,4)=0;
            lhs(0,5)=0;
            lhs(0,6)=0;
            lhs(0,7)=0;
            lhs(0,8)=clhs1;
            lhs(0,9)=0;
            lhs(0,10)=clhs2;
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
            lhs(1,9)=clhs1;
            lhs(1,10)=0;
            lhs(1,11)=clhs2;
            lhs(2,0)=0;
            lhs(2,1)=0;
            lhs(2,2)=0;
            lhs(2,3)=0;
            lhs(2,4)=0;
            lhs(2,5)=0;
            lhs(2,6)=0;
            lhs(2,7)=0;
            lhs(2,8)=clhs4;
            lhs(2,9)=0;
            lhs(2,10)=clhs5;
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
            lhs(3,9)=clhs4;
            lhs(3,10)=0;
            lhs(3,11)=clhs5;
            lhs(4,0)=0;
            lhs(4,1)=0;
            lhs(4,2)=0;
            lhs(4,3)=0;
            lhs(4,4)=0;
            lhs(4,5)=0;
            lhs(4,6)=0;
            lhs(4,7)=0;
            lhs(4,8)=clhs7;
            lhs(4,9)=0;
            lhs(4,10)=clhs8;
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
            lhs(5,9)=clhs7;
            lhs(5,10)=0;
            lhs(5,11)=clhs8;
            lhs(6,0)=0;
            lhs(6,1)=0;
            lhs(6,2)=0;
            lhs(6,3)=0;
            lhs(6,4)=0;
            lhs(6,5)=0;
            lhs(6,6)=0;
            lhs(6,7)=0;
            lhs(6,8)=clhs10;
            lhs(6,9)=0;
            lhs(6,10)=clhs11;
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
            lhs(7,9)=clhs10;
            lhs(7,10)=0;
            lhs(7,11)=clhs11;
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
        const Matrix DN1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix masternormal  = rContactData.NormalsMaster;
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const double Dt           = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double clhs0 =             N2[0]*Phi[0]*detJ*normalslave(0,0);
const double clhs1 =             N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0);
const double clhs2 =             N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1);
const double clhs3 =             N2[1]*Phi[0]*detJ*normalslave(0,0);
const double clhs4 =             N1[0]*Phi[0]*detJ*normalslave(0,0);
const double clhs5 =             N1[1]*Phi[0]*detJ*normalslave(0,0);
const double clhs6 =             N2[0]*Phi[0]*detJ*normalslave(0,1);
const double clhs7 =             N2[1]*Phi[0]*detJ*normalslave(0,1);
const double clhs8 =             N1[0]*Phi[0]*detJ*normalslave(0,1);
const double clhs9 =             N1[1]*Phi[0]*detJ*normalslave(0,1);
const double clhs10 =             N2[0]*Phi[1]*detJ*normalslave(1,0);
const double clhs11 =             N2[1]*Phi[1]*detJ*normalslave(1,0);
const double clhs12 =             N1[0]*Phi[1]*detJ*normalslave(1,0);
const double clhs13 =             N1[1]*Phi[1]*detJ*normalslave(1,0);
const double clhs14 =             N2[0]*Phi[1]*detJ*normalslave(1,1);
const double clhs15 =             N2[1]*Phi[1]*detJ*normalslave(1,1);
const double clhs16 =             N1[0]*Phi[1]*detJ*normalslave(1,1);
const double clhs17 =             N1[1]*Phi[1]*detJ*normalslave(1,1);
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
            lhs(8,0)=clhs0*clhs1;
            lhs(8,1)=clhs0*clhs2;
            lhs(8,2)=clhs1*clhs3;
            lhs(8,3)=clhs2*clhs3;
            lhs(8,4)=-clhs1*clhs4;
            lhs(8,5)=-clhs2*clhs4;
            lhs(8,6)=-clhs1*clhs5;
            lhs(8,7)=-clhs2*clhs5;
            lhs(8,8)=0;
            lhs(8,9)=0;
            lhs(8,10)=0;
            lhs(8,11)=0;
            lhs(9,0)=clhs1*clhs6;
            lhs(9,1)=clhs2*clhs6;
            lhs(9,2)=clhs1*clhs7;
            lhs(9,3)=clhs2*clhs7;
            lhs(9,4)=-clhs1*clhs8;
            lhs(9,5)=-clhs2*clhs8;
            lhs(9,6)=-clhs1*clhs9;
            lhs(9,7)=-clhs2*clhs9;
            lhs(9,8)=0;
            lhs(9,9)=0;
            lhs(9,10)=0;
            lhs(9,11)=0;
            lhs(10,0)=clhs1*clhs10;
            lhs(10,1)=clhs10*clhs2;
            lhs(10,2)=clhs1*clhs11;
            lhs(10,3)=clhs11*clhs2;
            lhs(10,4)=-clhs1*clhs12;
            lhs(10,5)=-clhs12*clhs2;
            lhs(10,6)=-clhs1*clhs13;
            lhs(10,7)=-clhs13*clhs2;
            lhs(10,8)=0;
            lhs(10,9)=0;
            lhs(10,10)=0;
            lhs(10,11)=0;
            lhs(11,0)=clhs1*clhs14;
            lhs(11,1)=clhs14*clhs2;
            lhs(11,2)=clhs1*clhs15;
            lhs(11,3)=clhs15*clhs2;
            lhs(11,4)=-clhs1*clhs16;
            lhs(11,5)=-clhs16*clhs2;
            lhs(11,6)=-clhs1*clhs17;
            lhs(11,7)=-clhs17*clhs2;
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
        const Matrix DN1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix masternormal  = rContactData.NormalsMaster;
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const double Dt           = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

    const double clhs0 =             N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0);
const double clhs1 =             1.0/Dt;
const double clhs2 =             N2[0]*Phi[0]*clhs1*detJ*tan1slave(0,0);
const double clhs3 =             N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1);
const double clhs4 =             N2[1]*Phi[0]*clhs1*detJ*tan1slave(0,0);
const double clhs5 =             N1[0]*Phi[0]*clhs1*detJ*tan1slave(0,0);
const double clhs6 =             N1[1]*Phi[0]*clhs1*detJ*tan1slave(0,0);
const double clhs7 =             N2[0]*Phi[0]*clhs1*detJ*tan1slave(0,1);
const double clhs8 =             N2[1]*Phi[0]*clhs1*detJ*tan1slave(0,1);
const double clhs9 =             N1[0]*Phi[0]*clhs1*detJ*tan1slave(0,1);
const double clhs10 =             N1[1]*Phi[0]*clhs1*detJ*tan1slave(0,1);
const double clhs11 =             N2[0]*Phi[1]*clhs1*detJ*tan1slave(1,0);
const double clhs12 =             N2[1]*Phi[1]*clhs1*detJ*tan1slave(1,0);
const double clhs13 =             N1[0]*Phi[1]*clhs1*detJ*tan1slave(1,0);
const double clhs14 =             N1[1]*Phi[1]*clhs1*detJ*tan1slave(1,0);
const double clhs15 =             N2[0]*Phi[1]*clhs1*detJ*tan1slave(1,1);
const double clhs16 =             N2[1]*Phi[1]*clhs1*detJ*tan1slave(1,1);
const double clhs17 =             N1[0]*Phi[1]*clhs1*detJ*tan1slave(1,1);
const double clhs18 =             N1[1]*Phi[1]*clhs1*detJ*tan1slave(1,1);
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
            lhs(8,0)=clhs0*clhs2;
            lhs(8,1)=clhs2*clhs3;
            lhs(8,2)=clhs0*clhs4;
            lhs(8,3)=clhs3*clhs4;
            lhs(8,4)=-clhs0*clhs5;
            lhs(8,5)=-clhs3*clhs5;
            lhs(8,6)=-clhs0*clhs6;
            lhs(8,7)=-clhs3*clhs6;
            lhs(8,8)=0;
            lhs(8,9)=0;
            lhs(8,10)=0;
            lhs(8,11)=0;
            lhs(9,0)=clhs0*clhs7;
            lhs(9,1)=clhs3*clhs7;
            lhs(9,2)=clhs0*clhs8;
            lhs(9,3)=clhs3*clhs8;
            lhs(9,4)=-clhs0*clhs9;
            lhs(9,5)=-clhs3*clhs9;
            lhs(9,6)=-clhs0*clhs10;
            lhs(9,7)=-clhs10*clhs3;
            lhs(9,8)=0;
            lhs(9,9)=0;
            lhs(9,10)=0;
            lhs(9,11)=0;
            lhs(10,0)=clhs0*clhs11;
            lhs(10,1)=clhs11*clhs3;
            lhs(10,2)=clhs0*clhs12;
            lhs(10,3)=clhs12*clhs3;
            lhs(10,4)=-clhs0*clhs13;
            lhs(10,5)=-clhs13*clhs3;
            lhs(10,6)=-clhs0*clhs14;
            lhs(10,7)=-clhs14*clhs3;
            lhs(10,8)=0;
            lhs(10,9)=0;
            lhs(10,10)=0;
            lhs(10,11)=0;
            lhs(11,0)=clhs0*clhs15;
            lhs(11,1)=clhs15*clhs3;
            lhs(11,2)=clhs0*clhs16;
            lhs(11,3)=clhs16*clhs3;
            lhs(11,4)=-clhs0*clhs17;
            lhs(11,5)=-clhs17*clhs3;
            lhs(11,6)=-clhs0*clhs18;
            lhs(11,7)=-clhs18*clhs3;
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
            rhs[0]=-crhs0*crhs1;
            rhs[1]=-crhs0*crhs2;
            rhs[2]=-crhs1*crhs3;
            rhs[3]=-crhs2*crhs3;
            rhs[4]=crhs1*crhs4;
            rhs[5]=crhs2*crhs4;
            rhs[6]=crhs1*crhs5;
            rhs[7]=crhs2*crhs5;
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