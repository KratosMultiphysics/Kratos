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
        const double J, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;

    const double clhs0 =             J*N2[0];
    const double clhs1 =             Phi[0]*clhs0;
    const double clhs2 =             Phi[1]*clhs0;
    const double clhs3 =             J*N2[1];
    const double clhs4 =             Phi[0]*clhs3;
    const double clhs5 =             Phi[1]*clhs3;
    const double clhs6 =             J*N1[0];
    const double clhs7 =             -Phi[0]*clhs6;
    const double clhs8 =             -Phi[1]*clhs6;
    const double clhs9 =             J*N1[1];
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
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHSContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double J, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster = rContactData.NormalsMaster;
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const double Dt           = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData, rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData, rContactData.MasterGeometry, false);
    const Matrix x1 = GetCoordinates(rContactData, rContactData.SlaveGeometry, true);
    const Matrix x2 = GetCoordinates(rContactData, rContactData.MasterGeometry, true);
    const Matrix u1 = GetVariableMatrix(rContactData, DISPLACEMENT, 0, rContactData.SlaveGeometry);
    const Matrix u2 = GetVariableMatrix(rContactData, DISPLACEMENT, 0, rContactData.MasterGeometry);
    const Matrix u1previous = GetVariableMatrix(rContactData, DISPLACEMENT, 1, rContactData.SlaveGeometry);
    const Matrix u2previous = GetVariableMatrix(rContactData, DISPLACEMENT, 1, rContactData.MasterGeometry);
    const Matrix v1 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.SlaveGeometry); // TODO: Replace with the Bossak velocity
    const Matrix v2 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.MasterGeometry);

    const double clhs0 =             1.0/Dt;
    const double clhs1 =             J*N2[0]*Phi[0]*clhs0;
    const double clhs2 =             N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0);
    const double clhs3 =             Dt*normalslave(0,0);
    const double clhs4 =             N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0);
    const double clhs5 =             clhs2*tan1slave(0,0) + clhs3*clhs4;
    const double clhs6 =             N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1);
    const double clhs7 =             N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1);
    const double clhs8 =             clhs3*clhs7 + clhs6*tan1slave(0,0);
    const double clhs9 =             J*N2[1]*Phi[0]*clhs0;
    const double clhs10 =             J*N1[0]*Phi[0]*clhs0;
    const double clhs11 =             J*N1[1]*Phi[0]*clhs0;
    const double clhs12 =             Dt*normalslave(0,1);
    const double clhs13 =             clhs12*clhs4 + clhs2*tan1slave(0,1);
    const double clhs14 =             clhs12*clhs7 + clhs6*tan1slave(0,1);
    const double clhs15 =             J*N2[0]*Phi[1]*clhs0;
    const double clhs16 =             Dt*normalslave(1,0);
    const double clhs17 =             clhs16*clhs4 + clhs2*tan1slave(1,0);
    const double clhs18 =             clhs16*clhs7 + clhs6*tan1slave(1,0);
    const double clhs19 =             J*N2[1]*Phi[1]*clhs0;
    const double clhs20 =             J*N1[0]*Phi[1]*clhs0;
    const double clhs21 =             J*N1[1]*Phi[1]*clhs0;
    const double clhs22 =             Dt*normalslave(1,1);
    const double clhs23 =             clhs2*tan1slave(1,1) + clhs22*clhs4;
    const double clhs24 =             clhs22*clhs7 + clhs6*tan1slave(1,1);
    
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
    lhs(8,0)=-clhs1*clhs5;
    lhs(8,1)=-clhs1*clhs8;
    lhs(8,2)=-clhs5*clhs9;
    lhs(8,3)=-clhs8*clhs9;
    lhs(8,4)=clhs10*clhs5;
    lhs(8,5)=clhs10*clhs8;
    lhs(8,6)=clhs11*clhs5;
    lhs(8,7)=clhs11*clhs8;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs1*clhs13;
    lhs(9,1)=-clhs1*clhs14;
    lhs(9,2)=-clhs13*clhs9;
    lhs(9,3)=-clhs14*clhs9;
    lhs(9,4)=clhs10*clhs13;
    lhs(9,5)=clhs10*clhs14;
    lhs(9,6)=clhs11*clhs13;
    lhs(9,7)=clhs11*clhs14;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs15*clhs17;
    lhs(10,1)=-clhs15*clhs18;
    lhs(10,2)=-clhs17*clhs19;
    lhs(10,3)=-clhs18*clhs19;
    lhs(10,4)=clhs17*clhs20;
    lhs(10,5)=clhs18*clhs20;
    lhs(10,6)=clhs17*clhs21;
    lhs(10,7)=clhs18*clhs21;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs15*clhs23;
    lhs(11,1)=-clhs15*clhs24;
    lhs(11,2)=-clhs19*clhs23;
    lhs(11,3)=-clhs19*clhs24;
    lhs(11,4)=clhs20*clhs23;
    lhs(11,5)=clhs20*clhs24;
    lhs(11,6)=clhs21*clhs23;
    lhs(11,7)=clhs21*clhs24;
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
        const double J, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;

    const Matrix lm           = rContactData.LagrangeMultipliers;

    const double crhs0 =             J*N2[0];
    const double crhs1 =             Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double crhs2 =             Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double crhs3 =             J*N2[1];
    const double crhs4 =             J*N1[0];
    const double crhs5 =             J*N1[1];

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
    
    static inline array_1d<double,12> ComputeGaussPointRHSContactContribution(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double J, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const Matrix normalmaster = rContactData.NormalsMaster;
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const double Dt           = rContactData.Dt;
    
    const Matrix X1 = GetCoordinates(rContactData, rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData, rContactData.MasterGeometry, false);
    const Matrix x1 = GetCoordinates(rContactData, rContactData.SlaveGeometry, true);
    const Matrix x2 = GetCoordinates(rContactData, rContactData.MasterGeometry, true);
    const Matrix u1 = GetVariableMatrix(rContactData, DISPLACEMENT, 0, rContactData.SlaveGeometry);
    const Matrix u2 = GetVariableMatrix(rContactData, DISPLACEMENT, 0, rContactData.MasterGeometry);
    const Matrix u1previous = GetVariableMatrix(rContactData, DISPLACEMENT, 1, rContactData.SlaveGeometry);
    const Matrix u2previous = GetVariableMatrix(rContactData, DISPLACEMENT, 1, rContactData.MasterGeometry);
    const Matrix v1 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.SlaveGeometry); // TODO: Replace with the Bossak velocity
    const Matrix v2 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.MasterGeometry);

    const double crhs0 =             1.0/Dt;
    const double crhs1 =             J*Phi[0]*crhs0;
    const double crhs2 =             Dt*((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N1[0]*(X1(0,0) + u1(0,0)) + N1[1]*(X1(1,0) + u1(1,0)) - N2[0]*(X2(0,0) + u2(0,0)) - N2[1]*(X2(1,0) + u2(1,0))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,1) + u1(1,1)) - N2[0]*(X2(0,1) + u2(0,1)) - N2[1]*(X2(1,1) + u2(1,1))));
    const double crhs3 =             (N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0))*(N1[0]*(u1(0,0) - u1previous(0,0)) + N1[1]*(u1(1,0) - u1previous(1,0)) - N2[0]*(u2(0,0) - u2previous(0,0)) - N2[1]*(u2(1,0) - u2previous(1,0))) + (N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1))*(N1[0]*(u1(0,1) - u1previous(0,1)) + N1[1]*(u1(1,1) - u1previous(1,1)) - N2[0]*(u2(0,1) - u2previous(0,1)) - N2[1]*(u2(1,1) - u2previous(1,1)));
    const double crhs4 =             J*Phi[1]*crhs0;
    
    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs1*(crhs2*normalslave(0,0) + crhs3*tan1slave(0,0));
    rhs[9]=-crhs1*(crhs2*normalslave(0,1) + crhs3*tan1slave(0,1));
    rhs[10]=-crhs4*(crhs2*normalslave(1,0) + crhs3*tan1slave(1,0));
    rhs[11]=-crhs4*(crhs2*normalslave(1,1) + crhs3*tan1slave(1,1));

    return rhs;
}

private:
    
    static inline Matrix GetCoordinates(
        const ContactData& rContactData,
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
        ContactData const & rContactData,
        const Variable<array_1d<double,3> >& rVarName,
        unsigned int step,
        const GeometryType& nodes
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