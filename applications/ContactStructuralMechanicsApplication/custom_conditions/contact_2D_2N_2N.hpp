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
    
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const double Dt           = rContactData.Dt;

    const double clhs0 =             J*N2[0]*Phi[0];
    const double clhs1 =             N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0);
    const double clhs2 =             clhs1*normalslave(0,0);
    const double clhs3 =             1.0/Dt;
    const double clhs4 =             clhs3*tan1slave(0,0);
    const double clhs5 =             N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0);
    const double clhs6 =             clhs4*clhs5;
    const double clhs7 =             clhs2 - clhs6;
    const double clhs8 =             N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1);
    const double clhs9 =             clhs8*normalslave(0,0);
    const double clhs10 =             N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1);
    const double clhs11 =             clhs10*clhs4;
    const double clhs12 =             -clhs11 + clhs9;
    const double clhs13 =             J*N2[1]*Phi[0];
    const double clhs14 =             J*N1[0]*Phi[0];
    const double clhs15 =             -clhs2 + clhs6;
    const double clhs16 =             clhs11 - clhs9;
    const double clhs17 =             J*N1[1]*Phi[0];
    const double clhs18 =             clhs1*normalslave(0,1);
    const double clhs19 =             clhs3*tan1slave(0,1);
    const double clhs20 =             clhs19*clhs5;
    const double clhs21 =             clhs18 - clhs20;
    const double clhs22 =             clhs8*normalslave(0,1);
    const double clhs23 =             clhs10*clhs19;
    const double clhs24 =             clhs22 - clhs23;
    const double clhs25 =             -clhs18 + clhs20;
    const double clhs26 =             -clhs22 + clhs23;
    const double clhs27 =             J*N2[0]*Phi[1];
    const double clhs28 =             clhs1*normalslave(1,0);
    const double clhs29 =             clhs3*tan1slave(1,0);
    const double clhs30 =             clhs29*clhs5;
    const double clhs31 =             clhs28 - clhs30;
    const double clhs32 =             clhs8*normalslave(1,0);
    const double clhs33 =             clhs10*clhs29;
    const double clhs34 =             clhs32 - clhs33;
    const double clhs35 =             J*N2[1]*Phi[1];
    const double clhs36 =             J*N1[0]*Phi[1];
    const double clhs37 =             -clhs28 + clhs30;
    const double clhs38 =             -clhs32 + clhs33;
    const double clhs39 =             J*N1[1]*Phi[1];
    const double clhs40 =             clhs1*normalslave(1,1);
    const double clhs41 =             clhs3*tan1slave(1,1);
    const double clhs42 =             clhs41*clhs5;
    const double clhs43 =             clhs40 - clhs42;
    const double clhs44 =             clhs8*normalslave(1,1);
    const double clhs45 =             clhs10*clhs41;
    const double clhs46 =             clhs44 - clhs45;
    const double clhs47 =             -clhs40 + clhs42;
    const double clhs48 =             -clhs44 + clhs45;

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
    lhs(8,0)=clhs0*clhs7;
    lhs(8,1)=clhs0*clhs12;
    lhs(8,2)=clhs13*clhs7;
    lhs(8,3)=clhs12*clhs13;
    lhs(8,4)=clhs14*clhs15;
    lhs(8,5)=clhs14*clhs16;
    lhs(8,6)=clhs15*clhs17;
    lhs(8,7)=clhs16*clhs17;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs0*clhs21;
    lhs(9,1)=clhs0*clhs24;
    lhs(9,2)=clhs13*clhs21;
    lhs(9,3)=clhs13*clhs24;
    lhs(9,4)=clhs14*clhs25;
    lhs(9,5)=clhs14*clhs26;
    lhs(9,6)=clhs17*clhs25;
    lhs(9,7)=clhs17*clhs26;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs27*clhs31;
    lhs(10,1)=clhs27*clhs34;
    lhs(10,2)=clhs31*clhs35;
    lhs(10,3)=clhs34*clhs35;
    lhs(10,4)=clhs36*clhs37;
    lhs(10,5)=clhs36*clhs38;
    lhs(10,6)=clhs37*clhs39;
    lhs(10,7)=clhs38*clhs39;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs27*clhs43;
    lhs(11,1)=clhs27*clhs46;
    lhs(11,2)=clhs35*clhs43;
    lhs(11,3)=clhs35*clhs46;
    lhs(11,4)=clhs36*clhs47;
    lhs(11,5)=clhs36*clhs48;
    lhs(11,6)=clhs39*clhs47;
    lhs(11,7)=clhs39*clhs48;
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
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Vector Gaps         = rContactData.Gaps;
    const double Dt           = rContactData.Dt;
    
    const Matrix v1 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.SlaveGeometry); 
    const Matrix v2 = GetVariableMatrix(rContactData, VELOCITY, 0, rContactData.MasterGeometry);

    const double crhs0 =             J*Phi[0];
    const double crhs1 =             - inner_prod(Gaps, N1);
    const double crhs2 =             ((N1[0]*tan1slave(0,0) + N1[1]*tan1slave(1,0))*(N1[0]*Dt*v1(0,0) + N1[1]*Dt*v1(1,0) - N2[0]*Dt*v2(0,0) - N2[1]*Dt*v2(1,0)) + (N1[0]*tan1slave(0,1) + N1[1]*tan1slave(1,1))*(N1[0]*Dt*v1(0,1) + N1[1]*Dt*v1(1,1) - N2[0]*Dt*v2(0,1) - N2[1]*Dt*v2(1,1)))/Dt;
    const double crhs3 =             J*Phi[1];
    
    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs0*(crhs1*normalslave(0,0) - crhs2*tan1slave(0,0));
    rhs[9]=crhs0*(crhs1*normalslave(0,1) - crhs2*tan1slave(0,1));
    rhs[10]=crhs3*(crhs1*normalslave(1,0) - crhs2*tan1slave(1,0));
    rhs[11]=crhs3*(crhs1*normalslave(1,1) - crhs2*tan1slave(1,1));

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