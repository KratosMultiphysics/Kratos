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

#if !defined(KRATOS_CONTACT_ACTIVE2D2N2N)
#define KRATOS_CONTACT_ACTIVE2D2N2N

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
        
class ContactActive2D2N2N
{
public:
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const ContactData& rContactData
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster  = rContactData.NormalsMaster;
    const Vector normalmasterg = prod(trans(normalmaster), N2);
    const Matrix normalslave   = rContactData.NormalsSlave;
    const Matrix tan1slave     = rContactData.Tangent1Slave;
    const Matrix lm            = rContactData.LagrangeMultipliers;
    const Vector Gaps         = rContactData.Gaps;
    const double Dt            = rContactData.Dt;
    const double epsilon     = rContactData.epsilon;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

//subs_debug
    
    const double Dtan1slave11u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave11u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave11u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave11u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave10u111 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave10u110 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave10u101 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave10u100 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave01u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave01u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave01u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave01u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave00u111 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave00u110 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave00u101 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dtan1slave00u100 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave11u111 =     -(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave11u110 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) - (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave11u101 =     -(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave11u100 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) - (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave10u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave10u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave10u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave10u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave01u111 =     -(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave01u110 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) - (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave01u101 =     -(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave01u100 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) - (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave00u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave00u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave00u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double Dnormalslave00u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
    const double DN21u211 =     1.0/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     1.0/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -(normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     -1/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     -1/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     (normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     (normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double Dgapgu211 =     (DN20u211*(X2(0,0) + u2(0,0)) + DN21u211*(X2(1,0) + u2(1,0)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(DN20u211*(X2(0,1) + u2(0,1)) + DN21u211*(X2(1,1) + u2(1,1)) + N2[1]);
    const double Dgapgu210 =     (DN20u210*(X2(0,1) + u2(0,1)) + DN21u210*(X2(1,1) + u2(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(DN20u210*(X2(0,0) + u2(0,0)) + DN21u210*(X2(1,0) + u2(1,0)) + N2[1]);
    const double Dgapgu201 =     (DN20u201*(X2(0,0) + u2(0,0)) + DN21u201*(X2(1,0) + u2(1,0)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(DN20u201*(X2(0,1) + u2(0,1)) + DN21u201*(X2(1,1) + u2(1,1)) + N2[0]);
    const double Dgapgu200 =     (DN20u200*(X2(0,1) + u2(0,1)) + DN21u200*(X2(1,1) + u2(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(DN20u200*(X2(0,0) + u2(0,0)) + DN21u200*(X2(1,0) + u2(1,0)) + N2[0]);
    const double Dgapgu111 =     (DN20u111*(X2(0,0) + u2(0,0)) + DN21u111*(X2(1,0) + u2(1,0)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(DN20u111*(X2(0,1) + u2(0,1)) + DN21u111*(X2(1,1) + u2(1,1)) - N1[1]);
    const double Dgapgu110 =     (DN20u110*(X2(0,1) + u2(0,1)) + DN21u110*(X2(1,1) + u2(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(DN20u110*(X2(0,0) + u2(0,0)) + DN21u110*(X2(1,0) + u2(1,0)) - N1[1]);
    const double Dgapgu101 =     (DN20u101*(X2(0,0) + u2(0,0)) + DN21u101*(X2(1,0) + u2(1,0)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(DN20u101*(X2(0,1) + u2(0,1)) + DN21u101*(X2(1,1) + u2(1,1)) - N1[0]);
    const double Dgapgu100 =     (DN20u100*(X2(0,1) + u2(0,1)) + DN21u100*(X2(1,1) + u2(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(DN20u100*(X2(0,0) + u2(0,0)) + DN21u100*(X2(1,0) + u2(1,0)) - N1[0]);
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));

    const double clhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs1 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double clhs4 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs6 =     N1[0]*clhs4 + N1[1]*clhs5;
    const double clhs7 =     clhs6*epsilon;
    const double clhs8 =     clhs2*(clhs3 - clhs7);
    const double clhs9 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs10 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs11 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs12 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs13 =     clhs12*clhs3;
    const double clhs14 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs15 =     clhs2*clhs3;
    const double clhs16 =     clhs12*clhs7;
    const double clhs17 =     clhs2*clhs6*epsilon;
    const double clhs18 =     clhs2*epsilon;
    const double clhs19 =     Dnormalslave00u100; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs20 =     Dnormalslave10u100; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs21 =     clhs18*(N1[0]*clhs19 + N1[1]*clhs20);
    const double clhs22 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs23 =     clhs22*clhs3;
    const double clhs24 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs25 =     clhs22*clhs7;
    const double clhs26 =     Dnormalslave00u101; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs27 =     Dnormalslave10u101; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs28 =     clhs18*(N1[0]*clhs26 + N1[1]*clhs27);
    const double clhs29 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs30 =     clhs29*clhs3;
    const double clhs31 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs32 =     clhs29*clhs7;
    const double clhs33 =     Dnormalslave00u110; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs34 =     Dnormalslave10u110; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs35 =     clhs18*(N1[0]*clhs33 + N1[1]*clhs34);
    const double clhs36 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs37 =     clhs3*clhs36;
    const double clhs38 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs39 =     clhs36*clhs7;
    const double clhs40 =     Dnormalslave00u111; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs41 =     Dnormalslave10u111; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs42 =     clhs18*(N1[0]*clhs40 + N1[1]*clhs41);
    const double clhs43 =     Phi[0]*clhs2;
    const double clhs44 =     clhs0*clhs43;
    const double clhs45 =     Phi[1]*clhs2;
    const double clhs46 =     clhs0*clhs45;
    const double clhs47 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double clhs48 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs49 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs50 =     N1[0]*clhs48 + N1[1]*clhs49;
    const double clhs51 =     clhs50*epsilon;
    const double clhs52 =     clhs2*(clhs47 - clhs51);
    const double clhs53 =     clhs12*clhs47;
    const double clhs54 =     clhs2*clhs47;
    const double clhs55 =     clhs12*clhs51;
    const double clhs56 =     clhs2*clhs50*epsilon;
    const double clhs57 =     Dnormalslave01u100; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs58 =     Dnormalslave11u100; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs59 =     clhs18*(N1[0]*clhs57 + N1[1]*clhs58);
    const double clhs60 =     clhs22*clhs47;
    const double clhs61 =     clhs22*clhs51;
    const double clhs62 =     Dnormalslave01u101; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs63 =     Dnormalslave11u101; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs64 =     clhs18*(N1[0]*clhs62 + N1[1]*clhs63);
    const double clhs65 =     clhs29*clhs47;
    const double clhs66 =     clhs29*clhs51;
    const double clhs67 =     Dnormalslave01u110; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     Dnormalslave11u110; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     clhs18*(N1[0]*clhs67 + N1[1]*clhs68);
    const double clhs70 =     clhs36*clhs47;
    const double clhs71 =     clhs36*clhs51;
    const double clhs72 =     Dnormalslave01u111; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs73 =     Dnormalslave11u111; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs74 =     clhs18*(N1[0]*clhs72 + N1[1]*clhs73);
    const double clhs75 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs76 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs77 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs78 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs79 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs80 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs81 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs82 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs83 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs84 =     clhs43*clhs75;
    const double clhs85 =     clhs45*clhs75;
    const double clhs86 =     -clhs13 + clhs16 + clhs21;
    const double clhs87 =     -clhs23 + clhs25 + clhs28;
    const double clhs88 =     -clhs30 + clhs32 + clhs35;
    const double clhs89 =     -clhs37 + clhs39 + clhs42;
    const double clhs90 =     N1[0]*clhs2;
    const double clhs91 =     -Phi[0]*clhs90;
    const double clhs92 =     -Phi[1]*clhs90;
    const double clhs93 =     -clhs53 + clhs55 + clhs59;
    const double clhs94 =     -clhs60 + clhs61 + clhs64;
    const double clhs95 =     -clhs65 + clhs66 + clhs69;
    const double clhs96 =     -clhs70 + clhs71 + clhs74;
    const double clhs97 =     N1[1]*clhs2;
    const double clhs98 =     -Phi[0]*clhs97;
    const double clhs99 =     -Phi[1]*clhs97;
    const double clhs100 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs101 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs102 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs103 =     1.0/Dt;
    const double clhs104 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs105 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs106 =     N1[0]*clhs104 + N1[1]*clhs105;
    const double clhs107 =     Dt*v2(0,1);
    const double clhs108 =     Dt*v2(1,1);
    const double clhs109 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs110 =     N1[0]*clhs102 + N1[1]*clhs109;
    const double clhs111 =     Dt*v2(0,0);
    const double clhs112 =     Dt*v2(1,0);
    const double clhs113 =     clhs103*(clhs106*(clhs1*clhs107 + clhs108*clhs76) + clhs110*(clhs0 + clhs1*clhs111 + clhs112*clhs76));
    const double clhs114 =     clhs101*clhs4 + clhs102*clhs113;
    const double clhs115 =     clhs101*clhs48 + clhs104*clhs113;
    const double clhs116 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs117 =     clhs103*(clhs106*(clhs0 + clhs107*clhs9 + clhs108*clhs77) + clhs110*(clhs111*clhs9 + clhs112*clhs77));
    const double clhs118 =     clhs102*clhs117 + clhs116*clhs4;
    const double clhs119 =     clhs104*clhs117 + clhs116*clhs48;
    const double clhs120 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs121 =     clhs103*(clhs106*(clhs10*clhs107 + clhs108*clhs78) + clhs110*(clhs10*clhs111 + clhs112*clhs78 + clhs75));
    const double clhs122 =     clhs102*clhs121 + clhs120*clhs4;
    const double clhs123 =     clhs104*clhs121 + clhs120*clhs48;
    const double clhs124 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs125 =     clhs103*(clhs106*(clhs107*clhs11 + clhs108*clhs79 + clhs75) + clhs110*(clhs11*clhs111 + clhs112*clhs79));
    const double clhs126 =     clhs102*clhs125 + clhs124*clhs4;
    const double clhs127 =     clhs104*clhs125 + clhs124*clhs48;
    const double clhs128 =     clhs100*clhs2;
    const double clhs129 =     clhs2*clhs4;
    const double clhs130 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs131 =     clhs100*clhs4;
    const double clhs132 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs133 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs0*clhs111 - clhs112*clhs75;
    const double clhs134 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs0*clhs107 - clhs108*clhs75;
    const double clhs135 =     clhs106*clhs134 + clhs110*clhs133;
    const double clhs136 =     clhs103*clhs135*clhs2;
    const double clhs137 =     clhs102*clhs103*clhs135;
    const double clhs138 =     -N1[0];
    const double clhs139 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs140 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs141 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs142 =     clhs103*clhs2*(clhs106*(clhs107*clhs14 + clhs108*clhs80) + clhs110*(clhs111*clhs14 + clhs112*clhs80 + clhs138) - clhs133*(N1[0]*clhs132 + N1[1]*clhs139) - clhs134*(N1[0]*clhs140 + N1[1]*clhs141));
    const double clhs143 =     -clhs102*clhs142 - clhs12*clhs131 + clhs12*clhs137 - clhs128*clhs19 - clhs129*clhs130 + clhs132*clhs136;
    const double clhs144 =     clhs2*clhs48;
    const double clhs145 =     clhs100*clhs48;
    const double clhs146 =     clhs103*clhs104*clhs135;
    const double clhs147 =     -clhs104*clhs142 - clhs12*clhs145 + clhs12*clhs146 - clhs128*clhs57 - clhs130*clhs144 + clhs136*clhs140;
    const double clhs148 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs149 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs150 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs151 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs152 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs153 =     clhs103*clhs2*(clhs106*(clhs107*clhs24 + clhs108*clhs81 + clhs138) + clhs110*(clhs111*clhs24 + clhs112*clhs81) - clhs133*(N1[0]*clhs149 + N1[1]*clhs150) - clhs134*(N1[0]*clhs151 + N1[1]*clhs152));
    const double clhs154 =     -clhs102*clhs153 - clhs128*clhs26 - clhs129*clhs148 - clhs131*clhs22 + clhs136*clhs149 + clhs137*clhs22;
    const double clhs155 =     -clhs104*clhs153 - clhs128*clhs62 + clhs136*clhs151 - clhs144*clhs148 - clhs145*clhs22 + clhs146*clhs22;
    const double clhs156 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs157 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs158 =     -N1[1];
    const double clhs159 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs160 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs161 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs162 =     clhs103*clhs2*(clhs106*(clhs107*clhs31 + clhs108*clhs82) + clhs110*(clhs111*clhs31 + clhs112*clhs82 + clhs158) - clhs133*(N1[0]*clhs157 + N1[1]*clhs159) - clhs134*(N1[0]*clhs160 + N1[1]*clhs161));
    const double clhs163 =     -clhs102*clhs162 - clhs128*clhs33 - clhs129*clhs156 - clhs131*clhs29 + clhs136*clhs157 + clhs137*clhs29;
    const double clhs164 =     -clhs104*clhs162 - clhs128*clhs67 + clhs136*clhs160 - clhs144*clhs156 - clhs145*clhs29 + clhs146*clhs29;
    const double clhs165 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs166 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs167 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs168 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs169 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs170 =     clhs103*clhs2*(clhs106*(clhs107*clhs38 + clhs108*clhs83 + clhs158) + clhs110*(clhs111*clhs38 + clhs112*clhs83) - clhs133*(N1[0]*clhs166 + N1[1]*clhs167) - clhs134*(N1[0]*clhs168 + N1[1]*clhs169));
    const double clhs171 =     -clhs102*clhs170 - clhs128*clhs40 - clhs129*clhs165 - clhs131*clhs36 + clhs136*clhs166 + clhs137*clhs36;
    const double clhs172 =     -clhs104*clhs170 - clhs128*clhs72 + clhs136*clhs168 - clhs144*clhs165 - clhs145*clhs36 + clhs146*clhs36;
    const double clhs173 =     clhs101*clhs5 + clhs109*clhs113;
    const double clhs174 =     clhs101*clhs49 + clhs105*clhs113;
    const double clhs175 =     clhs109*clhs117 + clhs116*clhs5;
    const double clhs176 =     clhs105*clhs117 + clhs116*clhs49;
    const double clhs177 =     clhs109*clhs121 + clhs120*clhs5;
    const double clhs178 =     clhs105*clhs121 + clhs120*clhs49;
    const double clhs179 =     clhs109*clhs125 + clhs124*clhs5;
    const double clhs180 =     clhs105*clhs125 + clhs124*clhs49;
    const double clhs181 =     clhs2*clhs5;
    const double clhs182 =     clhs100*clhs5;
    const double clhs183 =     clhs103*clhs109*clhs135;
    const double clhs184 =     -clhs109*clhs142 - clhs12*clhs182 + clhs12*clhs183 - clhs128*clhs20 - clhs130*clhs181 + clhs136*clhs139;
    const double clhs185 =     clhs2*clhs49;
    const double clhs186 =     clhs100*clhs49;
    const double clhs187 =     clhs103*clhs105*clhs135;
    const double clhs188 =     -clhs105*clhs142 - clhs12*clhs186 + clhs12*clhs187 - clhs128*clhs58 - clhs130*clhs185 + clhs136*clhs141;
    const double clhs189 =     -clhs109*clhs153 - clhs128*clhs27 + clhs136*clhs150 - clhs148*clhs181 - clhs182*clhs22 + clhs183*clhs22;
    const double clhs190 =     -clhs105*clhs153 - clhs128*clhs63 + clhs136*clhs152 - clhs148*clhs185 - clhs186*clhs22 + clhs187*clhs22;
    const double clhs191 =     -clhs109*clhs162 - clhs128*clhs34 + clhs136*clhs159 - clhs156*clhs181 - clhs182*clhs29 + clhs183*clhs29;
    const double clhs192 =     -clhs105*clhs162 - clhs128*clhs68 + clhs136*clhs161 - clhs156*clhs185 - clhs186*clhs29 + clhs187*clhs29;
    const double clhs193 =     -clhs109*clhs170 - clhs128*clhs41 + clhs136*clhs167 - clhs165*clhs181 - clhs182*clhs36 + clhs183*clhs36;
    const double clhs194 =     -clhs105*clhs170 - clhs128*clhs73 + clhs136*clhs169 - clhs165*clhs185 - clhs186*clhs36 + clhs187*clhs36;

    lhs(0,0)=clhs1*clhs8;
    lhs(0,1)=clhs8*clhs9;
    lhs(0,2)=clhs10*clhs8;
    lhs(0,3)=clhs11*clhs8;
    lhs(0,4)=clhs0*clhs13 - clhs0*clhs16 - clhs0*clhs21 + clhs14*clhs15 - clhs14*clhs17;
    lhs(0,5)=clhs0*clhs23 - clhs0*clhs25 - clhs0*clhs28 + clhs15*clhs24 - clhs17*clhs24;
    lhs(0,6)=clhs0*clhs30 - clhs0*clhs32 - clhs0*clhs35 + clhs15*clhs31 - clhs17*clhs31;
    lhs(0,7)=clhs0*clhs37 - clhs0*clhs39 - clhs0*clhs42 + clhs15*clhs38 - clhs17*clhs38;
    lhs(0,8)=clhs44;
    lhs(0,9)=0;
    lhs(0,10)=clhs46;
    lhs(0,11)=0;
    lhs(1,0)=clhs1*clhs52;
    lhs(1,1)=clhs52*clhs9;
    lhs(1,2)=clhs10*clhs52;
    lhs(1,3)=clhs11*clhs52;
    lhs(1,4)=clhs0*clhs53 - clhs0*clhs55 - clhs0*clhs59 + clhs14*clhs54 - clhs14*clhs56;
    lhs(1,5)=clhs0*clhs60 - clhs0*clhs61 - clhs0*clhs64 + clhs24*clhs54 - clhs24*clhs56;
    lhs(1,6)=clhs0*clhs65 - clhs0*clhs66 - clhs0*clhs69 + clhs31*clhs54 - clhs31*clhs56;
    lhs(1,7)=clhs0*clhs70 - clhs0*clhs71 - clhs0*clhs74 + clhs38*clhs54 - clhs38*clhs56;
    lhs(1,8)=0;
    lhs(1,9)=clhs44;
    lhs(1,10)=0;
    lhs(1,11)=clhs46;
    lhs(2,0)=clhs76*clhs8;
    lhs(2,1)=clhs77*clhs8;
    lhs(2,2)=clhs78*clhs8;
    lhs(2,3)=clhs79*clhs8;
    lhs(2,4)=clhs13*clhs75 + clhs15*clhs80 - clhs16*clhs75 - clhs17*clhs80 - clhs21*clhs75;
    lhs(2,5)=clhs15*clhs81 - clhs17*clhs81 + clhs23*clhs75 - clhs25*clhs75 - clhs28*clhs75;
    lhs(2,6)=clhs15*clhs82 - clhs17*clhs82 + clhs30*clhs75 - clhs32*clhs75 - clhs35*clhs75;
    lhs(2,7)=clhs15*clhs83 - clhs17*clhs83 + clhs37*clhs75 - clhs39*clhs75 - clhs42*clhs75;
    lhs(2,8)=clhs84;
    lhs(2,9)=0;
    lhs(2,10)=clhs85;
    lhs(2,11)=0;
    lhs(3,0)=clhs52*clhs76;
    lhs(3,1)=clhs52*clhs77;
    lhs(3,2)=clhs52*clhs78;
    lhs(3,3)=clhs52*clhs79;
    lhs(3,4)=clhs53*clhs75 + clhs54*clhs80 - clhs55*clhs75 - clhs56*clhs80 - clhs59*clhs75;
    lhs(3,5)=clhs54*clhs81 - clhs56*clhs81 + clhs60*clhs75 - clhs61*clhs75 - clhs64*clhs75;
    lhs(3,6)=clhs54*clhs82 - clhs56*clhs82 + clhs65*clhs75 - clhs66*clhs75 - clhs69*clhs75;
    lhs(3,7)=clhs54*clhs83 - clhs56*clhs83 + clhs70*clhs75 - clhs71*clhs75 - clhs74*clhs75;
    lhs(3,8)=0;
    lhs(3,9)=clhs84;
    lhs(3,10)=0;
    lhs(3,11)=clhs85;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=N1[0]*clhs86;
    lhs(4,5)=N1[0]*clhs87;
    lhs(4,6)=N1[0]*clhs88;
    lhs(4,7)=N1[0]*clhs89;
    lhs(4,8)=clhs91;
    lhs(4,9)=0;
    lhs(4,10)=clhs92;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=N1[0]*clhs93;
    lhs(5,5)=N1[0]*clhs94;
    lhs(5,6)=N1[0]*clhs95;
    lhs(5,7)=N1[0]*clhs96;
    lhs(5,8)=0;
    lhs(5,9)=clhs91;
    lhs(5,10)=0;
    lhs(5,11)=clhs92;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=N1[1]*clhs86;
    lhs(6,5)=N1[1]*clhs87;
    lhs(6,6)=N1[1]*clhs88;
    lhs(6,7)=N1[1]*clhs89;
    lhs(6,8)=clhs98;
    lhs(6,9)=0;
    lhs(6,10)=clhs99;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=N1[1]*clhs93;
    lhs(7,5)=N1[1]*clhs94;
    lhs(7,6)=N1[1]*clhs95;
    lhs(7,7)=N1[1]*clhs96;
    lhs(7,8)=0;
    lhs(7,9)=clhs98;
    lhs(7,10)=0;
    lhs(7,11)=clhs99;
    lhs(8,0)=-clhs43*(GPnormal[0]*clhs114 + GPnormal[1]*clhs115);
    lhs(8,1)=-clhs43*(GPnormal[0]*clhs118 + GPnormal[1]*clhs119);
    lhs(8,2)=-clhs43*(GPnormal[0]*clhs122 + GPnormal[1]*clhs123);
    lhs(8,3)=-clhs43*(GPnormal[0]*clhs126 + GPnormal[1]*clhs127);
    lhs(8,4)=Phi[0]*(GPnormal[0]*clhs143 + GPnormal[1]*clhs147);
    lhs(8,5)=Phi[0]*(GPnormal[0]*clhs154 + GPnormal[1]*clhs155);
    lhs(8,6)=Phi[0]*(GPnormal[0]*clhs163 + GPnormal[1]*clhs164);
    lhs(8,7)=Phi[0]*(GPnormal[0]*clhs171 + GPnormal[1]*clhs172);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs43*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs115);
    lhs(9,1)=-clhs43*(GPtangent1[0]*clhs118 + GPtangent1[1]*clhs119);
    lhs(9,2)=-clhs43*(GPtangent1[0]*clhs122 + GPtangent1[1]*clhs123);
    lhs(9,3)=-clhs43*(GPtangent1[0]*clhs126 + GPtangent1[1]*clhs127);
    lhs(9,4)=Phi[0]*(GPtangent1[0]*clhs143 + GPtangent1[1]*clhs147);
    lhs(9,5)=Phi[0]*(GPtangent1[0]*clhs154 + GPtangent1[1]*clhs155);
    lhs(9,6)=Phi[0]*(GPtangent1[0]*clhs163 + GPtangent1[1]*clhs164);
    lhs(9,7)=Phi[0]*(GPtangent1[0]*clhs171 + GPtangent1[1]*clhs172);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs45*(GPnormal[0]*clhs173 + GPnormal[1]*clhs174);
    lhs(10,1)=-clhs45*(GPnormal[0]*clhs175 + GPnormal[1]*clhs176);
    lhs(10,2)=-clhs45*(GPnormal[0]*clhs177 + GPnormal[1]*clhs178);
    lhs(10,3)=-clhs45*(GPnormal[0]*clhs179 + GPnormal[1]*clhs180);
    lhs(10,4)=Phi[1]*(GPnormal[0]*clhs184 + GPnormal[1]*clhs188);
    lhs(10,5)=Phi[1]*(GPnormal[0]*clhs189 + GPnormal[1]*clhs190);
    lhs(10,6)=Phi[1]*(GPnormal[0]*clhs191 + GPnormal[1]*clhs192);
    lhs(10,7)=Phi[1]*(GPnormal[0]*clhs193 + GPnormal[1]*clhs194);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs45*(GPtangent1[0]*clhs173 + GPtangent1[1]*clhs174);
    lhs(11,1)=-clhs45*(GPtangent1[0]*clhs175 + GPtangent1[1]*clhs176);
    lhs(11,2)=-clhs45*(GPtangent1[0]*clhs177 + GPtangent1[1]*clhs178);
    lhs(11,3)=-clhs45*(GPtangent1[0]*clhs179 + GPtangent1[1]*clhs180);
    lhs(11,4)=Phi[1]*(GPtangent1[0]*clhs184 + GPtangent1[1]*clhs188);
    lhs(11,5)=Phi[1]*(GPtangent1[0]*clhs189 + GPtangent1[1]*clhs190);
    lhs(11,6)=Phi[1]*(GPtangent1[0]*clhs191 + GPtangent1[1]*clhs192);
    lhs(11,7)=Phi[1]*(GPtangent1[0]*clhs193 + GPtangent1[1]*clhs194);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi,
        const double detJ, 
        const ContactData& rContactData
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave  = rContactData.NormalsSlave;
    const Matrix tan1slave    = rContactData.Tangent1Slave;
    const Vector Gaps         = rContactData.Gaps;
    const Matrix lm           = rContactData.LagrangeMultipliers;
    const double Dt           = rContactData.Dt;
    const double epsilon   = rContactData.epsilon;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]*lm(0,0);
    const double crhs3 =     Phi[1]*lm(1,0);
    const double crhs4 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs6 =     epsilon*(N1[0]*crhs4 + N1[1]*crhs5);
    const double crhs7 =     crhs1*(-crhs2 - crhs3 + crhs6);
    const double crhs8 =     Phi[0]*lm(0,1);
    const double crhs9 =     Phi[1]*lm(1,1);
    const double crhs10 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs11 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     epsilon*(N1[0]*crhs10 + N1[1]*crhs11);
    const double crhs13 =     crhs1*(crhs12 - crhs8 - crhs9);
    const double crhs14 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs15 =     N1[0]*crhs1;
    const double crhs16 =     crhs2 + crhs3 - crhs6;
    const double crhs17 =     -crhs12 + crhs8 + crhs9;
    const double crhs18 =     N1[1]*crhs1;
    const double crhs19 =     Phi[0]*crhs1;
    const double crhs20 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs21 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs22 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs23 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs24 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs25 =     ((N1[0]*crhs21 + N1[1]*crhs22)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs0*(Dt*v2(0,0)) - crhs14*(Dt*v2(1,0))) + (N1[0]*crhs23 + N1[1]*crhs24)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs0*(Dt*v2(0,1)) - crhs14*(Dt*v2(1,1))))/Dt;
    const double crhs26 =     -crhs20*crhs4 + crhs21*crhs25;
    const double crhs27 =     -crhs10*crhs20 + crhs23*crhs25;
    const double crhs28 =     Phi[1]*crhs1;
    const double crhs29 =     -crhs20*crhs5 + crhs22*crhs25;
    const double crhs30 =     -crhs11*crhs20 + crhs24*crhs25;

    rhs[0]=crhs0*crhs7;
    rhs[1]=crhs0*crhs13;
    rhs[2]=crhs14*crhs7;
    rhs[3]=crhs13*crhs14;
    rhs[4]=crhs15*crhs16;
    rhs[5]=crhs15*crhs17;
    rhs[6]=crhs16*crhs18;
    rhs[7]=crhs17*crhs18;
    rhs[8]=-crhs19*(GPnormal[0]*crhs26 + GPnormal[1]*crhs27);
    rhs[9]=-crhs19*(GPtangent1[0]*crhs26 + GPtangent1[1]*crhs27);
    rhs[10]=-crhs28*(GPnormal[0]*crhs29 + GPnormal[1]*crhs30);
    rhs[11]=-crhs28*(GPtangent1[0]*crhs29 + GPtangent1[1]*crhs30);

    
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

};// class ContactActive2D2N2N
}
#endif /* KRATOS_CONTACT_ACTIVE2D2N2N defined */