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

#if !defined(KRATOS_CONTACT_INACTIVE2D2N2N)
#define KRATOS_CONTACT_INACTIVE2D2N2N

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
        
class ContactInactive2D2N2N
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

    const double clhs0 =     1.0/epsilon;
    const double clhs1 =     Phi[0]*clhs0;
    const double clhs2 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     Dnormalslave00u100; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs4 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs6 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs8 =     Phi[0]*(clhs2*lm(0,0) + clhs5*lm(0,1)) + Phi[1]*(clhs6*lm(1,0) + clhs7*lm(1,1));
    const double clhs9 =     clhs4*clhs8;
    const double clhs10 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs11 =     clhs2*clhs8;
    const double clhs12 =     Dnormalslave01u100; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs13 =     Dnormalslave10u100; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs14 =     Dnormalslave11u100; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs15 =     clhs4*(Phi[0]*(clhs12*lm(0,1) + clhs3*lm(0,0)) + Phi[1]*(clhs13*lm(1,0) + clhs14*lm(1,1)));
    const double clhs16 =     clhs10*clhs11 + clhs15*clhs2 + clhs3*clhs9;
    const double clhs17 =     clhs5*clhs8;
    const double clhs18 =     clhs10*clhs17 + clhs12*clhs9 + clhs15*clhs5;
    const double clhs19 =     Dnormalslave00u101; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs20 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs21 =     Dnormalslave01u101; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs22 =     Dnormalslave10u101; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs23 =     Dnormalslave11u101; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs24 =     clhs4*(Phi[0]*(clhs19*lm(0,0) + clhs21*lm(0,1)) + Phi[1]*(clhs22*lm(1,0) + clhs23*lm(1,1)));
    const double clhs25 =     clhs11*clhs20 + clhs19*clhs9 + clhs2*clhs24;
    const double clhs26 =     clhs17*clhs20 + clhs21*clhs9 + clhs24*clhs5;
    const double clhs27 =     Dnormalslave00u110; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs28 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs29 =     Dnormalslave01u110; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs30 =     Dnormalslave10u110; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs31 =     Dnormalslave11u110; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs32 =     clhs4*(Phi[0]*(clhs27*lm(0,0) + clhs29*lm(0,1)) + Phi[1]*(clhs30*lm(1,0) + clhs31*lm(1,1)));
    const double clhs33 =     clhs11*clhs28 + clhs2*clhs32 + clhs27*clhs9;
    const double clhs34 =     clhs17*clhs28 + clhs29*clhs9 + clhs32*clhs5;
    const double clhs35 =     Dnormalslave00u111; // DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs36 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs37 =     Dnormalslave01u111; // DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs38 =     Dnormalslave10u111; // DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs39 =     Dnormalslave11u111; // DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs40 =     clhs4*(Phi[0]*(clhs35*lm(0,0) + clhs37*lm(0,1)) + Phi[1]*(clhs38*lm(1,0) + clhs39*lm(1,1)));
    const double clhs41 =     clhs11*clhs36 + clhs2*clhs40 + clhs35*clhs9;
    const double clhs42 =     clhs17*clhs36 + clhs37*clhs9 + clhs40*clhs5;
    const double clhs43 =     std::pow(Phi[0], 2)*clhs0*clhs4;
    const double clhs44 =     -clhs2*clhs43*clhs5;
    const double clhs45 =     Phi[0]*Phi[1]*clhs0*clhs2*clhs4;
    const double clhs46 =     -clhs45*clhs6;
    const double clhs47 =     -clhs45*clhs7;
    const double clhs48 =     Phi[0]*Phi[1]*clhs0*clhs4*clhs5;
    const double clhs49 =     -clhs48*clhs6;
    const double clhs50 =     -clhs48*clhs7;
    const double clhs51 =     Phi[1]*clhs0;
    const double clhs52 =     clhs6*clhs8;
    const double clhs53 =     clhs10*clhs52 + clhs13*clhs9 + clhs15*clhs6;
    const double clhs54 =     clhs7*clhs8;
    const double clhs55 =     clhs10*clhs54 + clhs14*clhs9 + clhs15*clhs7;
    const double clhs56 =     clhs20*clhs52 + clhs22*clhs9 + clhs24*clhs6;
    const double clhs57 =     clhs20*clhs54 + clhs23*clhs9 + clhs24*clhs7;
    const double clhs58 =     clhs28*clhs52 + clhs30*clhs9 + clhs32*clhs6;
    const double clhs59 =     clhs28*clhs54 + clhs31*clhs9 + clhs32*clhs7;
    const double clhs60 =     clhs36*clhs52 + clhs38*clhs9 + clhs40*clhs6;
    const double clhs61 =     clhs36*clhs54 + clhs39*clhs9 + clhs40*clhs7;
    const double clhs62 =     std::pow(Phi[1], 2)*clhs0*clhs4;
    const double clhs63 =     -clhs6*clhs62*clhs7;

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
    lhs(8,0)=0;
    lhs(8,1)=0;
    lhs(8,2)=0;
    lhs(8,3)=0;
    lhs(8,4)=-clhs1*(GPnormal[0]*clhs16 + GPnormal[1]*clhs18);
    lhs(8,5)=-clhs1*(GPnormal[0]*clhs25 + GPnormal[1]*clhs26);
    lhs(8,6)=-clhs1*(GPnormal[0]*clhs33 + GPnormal[1]*clhs34);
    lhs(8,7)=-clhs1*(GPnormal[0]*clhs41 + GPnormal[1]*clhs42);
    lhs(8,8)=-std::pow(clhs2, 2)*clhs43;
    lhs(8,9)=clhs44;
    lhs(8,10)=clhs46;
    lhs(8,11)=clhs47;
    lhs(9,0)=0;
    lhs(9,1)=0;
    lhs(9,2)=0;
    lhs(9,3)=0;
    lhs(9,4)=-clhs1*(GPtangent1[0]*clhs16 + GPtangent1[1]*clhs18);
    lhs(9,5)=-clhs1*(GPtangent1[0]*clhs25 + GPtangent1[1]*clhs26);
    lhs(9,6)=-clhs1*(GPtangent1[0]*clhs33 + GPtangent1[1]*clhs34);
    lhs(9,7)=-clhs1*(GPtangent1[0]*clhs41 + GPtangent1[1]*clhs42);
    lhs(9,8)=clhs44;
    lhs(9,9)=-clhs43*std::pow(clhs5, 2);
    lhs(9,10)=clhs49;
    lhs(9,11)=clhs50;
    lhs(10,0)=0;
    lhs(10,1)=0;
    lhs(10,2)=0;
    lhs(10,3)=0;
    lhs(10,4)=-clhs51*(GPnormal[0]*clhs53 + GPnormal[1]*clhs55);
    lhs(10,5)=-clhs51*(GPnormal[0]*clhs56 + GPnormal[1]*clhs57);
    lhs(10,6)=-clhs51*(GPnormal[0]*clhs58 + GPnormal[1]*clhs59);
    lhs(10,7)=-clhs51*(GPnormal[0]*clhs60 + GPnormal[1]*clhs61);
    lhs(10,8)=clhs46;
    lhs(10,9)=clhs49;
    lhs(10,10)=-std::pow(clhs6, 2)*clhs62;
    lhs(10,11)=clhs63;
    lhs(11,0)=0;
    lhs(11,1)=0;
    lhs(11,2)=0;
    lhs(11,3)=0;
    lhs(11,4)=-clhs51*(GPtangent1[0]*clhs53 + GPtangent1[1]*clhs55);
    lhs(11,5)=-clhs51*(GPtangent1[0]*clhs56 + GPtangent1[1]*clhs57);
    lhs(11,6)=-clhs51*(GPtangent1[0]*clhs58 + GPtangent1[1]*clhs59);
    lhs(11,7)=-clhs51*(GPtangent1[0]*clhs60 + GPtangent1[1]*clhs61);
    lhs(11,8)=clhs47;
    lhs(11,9)=clhs50;
    lhs(11,10)=clhs63;
    lhs(11,11)=-clhs62*std::pow(clhs7, 2);

    
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
    
    const double crhs0 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/epsilon;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs6 =     Phi[0]*(crhs0*lm(0,0) + crhs1*lm(0,1)) + Phi[1]*(crhs4*lm(1,0) + crhs5*lm(1,1));
    const double crhs7 =     Phi[0]*crhs2*crhs3*crhs6;
    const double crhs8 =     Phi[1]*crhs2*crhs3*crhs6;

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs7*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    rhs[9]=crhs7*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    rhs[10]=crhs8*(GPnormal[0]*crhs4 + GPnormal[1]*crhs5);
    rhs[11]=crhs8*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs5);

    
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

};// class ContactInactive2D2N2N
}
#endif /* KRATOS_CONTACT_INACTIVE2D2N2N defined */