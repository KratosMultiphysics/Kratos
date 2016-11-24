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
#include <boost/math/special_functions/sign.hpp>

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
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointActiveLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster   = rContactData.NormalsMaster;
    const Vector normalmasterg  = prod(trans(normalmaster), N2);
    const Matrix normalslave    = rContactData.NormalsSlave;
    const Matrix tan1slave      = rContactData.Tangent1Slave;
    const Matrix lm             = rContactData.LagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
//     const double epsilon_normal = rContactData.epsilon_normal;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

    const double Dnormalmasterg1u211 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u210 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u201 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u200 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u211 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u210 =     N2[0]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u201 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u200 =     N2[0]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave11u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave11u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave11u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave11u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave10u111 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave10u110 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave10u101 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave10u100 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave01u111 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave01u110 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave01u101 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave01u100 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))*(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave00u111 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave00u110 =     0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave00u101 =     (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
//     const double Dtan1slave00u100 =     -0.5/std::sqrt(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2)) + (-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0))*(-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::pow(std::pow(-0.5*X1(0,0) + 0.5*X1(1,0) - 0.5*u1(0,0) + 0.5*u1(1,0), 2) + std::pow(-0.5*X1(0,1) + 0.5*X1(1,1) - 0.5*u1(0,1) + 0.5*u1(1,1), 2), 3.0L/2.0L);
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
    const double DN21u211 =     -((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
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
    const double clhs4 =     clhs2*clhs3;
    const double clhs5 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs6 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs7 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs8 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs9 =     clhs8*N2[0] + detJ*DN20u100; // CLHS8*N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs10 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs11 =     clhs10*N2[0] + detJ*DN20u101; // CLHS10*N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs12 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs13 =     clhs12*N2[0] + detJ*DN20u110; // CLHS12*N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs14 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs15 =     clhs14*N2[0] + detJ*DN20u111; // CLHS14*N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs16 =     Phi[0]*clhs2;
    const double clhs17 =     clhs0*clhs16;
    const double clhs18 =     Phi[1]*clhs2;
    const double clhs19 =     clhs0*clhs18;
    const double clhs20 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double clhs21 =     clhs2*clhs20;
    const double clhs22 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs23 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs24 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs25 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs26 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs27 =     clhs8*N2[1] + detJ*DN21u100; // CLHS8*N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs28 =     clhs10*N2[1] + detJ*DN21u101; // CLHS10*N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs29 =     clhs12*N2[1] + detJ*DN21u110; // CLHS12*N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs30 =     clhs14*N2[1] + detJ*DN21u111; // CLHS14*N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)) + DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs31 =     clhs16*clhs22;
    const double clhs32 =     clhs18*clhs22;
    const double clhs33 =     N1[0]*clhs3;
    const double clhs34 =     N1[0]*clhs2;
    const double clhs35 =     -Phi[0]*clhs34;
    const double clhs36 =     -Phi[1]*clhs34;
    const double clhs37 =     N1[0]*clhs20;
    const double clhs38 =     N1[1]*clhs3;
    const double clhs39 =     N1[1]*clhs2;
    const double clhs40 =     -Phi[0]*clhs39;
    const double clhs41 =     -Phi[1]*clhs39;
    const double clhs42 =     N1[1]*clhs20;
    const double clhs43 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs44 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs45 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs46 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs47 =     Phi[0]*clhs2*(GPnormal[0]*clhs45 + GPnormal[1]*clhs46);
    const double clhs48 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs49 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs50 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs51 =     clhs2*clhs43;
    const double clhs52 =     clhs2*clhs45;
    const double clhs53 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs54 =     clhs43*clhs45;
    const double clhs55 =     clhs51*Dnormalslave00u100 + clhs52*clhs53 + clhs54*clhs8; // CLHS51*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS52*CLHS53 + CLHS54*CLHS8
    const double clhs56 =     clhs2*clhs46;
    const double clhs57 =     clhs43*clhs46;
    const double clhs58 =     clhs51*Dnormalslave01u100 + clhs53*clhs56 + clhs57*clhs8; // CLHS51*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS53*CLHS56 + CLHS57*CLHS8
    const double clhs59 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs60 =     clhs10*clhs54 + clhs51*Dnormalslave00u101 + clhs52*clhs59; // CLHS10*CLHS54 + CLHS51*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS52*CLHS59
    const double clhs61 =     clhs10*clhs57 + clhs51*Dnormalslave01u101 + clhs56*clhs59; // CLHS10*CLHS57 + CLHS51*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS56*CLHS59
    const double clhs62 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs63 =     clhs12*clhs54 + clhs51*Dnormalslave00u110 + clhs52*clhs62; // CLHS12*CLHS54 + CLHS51*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS52*CLHS62
    const double clhs64 =     clhs12*clhs57 + clhs51*Dnormalslave01u110 + clhs56*clhs62; // CLHS12*CLHS57 + CLHS51*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS56*CLHS62
    const double clhs65 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs66 =     clhs14*clhs54 + clhs51*Dnormalslave00u111 + clhs52*clhs65; // CLHS14*CLHS54 + CLHS51*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS52*CLHS65
    const double clhs67 =     clhs14*clhs57 + clhs51*Dnormalslave01u111 + clhs56*clhs65; // CLHS14*CLHS57 + CLHS51*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS56*CLHS65
    const double clhs68 =     Phi[0]*clhs2*(GPtangent1[0]*clhs45 + GPtangent1[1]*clhs46);
    const double clhs69 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs70 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs71 =     Phi[1]*clhs2*(GPnormal[0]*clhs69 + GPnormal[1]*clhs70);
    const double clhs72 =     clhs2*clhs69;
    const double clhs73 =     clhs43*clhs69;
    const double clhs74 =     clhs51*Dnormalslave10u100 + clhs53*clhs72 + clhs73*clhs8; // CLHS51*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS53*CLHS72 + CLHS73*CLHS8
    const double clhs75 =     clhs2*clhs70;
    const double clhs76 =     clhs43*clhs70;
    const double clhs77 =     clhs51*Dnormalslave11u100 + clhs53*clhs75 + clhs76*clhs8; // CLHS51*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS53*CLHS75 + CLHS76*CLHS8
    const double clhs78 =     clhs10*clhs73 + clhs51*Dnormalslave10u101 + clhs59*clhs72; // CLHS10*CLHS73 + CLHS51*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS59*CLHS72
    const double clhs79 =     clhs10*clhs76 + clhs51*Dnormalslave11u101 + clhs59*clhs75; // CLHS10*CLHS76 + CLHS51*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS59*CLHS75
    const double clhs80 =     clhs12*clhs73 + clhs51*Dnormalslave10u110 + clhs62*clhs72; // CLHS12*CLHS73 + CLHS51*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS62*CLHS72
    const double clhs81 =     clhs12*clhs76 + clhs51*Dnormalslave11u110 + clhs62*clhs75; // CLHS12*CLHS76 + CLHS51*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS62*CLHS75
    const double clhs82 =     clhs14*clhs73 + clhs51*Dnormalslave10u111 + clhs65*clhs72; // CLHS14*CLHS73 + CLHS51*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS65*CLHS72
    const double clhs83 =     clhs14*clhs76 + clhs51*Dnormalslave11u111 + clhs65*clhs75; // CLHS14*CLHS76 + CLHS51*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS65*CLHS75
    const double clhs84 =     Phi[1]*clhs2*(GPtangent1[0]*clhs69 + GPtangent1[1]*clhs70);

    lhs(0,0)=clhs1*clhs4;
    lhs(0,1)=clhs4*clhs5;
    lhs(0,2)=clhs4*clhs6;
    lhs(0,3)=clhs4*clhs7;
    lhs(0,4)=clhs3*clhs9;
    lhs(0,5)=clhs11*clhs3;
    lhs(0,6)=clhs13*clhs3;
    lhs(0,7)=clhs15*clhs3;
    lhs(0,8)=clhs17;
    lhs(0,9)=0;
    lhs(0,10)=clhs19;
    lhs(0,11)=0;
    lhs(1,0)=clhs1*clhs21;
    lhs(1,1)=clhs21*clhs5;
    lhs(1,2)=clhs21*clhs6;
    lhs(1,3)=clhs21*clhs7;
    lhs(1,4)=clhs20*clhs9;
    lhs(1,5)=clhs11*clhs20;
    lhs(1,6)=clhs13*clhs20;
    lhs(1,7)=clhs15*clhs20;
    lhs(1,8)=0;
    lhs(1,9)=clhs17;
    lhs(1,10)=0;
    lhs(1,11)=clhs19;
    lhs(2,0)=clhs23*clhs4;
    lhs(2,1)=clhs24*clhs4;
    lhs(2,2)=clhs25*clhs4;
    lhs(2,3)=clhs26*clhs4;
    lhs(2,4)=clhs27*clhs3;
    lhs(2,5)=clhs28*clhs3;
    lhs(2,6)=clhs29*clhs3;
    lhs(2,7)=clhs3*clhs30;
    lhs(2,8)=clhs31;
    lhs(2,9)=0;
    lhs(2,10)=clhs32;
    lhs(2,11)=0;
    lhs(3,0)=clhs21*clhs23;
    lhs(3,1)=clhs21*clhs24;
    lhs(3,2)=clhs21*clhs25;
    lhs(3,3)=clhs21*clhs26;
    lhs(3,4)=clhs20*clhs27;
    lhs(3,5)=clhs20*clhs28;
    lhs(3,6)=clhs20*clhs29;
    lhs(3,7)=clhs20*clhs30;
    lhs(3,8)=0;
    lhs(3,9)=clhs31;
    lhs(3,10)=0;
    lhs(3,11)=clhs32;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=-clhs33*clhs8;
    lhs(4,5)=-clhs10*clhs33;
    lhs(4,6)=-clhs12*clhs33;
    lhs(4,7)=-clhs14*clhs33;
    lhs(4,8)=clhs35;
    lhs(4,9)=0;
    lhs(4,10)=clhs36;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=-clhs37*clhs8;
    lhs(5,5)=-clhs10*clhs37;
    lhs(5,6)=-clhs12*clhs37;
    lhs(5,7)=-clhs14*clhs37;
    lhs(5,8)=0;
    lhs(5,9)=clhs35;
    lhs(5,10)=0;
    lhs(5,11)=clhs36;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=-clhs38*clhs8;
    lhs(6,5)=-clhs10*clhs38;
    lhs(6,6)=-clhs12*clhs38;
    lhs(6,7)=-clhs14*clhs38;
    lhs(6,8)=clhs40;
    lhs(6,9)=0;
    lhs(6,10)=clhs41;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=-clhs42*clhs8;
    lhs(7,5)=-clhs10*clhs42;
    lhs(7,6)=-clhs12*clhs42;
    lhs(7,7)=-clhs14*clhs42;
    lhs(7,8)=0;
    lhs(7,9)=clhs40;
    lhs(7,10)=0;
    lhs(7,11)=clhs41;
    lhs(8,0)=-clhs44*clhs47;
    lhs(8,1)=-clhs47*clhs48;
    lhs(8,2)=-clhs47*clhs49;
    lhs(8,3)=-clhs47*clhs50;
    lhs(8,4)=-Phi[0]*(GPnormal[0]*clhs55 + GPnormal[1]*clhs58);
    lhs(8,5)=-Phi[0]*(GPnormal[0]*clhs60 + GPnormal[1]*clhs61);
    lhs(8,6)=-Phi[0]*(GPnormal[0]*clhs63 + GPnormal[1]*clhs64);
    lhs(8,7)=-Phi[0]*(GPnormal[0]*clhs66 + GPnormal[1]*clhs67);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs44*clhs68;
    lhs(9,1)=-clhs48*clhs68;
    lhs(9,2)=-clhs49*clhs68;
    lhs(9,3)=-clhs50*clhs68;
    lhs(9,4)=-Phi[0]*(GPtangent1[0]*clhs55 + GPtangent1[1]*clhs58);
    lhs(9,5)=-Phi[0]*(GPtangent1[0]*clhs60 + GPtangent1[1]*clhs61);
    lhs(9,6)=-Phi[0]*(GPtangent1[0]*clhs63 + GPtangent1[1]*clhs64);
    lhs(9,7)=-Phi[0]*(GPtangent1[0]*clhs66 + GPtangent1[1]*clhs67);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs44*clhs71;
    lhs(10,1)=-clhs48*clhs71;
    lhs(10,2)=-clhs49*clhs71;
    lhs(10,3)=-clhs50*clhs71;
    lhs(10,4)=-Phi[1]*(GPnormal[0]*clhs74 + GPnormal[1]*clhs77);
    lhs(10,5)=-Phi[1]*(GPnormal[0]*clhs78 + GPnormal[1]*clhs79);
    lhs(10,6)=-Phi[1]*(GPnormal[0]*clhs80 + GPnormal[1]*clhs81);
    lhs(10,7)=-Phi[1]*(GPnormal[0]*clhs82 + GPnormal[1]*clhs83);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs44*clhs84;
    lhs(11,1)=-clhs48*clhs84;
    lhs(11,2)=-clhs49*clhs84;
    lhs(11,3)=-clhs50*clhs84;
    lhs(11,4)=-Phi[1]*(GPtangent1[0]*clhs74 + GPtangent1[1]*clhs77);
    lhs(11,5)=-Phi[1]*(GPtangent1[0]*clhs78 + GPtangent1[1]*clhs79);
    lhs(11,6)=-Phi[1]*(GPtangent1[0]*clhs80 + GPtangent1[1]*clhs81);
    lhs(11,7)=-Phi[1]*(GPtangent1[0]*clhs82 + GPtangent1[1]*clhs83);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointStickLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Vector normalmasterg   = prod(trans(normalmaster), N2);
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

    const double Dnormalmasterg1u211 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u210 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u201 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u200 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u211 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u210 =     N2[0]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u201 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u200 =     N2[0]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double DN21u211 =     -((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     Phi[0]*clhs0*clhs1*(GPnormal[0]*clhs2 + GPnormal[1]*clhs3);
    const double clhs5 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs6 =     N1[0]*clhs3 + N1[1]*clhs5;
    const double clhs7 =     Dt*v2(0,1);
    const double clhs8 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs9 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs10 =     Dt*v2(1,1);
    const double clhs11 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs12 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs13 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs14 =     N1[0]*clhs2 + N1[1]*clhs13;
    const double clhs15 =     Dt*v2(0,0);
    const double clhs16 =     Dt*v2(1,0);
    const double clhs17 =     clhs14*(clhs12*clhs16 + clhs15*clhs9 + clhs8) + clhs6*(clhs10*clhs12 + clhs7*clhs9);
    const double clhs18 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs19 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     clhs14*(clhs15*clhs18 + clhs16*clhs19) + clhs6*(clhs10*clhs19 + clhs18*clhs7 + clhs8);
    const double clhs21 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs22 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs23 =     clhs14*(clhs11 + clhs15*clhs21 + clhs16*clhs22) + clhs6*(clhs10*clhs22 + clhs21*clhs7);
    const double clhs24 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs25 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs26 =     clhs14*(clhs15*clhs24 + clhs16*clhs25) + clhs6*(clhs10*clhs25 + clhs11 + clhs24*clhs7);
    const double clhs27 =     Phi[0]*clhs0;
    const double clhs28 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs29 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs11*clhs16 - clhs15*clhs8;
    const double clhs30 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs10*clhs11 - clhs7*clhs8;
    const double clhs31 =     clhs14*clhs29 + clhs30*clhs6;
    const double clhs32 =     clhs1*clhs31;
    const double clhs33 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs34 =     clhs2*clhs31;
    const double clhs35 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs36 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs37 =     -N1[0];
    const double clhs38 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs39 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs40 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs41 =     clhs1*(clhs14*(clhs15*clhs35 + clhs16*clhs36 + clhs37) - clhs29*(N1[0]*clhs28 + N1[1]*clhs38) - clhs30*(N1[0]*clhs39 + N1[1]*clhs40) + clhs6*(clhs10*clhs36 + clhs35*clhs7));
    const double clhs42 =     -clhs2*clhs41 + clhs28*clhs32 + clhs33*clhs34;
    const double clhs43 =     clhs3*clhs31;
    const double clhs44 =     -clhs3*clhs41 + clhs32*clhs39 + clhs33*clhs43;
    const double clhs45 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs46 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs47 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs48 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs49 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs50 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     clhs1*(clhs14*(clhs15*clhs47 + clhs16*clhs48) - clhs29*(N1[0]*clhs45 + N1[1]*clhs49) - clhs30*(N1[0]*clhs50 + N1[1]*clhs51) + clhs6*(clhs10*clhs48 + clhs37 + clhs47*clhs7));
    const double clhs53 =     -clhs2*clhs52 + clhs32*clhs45 + clhs34*clhs46;
    const double clhs54 =     -clhs3*clhs52 + clhs32*clhs50 + clhs43*clhs46;
    const double clhs55 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs56 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs57 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs58 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs59 =     -N1[1];
    const double clhs60 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs61 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     clhs1*(clhs14*(clhs15*clhs57 + clhs16*clhs58 + clhs59) - clhs29*(N1[0]*clhs55 + N1[1]*clhs60) - clhs30*(N1[0]*clhs61 + N1[1]*clhs62) + clhs6*(clhs10*clhs58 + clhs57*clhs7));
    const double clhs64 =     -clhs2*clhs63 + clhs32*clhs55 + clhs34*clhs56;
    const double clhs65 =     -clhs3*clhs63 + clhs32*clhs61 + clhs43*clhs56;
    const double clhs66 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs67 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs68 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs69 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs70 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs71 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs72 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs73 =     clhs1*(clhs14*(clhs15*clhs68 + clhs16*clhs69) - clhs29*(N1[0]*clhs66 + N1[1]*clhs70) - clhs30*(N1[0]*clhs71 + N1[1]*clhs72) + clhs6*(clhs10*clhs69 + clhs59 + clhs68*clhs7));
    const double clhs74 =     -clhs2*clhs73 + clhs32*clhs66 + clhs34*clhs67;
    const double clhs75 =     -clhs3*clhs73 + clhs32*clhs71 + clhs43*clhs67;
    const double clhs76 =     Phi[0]*clhs0*clhs1*(GPtangent1[0]*clhs2 + GPtangent1[1]*clhs3);
    const double clhs77 =     Phi[1]*clhs0*clhs1*(GPnormal[0]*clhs13 + GPnormal[1]*clhs5);
    const double clhs78 =     Phi[1]*clhs0;
    const double clhs79 =     clhs13*clhs31;
    const double clhs80 =     -clhs13*clhs41 + clhs32*clhs38 + clhs33*clhs79;
    const double clhs81 =     clhs31*clhs5;
    const double clhs82 =     clhs32*clhs40 + clhs33*clhs81 - clhs41*clhs5;
    const double clhs83 =     -clhs13*clhs52 + clhs32*clhs49 + clhs46*clhs79;
    const double clhs84 =     clhs32*clhs51 + clhs46*clhs81 - clhs5*clhs52;
    const double clhs85 =     -clhs13*clhs63 + clhs32*clhs60 + clhs56*clhs79;
    const double clhs86 =     clhs32*clhs62 - clhs5*clhs63 + clhs56*clhs81;
    const double clhs87 =     -clhs13*clhs73 + clhs32*clhs70 + clhs67*clhs79;
    const double clhs88 =     clhs32*clhs72 - clhs5*clhs73 + clhs67*clhs81;
    const double clhs89 =     Phi[1]*clhs0*clhs1*(GPtangent1[0]*clhs13 + GPtangent1[1]*clhs5);

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
    lhs(8,0)=-clhs17*clhs4;
    lhs(8,1)=-clhs20*clhs4;
    lhs(8,2)=-clhs23*clhs4;
    lhs(8,3)=-clhs26*clhs4;
    lhs(8,4)=clhs27*(GPnormal[0]*clhs42 + GPnormal[1]*clhs44);
    lhs(8,5)=clhs27*(GPnormal[0]*clhs53 + GPnormal[1]*clhs54);
    lhs(8,6)=clhs27*(GPnormal[0]*clhs64 + GPnormal[1]*clhs65);
    lhs(8,7)=clhs27*(GPnormal[0]*clhs74 + GPnormal[1]*clhs75);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs17*clhs76;
    lhs(9,1)=-clhs20*clhs76;
    lhs(9,2)=-clhs23*clhs76;
    lhs(9,3)=-clhs26*clhs76;
    lhs(9,4)=clhs27*(GPtangent1[0]*clhs42 + GPtangent1[1]*clhs44);
    lhs(9,5)=clhs27*(GPtangent1[0]*clhs53 + GPtangent1[1]*clhs54);
    lhs(9,6)=clhs27*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs65);
    lhs(9,7)=clhs27*(GPtangent1[0]*clhs74 + GPtangent1[1]*clhs75);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs17*clhs77;
    lhs(10,1)=-clhs20*clhs77;
    lhs(10,2)=-clhs23*clhs77;
    lhs(10,3)=-clhs26*clhs77;
    lhs(10,4)=clhs78*(GPnormal[0]*clhs80 + GPnormal[1]*clhs82);
    lhs(10,5)=clhs78*(GPnormal[0]*clhs83 + GPnormal[1]*clhs84);
    lhs(10,6)=clhs78*(GPnormal[0]*clhs85 + GPnormal[1]*clhs86);
    lhs(10,7)=clhs78*(GPnormal[0]*clhs87 + GPnormal[1]*clhs88);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs17*clhs89;
    lhs(11,1)=-clhs20*clhs89;
    lhs(11,2)=-clhs23*clhs89;
    lhs(11,3)=-clhs26*clhs89;
    lhs(11,4)=clhs78*(GPtangent1[0]*clhs80 + GPtangent1[1]*clhs82);
    lhs(11,5)=clhs78*(GPtangent1[0]*clhs83 + GPtangent1[1]*clhs84);
    lhs(11,6)=clhs78*(GPtangent1[0]*clhs85 + GPtangent1[1]*clhs86);
    lhs(11,7)=clhs78*(GPtangent1[0]*clhs87 + GPtangent1[1]*clhs88);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointSlipLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Vector normalmasterg   = prod(trans(normalmaster), N2);
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

    const double Dnormalmasterg1u211 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u210 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u201 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u200 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u211 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u210 =     N2[0]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u201 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u200 =     N2[0]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double DN21u211 =     -((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u211*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u211*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u211*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u211*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u210*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u210*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalmasterg0u210*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u210*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u201*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u201*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u201*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u201*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[1])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (-Dnormalmasterg0u200*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - Dnormalmasterg1u200*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(Dnormalmasterg0u200*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + Dnormalmasterg1u200*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)) + normalmasterg[0])/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1]) - normalmasterg[1]*(Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[1]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1]) - normalmasterg[1]*(Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1]) - normalmasterg[1]*(Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) - N1[0]*normalmasterg[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-normalmasterg[0]*(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1]) - normalmasterg[1]*(Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1]))*(normalmasterg[0]*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + normalmasterg[1]*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow(normalmasterg[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + normalmasterg[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     Phi[0]*clhs0*clhs1*(GPnormal[0]*clhs2 + GPnormal[1]*clhs3);
    const double clhs5 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs6 =     N1[0]*clhs3 + N1[1]*clhs5;
    const double clhs7 =     Dt*v2(0,1);
    const double clhs8 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs9 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs10 =     Dt*v2(1,1);
    const double clhs11 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs12 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs13 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs14 =     N1[0]*clhs2 + N1[1]*clhs13;
    const double clhs15 =     Dt*v2(0,0);
    const double clhs16 =     Dt*v2(1,0);
    const double clhs17 =     clhs14*(clhs12*clhs16 + clhs15*clhs9 + clhs8) + clhs6*(clhs10*clhs12 + clhs7*clhs9);
    const double clhs18 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs19 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     clhs14*(clhs15*clhs18 + clhs16*clhs19) + clhs6*(clhs10*clhs19 + clhs18*clhs7 + clhs8);
    const double clhs21 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs22 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs23 =     clhs14*(clhs11 + clhs15*clhs21 + clhs16*clhs22) + clhs6*(clhs10*clhs22 + clhs21*clhs7);
    const double clhs24 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs25 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs26 =     clhs14*(clhs15*clhs24 + clhs16*clhs25) + clhs6*(clhs10*clhs25 + clhs11 + clhs24*clhs7);
    const double clhs27 =     Phi[0]*clhs0;
    const double clhs28 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs29 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs11*clhs16 - clhs15*clhs8;
    const double clhs30 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs10*clhs11 - clhs7*clhs8;
    const double clhs31 =     clhs14*clhs29 + clhs30*clhs6;
    const double clhs32 =     clhs1*clhs31;
    const double clhs33 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs34 =     clhs2*clhs31;
    const double clhs35 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs36 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs37 =     -N1[0];
    const double clhs38 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs39 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs40 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs41 =     clhs1*(clhs14*(clhs15*clhs35 + clhs16*clhs36 + clhs37) - clhs29*(N1[0]*clhs28 + N1[1]*clhs38) - clhs30*(N1[0]*clhs39 + N1[1]*clhs40) + clhs6*(clhs10*clhs36 + clhs35*clhs7));
    const double clhs42 =     -clhs2*clhs41 + clhs28*clhs32 + clhs33*clhs34;
    const double clhs43 =     clhs3*clhs31;
    const double clhs44 =     -clhs3*clhs41 + clhs32*clhs39 + clhs33*clhs43;
    const double clhs45 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs46 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs47 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs48 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs49 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs50 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     clhs1*(clhs14*(clhs15*clhs47 + clhs16*clhs48) - clhs29*(N1[0]*clhs45 + N1[1]*clhs49) - clhs30*(N1[0]*clhs50 + N1[1]*clhs51) + clhs6*(clhs10*clhs48 + clhs37 + clhs47*clhs7));
    const double clhs53 =     -clhs2*clhs52 + clhs32*clhs45 + clhs34*clhs46;
    const double clhs54 =     -clhs3*clhs52 + clhs32*clhs50 + clhs43*clhs46;
    const double clhs55 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs56 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs57 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs58 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs59 =     -N1[1];
    const double clhs60 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs61 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     clhs1*(clhs14*(clhs15*clhs57 + clhs16*clhs58 + clhs59) - clhs29*(N1[0]*clhs55 + N1[1]*clhs60) - clhs30*(N1[0]*clhs61 + N1[1]*clhs62) + clhs6*(clhs10*clhs58 + clhs57*clhs7));
    const double clhs64 =     -clhs2*clhs63 + clhs32*clhs55 + clhs34*clhs56;
    const double clhs65 =     -clhs3*clhs63 + clhs32*clhs61 + clhs43*clhs56;
    const double clhs66 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs67 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs68 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs69 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs70 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs71 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs72 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs73 =     clhs1*(clhs14*(clhs15*clhs68 + clhs16*clhs69) - clhs29*(N1[0]*clhs66 + N1[1]*clhs70) - clhs30*(N1[0]*clhs71 + N1[1]*clhs72) + clhs6*(clhs10*clhs69 + clhs59 + clhs68*clhs7));
    const double clhs74 =     -clhs2*clhs73 + clhs32*clhs66 + clhs34*clhs67;
    const double clhs75 =     -clhs3*clhs73 + clhs32*clhs71 + clhs43*clhs67;
    const double clhs76 =     Phi[0]*clhs0*clhs1*(GPtangent1[0]*clhs2 + GPtangent1[1]*clhs3);
    const double clhs77 =     Phi[1]*clhs0*clhs1*(GPnormal[0]*clhs13 + GPnormal[1]*clhs5);
    const double clhs78 =     Phi[1]*clhs0;
    const double clhs79 =     clhs13*clhs31;
    const double clhs80 =     -clhs13*clhs41 + clhs32*clhs38 + clhs33*clhs79;
    const double clhs81 =     clhs31*clhs5;
    const double clhs82 =     clhs32*clhs40 + clhs33*clhs81 - clhs41*clhs5;
    const double clhs83 =     -clhs13*clhs52 + clhs32*clhs49 + clhs46*clhs79;
    const double clhs84 =     clhs32*clhs51 + clhs46*clhs81 - clhs5*clhs52;
    const double clhs85 =     -clhs13*clhs63 + clhs32*clhs60 + clhs56*clhs79;
    const double clhs86 =     clhs32*clhs62 - clhs5*clhs63 + clhs56*clhs81;
    const double clhs87 =     -clhs13*clhs73 + clhs32*clhs70 + clhs67*clhs79;
    const double clhs88 =     clhs32*clhs72 - clhs5*clhs73 + clhs67*clhs81;
    const double clhs89 =     Phi[1]*clhs0*clhs1*(GPtangent1[0]*clhs13 + GPtangent1[1]*clhs5);

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
    lhs(8,0)=-clhs17*clhs4;
    lhs(8,1)=-clhs20*clhs4;
    lhs(8,2)=-clhs23*clhs4;
    lhs(8,3)=-clhs26*clhs4;
    lhs(8,4)=clhs27*(GPnormal[0]*clhs42 + GPnormal[1]*clhs44);
    lhs(8,5)=clhs27*(GPnormal[0]*clhs53 + GPnormal[1]*clhs54);
    lhs(8,6)=clhs27*(GPnormal[0]*clhs64 + GPnormal[1]*clhs65);
    lhs(8,7)=clhs27*(GPnormal[0]*clhs74 + GPnormal[1]*clhs75);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs17*clhs76;
    lhs(9,1)=-clhs20*clhs76;
    lhs(9,2)=-clhs23*clhs76;
    lhs(9,3)=-clhs26*clhs76;
    lhs(9,4)=clhs27*(GPtangent1[0]*clhs42 + GPtangent1[1]*clhs44);
    lhs(9,5)=clhs27*(GPtangent1[0]*clhs53 + GPtangent1[1]*clhs54);
    lhs(9,6)=clhs27*(GPtangent1[0]*clhs64 + GPtangent1[1]*clhs65);
    lhs(9,7)=clhs27*(GPtangent1[0]*clhs74 + GPtangent1[1]*clhs75);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs17*clhs77;
    lhs(10,1)=-clhs20*clhs77;
    lhs(10,2)=-clhs23*clhs77;
    lhs(10,3)=-clhs26*clhs77;
    lhs(10,4)=clhs78*(GPnormal[0]*clhs80 + GPnormal[1]*clhs82);
    lhs(10,5)=clhs78*(GPnormal[0]*clhs83 + GPnormal[1]*clhs84);
    lhs(10,6)=clhs78*(GPnormal[0]*clhs85 + GPnormal[1]*clhs86);
    lhs(10,7)=clhs78*(GPnormal[0]*clhs87 + GPnormal[1]*clhs88);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs17*clhs89;
    lhs(11,1)=-clhs20*clhs89;
    lhs(11,2)=-clhs23*clhs89;
    lhs(11,3)=-clhs26*clhs89;
    lhs(11,4)=clhs78*(GPtangent1[0]*clhs80 + GPtangent1[1]*clhs82);
    lhs(11,5)=clhs78*(GPtangent1[0]*clhs83 + GPtangent1[1]*clhs84);
    lhs(11,6)=clhs78*(GPtangent1[0]*clhs85 + GPtangent1[1]*clhs86);
    lhs(11,7)=clhs78*(GPtangent1[0]*clhs87 + GPtangent1[1]*clhs88);
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,12,12> ComputeGaussPointInactiveLHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
//         const Matrix DPhi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,12,12> lhs;
    
//     const Matrix normalmaster   = rContactData.NormalsMaster;
//     const Vector normalmasterg  = prod(trans(normalmaster), N2);
//     const Matrix normalslave    = rContactData.NormalsSlave;
//     const Matrix tan1slave      = rContactData.Tangent1Slave;
//     const Matrix lm             = rContactData.LagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
//     const double epsilon_normal = rContactData.epsilon_normal;
//     
//     const Vector GPnormal     = prod(trans(normalslave), N1);
//     const Vector GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const Matrix X1 = rContactData.X1;
//     const Matrix X2 = rContactData.X2;
//     const Matrix u1 = rContactData.u1;
//     const Matrix u2 = rContactData.u2;
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
// 
//substitute_derivatives_variables_inactive 
//substitute_inactive_lhs
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointActiveRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave    = rContactData.NormalsSlave;
    const Matrix tan1slave      = rContactData.Tangent1Slave;
    const Matrix lm             = rContactData.LagrangeMultipliers;
//     const double Dt             = rContactData.Dt;
//     const double epsilon_normal = rContactData.epsilon_normal;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double crhs3 =     crhs1*crhs2;
    const double crhs4 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double crhs5 =     crhs1*crhs4;
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N1[0]*crhs1;
    const double crhs8 =     N1[1]*crhs1;
    const double crhs9 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs11 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs12 =     Phi[0]*crhs1*crhs11;
    const double crhs13 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs14 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     Phi[1]*crhs1*crhs11;

    rhs[0]=-crhs0*crhs3;
    rhs[1]=-crhs0*crhs5;
    rhs[2]=-crhs3*crhs6;
    rhs[3]=-crhs5*crhs6;
    rhs[4]=crhs2*crhs7;
    rhs[5]=crhs4*crhs7;
    rhs[6]=crhs2*crhs8;
    rhs[7]=crhs4*crhs8;
    rhs[8]=crhs12*(GPnormal[0]*crhs9 + GPnormal[1]*crhs10);
    rhs[9]=crhs12*(GPtangent1[0]*crhs9 + GPtangent1[1]*crhs10);
    rhs[10]=crhs15*(GPnormal[0]*crhs13 + GPnormal[1]*crhs14);
    rhs[11]=crhs15*(GPtangent1[0]*crhs13 + GPtangent1[1]*crhs14);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointStickRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
//     const Matrix lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs8 =     (N1[0]*crhs0 + N1[1]*crhs4)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs5*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0))) + (N1[0]*crhs1 + N1[1]*crhs7)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs5*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1)));
    const double crhs9 =     Phi[0]*crhs2*crhs3*crhs8;
    const double crhs10 =     Phi[1]*crhs2*crhs3*crhs8;

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    rhs[9]=-crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    rhs[10]=-crhs10*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    rhs[11]=-crhs10*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointSlipRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const double mu, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,12> rhs;
    
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
//     const Matrix lm              = rContactData.LagrangeMultipliers;
    const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
//     const double sign_tangpress = boost::math::sign(augmented_tangent_lm);
    
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;
    
    const double crhs0 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     1.0/Dt;
    const double crhs3 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs5 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs8 =     (N1[0]*crhs0 + N1[1]*crhs4)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs5*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0))) + (N1[0]*crhs1 + N1[1]*crhs7)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs5*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1)));
    const double crhs9 =     Phi[0]*crhs2*crhs3*crhs8;
    const double crhs10 =     Phi[1]*crhs2*crhs3*crhs8;

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=-crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    rhs[9]=-crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    rhs[10]=-crhs10*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    rhs[11]=-crhs10*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,12> ComputeGaussPointInactiveRHS(
        const Vector N1, 
        const Vector N2, 
        const Vector Phi, 
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    array_1d<double,12> rhs;
    
//     const Matrix normalslave     = rContactData.NormalsSlave;
//     const Matrix tan1slave       = rContactData.Tangent1Slave;
//     const Matrix lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon_normal  = rContactData.epsilon_normal;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
//     
//     const Vector GPnormal     = prod(trans(normalslave), N1);
//     const Vector GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
    
//substitute_inactive_rhs
    
    return rhs;
}

private:
};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */