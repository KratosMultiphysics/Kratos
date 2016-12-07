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

#if !defined(KRATOS_CONTACT2D2N2NDLM)
#define KRATOS_CONTACT2D2N2NDLM

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
        
class Contact2D2N2NDLM
{
public:
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointActiveLHS(
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
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Vector normalmasterg   = prod(trans(normalmaster), N2);
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
    const double epsilon         = rContactData.epsilon;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

//substitute_constants_derivatives_variables_active 
    const double Dnormalmasterg1u211 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u210 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u201 =     -N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg1u200 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[0]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - N2[1]*(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u211 =     0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u210 =     N2[0]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u201 =     -0.5*N2[0]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[0]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) - 0.5*N2[1]/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + N2[1]*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmasterg0u200 =     N2[0]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L) + N2[1]*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double clhs3 =     Phi[0]*dlm(0,0);
    const double clhs4 =     Phi[1]*dlm(1,0);
    const double clhs5 =     clhs3 + clhs4;
    const double clhs6 =     Phi[0]*lm(0,0);
    const double clhs7 =     Phi[1]*lm(1,0);
    const double clhs8 =     clhs6 + clhs7;
    const double clhs9 =     clhs5 + clhs8;
    const double clhs10 =     clhs2*clhs9;
    const double clhs11 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs12 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs13 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs14 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs15 =     clhs0*clhs5;
    const double clhs16 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs17 =     clhs2*clhs5;
    const double clhs18 =     clhs0*clhs8;
    const double clhs19 =     clhs2*clhs8;
    const double clhs20 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs21 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs22 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs23 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs24 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs25 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs26 =     Phi[0]*clhs2;
    const double clhs27 =     clhs0*clhs26;
    const double clhs28 =     Phi[1]*clhs2;
    const double clhs29 =     clhs0*clhs28;
    const double clhs30 =     Phi[0]*dlm(0,1);
    const double clhs31 =     Phi[1]*dlm(1,1);
    const double clhs32 =     clhs30 + clhs31;
    const double clhs33 =     Phi[0]*lm(0,1);
    const double clhs34 =     Phi[1]*lm(1,1);
    const double clhs35 =     clhs33 + clhs34;
    const double clhs36 =     clhs32 + clhs35;
    const double clhs37 =     clhs2*clhs36;
    const double clhs38 =     clhs0*clhs32;
    const double clhs39 =     clhs2*clhs32;
    const double clhs40 =     clhs0*clhs35;
    const double clhs41 =     clhs2*clhs35;
    const double clhs42 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs43 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs44 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs45 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs46 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs47 =     clhs42*clhs5;
    const double clhs48 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs49 =     clhs42*clhs8;
    const double clhs50 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs51 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs52 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs53 =     clhs26*clhs42;
    const double clhs54 =     clhs28*clhs42;
    const double clhs55 =     clhs32*clhs42;
    const double clhs56 =     clhs35*clhs42;
    const double clhs57 =     N1[0]*clhs9;
    const double clhs58 =     N1[0]*clhs2;
    const double clhs59 =     -Phi[0]*clhs58;
    const double clhs60 =     -Phi[1]*clhs58;
    const double clhs61 =     N1[0]*clhs36;
    const double clhs62 =     N1[1]*clhs9;
    const double clhs63 =     N1[1]*clhs2;
    const double clhs64 =     -Phi[0]*clhs63;
    const double clhs65 =     -Phi[1]*clhs63;
    const double clhs66 =     N1[1]*clhs36;
    const double clhs67 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs68 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs69 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs70 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs71 =     Phi[0]*clhs2*(GPnormal[0]*clhs69 + GPnormal[1]*clhs70);
    const double clhs72 =     -clhs68*clhs71;
    const double clhs73 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs74 =     -clhs71*clhs73;
    const double clhs75 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs76 =     -clhs71*clhs75;
    const double clhs77 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs78 =     -clhs71*clhs77;
    const double clhs79 =     epsilon*(-clhs3 - clhs4 + clhs8);
    const double clhs80 =     clhs14*clhs79;
    const double clhs81 =     clhs2*clhs67;
    const double clhs82 =     clhs2*clhs69;
    const double clhs83 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs84 =     clhs67*clhs69;
    const double clhs85 =     clhs14*clhs84 + clhs81*Dnormalslave00u100 + clhs82*clhs83; // CLHS14*CLHS84 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS82*CLHS83
    const double clhs86 =     clhs80 + clhs85;
    const double clhs87 =     epsilon*(-clhs30 - clhs31 + clhs35);
    const double clhs88 =     clhs14*clhs87;
    const double clhs89 =     clhs2*clhs70;
    const double clhs90 =     clhs67*clhs70;
    const double clhs91 =     clhs14*clhs90 + clhs81*Dnormalslave01u100 + clhs83*clhs89; // CLHS14*CLHS90 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS83*CLHS89
    const double clhs92 =     clhs88 + clhs91;
    const double clhs93 =     clhs20*clhs79;
    const double clhs94 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs95 =     clhs20*clhs84 + clhs81*Dnormalslave00u101 + clhs82*clhs94; // CLHS20*CLHS84 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS82*CLHS94
    const double clhs96 =     clhs93 + clhs95;
    const double clhs97 =     clhs20*clhs87;
    const double clhs98 =     clhs20*clhs90 + clhs81*Dnormalslave01u101 + clhs89*clhs94; // CLHS20*CLHS90 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS89*CLHS94
    const double clhs99 =     clhs97 + clhs98;
    const double clhs100 =     clhs22*clhs79;
    const double clhs101 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs102 =     clhs101*clhs82 + clhs22*clhs84 + clhs81*Dnormalslave00u110; // CLHS101*CLHS82 + CLHS22*CLHS84 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs103 =     clhs100 + clhs102;
    const double clhs104 =     clhs22*clhs87;
    const double clhs105 =     clhs101*clhs89 + clhs22*clhs90 + clhs81*Dnormalslave01u110; // CLHS101*CLHS89 + CLHS22*CLHS90 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs106 =     clhs104 + clhs105;
    const double clhs107 =     clhs24*clhs79;
    const double clhs108 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs109 =     clhs108*clhs82 + clhs24*clhs84 + clhs81*Dnormalslave00u111; // CLHS108*CLHS82 + CLHS24*CLHS84 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs110 =     clhs107 + clhs109;
    const double clhs111 =     clhs24*clhs87;
    const double clhs112 =     clhs108*clhs89 + clhs24*clhs90 + clhs81*Dnormalslave01u111; // CLHS108*CLHS89 + CLHS24*CLHS90 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs113 =     clhs111 + clhs112;
    const double clhs114 =     std::pow(Phi[0], 2);
    const double clhs115 =     GPnormal[0]*clhs2*epsilon;
    const double clhs116 =     clhs114*clhs115;
    const double clhs117 =     -clhs116;
    const double clhs118 =     GPnormal[1]*clhs2*epsilon;
    const double clhs119 =     clhs114*clhs118;
    const double clhs120 =     -clhs119;
    const double clhs121 =     Phi[0]*Phi[1]*clhs2*epsilon;
    const double clhs122 =     GPnormal[0]*clhs121;
    const double clhs123 =     -clhs122;
    const double clhs124 =     GPnormal[1]*clhs121;
    const double clhs125 =     -clhs124;
    const double clhs126 =     Phi[0]*clhs2*(GPtangent1[0]*clhs69 + GPtangent1[1]*clhs70);
    const double clhs127 =     -clhs126*clhs68;
    const double clhs128 =     -clhs126*clhs73;
    const double clhs129 =     -clhs126*clhs75;
    const double clhs130 =     -clhs126*clhs77;
    const double clhs131 =     GPtangent1[0]*clhs2*epsilon;
    const double clhs132 =     clhs114*clhs131;
    const double clhs133 =     -clhs132;
    const double clhs134 =     GPtangent1[1]*clhs2*epsilon;
    const double clhs135 =     clhs114*clhs134;
    const double clhs136 =     -clhs135;
    const double clhs137 =     GPtangent1[0]*clhs121;
    const double clhs138 =     -clhs137;
    const double clhs139 =     GPtangent1[1]*clhs121;
    const double clhs140 =     -clhs139;
    const double clhs141 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs142 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs143 =     Phi[1]*clhs2*(GPnormal[0]*clhs141 + GPnormal[1]*clhs142);
    const double clhs144 =     -clhs143*clhs68;
    const double clhs145 =     -clhs143*clhs73;
    const double clhs146 =     -clhs143*clhs75;
    const double clhs147 =     -clhs143*clhs77;
    const double clhs148 =     clhs141*clhs2;
    const double clhs149 =     clhs141*clhs67;
    const double clhs150 =     clhs14*clhs149 + clhs148*clhs83 + clhs81*Dnormalslave10u100; // CLHS14*CLHS149 + CLHS148*CLHS83 + CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs151 =     clhs150 + clhs80;
    const double clhs152 =     clhs142*clhs2;
    const double clhs153 =     clhs142*clhs67;
    const double clhs154 =     clhs14*clhs153 + clhs152*clhs83 + clhs81*Dnormalslave11u100; // CLHS14*CLHS153 + CLHS152*CLHS83 + CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs155 =     clhs154 + clhs88;
    const double clhs156 =     clhs148*clhs94 + clhs149*clhs20 + clhs81*Dnormalslave10u101; // CLHS148*CLHS94 + CLHS149*CLHS20 + CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs157 =     clhs156 + clhs93;
    const double clhs158 =     clhs152*clhs94 + clhs153*clhs20 + clhs81*Dnormalslave11u101; // CLHS152*CLHS94 + CLHS153*CLHS20 + CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs159 =     clhs158 + clhs97;
    const double clhs160 =     clhs101*clhs148 + clhs149*clhs22 + clhs81*Dnormalslave10u110; // CLHS101*CLHS148 + CLHS149*CLHS22 + CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs161 =     clhs100 + clhs160;
    const double clhs162 =     clhs101*clhs152 + clhs153*clhs22 + clhs81*Dnormalslave11u110; // CLHS101*CLHS152 + CLHS153*CLHS22 + CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs163 =     clhs104 + clhs162;
    const double clhs164 =     clhs108*clhs148 + clhs149*clhs24 + clhs81*Dnormalslave10u111; // CLHS108*CLHS148 + CLHS149*CLHS24 + CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs165 =     clhs107 + clhs164;
    const double clhs166 =     clhs108*clhs152 + clhs153*clhs24 + clhs81*Dnormalslave11u111; // CLHS108*CLHS152 + CLHS153*CLHS24 + CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs167 =     clhs111 + clhs166;
    const double clhs168 =     std::pow(Phi[1], 2);
    const double clhs169 =     clhs115*clhs168;
    const double clhs170 =     -clhs169;
    const double clhs171 =     clhs118*clhs168;
    const double clhs172 =     -clhs171;
    const double clhs173 =     Phi[1]*clhs2*(GPtangent1[0]*clhs141 + GPtangent1[1]*clhs142);
    const double clhs174 =     -clhs173*clhs68;
    const double clhs175 =     -clhs173*clhs73;
    const double clhs176 =     -clhs173*clhs75;
    const double clhs177 =     -clhs173*clhs77;
    const double clhs178 =     clhs131*clhs168;
    const double clhs179 =     -clhs178;
    const double clhs180 =     clhs134*clhs168;
    const double clhs181 =     -clhs180;
    const double clhs182 =     epsilon*(clhs5 - clhs6 - clhs7);
    const double clhs183 =     clhs14*clhs182;
    const double clhs184 =     clhs183 + clhs85;
    const double clhs185 =     epsilon*(clhs32 - clhs33 - clhs34);
    const double clhs186 =     clhs14*clhs185;
    const double clhs187 =     clhs186 + clhs91;
    const double clhs188 =     clhs182*clhs20;
    const double clhs189 =     clhs188 + clhs95;
    const double clhs190 =     clhs185*clhs20;
    const double clhs191 =     clhs190 + clhs98;
    const double clhs192 =     clhs182*clhs22;
    const double clhs193 =     clhs102 + clhs192;
    const double clhs194 =     clhs185*clhs22;
    const double clhs195 =     clhs105 + clhs194;
    const double clhs196 =     clhs182*clhs24;
    const double clhs197 =     clhs109 + clhs196;
    const double clhs198 =     clhs185*clhs24;
    const double clhs199 =     clhs112 + clhs198;
    const double clhs200 =     clhs150 + clhs183;
    const double clhs201 =     clhs154 + clhs186;
    const double clhs202 =     clhs156 + clhs188;
    const double clhs203 =     clhs158 + clhs190;
    const double clhs204 =     clhs160 + clhs192;
    const double clhs205 =     clhs162 + clhs194;
    const double clhs206 =     clhs164 + clhs196;
    const double clhs207 =     clhs166 + clhs198;

    lhs(0,0)=clhs1*clhs10;
    lhs(0,1)=clhs10*clhs11;
    lhs(0,2)=clhs10*clhs12;
    lhs(0,3)=clhs10*clhs13;
    lhs(0,4)=clhs14*clhs15 + clhs14*clhs18 + clhs16*clhs17 + clhs16*clhs19;
    lhs(0,5)=clhs15*clhs20 + clhs17*clhs21 + clhs18*clhs20 + clhs19*clhs21;
    lhs(0,6)=clhs15*clhs22 + clhs17*clhs23 + clhs18*clhs22 + clhs19*clhs23;
    lhs(0,7)=clhs15*clhs24 + clhs17*clhs25 + clhs18*clhs24 + clhs19*clhs25;
    lhs(0,8)=clhs27;
    lhs(0,9)=0;
    lhs(0,10)=clhs29;
    lhs(0,11)=0;
    lhs(0,12)=clhs27;
    lhs(0,13)=0;
    lhs(0,14)=clhs29;
    lhs(0,15)=0;
    lhs(1,0)=clhs1*clhs37;
    lhs(1,1)=clhs11*clhs37;
    lhs(1,2)=clhs12*clhs37;
    lhs(1,3)=clhs13*clhs37;
    lhs(1,4)=clhs14*clhs38 + clhs14*clhs40 + clhs16*clhs39 + clhs16*clhs41;
    lhs(1,5)=clhs20*clhs38 + clhs20*clhs40 + clhs21*clhs39 + clhs21*clhs41;
    lhs(1,6)=clhs22*clhs38 + clhs22*clhs40 + clhs23*clhs39 + clhs23*clhs41;
    lhs(1,7)=clhs24*clhs38 + clhs24*clhs40 + clhs25*clhs39 + clhs25*clhs41;
    lhs(1,8)=0;
    lhs(1,9)=clhs27;
    lhs(1,10)=0;
    lhs(1,11)=clhs29;
    lhs(1,12)=0;
    lhs(1,13)=clhs27;
    lhs(1,14)=0;
    lhs(1,15)=clhs29;
    lhs(2,0)=clhs10*clhs43;
    lhs(2,1)=clhs10*clhs44;
    lhs(2,2)=clhs10*clhs45;
    lhs(2,3)=clhs10*clhs46;
    lhs(2,4)=clhs14*clhs47 + clhs14*clhs49 + clhs17*clhs48 + clhs19*clhs48;
    lhs(2,5)=clhs17*clhs50 + clhs19*clhs50 + clhs20*clhs47 + clhs20*clhs49;
    lhs(2,6)=clhs17*clhs51 + clhs19*clhs51 + clhs22*clhs47 + clhs22*clhs49;
    lhs(2,7)=clhs17*clhs52 + clhs19*clhs52 + clhs24*clhs47 + clhs24*clhs49;
    lhs(2,8)=clhs53;
    lhs(2,9)=0;
    lhs(2,10)=clhs54;
    lhs(2,11)=0;
    lhs(2,12)=clhs53;
    lhs(2,13)=0;
    lhs(2,14)=clhs54;
    lhs(2,15)=0;
    lhs(3,0)=clhs37*clhs43;
    lhs(3,1)=clhs37*clhs44;
    lhs(3,2)=clhs37*clhs45;
    lhs(3,3)=clhs37*clhs46;
    lhs(3,4)=clhs14*clhs55 + clhs14*clhs56 + clhs39*clhs48 + clhs41*clhs48;
    lhs(3,5)=clhs20*clhs55 + clhs20*clhs56 + clhs39*clhs50 + clhs41*clhs50;
    lhs(3,6)=clhs22*clhs55 + clhs22*clhs56 + clhs39*clhs51 + clhs41*clhs51;
    lhs(3,7)=clhs24*clhs55 + clhs24*clhs56 + clhs39*clhs52 + clhs41*clhs52;
    lhs(3,8)=0;
    lhs(3,9)=clhs53;
    lhs(3,10)=0;
    lhs(3,11)=clhs54;
    lhs(3,12)=0;
    lhs(3,13)=clhs53;
    lhs(3,14)=0;
    lhs(3,15)=clhs54;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=-clhs14*clhs57;
    lhs(4,5)=-clhs20*clhs57;
    lhs(4,6)=-clhs22*clhs57;
    lhs(4,7)=-clhs24*clhs57;
    lhs(4,8)=clhs59;
    lhs(4,9)=0;
    lhs(4,10)=clhs60;
    lhs(4,11)=0;
    lhs(4,12)=clhs59;
    lhs(4,13)=0;
    lhs(4,14)=clhs60;
    lhs(4,15)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=-clhs14*clhs61;
    lhs(5,5)=-clhs20*clhs61;
    lhs(5,6)=-clhs22*clhs61;
    lhs(5,7)=-clhs24*clhs61;
    lhs(5,8)=0;
    lhs(5,9)=clhs59;
    lhs(5,10)=0;
    lhs(5,11)=clhs60;
    lhs(5,12)=0;
    lhs(5,13)=clhs59;
    lhs(5,14)=0;
    lhs(5,15)=clhs60;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=-clhs14*clhs62;
    lhs(6,5)=-clhs20*clhs62;
    lhs(6,6)=-clhs22*clhs62;
    lhs(6,7)=-clhs24*clhs62;
    lhs(6,8)=clhs64;
    lhs(6,9)=0;
    lhs(6,10)=clhs65;
    lhs(6,11)=0;
    lhs(6,12)=clhs64;
    lhs(6,13)=0;
    lhs(6,14)=clhs65;
    lhs(6,15)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=-clhs14*clhs66;
    lhs(7,5)=-clhs20*clhs66;
    lhs(7,6)=-clhs22*clhs66;
    lhs(7,7)=-clhs24*clhs66;
    lhs(7,8)=0;
    lhs(7,9)=clhs64;
    lhs(7,10)=0;
    lhs(7,11)=clhs65;
    lhs(7,12)=0;
    lhs(7,13)=clhs64;
    lhs(7,14)=0;
    lhs(7,15)=clhs65;
    lhs(8,0)=clhs72;
    lhs(8,1)=clhs74;
    lhs(8,2)=clhs76;
    lhs(8,3)=clhs78;
    lhs(8,4)=-Phi[0]*(GPnormal[0]*clhs86 + GPnormal[1]*clhs92);
    lhs(8,5)=-Phi[0]*(GPnormal[0]*clhs96 + GPnormal[1]*clhs99);
    lhs(8,6)=-Phi[0]*(GPnormal[0]*clhs103 + GPnormal[1]*clhs106);
    lhs(8,7)=-Phi[0]*(GPnormal[0]*clhs110 + GPnormal[1]*clhs113);
    lhs(8,8)=clhs117;
    lhs(8,9)=clhs120;
    lhs(8,10)=clhs123;
    lhs(8,11)=clhs125;
    lhs(8,12)=clhs116;
    lhs(8,13)=clhs119;
    lhs(8,14)=clhs122;
    lhs(8,15)=clhs124;
    lhs(9,0)=clhs127;
    lhs(9,1)=clhs128;
    lhs(9,2)=clhs129;
    lhs(9,3)=clhs130;
    lhs(9,4)=-Phi[0]*(GPtangent1[0]*clhs86 + GPtangent1[1]*clhs92);
    lhs(9,5)=-Phi[0]*(GPtangent1[0]*clhs96 + GPtangent1[1]*clhs99);
    lhs(9,6)=-Phi[0]*(GPtangent1[0]*clhs103 + GPtangent1[1]*clhs106);
    lhs(9,7)=-Phi[0]*(GPtangent1[0]*clhs110 + GPtangent1[1]*clhs113);
    lhs(9,8)=clhs133;
    lhs(9,9)=clhs136;
    lhs(9,10)=clhs138;
    lhs(9,11)=clhs140;
    lhs(9,12)=clhs132;
    lhs(9,13)=clhs135;
    lhs(9,14)=clhs137;
    lhs(9,15)=clhs139;
    lhs(10,0)=clhs144;
    lhs(10,1)=clhs145;
    lhs(10,2)=clhs146;
    lhs(10,3)=clhs147;
    lhs(10,4)=-Phi[1]*(GPnormal[0]*clhs151 + GPnormal[1]*clhs155);
    lhs(10,5)=-Phi[1]*(GPnormal[0]*clhs157 + GPnormal[1]*clhs159);
    lhs(10,6)=-Phi[1]*(GPnormal[0]*clhs161 + GPnormal[1]*clhs163);
    lhs(10,7)=-Phi[1]*(GPnormal[0]*clhs165 + GPnormal[1]*clhs167);
    lhs(10,8)=clhs123;
    lhs(10,9)=clhs125;
    lhs(10,10)=clhs170;
    lhs(10,11)=clhs172;
    lhs(10,12)=clhs122;
    lhs(10,13)=clhs124;
    lhs(10,14)=clhs169;
    lhs(10,15)=clhs171;
    lhs(11,0)=clhs174;
    lhs(11,1)=clhs175;
    lhs(11,2)=clhs176;
    lhs(11,3)=clhs177;
    lhs(11,4)=-Phi[1]*(GPtangent1[0]*clhs151 + GPtangent1[1]*clhs155);
    lhs(11,5)=-Phi[1]*(GPtangent1[0]*clhs157 + GPtangent1[1]*clhs159);
    lhs(11,6)=-Phi[1]*(GPtangent1[0]*clhs161 + GPtangent1[1]*clhs163);
    lhs(11,7)=-Phi[1]*(GPtangent1[0]*clhs165 + GPtangent1[1]*clhs167);
    lhs(11,8)=clhs138;
    lhs(11,9)=clhs140;
    lhs(11,10)=clhs179;
    lhs(11,11)=clhs181;
    lhs(11,12)=clhs137;
    lhs(11,13)=clhs139;
    lhs(11,14)=clhs178;
    lhs(11,15)=clhs180;
    lhs(12,0)=clhs72;
    lhs(12,1)=clhs74;
    lhs(12,2)=clhs76;
    lhs(12,3)=clhs78;
    lhs(12,4)=-Phi[0]*(GPnormal[0]*clhs184 + GPnormal[1]*clhs187);
    lhs(12,5)=-Phi[0]*(GPnormal[0]*clhs189 + GPnormal[1]*clhs191);
    lhs(12,6)=-Phi[0]*(GPnormal[0]*clhs193 + GPnormal[1]*clhs195);
    lhs(12,7)=-Phi[0]*(GPnormal[0]*clhs197 + GPnormal[1]*clhs199);
    lhs(12,8)=clhs116;
    lhs(12,9)=clhs119;
    lhs(12,10)=clhs122;
    lhs(12,11)=clhs124;
    lhs(12,12)=clhs117;
    lhs(12,13)=clhs120;
    lhs(12,14)=clhs123;
    lhs(12,15)=clhs125;
    lhs(13,0)=clhs127;
    lhs(13,1)=clhs128;
    lhs(13,2)=clhs129;
    lhs(13,3)=clhs130;
    lhs(13,4)=-Phi[0]*(GPtangent1[0]*clhs184 + GPtangent1[1]*clhs187);
    lhs(13,5)=-Phi[0]*(GPtangent1[0]*clhs189 + GPtangent1[1]*clhs191);
    lhs(13,6)=-Phi[0]*(GPtangent1[0]*clhs193 + GPtangent1[1]*clhs195);
    lhs(13,7)=-Phi[0]*(GPtangent1[0]*clhs197 + GPtangent1[1]*clhs199);
    lhs(13,8)=clhs132;
    lhs(13,9)=clhs135;
    lhs(13,10)=clhs137;
    lhs(13,11)=clhs139;
    lhs(13,12)=clhs133;
    lhs(13,13)=clhs136;
    lhs(13,14)=clhs138;
    lhs(13,15)=clhs140;
    lhs(14,0)=clhs144;
    lhs(14,1)=clhs145;
    lhs(14,2)=clhs146;
    lhs(14,3)=clhs147;
    lhs(14,4)=-Phi[1]*(GPnormal[0]*clhs200 + GPnormal[1]*clhs201);
    lhs(14,5)=-Phi[1]*(GPnormal[0]*clhs202 + GPnormal[1]*clhs203);
    lhs(14,6)=-Phi[1]*(GPnormal[0]*clhs204 + GPnormal[1]*clhs205);
    lhs(14,7)=-Phi[1]*(GPnormal[0]*clhs206 + GPnormal[1]*clhs207);
    lhs(14,8)=clhs122;
    lhs(14,9)=clhs124;
    lhs(14,10)=clhs169;
    lhs(14,11)=clhs171;
    lhs(14,12)=clhs123;
    lhs(14,13)=clhs125;
    lhs(14,14)=clhs170;
    lhs(14,15)=clhs172;
    lhs(15,0)=clhs174;
    lhs(15,1)=clhs175;
    lhs(15,2)=clhs176;
    lhs(15,3)=clhs177;
    lhs(15,4)=-Phi[1]*(GPtangent1[0]*clhs200 + GPtangent1[1]*clhs201);
    lhs(15,5)=-Phi[1]*(GPtangent1[0]*clhs202 + GPtangent1[1]*clhs203);
    lhs(15,6)=-Phi[1]*(GPtangent1[0]*clhs204 + GPtangent1[1]*clhs205);
    lhs(15,7)=-Phi[1]*(GPtangent1[0]*clhs206 + GPtangent1[1]*clhs207);
    lhs(15,8)=clhs137;
    lhs(15,9)=clhs139;
    lhs(15,10)=clhs178;
    lhs(15,11)=clhs180;
    lhs(15,12)=clhs138;
    lhs(15,13)=clhs140;
    lhs(15,14)=clhs179;
    lhs(15,15)=clhs181;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointStickLHS(
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
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Vector normalmasterg   = prod(trans(normalmaster), N2);
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

//substitute_constants_derivatives_variables_stick 
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
    const double clhs18 =     -clhs17*clhs4;
    const double clhs19 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     clhs14*(clhs15*clhs19 + clhs16*clhs20) + clhs6*(clhs10*clhs20 + clhs19*clhs7 + clhs8);
    const double clhs22 =     -clhs21*clhs4;
    const double clhs23 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs24 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs25 =     clhs14*(clhs11 + clhs15*clhs23 + clhs16*clhs24) + clhs6*(clhs10*clhs24 + clhs23*clhs7);
    const double clhs26 =     -clhs25*clhs4;
    const double clhs27 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs28 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs29 =     clhs14*(clhs15*clhs27 + clhs16*clhs28) + clhs6*(clhs10*clhs28 + clhs11 + clhs27*clhs7);
    const double clhs30 =     -clhs29*clhs4;
    const double clhs31 =     Phi[0]*clhs0;
    const double clhs32 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs33 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs11*clhs16 - clhs15*clhs8;
    const double clhs34 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs10*clhs11 - clhs7*clhs8;
    const double clhs35 =     clhs14*clhs33 + clhs34*clhs6;
    const double clhs36 =     clhs1*clhs35;
    const double clhs37 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs38 =     clhs2*clhs35;
    const double clhs39 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs40 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs41 =     -N1[0];
    const double clhs42 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs43 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs44 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     clhs1*(clhs14*(clhs15*clhs39 + clhs16*clhs40 + clhs41) - clhs33*(N1[0]*clhs32 + N1[1]*clhs42) - clhs34*(N1[0]*clhs43 + N1[1]*clhs44) + clhs6*(clhs10*clhs40 + clhs39*clhs7));
    const double clhs46 =     -clhs2*clhs45 + clhs32*clhs36 + clhs37*clhs38;
    const double clhs47 =     clhs3*clhs35;
    const double clhs48 =     -clhs3*clhs45 + clhs36*clhs43 + clhs37*clhs47;
    const double clhs49 =     clhs31*(GPnormal[0]*clhs46 + GPnormal[1]*clhs48);
    const double clhs50 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs53 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs54 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     clhs1*(clhs14*(clhs15*clhs52 + clhs16*clhs53) - clhs33*(N1[0]*clhs50 + N1[1]*clhs54) - clhs34*(N1[0]*clhs55 + N1[1]*clhs56) + clhs6*(clhs10*clhs53 + clhs41 + clhs52*clhs7));
    const double clhs58 =     -clhs2*clhs57 + clhs36*clhs50 + clhs38*clhs51;
    const double clhs59 =     -clhs3*clhs57 + clhs36*clhs55 + clhs47*clhs51;
    const double clhs60 =     clhs31*(GPnormal[0]*clhs58 + GPnormal[1]*clhs59);
    const double clhs61 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs64 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs65 =     -N1[1];
    const double clhs66 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs67 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     clhs1*(clhs14*(clhs15*clhs63 + clhs16*clhs64 + clhs65) - clhs33*(N1[0]*clhs61 + N1[1]*clhs66) - clhs34*(N1[0]*clhs67 + N1[1]*clhs68) + clhs6*(clhs10*clhs64 + clhs63*clhs7));
    const double clhs70 =     -clhs2*clhs69 + clhs36*clhs61 + clhs38*clhs62;
    const double clhs71 =     -clhs3*clhs69 + clhs36*clhs67 + clhs47*clhs62;
    const double clhs72 =     clhs31*(GPnormal[0]*clhs70 + GPnormal[1]*clhs71);
    const double clhs73 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs74 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs75 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs76 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs77 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs78 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs79 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs80 =     clhs1*(clhs14*(clhs15*clhs75 + clhs16*clhs76) - clhs33*(N1[0]*clhs73 + N1[1]*clhs77) - clhs34*(N1[0]*clhs78 + N1[1]*clhs79) + clhs6*(clhs10*clhs76 + clhs65 + clhs7*clhs75));
    const double clhs81 =     -clhs2*clhs80 + clhs36*clhs73 + clhs38*clhs74;
    const double clhs82 =     -clhs3*clhs80 + clhs36*clhs78 + clhs47*clhs74;
    const double clhs83 =     clhs31*(GPnormal[0]*clhs81 + GPnormal[1]*clhs82);
    const double clhs84 =     Phi[0]*clhs0*clhs1*(GPtangent1[0]*clhs2 + GPtangent1[1]*clhs3);
    const double clhs85 =     -clhs17*clhs84;
    const double clhs86 =     -clhs21*clhs84;
    const double clhs87 =     -clhs25*clhs84;
    const double clhs88 =     -clhs29*clhs84;
    const double clhs89 =     clhs31*(GPtangent1[0]*clhs46 + GPtangent1[1]*clhs48);
    const double clhs90 =     clhs31*(GPtangent1[0]*clhs58 + GPtangent1[1]*clhs59);
    const double clhs91 =     clhs31*(GPtangent1[0]*clhs70 + GPtangent1[1]*clhs71);
    const double clhs92 =     clhs31*(GPtangent1[0]*clhs81 + GPtangent1[1]*clhs82);
    const double clhs93 =     Phi[1]*clhs0*clhs1*(GPnormal[0]*clhs13 + GPnormal[1]*clhs5);
    const double clhs94 =     -clhs17*clhs93;
    const double clhs95 =     -clhs21*clhs93;
    const double clhs96 =     -clhs25*clhs93;
    const double clhs97 =     -clhs29*clhs93;
    const double clhs98 =     Phi[1]*clhs0;
    const double clhs99 =     clhs13*clhs35;
    const double clhs100 =     -clhs13*clhs45 + clhs36*clhs42 + clhs37*clhs99;
    const double clhs101 =     clhs35*clhs5;
    const double clhs102 =     clhs101*clhs37 + clhs36*clhs44 - clhs45*clhs5;
    const double clhs103 =     clhs98*(GPnormal[0]*clhs100 + GPnormal[1]*clhs102);
    const double clhs104 =     -clhs13*clhs57 + clhs36*clhs54 + clhs51*clhs99;
    const double clhs105 =     clhs101*clhs51 + clhs36*clhs56 - clhs5*clhs57;
    const double clhs106 =     clhs98*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    const double clhs107 =     -clhs13*clhs69 + clhs36*clhs66 + clhs62*clhs99;
    const double clhs108 =     clhs101*clhs62 + clhs36*clhs68 - clhs5*clhs69;
    const double clhs109 =     clhs98*(GPnormal[0]*clhs107 + GPnormal[1]*clhs108);
    const double clhs110 =     -clhs13*clhs80 + clhs36*clhs77 + clhs74*clhs99;
    const double clhs111 =     clhs101*clhs74 + clhs36*clhs79 - clhs5*clhs80;
    const double clhs112 =     clhs98*(GPnormal[0]*clhs110 + GPnormal[1]*clhs111);
    const double clhs113 =     Phi[1]*clhs0*clhs1*(GPtangent1[0]*clhs13 + GPtangent1[1]*clhs5);
    const double clhs114 =     -clhs113*clhs17;
    const double clhs115 =     -clhs113*clhs21;
    const double clhs116 =     -clhs113*clhs25;
    const double clhs117 =     -clhs113*clhs29;
    const double clhs118 =     clhs98*(GPtangent1[0]*clhs100 + GPtangent1[1]*clhs102);
    const double clhs119 =     clhs98*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    const double clhs120 =     clhs98*(GPtangent1[0]*clhs107 + GPtangent1[1]*clhs108);
    const double clhs121 =     clhs98*(GPtangent1[0]*clhs110 + GPtangent1[1]*clhs111);

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
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
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
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
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
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
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
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
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
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
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
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
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
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
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
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(8,0)=clhs18;
    lhs(8,1)=clhs22;
    lhs(8,2)=clhs26;
    lhs(8,3)=clhs30;
    lhs(8,4)=clhs49;
    lhs(8,5)=clhs60;
    lhs(8,6)=clhs72;
    lhs(8,7)=clhs83;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs85;
    lhs(9,1)=clhs86;
    lhs(9,2)=clhs87;
    lhs(9,3)=clhs88;
    lhs(9,4)=clhs89;
    lhs(9,5)=clhs90;
    lhs(9,6)=clhs91;
    lhs(9,7)=clhs92;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs94;
    lhs(10,1)=clhs95;
    lhs(10,2)=clhs96;
    lhs(10,3)=clhs97;
    lhs(10,4)=clhs103;
    lhs(10,5)=clhs106;
    lhs(10,6)=clhs109;
    lhs(10,7)=clhs112;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs114;
    lhs(11,1)=clhs115;
    lhs(11,2)=clhs116;
    lhs(11,3)=clhs117;
    lhs(11,4)=clhs118;
    lhs(11,5)=clhs119;
    lhs(11,6)=clhs120;
    lhs(11,7)=clhs121;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs18;
    lhs(12,1)=clhs22;
    lhs(12,2)=clhs26;
    lhs(12,3)=clhs30;
    lhs(12,4)=clhs49;
    lhs(12,5)=clhs60;
    lhs(12,6)=clhs72;
    lhs(12,7)=clhs83;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs85;
    lhs(13,1)=clhs86;
    lhs(13,2)=clhs87;
    lhs(13,3)=clhs88;
    lhs(13,4)=clhs89;
    lhs(13,5)=clhs90;
    lhs(13,6)=clhs91;
    lhs(13,7)=clhs92;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs94;
    lhs(14,1)=clhs95;
    lhs(14,2)=clhs96;
    lhs(14,3)=clhs97;
    lhs(14,4)=clhs103;
    lhs(14,5)=clhs106;
    lhs(14,6)=clhs109;
    lhs(14,7)=clhs112;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs114;
    lhs(15,1)=clhs115;
    lhs(15,2)=clhs116;
    lhs(15,3)=clhs117;
    lhs(15,4)=clhs118;
    lhs(15,5)=clhs119;
    lhs(15,6)=clhs120;
    lhs(15,7)=clhs121;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;

    
    return lhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointSlipLHS(
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
    bounded_matrix<double,16,16> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Vector normalmasterg   = prod(trans(normalmaster), N2);
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = rContactData.X1;
    const Matrix X2 = rContactData.X2;
    const Matrix u1 = rContactData.u1;
    const Matrix u2 = rContactData.u2;
    const Matrix v1 = rContactData.v1;
    const Matrix v2 = rContactData.v2;

//substitute_constants_derivatives_variables_slip 
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
    const double clhs18 =     -clhs17*clhs4;
    const double clhs19 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     clhs14*(clhs15*clhs19 + clhs16*clhs20) + clhs6*(clhs10*clhs20 + clhs19*clhs7 + clhs8);
    const double clhs22 =     -clhs21*clhs4;
    const double clhs23 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs24 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs25 =     clhs14*(clhs11 + clhs15*clhs23 + clhs16*clhs24) + clhs6*(clhs10*clhs24 + clhs23*clhs7);
    const double clhs26 =     -clhs25*clhs4;
    const double clhs27 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs28 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs29 =     clhs14*(clhs15*clhs27 + clhs16*clhs28) + clhs6*(clhs10*clhs28 + clhs11 + clhs27*clhs7);
    const double clhs30 =     -clhs29*clhs4;
    const double clhs31 =     Phi[0]*clhs0;
    const double clhs32 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs33 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs11*clhs16 - clhs15*clhs8;
    const double clhs34 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs10*clhs11 - clhs7*clhs8;
    const double clhs35 =     clhs14*clhs33 + clhs34*clhs6;
    const double clhs36 =     clhs1*clhs35;
    const double clhs37 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs38 =     clhs2*clhs35;
    const double clhs39 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs40 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs41 =     -N1[0];
    const double clhs42 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs43 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs44 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs45 =     clhs1*(clhs14*(clhs15*clhs39 + clhs16*clhs40 + clhs41) - clhs33*(N1[0]*clhs32 + N1[1]*clhs42) - clhs34*(N1[0]*clhs43 + N1[1]*clhs44) + clhs6*(clhs10*clhs40 + clhs39*clhs7));
    const double clhs46 =     -clhs2*clhs45 + clhs32*clhs36 + clhs37*clhs38;
    const double clhs47 =     clhs3*clhs35;
    const double clhs48 =     -clhs3*clhs45 + clhs36*clhs43 + clhs37*clhs47;
    const double clhs49 =     clhs31*(GPnormal[0]*clhs46 + GPnormal[1]*clhs48);
    const double clhs50 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs53 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs54 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     clhs1*(clhs14*(clhs15*clhs52 + clhs16*clhs53) - clhs33*(N1[0]*clhs50 + N1[1]*clhs54) - clhs34*(N1[0]*clhs55 + N1[1]*clhs56) + clhs6*(clhs10*clhs53 + clhs41 + clhs52*clhs7));
    const double clhs58 =     -clhs2*clhs57 + clhs36*clhs50 + clhs38*clhs51;
    const double clhs59 =     -clhs3*clhs57 + clhs36*clhs55 + clhs47*clhs51;
    const double clhs60 =     clhs31*(GPnormal[0]*clhs58 + GPnormal[1]*clhs59);
    const double clhs61 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs64 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs65 =     -N1[1];
    const double clhs66 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs67 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     clhs1*(clhs14*(clhs15*clhs63 + clhs16*clhs64 + clhs65) - clhs33*(N1[0]*clhs61 + N1[1]*clhs66) - clhs34*(N1[0]*clhs67 + N1[1]*clhs68) + clhs6*(clhs10*clhs64 + clhs63*clhs7));
    const double clhs70 =     -clhs2*clhs69 + clhs36*clhs61 + clhs38*clhs62;
    const double clhs71 =     -clhs3*clhs69 + clhs36*clhs67 + clhs47*clhs62;
    const double clhs72 =     clhs31*(GPnormal[0]*clhs70 + GPnormal[1]*clhs71);
    const double clhs73 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs74 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs75 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs76 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs77 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs78 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs79 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs80 =     clhs1*(clhs14*(clhs15*clhs75 + clhs16*clhs76) - clhs33*(N1[0]*clhs73 + N1[1]*clhs77) - clhs34*(N1[0]*clhs78 + N1[1]*clhs79) + clhs6*(clhs10*clhs76 + clhs65 + clhs7*clhs75));
    const double clhs81 =     -clhs2*clhs80 + clhs36*clhs73 + clhs38*clhs74;
    const double clhs82 =     -clhs3*clhs80 + clhs36*clhs78 + clhs47*clhs74;
    const double clhs83 =     clhs31*(GPnormal[0]*clhs81 + GPnormal[1]*clhs82);
    const double clhs84 =     Phi[0]*clhs0*clhs1*(GPtangent1[0]*clhs2 + GPtangent1[1]*clhs3);
    const double clhs85 =     -clhs17*clhs84;
    const double clhs86 =     -clhs21*clhs84;
    const double clhs87 =     -clhs25*clhs84;
    const double clhs88 =     -clhs29*clhs84;
    const double clhs89 =     clhs31*(GPtangent1[0]*clhs46 + GPtangent1[1]*clhs48);
    const double clhs90 =     clhs31*(GPtangent1[0]*clhs58 + GPtangent1[1]*clhs59);
    const double clhs91 =     clhs31*(GPtangent1[0]*clhs70 + GPtangent1[1]*clhs71);
    const double clhs92 =     clhs31*(GPtangent1[0]*clhs81 + GPtangent1[1]*clhs82);
    const double clhs93 =     Phi[1]*clhs0*clhs1*(GPnormal[0]*clhs13 + GPnormal[1]*clhs5);
    const double clhs94 =     -clhs17*clhs93;
    const double clhs95 =     -clhs21*clhs93;
    const double clhs96 =     -clhs25*clhs93;
    const double clhs97 =     -clhs29*clhs93;
    const double clhs98 =     Phi[1]*clhs0;
    const double clhs99 =     clhs13*clhs35;
    const double clhs100 =     -clhs13*clhs45 + clhs36*clhs42 + clhs37*clhs99;
    const double clhs101 =     clhs35*clhs5;
    const double clhs102 =     clhs101*clhs37 + clhs36*clhs44 - clhs45*clhs5;
    const double clhs103 =     clhs98*(GPnormal[0]*clhs100 + GPnormal[1]*clhs102);
    const double clhs104 =     -clhs13*clhs57 + clhs36*clhs54 + clhs51*clhs99;
    const double clhs105 =     clhs101*clhs51 + clhs36*clhs56 - clhs5*clhs57;
    const double clhs106 =     clhs98*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    const double clhs107 =     -clhs13*clhs69 + clhs36*clhs66 + clhs62*clhs99;
    const double clhs108 =     clhs101*clhs62 + clhs36*clhs68 - clhs5*clhs69;
    const double clhs109 =     clhs98*(GPnormal[0]*clhs107 + GPnormal[1]*clhs108);
    const double clhs110 =     -clhs13*clhs80 + clhs36*clhs77 + clhs74*clhs99;
    const double clhs111 =     clhs101*clhs74 + clhs36*clhs79 - clhs5*clhs80;
    const double clhs112 =     clhs98*(GPnormal[0]*clhs110 + GPnormal[1]*clhs111);
    const double clhs113 =     Phi[1]*clhs0*clhs1*(GPtangent1[0]*clhs13 + GPtangent1[1]*clhs5);
    const double clhs114 =     -clhs113*clhs17;
    const double clhs115 =     -clhs113*clhs21;
    const double clhs116 =     -clhs113*clhs25;
    const double clhs117 =     -clhs113*clhs29;
    const double clhs118 =     clhs98*(GPtangent1[0]*clhs100 + GPtangent1[1]*clhs102);
    const double clhs119 =     clhs98*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    const double clhs120 =     clhs98*(GPtangent1[0]*clhs107 + GPtangent1[1]*clhs108);
    const double clhs121 =     clhs98*(GPtangent1[0]*clhs110 + GPtangent1[1]*clhs111);

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
    lhs(0,12)=0;
    lhs(0,13)=0;
    lhs(0,14)=0;
    lhs(0,15)=0;
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
    lhs(1,12)=0;
    lhs(1,13)=0;
    lhs(1,14)=0;
    lhs(1,15)=0;
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
    lhs(2,12)=0;
    lhs(2,13)=0;
    lhs(2,14)=0;
    lhs(2,15)=0;
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
    lhs(3,12)=0;
    lhs(3,13)=0;
    lhs(3,14)=0;
    lhs(3,15)=0;
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
    lhs(4,12)=0;
    lhs(4,13)=0;
    lhs(4,14)=0;
    lhs(4,15)=0;
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
    lhs(5,12)=0;
    lhs(5,13)=0;
    lhs(5,14)=0;
    lhs(5,15)=0;
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
    lhs(6,12)=0;
    lhs(6,13)=0;
    lhs(6,14)=0;
    lhs(6,15)=0;
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
    lhs(7,12)=0;
    lhs(7,13)=0;
    lhs(7,14)=0;
    lhs(7,15)=0;
    lhs(8,0)=clhs18;
    lhs(8,1)=clhs22;
    lhs(8,2)=clhs26;
    lhs(8,3)=clhs30;
    lhs(8,4)=clhs49;
    lhs(8,5)=clhs60;
    lhs(8,6)=clhs72;
    lhs(8,7)=clhs83;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(8,12)=0;
    lhs(8,13)=0;
    lhs(8,14)=0;
    lhs(8,15)=0;
    lhs(9,0)=clhs85;
    lhs(9,1)=clhs86;
    lhs(9,2)=clhs87;
    lhs(9,3)=clhs88;
    lhs(9,4)=clhs89;
    lhs(9,5)=clhs90;
    lhs(9,6)=clhs91;
    lhs(9,7)=clhs92;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(9,12)=0;
    lhs(9,13)=0;
    lhs(9,14)=0;
    lhs(9,15)=0;
    lhs(10,0)=clhs94;
    lhs(10,1)=clhs95;
    lhs(10,2)=clhs96;
    lhs(10,3)=clhs97;
    lhs(10,4)=clhs103;
    lhs(10,5)=clhs106;
    lhs(10,6)=clhs109;
    lhs(10,7)=clhs112;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(10,12)=0;
    lhs(10,13)=0;
    lhs(10,14)=0;
    lhs(10,15)=0;
    lhs(11,0)=clhs114;
    lhs(11,1)=clhs115;
    lhs(11,2)=clhs116;
    lhs(11,3)=clhs117;
    lhs(11,4)=clhs118;
    lhs(11,5)=clhs119;
    lhs(11,6)=clhs120;
    lhs(11,7)=clhs121;
    lhs(11,8)=0;
    lhs(11,9)=0;
    lhs(11,10)=0;
    lhs(11,11)=0;
    lhs(11,12)=0;
    lhs(11,13)=0;
    lhs(11,14)=0;
    lhs(11,15)=0;
    lhs(12,0)=clhs18;
    lhs(12,1)=clhs22;
    lhs(12,2)=clhs26;
    lhs(12,3)=clhs30;
    lhs(12,4)=clhs49;
    lhs(12,5)=clhs60;
    lhs(12,6)=clhs72;
    lhs(12,7)=clhs83;
    lhs(12,8)=0;
    lhs(12,9)=0;
    lhs(12,10)=0;
    lhs(12,11)=0;
    lhs(12,12)=0;
    lhs(12,13)=0;
    lhs(12,14)=0;
    lhs(12,15)=0;
    lhs(13,0)=clhs85;
    lhs(13,1)=clhs86;
    lhs(13,2)=clhs87;
    lhs(13,3)=clhs88;
    lhs(13,4)=clhs89;
    lhs(13,5)=clhs90;
    lhs(13,6)=clhs91;
    lhs(13,7)=clhs92;
    lhs(13,8)=0;
    lhs(13,9)=0;
    lhs(13,10)=0;
    lhs(13,11)=0;
    lhs(13,12)=0;
    lhs(13,13)=0;
    lhs(13,14)=0;
    lhs(13,15)=0;
    lhs(14,0)=clhs94;
    lhs(14,1)=clhs95;
    lhs(14,2)=clhs96;
    lhs(14,3)=clhs97;
    lhs(14,4)=clhs103;
    lhs(14,5)=clhs106;
    lhs(14,6)=clhs109;
    lhs(14,7)=clhs112;
    lhs(14,8)=0;
    lhs(14,9)=0;
    lhs(14,10)=0;
    lhs(14,11)=0;
    lhs(14,12)=0;
    lhs(14,13)=0;
    lhs(14,14)=0;
    lhs(14,15)=0;
    lhs(15,0)=clhs114;
    lhs(15,1)=clhs115;
    lhs(15,2)=clhs116;
    lhs(15,3)=clhs117;
    lhs(15,4)=clhs118;
    lhs(15,5)=clhs119;
    lhs(15,6)=clhs120;
    lhs(15,7)=clhs121;
    lhs(15,8)=0;
    lhs(15,9)=0;
    lhs(15,10)=0;
    lhs(15,11)=0;
    lhs(15,12)=0;
    lhs(15,13)=0;
    lhs(15,14)=0;
    lhs(15,15)=0;

    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,16,16> ComputeGaussPointInactiveLHS(
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
    bounded_matrix<double,16,16> lhs;
    
//     const Matrix normalmaster    = rContactData.NormalsMaster;
//     const Vector normalmasterg   = prod(trans(normalmaster), N2);
//     const Matrix normalslave     = rContactData.NormalsSlave;
//     const Matrix tan1slave       = rContactData.Tangent1Slave;
//     const Matrix lm              = rContactData.LagrangeMultipliers;
//     const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
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
//substitute_constants_derivatives_variables_inactive 
//substitute_derivatives_variables_inactive 
//substitute_inactive_lhs
    
    return lhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointActiveRHS(
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
    array_1d<double,16> rhs;
    
    const Matrix normalslave    = rContactData.NormalsSlave;
    const Matrix tan1slave      = rContactData.Tangent1Slave;
    const Matrix lm             = rContactData.LagrangeMultipliers;
    const Matrix dlm            = rContactData.DoubleLagrangeMultipliers;
    const double epsilon        = rContactData.epsilon;
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]*dlm(0,0) + Phi[1]*dlm(1,0);
    const double crhs3 =     Phi[0]*lm(0,0);
    const double crhs4 =     Phi[1]*lm(1,0);
    const double crhs5 =     crhs2 + crhs3 + crhs4;
    const double crhs6 =     crhs1*crhs5;
    const double crhs7 =     Phi[0]*dlm(0,1) + Phi[1]*dlm(1,1);
    const double crhs8 =     Phi[0]*lm(0,1);
    const double crhs9 =     Phi[1]*lm(1,1);
    const double crhs10 =     crhs7 + crhs8 + crhs9;
    const double crhs11 =     crhs1*crhs10;
    const double crhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs13 =     N1[0]*crhs1;
    const double crhs14 =     N1[1]*crhs1;
    const double crhs15 =     Phi[0]*crhs1;
    const double crhs16 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs17 =     crhs16*normalslave(0,0); // CRHS16*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs18 =     epsilon*(crhs2 - crhs3 - crhs4);
    const double crhs19 =     -crhs17 + crhs18;
    const double crhs20 =     crhs16*normalslave(0,1); // CRHS16*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs21 =     epsilon*(crhs7 - crhs8 - crhs9);
    const double crhs22 =     -crhs20 + crhs21;
    const double crhs23 =     Phi[1]*crhs1;
    const double crhs24 =     crhs16*normalslave(1,0); // CRHS16*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs25 =     crhs18 - crhs24;
    const double crhs26 =     crhs16*normalslave(1,1); // CRHS16*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs27 =     crhs21 - crhs26;
    const double crhs28 =     crhs17 + crhs18;
    const double crhs29 =     crhs20 + crhs21;
    const double crhs30 =     crhs18 + crhs24;
    const double crhs31 =     crhs21 + crhs26;

    rhs[0]=-crhs0*crhs6;
    rhs[1]=-crhs0*crhs11;
    rhs[2]=-crhs12*crhs6;
    rhs[3]=-crhs11*crhs12;
    rhs[4]=crhs13*crhs5;
    rhs[5]=crhs10*crhs13;
    rhs[6]=crhs14*crhs5;
    rhs[7]=crhs10*crhs14;
    rhs[8]=-crhs15*(GPnormal[0]*crhs19 + GPnormal[1]*crhs22);
    rhs[9]=-crhs15*(GPtangent1[0]*crhs19 + GPtangent1[1]*crhs22);
    rhs[10]=-crhs23*(GPnormal[0]*crhs25 + GPnormal[1]*crhs27);
    rhs[11]=-crhs23*(GPtangent1[0]*crhs25 + GPtangent1[1]*crhs27);
    rhs[12]=crhs15*(GPnormal[0]*crhs28 + GPnormal[1]*crhs29);
    rhs[13]=crhs15*(GPtangent1[0]*crhs28 + GPtangent1[1]*crhs29);
    rhs[14]=crhs23*(GPnormal[0]*crhs30 + GPnormal[1]*crhs31);
    rhs[15]=crhs23*(GPtangent1[0]*crhs30 + GPtangent1[1]*crhs31);

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointStickRHS(
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
    array_1d<double,16> rhs;
    
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
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
    const double crhs10 =     -crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs11 =     -crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs12 =     Phi[1]*crhs2*crhs3*crhs8;
    const double crhs13 =     -crhs12*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    const double crhs14 =     -crhs12*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs10;
    rhs[9]=crhs11;
    rhs[10]=crhs13;
    rhs[11]=crhs14;
    rhs[12]=crhs10;
    rhs[13]=crhs11;
    rhs[14]=crhs13;
    rhs[15]=crhs14;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointSlipRHS(
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
    array_1d<double,16> rhs;
    
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
    const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
    const double Dt              = rContactData.Dt;
    
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
    const double crhs10 =     -crhs9*(GPnormal[0]*crhs0 + GPnormal[1]*crhs1);
    const double crhs11 =     -crhs9*(GPtangent1[0]*crhs0 + GPtangent1[1]*crhs1);
    const double crhs12 =     Phi[1]*crhs2*crhs3*crhs8;
    const double crhs13 =     -crhs12*(GPnormal[0]*crhs4 + GPnormal[1]*crhs7);
    const double crhs14 =     -crhs12*(GPtangent1[0]*crhs4 + GPtangent1[1]*crhs7);

    rhs[0]=0;
    rhs[1]=0;
    rhs[2]=0;
    rhs[3]=0;
    rhs[4]=0;
    rhs[5]=0;
    rhs[6]=0;
    rhs[7]=0;
    rhs[8]=crhs10;
    rhs[9]=crhs11;
    rhs[10]=crhs13;
    rhs[11]=crhs14;
    rhs[12]=crhs10;
    rhs[13]=crhs11;
    rhs[14]=crhs13;
    rhs[15]=crhs14;

    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline array_1d<double,16> ComputeGaussPointInactiveRHS(
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
    array_1d<double,16> rhs;
    
//     const Matrix normalslave     = rContactData.NormalsSlave;
//     const Matrix tan1slave       = rContactData.Tangent1Slave;
//     const Matrix lm              = rContactData.LagrangeMultipliers;
//     const Matrix dlm             = rContactData.DoubleLagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
//     const double epsilon  = rContactData.epsilon;
//     const double epsilon_tangent = rContactData.epsilon_tangent;
//     
//     const Vector GPnormal     = prod(trans(normalslave), N1);
//     const Vector GPtangent1   = prod(trans(tan1slave), N1);
//     
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
//     
//substitute_inactive_rhs
    
    return rhs;
}

    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,2,2> ComputeDeltaDe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,2,2> DeltaDe;
    
    const Matrix X1 = rContactData.X1;
    const Matrix u1 = rContactData.u1;

    double DeltaDetJ = 0.0;
    
    if (derivative_index == 0)
    {
        DeltaDetJ = (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 1)
    {
        DeltaDetJ = (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 2)
    {
        DeltaDetJ = (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 3)
    {
        DeltaDetJ = (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    
    DeltaDe(0,0) = DeltaDetJ * N1[0];
    DeltaDe(0,1) = 0;
    DeltaDe(1,0) = 0;
    DeltaDe(1,1) = DeltaDetJ * N1[1];
    
    return DeltaDe;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,2,2> ComputeDeltaMe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,2,2> DeltaMe;
    
    const Matrix X1 = rContactData.X1;
    const Matrix u1 = rContactData.u1;

    double DeltaDetJ = 0.0;
    
    if (derivative_index == 0)
    {
        DeltaDetJ = (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 1)
    {
        DeltaDetJ = (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 2)
    {
        DeltaDetJ = (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    else if (derivative_index == 3)
    {
        DeltaDetJ = (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    }
    
    DeltaMe(0,0) = DeltaDetJ * N1[0] * N1[0];
    DeltaMe(0,1) = DeltaDetJ * N1[0] * N1[1];
    DeltaMe(1,0) = DeltaDetJ * N1[1] * N1[0];
    DeltaMe(1,1) = DeltaDetJ * N1[1] * N1[1];
    
    return DeltaMe;
}

private:
};// class Contact2D2N2NDLM
}
#endif /* KRATOS_CONTACT2D2N2NDLM defined */
