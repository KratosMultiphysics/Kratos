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
        const double detJ, 
        const ContactData& rContactData,
        const double& augmented_normal_lm,
        const double& augmented_tangent_lm,
        const double& integration_point_gap,
        const double& integration_point_slip
        )
{
    bounded_matrix<double,12,12> lhs;
    
    const Matrix normalmaster    = rContactData.NormalsMaster;
    const Matrix normalslave     = rContactData.NormalsSlave;
    const Matrix tan1slave       = rContactData.Tangent1Slave;
    const Matrix lm              = rContactData.LagrangeMultipliers;
//     const double Dt              = rContactData.Dt;
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
    
    const std::vector<Matrix> DeltaAe = rContactData.DeltaAe;
    Vector DerivativePhiu00 = prod(DeltaAe[0],N1);
    Vector DerivativePhiu01 = prod(DeltaAe[1],N1);
    Vector DerivativePhiu10 = prod(DeltaAe[2],N1);
    Vector DerivativePhiu11 = prod(DeltaAe[3],N1);

//substitute_constants_derivatives_variables_active 
    const double Dnormalmaster11u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double DPhi1u111 =     DerivativePhiu11[1];
    const double DPhi1u110 =     DerivativePhiu10[1];
    const double DPhi1u101 =     DerivativePhiu01[1];
    const double DPhi1u100 =     DerivativePhiu00[1];
    const double DPhi0u111 =     DerivativePhiu11[0];
    const double DPhi0u110 =     DerivativePhiu10[0];
    const double DPhi0u101 =     DerivativePhiu01[0];
    const double DPhi0u100 =     DerivativePhiu00[0];
    const double DN21u211 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
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
    const double clhs3 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs3*lm(0,0) + clhs4*lm(1,0);
    const double clhs6 =     clhs2*clhs5;
    const double clhs7 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs8 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs9 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs10 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs11 =     clhs10*clhs5;
    const double clhs12 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs13 =     clhs0*clhs2;
    const double clhs14 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs15 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs16 =     clhs14*lm(0,0) + clhs15*lm(1,0);
    const double clhs17 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs18 =     clhs17*clhs5;
    const double clhs19 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs20 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs21 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs22 =     clhs20*lm(0,0) + clhs21*lm(1,0);
    const double clhs23 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs24 =     clhs23*clhs5;
    const double clhs25 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs26 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs27 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs28 =     clhs26*lm(0,0) + clhs27*lm(1,0);
    const double clhs29 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs30 =     clhs29*clhs5;
    const double clhs31 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs32 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs33 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs34 =     clhs32*lm(0,0) + clhs33*lm(1,0);
    const double clhs35 =     clhs13*clhs3;
    const double clhs36 =     clhs13*clhs4;
    const double clhs37 =     clhs3*lm(0,1) + clhs4*lm(1,1);
    const double clhs38 =     clhs2*clhs37;
    const double clhs39 =     clhs10*clhs37;
    const double clhs40 =     clhs14*lm(0,1) + clhs15*lm(1,1);
    const double clhs41 =     clhs17*clhs37;
    const double clhs42 =     clhs20*lm(0,1) + clhs21*lm(1,1);
    const double clhs43 =     clhs23*clhs37;
    const double clhs44 =     clhs26*lm(0,1) + clhs27*lm(1,1);
    const double clhs45 =     clhs29*clhs37;
    const double clhs46 =     clhs32*lm(0,1) + clhs33*lm(1,1);
    const double clhs47 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs48 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs49 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs50 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs51 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs52 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs53 =     clhs2*clhs47;
    const double clhs54 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs55 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs56 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs57 =     clhs3*clhs53;
    const double clhs58 =     clhs4*clhs53;
    const double clhs59 =     clhs11 + clhs16*clhs2;
    const double clhs60 =     clhs18 + clhs2*clhs22;
    const double clhs61 =     clhs2*clhs28 + clhs24;
    const double clhs62 =     clhs2*clhs34 + clhs30;
    const double clhs63 =     N1[0]*clhs2;
    const double clhs64 =     -clhs3*clhs63;
    const double clhs65 =     -clhs4*clhs63;
    const double clhs66 =     clhs2*clhs40 + clhs39;
    const double clhs67 =     clhs2*clhs42 + clhs41;
    const double clhs68 =     clhs2*clhs44 + clhs43;
    const double clhs69 =     clhs2*clhs46 + clhs45;
    const double clhs70 =     N1[1]*clhs2;
    const double clhs71 =     -clhs3*clhs70;
    const double clhs72 =     -clhs4*clhs70;
    const double clhs73 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs74 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs75 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs76 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs77 =     clhs2*clhs3*(GPnormal[0]*clhs75 + GPnormal[1]*clhs76);
    const double clhs78 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs79 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs80 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs81 =     clhs2*clhs3*clhs73;
    const double clhs82 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs83 =     clhs2*clhs3*clhs75;
    const double clhs84 =     clhs3*clhs73*clhs75;
    const double clhs85 =     clhs2*clhs73*clhs75;
    const double clhs86 =     clhs10*clhs84 + clhs14*clhs85 + clhs81*Dnormalslave00u100 + clhs82*clhs83; // CLHS10*CLHS84 + CLHS14*CLHS85 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS82*CLHS83
    const double clhs87 =     clhs2*clhs3*clhs76;
    const double clhs88 =     clhs3*clhs73*clhs76;
    const double clhs89 =     clhs2*clhs73*clhs76;
    const double clhs90 =     clhs10*clhs88 + clhs14*clhs89 + clhs81*Dnormalslave01u100 + clhs82*clhs87; // CLHS10*CLHS88 + CLHS14*CLHS89 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS82*CLHS87
    const double clhs91 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs92 =     clhs17*clhs84 + clhs20*clhs85 + clhs81*Dnormalslave00u101 + clhs83*clhs91; // CLHS17*CLHS84 + CLHS20*CLHS85 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS83*CLHS91
    const double clhs93 =     clhs17*clhs88 + clhs20*clhs89 + clhs81*Dnormalslave01u101 + clhs87*clhs91; // CLHS17*CLHS88 + CLHS20*CLHS89 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS87*CLHS91
    const double clhs94 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs95 =     clhs23*clhs84 + clhs26*clhs85 + clhs81*Dnormalslave00u110 + clhs83*clhs94; // CLHS23*CLHS84 + CLHS26*CLHS85 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS83*CLHS94
    const double clhs96 =     clhs23*clhs88 + clhs26*clhs89 + clhs81*Dnormalslave01u110 + clhs87*clhs94; // CLHS23*CLHS88 + CLHS26*CLHS89 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS87*CLHS94
    const double clhs97 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs98 =     clhs29*clhs84 + clhs32*clhs85 + clhs81*Dnormalslave00u111 + clhs83*clhs97; // CLHS29*CLHS84 + CLHS32*CLHS85 + CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS83*CLHS97
    const double clhs99 =     clhs29*clhs88 + clhs32*clhs89 + clhs81*Dnormalslave01u111 + clhs87*clhs97; // CLHS29*CLHS88 + CLHS32*CLHS89 + CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS87*CLHS97
    const double clhs100 =     clhs2*clhs3*(GPtangent1[0]*clhs75 + GPtangent1[1]*clhs76);
    const double clhs101 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs102 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs103 =     clhs2*clhs4*(GPnormal[0]*clhs101 + GPnormal[1]*clhs102);
    const double clhs104 =     clhs2*clhs4*clhs73;
    const double clhs105 =     clhs101*clhs2*clhs4;
    const double clhs106 =     clhs101*clhs4*clhs73;
    const double clhs107 =     clhs101*clhs2*clhs73;
    const double clhs108 =     clhs10*clhs106 + clhs104*Dnormalslave10u100 + clhs105*clhs82 + clhs107*clhs15; // CLHS10*CLHS106 + CLHS104*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS105*CLHS82 + CLHS107*CLHS15
    const double clhs109 =     clhs102*clhs2*clhs4;
    const double clhs110 =     clhs102*clhs4*clhs73;
    const double clhs111 =     clhs102*clhs2*clhs73;
    const double clhs112 =     clhs10*clhs110 + clhs104*Dnormalslave11u100 + clhs109*clhs82 + clhs111*clhs15; // CLHS10*CLHS110 + CLHS104*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS109*CLHS82 + CLHS111*CLHS15
    const double clhs113 =     clhs104*Dnormalslave10u101 + clhs105*clhs91 + clhs106*clhs17 + clhs107*clhs21; // CLHS104*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS105*CLHS91 + CLHS106*CLHS17 + CLHS107*CLHS21
    const double clhs114 =     clhs104*Dnormalslave11u101 + clhs109*clhs91 + clhs110*clhs17 + clhs111*clhs21; // CLHS104*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS109*CLHS91 + CLHS110*CLHS17 + CLHS111*CLHS21
    const double clhs115 =     clhs104*Dnormalslave10u110 + clhs105*clhs94 + clhs106*clhs23 + clhs107*clhs27; // CLHS104*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS105*CLHS94 + CLHS106*CLHS23 + CLHS107*CLHS27
    const double clhs116 =     clhs104*Dnormalslave11u110 + clhs109*clhs94 + clhs110*clhs23 + clhs111*clhs27; // CLHS104*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS109*CLHS94 + CLHS110*CLHS23 + CLHS111*CLHS27
    const double clhs117 =     clhs104*Dnormalslave10u111 + clhs105*clhs97 + clhs106*clhs29 + clhs107*clhs33; // CLHS104*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS105*CLHS97 + CLHS106*CLHS29 + CLHS107*CLHS33
    const double clhs118 =     clhs104*Dnormalslave11u111 + clhs109*clhs97 + clhs110*clhs29 + clhs111*clhs33; // CLHS104*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS109*CLHS97 + CLHS110*CLHS29 + CLHS111*CLHS33
    const double clhs119 =     clhs2*clhs4*(GPtangent1[0]*clhs101 + GPtangent1[1]*clhs102);

    lhs(0,0)=clhs1*clhs6;
    lhs(0,1)=clhs6*clhs7;
    lhs(0,2)=clhs6*clhs8;
    lhs(0,3)=clhs6*clhs9;
    lhs(0,4)=clhs0*clhs11 + clhs12*clhs6 + clhs13*clhs16;
    lhs(0,5)=clhs0*clhs18 + clhs13*clhs22 + clhs19*clhs6;
    lhs(0,6)=clhs0*clhs24 + clhs13*clhs28 + clhs25*clhs6;
    lhs(0,7)=clhs0*clhs30 + clhs13*clhs34 + clhs31*clhs6;
    lhs(0,8)=clhs35;
    lhs(0,9)=0;
    lhs(0,10)=clhs36;
    lhs(0,11)=0;
    lhs(1,0)=clhs1*clhs38;
    lhs(1,1)=clhs38*clhs7;
    lhs(1,2)=clhs38*clhs8;
    lhs(1,3)=clhs38*clhs9;
    lhs(1,4)=clhs0*clhs39 + clhs12*clhs38 + clhs13*clhs40;
    lhs(1,5)=clhs0*clhs41 + clhs13*clhs42 + clhs19*clhs38;
    lhs(1,6)=clhs0*clhs43 + clhs13*clhs44 + clhs25*clhs38;
    lhs(1,7)=clhs0*clhs45 + clhs13*clhs46 + clhs31*clhs38;
    lhs(1,8)=0;
    lhs(1,9)=clhs35;
    lhs(1,10)=0;
    lhs(1,11)=clhs36;
    lhs(2,0)=clhs48*clhs6;
    lhs(2,1)=clhs49*clhs6;
    lhs(2,2)=clhs50*clhs6;
    lhs(2,3)=clhs51*clhs6;
    lhs(2,4)=clhs11*clhs47 + clhs16*clhs53 + clhs52*clhs6;
    lhs(2,5)=clhs18*clhs47 + clhs22*clhs53 + clhs54*clhs6;
    lhs(2,6)=clhs24*clhs47 + clhs28*clhs53 + clhs55*clhs6;
    lhs(2,7)=clhs30*clhs47 + clhs34*clhs53 + clhs56*clhs6;
    lhs(2,8)=clhs57;
    lhs(2,9)=0;
    lhs(2,10)=clhs58;
    lhs(2,11)=0;
    lhs(3,0)=clhs38*clhs48;
    lhs(3,1)=clhs38*clhs49;
    lhs(3,2)=clhs38*clhs50;
    lhs(3,3)=clhs38*clhs51;
    lhs(3,4)=clhs38*clhs52 + clhs39*clhs47 + clhs40*clhs53;
    lhs(3,5)=clhs38*clhs54 + clhs41*clhs47 + clhs42*clhs53;
    lhs(3,6)=clhs38*clhs55 + clhs43*clhs47 + clhs44*clhs53;
    lhs(3,7)=clhs38*clhs56 + clhs45*clhs47 + clhs46*clhs53;
    lhs(3,8)=0;
    lhs(3,9)=clhs57;
    lhs(3,10)=0;
    lhs(3,11)=clhs58;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=-N1[0]*clhs59;
    lhs(4,5)=-N1[0]*clhs60;
    lhs(4,6)=-N1[0]*clhs61;
    lhs(4,7)=-N1[0]*clhs62;
    lhs(4,8)=clhs64;
    lhs(4,9)=0;
    lhs(4,10)=clhs65;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=-N1[0]*clhs66;
    lhs(5,5)=-N1[0]*clhs67;
    lhs(5,6)=-N1[0]*clhs68;
    lhs(5,7)=-N1[0]*clhs69;
    lhs(5,8)=0;
    lhs(5,9)=clhs64;
    lhs(5,10)=0;
    lhs(5,11)=clhs65;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=-N1[1]*clhs59;
    lhs(6,5)=-N1[1]*clhs60;
    lhs(6,6)=-N1[1]*clhs61;
    lhs(6,7)=-N1[1]*clhs62;
    lhs(6,8)=clhs71;
    lhs(6,9)=0;
    lhs(6,10)=clhs72;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=-N1[1]*clhs66;
    lhs(7,5)=-N1[1]*clhs67;
    lhs(7,6)=-N1[1]*clhs68;
    lhs(7,7)=-N1[1]*clhs69;
    lhs(7,8)=0;
    lhs(7,9)=clhs71;
    lhs(7,10)=0;
    lhs(7,11)=clhs72;
    lhs(8,0)=-clhs74*clhs77;
    lhs(8,1)=-clhs77*clhs78;
    lhs(8,2)=-clhs77*clhs79;
    lhs(8,3)=-clhs77*clhs80;
    lhs(8,4)=-GPnormal[0]*clhs86 - GPnormal[1]*clhs90;
    lhs(8,5)=-GPnormal[0]*clhs92 - GPnormal[1]*clhs93;
    lhs(8,6)=-GPnormal[0]*clhs95 - GPnormal[1]*clhs96;
    lhs(8,7)=-GPnormal[0]*clhs98 - GPnormal[1]*clhs99;
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs100*clhs74;
    lhs(9,1)=-clhs100*clhs78;
    lhs(9,2)=-clhs100*clhs79;
    lhs(9,3)=-clhs100*clhs80;
    lhs(9,4)=-GPtangent1[0]*clhs86 - GPtangent1[1]*clhs90;
    lhs(9,5)=-GPtangent1[0]*clhs92 - GPtangent1[1]*clhs93;
    lhs(9,6)=-GPtangent1[0]*clhs95 - GPtangent1[1]*clhs96;
    lhs(9,7)=-GPtangent1[0]*clhs98 - GPtangent1[1]*clhs99;
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs103*clhs74;
    lhs(10,1)=-clhs103*clhs78;
    lhs(10,2)=-clhs103*clhs79;
    lhs(10,3)=-clhs103*clhs80;
    lhs(10,4)=-GPnormal[0]*clhs108 - GPnormal[1]*clhs112;
    lhs(10,5)=-GPnormal[0]*clhs113 - GPnormal[1]*clhs114;
    lhs(10,6)=-GPnormal[0]*clhs115 - GPnormal[1]*clhs116;
    lhs(10,7)=-GPnormal[0]*clhs117 - GPnormal[1]*clhs118;
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs119*clhs74;
    lhs(11,1)=-clhs119*clhs78;
    lhs(11,2)=-clhs119*clhs79;
    lhs(11,3)=-clhs119*clhs80;
    lhs(11,4)=-GPtangent1[0]*clhs108 - GPtangent1[1]*clhs112;
    lhs(11,5)=-GPtangent1[0]*clhs113 - GPtangent1[1]*clhs114;
    lhs(11,6)=-GPtangent1[0]*clhs115 - GPtangent1[1]*clhs116;
    lhs(11,7)=-GPtangent1[0]*clhs117 - GPtangent1[1]*clhs118;
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
    
    const std::vector<Matrix> DeltaAe = rContactData.DeltaAe;
    Vector DerivativePhiu00 = prod(DeltaAe[0],N1);
    Vector DerivativePhiu01 = prod(DeltaAe[1],N1);
    Vector DerivativePhiu10 = prod(DeltaAe[2],N1);
    Vector DerivativePhiu11 = prod(DeltaAe[3],N1);
    
//substitute_constants_derivatives_variables_stick 
    const double Dnormalmaster11u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double DPhi1u111 =     DerivativePhiu11[1];
    const double DPhi1u110 =     DerivativePhiu10[1];
    const double DPhi1u101 =     DerivativePhiu01[1];
    const double DPhi1u100 =     DerivativePhiu00[1];
    const double DPhi0u111 =     DerivativePhiu11[0];
    const double DPhi0u110 =     DerivativePhiu10[0];
    const double DPhi0u101 =     DerivativePhiu01[0];
    const double DPhi0u100 =     DerivativePhiu00[0];
    const double DN21u211 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[0]*clhs4 + N1[1]*clhs6;
    const double clhs8 =     Dt*v2(0,1);
    const double clhs9 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v2(1,1);
    const double clhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs14 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs15 =     N1[0]*clhs3 + N1[1]*clhs14;
    const double clhs16 =     Dt*v2(0,0);
    const double clhs17 =     Dt*v2(1,0);
    const double clhs18 =     clhs15*(clhs10*clhs16 + clhs13*clhs17 + clhs9) + clhs7*(clhs10*clhs8 + clhs11*clhs13);
    const double clhs19 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     clhs15*(clhs16*clhs19 + clhs17*clhs20) + clhs7*(clhs11*clhs20 + clhs19*clhs8 + clhs9);
    const double clhs22 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs23 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs24 =     clhs15*(clhs12 + clhs16*clhs22 + clhs17*clhs23) + clhs7*(clhs11*clhs23 + clhs22*clhs8);
    const double clhs25 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs26 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs27 =     clhs15*(clhs16*clhs25 + clhs17*clhs26) + clhs7*(clhs11*clhs26 + clhs12 + clhs25*clhs8);
    const double clhs28 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs29 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs12*clhs17 - clhs16*clhs9;
    const double clhs30 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs11*clhs12 - clhs8*clhs9;
    const double clhs31 =     clhs15*clhs29 + clhs30*clhs7;
    const double clhs32 =     clhs1*clhs2*clhs31;
    const double clhs33 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs34 =     clhs1*clhs3*clhs31;
    const double clhs35 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs36 =     clhs2*clhs3*clhs31;
    const double clhs37 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs38 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs39 =     -N1[0];
    const double clhs40 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs41 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs42 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs43 =     clhs15*(clhs16*clhs37 + clhs17*clhs38 + clhs39) - clhs29*(N1[0]*clhs28 + N1[1]*clhs40) - clhs30*(N1[0]*clhs41 + N1[1]*clhs42) + clhs7*(clhs11*clhs38 + clhs37*clhs8);
    const double clhs44 =     clhs1*clhs2*clhs43;
    const double clhs45 =     clhs28*clhs32 - clhs3*clhs44 + clhs33*clhs34 + clhs35*clhs36;
    const double clhs46 =     clhs1*clhs31*clhs4;
    const double clhs47 =     clhs2*clhs31*clhs4;
    const double clhs48 =     clhs32*clhs41 + clhs33*clhs46 + clhs35*clhs47 - clhs4*clhs44;
    const double clhs49 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs50 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs53 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs54 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     clhs15*(clhs16*clhs52 + clhs17*clhs53) - clhs29*(N1[0]*clhs49 + N1[1]*clhs54) - clhs30*(N1[0]*clhs55 + N1[1]*clhs56) + clhs7*(clhs11*clhs53 + clhs39 + clhs52*clhs8);
    const double clhs58 =     clhs1*clhs2*clhs57;
    const double clhs59 =     -clhs3*clhs58 + clhs32*clhs49 + clhs34*clhs50 + clhs36*clhs51;
    const double clhs60 =     clhs32*clhs55 - clhs4*clhs58 + clhs46*clhs50 + clhs47*clhs51;
    const double clhs61 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs64 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs65 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs66 =     -N1[1];
    const double clhs67 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs70 =     clhs15*(clhs16*clhs64 + clhs17*clhs65 + clhs66) - clhs29*(N1[0]*clhs61 + N1[1]*clhs67) - clhs30*(N1[0]*clhs68 + N1[1]*clhs69) + clhs7*(clhs11*clhs65 + clhs64*clhs8);
    const double clhs71 =     clhs1*clhs2*clhs70;
    const double clhs72 =     -clhs3*clhs71 + clhs32*clhs61 + clhs34*clhs62 + clhs36*clhs63;
    const double clhs73 =     clhs32*clhs68 - clhs4*clhs71 + clhs46*clhs62 + clhs47*clhs63;
    const double clhs74 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs75 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs76 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs77 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs78 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs79 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs80 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs81 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs82 =     clhs15*(clhs16*clhs77 + clhs17*clhs78) - clhs29*(N1[0]*clhs74 + N1[1]*clhs79) - clhs30*(N1[0]*clhs80 + N1[1]*clhs81) + clhs7*(clhs11*clhs78 + clhs66 + clhs77*clhs8);
    const double clhs83 =     clhs1*clhs2*clhs82;
    const double clhs84 =     -clhs3*clhs83 + clhs32*clhs74 + clhs34*clhs75 + clhs36*clhs76;
    const double clhs85 =     clhs32*clhs80 - clhs4*clhs83 + clhs46*clhs75 + clhs47*clhs76;
    const double clhs86 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs87 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs88 =     clhs0*clhs2*clhs87*(GPnormal[0]*clhs14 + GPnormal[1]*clhs6);
    const double clhs89 =     clhs2*clhs31*clhs87;
    const double clhs90 =     clhs14*clhs31*clhs87;
    const double clhs91 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs92 =     clhs14*clhs2*clhs31;
    const double clhs93 =     clhs2*clhs43*clhs87;
    const double clhs94 =     -clhs14*clhs93 + clhs33*clhs90 + clhs40*clhs89 + clhs91*clhs92;
    const double clhs95 =     clhs31*clhs6*clhs87;
    const double clhs96 =     clhs2*clhs31*clhs6;
    const double clhs97 =     clhs33*clhs95 + clhs42*clhs89 - clhs6*clhs93 + clhs91*clhs96;
    const double clhs98 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs99 =     clhs2*clhs57*clhs87;
    const double clhs100 =     -clhs14*clhs99 + clhs50*clhs90 + clhs54*clhs89 + clhs92*clhs98;
    const double clhs101 =     clhs50*clhs95 + clhs56*clhs89 - clhs6*clhs99 + clhs96*clhs98;
    const double clhs102 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs103 =     clhs2*clhs70*clhs87;
    const double clhs104 =     clhs102*clhs92 - clhs103*clhs14 + clhs62*clhs90 + clhs67*clhs89;
    const double clhs105 =     clhs102*clhs96 - clhs103*clhs6 + clhs62*clhs95 + clhs69*clhs89;
    const double clhs106 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs107 =     clhs2*clhs82*clhs87;
    const double clhs108 =     clhs106*clhs92 - clhs107*clhs14 + clhs75*clhs90 + clhs79*clhs89;
    const double clhs109 =     clhs106*clhs96 - clhs107*clhs6 + clhs75*clhs95 + clhs81*clhs89;
    const double clhs110 =     clhs0*clhs2*clhs87*(GPtangent1[0]*clhs14 + GPtangent1[1]*clhs6);

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
    lhs(8,0)=-clhs18*clhs5;
    lhs(8,1)=-clhs21*clhs5;
    lhs(8,2)=-clhs24*clhs5;
    lhs(8,3)=-clhs27*clhs5;
    lhs(8,4)=clhs0*(GPnormal[0]*clhs45 + GPnormal[1]*clhs48);
    lhs(8,5)=clhs0*(GPnormal[0]*clhs59 + GPnormal[1]*clhs60);
    lhs(8,6)=clhs0*(GPnormal[0]*clhs72 + GPnormal[1]*clhs73);
    lhs(8,7)=clhs0*(GPnormal[0]*clhs84 + GPnormal[1]*clhs85);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs18*clhs86;
    lhs(9,1)=-clhs21*clhs86;
    lhs(9,2)=-clhs24*clhs86;
    lhs(9,3)=-clhs27*clhs86;
    lhs(9,4)=clhs0*(GPtangent1[0]*clhs45 + GPtangent1[1]*clhs48);
    lhs(9,5)=clhs0*(GPtangent1[0]*clhs59 + GPtangent1[1]*clhs60);
    lhs(9,6)=clhs0*(GPtangent1[0]*clhs72 + GPtangent1[1]*clhs73);
    lhs(9,7)=clhs0*(GPtangent1[0]*clhs84 + GPtangent1[1]*clhs85);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs18*clhs88;
    lhs(10,1)=-clhs21*clhs88;
    lhs(10,2)=-clhs24*clhs88;
    lhs(10,3)=-clhs27*clhs88;
    lhs(10,4)=clhs0*(GPnormal[0]*clhs94 + GPnormal[1]*clhs97);
    lhs(10,5)=clhs0*(GPnormal[0]*clhs100 + GPnormal[1]*clhs101);
    lhs(10,6)=clhs0*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    lhs(10,7)=clhs0*(GPnormal[0]*clhs108 + GPnormal[1]*clhs109);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs110*clhs18;
    lhs(11,1)=-clhs110*clhs21;
    lhs(11,2)=-clhs110*clhs24;
    lhs(11,3)=-clhs110*clhs27;
    lhs(11,4)=clhs0*(GPtangent1[0]*clhs94 + GPtangent1[1]*clhs97);
    lhs(11,5)=clhs0*(GPtangent1[0]*clhs100 + GPtangent1[1]*clhs101);
    lhs(11,6)=clhs0*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    lhs(11,7)=clhs0*(GPtangent1[0]*clhs108 + GPtangent1[1]*clhs109);
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
    
    const std::vector<Matrix> DeltaAe = rContactData.DeltaAe;
    Vector DerivativePhiu00 = prod(DeltaAe[0],N1);
    Vector DerivativePhiu01 = prod(DeltaAe[1],N1);
    Vector DerivativePhiu10 = prod(DeltaAe[2],N1);
    Vector DerivativePhiu11 = prod(DeltaAe[3],N1);

//substitute_constants_derivatives_variables_slip 
    const double Dnormalmaster11u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster11u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster10u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u211 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u210 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u201 =     -(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster01u200 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) - (-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0))*(-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u211 =     0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(0.25*X2(0,1) - 0.25*X2(1,1) + 0.25*u2(0,1) - 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u210 =     (0.25*X2(0,0) - 0.25*X2(1,0) + 0.25*u2(0,0) - 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u201 =     -0.5/std::sqrt(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2)) + (-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))*(-0.25*X2(0,1) + 0.25*X2(1,1) - 0.25*u2(0,1) + 0.25*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
    const double Dnormalmaster00u200 =     (-0.25*X2(0,0) + 0.25*X2(1,0) - 0.25*u2(0,0) + 0.25*u2(1,0))*(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1))/std::pow(std::pow(-0.5*X2(0,0) + 0.5*X2(1,0) - 0.5*u2(0,0) + 0.5*u2(1,0), 2) + std::pow(-0.5*X2(0,1) + 0.5*X2(1,1) - 0.5*u2(0,1) + 0.5*u2(1,1), 2), 3.0L/2.0L);
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
    const double DPhi1u111 =     DerivativePhiu11[1];
    const double DPhi1u110 =     DerivativePhiu10[1];
    const double DPhi1u101 =     DerivativePhiu01[1];
    const double DPhi1u100 =     DerivativePhiu00[1];
    const double DPhi0u111 =     DerivativePhiu11[0];
    const double DPhi0u110 =     DerivativePhiu10[0];
    const double DPhi0u101 =     DerivativePhiu01[0];
    const double DPhi0u100 =     DerivativePhiu00[0];
    const double DN21u211 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u210 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u201 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u200 =     -((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN21u111 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u110 =     -(-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u101 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN21u100 =     -(-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u211 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u211*N2[0] + Dnormalmaster10u211*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u211*N2[0] + Dnormalmaster11u211*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u210 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((Dnormalmaster00u210*N2[0] + Dnormalmaster10u210*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u210*N2[0] + Dnormalmaster11u210*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - 1)/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) + (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u201 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1) + (Dnormalmaster00u201*N2[0] + Dnormalmaster10u201*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u201*N2[0] + Dnormalmaster11u201*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u200 =     ((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) - (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0) + (Dnormalmaster00u200*N2[0] + Dnormalmaster10u200*N2[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (Dnormalmaster01u200*N2[0] + Dnormalmaster11u200*N2[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1)) - (N1[0]*(X1(0,0) + u1(0,0)) + N1[0]*(X1(0,1) + u1(0,1)) + N1[1]*(X1(1,0) + u1(1,0)) + N1[1]*(X1(1,1) + u1(1,1)) - X2(1,0) - X2(1,1) - u2(1,0) - u2(1,1) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))))/std::pow(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1), 2);
    const double DN20u111 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u110 =     (-N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[1] + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u101 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DN20u100 =     (-N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) - N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + N1[0] + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))) + (N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(-(Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) - (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)))*((N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0))*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + X2(0,0) + u2(0,0)) + (N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1))*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + X2(0,1) + u2(0,1)))/std::pow((N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0))*(N2[0]*normalmaster(0,0) + N2[1]*normalmaster(1,0)) + (N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1))*(N2[0]*normalmaster(0,1) + N2[1]*normalmaster(1,1)), 2))/(X2(0,0) + X2(0,1) - X2(1,0) - X2(1,1) + u2(0,0) + u2(0,1) - u2(1,0) - u2(1,1));
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
 
    const double clhs0 =     1.0/Dt;
    const double clhs1 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs2 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs3 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs4 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs5 =     clhs0*clhs1*clhs2*(GPnormal[0]*clhs3 + GPnormal[1]*clhs4);
    const double clhs6 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs7 =     N1[0]*clhs4 + N1[1]*clhs6;
    const double clhs8 =     Dt*v2(0,1);
    const double clhs9 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs10 =     DN20u200; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs11 =     Dt*v2(1,1);
    const double clhs12 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs13 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs14 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs15 =     N1[0]*clhs3 + N1[1]*clhs14;
    const double clhs16 =     Dt*v2(0,0);
    const double clhs17 =     Dt*v2(1,0);
    const double clhs18 =     clhs15*(clhs10*clhs16 + clhs13*clhs17 + clhs9) + clhs7*(clhs10*clhs8 + clhs11*clhs13);
    const double clhs19 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs20 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs21 =     clhs15*(clhs16*clhs19 + clhs17*clhs20) + clhs7*(clhs11*clhs20 + clhs19*clhs8 + clhs9);
    const double clhs22 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs23 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs24 =     clhs15*(clhs12 + clhs16*clhs22 + clhs17*clhs23) + clhs7*(clhs11*clhs23 + clhs22*clhs8);
    const double clhs25 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs26 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs27 =     clhs15*(clhs16*clhs25 + clhs17*clhs26) + clhs7*(clhs11*clhs26 + clhs12 + clhs25*clhs8);
    const double clhs28 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs29 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs12*clhs17 - clhs16*clhs9;
    const double clhs30 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs11*clhs12 - clhs8*clhs9;
    const double clhs31 =     clhs15*clhs29 + clhs30*clhs7;
    const double clhs32 =     clhs1*clhs2*clhs31;
    const double clhs33 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs34 =     clhs1*clhs3*clhs31;
    const double clhs35 =     DPhi0u100; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs36 =     clhs2*clhs3*clhs31;
    const double clhs37 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs38 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs39 =     -N1[0];
    const double clhs40 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs41 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs42 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs43 =     clhs15*(clhs16*clhs37 + clhs17*clhs38 + clhs39) - clhs29*(N1[0]*clhs28 + N1[1]*clhs40) - clhs30*(N1[0]*clhs41 + N1[1]*clhs42) + clhs7*(clhs11*clhs38 + clhs37*clhs8);
    const double clhs44 =     clhs1*clhs2*clhs43;
    const double clhs45 =     clhs28*clhs32 - clhs3*clhs44 + clhs33*clhs34 + clhs35*clhs36;
    const double clhs46 =     clhs1*clhs31*clhs4;
    const double clhs47 =     clhs2*clhs31*clhs4;
    const double clhs48 =     clhs32*clhs41 + clhs33*clhs46 + clhs35*clhs47 - clhs4*clhs44;
    const double clhs49 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs50 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs51 =     DPhi0u101; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs52 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs53 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs54 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs55 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs56 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs57 =     clhs15*(clhs16*clhs52 + clhs17*clhs53) - clhs29*(N1[0]*clhs49 + N1[1]*clhs54) - clhs30*(N1[0]*clhs55 + N1[1]*clhs56) + clhs7*(clhs11*clhs53 + clhs39 + clhs52*clhs8);
    const double clhs58 =     clhs1*clhs2*clhs57;
    const double clhs59 =     -clhs3*clhs58 + clhs32*clhs49 + clhs34*clhs50 + clhs36*clhs51;
    const double clhs60 =     clhs32*clhs55 - clhs4*clhs58 + clhs46*clhs50 + clhs47*clhs51;
    const double clhs61 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs62 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs63 =     DPhi0u110; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs64 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs65 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs66 =     -N1[1];
    const double clhs67 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs68 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs69 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs70 =     clhs15*(clhs16*clhs64 + clhs17*clhs65 + clhs66) - clhs29*(N1[0]*clhs61 + N1[1]*clhs67) - clhs30*(N1[0]*clhs68 + N1[1]*clhs69) + clhs7*(clhs11*clhs65 + clhs64*clhs8);
    const double clhs71 =     clhs1*clhs2*clhs70;
    const double clhs72 =     -clhs3*clhs71 + clhs32*clhs61 + clhs34*clhs62 + clhs36*clhs63;
    const double clhs73 =     clhs32*clhs68 - clhs4*clhs71 + clhs46*clhs62 + clhs47*clhs63;
    const double clhs74 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs75 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs76 =     DPhi0u111; // DERIVATIVE(PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs77 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs78 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs79 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs80 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs81 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs82 =     clhs15*(clhs16*clhs77 + clhs17*clhs78) - clhs29*(N1[0]*clhs74 + N1[1]*clhs79) - clhs30*(N1[0]*clhs80 + N1[1]*clhs81) + clhs7*(clhs11*clhs78 + clhs66 + clhs77*clhs8);
    const double clhs83 =     clhs1*clhs2*clhs82;
    const double clhs84 =     -clhs3*clhs83 + clhs32*clhs74 + clhs34*clhs75 + clhs36*clhs76;
    const double clhs85 =     clhs32*clhs80 - clhs4*clhs83 + clhs46*clhs75 + clhs47*clhs76;
    const double clhs86 =     clhs0*clhs1*clhs2*(GPtangent1[0]*clhs3 + GPtangent1[1]*clhs4);
    const double clhs87 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs88 =     clhs0*clhs2*clhs87*(GPnormal[0]*clhs14 + GPnormal[1]*clhs6);
    const double clhs89 =     clhs2*clhs31*clhs87;
    const double clhs90 =     clhs14*clhs31*clhs87;
    const double clhs91 =     DPhi1u100; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs92 =     clhs14*clhs2*clhs31;
    const double clhs93 =     clhs2*clhs43*clhs87;
    const double clhs94 =     -clhs14*clhs93 + clhs33*clhs90 + clhs40*clhs89 + clhs91*clhs92;
    const double clhs95 =     clhs31*clhs6*clhs87;
    const double clhs96 =     clhs2*clhs31*clhs6;
    const double clhs97 =     clhs33*clhs95 + clhs42*clhs89 - clhs6*clhs93 + clhs91*clhs96;
    const double clhs98 =     DPhi1u101; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs99 =     clhs2*clhs57*clhs87;
    const double clhs100 =     -clhs14*clhs99 + clhs50*clhs90 + clhs54*clhs89 + clhs92*clhs98;
    const double clhs101 =     clhs50*clhs95 + clhs56*clhs89 - clhs6*clhs99 + clhs96*clhs98;
    const double clhs102 =     DPhi1u110; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs103 =     clhs2*clhs70*clhs87;
    const double clhs104 =     clhs102*clhs92 - clhs103*clhs14 + clhs62*clhs90 + clhs67*clhs89;
    const double clhs105 =     clhs102*clhs96 - clhs103*clhs6 + clhs62*clhs95 + clhs69*clhs89;
    const double clhs106 =     DPhi1u111; // DERIVATIVE(PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs107 =     clhs2*clhs82*clhs87;
    const double clhs108 =     clhs106*clhs92 - clhs107*clhs14 + clhs75*clhs90 + clhs79*clhs89;
    const double clhs109 =     clhs106*clhs96 - clhs107*clhs6 + clhs75*clhs95 + clhs81*clhs89;
    const double clhs110 =     clhs0*clhs2*clhs87*(GPtangent1[0]*clhs14 + GPtangent1[1]*clhs6);

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
    lhs(8,0)=-clhs18*clhs5;
    lhs(8,1)=-clhs21*clhs5;
    lhs(8,2)=-clhs24*clhs5;
    lhs(8,3)=-clhs27*clhs5;
    lhs(8,4)=clhs0*(GPnormal[0]*clhs45 + GPnormal[1]*clhs48);
    lhs(8,5)=clhs0*(GPnormal[0]*clhs59 + GPnormal[1]*clhs60);
    lhs(8,6)=clhs0*(GPnormal[0]*clhs72 + GPnormal[1]*clhs73);
    lhs(8,7)=clhs0*(GPnormal[0]*clhs84 + GPnormal[1]*clhs85);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=-clhs18*clhs86;
    lhs(9,1)=-clhs21*clhs86;
    lhs(9,2)=-clhs24*clhs86;
    lhs(9,3)=-clhs27*clhs86;
    lhs(9,4)=clhs0*(GPtangent1[0]*clhs45 + GPtangent1[1]*clhs48);
    lhs(9,5)=clhs0*(GPtangent1[0]*clhs59 + GPtangent1[1]*clhs60);
    lhs(9,6)=clhs0*(GPtangent1[0]*clhs72 + GPtangent1[1]*clhs73);
    lhs(9,7)=clhs0*(GPtangent1[0]*clhs84 + GPtangent1[1]*clhs85);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=-clhs18*clhs88;
    lhs(10,1)=-clhs21*clhs88;
    lhs(10,2)=-clhs24*clhs88;
    lhs(10,3)=-clhs27*clhs88;
    lhs(10,4)=clhs0*(GPnormal[0]*clhs94 + GPnormal[1]*clhs97);
    lhs(10,5)=clhs0*(GPnormal[0]*clhs100 + GPnormal[1]*clhs101);
    lhs(10,6)=clhs0*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    lhs(10,7)=clhs0*(GPnormal[0]*clhs108 + GPnormal[1]*clhs109);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=-clhs110*clhs18;
    lhs(11,1)=-clhs110*clhs21;
    lhs(11,2)=-clhs110*clhs24;
    lhs(11,3)=-clhs110*clhs27;
    lhs(11,4)=clhs0*(GPtangent1[0]*clhs94 + GPtangent1[1]*clhs97);
    lhs(11,5)=clhs0*(GPtangent1[0]*clhs100 + GPtangent1[1]*clhs101);
    lhs(11,6)=clhs0*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    lhs(11,7)=clhs0*(GPtangent1[0]*clhs108 + GPtangent1[1]*clhs109);
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
    
//     const Matrix normalmaster    = rContactData.NormalsMaster;
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
//     const Matrix X1 = rContactData.X1;
//     const Matrix X2 = rContactData.X2;
//     const Matrix u1 = rContactData.u1;
//     const Matrix u2 = rContactData.u2;
//     const Matrix v1 = rContactData.v1;
//     const Matrix v2 = rContactData.v2;
//     
//     const std::vector<Matrix> DeltaAe = rContactData.DeltaAe;
//     Vector DerivativePhiu00 = prod(DeltaAe[0],N1);
//     Vector DerivativePhiu01 = prod(DeltaAe[1],N1);
//     Vector DerivativePhiu10 = prod(DeltaAe[2],N1);
//     Vector DerivativePhiu11 = prod(DeltaAe[3],N1);
// 
//substitute_constants_derivatives_variables_inactive 
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
    const double crhs2 =     Phi[0]; // PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs3 =     Phi[1]; // PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs4 =     crhs2*lm(0,0) + crhs3*lm(1,0);
    const double crhs5 =     crhs1*crhs4;
    const double crhs6 =     crhs2*lm(0,1) + crhs3*lm(1,1);
    const double crhs7 =     crhs1*crhs6;
    const double crhs8 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     N1[0]*crhs1;
    const double crhs10 =     N1[1]*crhs1;
    const double crhs11 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     integration_point_gap; // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs14 =     crhs1*crhs13*crhs2;
    const double crhs15 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs16 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs17 =     crhs1*crhs13*crhs3;

    rhs[0]=-crhs0*crhs5;
    rhs[1]=-crhs0*crhs7;
    rhs[2]=-crhs5*crhs8;
    rhs[3]=-crhs7*crhs8;
    rhs[4]=crhs4*crhs9;
    rhs[5]=crhs6*crhs9;
    rhs[6]=crhs10*crhs4;
    rhs[7]=crhs10*crhs6;
    rhs[8]=crhs14*(GPnormal[0]*crhs11 + GPnormal[1]*crhs12);
    rhs[9]=crhs14*(GPtangent1[0]*crhs11 + GPtangent1[1]*crhs12);
    rhs[10]=crhs17*(GPnormal[0]*crhs15 + GPnormal[1]*crhs16);
    rhs[11]=crhs17*(GPtangent1[0]*crhs15 + GPtangent1[1]*crhs16);

    
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
    const Matrix lm              = rContactData.LagrangeMultipliers;
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
    const double crhs9 =     crhs2*crhs3*crhs8*Phi[0]; // CRHS2*CRHS3*CRHS8*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     crhs2*crhs3*crhs8*Phi[1]; // CRHS2*CRHS3*CRHS8*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

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
    const Matrix lm              = rContactData.LagrangeMultipliers;
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
    const double crhs9 =     crhs2*crhs3*crhs8*Phi[0]; // CRHS2*CRHS3*CRHS8*PHI[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     crhs2*crhs3*crhs8*Phi[1]; // CRHS2*CRHS3*CRHS8*PHI[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1))

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
//     
//substitute_inactive_rhs
    
    return rhs;
}
    /***********************************************************************************/
    /***********************************************************************************/
    
    static inline bounded_matrix<double,4,4> ComputeDeltaDe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,4,4> DeltaDe;
    
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
    
    static inline bounded_matrix<double,4,4> ComputeDeltaMe(
        const Vector N1, 
        const ContactData& rContactData,
        const unsigned int derivative_index
        )
{
    bounded_matrix<double,4,4> DeltaMe;
    
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
};// class Contact2D2N2N
}
#endif /* KRATOS_CONTACT2D2N2N defined */
