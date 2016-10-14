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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix X1 = GetCoordinates(rContactData.SlaveGeometry, false);
    const Matrix X2 = GetCoordinates(rContactData.MasterGeometry, false);
    const Matrix u1 = GetVariableMatrix(rContactData.SlaveGeometry, DISPLACEMENT, 0);
    const Matrix u2 = GetVariableMatrix(rContactData.MasterGeometry, DISPLACEMENT, 0);
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);

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
    const double clhs4 =     clhs2*clhs3;
    const double clhs5 =     DN20u201; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs6 =     DN20u210; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs7 =     DN20u211; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs8 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs9 =     DN20u100; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs10 =     clhs0*clhs8 + clhs2*clhs9;
    const double clhs11 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs12 =     DN20u101; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs13 =     clhs0*clhs11 + clhs12*clhs2;
    const double clhs14 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs15 =     DN20u110; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs16 =     clhs0*clhs14 + clhs15*clhs2;
    const double clhs17 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs18 =     DN20u111; // DERIVATIVE(N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs19 =     clhs0*clhs17 + clhs18*clhs2;
    const double clhs20 =     Phi[0]*clhs2;
    const double clhs21 =     -clhs0*clhs20;
    const double clhs22 =     Phi[1]*clhs2;
    const double clhs23 =     -clhs0*clhs22;
    const double clhs24 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double clhs25 =     clhs2*clhs24;
    const double clhs26 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs27 =     DN21u200; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs28 =     DN21u201; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs29 =     DN21u210; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs30 =     DN21u211; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs31 =     DN21u100; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs32 =     clhs2*clhs31 + clhs26*clhs8;
    const double clhs33 =     DN21u101; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs34 =     clhs11*clhs26 + clhs2*clhs33;
    const double clhs35 =     DN21u110; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs36 =     clhs14*clhs26 + clhs2*clhs35;
    const double clhs37 =     DN21u111; // DERIVATIVE(N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs38 =     clhs17*clhs26 + clhs2*clhs37;
    const double clhs39 =     -clhs20*clhs26;
    const double clhs40 =     -clhs22*clhs26;
    const double clhs41 =     N1[0]*clhs3;
    const double clhs42 =     N1[0]*clhs2;
    const double clhs43 =     Phi[0]*clhs42;
    const double clhs44 =     Phi[1]*clhs42;
    const double clhs45 =     N1[0]*clhs24;
    const double clhs46 =     N1[1]*clhs3;
    const double clhs47 =     N1[1]*clhs2;
    const double clhs48 =     Phi[0]*clhs47;
    const double clhs49 =     Phi[1]*clhs47;
    const double clhs50 =     N1[1]*clhs24;
    const double clhs51 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs52 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs53 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs54 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs55 =     1.0/Dt;
    const double clhs56 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs57 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs58 =     N1[0]*clhs56 + N1[1]*clhs57;
    const double clhs59 =     Dt*v2(0,1);
    const double clhs60 =     Dt*v2(1,1);
    const double clhs61 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs62 =     N1[0]*clhs54 + N1[1]*clhs61;
    const double clhs63 =     Dt*v2(0,0);
    const double clhs64 =     Dt*v2(1,0);
    const double clhs65 =     clhs55*(clhs58*(clhs1*clhs59 + clhs27*clhs60) + clhs62*(clhs0 + clhs1*clhs63 + clhs27*clhs64));
    const double clhs66 =     clhs51*clhs53 + clhs54*clhs65;
    const double clhs67 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs68 =     clhs53*clhs67 + clhs56*clhs65;
    const double clhs69 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs70 =     clhs55*(clhs58*(clhs0 + clhs28*clhs60 + clhs5*clhs59) + clhs62*(clhs28*clhs64 + clhs5*clhs63));
    const double clhs71 =     clhs51*clhs69 + clhs54*clhs70;
    const double clhs72 =     clhs56*clhs70 + clhs67*clhs69;
    const double clhs73 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs74 =     clhs55*(clhs58*(clhs29*clhs60 + clhs59*clhs6) + clhs62*(clhs26 + clhs29*clhs64 + clhs6*clhs63));
    const double clhs75 =     clhs51*clhs73 + clhs54*clhs74;
    const double clhs76 =     clhs56*clhs74 + clhs67*clhs73;
    const double clhs77 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs78 =     clhs55*(clhs58*(clhs26 + clhs30*clhs60 + clhs59*clhs7) + clhs62*(clhs30*clhs64 + clhs63*clhs7));
    const double clhs79 =     clhs51*clhs77 + clhs54*clhs78;
    const double clhs80 =     clhs56*clhs78 + clhs67*clhs77;
    const double clhs81 =     clhs2*clhs52;
    const double clhs82 =     clhs2*clhs51;
    const double clhs83 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs84 =     clhs51*clhs52;
    const double clhs85 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs86 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - clhs0*clhs63 - clhs26*clhs64;
    const double clhs87 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - clhs0*clhs59 - clhs26*clhs60;
    const double clhs88 =     clhs58*clhs87 + clhs62*clhs86;
    const double clhs89 =     clhs2*clhs55*clhs88;
    const double clhs90 =     clhs54*clhs55*clhs88;
    const double clhs91 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs92 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs93 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs94 =     clhs2*clhs55*(-clhs58*(clhs31*clhs60 + clhs59*clhs9) + clhs62*(N1[0] - clhs31*clhs64 - clhs63*clhs9) + clhs86*(N1[0]*clhs85 + N1[1]*clhs91) + clhs87*(N1[0]*clhs92 + N1[1]*clhs93));
    const double clhs95 =     -clhs8*clhs84 + clhs8*clhs90 - clhs81*Dnormalslave00u100 - clhs82*clhs83 + clhs85*clhs89 + clhs94*tan1slave(0,0); // -CLHS8*CLHS84 + CLHS8*CLHS90 - CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS82*CLHS83 + CLHS85*CLHS89 + CLHS94*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs96 =     clhs2*clhs67;
    const double clhs97 =     clhs52*clhs67;
    const double clhs98 =     clhs55*clhs56*clhs88;
    const double clhs99 =     -clhs8*clhs97 + clhs8*clhs98 - clhs81*Dnormalslave01u100 - clhs83*clhs96 + clhs89*clhs92 + clhs94*tan1slave(0,1); // -CLHS8*CLHS97 + CLHS8*CLHS98 - CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS83*CLHS96 + CLHS89*CLHS92 + CLHS94*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs100 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs101 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs102 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs103 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs104 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs105 =     clhs2*clhs55*(clhs58*(N1[0] - clhs12*clhs59 - clhs33*clhs60) - clhs62*(clhs12*clhs63 + clhs33*clhs64) + clhs86*(N1[0]*clhs101 + N1[1]*clhs102) + clhs87*(N1[0]*clhs103 + N1[1]*clhs104));
    const double clhs106 =     -clhs100*clhs82 + clhs101*clhs89 + clhs105*tan1slave(0,0) - clhs11*clhs84 + clhs11*clhs90 - clhs81*Dnormalslave00u101; // -CLHS100*CLHS82 + CLHS101*CLHS89 + CLHS105*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS11*CLHS84 + CLHS11*CLHS90 - CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs107 =     -clhs100*clhs96 + clhs103*clhs89 + clhs105*tan1slave(0,1) - clhs11*clhs97 + clhs11*clhs98 - clhs81*Dnormalslave01u101; // -CLHS100*CLHS96 + CLHS103*CLHS89 + CLHS105*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS11*CLHS97 + CLHS11*CLHS98 - CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs108 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs109 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs110 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs111 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs112 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs113 =     clhs2*clhs55*(-clhs58*(clhs15*clhs59 + clhs35*clhs60) + clhs62*(N1[1] - clhs15*clhs63 - clhs35*clhs64) + clhs86*(N1[0]*clhs109 + N1[1]*clhs110) + clhs87*(N1[0]*clhs111 + N1[1]*clhs112));
    const double clhs114 =     -clhs108*clhs82 + clhs109*clhs89 + clhs113*tan1slave(0,0) - clhs14*clhs84 + clhs14*clhs90 - clhs81*Dnormalslave00u110; // -CLHS108*CLHS82 + CLHS109*CLHS89 + CLHS113*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS14*CLHS84 + CLHS14*CLHS90 - CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs115 =     -clhs108*clhs96 + clhs111*clhs89 + clhs113*tan1slave(0,1) - clhs14*clhs97 + clhs14*clhs98 - clhs81*Dnormalslave01u110; // -CLHS108*CLHS96 + CLHS111*CLHS89 + CLHS113*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS14*CLHS97 + CLHS14*CLHS98 - CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs116 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs117 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs118 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs119 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs120 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs121 =     clhs2*clhs55*(clhs58*(N1[1] - clhs18*clhs59 - clhs37*clhs60) - clhs62*(clhs18*clhs63 + clhs37*clhs64) + clhs86*(N1[0]*clhs117 + N1[1]*clhs118) + clhs87*(N1[0]*clhs119 + N1[1]*clhs120));
    const double clhs122 =     -clhs116*clhs82 + clhs117*clhs89 + clhs121*tan1slave(0,0) - clhs17*clhs84 + clhs17*clhs90 - clhs81*Dnormalslave00u111; // -CLHS116*CLHS82 + CLHS117*CLHS89 + CLHS121*TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS17*CLHS84 + CLHS17*CLHS90 - CLHS81*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs123 =     -clhs116*clhs96 + clhs119*clhs89 + clhs121*tan1slave(0,1) - clhs17*clhs97 + clhs17*clhs98 - clhs81*Dnormalslave01u111; // -CLHS116*CLHS96 + CLHS119*CLHS89 + CLHS121*TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS17*CLHS97 + CLHS17*CLHS98 - CLHS81*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs124 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs125 =     clhs124*clhs53 + clhs61*clhs65;
    const double clhs126 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs127 =     clhs126*clhs53 + clhs57*clhs65;
    const double clhs128 =     clhs124*clhs69 + clhs61*clhs70;
    const double clhs129 =     clhs126*clhs69 + clhs57*clhs70;
    const double clhs130 =     clhs124*clhs73 + clhs61*clhs74;
    const double clhs131 =     clhs126*clhs73 + clhs57*clhs74;
    const double clhs132 =     clhs124*clhs77 + clhs61*clhs78;
    const double clhs133 =     clhs126*clhs77 + clhs57*clhs78;
    const double clhs134 =     clhs124*clhs2;
    const double clhs135 =     clhs124*clhs52;
    const double clhs136 =     clhs55*clhs61*clhs88;
    const double clhs137 =     -clhs134*clhs83 - clhs135*clhs8 + clhs136*clhs8 - clhs81*Dnormalslave10u100 + clhs89*clhs91 + clhs94*tan1slave(1,0); // -CLHS134*CLHS83 - CLHS135*CLHS8 + CLHS136*CLHS8 - CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS89*CLHS91 + CLHS94*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs138 =     clhs126*clhs2;
    const double clhs139 =     clhs126*clhs52;
    const double clhs140 =     clhs55*clhs57*clhs88;
    const double clhs141 =     -clhs138*clhs83 - clhs139*clhs8 + clhs140*clhs8 - clhs81*Dnormalslave11u100 + clhs89*clhs93 + clhs94*tan1slave(1,1); // -CLHS138*CLHS83 - CLHS139*CLHS8 + CLHS140*CLHS8 - CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS89*CLHS93 + CLHS94*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs142 =     -clhs100*clhs134 + clhs102*clhs89 + clhs105*tan1slave(1,0) - clhs11*clhs135 + clhs11*clhs136 - clhs81*Dnormalslave10u101; // -CLHS100*CLHS134 + CLHS102*CLHS89 + CLHS105*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS11*CLHS135 + CLHS11*CLHS136 - CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs143 =     -clhs100*clhs138 + clhs104*clhs89 + clhs105*tan1slave(1,1) - clhs11*clhs139 + clhs11*clhs140 - clhs81*Dnormalslave11u101; // -CLHS100*CLHS138 + CLHS104*CLHS89 + CLHS105*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS11*CLHS139 + CLHS11*CLHS140 - CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs144 =     -clhs108*clhs134 + clhs110*clhs89 + clhs113*tan1slave(1,0) - clhs135*clhs14 + clhs136*clhs14 - clhs81*Dnormalslave10u110; // -CLHS108*CLHS134 + CLHS110*CLHS89 + CLHS113*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS135*CLHS14 + CLHS136*CLHS14 - CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs145 =     -clhs108*clhs138 + clhs112*clhs89 + clhs113*tan1slave(1,1) - clhs139*clhs14 + clhs14*clhs140 - clhs81*Dnormalslave11u110; // -CLHS108*CLHS138 + CLHS112*CLHS89 + CLHS113*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS139*CLHS14 + CLHS14*CLHS140 - CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs146 =     -clhs116*clhs134 + clhs118*clhs89 + clhs121*tan1slave(1,0) - clhs135*clhs17 + clhs136*clhs17 - clhs81*Dnormalslave10u111; // -CLHS116*CLHS134 + CLHS118*CLHS89 + CLHS121*TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS135*CLHS17 + CLHS136*CLHS17 - CLHS81*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs147 =     -clhs116*clhs138 + clhs120*clhs89 + clhs121*tan1slave(1,1) - clhs139*clhs17 + clhs140*clhs17 - clhs81*Dnormalslave11u111; // -CLHS116*CLHS138 + CLHS120*CLHS89 + CLHS121*TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) - CLHS139*CLHS17 + CLHS140*CLHS17 - CLHS81*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))

    lhs(0,0)=-clhs1*clhs4;
    lhs(0,1)=-clhs4*clhs5;
    lhs(0,2)=-clhs4*clhs6;
    lhs(0,3)=-clhs4*clhs7;
    lhs(0,4)=-clhs10*clhs3;
    lhs(0,5)=-clhs13*clhs3;
    lhs(0,6)=-clhs16*clhs3;
    lhs(0,7)=-clhs19*clhs3;
    lhs(0,8)=clhs21;
    lhs(0,9)=0;
    lhs(0,10)=clhs23;
    lhs(0,11)=0;
    lhs(1,0)=-clhs1*clhs25;
    lhs(1,1)=-clhs25*clhs5;
    lhs(1,2)=-clhs25*clhs6;
    lhs(1,3)=-clhs25*clhs7;
    lhs(1,4)=-clhs10*clhs24;
    lhs(1,5)=-clhs13*clhs24;
    lhs(1,6)=-clhs16*clhs24;
    lhs(1,7)=-clhs19*clhs24;
    lhs(1,8)=0;
    lhs(1,9)=clhs21;
    lhs(1,10)=0;
    lhs(1,11)=clhs23;
    lhs(2,0)=-clhs27*clhs4;
    lhs(2,1)=-clhs28*clhs4;
    lhs(2,2)=-clhs29*clhs4;
    lhs(2,3)=-clhs30*clhs4;
    lhs(2,4)=-clhs3*clhs32;
    lhs(2,5)=-clhs3*clhs34;
    lhs(2,6)=-clhs3*clhs36;
    lhs(2,7)=-clhs3*clhs38;
    lhs(2,8)=clhs39;
    lhs(2,9)=0;
    lhs(2,10)=clhs40;
    lhs(2,11)=0;
    lhs(3,0)=-clhs25*clhs27;
    lhs(3,1)=-clhs25*clhs28;
    lhs(3,2)=-clhs25*clhs29;
    lhs(3,3)=-clhs25*clhs30;
    lhs(3,4)=-clhs24*clhs32;
    lhs(3,5)=-clhs24*clhs34;
    lhs(3,6)=-clhs24*clhs36;
    lhs(3,7)=-clhs24*clhs38;
    lhs(3,8)=0;
    lhs(3,9)=clhs39;
    lhs(3,10)=0;
    lhs(3,11)=clhs40;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=clhs41*clhs8;
    lhs(4,5)=clhs11*clhs41;
    lhs(4,6)=clhs14*clhs41;
    lhs(4,7)=clhs17*clhs41;
    lhs(4,8)=clhs43;
    lhs(4,9)=0;
    lhs(4,10)=clhs44;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=clhs45*clhs8;
    lhs(5,5)=clhs11*clhs45;
    lhs(5,6)=clhs14*clhs45;
    lhs(5,7)=clhs17*clhs45;
    lhs(5,8)=0;
    lhs(5,9)=clhs43;
    lhs(5,10)=0;
    lhs(5,11)=clhs44;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=clhs46*clhs8;
    lhs(6,5)=clhs11*clhs46;
    lhs(6,6)=clhs14*clhs46;
    lhs(6,7)=clhs17*clhs46;
    lhs(6,8)=clhs48;
    lhs(6,9)=0;
    lhs(6,10)=clhs49;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=clhs50*clhs8;
    lhs(7,5)=clhs11*clhs50;
    lhs(7,6)=clhs14*clhs50;
    lhs(7,7)=clhs17*clhs50;
    lhs(7,8)=0;
    lhs(7,9)=clhs48;
    lhs(7,10)=0;
    lhs(7,11)=clhs49;
    lhs(8,0)=clhs20*(GPnormal[0]*clhs66 + GPnormal[1]*clhs68);
    lhs(8,1)=clhs20*(GPnormal[0]*clhs71 + GPnormal[1]*clhs72);
    lhs(8,2)=clhs20*(GPnormal[0]*clhs75 + GPnormal[1]*clhs76);
    lhs(8,3)=clhs20*(GPnormal[0]*clhs79 + GPnormal[1]*clhs80);
    lhs(8,4)=-Phi[0]*(GPnormal[0]*clhs95 + GPnormal[1]*clhs99);
    lhs(8,5)=-Phi[0]*(GPnormal[0]*clhs106 + GPnormal[1]*clhs107);
    lhs(8,6)=-Phi[0]*(GPnormal[0]*clhs114 + GPnormal[1]*clhs115);
    lhs(8,7)=-Phi[0]*(GPnormal[0]*clhs122 + GPnormal[1]*clhs123);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs20*(GPtangent1[0]*clhs66 + GPtangent1[1]*clhs68);
    lhs(9,1)=clhs20*(GPtangent1[0]*clhs71 + GPtangent1[1]*clhs72);
    lhs(9,2)=clhs20*(GPtangent1[0]*clhs75 + GPtangent1[1]*clhs76);
    lhs(9,3)=clhs20*(GPtangent1[0]*clhs79 + GPtangent1[1]*clhs80);
    lhs(9,4)=-Phi[0]*(GPtangent1[0]*clhs95 + GPtangent1[1]*clhs99);
    lhs(9,5)=-Phi[0]*(GPtangent1[0]*clhs106 + GPtangent1[1]*clhs107);
    lhs(9,6)=-Phi[0]*(GPtangent1[0]*clhs114 + GPtangent1[1]*clhs115);
    lhs(9,7)=-Phi[0]*(GPtangent1[0]*clhs122 + GPtangent1[1]*clhs123);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs22*(GPnormal[0]*clhs125 + GPnormal[1]*clhs127);
    lhs(10,1)=clhs22*(GPnormal[0]*clhs128 + GPnormal[1]*clhs129);
    lhs(10,2)=clhs22*(GPnormal[0]*clhs130 + GPnormal[1]*clhs131);
    lhs(10,3)=clhs22*(GPnormal[0]*clhs132 + GPnormal[1]*clhs133);
    lhs(10,4)=-Phi[1]*(GPnormal[0]*clhs137 + GPnormal[1]*clhs141);
    lhs(10,5)=-Phi[1]*(GPnormal[0]*clhs142 + GPnormal[1]*clhs143);
    lhs(10,6)=-Phi[1]*(GPnormal[0]*clhs144 + GPnormal[1]*clhs145);
    lhs(10,7)=-Phi[1]*(GPnormal[0]*clhs146 + GPnormal[1]*clhs147);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs22*(GPtangent1[0]*clhs125 + GPtangent1[1]*clhs127);
    lhs(11,1)=clhs22*(GPtangent1[0]*clhs128 + GPtangent1[1]*clhs129);
    lhs(11,2)=clhs22*(GPtangent1[0]*clhs130 + GPtangent1[1]*clhs131);
    lhs(11,3)=clhs22*(GPtangent1[0]*clhs132 + GPtangent1[1]*clhs133);
    lhs(11,4)=-Phi[1]*(GPtangent1[0]*clhs137 + GPtangent1[1]*clhs141);
    lhs(11,5)=-Phi[1]*(GPtangent1[0]*clhs142 + GPtangent1[1]*clhs143);
    lhs(11,6)=-Phi[1]*(GPtangent1[0]*clhs144 + GPtangent1[1]*clhs145);
    lhs(11,7)=-Phi[1]*(GPtangent1[0]*clhs146 + GPtangent1[1]*clhs147);
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
    
    const Vector GPnormal     = prod(trans(normalslave), N1);
    const Vector GPtangent1   = prod(trans(tan1slave), N1);
    
    const Matrix v1 = GetVariableMatrix(rContactData.SlaveGeometry, VELOCITY, 0); 
    const Matrix v2 = GetVariableMatrix(rContactData.MasterGeometry, VELOCITY, 0);
    
    const double crhs0 =     N2[0]; // N2[0](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs1 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs2 =     Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double crhs3 =     crhs1*crhs2;
    const double crhs4 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double crhs5 =     crhs1*crhs4;
    const double crhs6 =     N2[1]; // N2[1](U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs7 =     N1[0]*crhs1;
    const double crhs8 =     N1[1]*crhs1;
    const double crhs9 =     Phi[0]*crhs1;
    const double crhs10 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs11 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs14 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     ((N1[0]*crhs11 + N1[1]*crhs12)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - crhs0*(Dt*v2(0,0)) - crhs6*(Dt*v2(1,0))) + (N1[0]*crhs13 + N1[1]*crhs14)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - crhs0*(Dt*v2(0,1)) - crhs6*(Dt*v2(1,1))))/Dt;
    const double crhs16 =     -crhs10*normalslave(0,0) + crhs11*crhs15; // -CRHS10*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) + CRHS11*CRHS15
    const double crhs17 =     -crhs10*normalslave(0,1) + crhs13*crhs15; // -CRHS10*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) + CRHS13*CRHS15
    const double crhs18 =     Phi[1]*crhs1;
    const double crhs19 =     -crhs10*normalslave(1,0) + crhs12*crhs15; // -CRHS10*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) + CRHS12*CRHS15
    const double crhs20 =     -crhs10*normalslave(1,1) + crhs14*crhs15; // -CRHS10*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)) + CRHS14*CRHS15

    rhs[0]=crhs0*crhs3;
    rhs[1]=crhs0*crhs5;
    rhs[2]=crhs3*crhs6;
    rhs[3]=crhs5*crhs6;
    rhs[4]=-crhs2*crhs7;
    rhs[5]=-crhs4*crhs7;
    rhs[6]=-crhs2*crhs8;
    rhs[7]=-crhs4*crhs8;
    rhs[8]=crhs9*(GPnormal[0]*crhs16 + GPnormal[1]*crhs17);
    rhs[9]=crhs9*(GPtangent1[0]*crhs16 + GPtangent1[1]*crhs17);
    rhs[10]=crhs18*(GPnormal[0]*crhs19 + GPnormal[1]*crhs20);
    rhs[11]=crhs18*(GPtangent1[0]*crhs19 + GPtangent1[1]*crhs20);

    
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