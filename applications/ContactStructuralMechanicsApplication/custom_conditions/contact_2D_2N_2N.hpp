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
    
    const Matrix NormalsMaster = rContactData.NormalsMaster;
    const Vector masternormal  = prod(NormalsMaster, N2);
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
    const double Dgapgu211 =     N2[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1));
    const double Dgapgu210 =     N2[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0));
    const double Dgapgu201 =     N2[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1));
    const double Dgapgu200 =     N2[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0));
    const double Dgapgu111 =     -N1[1]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (Dnormalslave00u111*N1[0] + Dnormalslave10u111*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u111*N1[0] + Dnormalslave11u111*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1)));
    const double Dgapgu110 =     -N1[1]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (Dnormalslave00u110*N1[0] + Dnormalslave10u110*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u110*N1[0] + Dnormalslave11u110*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1)));
    const double Dgapgu101 =     -N1[0]*(N1[0]*normalslave(0,1) + N1[1]*normalslave(1,1)) + (Dnormalslave00u101*N1[0] + Dnormalslave10u101*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u101*N1[0] + Dnormalslave11u101*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1)));
    const double Dgapgu100 =     -N1[0]*(N1[0]*normalslave(0,0) + N1[1]*normalslave(1,0)) + (Dnormalslave00u100*N1[0] + Dnormalslave10u100*N1[1])*(-N1[0]*(X1(0,0) + u1(0,0)) - N1[1]*(X1(1,0) + u1(1,0)) + N2[0]*(X2(0,0) + u2(0,0)) + N2[1]*(X2(1,0) + u2(1,0))) + (Dnormalslave01u100*N1[0] + Dnormalslave11u100*N1[1])*(-N1[0]*(X1(0,1) + u1(0,1)) - N1[1]*(X1(1,1) + u1(1,1)) + N2[0]*(X2(0,1) + u2(0,1)) + N2[1]*(X2(1,1) + u2(1,1)));
    const double DdetJu111 =     (-0.25*X1(0,1) + 0.25*X1(1,1) - 0.25*u1(0,1) + 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu110 =     (-0.25*X1(0,0) + 0.25*X1(1,0) - 0.25*u1(0,0) + 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu101 =     (0.25*X1(0,1) - 0.25*X1(1,1) + 0.25*u1(0,1) - 0.25*u1(1,1))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));
    const double DdetJu100 =     (0.25*X1(0,0) - 0.25*X1(1,0) + 0.25*u1(0,0) - 0.25*u1(1,0))/std::sqrt(0.25*std::pow(X1(0,0), 2) - 0.5*X1(0,0)*X1(1,0) + 0.5*X1(0,0)*u1(0,0) - 0.5*X1(0,0)*u1(1,0) + 0.25*std::pow(X1(0,1), 2) - 0.5*X1(0,1)*X1(1,1) + 0.5*X1(0,1)*u1(0,1) - 0.5*X1(0,1)*u1(1,1) + 0.25*std::pow(X1(1,0), 2) - 0.5*X1(1,0)*u1(0,0) + 0.5*X1(1,0)*u1(1,0) + 0.25*std::pow(X1(1,1), 2) - 0.5*X1(1,1)*u1(0,1) + 0.5*X1(1,1)*u1(1,1) + 0.25*std::pow(u1(0,0), 2) - 0.5*u1(0,0)*u1(1,0) + 0.25*std::pow(u1(0,1), 2) - 0.5*u1(0,1)*u1(1,1) + 0.25*std::pow(u1(1,0), 2) + 0.25*std::pow(u1(1,1), 2));

    const double clhs0 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs1 =     DdetJu100; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs2 =     Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double clhs3 =     N2[0]*clhs2;
    const double clhs4 =     DdetJu101; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs5 =     DdetJu110; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs6 =     DdetJu111; // DERIVATIVE(DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs7 =     N2[0]*clhs0;
    const double clhs8 =     -Phi[0]*clhs7;
    const double clhs9 =     -Phi[1]*clhs7;
    const double clhs10 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double clhs11 =     N2[0]*clhs10;
    const double clhs12 =     N2[1]*clhs2;
    const double clhs13 =     N2[1]*clhs0;
    const double clhs14 =     -Phi[0]*clhs13;
    const double clhs15 =     -Phi[1]*clhs13;
    const double clhs16 =     N2[1]*clhs10;
    const double clhs17 =     N1[0]*clhs2;
    const double clhs18 =     N1[0]*clhs0;
    const double clhs19 =     Phi[0]*clhs18;
    const double clhs20 =     Phi[1]*clhs18;
    const double clhs21 =     N1[0]*clhs10;
    const double clhs22 =     N1[1]*clhs2;
    const double clhs23 =     N1[1]*clhs0;
    const double clhs24 =     Phi[0]*clhs23;
    const double clhs25 =     Phi[1]*clhs23;
    const double clhs26 =     N1[1]*clhs10;
    const double clhs27 =     Phi[0]*clhs0;
    const double clhs28 =     normalslave(0,0); // NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs29 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double clhs30 =     Dgapgu200; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,0))
    const double clhs31 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs32 =     1.0/Dt;
    const double clhs33 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs34 =     N1[0]*clhs31 + N1[1]*clhs33;
    const double clhs35 =     N2[0]*clhs32*clhs34;
    const double clhs36 =     clhs28*clhs30 + clhs31*clhs35;
    const double clhs37 =     normalslave(0,1); // NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs38 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs39 =     clhs30*clhs37 + clhs35*clhs38;
    const double clhs40 =     Dgapgu201; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(0,1))
    const double clhs41 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs42 =     N1[0]*clhs38 + N1[1]*clhs41;
    const double clhs43 =     N2[0]*clhs32*clhs42;
    const double clhs44 =     clhs28*clhs40 + clhs31*clhs43;
    const double clhs45 =     clhs37*clhs40 + clhs38*clhs43;
    const double clhs46 =     Dgapgu210; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,0))
    const double clhs47 =     N2[1]*clhs32*clhs34;
    const double clhs48 =     clhs28*clhs46 + clhs31*clhs47;
    const double clhs49 =     clhs37*clhs46 + clhs38*clhs47;
    const double clhs50 =     Dgapgu211; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U2(1,1))
    const double clhs51 =     N2[1]*clhs32*clhs42;
    const double clhs52 =     clhs28*clhs50 + clhs31*clhs51;
    const double clhs53 =     clhs37*clhs50 + clhs38*clhs51;
    const double clhs54 =     clhs0*clhs29;
    const double clhs55 =     clhs0*clhs28;
    const double clhs56 =     Dgapgu100; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,0))
    const double clhs57 =     clhs28*clhs29;
    const double clhs58 =     Dtan1slave00u100; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs59 =     N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - N2[0]*(Dt*v2(0,0)) - N2[1]*(Dt*v2(1,0));
    const double clhs60 =     N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - N2[0]*(Dt*v2(0,1)) - N2[1]*(Dt*v2(1,1));
    const double clhs61 =     clhs34*clhs59 + clhs42*clhs60;
    const double clhs62 =     clhs0*clhs32*clhs61;
    const double clhs63 =     clhs31*clhs32*clhs61;
    const double clhs64 =     Dtan1slave10u100; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs65 =     Dtan1slave01u100; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs66 =     Dtan1slave11u100; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0))
    const double clhs67 =     clhs0*clhs32*(N1[0]*clhs34 + clhs59*(N1[0]*clhs58 + N1[1]*clhs64) + clhs60*(N1[0]*clhs65 + N1[1]*clhs66));
    const double clhs68 =     -clhs1*clhs57 + clhs1*clhs63 + tan1slave(0,0)*clhs67 - clhs54*Dnormalslave00u100 - clhs55*clhs56 + clhs58*clhs62; // -CLHS1*CLHS57 + CLHS1*CLHS63 + TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS67 - CLHS54*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS55*CLHS56 + CLHS58*CLHS62
    const double clhs69 =     clhs0*clhs37;
    const double clhs70 =     clhs29*clhs37;
    const double clhs71 =     clhs32*clhs38*clhs61;
    const double clhs72 =     -clhs1*clhs70 + clhs1*clhs71 + tan1slave(0,1)*clhs67 - clhs54*Dnormalslave01u100 - clhs56*clhs69 + clhs62*clhs65; // -CLHS1*CLHS70 + CLHS1*CLHS71 + TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS67 - CLHS54*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) - CLHS56*CLHS69 + CLHS62*CLHS65
    const double clhs73 =     Dgapgu101; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(0,1))
    const double clhs74 =     Dtan1slave00u101; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs75 =     Dtan1slave10u101; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs76 =     Dtan1slave01u101; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs77 =     Dtan1slave11u101; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1))
    const double clhs78 =     clhs0*clhs32*(N1[0]*clhs42 + clhs59*(N1[0]*clhs74 + N1[1]*clhs75) + clhs60*(N1[0]*clhs76 + N1[1]*clhs77));
    const double clhs79 =     tan1slave(0,0)*clhs78 - clhs4*clhs57 + clhs4*clhs63 - clhs54*Dnormalslave00u101 - clhs55*clhs73 + clhs62*clhs74; // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS78 - CLHS4*CLHS57 + CLHS4*CLHS63 - CLHS54*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) - CLHS55*CLHS73 + CLHS62*CLHS74
    const double clhs80 =     tan1slave(0,1)*clhs78 - clhs4*clhs70 + clhs4*clhs71 - clhs54*Dnormalslave01u101 + clhs62*clhs76 - clhs69*clhs73; // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS78 - CLHS4*CLHS70 + CLHS4*CLHS71 - CLHS54*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS62*CLHS76 - CLHS69*CLHS73
    const double clhs81 =     Dgapgu110; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,0))
    const double clhs82 =     Dtan1slave00u110; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs83 =     Dtan1slave10u110; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs84 =     Dtan1slave01u110; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs85 =     Dtan1slave11u110; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0))
    const double clhs86 =     clhs0*clhs32*(N1[1]*clhs34 + clhs59*(N1[0]*clhs82 + N1[1]*clhs83) + clhs60*(N1[0]*clhs84 + N1[1]*clhs85));
    const double clhs87 =     tan1slave(0,0)*clhs86 - clhs5*clhs57 + clhs5*clhs63 - clhs54*Dnormalslave00u110 - clhs55*clhs81 + clhs62*clhs82; // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS86 - CLHS5*CLHS57 + CLHS5*CLHS63 - CLHS54*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) - CLHS55*CLHS81 + CLHS62*CLHS82
    const double clhs88 =     tan1slave(0,1)*clhs86 - clhs5*clhs70 + clhs5*clhs71 - clhs54*Dnormalslave01u110 + clhs62*clhs84 - clhs69*clhs81; // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS86 - CLHS5*CLHS70 + CLHS5*CLHS71 - CLHS54*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS62*CLHS84 - CLHS69*CLHS81
    const double clhs89 =     Dgapgu111; // DERIVATIVE(GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1)), U1(1,1))
    const double clhs90 =     Dtan1slave00u111; // DERIVATIVE(TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs91 =     Dtan1slave10u111; // DERIVATIVE(TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs92 =     Dtan1slave01u111; // DERIVATIVE(TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs93 =     Dtan1slave11u111; // DERIVATIVE(TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1))
    const double clhs94 =     clhs0*clhs32*(N1[1]*clhs42 + clhs59*(N1[0]*clhs90 + N1[1]*clhs91) + clhs60*(N1[0]*clhs92 + N1[1]*clhs93));
    const double clhs95 =     tan1slave(0,0)*clhs94 - clhs54*Dnormalslave00u111 - clhs55*clhs89 - clhs57*clhs6 + clhs6*clhs63 + clhs62*clhs90; // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS94 - CLHS54*DERIVATIVE(NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS55*CLHS89 - CLHS57*CLHS6 + CLHS6*CLHS63 + CLHS62*CLHS90
    const double clhs96 =     tan1slave(0,1)*clhs94 - clhs54*Dnormalslave01u111 - clhs6*clhs70 + clhs6*clhs71 + clhs62*clhs92 - clhs69*clhs89; // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS94 - CLHS54*DERIVATIVE(NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) - CLHS6*CLHS70 + CLHS6*CLHS71 + CLHS62*CLHS92 - CLHS69*CLHS89
    const double clhs97 =     Phi[1]*clhs0;
    const double clhs98 =     normalslave(1,0); // NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs99 =     clhs30*clhs98 + clhs33*clhs35;
    const double clhs100 =     normalslave(1,1); // NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double clhs101 =     clhs100*clhs30 + clhs35*clhs41;
    const double clhs102 =     clhs33*clhs43 + clhs40*clhs98;
    const double clhs103 =     clhs100*clhs40 + clhs41*clhs43;
    const double clhs104 =     clhs33*clhs47 + clhs46*clhs98;
    const double clhs105 =     clhs100*clhs46 + clhs41*clhs47;
    const double clhs106 =     clhs33*clhs51 + clhs50*clhs98;
    const double clhs107 =     clhs100*clhs50 + clhs41*clhs51;
    const double clhs108 =     clhs0*clhs98;
    const double clhs109 =     clhs29*clhs98;
    const double clhs110 =     clhs32*clhs33*clhs61;
    const double clhs111 =     -clhs1*clhs109 + clhs1*clhs110 - clhs108*clhs56 + tan1slave(1,0)*clhs67 - clhs54*Dnormalslave10u100 + clhs62*clhs64; // -CLHS1*CLHS109 + CLHS1*CLHS110 - CLHS108*CLHS56 + TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS67 - CLHS54*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS62*CLHS64
    const double clhs112 =     clhs0*clhs100;
    const double clhs113 =     clhs100*clhs29;
    const double clhs114 =     clhs32*clhs41*clhs61;
    const double clhs115 =     -clhs1*clhs113 + clhs1*clhs114 - clhs112*clhs56 + tan1slave(1,1)*clhs67 - clhs54*Dnormalslave11u100 + clhs62*clhs66; // -CLHS1*CLHS113 + CLHS1*CLHS114 - CLHS112*CLHS56 + TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS67 - CLHS54*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,0)) + CLHS62*CLHS66
    const double clhs116 =     -clhs108*clhs73 - clhs109*clhs4 + clhs110*clhs4 + tan1slave(1,0)*clhs78 - clhs54*Dnormalslave10u101 + clhs62*clhs75; // -CLHS108*CLHS73 - CLHS109*CLHS4 + CLHS110*CLHS4 + TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS78 - CLHS54*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS62*CLHS75
    const double clhs117 =     -clhs112*clhs73 - clhs113*clhs4 + clhs114*clhs4 + tan1slave(1,1)*clhs78 - clhs54*Dnormalslave11u101 + clhs62*clhs77; // -CLHS112*CLHS73 - CLHS113*CLHS4 + CLHS114*CLHS4 + TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS78 - CLHS54*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(0,1)) + CLHS62*CLHS77
    const double clhs118 =     -clhs108*clhs81 - clhs109*clhs5 + clhs110*clhs5 + tan1slave(1,0)*clhs86 - clhs54*Dnormalslave10u110 + clhs62*clhs83; // -CLHS108*CLHS81 - CLHS109*CLHS5 + CLHS110*CLHS5 + TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS86 - CLHS54*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS62*CLHS83
    const double clhs119 =     -clhs112*clhs81 - clhs113*clhs5 + clhs114*clhs5 + tan1slave(1,1)*clhs86 - clhs54*Dnormalslave11u110 + clhs62*clhs85; // -CLHS112*CLHS81 - CLHS113*CLHS5 + CLHS114*CLHS5 + TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS86 - CLHS54*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,0)) + CLHS62*CLHS85
    const double clhs120 =     -clhs108*clhs89 - clhs109*clhs6 + clhs110*clhs6 + tan1slave(1,0)*clhs94 - clhs54*Dnormalslave10u111 + clhs62*clhs91; // -CLHS108*CLHS89 - CLHS109*CLHS6 + CLHS110*CLHS6 + TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS94 - CLHS54*DERIVATIVE(NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS62*CLHS91
    const double clhs121 =     -clhs112*clhs89 - clhs113*clhs6 + clhs114*clhs6 + tan1slave(1,1)*clhs94 - clhs54*Dnormalslave11u111 + clhs62*clhs93; // -CLHS112*CLHS89 - CLHS113*CLHS6 + CLHS114*CLHS6 + TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))*CLHS94 - CLHS54*DERIVATIVE(NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1)), U1(1,1)) + CLHS62*CLHS93

    lhs(0,0)=0;
    lhs(0,1)=0;
    lhs(0,2)=0;
    lhs(0,3)=0;
    lhs(0,4)=-clhs1*clhs3;
    lhs(0,5)=-clhs3*clhs4;
    lhs(0,6)=-clhs3*clhs5;
    lhs(0,7)=-clhs3*clhs6;
    lhs(0,8)=clhs8;
    lhs(0,9)=0;
    lhs(0,10)=clhs9;
    lhs(0,11)=0;
    lhs(1,0)=0;
    lhs(1,1)=0;
    lhs(1,2)=0;
    lhs(1,3)=0;
    lhs(1,4)=-clhs1*clhs11;
    lhs(1,5)=-clhs11*clhs4;
    lhs(1,6)=-clhs11*clhs5;
    lhs(1,7)=-clhs11*clhs6;
    lhs(1,8)=0;
    lhs(1,9)=clhs8;
    lhs(1,10)=0;
    lhs(1,11)=clhs9;
    lhs(2,0)=0;
    lhs(2,1)=0;
    lhs(2,2)=0;
    lhs(2,3)=0;
    lhs(2,4)=-clhs1*clhs12;
    lhs(2,5)=-clhs12*clhs4;
    lhs(2,6)=-clhs12*clhs5;
    lhs(2,7)=-clhs12*clhs6;
    lhs(2,8)=clhs14;
    lhs(2,9)=0;
    lhs(2,10)=clhs15;
    lhs(2,11)=0;
    lhs(3,0)=0;
    lhs(3,1)=0;
    lhs(3,2)=0;
    lhs(3,3)=0;
    lhs(3,4)=-clhs1*clhs16;
    lhs(3,5)=-clhs16*clhs4;
    lhs(3,6)=-clhs16*clhs5;
    lhs(3,7)=-clhs16*clhs6;
    lhs(3,8)=0;
    lhs(3,9)=clhs14;
    lhs(3,10)=0;
    lhs(3,11)=clhs15;
    lhs(4,0)=0;
    lhs(4,1)=0;
    lhs(4,2)=0;
    lhs(4,3)=0;
    lhs(4,4)=clhs1*clhs17;
    lhs(4,5)=clhs17*clhs4;
    lhs(4,6)=clhs17*clhs5;
    lhs(4,7)=clhs17*clhs6;
    lhs(4,8)=clhs19;
    lhs(4,9)=0;
    lhs(4,10)=clhs20;
    lhs(4,11)=0;
    lhs(5,0)=0;
    lhs(5,1)=0;
    lhs(5,2)=0;
    lhs(5,3)=0;
    lhs(5,4)=clhs1*clhs21;
    lhs(5,5)=clhs21*clhs4;
    lhs(5,6)=clhs21*clhs5;
    lhs(5,7)=clhs21*clhs6;
    lhs(5,8)=0;
    lhs(5,9)=clhs19;
    lhs(5,10)=0;
    lhs(5,11)=clhs20;
    lhs(6,0)=0;
    lhs(6,1)=0;
    lhs(6,2)=0;
    lhs(6,3)=0;
    lhs(6,4)=clhs1*clhs22;
    lhs(6,5)=clhs22*clhs4;
    lhs(6,6)=clhs22*clhs5;
    lhs(6,7)=clhs22*clhs6;
    lhs(6,8)=clhs24;
    lhs(6,9)=0;
    lhs(6,10)=clhs25;
    lhs(6,11)=0;
    lhs(7,0)=0;
    lhs(7,1)=0;
    lhs(7,2)=0;
    lhs(7,3)=0;
    lhs(7,4)=clhs1*clhs26;
    lhs(7,5)=clhs26*clhs4;
    lhs(7,6)=clhs26*clhs5;
    lhs(7,7)=clhs26*clhs6;
    lhs(7,8)=0;
    lhs(7,9)=clhs24;
    lhs(7,10)=0;
    lhs(7,11)=clhs25;
    lhs(8,0)=clhs27*(GPnormal[0]*clhs36 + GPnormal[1]*clhs39);
    lhs(8,1)=clhs27*(GPnormal[0]*clhs44 + GPnormal[1]*clhs45);
    lhs(8,2)=clhs27*(GPnormal[0]*clhs48 + GPnormal[1]*clhs49);
    lhs(8,3)=clhs27*(GPnormal[0]*clhs52 + GPnormal[1]*clhs53);
    lhs(8,4)=-Phi[0]*(GPnormal[0]*clhs68 + GPnormal[1]*clhs72);
    lhs(8,5)=-Phi[0]*(GPnormal[0]*clhs79 + GPnormal[1]*clhs80);
    lhs(8,6)=-Phi[0]*(GPnormal[0]*clhs87 + GPnormal[1]*clhs88);
    lhs(8,7)=-Phi[0]*(GPnormal[0]*clhs95 + GPnormal[1]*clhs96);
    lhs(8,8)=0;
    lhs(8,9)=0;
    lhs(8,10)=0;
    lhs(8,11)=0;
    lhs(9,0)=clhs27*(GPtangent1[0]*clhs36 + GPtangent1[1]*clhs39);
    lhs(9,1)=clhs27*(GPtangent1[0]*clhs44 + GPtangent1[1]*clhs45);
    lhs(9,2)=clhs27*(GPtangent1[0]*clhs48 + GPtangent1[1]*clhs49);
    lhs(9,3)=clhs27*(GPtangent1[0]*clhs52 + GPtangent1[1]*clhs53);
    lhs(9,4)=-Phi[0]*(GPtangent1[0]*clhs68 + GPtangent1[1]*clhs72);
    lhs(9,5)=-Phi[0]*(GPtangent1[0]*clhs79 + GPtangent1[1]*clhs80);
    lhs(9,6)=-Phi[0]*(GPtangent1[0]*clhs87 + GPtangent1[1]*clhs88);
    lhs(9,7)=-Phi[0]*(GPtangent1[0]*clhs95 + GPtangent1[1]*clhs96);
    lhs(9,8)=0;
    lhs(9,9)=0;
    lhs(9,10)=0;
    lhs(9,11)=0;
    lhs(10,0)=clhs97*(GPnormal[0]*clhs99 + GPnormal[1]*clhs101);
    lhs(10,1)=clhs97*(GPnormal[0]*clhs102 + GPnormal[1]*clhs103);
    lhs(10,2)=clhs97*(GPnormal[0]*clhs104 + GPnormal[1]*clhs105);
    lhs(10,3)=clhs97*(GPnormal[0]*clhs106 + GPnormal[1]*clhs107);
    lhs(10,4)=-Phi[1]*(GPnormal[0]*clhs111 + GPnormal[1]*clhs115);
    lhs(10,5)=-Phi[1]*(GPnormal[0]*clhs116 + GPnormal[1]*clhs117);
    lhs(10,6)=-Phi[1]*(GPnormal[0]*clhs118 + GPnormal[1]*clhs119);
    lhs(10,7)=-Phi[1]*(GPnormal[0]*clhs120 + GPnormal[1]*clhs121);
    lhs(10,8)=0;
    lhs(10,9)=0;
    lhs(10,10)=0;
    lhs(10,11)=0;
    lhs(11,0)=clhs97*(GPtangent1[0]*clhs99 + GPtangent1[1]*clhs101);
    lhs(11,1)=clhs97*(GPtangent1[0]*clhs102 + GPtangent1[1]*clhs103);
    lhs(11,2)=clhs97*(GPtangent1[0]*clhs104 + GPtangent1[1]*clhs105);
    lhs(11,3)=clhs97*(GPtangent1[0]*clhs106 + GPtangent1[1]*clhs107);
    lhs(11,4)=-Phi[1]*(GPtangent1[0]*clhs111 + GPtangent1[1]*clhs115);
    lhs(11,5)=-Phi[1]*(GPtangent1[0]*clhs116 + GPtangent1[1]*clhs117);
    lhs(11,6)=-Phi[1]*(GPtangent1[0]*clhs118 + GPtangent1[1]*clhs119);
    lhs(11,7)=-Phi[1]*(GPtangent1[0]*clhs120 + GPtangent1[1]*clhs121);
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
    
    const double crhs0 =     detJ; // DETJ(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs1 =     N2[0]*crhs0;
    const double crhs2 =     Phi[0]*lm(0,0) + Phi[1]*lm(1,0);
    const double crhs3 =     Phi[0]*lm(0,1) + Phi[1]*lm(1,1);
    const double crhs4 =     N2[1]*crhs0;
    const double crhs5 =     N1[0]*crhs0;
    const double crhs6 =     N1[1]*crhs0;
    const double crhs7 =     Phi[0]*crhs0;
    const double crhs8 =     inner_prod(Gaps, N1); // GAPG(U1(0,0), U1(0,1), U1(1,0), U1(1,1), U2(0,0), U2(0,1), U2(1,0), U2(1,1))
    const double crhs9 =     tan1slave(0,0); // TAN1SLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs10 =     tan1slave(1,0); // TAN1SLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs11 =     tan1slave(0,1); // TAN1SLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs12 =     tan1slave(1,1); // TAN1SLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs13 =     ((N1[0]*crhs11 + N1[1]*crhs12)*(N1[0]*(Dt*v1(0,1)) + N1[1]*(Dt*v1(1,1)) - N2[0]*(Dt*v2(0,1)) - N2[1]*(Dt*v2(1,1))) + (N1[0]*crhs9 + N1[1]*crhs10)*(N1[0]*(Dt*v1(0,0)) + N1[1]*(Dt*v1(1,0)) - N2[0]*(Dt*v2(0,0)) - N2[1]*(Dt*v2(1,0))))/Dt;
    const double crhs14 =     crhs13*crhs9 - crhs8*normalslave(0,0); // CRHS13*CRHS9 - CRHS8*NORMALSLAVE(0,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs15 =     crhs11*crhs13 - crhs8*normalslave(0,1); // CRHS11*CRHS13 - CRHS8*NORMALSLAVE(0,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs16 =     Phi[1]*crhs0;
    const double crhs17 =     crhs10*crhs13 - crhs8*normalslave(1,0); // CRHS10*CRHS13 - CRHS8*NORMALSLAVE(1,0)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))
    const double crhs18 =     crhs12*crhs13 - crhs8*normalslave(1,1); // CRHS12*CRHS13 - CRHS8*NORMALSLAVE(1,1)(U1(0,0), U1(0,1), U1(1,0), U1(1,1))

    rhs[0]=crhs1*crhs2;
    rhs[1]=crhs1*crhs3;
    rhs[2]=crhs2*crhs4;
    rhs[3]=crhs3*crhs4;
    rhs[4]=-crhs2*crhs5;
    rhs[5]=-crhs3*crhs5;
    rhs[6]=-crhs2*crhs6;
    rhs[7]=-crhs3*crhs6;
    rhs[8]=crhs7*(GPnormal[0]*crhs14 + GPnormal[1]*crhs15);
    rhs[9]=crhs7*(GPtangent1[0]*crhs14 + GPtangent1[1]*crhs15);
    rhs[10]=crhs16*(GPnormal[0]*crhs17 + GPnormal[1]*crhs18);
    rhs[11]=crhs16*(GPtangent1[0]*crhs17 + GPtangent1[1]*crhs18);

    
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