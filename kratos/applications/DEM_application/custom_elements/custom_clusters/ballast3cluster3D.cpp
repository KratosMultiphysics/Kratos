//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "ballast3cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast3Cluster3D::Ballast3Cluster3D() : Cluster3D() {}
            
      
    Ballast3Cluster3D::Ballast3Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast3Cluster3D::Ballast3Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast3Cluster3D::Ballast3Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast3Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast3Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast3Cluster3D::~Ballast3Cluster3D() {}
      
    
    void Ballast3Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 25;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.725 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.5797,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.5797;
                
//         double a, b, c;
//         
//         a = 0.572749; // this is the original semi-axis (1.1455/2)
//         b = a; c = b;
  
        mListOfCoordinates[ 0][0]= 0.604999; mListOfCoordinates[ 0][1]=-0.181178; mListOfCoordinates[ 0][2]=-0.535269;
        mListOfCoordinates[ 1][0]= 0.353815; mListOfCoordinates[ 1][1]=-0.331701; mListOfCoordinates[ 1][2]=-0.357472;
        mListOfCoordinates[ 2][0]=-0.400212; mListOfCoordinates[ 2][1]=-0.180159; mListOfCoordinates[ 2][2]=-0.501756;
        mListOfCoordinates[ 3][0]=-0.320302; mListOfCoordinates[ 3][1]=-0.039794; mListOfCoordinates[ 3][2]=-0.421962;
        mListOfCoordinates[ 4][0]= 0.632890; mListOfCoordinates[ 4][1]= 0.159042; mListOfCoordinates[ 4][2]=-0.590095;
        mListOfCoordinates[ 5][0]=-0.172050; mListOfCoordinates[ 5][1]= 0.505000; mListOfCoordinates[ 5][2]=-0.234401;
        mListOfCoordinates[ 6][0]=-0.113808; mListOfCoordinates[ 6][1]= 0.123756; mListOfCoordinates[ 6][2]=-0.246717;
        mListOfCoordinates[ 7][0]= 0.099478; mListOfCoordinates[ 7][1]= 0.520025; mListOfCoordinates[ 7][2]=-0.331344;
        mListOfCoordinates[ 8][0]=-0.134570; mListOfCoordinates[ 8][1]=-0.054170; mListOfCoordinates[ 8][2]=-0.229200;
        mListOfCoordinates[ 9][0]=-0.207965; mListOfCoordinates[ 9][1]= 0.060653; mListOfCoordinates[ 9][2]=-0.535985;
        mListOfCoordinates[10][0]=-0.371236; mListOfCoordinates[10][1]= 0.376690; mListOfCoordinates[10][2]=-0.482355;
        mListOfCoordinates[11][0]= 0.063781; mListOfCoordinates[11][1]=-0.110294; mListOfCoordinates[11][2]=-0.604237;
        mListOfCoordinates[12][0]=-0.409960; mListOfCoordinates[12][1]=-0.372155; mListOfCoordinates[12][2]=-0.366815;
        mListOfCoordinates[13][0]= 0.364803; mListOfCoordinates[13][1]= 0.302742; mListOfCoordinates[13][2]=-0.434199;
        mListOfCoordinates[14][0]=-0.195305; mListOfCoordinates[14][1]=-0.636214; mListOfCoordinates[14][2]=-0.092966;
        mListOfCoordinates[15][0]= 0.075983; mListOfCoordinates[15][1]=-0.450965; mListOfCoordinates[15][2]=-0.162783;
        mListOfCoordinates[16][0]=-0.245956; mListOfCoordinates[16][1]=-0.304217; mListOfCoordinates[16][2]=-0.132043;
        mListOfCoordinates[17][0]=-0.246510; mListOfCoordinates[17][1]= 0.155701; mListOfCoordinates[17][2]=-0.449807;
        mListOfCoordinates[18][0]=-0.426278; mListOfCoordinates[18][1]=-0.285532; mListOfCoordinates[18][2]=-0.400271;
        mListOfCoordinates[19][0]= 0.523805; mListOfCoordinates[19][1]= 0.033851; mListOfCoordinates[19][2]=-0.531721;
        mListOfCoordinates[20][0]= 0.072404; mListOfCoordinates[20][1]= 0.162301; mListOfCoordinates[20][2]=-0.308841;
        mListOfCoordinates[21][0]=-0.077444; mListOfCoordinates[21][1]= 0.339949; mListOfCoordinates[21][2]=-0.251502;
        mListOfCoordinates[22][0]= 0.280536; mListOfCoordinates[22][1]=-0.037411; mListOfCoordinates[22][2]=-0.416817;
        mListOfCoordinates[23][0]=-0.354080; mListOfCoordinates[23][1]=-0.224467; mListOfCoordinates[23][2]=-0.322169;
        mListOfCoordinates[24][0]=-0.090110; mListOfCoordinates[24][1]= 0.012351; mListOfCoordinates[24][2]=-0.262796;

        
        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0] = 0.414595;
        mListOfRadii[ 1] = 0.291621;
        mListOfRadii[ 2] = 0.449349;
        mListOfRadii[ 3] = 0.472670;
        mListOfRadii[ 4] = 0.349904;
        mListOfRadii[ 5] = 0.403162;
        mListOfRadii[ 6] = 0.479166;
        mListOfRadii[ 7] = 0.351946;
        mListOfRadii[ 8] = 0.510088;
        mListOfRadii[ 9] = 0.483440;
        mListOfRadii[10] = 0.330173;
        mListOfRadii[11] = 0.387938;
        mListOfRadii[12] = 0.340417;
        mListOfRadii[13] = 0.336720;
        mListOfRadii[14] = 0.289903;
        mListOfRadii[15] = 0.309353;
        mListOfRadii[16] = 0.307321;
        mListOfRadii[17] = 0.443707;
        mListOfRadii[18] = 0.341255;
        mListOfRadii[19] = 0.332425;
        mListOfRadii[20] = 0.401344;
        mListOfRadii[21] = 0.427726;
        mListOfRadii[22] = 0.422835;
        mListOfRadii[23] = 0.357376;
        mListOfRadii[24] = 0.487111;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 1.910646 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.230156973;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.261441594;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.311897872;
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast3Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast3Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast3Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

