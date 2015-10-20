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
#include "ballast2cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast2Cluster3D::Ballast2Cluster3D() : Cluster3D() {}
            
      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast2Cluster3D::Ballast2Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast2Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast2Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast2Cluster3D::~Ballast2Cluster3D() {}
      
    
    void Ballast2Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 28;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.241 (in meters) was the medium diameter of the rock in GiD (Rock3_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.8058,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.8058;
     
        mListOfCoordinates[ 0][0]= 0.031233; mListOfCoordinates[ 0][1]= 0.039652; mListOfCoordinates[ 0][2]=-0.071418;
        mListOfCoordinates[ 1][0]=-0.500669; mListOfCoordinates[ 1][1]=-0.021090; mListOfCoordinates[ 1][2]= 0.016663;
        mListOfCoordinates[ 2][0]= 0.043345; mListOfCoordinates[ 2][1]= 0.097122; mListOfCoordinates[ 2][2]=-0.064253;
        mListOfCoordinates[ 3][0]=-0.011465; mListOfCoordinates[ 3][1]=-0.039297; mListOfCoordinates[ 3][2]=-0.035967;
        mListOfCoordinates[ 4][0]=-0.431808; mListOfCoordinates[ 4][1]=-0.265948; mListOfCoordinates[ 4][2]= 0.090583;
        mListOfCoordinates[ 5][0]=-0.028961; mListOfCoordinates[ 5][1]= 0.196849; mListOfCoordinates[ 5][2]=-0.114789;
        mListOfCoordinates[ 6][0]=-0.331251; mListOfCoordinates[ 6][1]=-0.068507; mListOfCoordinates[ 6][2]= 0.057451;
        mListOfCoordinates[ 7][0]= 0.222432; mListOfCoordinates[ 7][1]=-0.056498; mListOfCoordinates[ 7][2]=-0.011545;
        mListOfCoordinates[ 8][0]=-0.477574; mListOfCoordinates[ 8][1]=-0.027498; mListOfCoordinates[ 8][2]= 0.012307;
        mListOfCoordinates[ 9][0]=-0.473843; mListOfCoordinates[ 9][1]=-0.171022; mListOfCoordinates[ 9][2]= 0.066582;
        mListOfCoordinates[10][0]= 0.282513; mListOfCoordinates[10][1]= 0.176342; mListOfCoordinates[10][2]=-0.110049;
        mListOfCoordinates[11][0]= 0.422548; mListOfCoordinates[11][1]=-0.235309; mListOfCoordinates[11][2]= 0.080327;
        mListOfCoordinates[12][0]= 0.145562; mListOfCoordinates[12][1]= 0.308839; mListOfCoordinates[12][2]=-0.148519;
        mListOfCoordinates[13][0]= 0.058483; mListOfCoordinates[13][1]= 0.084395; mListOfCoordinates[13][2]=-0.062576;
        mListOfCoordinates[14][0]=-0.205940; mListOfCoordinates[14][1]=-0.096927; mListOfCoordinates[14][2]= 0.005223;
        mListOfCoordinates[15][0]= 0.139789; mListOfCoordinates[15][1]=-0.180467; mListOfCoordinates[15][2]= 0.035945;
        mListOfCoordinates[16][0]=-0.273070; mListOfCoordinates[16][1]=-0.304840; mListOfCoordinates[16][2]= 0.078110;
        mListOfCoordinates[17][0]= 0.398024; mListOfCoordinates[17][1]=-0.042855; mListOfCoordinates[17][2]=-0.014173;
        mListOfCoordinates[18][0]=-0.290355; mListOfCoordinates[18][1]=-0.272346; mListOfCoordinates[18][2]= 0.075184;
        mListOfCoordinates[19][0]= 0.030095; mListOfCoordinates[19][1]= 0.335917; mListOfCoordinates[19][2]=-0.164594;
        mListOfCoordinates[20][0]= 0.113366; mListOfCoordinates[20][1]= 0.276485; mListOfCoordinates[20][2]=-0.133785;
        mListOfCoordinates[21][0]= 0.255818; mListOfCoordinates[21][1]=-0.133544; mListOfCoordinates[21][2]= 0.048835;
        mListOfCoordinates[22][0]= 0.145044; mListOfCoordinates[22][1]=-0.006133; mListOfCoordinates[22][2]=-0.019225;
        mListOfCoordinates[23][0]= 0.093352; mListOfCoordinates[23][1]=-0.101307; mListOfCoordinates[23][2]=-0.005626;
        mListOfCoordinates[24][0]=-0.076947; mListOfCoordinates[24][1]= 0.123827; mListOfCoordinates[24][2]=-0.086978;
        mListOfCoordinates[25][0]=-0.207412; mListOfCoordinates[25][1]= 0.111712; mListOfCoordinates[25][2]=-0.072115;
        mListOfCoordinates[26][0]=-0.021493; mListOfCoordinates[26][1]= 0.290421; mListOfCoordinates[26][2]=-0.149556;
        mListOfCoordinates[27][0]=-0.024671; mListOfCoordinates[27][1]= 0.093951; mListOfCoordinates[27][2]=-0.073355;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0] = 0.339049;
        mListOfRadii[ 1] = 0.233528;
        mListOfRadii[ 2] = 0.388827;
        mListOfRadii[ 3] = 0.331505;
        mListOfRadii[ 4] = 0.227643;
        mListOfRadii[ 5] = 0.374958;
        mListOfRadii[ 6] = 0.262937;
        mListOfRadii[ 7] = 0.291307;
        mListOfRadii[ 8] = 0.245471;
        mListOfRadii[ 9] = 0.253552;
        mListOfRadii[10] = 0.265224;
        mListOfRadii[11] = 0.269778;
        mListOfRadii[12] = 0.288632;
        mListOfRadii[13] = 0.380106;
        mListOfRadii[14] = 0.300409;
        mListOfRadii[15] = 0.333144;
        mListOfRadii[16] = 0.110335;
        mListOfRadii[17] = 0.268086;
        mListOfRadii[18] = 0.202354;
        mListOfRadii[19] = 0.293389;
        mListOfRadii[20] = 0.284530;
        mListOfRadii[21] = 0.321070;
        mListOfRadii[22] = 0.365016;
        mListOfRadii[23] = 0.329093;
        mListOfRadii[24] = 0.345479;
        mListOfRadii[25] = 0.319487;
        mListOfRadii[26] = 0.268785;
        mListOfRadii[27] = 0.387152;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.632108 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.096309691;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.133694924;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.165530265; 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast2Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast2Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast2Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

