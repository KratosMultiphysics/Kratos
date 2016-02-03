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
        
        int number_of_spheres = 43;
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
 
        mListOfCoordinates[ 0][0]= 0.282367; mListOfCoordinates[ 0][1]=-0.191849; mListOfCoordinates[ 0][2]=-0.076992;
        mListOfCoordinates[ 1][0]=-0.326788; mListOfCoordinates[ 1][1]= 0.174716; mListOfCoordinates[ 1][2]= 0.211693;
        mListOfCoordinates[ 2][0]= 0.618972; mListOfCoordinates[ 2][1]= 0.184459; mListOfCoordinates[ 2][2]= 0.056258;
        mListOfCoordinates[ 3][0]=-0.165500; mListOfCoordinates[ 3][1]=-0.040034; mListOfCoordinates[ 3][2]=-0.279015;
        mListOfCoordinates[ 4][0]= 0.161947; mListOfCoordinates[ 4][1]=-0.681441; mListOfCoordinates[ 4][2]= 0.175952;
        mListOfCoordinates[ 5][0]=-0.523620; mListOfCoordinates[ 5][1]= 0.053050; mListOfCoordinates[ 5][2]= 0.227297;
        mListOfCoordinates[ 6][0]=-0.440551; mListOfCoordinates[ 6][1]=-0.200588; mListOfCoordinates[ 6][2]=-0.201722;
        mListOfCoordinates[ 7][0]=-0.134893; mListOfCoordinates[ 7][1]= 0.387507; mListOfCoordinates[ 7][2]= 0.215766;
        mListOfCoordinates[ 8][0]= 0.168112; mListOfCoordinates[ 8][1]= 0.037289; mListOfCoordinates[ 8][2]= 0.086013;
        mListOfCoordinates[ 9][0]=-0.579941; mListOfCoordinates[ 9][1]=-0.233373; mListOfCoordinates[ 9][2]= 0.257002;
        mListOfCoordinates[10][0]=-0.129074; mListOfCoordinates[10][1]=-0.201429; mListOfCoordinates[10][2]= 0.205499;
        mListOfCoordinates[11][0]=-0.652185; mListOfCoordinates[11][1]= 0.100605; mListOfCoordinates[11][2]= 0.257928;
        mListOfCoordinates[12][0]= 0.240681; mListOfCoordinates[12][1]= 0.466673; mListOfCoordinates[12][2]=-0.336280;
        mListOfCoordinates[13][0]=-0.145779; mListOfCoordinates[13][1]= 0.159040; mListOfCoordinates[13][2]= 0.201385;
        mListOfCoordinates[14][0]= 0.212846; mListOfCoordinates[14][1]= 0.181019; mListOfCoordinates[14][2]=-0.278521;
        mListOfCoordinates[15][0]=-0.369727; mListOfCoordinates[15][1]= 0.387710; mListOfCoordinates[15][2]= 0.255948;
        mListOfCoordinates[16][0]= 0.246971; mListOfCoordinates[16][1]= 0.667141; mListOfCoordinates[16][2]= 0.130134;
        mListOfCoordinates[17][0]= 0.262716; mListOfCoordinates[17][1]=-0.174559; mListOfCoordinates[17][2]=-0.319388;
        mListOfCoordinates[18][0]= 0.171921; mListOfCoordinates[18][1]= 0.294305; mListOfCoordinates[18][2]= 0.124303;
        mListOfCoordinates[19][0]= 0.120008; mListOfCoordinates[19][1]=-0.357930; mListOfCoordinates[19][2]= 0.195817;
        mListOfCoordinates[20][0]= 0.007487; mListOfCoordinates[20][1]=-0.576526; mListOfCoordinates[20][2]= 0.252003;
        mListOfCoordinates[21][0]= 0.405718; mListOfCoordinates[21][1]= 0.283885; mListOfCoordinates[21][2]= 0.079109;
        mListOfCoordinates[22][0]=-0.708690; mListOfCoordinates[22][1]=-0.056333; mListOfCoordinates[22][2]= 0.261727;
        mListOfCoordinates[23][0]= 0.300941; mListOfCoordinates[23][1]=-0.463648; mListOfCoordinates[23][2]= 0.075084;
        mListOfCoordinates[24][0]=-0.412075; mListOfCoordinates[24][1]=-0.101999; mListOfCoordinates[24][2]= 0.206846;
        mListOfCoordinates[25][0]= 0.169464; mListOfCoordinates[25][1]=-0.110553; mListOfCoordinates[25][2]=-0.515637;
        mListOfCoordinates[26][0]= 0.280113; mListOfCoordinates[26][1]=-0.401824; mListOfCoordinates[26][2]=-0.396139;
        mListOfCoordinates[27][0]= 0.143246; mListOfCoordinates[27][1]= 0.241375; mListOfCoordinates[27][2]=-0.233068;
        mListOfCoordinates[28][0]=-0.435998; mListOfCoordinates[28][1]=-0.376401; mListOfCoordinates[28][2]= 0.247957;
        mListOfCoordinates[29][0]=-0.148025; mListOfCoordinates[29][1]= 0.060517; mListOfCoordinates[29][2]= 0.174888;
        mListOfCoordinates[30][0]= 0.089942; mListOfCoordinates[30][1]= 0.515272; mListOfCoordinates[30][2]= 0.168312;
        mListOfCoordinates[31][0]=-0.126219; mListOfCoordinates[31][1]=-0.056626; mListOfCoordinates[31][2]=-0.134893;
        mListOfCoordinates[32][0]= 0.304421; mListOfCoordinates[32][1]=-0.749720; mListOfCoordinates[32][2]= 0.099991;
        mListOfCoordinates[33][0]= 0.341687; mListOfCoordinates[33][1]= 0.518540; mListOfCoordinates[33][2]= 0.080518;
        mListOfCoordinates[34][0]= 0.284475; mListOfCoordinates[34][1]=-0.370066; mListOfCoordinates[34][2]=-0.254208;
        mListOfCoordinates[35][0]= 0.106085; mListOfCoordinates[35][1]= 0.069167; mListOfCoordinates[35][2]=-0.528590;
        mListOfCoordinates[36][0]=-0.531031; mListOfCoordinates[36][1]= 0.220978; mListOfCoordinates[36][2]= 0.248842;
        mListOfCoordinates[37][0]= 0.245845; mListOfCoordinates[37][1]= 0.486141; mListOfCoordinates[37][2]=-0.109937;
        mListOfCoordinates[38][0]= 0.251911; mListOfCoordinates[38][1]= 0.361353; mListOfCoordinates[38][2]=-0.379873;
        mListOfCoordinates[39][0]= 0.233034; mListOfCoordinates[39][1]=-0.194649; mListOfCoordinates[39][2]= 0.109307;
        mListOfCoordinates[40][0]= 0.436142; mListOfCoordinates[40][1]= 0.123600; mListOfCoordinates[40][2]=-0.369742;
        mListOfCoordinates[41][0]=-0.228657; mListOfCoordinates[41][1]=-0.565929; mListOfCoordinates[41][2]= 0.246015;
        mListOfCoordinates[42][0]=-0.028287; mListOfCoordinates[42][1]= 0.131147; mListOfCoordinates[42][2]=-0.437569;

        
        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.423948;
        mListOfRadii[ 1]= 0.357812;
        mListOfRadii[ 2]= 0.041519;
        mListOfRadii[ 3]= 0.363435;
        mListOfRadii[ 4]= 0.162124;
        mListOfRadii[ 5]= 0.284758;
        mListOfRadii[ 6]= 0.157085;
        mListOfRadii[ 7]= 0.317212;
        mListOfRadii[ 8]= 0.505934;
        mListOfRadii[ 9]= 0.246213;
        mListOfRadii[10]= 0.328157;
        mListOfRadii[11]= 0.249119;
        mListOfRadii[12]= 0.249428;
        mListOfRadii[13]= 0.329662;
        mListOfRadii[14]= 0.381627;
        mListOfRadii[15]= 0.233132;
        mListOfRadii[16]= 0.211288;
        mListOfRadii[17]= 0.417010;
        mListOfRadii[18]= 0.418231;
        mListOfRadii[19]= 0.289032;
        mListOfRadii[20]= 0.169991;
        mListOfRadii[21]= 0.255457;
        mListOfRadii[22]= 0.249474;
        mListOfRadii[23]= 0.325893;
        mListOfRadii[24]= 0.332234;
        mListOfRadii[25]= 0.299024;
        mListOfRadii[26]= 0.289954;
        mListOfRadii[27]= 0.400860;
        mListOfRadii[28]= 0.236031;
        mListOfRadii[29]= 0.345723;
        mListOfRadii[30]= 0.332815;
        mListOfRadii[31]= 0.348808;
        mListOfRadii[32]= 0.220961;
        mListOfRadii[33]= 0.254481;
        mListOfRadii[34]= 0.318738;
        mListOfRadii[35]= 0.295372;
        mListOfRadii[36]= 0.259937;
        mListOfRadii[37]= 0.241836;
        mListOfRadii[38]= 0.249171;
        mListOfRadii[39]= 0.360297;
        mListOfRadii[40]= 0.181715;
        mListOfRadii[41]= 0.135446;
        mListOfRadii[42]= 0.243484;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 1.299807 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.180015;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.205021;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.243482;
  
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

