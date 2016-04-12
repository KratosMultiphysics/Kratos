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
      
    
    void Ballast3Cluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres = 43;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.446688 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.691234,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.691234;     
 
        mListOfCoordinates[ 0][0]= 0.250209; mListOfCoordinates[ 0][1]=-0.191917; mListOfCoordinates[ 0][2]=-0.066919;
        mListOfCoordinates[ 1][0]=-0.358946; mListOfCoordinates[ 1][1]= 0.174648; mListOfCoordinates[ 1][2]= 0.221766;
        mListOfCoordinates[ 2][0]= 0.586814; mListOfCoordinates[ 2][1]= 0.184391; mListOfCoordinates[ 2][2]= 0.066331;
        mListOfCoordinates[ 3][0]=-0.197658; mListOfCoordinates[ 3][1]=-0.040102; mListOfCoordinates[ 3][2]=-0.268942;
        mListOfCoordinates[ 4][0]= 0.129789; mListOfCoordinates[ 4][1]=-0.681509; mListOfCoordinates[ 4][2]= 0.186025;
        mListOfCoordinates[ 5][0]=-0.555778; mListOfCoordinates[ 5][1]= 0.052982; mListOfCoordinates[ 5][2]= 0.237370;
        mListOfCoordinates[ 6][0]=-0.472709; mListOfCoordinates[ 6][1]=-0.200656; mListOfCoordinates[ 6][2]=-0.191649;
        mListOfCoordinates[ 7][0]=-0.167051; mListOfCoordinates[ 7][1]= 0.387439; mListOfCoordinates[ 7][2]= 0.225839;
        mListOfCoordinates[ 8][0]= 0.135954; mListOfCoordinates[ 8][1]= 0.037221; mListOfCoordinates[ 8][2]= 0.096086;
        mListOfCoordinates[ 9][0]=-0.612099; mListOfCoordinates[ 9][1]=-0.233441; mListOfCoordinates[ 9][2]= 0.267075;
        mListOfCoordinates[10][0]=-0.161232; mListOfCoordinates[10][1]=-0.201497; mListOfCoordinates[10][2]= 0.215572;
        mListOfCoordinates[11][0]=-0.684343; mListOfCoordinates[11][1]= 0.100537; mListOfCoordinates[11][2]= 0.268001;
        mListOfCoordinates[12][0]= 0.208523; mListOfCoordinates[12][1]= 0.466605; mListOfCoordinates[12][2]=-0.326207;
        mListOfCoordinates[13][0]=-0.177937; mListOfCoordinates[13][1]= 0.158972; mListOfCoordinates[13][2]= 0.211458;
        mListOfCoordinates[14][0]= 0.180688; mListOfCoordinates[14][1]= 0.180951; mListOfCoordinates[14][2]=-0.268448;
        mListOfCoordinates[15][0]=-0.401885; mListOfCoordinates[15][1]= 0.387642; mListOfCoordinates[15][2]= 0.266021;
        mListOfCoordinates[16][0]= 0.214813; mListOfCoordinates[16][1]= 0.667073; mListOfCoordinates[16][2]= 0.140207;
        mListOfCoordinates[17][0]= 0.230558; mListOfCoordinates[17][1]=-0.174627; mListOfCoordinates[17][2]=-0.309315;
        mListOfCoordinates[18][0]= 0.139763; mListOfCoordinates[18][1]= 0.294237; mListOfCoordinates[18][2]= 0.134376;
        mListOfCoordinates[19][0]= 0.087850; mListOfCoordinates[19][1]=-0.357998; mListOfCoordinates[19][2]= 0.205890;
        mListOfCoordinates[20][0]=-0.024671; mListOfCoordinates[20][1]=-0.576594; mListOfCoordinates[20][2]= 0.262076;
        mListOfCoordinates[21][0]= 0.373560; mListOfCoordinates[21][1]= 0.283817; mListOfCoordinates[21][2]= 0.089182;
        mListOfCoordinates[22][0]=-0.740848; mListOfCoordinates[22][1]=-0.056401; mListOfCoordinates[22][2]= 0.271800;
        mListOfCoordinates[23][0]= 0.268783; mListOfCoordinates[23][1]=-0.463716; mListOfCoordinates[23][2]= 0.085157;
        mListOfCoordinates[24][0]=-0.444233; mListOfCoordinates[24][1]=-0.102067; mListOfCoordinates[24][2]= 0.216919;
        mListOfCoordinates[25][0]= 0.137306; mListOfCoordinates[25][1]=-0.110621; mListOfCoordinates[25][2]=-0.505564;
        mListOfCoordinates[26][0]= 0.247955; mListOfCoordinates[26][1]=-0.401892; mListOfCoordinates[26][2]=-0.386066;
        mListOfCoordinates[27][0]= 0.111088; mListOfCoordinates[27][1]= 0.241307; mListOfCoordinates[27][2]=-0.222995;
        mListOfCoordinates[28][0]=-0.468156; mListOfCoordinates[28][1]=-0.376469; mListOfCoordinates[28][2]= 0.258030;
        mListOfCoordinates[29][0]=-0.180183; mListOfCoordinates[29][1]= 0.060449; mListOfCoordinates[29][2]= 0.184961;
        mListOfCoordinates[30][0]= 0.057784; mListOfCoordinates[30][1]= 0.515204; mListOfCoordinates[30][2]= 0.178385;
        mListOfCoordinates[31][0]=-0.158377; mListOfCoordinates[31][1]=-0.056694; mListOfCoordinates[31][2]=-0.124820;
        mListOfCoordinates[32][0]= 0.272263; mListOfCoordinates[32][1]=-0.749788; mListOfCoordinates[32][2]= 0.110064;
        mListOfCoordinates[33][0]= 0.309529; mListOfCoordinates[33][1]= 0.518472; mListOfCoordinates[33][2]= 0.090591;
        mListOfCoordinates[34][0]= 0.252317; mListOfCoordinates[34][1]=-0.370134; mListOfCoordinates[34][2]=-0.244135;
        mListOfCoordinates[35][0]= 0.073927; mListOfCoordinates[35][1]= 0.069099; mListOfCoordinates[35][2]=-0.518517;
        mListOfCoordinates[36][0]=-0.563189; mListOfCoordinates[36][1]= 0.220910; mListOfCoordinates[36][2]= 0.258915;
        mListOfCoordinates[37][0]= 0.213687; mListOfCoordinates[37][1]= 0.486073; mListOfCoordinates[37][2]=-0.099864;
        mListOfCoordinates[38][0]= 0.219753; mListOfCoordinates[38][1]= 0.361285; mListOfCoordinates[38][2]=-0.369800;
        mListOfCoordinates[39][0]= 0.200876; mListOfCoordinates[39][1]=-0.194717; mListOfCoordinates[39][2]= 0.119380;
        mListOfCoordinates[40][0]= 0.403984; mListOfCoordinates[40][1]= 0.123532; mListOfCoordinates[40][2]=-0.359669;
        mListOfCoordinates[41][0]=-0.260815; mListOfCoordinates[41][1]=-0.565997; mListOfCoordinates[41][2]= 0.256088;
        mListOfCoordinates[42][0]=-0.060445; mListOfCoordinates[42][1]= 0.131079; mListOfCoordinates[42][2]=-0.427496;

        
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
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 1.585343 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.215133;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.246313;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.283201;
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast3Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast3Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Ballast3Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Ballast3Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

