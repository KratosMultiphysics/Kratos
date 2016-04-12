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
      
    
    void Ballast2Cluster3D::CustomInitialize(ProcessInfo& r_process_info) {
        
        int number_of_spheres = 33;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.038093 (in meters) was the medium diameter of the rock in GiD (Rock3_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.963304,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.963304438;
     
        mListOfCoordinates[ 0][0]=-0.530356; mListOfCoordinates[ 0][1]=-0.258546; mListOfCoordinates[ 0][2]=-0.129540;
        mListOfCoordinates[ 1][0]= 0.099702; mListOfCoordinates[ 1][1]=-0.009735; mListOfCoordinates[ 1][2]=-0.025147;
        mListOfCoordinates[ 2][0]= 0.063848; mListOfCoordinates[ 2][1]= 0.452319; mListOfCoordinates[ 2][2]=-0.133395;
        mListOfCoordinates[ 3][0]= 0.073609; mListOfCoordinates[ 3][1]= 0.240192; mListOfCoordinates[ 3][2]=-0.067885;
        mListOfCoordinates[ 4][0]=-0.484912; mListOfCoordinates[ 4][1]=-0.344937; mListOfCoordinates[ 4][2]=-0.078579;
        mListOfCoordinates[ 5][0]=-0.345917; mListOfCoordinates[ 5][1]=-0.177503; mListOfCoordinates[ 5][2]=-0.043127;
        mListOfCoordinates[ 6][0]= 0.129778; mListOfCoordinates[ 6][1]= 0.523061; mListOfCoordinates[ 6][2]=-0.063640;
        mListOfCoordinates[ 7][0]=-0.041941; mListOfCoordinates[ 7][1]= 0.131590; mListOfCoordinates[ 7][2]= 0.224032;
        mListOfCoordinates[ 8][0]= 0.171637; mListOfCoordinates[ 8][1]=-0.332287; mListOfCoordinates[ 8][2]= 0.226044;
        mListOfCoordinates[ 9][0]= 0.110953; mListOfCoordinates[ 9][1]=-0.228321; mListOfCoordinates[ 9][2]= 0.169516;
        mListOfCoordinates[10][0]=-0.145833; mListOfCoordinates[10][1]= 0.035002; mListOfCoordinates[10][2]=-0.215776;
        mListOfCoordinates[11][0]= 0.045387; mListOfCoordinates[11][1]=-0.046993; mListOfCoordinates[11][2]=-0.048467;
        mListOfCoordinates[12][0]= 0.188427; mListOfCoordinates[12][1]= 0.430190; mListOfCoordinates[12][2]= 0.007965;
        mListOfCoordinates[13][0]= 0.162006; mListOfCoordinates[13][1]=-0.068790; mListOfCoordinates[13][2]= 0.073529;
        mListOfCoordinates[14][0]=-0.066027; mListOfCoordinates[14][1]=-0.095492; mListOfCoordinates[14][2]=-0.072093;
        mListOfCoordinates[15][0]= 0.035334; mListOfCoordinates[15][1]=-0.078797; mListOfCoordinates[15][2]= 0.060776;
        mListOfCoordinates[16][0]=-0.032452; mListOfCoordinates[16][1]=-0.313982; mListOfCoordinates[16][2]= 0.064984;
        mListOfCoordinates[17][0]=-0.240370; mListOfCoordinates[17][1]=-0.377642; mListOfCoordinates[17][2]=-0.007623;
        mListOfCoordinates[18][0]=-0.229537; mListOfCoordinates[18][1]=-0.164233; mListOfCoordinates[18][2]=-0.037817;
        mListOfCoordinates[19][0]= 0.059738; mListOfCoordinates[19][1]= 0.367246; mListOfCoordinates[19][2]= 0.046050;
        mListOfCoordinates[20][0]= 0.224779; mListOfCoordinates[20][1]=-0.088784; mListOfCoordinates[20][2]= 0.117009;
        mListOfCoordinates[21][0]=-0.043973; mListOfCoordinates[21][1]= 0.125961; mListOfCoordinates[21][2]=-0.200320;
        mListOfCoordinates[22][0]=-0.197870; mListOfCoordinates[22][1]=-0.059610; mListOfCoordinates[22][2]=-0.139197;
        mListOfCoordinates[23][0]=-0.003357; mListOfCoordinates[23][1]=-0.000498; mListOfCoordinates[23][2]=-0.128394;
        mListOfCoordinates[24][0]= 0.255713; mListOfCoordinates[24][1]= 0.037975; mListOfCoordinates[24][2]= 0.065858;
        mListOfCoordinates[25][0]=-0.047541; mListOfCoordinates[25][1]= 0.338212; mListOfCoordinates[25][2]=-0.273544;
        mListOfCoordinates[26][0]=-0.379316; mListOfCoordinates[26][1]=-0.347149; mListOfCoordinates[26][2]=-0.059194;
        mListOfCoordinates[27][0]=-0.339100; mListOfCoordinates[27][1]=-0.058273; mListOfCoordinates[27][2]=-0.145909;
        mListOfCoordinates[28][0]=-0.029087; mListOfCoordinates[28][1]= 0.471523; mListOfCoordinates[28][2]=-0.223240;
        mListOfCoordinates[29][0]=-0.165543; mListOfCoordinates[29][1]=-0.068171; mListOfCoordinates[29][2]= 0.068490;
        mListOfCoordinates[30][0]= 0.260955; mListOfCoordinates[30][1]=-0.261653; mListOfCoordinates[30][2]= 0.231206;
        mListOfCoordinates[31][0]= 0.150647; mListOfCoordinates[31][1]= 0.188584; mListOfCoordinates[31][2]= 0.013834;
        mListOfCoordinates[32][0]=-0.055068; mListOfCoordinates[32][1]= 0.452083; mListOfCoordinates[32][2]= 0.259645;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.148733;
        mListOfRadii[ 1]= 0.344678;
        mListOfRadii[ 2]= 0.255659;
        mListOfRadii[ 3]= 0.281480;
        mListOfRadii[ 4]= 0.157250;
        mListOfRadii[ 5]= 0.297355;
        mListOfRadii[ 6]= 0.217683;
        mListOfRadii[ 7]= 0.125977;
        mListOfRadii[ 8]= 0.190703;
        mListOfRadii[ 9]= 0.280329;
        mListOfRadii[10]= 0.215576;
        mListOfRadii[11]= 0.372309;
        mListOfRadii[12]= 0.259056;
        mListOfRadii[13]= 0.374545;
        mListOfRadii[14]= 0.339278;
        mListOfRadii[15]= 0.384481;
        mListOfRadii[16]= 0.274594;
        mListOfRadii[17]= 0.226771;
        mListOfRadii[18]= 0.311030;
        mListOfRadii[19]= 0.275784;
        mListOfRadii[20]= 0.309647;
        mListOfRadii[21]= 0.217070;
        mListOfRadii[22]= 0.260999;
        mListOfRadii[23]= 0.300888;
        mListOfRadii[24]= 0.255980;
        mListOfRadii[25]= 0.146572;
        mListOfRadii[26]= 0.195571;
        mListOfRadii[27]= 0.192388;
        mListOfRadii[28]= 0.165185;
        mListOfRadii[29]= 0.275272;
        mListOfRadii[30]= 0.193322;
        mListOfRadii[31]= 0.293912;
        mListOfRadii[32]= 0.039266;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.585744 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.100814;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.224904;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.250009;

    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast2Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast2Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Ballast2Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Ballast2Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

