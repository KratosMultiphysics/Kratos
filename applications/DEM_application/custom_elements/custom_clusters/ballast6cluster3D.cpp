//
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Joaquín $
//   Date:                $Date: 2016-01-26 15:00:00 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "ballast6cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast6Cluster3D::Ballast6Cluster3D() : Cluster3D() {}
            
      
    Ballast6Cluster3D::Ballast6Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast6Cluster3D::Ballast6Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast6Cluster3D::Ballast6Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast6Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast6Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast6Cluster3D::~Ballast6Cluster3D() {}
      
    
    void Ballast6Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 44;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 0.989332 (in meters) was the medium diameter of the rock in GiD (grain_0452_geom.gid file)
        // so we should multiply every size that follow by the inverse of that number, 1.010783,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 1.010783;

        mListOfCoordinates[ 0][0]= 0.154526; mListOfCoordinates[ 0][1]=-0.232991; mListOfCoordinates[ 0][2]=-0.357640;
        mListOfCoordinates[ 1][0]= 0.057113; mListOfCoordinates[ 1][1]= 0.176699; mListOfCoordinates[ 1][2]=-0.104017;
        mListOfCoordinates[ 2][0]= 0.122077; mListOfCoordinates[ 2][1]= 0.172146; mListOfCoordinates[ 2][2]=-0.288484;
        mListOfCoordinates[ 3][0]=-0.279310; mListOfCoordinates[ 3][1]= 0.184894; mListOfCoordinates[ 3][2]= 0.335374;
        mListOfCoordinates[ 4][0]= 0.023532; mListOfCoordinates[ 4][1]= 0.000773; mListOfCoordinates[ 4][2]= 0.140640;
        mListOfCoordinates[ 5][0]= 0.043812; mListOfCoordinates[ 5][1]=-0.161539; mListOfCoordinates[ 5][2]= 0.123590;
        mListOfCoordinates[ 6][0]=-0.026180; mListOfCoordinates[ 6][1]= 0.036596; mListOfCoordinates[ 6][2]= 0.142566;
        mListOfCoordinates[ 7][0]= 0.061040; mListOfCoordinates[ 7][1]=-0.128848; mListOfCoordinates[ 7][2]= 0.020946;
        mListOfCoordinates[ 8][0]=-0.228922; mListOfCoordinates[ 8][1]= 0.390636; mListOfCoordinates[ 8][2]= 0.445848;
        mListOfCoordinates[ 9][0]= 0.016493; mListOfCoordinates[ 9][1]= 0.045919; mListOfCoordinates[ 9][2]=-0.468434;
        mListOfCoordinates[10][0]=-0.245628; mListOfCoordinates[10][1]= 0.211156; mListOfCoordinates[10][2]= 0.372786;
        mListOfCoordinates[11][0]= 0.283690; mListOfCoordinates[11][1]=-0.197271; mListOfCoordinates[11][2]=-0.132999;
        mListOfCoordinates[12][0]= 0.109943; mListOfCoordinates[12][1]=-0.169405; mListOfCoordinates[12][2]=-0.446046;
        mListOfCoordinates[13][0]= 0.137166; mListOfCoordinates[13][1]=-0.304142; mListOfCoordinates[13][2]=-0.185309;
        mListOfCoordinates[14][0]=-0.070918; mListOfCoordinates[14][1]= 0.148857; mListOfCoordinates[14][2]= 0.342722;
        mListOfCoordinates[15][0]=-0.246180; mListOfCoordinates[15][1]= 0.424930; mListOfCoordinates[15][2]= 0.309099;
        mListOfCoordinates[16][0]= 0.142235; mListOfCoordinates[16][1]=-0.321863; mListOfCoordinates[16][2]=-0.294089;
        mListOfCoordinates[17][0]= 0.070753; mListOfCoordinates[17][1]=-0.038127; mListOfCoordinates[17][2]=-0.201478;
        mListOfCoordinates[18][0]=-0.201491; mListOfCoordinates[18][1]=-0.028836; mListOfCoordinates[18][2]= 0.201593;
        mListOfCoordinates[19][0]= 0.142080; mListOfCoordinates[19][1]=-0.149696; mListOfCoordinates[19][2]=-0.198547;
        mListOfCoordinates[20][0]=-0.004741; mListOfCoordinates[20][1]=-0.022268; mListOfCoordinates[20][2]=-0.036651;
        mListOfCoordinates[21][0]= 0.109159; mListOfCoordinates[21][1]=-0.255812; mListOfCoordinates[21][2]= 0.120899;
        mListOfCoordinates[22][0]=-0.397062; mListOfCoordinates[22][1]= 0.323006; mListOfCoordinates[22][2]= 0.291901;
        mListOfCoordinates[23][0]=-0.191167; mListOfCoordinates[23][1]=-0.024045; mListOfCoordinates[23][2]=-0.281356;
        mListOfCoordinates[24][0]=-0.231739; mListOfCoordinates[24][1]= 0.083142; mListOfCoordinates[24][2]= 0.269277;
        mListOfCoordinates[25][0]= 0.127492; mListOfCoordinates[25][1]= 0.226223; mListOfCoordinates[25][2]=-0.209853;
        mListOfCoordinates[26][0]= 0.107985; mListOfCoordinates[26][1]=-0.311237; mListOfCoordinates[26][2]=-0.036783;
        mListOfCoordinates[27][0]=-0.176574; mListOfCoordinates[27][1]=-0.230671; mListOfCoordinates[27][2]= 0.189773;
        mListOfCoordinates[28][0]=-0.095595; mListOfCoordinates[28][1]= 0.221152; mListOfCoordinates[28][2]= 0.463295;
        mListOfCoordinates[29][0]=-0.161136; mListOfCoordinates[29][1]= 0.377894; mListOfCoordinates[29][2]= 0.227354;
        mListOfCoordinates[30][0]=-0.119950; mListOfCoordinates[30][1]= 0.001862; mListOfCoordinates[30][2]=-0.402085;
        mListOfCoordinates[31][0]= 0.015778; mListOfCoordinates[31][1]= 0.322661; mListOfCoordinates[31][2]= 0.005045;
        mListOfCoordinates[32][0]=-0.204440; mListOfCoordinates[32][1]= 0.300092; mListOfCoordinates[32][2]= 0.477911;
        mListOfCoordinates[33][0]=-0.213493; mListOfCoordinates[33][1]=-0.055403; mListOfCoordinates[33][2]=-0.085896;
        mListOfCoordinates[34][0]= 0.042008; mListOfCoordinates[34][1]=-0.016415; mListOfCoordinates[34][2]=-0.255337;
        mListOfCoordinates[35][0]= 0.036951; mListOfCoordinates[35][1]=-0.060371; mListOfCoordinates[35][2]=-0.481100;
        mListOfCoordinates[36][0]=-0.258692; mListOfCoordinates[36][1]=-0.181086; mListOfCoordinates[36][2]= 0.214894;
        mListOfCoordinates[37][0]=-0.295765; mListOfCoordinates[37][1]=-0.049771; mListOfCoordinates[37][2]= 0.076606;
        mListOfCoordinates[38][0]= 0.118710; mListOfCoordinates[38][1]= 0.179190; mListOfCoordinates[38][2]=-0.454814;
        mListOfCoordinates[39][0]=-0.067306; mListOfCoordinates[39][1]=-0.318055; mListOfCoordinates[39][2]= 0.171100;
        mListOfCoordinates[40][0]=-0.063380; mListOfCoordinates[40][1]= 0.320196; mListOfCoordinates[40][2]= 0.133387;
        mListOfCoordinates[41][0]=-0.237025; mListOfCoordinates[41][1]= 0.267292; mListOfCoordinates[41][2]= 0.302601;
        mListOfCoordinates[42][0]= 0.255741; mListOfCoordinates[42][1]=-0.230471; mListOfCoordinates[42][2]=-0.005797;
        mListOfCoordinates[43][0]=-0.280429; mListOfCoordinates[43][1]=-0.036772; mListOfCoordinates[43][2]=-0.225366; 

        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.187781;
        mListOfRadii[ 1]= 0.264401;
        mListOfRadii[ 2]= 0.207506;
        mListOfRadii[ 3]= 0.179643;
        mListOfRadii[ 4]= 0.343981;
        mListOfRadii[ 5]= 0.270950;
        mListOfRadii[ 6]= 0.377743;
        mListOfRadii[ 7]= 0.326901;
        mListOfRadii[ 8]= 0.065563;
        mListOfRadii[ 9]= 0.226106;
        mListOfRadii[10]= 0.186755;
        mListOfRadii[11]= 0.161458;
        mListOfRadii[12]= 0.165767;
        mListOfRadii[13]= 0.163999;
        mListOfRadii[14]= 0.234043;
        mListOfRadii[15]= 0.137725;
        mListOfRadii[16]= 0.123157;
        mListOfRadii[17]= 0.315835;
        mListOfRadii[18]= 0.242484;
        mListOfRadii[19]= 0.274473;
        mListOfRadii[20]= 0.347929;
        mListOfRadii[21]= 0.271708;
        mListOfRadii[22]= 0.040709;
        mListOfRadii[23]= 0.130733;
        mListOfRadii[24]= 0.211292;
        mListOfRadii[25]= 0.206535;
        mListOfRadii[26]= 0.201208;
        mListOfRadii[27]= 0.169170;
        mListOfRadii[28]= 0.135061;
        mListOfRadii[29]= 0.155623;
        mListOfRadii[30]= 0.153920;
        mListOfRadii[31]= 0.140006;
        mListOfRadii[32]= 0.123088;
        mListOfRadii[33]= 0.201459;
        mListOfRadii[34]= 0.289651;
        mListOfRadii[35]= 0.192456;
        mListOfRadii[36]= 0.144829;
        mListOfRadii[37]= 0.139493;
        mListOfRadii[38]= 0.133983;
        mListOfRadii[39]= 0.164085;
        mListOfRadii[40]= 0.171098;
        mListOfRadii[41]= 0.212894;
        mListOfRadii[42]= 0.168465;
        mListOfRadii[43]= 0.060335;
        
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.507020 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.092287;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.191426;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.223336;

    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast6Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Ballast6Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Ballast6Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

