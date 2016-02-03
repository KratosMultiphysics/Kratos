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
#include "includes/define.h"
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
        
        // 1.33944 (in meters) was the medium diameter of the rock in GiD (grain_0452_geom.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.746581,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.746581;

        mListOfCoordinates[ 0][0]= 0.202681; mListOfCoordinates[ 0][1]=-0.253223; mListOfCoordinates[ 0][2]=-0.362802;
        mListOfCoordinates[ 1][0]= 0.105268; mListOfCoordinates[ 1][1]= 0.156467; mListOfCoordinates[ 1][2]=-0.109179;
        mListOfCoordinates[ 2][0]= 0.170232; mListOfCoordinates[ 2][1]= 0.151914; mListOfCoordinates[ 2][2]=-0.293646;
        mListOfCoordinates[ 3][0]=-0.231155; mListOfCoordinates[ 3][1]= 0.164662; mListOfCoordinates[ 3][2]= 0.330212;
        mListOfCoordinates[ 4][0]= 0.071687; mListOfCoordinates[ 4][1]=-0.019459; mListOfCoordinates[ 4][2]= 0.135478;
        mListOfCoordinates[ 5][0]= 0.091967; mListOfCoordinates[ 5][1]=-0.181771; mListOfCoordinates[ 5][2]= 0.118428;
        mListOfCoordinates[ 6][0]= 0.021975; mListOfCoordinates[ 6][1]= 0.016364; mListOfCoordinates[ 6][2]= 0.137404;
        mListOfCoordinates[ 7][0]= 0.109195; mListOfCoordinates[ 7][1]=-0.149080; mListOfCoordinates[ 7][2]= 0.015784;
        mListOfCoordinates[ 8][0]=-0.180767; mListOfCoordinates[ 8][1]= 0.370404; mListOfCoordinates[ 8][2]= 0.440686;
        mListOfCoordinates[ 9][0]= 0.064648; mListOfCoordinates[ 9][1]= 0.025687; mListOfCoordinates[ 9][2]=-0.473596;
        mListOfCoordinates[10][0]=-0.197473; mListOfCoordinates[10][1]= 0.190924; mListOfCoordinates[10][2]= 0.367624;
        mListOfCoordinates[11][0]= 0.331845; mListOfCoordinates[11][1]=-0.217503; mListOfCoordinates[11][2]=-0.138161;
        mListOfCoordinates[12][0]= 0.158098; mListOfCoordinates[12][1]=-0.189637; mListOfCoordinates[12][2]=-0.451208;
        mListOfCoordinates[13][0]= 0.185321; mListOfCoordinates[13][1]=-0.324374; mListOfCoordinates[13][2]=-0.190471;
        mListOfCoordinates[14][0]=-0.022763; mListOfCoordinates[14][1]= 0.128625; mListOfCoordinates[14][2]= 0.337560;
        mListOfCoordinates[15][0]=-0.198025; mListOfCoordinates[15][1]= 0.404698; mListOfCoordinates[15][2]= 0.303937;
        mListOfCoordinates[16][0]= 0.190390; mListOfCoordinates[16][1]=-0.342095; mListOfCoordinates[16][2]=-0.299251;
        mListOfCoordinates[17][0]= 0.118908; mListOfCoordinates[17][1]=-0.058359; mListOfCoordinates[17][2]=-0.206640;
        mListOfCoordinates[18][0]=-0.153336; mListOfCoordinates[18][1]=-0.049068; mListOfCoordinates[18][2]= 0.196431;
        mListOfCoordinates[19][0]= 0.190235; mListOfCoordinates[19][1]=-0.169928; mListOfCoordinates[19][2]=-0.203709;
        mListOfCoordinates[20][0]= 0.043414; mListOfCoordinates[20][1]=-0.042500; mListOfCoordinates[20][2]=-0.041813;
        mListOfCoordinates[21][0]= 0.157314; mListOfCoordinates[21][1]=-0.276044; mListOfCoordinates[21][2]= 0.115737;
        mListOfCoordinates[22][0]=-0.348907; mListOfCoordinates[22][1]= 0.302774; mListOfCoordinates[22][2]= 0.286739;
        mListOfCoordinates[23][0]=-0.143012; mListOfCoordinates[23][1]=-0.044277; mListOfCoordinates[23][2]=-0.286518;
        mListOfCoordinates[24][0]=-0.183584; mListOfCoordinates[24][1]= 0.062910; mListOfCoordinates[24][2]= 0.264115;
        mListOfCoordinates[25][0]= 0.175647; mListOfCoordinates[25][1]= 0.205991; mListOfCoordinates[25][2]=-0.215015;
        mListOfCoordinates[26][0]= 0.156140; mListOfCoordinates[26][1]=-0.331469; mListOfCoordinates[26][2]=-0.041945;
        mListOfCoordinates[27][0]=-0.128419; mListOfCoordinates[27][1]=-0.250903; mListOfCoordinates[27][2]= 0.184611;
        mListOfCoordinates[28][0]=-0.047440; mListOfCoordinates[28][1]= 0.200920; mListOfCoordinates[28][2]= 0.458133;
        mListOfCoordinates[29][0]=-0.112981; mListOfCoordinates[29][1]= 0.357662; mListOfCoordinates[29][2]= 0.222192;
        mListOfCoordinates[30][0]=-0.071795; mListOfCoordinates[30][1]=-0.018370; mListOfCoordinates[30][2]=-0.407247;
        mListOfCoordinates[31][0]= 0.063933; mListOfCoordinates[31][1]= 0.302429; mListOfCoordinates[31][2]=-0.000117;
        mListOfCoordinates[32][0]=-0.156285; mListOfCoordinates[32][1]= 0.279860; mListOfCoordinates[32][2]= 0.472749;
        mListOfCoordinates[33][0]=-0.165338; mListOfCoordinates[33][1]=-0.075635; mListOfCoordinates[33][2]=-0.091058;
        mListOfCoordinates[34][0]= 0.090163; mListOfCoordinates[34][1]=-0.036647; mListOfCoordinates[34][2]=-0.260499;
        mListOfCoordinates[35][0]= 0.085106; mListOfCoordinates[35][1]=-0.080603; mListOfCoordinates[35][2]=-0.486262;
        mListOfCoordinates[36][0]=-0.210537; mListOfCoordinates[36][1]=-0.201318; mListOfCoordinates[36][2]= 0.209732;
        mListOfCoordinates[37][0]=-0.247610; mListOfCoordinates[37][1]=-0.070003; mListOfCoordinates[37][2]= 0.071444;
        mListOfCoordinates[38][0]= 0.166865; mListOfCoordinates[38][1]= 0.158958; mListOfCoordinates[38][2]=-0.459976;
        mListOfCoordinates[39][0]=-0.019151; mListOfCoordinates[39][1]=-0.338287; mListOfCoordinates[39][2]= 0.165938;
        mListOfCoordinates[40][0]=-0.015225; mListOfCoordinates[40][1]= 0.299964; mListOfCoordinates[40][2]= 0.128225;
        mListOfCoordinates[41][0]=-0.188870; mListOfCoordinates[41][1]= 0.247060; mListOfCoordinates[41][2]= 0.297439;
        mListOfCoordinates[42][0]= 0.303896; mListOfCoordinates[42][1]=-0.250703; mListOfCoordinates[42][2]=-0.010959;
        mListOfCoordinates[43][0]=-0.232274; mListOfCoordinates[43][1]=-0.057004; mListOfCoordinates[43][2]=-0.230528; 

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
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.42385 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.070320;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.109187;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.127069;

    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast6Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast6Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast6Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast6Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

