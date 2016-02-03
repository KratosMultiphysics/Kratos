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
#include "ballast4cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast4Cluster3D::Ballast4Cluster3D() : Cluster3D() {}
            
      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast4Cluster3D::Ballast4Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast4Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast4Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast4Cluster3D::~Ballast4Cluster3D() {}
      
    
    void Ballast4Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 29;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.183 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.8453,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.8453;
        
        mListOfCoordinates[ 0][0]=-0.097904; mListOfCoordinates[ 0][1]=-0.026805; mListOfCoordinates[ 0][2]= 0.223489;
        mListOfCoordinates[ 1][0]=-0.054885; mListOfCoordinates[ 1][1]=-0.052035; mListOfCoordinates[ 1][2]=-0.141507;
        mListOfCoordinates[ 2][0]=-0.178429; mListOfCoordinates[ 2][1]=-0.168440; mListOfCoordinates[ 2][2]=-0.044703;
        mListOfCoordinates[ 3][0]= 0.091604; mListOfCoordinates[ 3][1]=-0.282024; mListOfCoordinates[ 3][2]= 0.009452;
        mListOfCoordinates[ 4][0]=-0.025909; mListOfCoordinates[ 4][1]= 0.055118; mListOfCoordinates[ 4][2]= 0.331539;
        mListOfCoordinates[ 5][0]= 0.177101; mListOfCoordinates[ 5][1]= 0.243205; mListOfCoordinates[ 5][2]=-0.006355;
        mListOfCoordinates[ 6][0]=-0.038003; mListOfCoordinates[ 6][1]= 0.070991; mListOfCoordinates[ 6][2]= 0.550620;
        mListOfCoordinates[ 7][0]=-0.329403; mListOfCoordinates[ 7][1]=-0.246177; mListOfCoordinates[ 7][2]=-0.189402;
        mListOfCoordinates[ 8][0]=-0.044281; mListOfCoordinates[ 8][1]=-0.132155; mListOfCoordinates[ 8][2]=-0.368752;
        mListOfCoordinates[ 9][0]= 0.107029; mListOfCoordinates[ 9][1]= 0.075021; mListOfCoordinates[ 9][2]=-0.080891;
        mListOfCoordinates[10][0]=-0.145657; mListOfCoordinates[10][1]= 0.001185; mListOfCoordinates[10][2]= 0.592708;
        mListOfCoordinates[11][0]=-0.169341; mListOfCoordinates[11][1]=-0.021463; mListOfCoordinates[11][2]= 0.492779;
        mListOfCoordinates[12][0]=-0.196454; mListOfCoordinates[12][1]=-0.168665; mListOfCoordinates[12][2]=-0.236774;
        mListOfCoordinates[13][0]=-0.033011; mListOfCoordinates[13][1]=-0.286510; mListOfCoordinates[13][2]= 0.019847;
        mListOfCoordinates[14][0]= 0.088169; mListOfCoordinates[14][1]= 0.243291; mListOfCoordinates[14][2]= 0.024653;
        mListOfCoordinates[15][0]=-0.007800; mListOfCoordinates[15][1]= 0.069914; mListOfCoordinates[15][2]= 0.025962;
        mListOfCoordinates[16][0]= 0.269431; mListOfCoordinates[16][1]= 0.521840; mListOfCoordinates[16][2]= 0.271082;
        mListOfCoordinates[17][0]= 0.242360; mListOfCoordinates[17][1]= 0.248042; mListOfCoordinates[17][2]=-0.251142;
        mListOfCoordinates[18][0]=-0.210488; mListOfCoordinates[18][1]=-0.161256; mListOfCoordinates[18][2]= 0.320421;
        mListOfCoordinates[19][0]=-0.343826; mListOfCoordinates[19][1]=-0.270268; mListOfCoordinates[19][2]=-0.116524;
        mListOfCoordinates[20][0]= 0.229303; mListOfCoordinates[20][1]= 0.185479; mListOfCoordinates[20][2]=-0.363897;
        mListOfCoordinates[21][0]= 0.218969; mListOfCoordinates[21][1]=-0.219529; mListOfCoordinates[21][2]= 0.002036;
        mListOfCoordinates[22][0]=-0.037602; mListOfCoordinates[22][1]=-0.033851; mListOfCoordinates[22][2]=-0.053065;
        mListOfCoordinates[23][0]=-0.259410; mListOfCoordinates[23][1]=-0.120814; mListOfCoordinates[23][2]= 0.435163;
        mListOfCoordinates[24][0]= 0.166656; mListOfCoordinates[24][1]= 0.204376; mListOfCoordinates[24][2]=-0.546135;
        mListOfCoordinates[25][0]= 0.154313; mListOfCoordinates[25][1]= 0.067902; mListOfCoordinates[25][2]=-0.315652;
        mListOfCoordinates[26][0]= 0.264740; mListOfCoordinates[26][1]= 0.288433; mListOfCoordinates[26][2]=-0.173476;
        mListOfCoordinates[27][0]= 0.074453; mListOfCoordinates[27][1]=-0.015172; mListOfCoordinates[27][2]=-0.085927;
        mListOfCoordinates[28][0]= 0.088261; mListOfCoordinates[28][1]=-0.069631; mListOfCoordinates[28][2]=-0.325565;

        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.343889;
        mListOfRadii[ 1]= 0.343504;
        mListOfRadii[ 2]= 0.375064;
        mListOfRadii[ 3]= 0.227838;
        mListOfRadii[ 4]= 0.301032;
        mListOfRadii[ 5]= 0.433485;
        mListOfRadii[ 6]= 0.232885;
        mListOfRadii[ 7]= 0.201071;
        mListOfRadii[ 8]= 0.175775;
        mListOfRadii[ 9]= 0.477571;
        mListOfRadii[10]= 0.209730;
        mListOfRadii[11]= 0.250653;
        mListOfRadii[12]= 0.222536;
        mListOfRadii[13]= 0.244465;
        mListOfRadii[14]= 0.421175;
        mListOfRadii[15]= 0.414231;
        mListOfRadii[16]= 0.047803;
        mListOfRadii[17]= 0.258038;
        mListOfRadii[18]= 0.263474;
        mListOfRadii[19]= 0.203770;
        mListOfRadii[20]= 0.206013;
        mListOfRadii[21]= 0.251268;
        mListOfRadii[22]= 0.392997;
        mListOfRadii[23]= 0.207835;
        mListOfRadii[24]= 0.033922;
        mListOfRadii[25]= 0.296245;
        mListOfRadii[26]= 0.285736;
        mListOfRadii[27]= 0.407137;
        mListOfRadii[28]= 0.237066;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.662181 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.106204;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.124578;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.155769; 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast4Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast4Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast4Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

