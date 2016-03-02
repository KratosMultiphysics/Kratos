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
        
        // 1.154643 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.866068,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.866068;
        
        mListOfCoordinates[ 0][0]=-0.125806; mListOfCoordinates[ 0][1]=-0.075738; mListOfCoordinates[ 0][2]= 0.197613;
        mListOfCoordinates[ 1][0]=-0.082787; mListOfCoordinates[ 1][1]=-0.100968; mListOfCoordinates[ 1][2]=-0.167383;
        mListOfCoordinates[ 2][0]=-0.206331; mListOfCoordinates[ 2][1]=-0.217373; mListOfCoordinates[ 2][2]=-0.070579;
        mListOfCoordinates[ 3][0]= 0.063702; mListOfCoordinates[ 3][1]=-0.330957; mListOfCoordinates[ 3][2]=-0.016424;
        mListOfCoordinates[ 4][0]=-0.053811; mListOfCoordinates[ 4][1]= 0.006185; mListOfCoordinates[ 4][2]= 0.305663;
        mListOfCoordinates[ 5][0]= 0.149199; mListOfCoordinates[ 5][1]= 0.194272; mListOfCoordinates[ 5][2]=-0.032231;
        mListOfCoordinates[ 6][0]=-0.065905; mListOfCoordinates[ 6][1]= 0.022058; mListOfCoordinates[ 6][2]= 0.524744;
        mListOfCoordinates[ 7][0]=-0.357305; mListOfCoordinates[ 7][1]=-0.295110; mListOfCoordinates[ 7][2]=-0.215278;
        mListOfCoordinates[ 8][0]=-0.072183; mListOfCoordinates[ 8][1]=-0.181088; mListOfCoordinates[ 8][2]=-0.394628;
        mListOfCoordinates[ 9][0]= 0.079127; mListOfCoordinates[ 9][1]= 0.026088; mListOfCoordinates[ 9][2]=-0.106767;
        mListOfCoordinates[10][0]=-0.173559; mListOfCoordinates[10][1]=-0.047748; mListOfCoordinates[10][2]= 0.566832;
        mListOfCoordinates[11][0]=-0.197243; mListOfCoordinates[11][1]=-0.070396; mListOfCoordinates[11][2]= 0.466903;
        mListOfCoordinates[12][0]=-0.224356; mListOfCoordinates[12][1]=-0.217598; mListOfCoordinates[12][2]=-0.262650;
        mListOfCoordinates[13][0]=-0.060913; mListOfCoordinates[13][1]=-0.335443; mListOfCoordinates[13][2]=-0.006029;
        mListOfCoordinates[14][0]= 0.060267; mListOfCoordinates[14][1]= 0.194358; mListOfCoordinates[14][2]=-0.001223;
        mListOfCoordinates[15][0]=-0.035702; mListOfCoordinates[15][1]= 0.020981; mListOfCoordinates[15][2]= 0.000086;
        mListOfCoordinates[16][0]= 0.241529; mListOfCoordinates[16][1]= 0.472907; mListOfCoordinates[16][2]= 0.245206;
        mListOfCoordinates[17][0]= 0.214458; mListOfCoordinates[17][1]= 0.199109; mListOfCoordinates[17][2]=-0.277018;
        mListOfCoordinates[18][0]=-0.238390; mListOfCoordinates[18][1]=-0.210189; mListOfCoordinates[18][2]= 0.294545;
        mListOfCoordinates[19][0]=-0.371728; mListOfCoordinates[19][1]=-0.319201; mListOfCoordinates[19][2]=-0.142400;
        mListOfCoordinates[20][0]= 0.201401; mListOfCoordinates[20][1]= 0.136546; mListOfCoordinates[20][2]=-0.389773;
        mListOfCoordinates[21][0]= 0.191067; mListOfCoordinates[21][1]=-0.268462; mListOfCoordinates[21][2]=-0.023840;
        mListOfCoordinates[22][0]=-0.065504; mListOfCoordinates[22][1]=-0.082784; mListOfCoordinates[22][2]=-0.078941;
        mListOfCoordinates[23][0]=-0.287312; mListOfCoordinates[23][1]=-0.169747; mListOfCoordinates[23][2]= 0.409287;
        mListOfCoordinates[24][0]= 0.138754; mListOfCoordinates[24][1]= 0.155443; mListOfCoordinates[24][2]=-0.572011;
        mListOfCoordinates[25][0]= 0.126411; mListOfCoordinates[25][1]= 0.018969; mListOfCoordinates[25][2]=-0.341528;
        mListOfCoordinates[26][0]= 0.236838; mListOfCoordinates[26][1]= 0.239500; mListOfCoordinates[26][2]=-0.199352;
        mListOfCoordinates[27][0]= 0.046551; mListOfCoordinates[27][1]=-0.064105; mListOfCoordinates[27][2]=-0.111803;
        mListOfCoordinates[28][0]= 0.060359; mListOfCoordinates[28][1]=-0.118564; mListOfCoordinates[28][2]=-0.351441;

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
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.806013 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.132580;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.166454;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.209359; 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::InitializeSolutionStep(const ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::FinalizeSolutionStep(const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast4Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast4Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Ballast4Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Ballast4Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

