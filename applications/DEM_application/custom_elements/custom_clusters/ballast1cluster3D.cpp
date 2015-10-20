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
#include "ballast1cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast1Cluster3D::Ballast1Cluster3D() : Cluster3D() {}
            
      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast1Cluster3D::Ballast1Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast1Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast1Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast1Cluster3D::~Ballast1Cluster3D() {}
      
    
    void Ballast1Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 34;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 2.5900227 (in meters) was the medium diameter of the rock in GiD (Rock1_01.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.386097,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.9879;

        mListOfCoordinates[ 0][0]=-0.239249; mListOfCoordinates[ 0][1]= 0.135066; mListOfCoordinates[ 0][2]= 0.008573;
        mListOfCoordinates[ 1][0]= 0.061018; mListOfCoordinates[ 1][1]= 0.231957; mListOfCoordinates[ 1][2]= 0.013032;
        mListOfCoordinates[ 2][0]= 0.034042; mListOfCoordinates[ 2][1]=-0.218863; mListOfCoordinates[ 2][2]=-0.005182;
        mListOfCoordinates[ 3][0]=-0.120326; mListOfCoordinates[ 3][1]= 0.294719; mListOfCoordinates[ 3][2]= 0.022409;
        mListOfCoordinates[ 4][0]= 0.197621; mListOfCoordinates[ 4][1]=-0.138049; mListOfCoordinates[ 4][2]=-0.006054;
        mListOfCoordinates[ 5][0]=-0.215914; mListOfCoordinates[ 5][1]= 0.046352; mListOfCoordinates[ 5][2]=-0.002806;
        mListOfCoordinates[ 6][0]=-0.398867; mListOfCoordinates[ 6][1]=-0.173338; mListOfCoordinates[ 6][2]=-0.001860;
        mListOfCoordinates[ 7][0]= 0.615750; mListOfCoordinates[ 7][1]=-0.016372; mListOfCoordinates[ 7][2]=-0.005601;
        mListOfCoordinates[ 8][0]=-0.375816; mListOfCoordinates[ 8][1]= 0.171554; mListOfCoordinates[ 8][2]= 0.011697;
        mListOfCoordinates[ 9][0]=-0.324130; mListOfCoordinates[ 9][1]= 0.304857; mListOfCoordinates[ 9][2]= 0.001506;
        mListOfCoordinates[10][0]=-0.563655; mListOfCoordinates[10][1]= 0.166050; mListOfCoordinates[10][2]=-0.013256;
        mListOfCoordinates[11][0]=-0.609904; mListOfCoordinates[11][1]=-0.051607; mListOfCoordinates[11][2]=-0.012257;
        mListOfCoordinates[12][0]= 0.555884; mListOfCoordinates[12][1]= 0.140029; mListOfCoordinates[12][2]= 0.004157;
        mListOfCoordinates[13][0]=-0.031005; mListOfCoordinates[13][1]= 0.265244; mListOfCoordinates[13][2]= 0.008116;
        mListOfCoordinates[14][0]= 0.422827; mListOfCoordinates[14][1]=-0.102716; mListOfCoordinates[14][2]=-0.008833;
        mListOfCoordinates[15][0]=-0.235162; mListOfCoordinates[15][1]=-0.239884; mListOfCoordinates[15][2]=-0.001202;
        mListOfCoordinates[16][0]= 0.305764; mListOfCoordinates[16][1]= 0.214467; mListOfCoordinates[16][2]= 0.024351;
        mListOfCoordinates[17][0]=-0.371884; mListOfCoordinates[17][1]=-0.041089; mListOfCoordinates[17][2]=-0.000486;
        mListOfCoordinates[18][0]= 0.158129; mListOfCoordinates[18][1]=-0.180230; mListOfCoordinates[18][2]= 0.001457;
        mListOfCoordinates[19][0]=-0.212650; mListOfCoordinates[19][1]= 0.257694; mListOfCoordinates[19][2]= 0.014582;
        mListOfCoordinates[20][0]=-0.050354; mListOfCoordinates[20][1]= 0.007217; mListOfCoordinates[20][2]=-0.012359;
        mListOfCoordinates[21][0]= 0.197770; mListOfCoordinates[21][1]= 0.055864; mListOfCoordinates[21][2]= 0.081300;
        mListOfCoordinates[22][0]=-0.508956; mListOfCoordinates[22][1]= 0.022799; mListOfCoordinates[22][2]= 0.010677;
        mListOfCoordinates[23][0]=-0.149272; mListOfCoordinates[23][1]=-0.061849; mListOfCoordinates[23][2]=-0.016144;
        mListOfCoordinates[24][0]=-0.421250; mListOfCoordinates[24][1]= 0.180733; mListOfCoordinates[24][2]= 0.011210;
        mListOfCoordinates[25][0]=-0.131319; mListOfCoordinates[25][1]=-0.207744; mListOfCoordinates[25][2]= 0.003377;
        mListOfCoordinates[26][0]=-0.067746; mListOfCoordinates[26][1]= 0.149392; mListOfCoordinates[26][2]= 0.003433;
        mListOfCoordinates[27][0]=-0.156173; mListOfCoordinates[27][1]=-0.210251; mListOfCoordinates[27][2]=-0.002048;
        mListOfCoordinates[28][0]=-0.064504; mListOfCoordinates[28][1]=-0.079852; mListOfCoordinates[28][2]=-0.015108;
        mListOfCoordinates[29][0]= 0.253230; mListOfCoordinates[29][1]= 0.126287; mListOfCoordinates[29][2]= 0.012017;
        mListOfCoordinates[30][0]= 0.006384; mListOfCoordinates[30][1]=-0.073943; mListOfCoordinates[30][2]=-0.111380;
        mListOfCoordinates[31][0]= 0.367258; mListOfCoordinates[31][1]=-0.032396; mListOfCoordinates[31][2]=-0.004700;
        mListOfCoordinates[32][0]= 0.070691; mListOfCoordinates[32][1]= 0.069217; mListOfCoordinates[32][2]= 0.012072;
        mListOfCoordinates[33][0]= 0.374490; mListOfCoordinates[33][1]= 0.076991; mListOfCoordinates[33][2]= 0.001279;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.199183;
        mListOfRadii[ 1]= 0.218745;
        mListOfRadii[ 2]= 0.255778;
        mListOfRadii[ 3]= 0.184401;
        mListOfRadii[ 4]= 0.268336;
        mListOfRadii[ 5]= 0.209218;
        mListOfRadii[ 6]= 0.206688;
        mListOfRadii[ 7]= 0.184994;
        mListOfRadii[ 8]= 0.183692;
        mListOfRadii[ 9]= 0.194576;
        mListOfRadii[10]= 0.166181;
        mListOfRadii[11]= 0.165078;
        mListOfRadii[12]= 0.218567;
        mListOfRadii[13]= 0.181078;
        mListOfRadii[14]= 0.212356;
        mListOfRadii[15]= 0.236673;
        mListOfRadii[16]= 0.142964;
        mListOfRadii[17]= 0.192695;
        mListOfRadii[18]= 0.197730;
        mListOfRadii[19]= 0.191680;
        mListOfRadii[20]= 0.202352;
        mListOfRadii[21]= 0.100967;
        mListOfRadii[22]= 0.169001;
        mListOfRadii[23]= 0.202381;
        mListOfRadii[24]= 0.200122;
        mListOfRadii[25]= 0.200207;
        mListOfRadii[26]= 0.212194;
        mListOfRadii[27]= 0.194278;
        mListOfRadii[28]= 0.192909;
        mListOfRadii[29]= 0.169839;
        mListOfRadii[30]= 0.103225;
        mListOfRadii[31]= 0.195763;
        mListOfRadii[32]= 0.207302;
        mListOfRadii[33]= 0.191825;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.386097 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.059488566;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.145587554;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.181080732;

    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast1Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast1Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast1Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast1Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

