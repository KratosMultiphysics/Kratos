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
#include "ballast5cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Ballast5Cluster3D::Ballast5Cluster3D() : Cluster3D() {}
            
      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Ballast5Cluster3D::Ballast5Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Ballast5Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Ballast5Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Ballast5Cluster3D::~Ballast5Cluster3D() {}
      
    
    void Ballast5Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 45;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 0.834368 (in meters) was the medium diameter of the rock in GiD (grain_0451_geom.gid file)
        // so we should multiply every size that follow by the inverse of that number, 1.198512,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 1.198512;

        mListOfCoordinates[ 0][0]= 0.122446; mListOfCoordinates[ 0][1]= 0.054697; mListOfCoordinates[ 0][2]= 0.314590;
        mListOfCoordinates[ 1][0]= 0.006550; mListOfCoordinates[ 1][1]=-0.006551; mListOfCoordinates[ 1][2]= 0.203957;
        mListOfCoordinates[ 2][0]=-0.335715; mListOfCoordinates[ 2][1]=-0.099698; mListOfCoordinates[ 2][2]= 0.222338;
        mListOfCoordinates[ 3][0]= 0.062857; mListOfCoordinates[ 3][1]= 0.087804; mListOfCoordinates[ 3][2]= 0.490972;
        mListOfCoordinates[ 4][0]= 0.384461; mListOfCoordinates[ 4][1]=-0.131279; mListOfCoordinates[ 4][2]= 0.098535;
        mListOfCoordinates[ 5][0]=-0.049615; mListOfCoordinates[ 5][1]=-0.004342; mListOfCoordinates[ 5][2]= 0.487724;
        mListOfCoordinates[ 6][0]=-0.199512; mListOfCoordinates[ 6][1]= 0.035878; mListOfCoordinates[ 6][2]=-0.473847;
        mListOfCoordinates[ 7][0]=-0.204450; mListOfCoordinates[ 7][1]=-0.049991; mListOfCoordinates[ 7][2]= 0.509979;
        mListOfCoordinates[ 8][0]=-0.104326; mListOfCoordinates[ 8][1]= 0.045843; mListOfCoordinates[ 8][2]=-0.580824;
        mListOfCoordinates[ 9][0]= 0.018122; mListOfCoordinates[ 9][1]= 0.045995; mListOfCoordinates[ 9][2]=-0.443339;
        mListOfCoordinates[10][0]= 0.111737; mListOfCoordinates[10][1]= 0.027518; mListOfCoordinates[10][2]=-0.195233;
        mListOfCoordinates[11][0]= 0.089208; mListOfCoordinates[11][1]= 0.071456; mListOfCoordinates[11][2]= 0.159663;
        mListOfCoordinates[12][0]=-0.032996; mListOfCoordinates[12][1]= 0.055237; mListOfCoordinates[12][2]=-0.253771;
        mListOfCoordinates[13][0]=-0.253628; mListOfCoordinates[13][1]=-0.012385; mListOfCoordinates[13][2]=-0.042384;
        mListOfCoordinates[14][0]= 0.394031; mListOfCoordinates[14][1]=-0.140282; mListOfCoordinates[14][2]=-0.036371;
        mListOfCoordinates[15][0]=-0.159027; mListOfCoordinates[15][1]=-0.046563; mListOfCoordinates[15][2]= 0.340039;
        mListOfCoordinates[16][0]=-0.022858; mListOfCoordinates[16][1]=-0.081712; mListOfCoordinates[16][2]= 0.318210;
        mListOfCoordinates[17][0]=-0.313193; mListOfCoordinates[17][1]=-0.099407; mListOfCoordinates[17][2]= 0.448613;
        mListOfCoordinates[18][0]=-0.155450; mListOfCoordinates[18][1]=-0.012320; mListOfCoordinates[18][2]= 0.076087;
        mListOfCoordinates[19][0]=-0.170394; mListOfCoordinates[19][1]=-0.047290; mListOfCoordinates[19][2]= 0.272506;
        mListOfCoordinates[20][0]=-0.358238; mListOfCoordinates[20][1]=-0.140202; mListOfCoordinates[20][2]= 0.313960;
        mListOfCoordinates[21][0]= 0.005269; mListOfCoordinates[21][1]=-0.032401; mListOfCoordinates[21][2]=-0.584301;
        mListOfCoordinates[22][0]= 0.299964; mListOfCoordinates[22][1]= 0.034193; mListOfCoordinates[22][2]= 0.413795;
        mListOfCoordinates[23][0]= 0.347250; mListOfCoordinates[23][1]=-0.159755; mListOfCoordinates[23][2]= 0.207182;
        mListOfCoordinates[24][0]=-0.150415; mListOfCoordinates[24][1]= 0.042863; mListOfCoordinates[24][2]=-0.380773;
        mListOfCoordinates[25][0]= 0.269520; mListOfCoordinates[25][1]=-0.073773; mListOfCoordinates[25][2]= 0.310565;
        mListOfCoordinates[26][0]= 0.022550; mListOfCoordinates[26][1]= 0.036801; mListOfCoordinates[26][2]=-0.048536;
        mListOfCoordinates[27][0]=-0.294539; mListOfCoordinates[27][1]=-0.039049; mListOfCoordinates[27][2]= 0.097164;
        mListOfCoordinates[28][0]= 0.266533; mListOfCoordinates[28][1]=-0.197350; mListOfCoordinates[28][2]= 0.273966;
        mListOfCoordinates[29][0]=-0.151871; mListOfCoordinates[29][1]= 0.051728; mListOfCoordinates[29][2]=-0.268199;
        mListOfCoordinates[30][0]= 0.226666; mListOfCoordinates[30][1]=-0.074558; mListOfCoordinates[30][2]= 0.064190;
        mListOfCoordinates[31][0]= 0.414736; mListOfCoordinates[31][1]=-0.165577; mListOfCoordinates[31][2]=-0.144017;
        mListOfCoordinates[32][0]= 0.283945; mListOfCoordinates[32][1]=-0.095950; mListOfCoordinates[32][2]=-0.097269;
        mListOfCoordinates[33][0]= 0.121527; mListOfCoordinates[33][1]= 0.040626; mListOfCoordinates[33][2]=-0.347463;
        mListOfCoordinates[34][0]= 0.050257; mListOfCoordinates[34][1]=-0.202142; mListOfCoordinates[34][2]= 0.333495;
        mListOfCoordinates[35][0]= 0.126594; mListOfCoordinates[35][1]= 0.018129; mListOfCoordinates[35][2]=-0.054145;
        mListOfCoordinates[36][0]= 0.219640; mListOfCoordinates[36][1]=-0.076416; mListOfCoordinates[36][2]= 0.200113;
        mListOfCoordinates[37][0]= 0.230678; mListOfCoordinates[37][1]=-0.063279; mListOfCoordinates[37][2]=-0.385387;
        mListOfCoordinates[38][0]= 0.351292; mListOfCoordinates[38][1]=-0.102385; mListOfCoordinates[38][2]=-0.323064;
        mListOfCoordinates[39][0]= 0.090228; mListOfCoordinates[39][1]=-0.051419; mListOfCoordinates[39][2]=-0.495200;
        mListOfCoordinates[40][0]= 0.267021; mListOfCoordinates[40][1]=-0.068394; mListOfCoordinates[40][2]=-0.274295;
        mListOfCoordinates[41][0]=-0.158768; mListOfCoordinates[41][1]= 0.023457; mListOfCoordinates[41][2]=-0.136838;
        mListOfCoordinates[42][0]= 0.108756; mListOfCoordinates[42][1]=-0.168985; mListOfCoordinates[42][2]= 0.281052;
        mListOfCoordinates[43][0]=-0.172632; mListOfCoordinates[43][1]=-0.026780; mListOfCoordinates[43][2]= 0.598262;
        mListOfCoordinates[44][0]= 0.388411; mListOfCoordinates[44][1]=-0.132593; mListOfCoordinates[44][2]=-0.233620;

        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.204898;
        mListOfRadii[ 1]= 0.248121;
        mListOfRadii[ 2]= 0.129997;
        mListOfRadii[ 3]= 0.110971;
        mListOfRadii[ 4]= 0.101141;
        mListOfRadii[ 5]= 0.183270;
        mListOfRadii[ 6]= 0.076334;
        mListOfRadii[ 7]= 0.126306;
        mListOfRadii[ 8]= 0.158541;
        mListOfRadii[ 9]= 0.203504;
        mListOfRadii[10]= 0.211509;
        mListOfRadii[11]= 0.207299;
        mListOfRadii[12]= 0.154789;
        mListOfRadii[13]= 0.084444;
        mListOfRadii[14]= 0.090622;
        mListOfRadii[15]= 0.205679;
        mListOfRadii[16]= 0.219614;
        mListOfRadii[17]= 0.096091;
        mListOfRadii[18]= 0.156978;
        mListOfRadii[19]= 0.233531;
        mListOfRadii[20]= 0.094059;
        mListOfRadii[21]= 0.101306;
        mListOfRadii[22]= 0.049709;
        mListOfRadii[23]= 0.116412;
        mListOfRadii[24]= 0.113252;
        mListOfRadii[25]= 0.162797;
        mListOfRadii[26]= 0.213357;
        mListOfRadii[27]= 0.124343;
        mListOfRadii[28]= 0.139851;
        mListOfRadii[29]= 0.097662;
        mListOfRadii[30]= 0.168403;
        mListOfRadii[31]= 0.077644;
        mListOfRadii[32]= 0.150284;
        mListOfRadii[33]= 0.199167;
        mListOfRadii[34]= 0.124836;
        mListOfRadii[35]= 0.209622;
        mListOfRadii[36]= 0.186548;
        mListOfRadii[37]= 0.138241;
        mListOfRadii[38]= 0.085285;
        mListOfRadii[39]= 0.130440;
        mListOfRadii[40]= 0.128303;
        mListOfRadii[41]= 0.112418;
        mListOfRadii[42]= 0.153510;
        mListOfRadii[43]= 0.114211;
        mListOfRadii[44]= 0.073588;
        
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.304139 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.133008;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.242195;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.299279;

    }
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Ballast5Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Ballast5Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
    void Ballast5Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}
    double Ballast5Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

