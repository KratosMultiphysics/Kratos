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
        
        int number_of_spheres = 37;
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
     
        mListOfCoordinates[ 0][0]= 0.093055; mListOfCoordinates[ 0][1]= 0.407802; mListOfCoordinates[ 0][2]=-0.079434;
        mListOfCoordinates[ 1][0]= 0.220836; mListOfCoordinates[ 1][1]=-0.141754; mListOfCoordinates[ 1][2]= 0.143634;
        mListOfCoordinates[ 2][0]= 0.260100; mListOfCoordinates[ 2][1]=-0.255984; mListOfCoordinates[ 2][2]= 0.223026;
        mListOfCoordinates[ 3][0]= 0.108624; mListOfCoordinates[ 3][1]=-0.049335; mListOfCoordinates[ 3][2]= 0.044379;
        mListOfCoordinates[ 4][0]= 0.106348; mListOfCoordinates[ 4][1]= 0.322341; mListOfCoordinates[ 4][2]= 0.011322;
        mListOfCoordinates[ 5][0]=-0.320713; mListOfCoordinates[ 5][1]=-0.211575; mListOfCoordinates[ 5][2]=-0.016293;
        mListOfCoordinates[ 6][0]=-0.266279; mListOfCoordinates[ 6][1]=-0.411144; mListOfCoordinates[ 6][2]=-0.013344;
        mListOfCoordinates[ 7][0]= 0.157703; mListOfCoordinates[ 7][1]=-0.123736; mListOfCoordinates[ 7][2]= 0.055287;
        mListOfCoordinates[ 8][0]= 0.088778; mListOfCoordinates[ 8][1]=-0.108746; mListOfCoordinates[ 8][2]=-0.108050;
        mListOfCoordinates[ 9][0]=-0.176388; mListOfCoordinates[ 9][1]=-0.355602; mListOfCoordinates[ 9][2]= 0.018067;
        mListOfCoordinates[10][0]= 0.401526; mListOfCoordinates[10][1]= 0.513011; mListOfCoordinates[10][2]=-0.070983;
        mListOfCoordinates[11][0]= 0.161663; mListOfCoordinates[11][1]=-0.292353; mListOfCoordinates[11][2]= 0.209007;
        mListOfCoordinates[12][0]=-0.312592; mListOfCoordinates[12][1]=-0.114346; mListOfCoordinates[12][2]=-0.127991;
        mListOfCoordinates[13][0]= 0.204706; mListOfCoordinates[13][1]= 0.398562; mListOfCoordinates[13][2]= 0.040206;
        mListOfCoordinates[14][0]=-0.016986; mListOfCoordinates[14][1]= 0.433223; mListOfCoordinates[14][2]=-0.218050;
        mListOfCoordinates[15][0]= 0.202056; mListOfCoordinates[15][1]= 0.462982; mListOfCoordinates[15][2]=-0.038994;
        mListOfCoordinates[16][0]=-0.156970; mListOfCoordinates[16][1]=-0.196157; mListOfCoordinates[16][2]=-0.014482;
        mListOfCoordinates[17][0]=-0.115325; mListOfCoordinates[17][1]= 0.010315; mListOfCoordinates[17][2]=-0.208098;
        mListOfCoordinates[18][0]= 0.021609; mListOfCoordinates[18][1]=-0.098776; mListOfCoordinates[18][2]= 0.106814;
        mListOfCoordinates[19][0]=-0.393505; mListOfCoordinates[19][1]=-0.344451; mListOfCoordinates[19][2]=-0.041832;
        mListOfCoordinates[20][0]=-0.217143; mListOfCoordinates[20][1]=-0.128369; mListOfCoordinates[20][2]=-0.114458;
        mListOfCoordinates[21][0]=-0.054566; mListOfCoordinates[21][1]=-0.362400; mListOfCoordinates[21][2]= 0.081134;
        mListOfCoordinates[22][0]=-0.054108; mListOfCoordinates[22][1]=-0.155989; mListOfCoordinates[22][2]=-0.042953;
        mListOfCoordinates[23][0]= 0.093046; mListOfCoordinates[23][1]= 0.229592; mListOfCoordinates[23][2]=-0.045959;
        mListOfCoordinates[24][0]=-0.008726; mListOfCoordinates[24][1]= 0.010088; mListOfCoordinates[24][2]=-0.153146;
        mListOfCoordinates[25][0]=-0.016238; mListOfCoordinates[25][1]= 0.285105; mListOfCoordinates[25][2]=-0.217501;
        mListOfCoordinates[26][0]= 0.412661; mListOfCoordinates[26][1]=-0.068805; mListOfCoordinates[26][2]= 0.092668;
        mListOfCoordinates[27][0]=-0.012104; mListOfCoordinates[27][1]= 0.119654; mListOfCoordinates[27][2]=-0.191545;
        mListOfCoordinates[28][0]=-0.490615; mListOfCoordinates[28][1]=-0.299183; mListOfCoordinates[28][2]=-0.104965;
        mListOfCoordinates[29][0]= 0.099773; mListOfCoordinates[29][1]=-0.049406; mListOfCoordinates[29][2]=-0.012620;
        mListOfCoordinates[30][0]=-0.004407; mListOfCoordinates[30][1]= 0.317069; mListOfCoordinates[30][2]= 0.209789;
        mListOfCoordinates[31][0]= 0.025055; mListOfCoordinates[31][1]=-0.338895; mListOfCoordinates[31][2]= 0.097088;
        mListOfCoordinates[32][0]=-0.023297; mListOfCoordinates[32][1]=-0.067427; mListOfCoordinates[32][2]=-0.117750;
        mListOfCoordinates[33][0]= 0.225288; mListOfCoordinates[33][1]= 0.075036; mListOfCoordinates[33][2]= 0.052082;
        mListOfCoordinates[34][0]=-0.001800; mListOfCoordinates[34][1]= 0.400404; mListOfCoordinates[34][2]= 0.212799;
        mListOfCoordinates[35][0]= 0.002649; mListOfCoordinates[35][1]= 0.217111; mListOfCoordinates[35][2]= 0.196087;
        mListOfCoordinates[36][0]=-0.243698; mListOfCoordinates[36][1]=-0.027877; mListOfCoordinates[36][2]= 0.145046;


        for (int i = 0; i < number_of_spheres; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        mListOfRadii[ 0]= 0.282575;
        mListOfRadii[ 1]= 0.317906;
        mListOfRadii[ 2]= 0.234842;
        mListOfRadii[ 3]= 0.368941;
        mListOfRadii[ 4]= 0.299773;
        mListOfRadii[ 5]= 0.299790;
        mListOfRadii[ 6]= 0.211190;
        mListOfRadii[ 7]= 0.358930;
        mListOfRadii[ 8]= 0.293120;
        mListOfRadii[ 9]= 0.255088;
        mListOfRadii[10]= 0.000000;
        mListOfRadii[11]= 0.275923;
        mListOfRadii[12]= 0.221766;
        mListOfRadii[13]= 0.251946;
        mListOfRadii[14]= 0.155318;
        mListOfRadii[15]= 0.192317;
        mListOfRadii[16]= 0.308834;
        mListOfRadii[17]= 0.213049;
        mListOfRadii[18]= 0.356396;
        mListOfRadii[19]= 0.239544;
        mListOfRadii[20]= 0.276506;
        mListOfRadii[21]= 0.254338;
        mListOfRadii[22]= 0.334751;
        mListOfRadii[23]= 0.286505;
        mListOfRadii[24]= 0.259457;
        mListOfRadii[25]= 0.184158;
        mListOfRadii[26]= 0.144672;
        mListOfRadii[27]= 0.203052;
        mListOfRadii[28]= 0.173081;
        mListOfRadii[29]= 0.360278;
        mListOfRadii[30]= 0.106554;
        mListOfRadii[31]= 0.254407;
        mListOfRadii[32]= 0.303123;
        mListOfRadii[33]= 0.284376;
        mListOfRadii[34]= 0.107573;
        mListOfRadii[35]= 0.120925;
        mListOfRadii[36]= 0.134072;
        
        for (int i = 0; i < number_of_spheres; i++) { mListOfRadii[i]= mListOfRadii[i] * cl; }
       
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.494047 * cl * cl * cl;
                
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = cluster_mass * cl * cl * 0.150066;
        base_principal_moments_of_inertia[1] = cluster_mass * cl * cl * 0.217885;
        base_principal_moments_of_inertia[2] = cluster_mass * cl * cl * 0.267374;

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

