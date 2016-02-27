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
#include "rock1cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Rock1Cluster3D::Rock1Cluster3D() : Cluster3D() {}
            
      
    Rock1Cluster3D::Rock1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Rock1Cluster3D::Rock1Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Rock1Cluster3D::Rock1Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Rock1Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Rock1Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Rock1Cluster3D::~Rock1Cluster3D() {}
      
    
    void Rock1Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 85;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.1455 (in meters) was the medium diameter of the rock in GiD (rock05.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.87298,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                
        cl *= 0.87298;
                
        double a, b, c;
        
        a = 0.572749; // this is the original semi-axis (1.1455/2)
        b = a; c = b;
        
        mListOfCoordinates[ 0][0] =-0.4059089300155639; mListOfCoordinates[ 0][1] =-0.111142139244079; mListOfCoordinates[ 0][2] = 0.17642473061084818;
        mListOfCoordinates[ 1][0] = 0.4734444112777709; mListOfCoordinates[ 1][1] = 0.073340840530395; mListOfCoordinates[ 1][2] =-0.02537642452716770;
        mListOfCoordinates[ 2][0] = 0.4013676259994507; mListOfCoordinates[ 2][1] = 0.061800733184814; mListOfCoordinates[ 2][2] =-0.19055143506526923;
        mListOfCoordinates[ 3][0] = 0.4351699476242065; mListOfCoordinates[ 3][1] = 0.232795501327514; mListOfCoordinates[ 3][2] =-0.05244942595958779;
        mListOfCoordinates[ 4][0] = 0.5036787439346313; mListOfCoordinates[ 4][1] = 0.221658722305297; mListOfCoordinates[ 4][2] = 0.07132912433147442;
        mListOfCoordinates[ 5][0] = 0.5443247541427613; mListOfCoordinates[ 5][1] = 0.043943148040771; mListOfCoordinates[ 5][2] = 0.14715090754032101;
        mListOfCoordinates[ 6][0] = 0.4602648550033572; mListOfCoordinates[ 6][1] =-0.082053426933288; mListOfCoordinates[ 6][2] = 0.01905415561199230;
        mListOfCoordinates[ 7][0] = 0.3347237546920778; mListOfCoordinates[ 7][1] = 0.182205056762695; mListOfCoordinates[ 7][2] =-0.23256184537410624;
        mListOfCoordinates[ 8][0] = 0.5646031087875365; mListOfCoordinates[ 8][1] =-0.096323633766174; mListOfCoordinates[ 8][2] = 0.16870258734226251;
        mListOfCoordinates[ 9][0] = 0.3287808710098267; mListOfCoordinates[ 9][1] =-0.181287511062622; mListOfCoordinates[ 9][2] =-0.13149494359493286;
        mListOfCoordinates[10][0] = 0.2688478113174437; mListOfCoordinates[10][1] = 0.046211501312255; mListOfCoordinates[10][2] =-0.33859746625423442;
        mListOfCoordinates[11][0] = 0.3480440586090087; mListOfCoordinates[11][1] = 0.354911084365844; mListOfCoordinates[11][2] =-0.17006706864833898;
        mListOfCoordinates[12][0] = 0.4267853559494019; mListOfCoordinates[12][1] = 0.376573128509521; mListOfCoordinates[12][2] = 0.06393966071605650;
        mListOfCoordinates[13][0] = 0.4360494642257688; mListOfCoordinates[13][1] = 0.292704973983764; mListOfCoordinates[13][2] = 0.17923974592685621;
        mListOfCoordinates[14][0] = 0.4824625547409053; mListOfCoordinates[14][1] = 0.049741461181640; mListOfCoordinates[14][2] = 0.29653935091495587;
        mListOfCoordinates[15][0] = 0.3834672643661500; mListOfCoordinates[15][1] =-0.249273169898986; mListOfCoordinates[15][2] = 0.12303228857517201;
        mListOfCoordinates[16][0] = 0.2679400602340698; mListOfCoordinates[16][1] =-0.130030882644653; mListOfCoordinates[16][2] =-0.28996946809291740;
        mListOfCoordinates[17][0] = 0.1949276391983032; mListOfCoordinates[17][1] = 0.211629934310913; mListOfCoordinates[17][2] =-0.32624310643672866;
        mListOfCoordinates[18][0] = 0.3209657163619994; mListOfCoordinates[18][1] = 0.472557723617553; mListOfCoordinates[18][2] =-0.06411295297145880;
        mListOfCoordinates[19][0] = 0.3708844800949096; mListOfCoordinates[19][1] = 0.167581144332885; mListOfCoordinates[19][2] = 0.32892922384738998;
        mListOfCoordinates[20][0] = 0.4377691190719603; mListOfCoordinates[20][1] =-0.140802462387085; mListOfCoordinates[20][2] = 0.35195349428653788;
        mListOfCoordinates[21][0] = 0.4323379541397096; mListOfCoordinates[21][1] =-0.257158136177063; mListOfCoordinates[21][2] = 0.25625881483554824;
        mListOfCoordinates[22][0] = 0.2482549795150756; mListOfCoordinates[22][1] =-0.325963994216919; mListOfCoordinates[22][2] = 0.02574222035408002;
        mListOfCoordinates[23][0] = 0.1967921144485471; mListOfCoordinates[23][1] =-0.281363657379150; mListOfCoordinates[23][2] =-0.11305880143642444;
        mListOfCoordinates[24][0] = 0.0651245376110074; mListOfCoordinates[24][1] = 0.158836550521850; mListOfCoordinates[24][2] =-0.38417096822261870;
        mListOfCoordinates[25][0] = 0.1459672737121582; mListOfCoordinates[25][1] = 0.447385758590698; mListOfCoordinates[25][2] =-0.27225562894344341;
        mListOfCoordinates[26][0] = 0.2992995994567870; mListOfCoordinates[26][1] = 0.446245985031127; mListOfCoordinates[26][2] = 0.09576077749729083;
        mListOfCoordinates[27][0] = 0.3447016206741333; mListOfCoordinates[27][1] = 0.001292739486694; mListOfCoordinates[27][2] = 0.41058149626255069;
        mListOfCoordinates[28][0] = 0.3102501558303834; mListOfCoordinates[28][1] =-0.380650337600708; mListOfCoordinates[28][2] = 0.28436096937656341;
        mListOfCoordinates[29][0] = 0.1993804918289183; mListOfCoordinates[29][1] =-0.392162785148620; mListOfCoordinates[29][2] = 0.16701459009647390;
        mListOfCoordinates[30][0] = 0.1289813469886780; mListOfCoordinates[30][1] =-0.162152777862549; mListOfCoordinates[30][2] =-0.27579260995387939;
        mListOfCoordinates[31][0] = 0.0968052815437314; mListOfCoordinates[31][1] =-0.004518916320800; mListOfCoordinates[31][2] =-0.35759917256832191;
        mListOfCoordinates[32][0] = 0.1291622600555417; mListOfCoordinates[32][1] = 0.507337607955932; mListOfCoordinates[32][2] =-0.12748357331752747;
        mListOfCoordinates[33][0] = 0.2226009794235228; mListOfCoordinates[33][1] = 0.360639935684204; mListOfCoordinates[33][2] = 0.27350471117496422;
        mListOfCoordinates[34][0] = 0.2215446350097656; mListOfCoordinates[34][1] = 0.167033445358276; mListOfCoordinates[34][2] = 0.40253746740818119;
        mListOfCoordinates[35][0] = 0.3044548410415648; mListOfCoordinates[35][1] =-0.198477912902832; mListOfCoordinates[35][2] = 0.43408275454044359;
        mListOfCoordinates[36][0] = 0.0716768550872801; mListOfCoordinates[36][1] =-0.451560339164733; mListOfCoordinates[36][2] = 0.09333994553089208;
        mListOfCoordinates[37][0] = 0.0512112783432006; mListOfCoordinates[37][1] =-0.375374312400817; mListOfCoordinates[37][2] =-0.07637091577053123;
        mListOfCoordinates[38][0] =-0.1096965712547303; mListOfCoordinates[38][1] = 0.148397088241577; mListOfCoordinates[38][2] =-0.38420449979305277;
        mListOfCoordinates[39][0] =-0.0675114569664000; mListOfCoordinates[39][1] = 0.302265419387817; mListOfCoordinates[39][2] =-0.35367502934932654;
        mListOfCoordinates[40][0] =-0.0331960923194886; mListOfCoordinates[40][1] = 0.424204599380493; mListOfCoordinates[40][2] =-0.25188080232143312;
        mListOfCoordinates[41][0] = 0.1031348651885986; mListOfCoordinates[41][1] = 0.502572926712036; mListOfCoordinates[41][2] = 0.04219804468154908;
        mListOfCoordinates[42][0] = 0.1930046960830689; mListOfCoordinates[42][1] =-0.031802120971679; mListOfCoordinates[42][2] = 0.45881696016788448;
        mListOfCoordinates[43][0] = 0.1777827748298646; mListOfCoordinates[43][1] =-0.357880552291870; mListOfCoordinates[43][2] = 0.38924115030765538;
        mListOfCoordinates[44][0] = 0.0916415321350096; mListOfCoordinates[44][1] =-0.405296241760254; mListOfCoordinates[44][2] = 0.26290571157932252;
        mListOfCoordinates[45][0] =-0.0855454576492310; mListOfCoordinates[45][1] =-0.385461668205261; mListOfCoordinates[45][2] =-0.13886199452876991;
        mListOfCoordinates[46][0] =-0.0179865767478942; mListOfCoordinates[46][1] =-0.275569766807556; mListOfCoordinates[46][2] =-0.22950497548580179;
        mListOfCoordinates[47][0] =-0.0791099072456359; mListOfCoordinates[47][1] =-0.153505855369568; mListOfCoordinates[47][2] =-0.32649596745967907;
        mListOfCoordinates[48][0] =-0.0519204635620117; mListOfCoordinates[48][1] = 0.467497504806518; mListOfCoordinates[48][2] =-0.10875714070796949;
        mListOfCoordinates[49][0] = 0.0364614850997924; mListOfCoordinates[49][1] = 0.418675238418579; mListOfCoordinates[49][2] = 0.13530640499591851;
        mListOfCoordinates[50][0] = 0.0416160520553587; mListOfCoordinates[50][1] = 0.336533784866333; mListOfCoordinates[50][2] = 0.29296086771488117;
        mListOfCoordinates[51][0] = 0.0415515874862668; mListOfCoordinates[51][1] = 0.126838487243652; mListOfCoordinates[51][2] = 0.42631612322330442;
        mListOfCoordinates[52][0] = 0.1214497115135191; mListOfCoordinates[52][1] =-0.218759203433990; mListOfCoordinates[52][2] = 0.44422787859439788;
        mListOfCoordinates[53][0] =-0.0743732855796813; mListOfCoordinates[53][1] =-0.447769972419738; mListOfCoordinates[53][2] = 0.21266080210208838;
        mListOfCoordinates[54][0] =-0.1236146716117859; mListOfCoordinates[54][1] =-0.456402153015136; mListOfCoordinates[54][2] = 0.05444371111393042;
        mListOfCoordinates[55][0] =-0.1744168872833252; mListOfCoordinates[55][1] =-0.052497046661377; mListOfCoordinates[55][2] =-0.35289335896968860;
        mListOfCoordinates[56][0] =-0.2886909051895142; mListOfCoordinates[56][1] = 0.130849282836913; mListOfCoordinates[56][2] =-0.34010083997249568;
        mListOfCoordinates[57][0] =-0.2329027673721314; mListOfCoordinates[57][1] = 0.300810918807983; mListOfCoordinates[57][2] =-0.24706867215633371;
        mListOfCoordinates[58][0] =-0.0879464755058289; mListOfCoordinates[58][1] = 0.399750239944458; mListOfCoordinates[58][2] = 0.04652738082408891;
        mListOfCoordinates[59][0] =-0.0455796985626221; mListOfCoordinates[59][1] = 0.236131317901611; mListOfCoordinates[59][2] = 0.37295285265445788;
        mListOfCoordinates[60][0] =-0.0100626006603240; mListOfCoordinates[60][1] =-0.089224510574341; mListOfCoordinates[60][2] = 0.44915112268924690;
        mListOfCoordinates[61][0] = 0.0004674743175506; mListOfCoordinates[61][1] =-0.325020887565612; mListOfCoordinates[61][2] = 0.36416195929050499;
        mListOfCoordinates[62][0] =-0.2723078096389771; mListOfCoordinates[62][1] =-0.409598627090454; mListOfCoordinates[62][2] =-0.09644285886287651;
        mListOfCoordinates[63][0] =-0.2275727468490601; mListOfCoordinates[63][1] =-0.248208055305481; mListOfCoordinates[63][2] =-0.24545365979671535;
        mListOfCoordinates[64][0] =-0.2126167919158935; mListOfCoordinates[64][1] = 0.277475944137573; mListOfCoordinates[64][2] =-0.10951729791164322;
        mListOfCoordinates[65][0] =-0.1259296492576599; mListOfCoordinates[65][1] = 0.303611579513549; mListOfCoordinates[65][2] = 0.16006849529743139;
        mListOfCoordinates[66][0] =-0.1229696261405944; mListOfCoordinates[66][1] = 0.068101969909667; mListOfCoordinates[66][2] = 0.42015373518467047;
        mListOfCoordinates[67][0] =-0.2398387802124023; mListOfCoordinates[67][1] =-0.415882673263549; mListOfCoordinates[67][2] = 0.18596135160923044;
        mListOfCoordinates[68][0] =-0.3492557676315307; mListOfCoordinates[68][1] =-0.121619580650329; mListOfCoordinates[68][2] =-0.31143167207240990;
        mListOfCoordinates[69][0] =-0.3845537672042848; mListOfCoordinates[69][1] = 0.011275652694701; mListOfCoordinates[69][2] =-0.25518171308040655;
        mListOfCoordinates[70][0] =-0.3695017900466920; mListOfCoordinates[70][1] = 0.185836410903930; mListOfCoordinates[70][2] =-0.10226518380641902;
        mListOfCoordinates[71][0] =-0.2650948740005493; mListOfCoordinates[71][1] = 0.233218668746948; mListOfCoordinates[71][2] = 0.10217342174053190;
        mListOfCoordinates[72][0] =-0.2116164916992188; mListOfCoordinates[72][1] = 0.256066311264038; mListOfCoordinates[72][2] = 0.28925195620060062;
        mListOfCoordinates[73][0] =-0.1532823032379150; mListOfCoordinates[73][1] =-0.134830444717407; mListOfCoordinates[73][2] = 0.34775432722568578;
        mListOfCoordinates[74][0] =-0.3378754974365235; mListOfCoordinates[74][1] =-0.420167843055725; mListOfCoordinates[74][2] = 0.05040523698329924;
        mListOfCoordinates[75][0] =-0.4267853559494017; mListOfCoordinates[75][1] =-0.376573128509521; mListOfCoordinates[75][2] =-0.06393966071605650;
        mListOfCoordinates[76][0] =-0.3716675497055054; mListOfCoordinates[76][1] =-0.272471319770813; mListOfCoordinates[76][2] =-0.18140236356258382;
        mListOfCoordinates[77][0] =-0.4312976259231565; mListOfCoordinates[77][1] = 0.039080370712280; mListOfCoordinates[77][2] =-0.12559765985012117;
        mListOfCoordinates[78][0] =-0.2487496326446533; mListOfCoordinates[78][1] = 0.055495398712157; mListOfCoordinates[78][2] = 0.32740528643131217;
        mListOfCoordinates[79][0] =-0.2870210729598999; mListOfCoordinates[79][1] =-0.200657190132141; mListOfCoordinates[79][2] = 0.21992153646945933;
        mListOfCoordinates[80][0] =-0.4622705526351929; mListOfCoordinates[80][1] =-0.153876851463318; mListOfCoordinates[80][2] =-0.17122945477962492;
        mListOfCoordinates[81][0] =-0.4278863508224484; mListOfCoordinates[81][1] = 0.060835263824462; mListOfCoordinates[81][2] = 0.01314950716495477;
        mListOfCoordinates[82][0] =-0.3568963682174680; mListOfCoordinates[82][1] = 0.101075247192382; mListOfCoordinates[82][2] = 0.15478396227359870;
        mListOfCoordinates[83][0] =-0.3929164709091187; mListOfCoordinates[83][1] =-0.240009811210632; mListOfCoordinates[83][2] = 0.07486742961406722;
        mListOfCoordinates[84][0] =-0.4530466222763061; mListOfCoordinates[84][1] =-0.092443091583252; mListOfCoordinates[84][2] = 0.00000000000000000;
             
        for (int i = 0; i < 85; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        for (int i = 0; i < 85; i++) { mListOfRadii[i]= 0.15 * cl; }
                        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
        //double cluster_volume = 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c * cl * cl * cl; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        //double cluster_mass = particle_density * cluster_volume;
        
        ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c * cl * cl * cl;
        
        //GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = cluster_mass;
        //GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = 0.2 * cluster_mass * (b * b * cl * cl + c * c * cl * cl);
        //GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = 0.2 * cluster_mass * (a * a * cl * cl + c * c * cl * cl);
        //GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = 0.2 * cluster_mass * (a * a * cl * cl + b * b * cl * cl);
        //array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.2 * cluster_mass * cl * cl * (b * b + c * c);
        base_principal_moments_of_inertia[1] = 0.2 * cluster_mass * cl * cl * (a * a + c * c);
        base_principal_moments_of_inertia[2] = 0.2 * cluster_mass * cl * cl * (a * a + b * b); 
  
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Rock1Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock1Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Rock1Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Rock1Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

