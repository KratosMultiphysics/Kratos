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
#include "rock2cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "custom_utilities/properties_proxies.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Rock2Cluster3D::Rock2Cluster3D() : Cluster3D() {}
            
      
    Rock2Cluster3D::Rock2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Cluster3D(NewId, pGeometry) {}
      
      
    Rock2Cluster3D::Rock2Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Cluster3D(NewId, pGeometry, pProperties) {}

      
    Rock2Cluster3D::Rock2Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Cluster3D(NewId, ThisNodes) {}

      
    Element::Pointer Rock2Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Cluster3D::Pointer(new Rock2Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Rock2Cluster3D::~Rock2Cluster3D() {}
      
    
    void Rock2Cluster3D::CustomInitialize() {
        
        int number_of_spheres = 76;
        mListOfRadii.resize(number_of_spheres);
        mListOfCoordinates.resize(number_of_spheres);
        mListOfSphericParticles.resize(number_of_spheres);
        
        double cl = GetGeometry()[0].FastGetSolutionStepValue(CHARACTERISTIC_LENGTH);
        
        // 1.06937 (in meters) was the medium diameter of the rock in GiD (rock08.gid file)
        // so we should multiply every size that follow by the inverse of that number, 0.93513,
        // to obtain a 'unity' rock.
        // We then have to multiply again everything by 'cl' to obtain the desired dimensions
        // to adjust to the characteristic length given
                       
        cl *= 0.93513;
                
        double a, b, c;
        
        a = 0.534685 * cl; // this is the original semi-axis (1.06937/2)
        b = a; c = b;
        
        mListOfCoordinates[ 0][0] =-0.16199315983425355; mListOfCoordinates[ 0][1] =-0.42102103032777205; mListOfCoordinates[ 0][2] =-0.02792826491079630;
        mListOfCoordinates[ 1][0] = 0.23403579895576743; mListOfCoordinates[ 1][1] = 0.40993540195966893; mListOfCoordinates[ 1][2] =-0.05424263978611415;
        mListOfCoordinates[ 2][0] = 0.08185167087129807; mListOfCoordinates[ 2][1] = 0.44659347329098542; mListOfCoordinates[ 2][2] = 0.03579734091859615;
        mListOfCoordinates[ 3][0] = 0.28672188104535412; mListOfCoordinates[ 3][1] = 0.39125427548705594; mListOfCoordinates[ 3][2] = 0.10227848687803827;
        mListOfCoordinates[ 4][0] = 0.38029896030677018; mListOfCoordinates[ 4][1] = 0.35242748785860212; mListOfCoordinates[ 4][2] =-0.10140436471834796;
        mListOfCoordinates[ 5][0] = 0.27057548158552291; mListOfCoordinates[ 5][1] = 0.37384491563152789; mListOfCoordinates[ 5][2] =-0.23229535829529863;
        mListOfCoordinates[ 6][0] = 0.06026267949569630; mListOfCoordinates[ 6][1] = 0.39334375378995590; mListOfCoordinates[ 6][2] =-0.16336963750757943;
        mListOfCoordinates[ 7][0] = 0.17558052265734558; mListOfCoordinates[ 7][1] = 0.40740073984879610; mListOfCoordinates[ 7][2] = 0.20757942452270034;
        mListOfCoordinates[ 8][0] = 0.41820590832626819; mListOfCoordinates[ 8][1] = 0.31786192721489381; mListOfCoordinates[ 8][2] = 0.08125249053753252;
        mListOfCoordinates[ 9][0] = 0.11000435784720264; mListOfCoordinates[ 9][1] = 0.31302484446149681; mListOfCoordinates[ 9][2] =-0.28394352800670380;
        mListOfCoordinates[10][0] =-0.06389455326369086; mListOfCoordinates[10][1] = 0.39414058170737920; mListOfCoordinates[10][2] =-0.06975019684356226;
        mListOfCoordinates[11][0] =-0.09052976924670896; mListOfCoordinates[11][1] = 0.42409058009134892; mListOfCoordinates[11][2] = 0.13326739624340606;
        mListOfCoordinates[12][0] = 0.34739715083749628; mListOfCoordinates[12][1] = 0.36161555796118672; mListOfCoordinates[12][2] = 0.24317355824851400;
        mListOfCoordinates[13][0] = 0.47271469954349171; mListOfCoordinates[13][1] = 0.21267449463211122; mListOfCoordinates[13][2] =-0.08450416932362686;
        mListOfCoordinates[14][0] = 0.36482420140036492; mListOfCoordinates[14][1] = 0.26314507618418814; mListOfCoordinates[14][2] =-0.32807736813744742;
        mListOfCoordinates[15][0] =-0.07171471078474145; mListOfCoordinates[15][1] = 0.30229088791435077; mListOfCoordinates[15][2] =-0.28859506039726096;
        mListOfCoordinates[16][0] =-0.14853938092710262; mListOfCoordinates[16][1] = 0.29050450632993585; mListOfCoordinates[16][2] =-0.16019949679078488;
        mListOfCoordinates[17][0] =-0.02542489557970373; mListOfCoordinates[17][1] = 0.35344154022360297; mListOfCoordinates[17][2] = 0.28737822850033728;
        mListOfCoordinates[18][0] = 0.49661412053635773; mListOfCoordinates[18][1] = 0.17309764645871875; mListOfCoordinates[18][2] = 0.09643547650315283;
        mListOfCoordinates[19][0] = 0.42776582826526705; mListOfCoordinates[19][1] = 0.13353934972410145; mListOfCoordinates[19][2] =-0.26846148257444524;
        mListOfCoordinates[20][0] = 0.21732845629814312; mListOfCoordinates[20][1] = 0.22020054201769748; mListOfCoordinates[20][2] =-0.44661907906112447;
        mListOfCoordinates[21][0] = 0.00227466143361387; mListOfCoordinates[21][1] = 0.17588435314899811; mListOfCoordinates[21][2] =-0.37746248358808404;
        mListOfCoordinates[22][0] =-0.22069257188043767; mListOfCoordinates[22][1] = 0.35598023039737847; mListOfCoordinates[22][2] = 0.02270157338409156;
        mListOfCoordinates[23][0] =-0.15698261097845456; mListOfCoordinates[23][1] = 0.28187637880590422; mListOfCoordinates[23][2] = 0.31109781776714945;
        mListOfCoordinates[24][0] = 0.11916423535235821; mListOfCoordinates[24][1] = 0.27781106957775970; mListOfCoordinates[24][2] = 0.47249059686794553;
        mListOfCoordinates[25][0] = 0.27044715234663019; mListOfCoordinates[25][1] = 0.25393401859219200; mListOfCoordinates[25][2] = 0.42999059184925031;
        mListOfCoordinates[26][0] = 0.38496870816153433; mListOfCoordinates[26][1] = 0.21019157583476825; mListOfCoordinates[26][2] = 0.31420864383226349;
        mListOfCoordinates[27][0] = 0.45771033038392434; mListOfCoordinates[27][1] = 0.03829059952304658; mListOfCoordinates[27][2] =-0.10655135601305871;
        mListOfCoordinates[28][0] = 0.33298767612776958; mListOfCoordinates[28][1] = 0.09090230076745489; mListOfCoordinates[28][2] =-0.41776438722166392;
        mListOfCoordinates[29][0] = 0.10749968063714899; mListOfCoordinates[29][1] = 0.10665328629063189; mListOfCoordinates[29][2] =-0.52401511353307484;
        mListOfCoordinates[30][0] =-0.16140317288389727; mListOfCoordinates[30][1] = 0.14140137194018157; mListOfCoordinates[30][2] =-0.35532097617180575;
        mListOfCoordinates[31][0] =-0.23745365826193410; mListOfCoordinates[31][1] = 0.16083813066793887; mListOfCoordinates[31][2] =-0.18851285381355382;
        mListOfCoordinates[32][0] =-0.07305987061386065; mListOfCoordinates[32][1] =-0.00065661797248794; mListOfCoordinates[32][2] =-0.42230747967272114;
        mListOfCoordinates[33][0] =-0.29032909724408495; mListOfCoordinates[33][1] = 0.22314227902176711; mListOfCoordinates[33][2] =-0.01838220933376660;
        mListOfCoordinates[34][0] =-0.24789769041569165; mListOfCoordinates[34][1] = 0.13934085881543193; mListOfCoordinates[34][2] = 0.35680862036434996;
        mListOfCoordinates[35][0] =-0.12061955198471008; mListOfCoordinates[35][1] = 0.15323531509912738; mListOfCoordinates[35][2] = 0.45678019294849675;
        mListOfCoordinates[36][0] = 0.03266165116733488; mListOfCoordinates[36][1] = 0.12019910722150953; mListOfCoordinates[36][2] = 0.49946806304845853;
        mListOfCoordinates[37][0] = 0.17909680947438067; mListOfCoordinates[37][1] = 0.10073283454276896; mListOfCoordinates[37][2] = 0.45937156731120360;
        mListOfCoordinates[38][0] = 0.33993835789668292; mListOfCoordinates[38][1] = 0.06036472277878357; mListOfCoordinates[38][2] = 0.38452910142819013;
        mListOfCoordinates[39][0] = 0.43761253940163242; mListOfCoordinates[39][1] = 0.03912205149072844; mListOfCoordinates[39][2] = 0.20913620028902430;
        mListOfCoordinates[40][0] = 0.45231596514306288; mListOfCoordinates[40][1] =-0.08315249851788833; mListOfCoordinates[40][2] = 0.02825358288538609;
        mListOfCoordinates[41][0] = 0.39446189312953678; mListOfCoordinates[41][1] =-0.13566781051313517; mListOfCoordinates[41][2] =-0.18171126897548090;
        mListOfCoordinates[42][0] = 0.32248483922429350; mListOfCoordinates[42][1] =-0.09805364914821548; mListOfCoordinates[42][2] =-0.36745594563327894;
        mListOfCoordinates[43][0] = 0.12609595851529198; mListOfCoordinates[43][1] =-0.06245909414775910; mListOfCoordinates[43][2] =-0.51301039927846626;
        mListOfCoordinates[44][0] =-0.21412168653071084; mListOfCoordinates[44][1] =-0.06855460016488648; mListOfCoordinates[44][2] =-0.41992516273846886;
        mListOfCoordinates[45][0] =-0.30793521679445501; mListOfCoordinates[45][1] = 0.01362528178963540; mListOfCoordinates[45][2] =-0.28162027298296433;
        mListOfCoordinates[46][0] =-0.34570147574808485; mListOfCoordinates[46][1] = 0.16357966544130811; mListOfCoordinates[46][2] = 0.18177106915333380;
        mListOfCoordinates[47][0] =-0.01134334795676805; mListOfCoordinates[47][1] =-0.14309462094312733; mListOfCoordinates[47][2] =-0.50078455793871712;
        mListOfCoordinates[48][0] =-0.39497424163743589; mListOfCoordinates[48][1] =-0.06133509162921349; mListOfCoordinates[48][2] =-0.18139860004143155;
        mListOfCoordinates[49][0] =-0.38473417538332388; mListOfCoordinates[49][1] = 0.06463926253976748; mListOfCoordinates[49][2] =-0.00073471688929366;
        mListOfCoordinates[50][0] =-0.31526220786715264; mListOfCoordinates[50][1] =-0.03193045762958082; mListOfCoordinates[50][2] = 0.37730710363917802;
        mListOfCoordinates[51][0] =-0.13046953957806384; mListOfCoordinates[51][1] =-0.04329495897309523; mListOfCoordinates[51][2] = 0.45716182251965615;
        mListOfCoordinates[52][0] = 0.12699029822006469; mListOfCoordinates[52][1] =-0.08168091814715886; mListOfCoordinates[52][2] = 0.40577376882046873;
        mListOfCoordinates[53][0] = 0.25149714925894451; mListOfCoordinates[53][1] =-0.11721157807385976; mListOfCoordinates[53][2] = 0.29514363827599138;
        mListOfCoordinates[54][0] = 0.34307932557566617; mListOfCoordinates[54][1] =-0.22208373635346193; mListOfCoordinates[54][2] = 0.16072496717091384;
        mListOfCoordinates[55][0] = 0.32522680582358710; mListOfCoordinates[55][1] =-0.28524461075232421; mListOfCoordinates[55][2] =-0.04165798897635891;
        mListOfCoordinates[56][0] = 0.22554984505659637; mListOfCoordinates[56][1] =-0.25554477852684443; mListOfCoordinates[56][2] =-0.26026942204122994;
        mListOfCoordinates[57][0] = 0.12964384468042478; mListOfCoordinates[57][1] =-0.21603413770964500; mListOfCoordinates[57][2] =-0.39526767866933771;
        mListOfCoordinates[58][0] =-0.20759588315941702; mListOfCoordinates[58][1] =-0.22906630763253735; mListOfCoordinates[58][2] =-0.40448294489852660;
        mListOfCoordinates[59][0] =-0.40599489658690391; mListOfCoordinates[59][1] =-0.00942981571463086; mListOfCoordinates[59][2] = 0.20627890451390091;
        mListOfCoordinates[60][0] = 0.04878234143394022; mListOfCoordinates[60][1] =-0.15778808174483375; mListOfCoordinates[60][2] = 0.29525069398788012;
        mListOfCoordinates[61][0] =-0.46428793501847443; mListOfCoordinates[61][1] =-0.21822762278094002; mListOfCoordinates[61][2] =-0.26992903911253985;
        mListOfCoordinates[62][0] =-0.49661412053635778; mListOfCoordinates[62][1] =-0.17309764645871875; mListOfCoordinates[62][2] =-0.09643547650315281;
        mListOfCoordinates[63][0] =-0.47231687919083731; mListOfCoordinates[63][1] =-0.09855357809311918; mListOfCoordinates[63][2] = 0.08724943482085417;
        mListOfCoordinates[64][0] =-0.37207470234201112; mListOfCoordinates[64][1] =-0.15691712039136543; mListOfCoordinates[64][2] = 0.26845181412883673;
        mListOfCoordinates[65][0] =-0.16440058380854305; mListOfCoordinates[65][1] =-0.23618234946990771; mListOfCoordinates[65][2] = 0.33678898701891941;
        mListOfCoordinates[66][0] = 0.12794893274616437; mListOfCoordinates[66][1] =-0.27179489029087478; mListOfCoordinates[66][2] = 0.22776769956847437;
        mListOfCoordinates[67][0] = 0.20174042423804217; mListOfCoordinates[67][1] =-0.38117529992105548; mListOfCoordinates[67][2] = 0.00989673768889937;
        mListOfCoordinates[68][0] = 0.06237751245203371; mListOfCoordinates[68][1] =-0.31312457737456179; mListOfCoordinates[68][2] =-0.20611731032989428;
        mListOfCoordinates[69][0] =-0.16128821575663410; mListOfCoordinates[69][1] =-0.33007404348023262; mListOfCoordinates[69][2] =-0.23588550296606992;
        mListOfCoordinates[70][0] =-0.31113773809602596; mListOfCoordinates[70][1] =-0.30998310840619292; mListOfCoordinates[70][2] =-0.24943126752758416;
        mListOfCoordinates[71][0] =-0.44627395477283305; mListOfCoordinates[71][1] =-0.26565760311971431; mListOfCoordinates[71][2] = 0.11668485990001749;
        mListOfCoordinates[72][0] =-0.28276970983594640; mListOfCoordinates[72][1] =-0.32476364028211357; mListOfCoordinates[72][2] = 0.15570652222850240;
        mListOfCoordinates[73][0] =-0.10487601160465088; mListOfCoordinates[73][1] =-0.38078919136031775; mListOfCoordinates[73][2] = 0.18070971348415327;
        mListOfCoordinates[74][0] =-0.00761306269617045; mListOfCoordinates[74][1] =-0.39913408641999293; mListOfCoordinates[74][2] = 0.06737254532340403;
        mListOfCoordinates[75][0] =-0.29599030492129608; mListOfCoordinates[75][1] =-0.34935882464685342; mListOfCoordinates[75][2] =-0.09116480215731476;
 
        for (int i = 0; i < 76; i++) { 
            mListOfCoordinates[i][0] *= cl;  mListOfCoordinates[i][1] *= cl;  mListOfCoordinates[i][2] *= cl;
        }
        
        for (int i = 0; i < 76; i++) { mListOfRadii[i]= 0.15 * cl; }
                        
        //double particle_density = this->SlowGetDensity(); /////////////////////////////USE FAST
         
        //double cluster_volume = 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c; ////APPROXIMATE VOLUME, CALCULATE MORE EXACTLY
        
        //double cluster_mass = particle_density * cluster_volume;
        
        double cluster_mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = (this->SlowGetDensity()) * 0.3333333333 * 4.0 * KRATOS_M_PI * a * b * c;
        
        array_1d<double, 3>& base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        base_principal_moments_of_inertia[0] = 0.2 * cluster_mass * (b * b + c * c);
        base_principal_moments_of_inertia[1] = 0.2 * cluster_mass * (a * a + c * c);
        base_principal_moments_of_inertia[2] = 0.2 * cluster_mass * (a * a + b * b);
          
    }     
    
      
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& r_process_info) {}
  
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::InitializeSolutionStep(ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        KRATOS_CATCH("")
        
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::FinalizeSolutionStep(ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************
    
    void Rock2Cluster3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
          
        KRATOS_TRY

        KRATOS_CATCH("")

    }// Calculate

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& r_process_info) {}

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void Rock2Cluster3D::Calculate(const Variable<Vector>& rVariable, Vector& Output, const ProcessInfo& r_process_info){}
    void Rock2Cluster3D::Calculate(const Variable<Matrix>& rVariable, Matrix& Output, const ProcessInfo& r_process_info){}
    double Rock2Cluster3D::SlowGetDensity()                                        { return GetProperties()[PARTICLE_DENSITY];}

}  // namespace Kratos.

