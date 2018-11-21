#include "catch.hpp"

#include <ANurbs/src/Curve.h>
#include <ANurbs/src/CurveOnSurface.h>
#include <ANurbs/src/PointOnCurveProjection.h>

TEST_CASE( "Project point on spatial curve", "[PointOnCurveProjection]" ) {
    using namespace ANurbs;

    int degree = 4;
    int nbPoles = 8;
    bool isRational = false;

    Pointer<CurveGeometry3D> curveGeometry = New<CurveGeometry3D>(degree,
        nbPoles, isRational);
    {
        curveGeometry->SetKnot( 0,  3.0 );
        curveGeometry->SetKnot( 1,  3.0 );
        curveGeometry->SetKnot( 2,  3.0 );
        curveGeometry->SetKnot( 3,  3.0 );
        curveGeometry->SetKnot( 4,  6.5 );
        curveGeometry->SetKnot( 5, 10.0 );
        curveGeometry->SetKnot( 6, 13.5);
        curveGeometry->SetKnot( 7, 17.0 );
        curveGeometry->SetKnot( 8, 17.0 );
        curveGeometry->SetKnot( 9, 17.0 );
        curveGeometry->SetKnot(10, 17.0 );

        curveGeometry->SetPole(0, {  0, -25, - 5});
        curveGeometry->SetPole(1, {-15, -15,   0});
        curveGeometry->SetPole(2, {  5, - 5, - 3});
        curveGeometry->SetPole(3, { 15, -15,   3});
        curveGeometry->SetPole(4, { 25,   0,   6});
        curveGeometry->SetPole(5, { 15,  15,   6});
        curveGeometry->SetPole(6, {- 5, - 5, - 3});
        curveGeometry->SetPole(7, {-25,  15,   4});
    }

    Pointer<Curve3D> curve = New<Curve3D>(curveGeometry, curveGeometry->Domain());

    PointOnCurveProjection3D projection(curve, 1e-7);

    SECTION("Point 1")
    {
        Point3D point = {-25.1331415843, -38.9256022819, -3.2989320128};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.3983282912) );
    }

    SECTION("Point 2")
    {
        Point3D point = {35.6464813397, 27.3703996918, -41.1153099924};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.3339477287) );
    }

    SECTION("Point 3")
    {
        Point3D point = {-40.3995502695, 45.1689836547, -1.7412051334};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 4")
    {
        Point3D point = {36.117074539, 44.5183648237, 47.2049699152};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.0827431544) );
    }

    SECTION("Point 5")
    {
        Point3D point = {36.8315563476, -48.7314244261, 46.3990433125};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.4334414008) );
    }

    SECTION("Point 6")
    {
        Point3D point = {-39.7935307537, 1.0082708909, -48.4975476742};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.7432603141) );
    }

    SECTION("Point 7")
    {
        Point3D point = {39.2152096095, -39.0656723124, -28.995046196};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 8")
    {
        Point3D point = {-11.5997738492, 6.2506795657, 41.5889377667};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.7208217289) );
    }

    SECTION("Point 9")
    {
        Point3D point = {-49.8732305131, -40.7106279818, 48.4922331285};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.8219161918) );
    }

    SECTION("Point 10")
    {
        Point3D point = {-0.5889005263, 15.2143434459, -2.7129801701};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(15.4833842886) );
    }

    SECTION("Point 11")
    {
        Point3D point = {48.969280533, 1.8857173398, -5.5880641358};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(11.4650237679) );
    }

    SECTION("Point 12")
    {
        Point3D point = {-8.4467404794, 45.5121414715, -45.3669887015};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.9062391608) );
    }

    SECTION("Point 13")
    {
        Point3D point = {30.4369597139, -1.7965056709, 48.9445074922};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(11.6561943636) );
    }

    SECTION("Point 14")
    {
        Point3D point = {-44.3057219006, 33.0192715316, 47.8292196048};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 15")
    {
        Point3D point = {-5.7458762805, -43.1324416274, 40.1634508698};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.248644251) );
    }

    SECTION("Point 16")
    {
        Point3D point = {-40.9041742286, 3.1722395463, 4.5642140576};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 17")
    {
        Point3D point = {-31.6555538129, -49.6355080975, -48.3113358721};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.0165966193) );
    }

    SECTION("Point 18")
    {
        Point3D point = {-7.825023475, 48.957493342, 43.3268881837};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 19")
    {
        Point3D point = {-44.9014713033, 48.0409349306, -44.2031802117};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 20")
    {
        Point3D point = {22.8517943401, -29.0174949817, 12.8639449658};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(8.9331640387) );
    }

    SECTION("Point 21")
    {
        Point3D point = {27.4416171375, 22.6609359834, 15.6104371723};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.0698217031) );
    }

    SECTION("Point 22")
    {
        Point3D point = {-30.2095402406, -18.2692825646, 24.9043642426};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(4.0479717921) );
    }

    SECTION("Point 23")
    {
        Point3D point = {48.586275195, 41.7056994008, -14.2714379655};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.0251600483) );
    }

    SECTION("Point 24")
    {
        Point3D point = {18.275234135, -3.5222361579, -22.7704009846};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(7.3009323651) );
    }

    SECTION("Point 25")
    {
        Point3D point = {6.3712748496, -41.5209055373, -17.2412156906};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 26")
    {
        Point3D point = {1.2251402024, -22.859654237, -49.4462563188};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 27")
    {
        Point3D point = {-0.2169936664, 45.897933932, 7.9948189473};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 28")
    {
        Point3D point = {-10.3139075266, 7.8029314325, -49.8249060008};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.1366226545) );
    }

    SECTION("Point 29")
    {
        Point3D point = {42.6518123563, -7.5629428763, -48.0427275868};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.5900401609) );
    }

    SECTION("Point 30")
    {
        Point3D point = {-45.737014057, -28.2994790833, -30.3337322922};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.6872676663) );
    }

    SECTION("Point 31")
    {
        Point3D point = {1.2162083533, -9.9415968917, 14.8779786028};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(5.8379685342) );
    }

    SECTION("Point 32")
    {
        Point3D point = {29.9975908268, 19.9978367751, -14.8495243233};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.1388148975) );
    }

    SECTION("Point 33")
    {
        Point3D point = {-16.2058498553, -12.1394114393, -24.9289664323};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.9213943448) );
    }

    SECTION("Point 34")
    {
        Point3D point = {3.4080482802, -48.8883231296, -43.8845983678};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 35")
    {
        Point3D point = {46.1560620908, -41.6277617643, 8.0593691012};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.3884196883) );
    }

    SECTION("Point 36")
    {
        Point3D point = {20.6848680837, 44.7835938049, -28.666853336};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.9169158965) );
    }

    SECTION("Point 37")
    {
        Point3D point = {-48.5924598754, 31.7137622655, -23.0120238722};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 38")
    {
        Point3D point = {42.6690408926, 21.1015188466, 39.1260346347};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(12.5371008351) );
    }

    SECTION("Point 39")
    {
        Point3D point = {-9.4118899942, 43.2968541949, -20.261988449};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 40")
    {
        Point3D point = {-49.8420792631, 22.5401981606, -49.1368522864};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 41")
    {
        Point3D point = {-26.8725741547, 28.7433622306, 19.6908688032};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 42")
    {
        Point3D point = {-47.2028514823, -47.8318335478, 15.5816407248};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.6148767131) );
    }

    SECTION("Point 43")
    {
        Point3D point = {16.8426310023, 22.1283477601, 43.479231416};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.2391795456) );
    }

    SECTION("Point 44")
    {
        Point3D point = {-46.7923978329, -1.5107623076, 43.8335186307};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 45")
    {
        Point3D point = {46.2800507882, -1.1699410394, 23.3208604033};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(11.4617729153) );
    }

    SECTION("Point 46")
    {
        Point3D point = {-25.1075640671, 16.0016334923, -20.8414799398};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.8481879785) );
    }

    SECTION("Point 47")
    {
        Point3D point = {-41.8020922652, 49.4673997161, 22.9006189261};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 48")
    {
        Point3D point = {-23.2073077342, -44.9300117301, 22.010030305};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.429348821) );
    }

    SECTION("Point 49")
    {
        Point3D point = {37.8241089582, -17.2122999407, -26.5997939168};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.2337126983) );
    }

    SECTION("Point 50")
    {
        Point3D point = {2.5119125622, 24.8735006316, -33.4853518212};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(15.4768984221) );
    }

    SECTION("Point 51")
    {
        Point3D point = {42.3360173555, -22.3200439812, 37.2103834};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(10.3336146389) );
    }

    SECTION("Point 52")
    {
        Point3D point = {24.6305152656, 47.4646406236, 24.1349146581};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.4305844796) );
    }

    SECTION("Point 53")
    {
        Point3D point = {10.5149867295, -15.3445231101, 39.6555222057};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.6202379101) );
    }

    SECTION("Point 54")
    {
        Point3D point = {-0.6580345103, 17.6498819923, 21.9638905823};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(14.6181118337) );
    }

    SECTION("Point 55")
    {
        Point3D point = {21.9565900378, 4.854384649, -46.3083175459};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(6.9412270659) );
    }

    SECTION("Point 56")
    {
        Point3D point = {47.2674666426, 49.1388321385, 13.4482732338};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.0568860076) );
    }

    SECTION("Point 57")
    {
        Point3D point = {-25.7504245153, 24.6689192833, -43.3493452116};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.8806165748) );
    }

    SECTION("Point 58")
    {
        Point3D point = {-30.1640459244, 6.0843163431, 26.2018722371};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 59")
    {
        Point3D point = {29.3484592714, 46.1486992408, -8.1712008725};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.495386228) );
    }

    SECTION("Point 60")
    {
        Point3D point = {48.3516445841, 45.3574198277, -48.7276976457};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.2416276365) );
    }

    SECTION("Point 61")
    {
        Point3D point = {-45.9047522377, -19.3977520193, 2.7042823158};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.947078163) );
    }

    SECTION("Point 62")
    {
        Point3D point = {-48.2935732223, 3.1715559089, -21.2307443243};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 63")
    {
        Point3D point = {19.803537554, 1.7730305678, 2.7095494572};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(11.9476022633) );
    }

    SECTION("Point 64")
    {
        Point3D point = {12.4297294125, -49.8548706993, 4.3646752156};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 65")
    {
        Point3D point = {-23.3992451212, -33.0813171263, -24.6736706582};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.3112300201) );
    }

    SECTION("Point 66")
    {
        Point3D point = {-18.7366216764, 11.0967950249, 8.6394815979};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.7971490178) );
    }

    SECTION("Point 67")
    {
        Point3D point = {-26.0479715076, -28.0749771642, 46.442157075};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.9431939675) );
    }

    SECTION("Point 68")
    {
        Point3D point = {4.5678507325, -22.0657207407, -2.8295629904};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 69")
    {
        Point3D point = {47.0529675004, -49.7124435844, 26.3974328415};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.4338494463) );
    }

    SECTION("Point 70")
    {
        Point3D point = {-26.2698389014, -14.3729289828, -44.3610589459};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.7178330567) );
    }

    SECTION("Point 71")
    {
        Point3D point = {20.170082464, 3.6481735081, 24.0622370383};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(12.3083751719) );
    }

    SECTION("Point 72")
    {
        Point3D point = {-34.3673068957, -16.577460741, -17.3887349513};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.9002353736) );
    }

    SECTION("Point 73")
    {
        Point3D point = {48.5242923249, 23.3141597702, -0.400653505};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(12.5578436829) );
    }

    SECTION("Point 74")
    {
        Point3D point = {46.056583941, -22.4939376919, -6.313336434};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.8354127278) );
    }

    SECTION("Point 75")
    {
        Point3D point = {-16.0577103338, 28.7644069077, 44.1796470406};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 76")
    {
        Point3D point = {-49.6267911045, 23.8918445883, 27.4761605437};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 77")
    {
        Point3D point = {-12.9306530873, 47.8111545079, 23.2428442097};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 78")
    {
        Point3D point = {-38.7378187332, 48.3537688844, -21.8518963418};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 79")
    {
        Point3D point = {-45.7704795039, -29.7681998367, 17.3147712682};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.8500988476) );
    }

    SECTION("Point 80")
    {
        Point3D point = {-3.9792966815, -33.0479217614, 17.9478132482};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.2136732553) );
    }

    SECTION("Point 81")
    {
        Point3D point = {13.1598295938, 48.6314966803, -46.5716411344};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(14.4196487184) );
    }

    SECTION("Point 82")
    {
        Point3D point = {-20.2005061182, -4.6676250895, -2.054497065};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(16.4492271092) );
    }

    SECTION("Point 83")
    {
        Point3D point = {34.9103078642, -42.3299725132, -10.2239740362};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 84")
    {
        Point3D point = {-20.2698560759, 36.5470800952, -1.4485126135};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 85")
    {
        Point3D point = {-41.7693495246, 21.7847013249, -0.4240519136};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 86")
    {
        Point3D point = {41.6287406309, 13.4751787239, -27.3240220627};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(12.3286650542) );
    }

    SECTION("Point 87")
    {
        Point3D point = {-28.3444108806, -48.405116982, 49.6433370279};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.5989999055) );
    }

    SECTION("Point 88")
    {
        Point3D point = {10.5399498532, -40.6306048812, 29.2711783104};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 89")
    {
        Point3D point = {44.1834146595, -21.1933777068, 13.3290060625};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(10.1362505804) );
    }

    SECTION("Point 90")
    {
        Point3D point = {21.063313806, -25.9638172462, -35.295762953};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 91")
    {
        Point3D point = {-41.1744665313, -49.737137556, -16.550619419};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.3804623931) );
    }

    SECTION("Point 92")
    {
        Point3D point = {9.8580917157, 16.7146294223, -20.1967504202};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(14.7825668458) );
    }

    SECTION("Point 93")
    {
        Point3D point = {25.3265625217, -13.2317370098, -7.9272767799};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(9.1345423503) );
    }

    SECTION("Point 94")
    {
        Point3D point = {34.0210880078, -45.3797400908, -47.5821475487};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3) );
    }

    SECTION("Point 95")
    {
        Point3D point = {-44.0322639393, -31.9711322347, -10.1224126109};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.7075048514) );
    }

    SECTION("Point 96")
    {
        Point3D point = {-5.9085791725, -21.4804987756, 32.4249613483};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(4.1234285025) );
    }

    SECTION("Point 97")
    {
        Point3D point = {7.0652345927, 38.8497738581, 43.4287495881};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(13.6988253951) );
    }

    SECTION("Point 98")
    {
        Point3D point = {-29.6910216006, 41.4709048306, 32.7122103342};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(17) );
    }

    SECTION("Point 99")
    {
        Point3D point = {-17.4587137846, -23.6425445292, 4.9992781389};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(3.7460748971) );
    }

    SECTION("Point 100")
    {
        Point3D point = {10.8703085039, 39.8229706054, -12.0107919266};
        projection.Compute(point);
        REQUIRE( projection.Parameter() == Approx(14.1342363212) );
    }
}

TEST_CASE( "Project point on CurveOnSurface", "[PointOnCurveProjection]" ) {
    using namespace ANurbs;

    Pointer<CurveGeometry2D> curveGeometry;
    {
        int degree = 2;
        int nbPoles = 3;
        bool isRational = false;

        curveGeometry = New<CurveGeometry2D>(degree, nbPoles, isRational);

        curveGeometry->SetKnot( 0, 7.0);
        curveGeometry->SetKnot( 1, 7.0);
        curveGeometry->SetKnot( 2, 9.0);
        curveGeometry->SetKnot( 3, 9.0);

        curveGeometry->SetPole(0, {3, 2});
        curveGeometry->SetPole(1, {1, 4});
        curveGeometry->SetPole(2, {2, 5});
    }

    Pointer<SurfaceGeometry3D> surfaceGeometry;
    {
        int degreeU = 2;
        int degreeV = 2;

        int nbPolesU = 3;
        int nbPolesV = 3;

        bool isRational = false;

        auto knotsU = {1, 1, 3, 3};

        surfaceGeometry = New<SurfaceGeometry3D>(degreeU, degreeV, nbPolesU,
            nbPolesV, isRational);

        surfaceGeometry->SetKnotU(0, 1);
        surfaceGeometry->SetKnotU(1, 1);
        surfaceGeometry->SetKnotU(2, 3);
        surfaceGeometry->SetKnotU(3, 3);

        surfaceGeometry->SetKnotV(0, 2);
        surfaceGeometry->SetKnotV(1, 2);
        surfaceGeometry->SetKnotV(2, 6);
        surfaceGeometry->SetKnotV(3, 6);

        surfaceGeometry->SetPole(0, 0, { 0,  0,  3});
        surfaceGeometry->SetPole(0, 1, { 0,  5,  0});
        surfaceGeometry->SetPole(0, 2, { 0, 10,  2});
        surfaceGeometry->SetPole(1, 0, { 5,  0,  5});
        surfaceGeometry->SetPole(1, 1, { 5,  5,  0});
        surfaceGeometry->SetPole(1, 2, { 5, 10,  3});
        surfaceGeometry->SetPole(2, 0, {10,  0,  1});
        surfaceGeometry->SetPole(2, 1, {10,  5, -1});
        surfaceGeometry->SetPole(2, 2, {10, 10,  0});
    }

    Pointer<CurveOnSurface3D> curve = New<CurveOnSurface3D>(curveGeometry,
        surfaceGeometry, curveGeometry->Domain());

    Point3D point = {0, 0, 0};

    PointOnCurveProjection3D projection(curve, 1e-7);
    projection.Compute(point);
    
    REQUIRE( projection.Parameter() == Approx(7.8785895558) );
}