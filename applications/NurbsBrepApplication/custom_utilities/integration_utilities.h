#ifndef INTEGRATION_UTILITIES
#define INTEGRATION_UTILITIES

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/math/special_functions/binomial.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/variables.h"
#include "includes/define.h"
//#include "nurbs_brep_application.h"
//#include "nurbs_brep_application_variables.h"

// ==============================================================================

namespace Kratos
{
  class IntegrationUtilities
  {
  public:
    /// Destructor.
    virtual ~IntegrationUtilities() { }

    /// Pointer definition of MapperUtilities
    KRATOS_CLASS_POINTER_DEFINITION(IntegrationUtilities);

    /**
    * @Author T.Teschemacher, M.Fischer (11/2008)
    * @date   March, 2017
    * @brief   returns location of gauss points in gaussian space from -1 to 1.
    *
    *
    * @param [in] number_of_points  Number of gauss points in a row.
    * @param [out] gaussian_coords  Coords in gaussian space from -1 to 1.
    * @param [out] weights  Weighting of each point respectiveley.
    */
    static void getGaussPointLocations(unsigned int number_of_points,
      Vector& gaussian_coords,
      Vector& weights)
    {
      gaussian_coords.resize(number_of_points);
      weights.resize(number_of_points);
      if (number_of_points == 1)
      {
        gaussian_coords[0] = 0.0;
        weights[0] = 2.0;
      }
      else if (number_of_points == 2)
      {
        gaussian_coords[0] = -1.0 / (std::sqrt(3.0));
        gaussian_coords[1] = 1.0 / (std::sqrt(3.0));

        weights[0] = 1.0;
        weights[1] = 1.0;
      }
      else if (number_of_points == 3)
      {
        gaussian_coords[0] = -sqrt(0.6);
        gaussian_coords[1] = 0.0;
        gaussian_coords[2] = sqrt(0.6);

        weights[0] = 5.0 / 9.0;
        weights[1] = 8.0 / 9.0;
        weights[2] = 5.0 / 9.0;
      }
      else if (number_of_points == 4)
      {
        gaussian_coords[0] = -0.86113631159405257524;
        gaussian_coords[1] = -0.33998104358485626481;
        gaussian_coords[2] = 0.33998104358485626481;
        gaussian_coords[3] = 0.86113631159405257524;

        weights[0] = 0.34785484513745385736;
        weights[1] = 0.65214515486254614264;
        weights[2] = 0.65214515486254614264;
        weights[3] = 0.34785484513745385736;
      }
      else if (number_of_points == 5)
      {
        gaussian_coords[0] = -0.90617984593866399282;
        gaussian_coords[1] = -0.53846931010568309105;
        gaussian_coords[2] = 0;
        gaussian_coords[3] = 0.53846931010568309105;
        gaussian_coords[4] = 0.90617984593866399282;

        weights[0] = 0.23692688505618908749;
        weights[1] = 0.47862867049936646808;
        weights[2] = 0.56888888888888888888;
        weights[3] = 0.47862867049936646808;
        weights[4] = 0.23692688505618908749;
      }
      else if (number_of_points == 6)
      {
        gaussian_coords[0] = -0.9324695142031520;
        gaussian_coords[1] = -0.6612093864662645;
        gaussian_coords[2] = -0.2386191860831969;
        gaussian_coords[3] = 0.2386191860831969;
        gaussian_coords[4] = 0.6612093864662645;
        gaussian_coords[5] = 0.9324695142031520;

        weights[0] = .1713244923791703;
        weights[1] = .3607615730481386;
        weights[2] = .4679139345726911;
        weights[3] = .4679139345726911;
        weights[4] = .3607615730481386;
        weights[5] = .1713244923791703;
      }
      else if (number_of_points == 7)
      {
        gaussian_coords[0] = -0.9491079123427585;
        gaussian_coords[1] = -0.7415311855993944;
        gaussian_coords[2] = -0.4058451513773972;
        gaussian_coords[3] = 0;
        gaussian_coords[4] = 0.4058451513773972;
        gaussian_coords[5] = 0.7415311855993944;
        gaussian_coords[6] = 0.9491079123427585;

        weights[0] = 0.1294849661688697;
        weights[1] = 0.2797053914892767;
        weights[2] = 0.3818300505051189;
        weights[3] = 0.4179591836734694;
        weights[4] = 0.3818300505051189;
        weights[5] = 0.2797053914892767;
        weights[6] = 0.1294849661688697;
      }
      else if (number_of_points == 8)
      {
        gaussian_coords[0] = -0.9602898564975362;
        gaussian_coords[1] = -0.7966664774136267;
        gaussian_coords[2] = -0.5255324099163290;
        gaussian_coords[3] = -0.1834346424956498;
        gaussian_coords[4] = 0.1834346424956498;
        gaussian_coords[5] = 0.5255324099163290;
        gaussian_coords[6] = 0.7966664774136267;
        gaussian_coords[7] = 0.9602898564975362;

        weights[0] = 0.1012285362903763;
        weights[1] = 0.2223810344533745;
        weights[2] = 0.3137066458778873;
        weights[3] = 0.3626837833783620;
        weights[4] = 0.3626837833783620;
        weights[5] = 0.3137066458778873;
        weights[6] = 0.2223810344533745;
        weights[7] = 0.1012285362903763;
      }
      else if (number_of_points == 9)
      {
        gaussian_coords[0] = -0.9681602395076261;
        gaussian_coords[1] = -0.8360311073266358;
        gaussian_coords[2] = -0.6133714327005904;
        gaussian_coords[3] = -0.3242534234038089;
        gaussian_coords[4] = 0.0;
        gaussian_coords[5] = 0.3242534234038089;
        gaussian_coords[6] = 0.6133714327005904;
        gaussian_coords[7] = 0.8360311073266358;
        gaussian_coords[8] = 0.9681602395076261;

        weights[0] = 0.0812743883615744;
        weights[1] = 0.1806481606948574;
        weights[2] = 0.2606106964029354;
        weights[3] = 0.3123470770400029;
        weights[4] = 0.3302393550012598;
        weights[5] = 0.3123470770400028;
        weights[6] = 0.2606106964029355;
        weights[7] = 0.1806481606948574;
        weights[8] = 0.0812743883615744;
      }
      else if (number_of_points == 10)
      {
        gaussian_coords[0] = -0.973906528517172;
        gaussian_coords[1] = -0.865063366688985;
        gaussian_coords[2] = -0.679409568299024;
        gaussian_coords[3] = -0.433395394129247;
        gaussian_coords[4] = -0.148874338981631;
        gaussian_coords[5] = 0.148874338981631;
        gaussian_coords[6] = 0.433395394129247;
        gaussian_coords[7] = 0.679409568299024;
        gaussian_coords[8] = 0.865063366688985;
        gaussian_coords[9] = 0.973906528517172;

        weights[0] = 0.066671344308688;
        weights[1] = 0.149451349150581;
        weights[2] = 0.219086362515982;
        weights[3] = 0.269266719309996;
        weights[4] = 0.295524224714753;
        weights[5] = 0.295524224714753;
        weights[6] = 0.269266719309996;
        weights[7] = 0.219086362515982;
        weights[8] = 0.149451349150581;
        weights[9] = 0.066671344308688;
      }
      else if (number_of_points == 11)   // added by R.Schmidt 09/2011
      {
        gaussian_coords[0] = 9.782286581460570e-01;
        gaussian_coords[1] = 8.870625997680953e-01;
        gaussian_coords[2] = 7.301520055740493e-01;
        gaussian_coords[3] = 5.190961292068118e-01;
        gaussian_coords[4] = 2.695431559523450e-01;
        gaussian_coords[5] = -0.000000000000000e-01;
        gaussian_coords[6] = -2.695431559523450e-01;
        gaussian_coords[7] = -5.190961292068118e-01;
        gaussian_coords[8] = -7.301520055740493e-01;
        gaussian_coords[9] = -8.870625997680953e-01;
        gaussian_coords[10] = -9.782286581460570e-01;

        weights[0] = 5.566856711617367e-02;
        weights[1] = 1.255803694649046e-01;
        weights[2] = 1.862902109277343e-01;
        weights[3] = 2.331937645919905e-01;
        weights[4] = 2.628045445102467e-01;
        weights[5] = 2.729250867779006e-01;
        weights[6] = 2.628045445102467e-01;
        weights[7] = 2.331937645919905e-01;
        weights[8] = 1.862902109277343e-01;
        weights[9] = 1.255803694649046e-01;
        weights[10] = 5.566856711617367e-02;
      }
      else if (number_of_points == 12)// added by R.Schmidt 09/2011
      {
        gaussian_coords[0] = 9.815606342467193e-01;
        gaussian_coords[1] = 9.041172563704749e-01;
        gaussian_coords[2] = 7.699026741943047e-01;
        gaussian_coords[3] = 5.873179542866174e-01;
        gaussian_coords[4] = 3.678314989981802e-01;
        gaussian_coords[5] = 1.252334085114689e-01;
        gaussian_coords[6] = -1.252334085114689e-01;
        gaussian_coords[7] = -3.678314989981802e-01;
        gaussian_coords[8] = -5.873179542866174e-01;
        gaussian_coords[9] = -7.699026741943047e-01;
        gaussian_coords[10] = -9.041172563704749e-01;
        gaussian_coords[11] = -9.815606342467193e-01;

        weights[0] = 4.717533638651183e-02;
        weights[1] = 1.069393259953184e-01;
        weights[2] = 1.600783285433462e-01;
        weights[3] = 2.031674267230659e-01;
        weights[4] = 2.334925365383548e-01;
        weights[5] = 2.491470458134028e-01;
        weights[6] = 2.491470458134028e-01;
        weights[7] = 2.334925365383548e-01;
        weights[8] = 2.031674267230659e-01;
        weights[9] = 1.600783285433462e-01;
        weights[10] = 1.069393259953184e-01;
        weights[11] = 4.717533638651183e-02;
      }
      else if (number_of_points == 13)// added by R.Schmidt 09/2011
      {
        gaussian_coords[0] = 9.841830547185881e-01;
        gaussian_coords[1] = 9.175983992229780e-01;
        gaussian_coords[2] = 8.015780907333099e-01;
        gaussian_coords[3] = 6.423493394403402e-01;
        gaussian_coords[4] = 4.484927510364469e-01;
        gaussian_coords[5] = 2.304583159551348e-01;
        gaussian_coords[6] = -0.000000000000000e-01;
        gaussian_coords[7] = -2.304583159551348e-01;
        gaussian_coords[8] = -4.484927510364469e-01;
        gaussian_coords[9] = -6.423493394403402e-01;
        gaussian_coords[10] = -8.015780907333099e-01;
        gaussian_coords[11] = -9.175983992229780e-01;
        gaussian_coords[12] = -9.841830547185881e-01;

        weights[0] = 4.048400476531588e-02;
        weights[1] = 9.212149983772845e-02;
        weights[2] = 1.388735102197872e-01;
        weights[3] = 1.781459807619457e-01;
        weights[4] = 2.078160475368885e-01;
        weights[5] = 2.262831802628972e-01;
        weights[6] = 2.325515532308739e-01;
        weights[7] = 2.262831802628972e-01;
        weights[8] = 2.078160475368885e-01;
        weights[9] = 1.781459807619457e-01;
        weights[10] = 1.388735102197872e-01;
        weights[11] = 9.212149983772845e-02;
        weights[12] = 4.048400476531588e-02;
      }
      else if (number_of_points == 14)// added by R.Schmidt 09/2011
      {
        gaussian_coords[0] = 9.862838086968123e-01;
        gaussian_coords[1] = 9.284348836635735e-01;
        gaussian_coords[2] = 8.272013150697650e-01;
        gaussian_coords[3] = 6.872929048116855e-01;
        gaussian_coords[4] = 5.152486363581541e-01;
        gaussian_coords[5] = 3.191123689278898e-01;
        gaussian_coords[6] = 1.080549487073437e-01;
        gaussian_coords[7] = -1.080549487073437e-01;
        gaussian_coords[8] = -3.191123689278898e-01;
        gaussian_coords[9] = -5.152486363581541e-01;
        gaussian_coords[10] = -6.872929048116855e-01;
        gaussian_coords[11] = -8.272013150697650e-01;
        gaussian_coords[12] = -9.284348836635735e-01;
        gaussian_coords[13] = -9.862838086968123e-01;

        weights[0] = 3.511946033175186e-02;
        weights[1] = 8.015808715976021e-02;
        weights[2] = 1.215185706879032e-01;
        weights[3] = 1.572031671581935e-01;
        weights[4] = 1.855383974779378e-01;
        weights[5] = 2.051984637212956e-01;
        weights[6] = 2.152638534631578e-01;
        weights[7] = 2.152638534631578e-01;
        weights[8] = 2.051984637212956e-01;
        weights[9] = 1.855383974779378e-01;
        weights[10] = 1.572031671581935e-01;
        weights[11] = 1.215185706879032e-01;
        weights[12] = 8.015808715976021e-02;
        weights[13] = 3.511946033175186e-02;
      }
      else if (number_of_points == 15) //from Hughes isogeometric code, not checked
      {
        gaussian_coords[0] = -0.9879925180204854;
        gaussian_coords[1] = -0.9372733924007059;
        gaussian_coords[2] = -0.8482065834104272;
        gaussian_coords[3] = -0.7244177313601700;
        gaussian_coords[4] = -0.5709721726085388;
        gaussian_coords[5] = -0.3941513470775634;
        gaussian_coords[6] = -0.2011940939974345;
        gaussian_coords[7] = 0;
        gaussian_coords[8] = 0.2011940939974345;
        gaussian_coords[9] = 0.3941513470775634;
        gaussian_coords[10] = 0.5709721726085388;
        gaussian_coords[11] = 0.7244177313601700;
        gaussian_coords[12] = 0.8482065834104272;
        gaussian_coords[13] = 0.9372733924007059;
        gaussian_coords[14] = 0.9879925180204854;

        weights[0] = 0.03075324199611807;
        weights[1] = 0.07036604748811134;
        weights[2] = 0.1071592204671351;
        weights[3] = 0.1395706779261761;
        weights[4] = 0.1662692058169852;
        weights[5] = 0.1861610000155741;
        weights[6] = 0.1984314853271374;
        weights[7] = 0.2025782419255562;
        weights[8] = 0.1984314853271374;
        weights[9] = 0.1861610000155741;
        weights[10] = 0.1662692058169852;
        weights[11] = 0.1395706779261761;
        weights[12] = 0.1071592204671351;
        weights[13] = 0.07036604748811134;
        weights[14] = 0.03075324199611807;
      }
      else if (number_of_points == 64) //from http://pomax.github.io/bezierinfo/legendre-gauss.html, not checked
      {
        gaussian_coords[0] = -0.9993050417357720;
        gaussian_coords[1] = -0.9963401167719550;
        gaussian_coords[2] = -0.9910133714767440;
        gaussian_coords[3] = -0.9833362538846260;
        gaussian_coords[4] = -0.9733268277899110;
        gaussian_coords[5] = -0.9610087996520530;
        gaussian_coords[6] = -0.9464113748584020;
        gaussian_coords[7] = -0.9295691721319390;
        gaussian_coords[8] = -0.9105221370785020;
        gaussian_coords[9] = -0.8893154459951140;
        gaussian_coords[10] = -0.8659993981540920;
        gaussian_coords[11] = -0.8406292962525800;
        gaussian_coords[12] = -0.8132653151227970;
        gaussian_coords[13] = -0.7839723589433410;
        gaussian_coords[14] = -0.7528199072605310;
        gaussian_coords[15] = -0.7198818501716100;
        gaussian_coords[16] = -0.6852363130542330;
        gaussian_coords[17] = -0.6489654712546570;
        gaussian_coords[18] = -0.6111553551723930;
        gaussian_coords[19] = -0.5718956462026340;
        gaussian_coords[20] = -0.5312794640198940;
        gaussian_coords[21] = -0.4894031457070530;
        gaussian_coords[22] = -0.4463660172534640;
        gaussian_coords[23] = -0.4022701579639910;
        gaussian_coords[24] = -0.3572201583376680;
        gaussian_coords[25] = -0.3113228719902110;
        gaussian_coords[26] = -0.2646871622087670;
        gaussian_coords[27] = -0.2174236437400070;
        gaussian_coords[28] = -0.1696444204239920;
        gaussian_coords[29] = -0.1214628192961200;
        gaussian_coords[30] = -0.0729931217877990;
        gaussian_coords[31] = -0.0243502926634244;
        gaussian_coords[32] = 0.0243502926634244;
        gaussian_coords[33] = 0.0729931217877990;
        gaussian_coords[34] = 0.1214628192961200;
        gaussian_coords[35] = 0.1696444204239920;
        gaussian_coords[36] = 0.2174236437400070;
        gaussian_coords[37] = 0.2646871622087670;
        gaussian_coords[38] = 0.3113228719902110;
        gaussian_coords[39] = 0.3572201583376680;
        gaussian_coords[40] = 0.4022701579639910;
        gaussian_coords[41] = 0.4463660172534640;
        gaussian_coords[42] = 0.4894031457070530;
        gaussian_coords[43] = 0.5312794640198940;
        gaussian_coords[44] = 0.5718956462026340;
        gaussian_coords[45] = 0.6111553551723930;
        gaussian_coords[46] = 0.6489654712546570;
        gaussian_coords[47] = 0.6852363130542330;
        gaussian_coords[48] = 0.7198818501716100;
        gaussian_coords[49] = 0.7528199072605310;
        gaussian_coords[50] = 0.7839723589433410;
        gaussian_coords[51] = 0.8132653151227970;
        gaussian_coords[52] = 0.8406292962525800;
        gaussian_coords[53] = 0.8659993981540920;
        gaussian_coords[54] = 0.8893154459951140;
        gaussian_coords[55] = 0.9105221370785020;
        gaussian_coords[56] = 0.9295691721319390;
        gaussian_coords[57] = 0.9464113748584020;
        gaussian_coords[58] = 0.9610087996520530;
        gaussian_coords[59] = 0.9733268277899110;
        gaussian_coords[60] = 0.9833362538846260;
        gaussian_coords[61] = 0.9910133714767440;
        gaussian_coords[62] = 0.9963401167719550;
        gaussian_coords[63] = 0.9993050417357720;

        weights[0] = 0.0017832807216964;
        weights[1] = 0.0041470332605625;
        weights[2] = 0.0065044579689784;
        weights[3] = 0.0088467598263639;
        weights[4] = 0.0111681394601311;
        weights[5] = 0.0134630478967186;
        weights[6] = 0.0157260304760247;
        weights[7] = 0.0179517157756973;
        weights[8] = 0.0201348231535302;
        weights[9] = 0.0222701738083833;
        weights[10] = 0.0243527025687109;
        weights[11] = 0.0263774697150547;
        weights[12] = 0.0283396726142595;
        weights[13] = 0.0302346570724025;
        weights[14] = 0.0320579283548516;
        weights[15] = 0.0338051618371416;
        weights[16] = 0.0354722132568824;
        weights[17] = 0.0370551285402400;
        weights[18] = 0.0385501531786156;
        weights[19] = 0.0399537411327203;
        weights[20] = 0.0412625632426235;
        weights[21] = 0.0424735151236536;
        weights[22] = 0.0435837245293235;
        weights[23] = 0.0445905581637566;
        weights[24] = 0.0454916279274181;
        weights[25] = 0.0462847965813144;
        weights[26] = 0.0469681828162100;
        weights[27] = 0.0475401657148303;
        weights[28] = 0.0479993885964583;
        weights[29] = 0.0483447622348030;
        weights[30] = 0.0485754674415034;
        weights[31] = 0.0486909570091397;
        weights[32] = 0.0486909570091397;
        weights[33] = 0.0485754674415034;
        weights[34] = 0.0483447622348030;
        weights[35] = 0.0479993885964583;
        weights[36] = 0.0475401657148303;
        weights[37] = 0.0469681828162100;
        weights[38] = 0.0462847965813144;
        weights[39] = 0.0454916279274181;
        weights[40] = 0.0445905581637566;
        weights[41] = 0.0435837245293235;
        weights[42] = 0.0424735151236536;
        weights[43] = 0.0412625632426235;
        weights[44] = 0.0399537411327203;
        weights[45] = 0.0385501531786156;
        weights[46] = 0.0370551285402400;
        weights[47] = 0.0354722132568824;
        weights[48] = 0.0338051618371416;
        weights[49] = 0.0320579283548516;
        weights[50] = 0.0302346570724025;
        weights[51] = 0.0283396726142595;
        weights[52] = 0.0263774697150547;
        weights[53] = 0.0243527025687109;
        weights[54] = 0.0222701738083833;
        weights[55] = 0.0201348231535302;
        weights[56] = 0.0179517157756973;
        weights[57] = 0.0157260304760247;
        weights[58] = 0.0134630478967186;
        weights[59] = 0.0111681394601311;
        weights[60] = 0.0088467598263639;
        weights[61] = 0.0065044579689784;
        weights[62] = 0.0041470332605625;
        weights[63] = 0.0017832807216964;
      }
      else
      {
        KRATOS_THROW_ERROR(std::logic_error, "IntegrationUtilities::GetFace: face_id does not exist", "");
      }
    }

    /**
    * @Author T.Oberbichler
    * @date   March, 2017
    * @brief   returns location of gauss points for triangles in gaussian space from -1 to 1.
    *
    *
    * @param [in] number_of_points  Number of gauss points in a row.
    * @param [out] gaussian_coords  Coords in gaussian space from -1 to 1.
    * @param [out] weights  Weighting of each point respectiveley.
    */
    static void getGaussPointLocationsTriangle(unsigned int degree,
      Matrix& gaussian_coords,
      Vector& weights)
    {
      std::cout << "degree: " << degree << std::endl;
      switch (degree) {
      case 1:
        gaussian_coords.resize(1, 2);
        weights.resize(1);

        gaussian_coords(0, 0) = 0.33333333333333;
        gaussian_coords(0, 1) = 0.33333333333333;
        weights[0] = 1.00000000000000;
        return;
      case 2:
        gaussian_coords.resize(3, 2);
        weights.resize(3);

        gaussian_coords(0, 0) = 0.16666666666667;
        gaussian_coords(0, 1) = 0.16666666666667;
        gaussian_coords(1, 0) = 0.16666666666667;
        gaussian_coords(1, 1) = 0.66666666666667;
        gaussian_coords(2, 0) = 0.66666666666667;
        gaussian_coords(2, 1) = 0.16666666666667;
        weights[0] = 0.33333333333333;
        weights[1] = 0.33333333333333;
        weights[2] = 0.33333333333333;
        return;
      case 3:
        gaussian_coords.resize(4, 2);
        weights.resize(4);

        gaussian_coords(0, 0) = 0.33333333333333;
        gaussian_coords(0, 1) = 0.33333333333333;
        gaussian_coords(1, 0) = 0.20000000000000;
        gaussian_coords(1, 1) = 0.20000000000000;
        gaussian_coords(2, 0) = 0.20000000000000;
        gaussian_coords(2, 1) = 0.60000000000000;
        gaussian_coords(3, 0) = 0.60000000000000;
        gaussian_coords(3, 1) = 0.20000000000000;

        weights[0] = -0.56250000000000;
        weights[1] = 0.52083333333333;
        weights[2] = 0.52083333333333;
        weights[3] = 0.52083333333333;
        return;
      //case 4:
      //  return{
      //    IntegrationPoint(Vector2d(0.44594849091597, 0.44594849091597), 0.22338158967801),
      //    IntegrationPoint(Vector2d(0.44594849091597, 0.10810301816807), 0.22338158967801),
      //    IntegrationPoint(Vector2d(0.10810301816807, 0.44594849091597), 0.22338158967801),
      //    IntegrationPoint(Vector2d(0.09157621350977, 0.09157621350977), 0.10995174365532),
      //    IntegrationPoint(Vector2d(0.09157621350977, 0.81684757298046), 0.10995174365532),
      //    IntegrationPoint(Vector2d(0.81684757298046, 0.09157621350977), 0.10995174365532)
      //  };
      //case 5:
      //  return{
      //    IntegrationPoint(Vector2d(0.33333333333333, 0.33333333333333), 0.22500000000000),
      //    IntegrationPoint(Vector2d(0.47014206410511, 0.47014206410511), 0.13239415278851),
      //    IntegrationPoint(Vector2d(0.47014206410511, 0.05971587178977), 0.13239415278851),
      //    IntegrationPoint(Vector2d(0.05971587178977, 0.47014206410511), 0.13239415278851),
      //    IntegrationPoint(Vector2d(0.10128650732346, 0.10128650732346), 0.12593918054483),
      //    IntegrationPoint(Vector2d(0.10128650732346, 0.79742698535309), 0.12593918054483),
      //    IntegrationPoint(Vector2d(0.79742698535309, 0.10128650732346), 0.12593918054483)
      //  };
      //case 6:
      //  return{
      //    IntegrationPoint(Vector2d(0.24928674517091, 0.24928674517091), 0.11678627572638),
      //    IntegrationPoint(Vector2d(0.24928674517091, 0.50142650965818), 0.11678627572638),
      //    IntegrationPoint(Vector2d(0.50142650965818, 0.24928674517091), 0.11678627572638),
      //    IntegrationPoint(Vector2d(0.06308901449150, 0.06308901449150), 0.05084490637021),
      //    IntegrationPoint(Vector2d(0.06308901449150, 0.87382197101700), 0.05084490637021),
      //    IntegrationPoint(Vector2d(0.87382197101700, 0.06308901449150), 0.05084490637021),
      //    IntegrationPoint(Vector2d(0.31035245103378, 0.63650249912140), 0.08285107561837),
      //    IntegrationPoint(Vector2d(0.63650249912140, 0.05314504984482), 0.08285107561837),
      //    IntegrationPoint(Vector2d(0.05314504984482, 0.31035245103378), 0.08285107561837),
      //    IntegrationPoint(Vector2d(0.63650249912140, 0.31035245103378), 0.08285107561837),
      //    IntegrationPoint(Vector2d(0.31035245103378, 0.05314504984482), 0.08285107561837),
      //    IntegrationPoint(Vector2d(0.05314504984482, 0.63650249912140), 0.08285107561837)
      //  };
      //case 7:
      //  return{
      //    IntegrationPoint(Vector2d(0.33333333333333, 0.33333333333333), -0.14957004446768),
      //    IntegrationPoint(Vector2d(0.26034596607904, 0.26034596607904), 0.17561525743321),
      //    IntegrationPoint(Vector2d(0.26034596607904, 0.47930806784192), 0.17561525743321),
      //    IntegrationPoint(Vector2d(0.47930806784192, 0.26034596607904), 0.17561525743321),
      //    IntegrationPoint(Vector2d(0.06513010290222, 0.06513010290222), 0.05334723560884),
      //    IntegrationPoint(Vector2d(0.06513010290222, 0.86973979419557), 0.05334723560884),
      //    IntegrationPoint(Vector2d(0.86973979419557, 0.06513010290222), 0.05334723560884),
      //    IntegrationPoint(Vector2d(0.31286549600487, 0.63844418856981), 0.07711376089026),
      //    IntegrationPoint(Vector2d(0.63844418856981, 0.04869031542532), 0.07711376089026),
      //    IntegrationPoint(Vector2d(0.04869031542532, 0.31286549600487), 0.07711376089026),
      //    IntegrationPoint(Vector2d(0.63844418856981, 0.31286549600487), 0.07711376089026),
      //    IntegrationPoint(Vector2d(0.31286549600487, 0.04869031542532), 0.07711376089026),
      //    IntegrationPoint(Vector2d(0.04869031542532, 0.63844418856981), 0.07711376089026)
      //  };
      //case 8:
      //  return{
      //    IntegrationPoint(Vector2d(0.33333333333333, 0.33333333333333), 0.14431560767779),
      //    IntegrationPoint(Vector2d(0.45929258829272, 0.45929258829272), 0.09509163426728),
      //    IntegrationPoint(Vector2d(0.45929258829272, 0.08141482341455), 0.09509163426728),
      //    IntegrationPoint(Vector2d(0.08141482341455, 0.45929258829272), 0.09509163426728),
      //    IntegrationPoint(Vector2d(0.17056930775176, 0.17056930775176), 0.10321737053472),
      //    IntegrationPoint(Vector2d(0.17056930775176, 0.65886138449648), 0.10321737053472),
      //    IntegrationPoint(Vector2d(0.65886138449648, 0.17056930775176), 0.10321737053472),
      //    IntegrationPoint(Vector2d(0.05054722831703, 0.05054722831703), 0.03245849762320),
      //    IntegrationPoint(Vector2d(0.05054722831703, 0.89890554336594), 0.03245849762320),
      //    IntegrationPoint(Vector2d(0.89890554336594, 0.05054722831703), 0.03245849762320),
      //    IntegrationPoint(Vector2d(0.26311282963464, 0.72849239295540), 0.02723031417443),
      //    IntegrationPoint(Vector2d(0.72849239295540, 0.00839477740996), 0.02723031417443),
      //    IntegrationPoint(Vector2d(0.00839477740996, 0.26311282963464), 0.02723031417443),
      //    IntegrationPoint(Vector2d(0.72849239295540, 0.26311282963464), 0.02723031417443),
      //    IntegrationPoint(Vector2d(0.26311282963464, 0.00839477740996), 0.02723031417443),
      //    IntegrationPoint(Vector2d(0.00839477740996, 0.72849239295540), 0.02723031417443)
      //  };
      //case 9:
      //  return{
      //    IntegrationPoint(Vector2d(0.33333333333333, 0.33333333333333), 0.09713579628280),
      //    IntegrationPoint(Vector2d(0.48968251919874, 0.48968251919874), 0.03133470022714),
      //    IntegrationPoint(Vector2d(0.48968251919874, 0.02063496160252), 0.03133470022714),
      //    IntegrationPoint(Vector2d(0.02063496160252, 0.48968251919874), 0.03133470022714),
      //    IntegrationPoint(Vector2d(0.43708959149294, 0.43708959149294), 0.07782754100477),
      //    IntegrationPoint(Vector2d(0.43708959149294, 0.12582081701413), 0.07782754100477),
      //    IntegrationPoint(Vector2d(0.12582081701413, 0.43708959149294), 0.07782754100477),
      //    IntegrationPoint(Vector2d(0.18820353561903, 0.18820353561903), 0.07964773892721),
      //    IntegrationPoint(Vector2d(0.18820353561903, 0.62359292876193), 0.07964773892721),
      //    IntegrationPoint(Vector2d(0.62359292876193, 0.18820353561903), 0.07964773892721),
      //    IntegrationPoint(Vector2d(0.04472951339445, 0.04472951339445), 0.02557767565870),
      //    IntegrationPoint(Vector2d(0.04472951339445, 0.91054097321109), 0.02557767565870),
      //    IntegrationPoint(Vector2d(0.91054097321109, 0.04472951339445), 0.02557767565870),
      //    IntegrationPoint(Vector2d(0.22196298916077, 0.74119859878450), 0.04328353937729),
      //    IntegrationPoint(Vector2d(0.74119859878450, 0.03683841205474), 0.04328353937729),
      //    IntegrationPoint(Vector2d(0.03683841205474, 0.22196298916077), 0.04328353937729),
      //    IntegrationPoint(Vector2d(0.74119859878450, 0.22196298916077), 0.04328353937729),
      //    IntegrationPoint(Vector2d(0.22196298916077, 0.03683841205474), 0.04328353937729),
      //    IntegrationPoint(Vector2d(0.03683841205474, 0.74119859878450), 0.04328353937729)
      //  };
      //case 10:
      //  return{
      //    IntegrationPoint(Vector2d(0.33333333333333, 0.33333333333333), 0.09081799038275),
      //    IntegrationPoint(Vector2d(0.48557763338366, 0.48557763338366), 0.03672595775647),
      //    IntegrationPoint(Vector2d(0.48557763338366, 0.02884473323269), 0.03672595775647),
      //    IntegrationPoint(Vector2d(0.02884473323269, 0.48557763338366), 0.03672595775647),
      //    IntegrationPoint(Vector2d(0.10948157548504, 0.10948157548504), 0.04532105943553),
      //    IntegrationPoint(Vector2d(0.10948157548504, 0.78103684902993), 0.04532105943553),
      //    IntegrationPoint(Vector2d(0.78103684902993, 0.10948157548504), 0.04532105943553),
      //    IntegrationPoint(Vector2d(0.30793983876412, 0.55035294182100), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.55035294182100, 0.14170721941488), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.14170721941488, 0.30793983876412), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.55035294182100, 0.30793983876412), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.30793983876412, 0.14170721941488), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.14170721941488, 0.55035294182100), 0.07275791684542),
      //    IntegrationPoint(Vector2d(0.24667256063990, 0.72832390459741), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.72832390459741, 0.02500353476269), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.02500353476269, 0.24667256063990), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.72832390459741, 0.24667256063990), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.24667256063990, 0.02500353476269), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.02500353476269, 0.72832390459741), 0.02832724253106),
      //    IntegrationPoint(Vector2d(0.06680325101220, 0.92365593358750), 0.00942166696373),
      //    IntegrationPoint(Vector2d(0.92365593358750, 0.00954081540030), 0.00942166696373),
      //    IntegrationPoint(Vector2d(0.00954081540030, 0.06680325101220), 0.00942166696373),
      //    IntegrationPoint(Vector2d(0.92365593358750, 0.06680325101220), 0.00942166696373),
      //    IntegrationPoint(Vector2d(0.06680325101220, 0.00954081540030), 0.00942166696373),
      //    IntegrationPoint(Vector2d(0.00954081540030, 0.92365593358750), 0.00942166696373)
      //  };
      //case 11:
      //  return{
      //    IntegrationPoint(Vector2d(0.53461104827076, 0.53461104827076), 0.00092700632896),
      //    IntegrationPoint(Vector2d(0.53461104827076, -0.06922209654152), 0.00092700632896),
      //    IntegrationPoint(Vector2d(-0.06922209654152, 0.53461104827076), 0.00092700632896),
      //    IntegrationPoint(Vector2d(0.39896930296585, 0.39896930296585), 0.07714953491481),
      //    IntegrationPoint(Vector2d(0.39896930296585, 0.20206139406829), 0.07714953491481),
      //    IntegrationPoint(Vector2d(0.20206139406829, 0.39896930296585), 0.07714953491481),
      //    IntegrationPoint(Vector2d(0.20330990043128, 0.20330990043128), 0.05932297738077),
      //    IntegrationPoint(Vector2d(0.20330990043128, 0.59338019913744), 0.05932297738077),
      //    IntegrationPoint(Vector2d(0.59338019913744, 0.20330990043128), 0.05932297738077),
      //    IntegrationPoint(Vector2d(0.11935091228258, 0.11935091228258), 0.03618454050342),
      //    IntegrationPoint(Vector2d(0.11935091228258, 0.76129817543484), 0.03618454050342),
      //    IntegrationPoint(Vector2d(0.76129817543484, 0.11935091228258), 0.03618454050342),
      //    IntegrationPoint(Vector2d(0.03236494811128, 0.03236494811128), 0.01365973100268),
      //    IntegrationPoint(Vector2d(0.03236494811128, 0.93527010377745), 0.01365973100268),
      //    IntegrationPoint(Vector2d(0.93527010377745, 0.03236494811128), 0.01365973100268),
      //    IntegrationPoint(Vector2d(0.35662064826129, 0.59320121342821), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.59320121342821, 0.05017813831050), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.05017813831050, 0.35662064826129), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.59320121342821, 0.35662064826129), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.35662064826129, 0.05017813831050), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.05017813831050, 0.59320121342821), 0.05233711196220),
      //    IntegrationPoint(Vector2d(0.17148898030404, 0.80748900315979), 0.02070765963914),
      //    IntegrationPoint(Vector2d(0.80748900315979, 0.02102201653617), 0.02070765963914),
      //    IntegrationPoint(Vector2d(0.02102201653617, 0.17148898030404), 0.02070765963914),
      //    IntegrationPoint(Vector2d(0.80748900315979, 0.17148898030404), 0.02070765963914),
      //    IntegrationPoint(Vector2d(0.17148898030404, 0.02102201653617), 0.02070765963914),
      //    IntegrationPoint(Vector2d(0.02102201653617, 0.80748900315979), 0.02070765963914)
      //  };
      //case 12:
      //  return{
      //    IntegrationPoint(Vector2d(0.48821738977381, 0.48821738977381), 0.02573106644045),
      //    IntegrationPoint(Vector2d(0.48821738977381, 0.02356522045239), 0.02573106644045),
      //    IntegrationPoint(Vector2d(0.02356522045239, 0.48821738977381), 0.02573106644045),
      //    IntegrationPoint(Vector2d(0.43972439229446, 0.43972439229446), 0.04369254453804),
      //    IntegrationPoint(Vector2d(0.43972439229446, 0.12055121541108), 0.04369254453804),
      //    IntegrationPoint(Vector2d(0.12055121541108, 0.43972439229446), 0.04369254453804),
      //    IntegrationPoint(Vector2d(0.27121038501212, 0.27121038501212), 0.06285822421789),
      //    IntegrationPoint(Vector2d(0.27121038501212, 0.45757922997577), 0.06285822421789),
      //    IntegrationPoint(Vector2d(0.45757922997577, 0.27121038501212), 0.06285822421789),
      //    IntegrationPoint(Vector2d(0.12757614554159, 0.12757614554159), 0.03479611293071),
      //    IntegrationPoint(Vector2d(0.12757614554159, 0.74484770891683), 0.03479611293071),
      //    IntegrationPoint(Vector2d(0.74484770891683, 0.12757614554159), 0.03479611293071),
      //    IntegrationPoint(Vector2d(0.02131735045321, 0.02131735045321), 0.00616626105156),
      //    IntegrationPoint(Vector2d(0.02131735045321, 0.95736529909358), 0.00616626105156),
      //    IntegrationPoint(Vector2d(0.95736529909358, 0.02131735045321), 0.00616626105156),
      //    IntegrationPoint(Vector2d(0.27571326968551, 0.60894323577979), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.60894323577979, 0.11534349453470), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.11534349453470, 0.27571326968551), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.60894323577979, 0.27571326968551), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.27571326968551, 0.11534349453470), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.11534349453470, 0.60894323577979), 0.04037155776638),
      //    IntegrationPoint(Vector2d(0.28132558098994, 0.69583608678780), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.69583608678780, 0.02283833222226), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.02283833222226, 0.28132558098994), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.69583608678780, 0.28132558098994), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.28132558098994, 0.02283833222226), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.02283833222226, 0.69583608678780), 0.02235677320230),
      //    IntegrationPoint(Vector2d(0.11625191590760, 0.85801403354407), 0.01731623110866),
      //    IntegrationPoint(Vector2d(0.85801403354407, 0.02573405054833), 0.01731623110866),
      //    IntegrationPoint(Vector2d(0.02573405054833, 0.11625191590760), 0.01731623110866),
      //    IntegrationPoint(Vector2d(0.85801403354407, 0.11625191590760), 0.01731623110866),
      //    IntegrationPoint(Vector2d(0.11625191590760, 0.02573405054833), 0.01731623110866),
      //    IntegrationPoint(Vector2d(0.02573405054833, 0.85801403354407), 0.01731623110866)
      //  };
      default:
        KRATOS_THROW_ERROR(std::logic_error, "IntegrationUtilities::GetFace: face_id does not exist", "");
      }
    }

  protected:

  private:
    /// Default constructor.
    IntegrationUtilities() { }
    /// Assignment operator.
    IntegrationUtilities& operator=(IntegrationUtilities const& rOther);
  };
} // namespace Kratos.

#endif // INTEGRATION_UTILITIES
