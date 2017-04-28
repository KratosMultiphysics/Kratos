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
    static int getGaussPointLocations(unsigned int number_of_points,
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

  protected:

  private:
    /// Default constructor.
    IntegrationUtilities() { }
    /// Assignment operator.
    IntegrationUtilities& operator=(IntegrationUtilities const& rOther);
  };
} // namespace Kratos.

#endif // INTEGRATION_UTILITIES
