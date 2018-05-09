//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas (gcasas@cimmne.upc.edu) $
//   Date:                $Date: 2016-5-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
#include "basset_force_tools.h"
namespace Kratos
{
void BassetForceTools::FillDaitcheVectors(const int N, const int order, const int n_steps_per_quad_step)
{
    std::cout << "\nFilling up vectors of coefficients for Daitche quadrature...\n";

    if (!N){
        return;
    }

    std::vector<double>& Ajs = SphericSwimmingParticle<SphericParticle>::mAjs;
    std::vector<double>& Bns = SphericSwimmingParticle<SphericParticle>::mBns;
    std::vector<double>& Cns = SphericSwimmingParticle<SphericParticle>::mCns;
    std::vector<double>& Dns = SphericSwimmingParticle<SphericParticle>::mDns;
    std::vector<double>& Ens = SphericSwimmingParticle<SphericParticle>::mEns;

    if (order == 1){
        long double FourThirds = 4.0 / 3;

        Ajs.resize(n_steps_per_quad_step * N);
        Bns.resize(n_steps_per_quad_step * N);

        for (int i = 0; i < n_steps_per_quad_step; ++i){
            long double alpha = static_cast<long double>(i + 1) / n_steps_per_quad_step;
            long double sqrt_alpha            = std::sqrt(alpha);
            long double sqrt_alpha_plus       = std::sqrt(alpha + 1);
            Ajs[i] = static_cast<double>(FourThirds * sqrt_alpha);
            Bns[i] = static_cast<double>(FourThirds * sqrt_alpha);

            Ajs[n_steps_per_quad_step + i] = static_cast<double>(FourThirds * ((1 + alpha) * (sqrt_alpha_plus - sqrt_alpha)));
            Bns[n_steps_per_quad_step + i] = static_cast<double>(FourThirds * ((1 + alpha) * (sqrt_alpha_plus - sqrt_alpha)));
        }

        for (int j = 2; j < N; ++j){
            for (int i = 0; i < n_steps_per_quad_step; ++i){
                long double alpha = static_cast<long double>(i + 1) / n_steps_per_quad_step;
                long double sqrt_j             = std::sqrt(j - 1 + alpha);
                long double sqrt_j_minus       = std::sqrt(j - 2 + alpha);
                long double sqrt_j_plus        = std::sqrt(j + alpha);
                long double sqrt_j_cubed       = SWIMMING_POW_3(sqrt_j);
                long double sqrt_j_minus_cubed = SWIMMING_POW_3(sqrt_j_minus);
                long double sqrt_j_plus_cubed  = SWIMMING_POW_3(sqrt_j_plus);

                Ajs[n_steps_per_quad_step * j + i] = static_cast<double>(FourThirds * (sqrt_j_minus_cubed + sqrt_j_plus_cubed - 2 * sqrt_j_cubed));
                Bns[n_steps_per_quad_step * j + i] = static_cast<double>(FourThirds * (sqrt_j_minus_cubed - sqrt_j_cubed + 1.5 * sqrt_j));
            }
        }
    }

    else if (order == 2){
        long double OneFifteenth = 1.0 / 15;
        Ajs.resize(n_steps_per_quad_step * N);
        Bns.resize(n_steps_per_quad_step * N);
        Cns.resize(n_steps_per_quad_step * N);

        for (int i = 0; i < n_steps_per_quad_step; ++i){
            long double alpha = static_cast<long double>(i + 1) / n_steps_per_quad_step;
            long double OneOverFifteenAlphas = OneFifteenth / alpha;
            long double sqrt_alpha_plus_1       = std::sqrt(alpha + 1);
            long double sqrt_alpha_plus_2       = std::sqrt(alpha + 2);
            long double sqrt_alpha_plus_3       = std::sqrt(alpha + 3);
            Ajs[i]                             = 4 * sqrt_alpha_plus_1 * (4 * alpha - 1) * OneOverFifteenAlphas;
            Ajs[n_steps_per_quad_step + i]     = 4 * SWIMMING_POW_5(sqrt_alpha_plus_1) * OneOverFifteenAlphas + 2 * (sqrt_alpha_plus_2 * (6 + 11.0 * alpha + 4 * alpha * alpha) - sqrt_alpha_plus_1 * (9 + 13 * alpha + 4 * alpha * alpha)) * OneFifteenth;
            Ajs[2 * n_steps_per_quad_step + i] = 2 * OneFifteenth * (2 * sqrt_alpha_plus_1 * (-2 +  7 * alpha + 4 * alpha * alpha)
                                                                   - 3 * sqrt_alpha_plus_2 * ( 6 + 11 * alpha + 4 * alpha * alpha)
                                                                       + sqrt_alpha_plus_3 * (21 + 19 * alpha + 4 * alpha * alpha));
        }

        for (int j = 3; j < N; ++j){
            for (int i = 0; i < n_steps_per_quad_step; ++i){
                long double alpha = static_cast<long double>(i + 1) / n_steps_per_quad_step;
                long double sqrt_j_minus_2 = std::sqrt(j + alpha - 3);
                long double sqrt_j_minus   = std::sqrt(j + alpha - 2);
                long double sqrt_j         = std::sqrt(j + alpha - 1);
                long double sqrt_j_plus    = std::sqrt(j + alpha);
                long double sqrt_j_plus_2  = std::sqrt(j + alpha + 1);

                Ajs[n_steps_per_quad_step * j + i] = static_cast<double>(
                          8 * OneFifteenth * (      SWIMMING_POW_5(sqrt_j_plus_2) - 3 * SWIMMING_POW_5(sqrt_j_plus) + 3 * SWIMMING_POW_5(sqrt_j) - SWIMMING_POW_5(sqrt_j_minus))
                       + 10 * OneFifteenth * (    - SWIMMING_POW_3(sqrt_j_plus_2) + 3 * SWIMMING_POW_3(sqrt_j_plus) - 3 * SWIMMING_POW_3(sqrt_j) + SWIMMING_POW_3(sqrt_j_minus)));
                Bns[n_steps_per_quad_step * j + i] = static_cast<double>(
                          8 * OneFifteenth * (- 2 * SWIMMING_POW_5(sqrt_j)        + 3 * SWIMMING_POW_5(sqrt_j_minus)    - SWIMMING_POW_5(sqrt_j_minus_2))
                       + 10 * OneFifteenth * (  4 * SWIMMING_POW_3(sqrt_j)        - 3 * SWIMMING_POW_3(sqrt_j_minus)    + SWIMMING_POW_3(sqrt_j_minus_2)));
                Cns[n_steps_per_quad_step * j + i] = static_cast<double>(
                          8 * OneFifteenth * (      SWIMMING_POW_5(sqrt_j)            - SWIMMING_POW_5(sqrt_j_minus))
                       + 10 * OneFifteenth * (- 3 * SWIMMING_POW_3(sqrt_j)            + SWIMMING_POW_3(sqrt_j_minus)) + 2 * sqrt_j);

            }
        }
    }

    else {
        long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
        long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
        long double sqrt_5_over_63  = std::sqrt(static_cast<long double>(5)) / 63;
        long double OneOver315 = 1.0 / 315;
        Ajs.resize(N);
        Bns.resize(N);
        Cns.resize(N);
        Dns.resize(N);
        Ens.resize(N);
        Ajs[0] = 244 * sqrt_2_over_315;
        Ajs[1] = 362 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
        Ajs[2] = 5584 * OneOver315 - 1448 * sqrt_3_over_105 + 1464 * sqrt_2_over_315;
        Ajs[3] = 1130 * sqrt_5_over_63 - 22336 * OneOver315 + 2172 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
        long double OneOver9 = 1.0 / 9;

        for (int j = 3; j < N; ++j){
            long double sqrt_j_minus_5 = std::sqrt(static_cast<long double>(j - 5));
            long double sqrt_j_minus_4 = std::sqrt(static_cast<long double>(j - 4));
            long double sqrt_j_minus_3 = std::sqrt(static_cast<long double>(j - 3));
            long double sqrt_j_minus_2 = std::sqrt(static_cast<long double>(j - 2));
            long double sqrt_j_minus   = std::sqrt(static_cast<long double>(j - 1));
            long double sqrt_j         = std::sqrt(static_cast<long double>(j));
            long double sqrt_j_plus    = std::sqrt(static_cast<long double>(j + 1));
            long double sqrt_j_plus_2  = std::sqrt(static_cast<long double>(j + 2));
            Ajs[j] = static_cast<double>(
                     48 * OneOver315 * (       SWIMMING_POW_7(sqrt_j_plus_2)      + SWIMMING_POW_7(sqrt_j_minus_2)  - 4 * SWIMMING_POW_7(sqrt_j_plus)    - 4 * SWIMMING_POW_7(sqrt_j_minus)   + 6 * SWIMMING_POW_7(sqrt_j))
                    + 2 * OneOver9   * (   4 * SWIMMING_POW_3(sqrt_j_plus)    + 4 * SWIMMING_POW_3(sqrt_j_minus)        - SWIMMING_POW_3(sqrt_j_plus_2)      - SWIMMING_POW_3(sqrt_j_minus_2) - 6 * SWIMMING_POW_3(sqrt_j)));
            Bns[j] = static_cast<double>(
                     48 * OneOver315 * (       SWIMMING_POW_7(sqrt_j)         - 4 * SWIMMING_POW_7(sqrt_j_minus_2)  + 6 * SWIMMING_POW_7(sqrt_j_minus_3) - 4 * SWIMMING_POW_7(sqrt_j_minus_4)     + SWIMMING_POW_7(sqrt_j_minus_5)) - 168 * OneOver315 * SWIMMING_POW_5(sqrt_j)
                        + OneOver9   * (   4 * SWIMMING_POW_3(sqrt_j)         + 8 * SWIMMING_POW_3(sqrt_j_minus_2) - 12 * SWIMMING_POW_3(sqrt_j_minus_3) + 8 * SWIMMING_POW_3(sqrt_j_minus_4) - 2 * SWIMMING_POW_3(sqrt_j_minus_5)));
            Cns[j] = static_cast<double>(
                    48 * OneOver315 *  (       SWIMMING_POW_7(sqrt_j_minus_4) - 4 * SWIMMING_POW_7(sqrt_j_minus_3)  + 6 * SWIMMING_POW_7(sqrt_j_minus_2) - 3 * SWIMMING_POW_7(sqrt_j)) + 672 * OneOver315 * SWIMMING_POW_5(sqrt_j)
                        + OneOver9  *  (- 18 * SWIMMING_POW_3(sqrt_j)        - 12 * SWIMMING_POW_3(sqrt_j_minus_2)  + 8 * SWIMMING_POW_3(sqrt_j_minus_3) - 2 * SWIMMING_POW_3(sqrt_j_minus_4)));
            Dns[j] = static_cast<double>(
                    48 * OneOver315 *  (   3 * SWIMMING_POW_7(sqrt_j)         - 4 * SWIMMING_POW_7(sqrt_j_minus_2)      + SWIMMING_POW_7(sqrt_j_minus_3)) - 24 * OneOver9 * SWIMMING_POW_5(sqrt_j)
                         + OneOver9 *  (  36 * SWIMMING_POW_3(sqrt_j)         + 8 * SWIMMING_POW_3(sqrt_j_minus_2)  - 2 * SWIMMING_POW_3(sqrt_j_minus_3)));
            Ens[j] = static_cast<double>(
                    48 * OneOver315 *  (       SWIMMING_POW_7(sqrt_j_minus_2)     - SWIMMING_POW_7(sqrt_j)) + 336 * OneOver315 * SWIMMING_POW_5(sqrt_j)
                         + OneOver9 *  (- 22 * SWIMMING_POW_3(sqrt_j)         - 2 * SWIMMING_POW_3(sqrt_j_minus_2)) + 2 * sqrt_j);
        }

    }

    // Filling up the static variables of other versions of swimming particle, since each template instance has their own variables

    std::vector<double>& Ajs_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mAjs;
    std::vector<double>& Bns_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mBns;
    std::vector<double>& Cns_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mCns;
    std::vector<double>& Dns_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mDns;
    std::vector<double>& Ens_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mEns;

    Ajs_analytic = Ajs;
    Bns_analytic = Bns;
    Cns_analytic = Cns;
    Dns_analytic = Dns;
    Ens_analytic = Ens;

    std::cout << "...Finished filling up vectors of coefficients.\n";
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This method precalculates values needed for the implementation of the method described by van Hinsberg et al. (2011)
// The method requires the determination of m numbers ai and m other numbers ti which determine m exponentials.The sum of these exponentials is unsed as an approximation to the tail of the Basset kernel
// In the following particular values of these numbers are presented up to m = 10. These numbers were obtained by Casas and Ferrer (2016)
void BassetForceTools::FillHinsbergVectors(ModelPart& r_model_part, const int m, const int n_quad_delta_times_window)
{
    KRATOS_TRY
    mNumberOfQuadratureStepsInWindow = n_quad_delta_times_window;
    mNumberOfExponentials = m;

    if (!m){
        return;
    }

    std::cout << "\nFilling up vectors of coefficients for Hinsberg method with m = " << m << " ...\n";
    double & t_win = SphericSwimmingParticle<SphericParticle>::mTimeWindow;
    std::vector<double>& As = SphericSwimmingParticle<SphericParticle>::mAs;
    std::vector<double>& Ts = SphericSwimmingParticle<SphericParticle>::mTs;
    std::vector<double>& Alphas = SphericSwimmingParticle<SphericParticle>::mAlphas;
    std::vector<double>& Betas = SphericSwimmingParticle<SphericParticle>::mBetas;

    mTimeWindow = mNumberOfQuadratureStepsInWindow * r_model_part.GetProcessInfo()[TIME_STEPS_PER_QUADRATURE_STEP] * r_model_part.GetProcessInfo()[DELTA_TIME];
    t_win = mTimeWindow;

    As.resize(m);
    Ts.resize(m);
    Alphas.resize(m);
    Betas.resize(m);
    std::string points_optimization_method = "AbsValue"; // toggle between types of optimization

    if (points_optimization_method == "AbsValue"){
        if (m == 1){
            As[0] = 1.046347992;
            Ts[0] = 1.581186674;
        }

        else if (m == 2){
            As[0] = 0.566192817;
            Ts[0] = 0.717656182;
            As[1] = 0.864298391;
            Ts[1] = 8.925153279;
        }

        else if (m == 3){
            As[0] = 0.440072204;
            Ts[0] = 0.482318894;
            As[1] = 0.538287204;
            Ts[1] = 3.324763126;
            As[2] = 0.807797346;
            Ts[2] = 38.928376132;
        }

        else if (m == 4){
            As[0] = 0.374397988;
            Ts[0] = 0.365083559;
            As[1] = 0.421322343;
            Ts[1] = 1.820334739;
            As[2] = 0.517872275;
            Ts[2] = 11.809488351;
            As[3] = 0.761539469;
            Ts[3] = 127.109159354;
        }

        else if (m == 5){
            As[0] = 0.3450551877;
            Ts[0] = 0.3227320427;
            As[1] = 0.3762685526;
            Ts[1] = 1.4017593843;
            As[2] = 0.4383511621;
            Ts[2] = 7.3543952717;
            As[3] = 0.5502868981;
            Ts[3] = 52.9058339347;
            As[4] = 0.7701813938;
            Ts[4] = 699.4337431732;
        }

        else if (m == 6){
            As[0] = 0.3227460255;
            Ts[0] = 0.2894856389;
            As[1] = 0.3446901326;
            Ts[1] = 1.1312690586;
            As[2] = 0.3924441164;
            Ts[2] = 5.1207861657;
            As[3] = 0.471576099;
            Ts[3] = 29.6345412934;
            As[4] = 0.5990063177;
            Ts[4] = 256.64908268;
            As[5] = 0.7695849793;
            Ts[5] = 4254.1241751139;
        }

        else if (m == 7){
            As[0] = 0.2931405176;
            Ts[0] = 0.2413624327;
            As[1] = 0.3053190176;
            Ts[1] = 0.8199848671;
            As[2] = 0.3394616674;
            Ts[2] = 3.0838532791;
            As[3] = 0.3924532926;
            Ts[3] = 13.8047974118;
            As[4] = 0.4794140412;
            Ts[4] = 80.9779742728;
            As[5] = 0.5546383969;
            Ts[5] = 696.8320792921;
            As[6] = 0.6207864425;
            Ts[6] = 6133.2449027098;
        }

        else if (m == 8){
            As[0] = 0.2718360249;
            Ts[0] = 0.2192620346;
            As[1] = 0.2685924185;
            Ts[1] = 0.662026818;
            As[2] = 0.2871214552;
            Ts[2] = 2.0706383247;
            As[3] = 0.3249589764;
            Ts[3] = 7.2825402363;
            As[4] = 0.3805886345;
            Ts[4] = 31.0062809826;
            As[5] = 0.4469592071;
            Ts[5] = 169.6857783353;
            As[6] = 0.5474439544;
            Ts[6] = 1226.001409491;
            As[7] = 0.7637048975;
            Ts[7] = 17271.9375778519;
        }

        else if (m == 9){
            As[0] = 0.2570818336;
            Ts[0] = 0.1878604572;
            As[1] = 0.2610118588;
            Ts[1] = 0.5420260992;
            As[2] = 0.2799238451;
            Ts[2] = 1.6534881587;
            As[3] = 0.3051985477;
            Ts[3] = 5.5204876302;
            As[4] = 0.3418149337;
            Ts[4] = 20.8847203692;
            As[5] = 0.3892337642;
            Ts[5] = 93.9005719593;
            As[6] = 0.4655655296;
            Ts[6] = 532.1532341216;
            As[7] = 0.6107696402;
            Ts[7] = 4683.3937018005;
            As[8] = 0.784623916;
            Ts[8] = 93277.7129340798;
        }

        else if (m == 10){
            As[0] = 0.2520642358;
            Ts[0] = 0.1878604572;
            As[1] = 0.254913066;
            Ts[1] = 0.5306382498;
            As[2] = 0.2638832071;
            Ts[2] = 1.5524873935;
            As[3] = 0.2666445191;
            Ts[3] = 4.6517443725;
            As[4] = 0.2806268115;
            Ts[4] = 14.2413555446;
            As[5] = 0.344914608;
            Ts[5] = 50.7413819742;
            As[6] = 0.4566204962;
            Ts[6] = 263.7561507819;
            As[7] = 0.5663046247;
            Ts[7] = 2146.211201895;
            As[8] = 0.6253574036;
            Ts[8] = 26744.590748687;
            As[9] = 0.6932526975;
            Ts[9] = 348322.670028861;
        }
    } // if (points_optimization_method == 'AbsValue')

    else if (points_optimization_method == "TNorm"){
        if (m == 1){
            As[0] = 1.046347992;
            Ts[0] = 1.581186674;
        }

        else if (m == 2){
            As[0] = 0.547059755163;
            Ts[0] = 0.666683527541;
            As[1] = 0.844976749056;
            Ts[1] = 8.342487240727;
        }

        else if (m == 3){
            As[0] = 0.430797004966;
            Ts[0] = 0.452146141420;
            As[1] = 0.531940201577;
            Ts[1] = 3.059709731091;
            As[2] = 0.804647149264;
            Ts[2] = 36.769402334968;
        }

        else if (m == 4){
            As[0] = 0.371405161268;
            Ts[0] = 0.350505616164;
            As[1] = 0.422130638640;
            Ts[1] = 1.752574133530;
            As[2] = 0.524882763770;
            Ts[2] = 11.652875613838;
            As[3] = 0.781431790230;
            Ts[3] = 136.886412459808;
        }

        else if (m == 5){
            As[0] = 0.333573629082;
            Ts[0] = 0.290461028927;
            As[1] = 0.362933117265;
            Ts[1] = 1.203691574020;
            As[2] = 0.419725251898;
            Ts[2] = 5.937032480632;
            As[3] = 0.520201697957;
            Ts[3] = 39.145059811464;
            As[4] = 0.766103870248;
            Ts[4] = 452.822622886880;
        }

        else if (m == 6){
            As[0] = 0.306592856340;
            Ts[0] = 0.250430971277;
            As[1] = 0.324348018693;
            Ts[1] = 0.910305675752;
            As[2] = 0.361593254514;
            Ts[2] = 3.720497699441;
            As[3] = 0.418122688955;
            Ts[3] = 18.272742276055;
            As[4] = 0.516808573539;
            Ts[4] = 119.760302387015;
            As[5] = 0.755114941277;
            Ts[5] = 1369.901637784433;
        }

        else if (m == 7){
            As[0] = 0.286966758381;
            Ts[0] = 0.222956735521;
            As[1] = 0.297777435995;
            Ts[1] = 0.737995019307;
            As[2] = 0.324921880389;
            Ts[2] = 2.658309910255;
            As[3] = 0.363168742275;
            Ts[3] = 10.923732152073;
            As[4] = 0.420482447032;
            Ts[4] = 54.149026920965;
            As[5] = 0.520771163430;
            Ts[5] = 360.637576912179;
            As[6] = 0.755431859453;
            Ts[6] = 4254.124341110529;
        }

        else if (m == 8){
            As[0] = 0.269592611484;
            Ts[0] = 0.199808472359;
            As[1] = 0.275162851282;
            Ts[1] = 0.609421758879;
            As[2] = 0.295415500694;
            Ts[2] = 1.974629286510;
            As[3] = 0.322815911082;
            Ts[3] = 7.052730788656;
            As[4] = 0.360270264199;
            Ts[4] = 28.694274517306;
            As[5] = 0.415929367307;
            Ts[5] = 140.189096127922;
            As[6] = 0.512156883919;
            Ts[6] = 911.255504581087;
            As[7] = 0.740228044588;
            Ts[7] = 10263.341976325066;
        }

        else if (m == 9){
            As[0] = 0.256076630317;
            Ts[0] = 0.182624446580;
            As[1] = 0.258081235942;
            Ts[1] = 0.523116633552;
            As[2] = 0.273953523605;
            Ts[2] = 1.566580998223;
            As[3] = 0.295097701391;
            Ts[3] = 5.063946393605;
            As[4] = 0.322494108162;
            Ts[4] = 18.066436391438;
            As[5] = 0.359800501707;
            Ts[5] = 73.405444944814;
            As[6] = 0.415133110865;
            Ts[6] = 357.949475288163;
            As[7] = 0.510476026521;
            Ts[7] = 2319.768464890347;
            As[8] = 0.734899701220;
            Ts[8] = 25980.611692219158;
        }

        else if (m == 10){
            As[0] = 0.246702083145;
            Ts[0] = 0.171137410203;
            As[1] = 0.246474944419;
            Ts[1] = 0.469538455725;
            As[2] = 0.259717868181;
            Ts[2] = 1.333604723390;
            As[3] = 0.277340588232;
            Ts[3] = 4.038729849045;
            As[4] = 0.299501001874;
            Ts[4] = 13.268683433850;
            As[5] = 0.328282204667;
            Ts[5] = 48.350555319658;
            As[6] = 0.367882181136;
            Ts[6] = 202.201304412825;
            As[7] = 0.427624033730;
            Ts[7] = 1029.089927961859;
            As[8] = 0.533580013899;
            Ts[8] = 7177.875290938712;
            As[9] = 0.765266538864;
            Ts[9] = 93277.737373373078;
        }
    }

    else {
        KRATOS_THROW_ERROR(std::invalid_argument, "van Hinsberg's method is only implemented up to a number of exponentials m = 10.", m);
    }

    const double e = std::exp(1);

    for (int i = 0; i < m; i++){
        Ts[i] *= mTimeWindow;  // here we recover the dimensional tis from the non-dimensional tis (which are idependent of the time_window)
        Alphas[i] = std::sqrt(e / Ts[i]);
        Betas[i] = - 0.5 / Ts[i];
    }

    // Filling up the static variables of other versions of swimming particle, since each template instance has their own variables

     double & t_win_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mTimeWindow;
     std::vector<double>& As_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mAs;
     std::vector<double>& Ts_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mTs;
     std::vector<double>& Alphas_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mAlphas;
     std::vector<double>& Betas_analytic = SphericSwimmingParticle<AnalyticSphericParticle>::mBetas;

     t_win_analytic = t_win;
     As_analytic = As;
     Ts_analytic = Ts;
     Alphas_analytic = Alphas;
     Betas_analytic = Betas;

    // Filling up the particles' individual vectors

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        DenseVector<double>& hinsberg_tail_contributions = inode->GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
        hinsberg_tail_contributions.resize(3 * m + 3); //  the extra contribution will contain the oldest integrand that has been discarded from the integrands vector
        for (int i = 0; i < 3 * m + 3; i++){
            hinsberg_tail_contributions[i] = 0.0;
        }
    }

    std::cout << "...Finished filling up vectors of coefficients.\n";
    KRATOS_CATCH("")
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void BassetForceTools::AppendIntegrands(ModelPart& r_model_part)
{
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    double time = r_process_info[TIME] + r_process_info[DELTA_TIME];
    r_process_info[LAST_TIME_APPENDING] = time;

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        DenseVector<double>& historic_integrands             = inode->GetValue(BASSET_HISTORIC_INTEGRANDS);
        const array_1d<double, 3>& fluid_vel_projected  = inode->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel         = inode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel                    = fluid_vel_projected - particle_vel;

        int n = historic_integrands.size();

        historic_integrands.resize(n + 3);
        historic_integrands.insert_element(n,     slip_vel[0]);
        historic_integrands.insert_element(n + 1, slip_vel[1]);
        historic_integrands.insert_element(n + 2, slip_vel[2]);
    }

}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void BassetForceTools::AppendIntegrandsImplicit(ModelPart& r_model_part)
{
    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    process_info[LAST_TIME_APPENDING] = r_model_part.GetProcessInfo()[TIME];

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        DenseVector<double>& historic_integrands             = inode->GetValue(BASSET_HISTORIC_INTEGRANDS);
        const array_1d<double, 3>& fluid_vel_projected  = inode->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel         = inode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel                    = fluid_vel_projected - particle_vel;
        int n = historic_integrands.size();

        if (mFirstTimeAppending){
            mFirstTimeAppending = false;
            historic_integrands.resize(n + 9);
            historic_integrands.insert_element(n    , slip_vel[0]);
            historic_integrands.insert_element(n + 1, slip_vel[1]);
            historic_integrands.insert_element(n + 2, slip_vel[2]);
            historic_integrands.insert_element(n + 3, particle_vel[0]);
            historic_integrands.insert_element(n + 4, particle_vel[1]);
            historic_integrands.insert_element(n + 5, particle_vel[2]);
            historic_integrands.insert_element(n + 6, particle_vel[0]);
            historic_integrands.insert_element(n + 7, particle_vel[1]);
            historic_integrands.insert_element(n + 8, particle_vel[2]);
        }
        else {
            array_1d<double, 3> old_particle_vel;
            // we are using the last positions to store the old particle velocities which are not stored anywhere else
            old_particle_vel[0] = historic_integrands[n - 3];
            old_particle_vel[1] = historic_integrands[n - 2];
            old_particle_vel[2] = historic_integrands[n - 1];
            historic_integrands.resize(n + 3);
            historic_integrands.insert_element(n - 6, slip_vel[0]);
            historic_integrands.insert_element(n - 5, slip_vel[1]);
            historic_integrands.insert_element(n - 4, slip_vel[2]);
            historic_integrands.insert_element(n - 3, old_particle_vel[0]);
            historic_integrands.insert_element(n - 2, old_particle_vel[1]);
            historic_integrands.insert_element(n - 1, old_particle_vel[2]);
            historic_integrands.insert_element(n,     particle_vel[0]);
            historic_integrands.insert_element(n + 1, particle_vel[1]);
            historic_integrands.insert_element(n + 2, particle_vel[2]);
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void BassetForceTools::AppendIntegrandsWindow(ModelPart& r_model_part)
{
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    double time = r_process_info[TIME] + r_process_info[DELTA_TIME];
    r_process_info[LAST_TIME_APPENDING] = time;

    if (r_process_info[BASSET_FORCE_TYPE] >= 3 && mNumberOfExponentials > 0){

        for (ElementIterator iparticle = r_model_part.ElementsBegin(); iparticle != r_model_part.ElementsEnd(); iparticle++){
            Node<3>& node = iparticle->GetGeometry()[0];
            if (node.IsNot(BLOCKED)){
                DenseVector<double>& historic_integrands = node.GetValue(BASSET_HISTORIC_INTEGRANDS);

                if (int(historic_integrands.size()) >= 3 * mNumberOfQuadratureStepsInWindow){
                    DenseVector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
                    hinsberg_tail_contributions.resize(3 * mNumberOfExponentials + 3); // in case there is an inlet and new particles with empty vectors come about
                    hinsberg_tail_contributions[3 * mNumberOfExponentials]     = historic_integrands[0];
                    hinsberg_tail_contributions[3 * mNumberOfExponentials + 1] = historic_integrands[1];
                    hinsberg_tail_contributions[3 * mNumberOfExponentials + 2] = historic_integrands[2];
                }
            }
        }
    }

    for (ElementIterator iparticle = r_model_part.ElementsBegin(); iparticle != r_model_part.ElementsEnd(); iparticle++){
        Node<3>& node = iparticle->GetGeometry()[0];
        if (node.IsNot(BLOCKED)){
            DenseVector<double>& historic_integrands             = node.GetValue(BASSET_HISTORIC_INTEGRANDS);
            const array_1d<double, 3>& fluid_vel_projected  = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3>& particle_vel         = node.FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3> slip_vel                    = fluid_vel_projected - particle_vel;
            int n = historic_integrands.size();

            if (n < 3 * mNumberOfQuadratureStepsInWindow){ // list of integrands still growing
                historic_integrands.resize(n + 3);
                historic_integrands.insert_element(n,     slip_vel[0]);
                historic_integrands.insert_element(n + 1, slip_vel[1]);
                historic_integrands.insert_element(n + 2, slip_vel[2]);
            }

            else { // forget oldest integrand
                for (int i = 0; i < n / 3 - 1; i++){
                    historic_integrands[3 * i]     = historic_integrands[3 * (i + 1)];
                    historic_integrands[3 * i + 1] = historic_integrands[3 * (i + 1) + 1];
                    historic_integrands[3 * i + 2] = historic_integrands[3 * (i + 1) + 2];
                }

                historic_integrands.insert_element(n - 3, slip_vel[0]);
                historic_integrands.insert_element(n - 2, slip_vel[1]);
                historic_integrands.insert_element(n - 1, slip_vel[2]);
            }
        }
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

} // namespace Kratos
