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
void BassetForceTools::FillDaitcheVectors(const int N, const int order)
{
    std::cout << "\nFilling up vectors of coefficients for Daitche quadrature...";

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

        Ajs.resize(N);
        Bns.resize(N);
        Ajs[0] = FourThirds;
        Bns[0] = FourThirds;

        for (int j = 1; j < N; ++j){
            long double sqrt_j       = std::sqrt(static_cast<long double>(j));
            long double sqrt_j_minus = std::sqrt(static_cast<long double>(j - 1));
            long double sqrt_j_plus  = std::sqrt(static_cast<long double>(j + 1));
            long double sqrt_j_cubed       = sqrt_j * sqrt_j * sqrt_j;
            long double sqrt_j_minus_cubed = sqrt_j_minus * sqrt_j_minus * sqrt_j_minus;
            long double sqrt_j_plus_cubed  = sqrt_j_plus * sqrt_j_plus * sqrt_j_plus;

            Ajs[j] = static_cast<double>(FourThirds * (sqrt_j_minus_cubed + sqrt_j_plus_cubed - 2 * sqrt_j_cubed));
            Bns[j] = static_cast<double>(FourThirds * (sqrt_j_minus_cubed - sqrt_j_cubed + 1.5 * sqrt_j));
        }
    }

    else if (order == 2){
        long double OneFifteenth = 1.0 / 15;
        long double sqrt_2_over_5 = std::sqrt(static_cast<long double>(2)) / 5;
        long double sqrt_3_over_5 = std::sqrt(static_cast<long double>(3)) / 5;
        Ajs.resize(N);
        Bns.resize(N);
        Cns.resize(N);
        Ajs[0] = 4 * sqrt_2_over_5;
        Ajs[1] = 14 * sqrt_3_over_5 - 12 * sqrt_2_over_5;
        Ajs[2] = 176 * OneFifteenth - 42 * sqrt_3_over_5 + 12 * sqrt_2_over_5;
//        Bns[1] = 48 * sqrt_2_over_5;
//        Bns[2] = - 8 * sqrt_3_over_5 + 12 * sqrt_2_over_5;
//        Cns[2] = 10 * OneFifteenth * sqrt_2_over_5;
//        Cns[3] = 4 * sqrt_3_over_5 - 4 * sqrt_2_over_5;

        for (int j = 3; j < N; ++j){
            long double sqrt_j_minus_2 = std::sqrt(static_cast<long double>(j - 2));
            long double sqrt_j_minus   = std::sqrt(static_cast<long double>(j - 1));
            long double sqrt_j         = std::sqrt(static_cast<long double>(j));
            long double sqrt_j_plus    = std::sqrt(static_cast<long double>(j + 1));
            long double sqrt_j_plus_2  = std::sqrt(static_cast<long double>(j + 2));

            Ajs[j] = static_cast<double>(
                      8 * OneFifteenth * (      SWIMMING_POW_5(sqrt_j_plus_2) - 3 * SWIMMING_POW_5(sqrt_j_plus) + 3 * SWIMMING_POW_5(sqrt_j) - SWIMMING_POW_5(sqrt_j_minus))
                   + 10 * OneFifteenth * (    - SWIMMING_POW_3(sqrt_j_plus_2) + 3 * SWIMMING_POW_3(sqrt_j_plus) - 3 * SWIMMING_POW_3(sqrt_j) + SWIMMING_POW_3(sqrt_j_minus)));
            Bns[j] = static_cast<double>(
                      8 * OneFifteenth * (- 2 * SWIMMING_POW_5(sqrt_j)        + 3 * SWIMMING_POW_5(sqrt_j_minus)    - SWIMMING_POW_5(sqrt_j_minus_2))
                   + 10 * OneFifteenth * (  4 * SWIMMING_POW_3(sqrt_j)        - 3 * SWIMMING_POW_3(sqrt_j_minus)    + SWIMMING_POW_3(sqrt_j_minus_2)));
            Cns[j] = static_cast<double>(
                      8 * OneFifteenth * (      SWIMMING_POW_5(sqrt_j)            - SWIMMING_POW_5(sqrt_j_minus))
                   + 10 * OneFifteenth * (- 3 * SWIMMING_POW_3(sqrt_j)            + SWIMMING_POW_3(sqrt_j_minus)) + 2 * sqrt_j);
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
                    48 * OneOver315 *  (       SWIMMING_POW_7(sqrt_j_minus_4)  - 4 * SWIMMING_POW_7(sqrt_j_minus_3)  + 6 * SWIMMING_POW_7(sqrt_j_minus_2) - 3 * SWIMMING_POW_7(sqrt_j)) + 672 * OneOver315 * SWIMMING_POW_5(sqrt_j)
                        + OneOver9  *  (- 18 * SWIMMING_POW_3(sqrt_j)        - 12 * SWIMMING_POW_3(sqrt_j_minus_2)  + 8 * SWIMMING_POW_3(sqrt_j_minus_3) - 2 * SWIMMING_POW_3(sqrt_j_minus_4)));
            Dns[j] = static_cast<double>(
                    48 * OneOver315 *  (   3 * SWIMMING_POW_7(sqrt_j)          - 4 * SWIMMING_POW_7(sqrt_j_minus_2)      + SWIMMING_POW_7(sqrt_j_minus_3)) - 24 * OneOver9 * SWIMMING_POW_5(sqrt_j)
                         + OneOver9 *  (  36 * SWIMMING_POW_3(sqrt_j)          + 8 * SWIMMING_POW_3(sqrt_j_minus_2)  - 2 * SWIMMING_POW_3(sqrt_j_minus_3)));
            Ens[j] = static_cast<double>(
                    48 * OneOver315 *  (       SWIMMING_POW_7(sqrt_j_minus_2)      - SWIMMING_POW_7(sqrt_j)) + 336 * OneOver315 * SWIMMING_POW_5(sqrt_j)
                         + OneOver9 *  (- 22 * SWIMMING_POW_3(sqrt_j)          - 2 * SWIMMING_POW_3(sqrt_j_minus_2)) + 2 * sqrt_j);
        }
    }

    std::cout << "...Finished filling up vectors of coefficients.\n";
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// This method precalculates values needed for the implementation of the method described by Von Hinsberg (2011)
// The method requires the determination of m numbers ai and another numbers ti which deteermine m exponentials.
// In the following particular values of these numbers are presented u to m = 8. These numbers were obtained by Casas and Ferrer (2016)
void BassetForceTools::FillHinsbergVectors(ModelPart& r_model_part, const int m, const double time_window)
{
    if (!m){
        return;
    }

    std::cout << "\nFilling up vectors of coefficients for Hinsberg method with m = " << m << " ...";
    double & t_win = SphericSwimmingParticle<SphericParticle>::mTimeWindow;
    std::vector<double>& As = SphericSwimmingParticle<SphericParticle>::mAs;
    std::vector<double>& Ts = SphericSwimmingParticle<SphericParticle>::mTs;
    std::vector<double>& Alphas = SphericSwimmingParticle<SphericParticle>::mAlphas;
    std::vector<double>& Betas = SphericSwimmingParticle<SphericParticle>::mBetas;

    t_win = time_window;
    As.resize(m);
    Ts.resize(m);
    Alphas.resize(m);
    Betas.resize(m);

    // Filling up common vectors

    if (m == 1){
        As[0] = 1.055699152;
        Ts[0] = 1.656571537;
    }

    else if (m == 2){
        As[0] = 0.595936548;
        Ts[0] = 0.758737731;
        As[1] = 0.765862627;
        Ts[1] = 8.130844515;
    }

    else if (m == 3){
        As[0] = 0.457076294;
        Ts[0] = 0.523735503;
        As[1] = 0.52049493;
        Ts[1] = 3.463557465;
        As[2] = 0.730234918;
        Ts[2] = 35.209010652;
    }

    else if (m == 4){
        As[0] = 0.378123792;
        Ts[0] = 0.377168054;
        As[1] = 0.420250984;
        Ts[1] = 1.883663548;
        As[2] = 0.515234662;
        Ts[2] = 12.085534613;
        As[3] = 0.747882647;
        Ts[3] = 127.522988007;
    }

    else if (m == 5){
        As[0] = 0.377585589;
        Ts[0] = 0.361079805;
        As[1] = 0.389837358;
        Ts[1] = 1.758926107;
        As[2] = 0.414949491;
        Ts[2] = 8.640539541;
        As[3] = 0.503856364;
        Ts[3] = 57.10122954;
        As[4] = 0.607332741;
        Ts[4] = 448.083463993;
    }

    else if (m == 6){
        As[0] = 0.338300743;
        Ts[0] = 0.346126312;
        As[1] = 0.345524197;
        Ts[1] = 1.386290002;
        As[2] = 0.368960284;
        Ts[2] = 5.934710427;
        As[3] = 0.368902685;
        Ts[3] = 27.706980453;
        As[4] = 0.432603065;
        Ts[4] = 132.567423265;
        As[5] = 0.771632072;
        Ts[5] = 1371.238854372;
    }

    else if (m == 7){
        As[0] = 0.28596607;
        Ts[0] = 0.222165932;
        As[1] = 0.296473585;
        Ts[1] = 0.730394698;
        As[2] = 0.323588913;
        Ts[2] = 2.617417995;
        As[3] = 0.360741831;
        Ts[3] = 10.658764953;
        As[4] = 0.417056856;
        Ts[4] = 52.200869695;
        As[5] = 0.514260513;
        Ts[5] = 340.581473772;
        As[6] = 0.746597256;
        Ts[6] = 3862.218173227;
    }

    else if (m == 8){
        As[0] = 0.2696042712;
        Ts[0] = 0.1998261277;
        As[1] = 0.275174292;
        Ts[1] = 0.6095050731;
        As[2] = 0.2954251215;
        Ts[2] = 1.9749864724;
        As[3] = 0.3228225133;
        Ts[3] = 7.0542458787;
        As[4] = 0.3602769943;
        Ts[4] = 28.7011577713;
        As[5] = 0.4159412411;
        Ts[5] = 140.2284577358;
        As[6] = 0.5121587415;
        Ts[6] = 911.5441865997;
        As[7] = 0.7402444628;
        Ts[7] = 10266.8313490027;
    }

    else {
        KRATOS_THROW_ERROR(std::invalid_argument, "Hinsberg's method is only implemented up to a number of exponentials m = 8.", m);
    }

    const double e = std::exp(1);

    for (int i = 0; i < m; i++){
        Alphas[i] = std::sqrt(e / Ts[i]);
        Betas[i] = - 0.5 / Ts[i];
    }

    // Filling up the particles' individual vectors

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        vector<double>& hinsberg_tail_contributions = inode->GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
        hinsberg_tail_contributions.resize(3 * m);
        for (int i = 0; i < 3 * m; i++){
            hinsberg_tail_contributions[i] = 0.0;
        }
    }

    std::cout << "...Finished filling up vectors of coefficients.\n";
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
double BassetForceTools::Phi(const double x)
{
    if (fabs(x) < 1e-10){
        return (std::exp(x) - 1) / x;
    }
    else {
        return 1 + 0.5 * x + 1.0 / 6 * x * x;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void BassetForceTools::AddFdi(const int order, array_1d<double, 3>& F, const double t_win, const double ti, const double beta, const double dt, const vector<double>& historic_integrands)
{
    double normalized_dt = 0.5 * dt / ti;

    if (order == 2){
        const double coeff = 2 * std::sqrt(ti) * std::exp(beta * t_win + 0.5);
        F[0] +=  coeff * (historic_integrands[0] * (1 - Phi(- normalized_dt)) + historic_integrands[3] * std::exp(- normalized_dt) * (Phi(normalized_dt) - 1));
        F[1] +=  coeff * (historic_integrands[1] * (1 - Phi(- normalized_dt)) + historic_integrands[4] * std::exp(- normalized_dt) * (Phi(normalized_dt) - 1));
        F[2] +=  coeff * (historic_integrands[2] * (1 - Phi(- normalized_dt)) + historic_integrands[5] * std::exp(- normalized_dt) * (Phi(normalized_dt) - 1));
    }

    else if (order == 1){
        const double coeff = 0.5 * std::sqrt(1.0 / ti) / (SWIMMING_POW_3(beta) * SWIMMING_POW_2(dt));
        const double exp_1 = exp((t_win + dt) * beta + 0.5);
        const double exp_2 = exp(t_win * beta + 0.5);
        const double f00 = historic_integrands[0];
        const double f01 = historic_integrands[1];
        const double f02 = historic_integrands[2];
        const double f10 = historic_integrands[3];
        const double f11 = historic_integrands[4];
        const double f12 = historic_integrands[5];
        const double f20 = historic_integrands[6];
        const double f21 = historic_integrands[7];
        const double f22 = historic_integrands[8];
        F[0] += coeff * (exp_1 * (2 * (2 * f10 + f20) + f00 * (2 + dt * beta) - dt * beta * (f20 + 2 * dt * f10 * beta))\
                       - exp_2 * (4 * f10 + 2 * f20 + dt * beta * (4 * f10 + f20) + f00 * (2 + dt * beta * (3 + 2 * dt * beta))));
        F[1] += coeff * (exp_1 * (2 * (2 * f11 + f21) + f01 * (2 + dt * beta) - dt * beta * (f21 + 2 * dt * f11 * beta))\
                       - exp_2 * (4 * f11 + 2 * f21 + dt * beta * (4 * f11 + f21) + f01 * (2 + dt * beta * (3 + 2 * dt * beta))));
        F[2] += coeff * (exp_1 * (2 * (2 * f12 + f22) + f02 * (2 + dt * beta) - dt * beta * (f22 + 2 * dt * f12 * beta))\
                       - exp_2 * (4 * f12 + 2 * f22 + dt * beta * (4 * f12 + f22) + f02 * (2 + dt * beta * (3 + 2 * dt * beta))));
    }

    else {
        return;
    }
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
void BassetForceTools::AddFre(array_1d<double, 3>& old_Fi, const double beta, const double dt)
{
    const double exp_coeff = std::exp(beta * dt);
    old_Fi[0] = exp_coeff * old_Fi[0];
    old_Fi[1] = exp_coeff * old_Fi[1];
    old_Fi[2] = exp_coeff * old_Fi[2];
}
//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void BassetForceTools::AppendIntegrands(ModelPart& r_model_part)
{
    ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[LAST_TIME_APPENDING] = r_model_part.GetProcessInfo()[TIME];

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        vector<double>& historic_integrands             = inode->GetValue(BASSET_HISTORIC_INTEGRANDS);
        const array_1d<double, 3>& fluid_vel_projected  = inode->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel         = inode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel;
        noalias(slip_vel)                               = fluid_vel_projected - particle_vel;
        int n = historic_integrands.size();

        historic_integrands.resize(n + 3);
        historic_integrands.insert_element(n,     slip_vel[0]);
        historic_integrands.insert_element(n + 1, slip_vel[1]);
        historic_integrands.insert_element(n + 2, 0.0);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void BassetForceTools::AppendIntegrandsImplicit(ModelPart& r_model_part)
{
    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    process_info[LAST_TIME_APPENDING] = r_model_part.GetProcessInfo()[TIME];

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        vector<double>& historic_integrands             = inode->GetValue(BASSET_HISTORIC_INTEGRANDS);
        const array_1d<double, 3>& fluid_vel_projected  = inode->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel         = inode->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel;
        noalias(slip_vel)                               = fluid_vel_projected - particle_vel;
        int n = historic_integrands.size();

        if (mFirstTimeAppending){
            mFirstTimeAppending = false;
            historic_integrands.resize(n + 9);
            historic_integrands.insert_element(n    , slip_vel[0]);
            historic_integrands.insert_element(n + 1, slip_vel[1]);
            historic_integrands.insert_element(n + 2, 0.0);
            historic_integrands.insert_element(n + 3, particle_vel[0]);
            historic_integrands.insert_element(n + 4, particle_vel[1]);
            historic_integrands.insert_element(n + 5, 0.0);
            historic_integrands.insert_element(n + 6, particle_vel[0]);
            historic_integrands.insert_element(n + 7, particle_vel[1]);
            historic_integrands.insert_element(n + 8, 0.0);
        }
        else {
            array_1d<double, 3> old_particle_vel;
            old_particle_vel[0] = historic_integrands[n - 3];
            old_particle_vel[1] = historic_integrands[n - 2];
            old_particle_vel[2] = historic_integrands[n - 1];
            historic_integrands.resize(n + 3);
            historic_integrands.insert_element(n - 6, slip_vel[0]);
            historic_integrands.insert_element(n - 5, slip_vel[1]);
            historic_integrands.insert_element(n - 4, 0.0);
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
    double time = r_process_info[TIME];
    const double delta_time = time - r_process_info[LAST_TIME_APPENDING];
    r_process_info[LAST_TIME_APPENDING] = time;
    const double t_win = SphericSwimmingParticle<SphericParticle>::mTimeWindow;

    if (r_process_info[BASSET_FORCE_TYPE] == 3){
        const std::vector<double>& Ts    = SphericSwimmingParticle<SphericParticle>::mTs;
        const std::vector<double>& Betas = SphericSwimmingParticle<SphericParticle>::mBetas;
        const double t_win               = SphericSwimmingParticle<SphericParticle>::mTimeWindow;

        const int order = r_process_info[QUADRATURE_ORDER];

        for (ElementIterator iparticle = r_model_part.ElementsBegin(); iparticle != r_model_part.ElementsEnd(); iparticle++){
            Node<3>& node = iparticle->GetGeometry()[0];
            double initial_time;
            iparticle->Calculate(TIME, initial_time, r_process_info);

            if (time - initial_time > t_win){
                vector<double>& historic_integrands         = node.GetValue(BASSET_HISTORIC_INTEGRANDS);
                vector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
                int m = hinsberg_tail_contributions.size() / 3;
                array_1d<double, 3> Fi = ZeroVector(3);

                for (int i = 0; i < m; i++){
                    double ti = Ts[i];
                    double beta = Betas[i];
                    Fi[0] = hinsberg_tail_contributions[3 * i];
                    Fi[1] = hinsberg_tail_contributions[3 * i + 1];
                    Fi[2] = hinsberg_tail_contributions[3 * i + 2];
                    AddFre(Fi, beta, delta_time);
                    AddFdi(order, Fi, t_win, ti, beta, delta_time, historic_integrands);
                    hinsberg_tail_contributions[3 * i]     = Fi[0];
                    hinsberg_tail_contributions[3 * i + 1] = Fi[1];
                    hinsberg_tail_contributions[3 * i + 2] = Fi[2];
                }
            }
        }
    }

    for (ElementIterator iparticle = r_model_part.ElementsBegin(); iparticle != r_model_part.ElementsEnd(); iparticle++){
        Node<3>& node = iparticle->GetGeometry()[0];
        vector<double>& historic_integrands             = node.GetValue(BASSET_HISTORIC_INTEGRANDS);
        const array_1d<double, 3>& fluid_vel_projected  = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        const array_1d<double, 3>& particle_vel         = node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> slip_vel;
        noalias(slip_vel)                               = fluid_vel_projected - particle_vel;
        int n = historic_integrands.size();
        double initial_time;
        iparticle->Calculate(TIME, initial_time, r_process_info);
        if (time - initial_time <= t_win){ // list of integrands still growing
            historic_integrands.resize(n + 3);
            historic_integrands.insert_element(n,     slip_vel[0]);
            historic_integrands.insert_element(n + 1, slip_vel[1]);
            historic_integrands.insert_element(n + 2, 0.0);
        }

        else { // forget oldest integrand
            for (int i = 0; i < n / 3 - 1; i++){
                historic_integrands[3 * i]     = historic_integrands[3 * (i + 1)];
                historic_integrands[3 * i + 1] = historic_integrands[3 * (i + 1) + 1];
                historic_integrands[3 * i + 2] = historic_integrands[3 * (i + 1) + 2];
            }

            historic_integrands.insert_element(n - 3, slip_vel[0]);
            historic_integrands.insert_element(n - 2, slip_vel[1]);
            historic_integrands.insert_element(n - 1, 0.0);
        }
    }
}
} // namespace Kratos
