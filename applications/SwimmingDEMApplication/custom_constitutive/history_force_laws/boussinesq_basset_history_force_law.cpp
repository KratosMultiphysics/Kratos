// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#include "swimming_DEM_application.h"
#include "boussinesq_basset_history_force_law.h"

namespace Kratos {

    HistoryForceLaw::Pointer BoussinesqBassetHistoryForceLaw::Clone() const {
        BoussinesqBassetHistoryForceLaw::Pointer p_clone(new BoussinesqBassetHistoryForceLaw(*this));
        return p_clone;
    }

    BoussinesqBassetHistoryForceLaw::BoussinesqBassetHistoryForceLaw(Parameters r_parameters)
    {
        Parameters default_parameters( R"(
            {
                "name":"BoussinesqBassetHistoryForceLaw",
                "quadrature_order": 2,
                "time_steps_per_quadrature_step": 1,
                "n_init_basset_steps": 0,
                "mae_parameters": {
                    "do_use_mae": false,
                    "m": 10,
                    "window_time_interval": 0.1,
                    "type":4
                }
            }
            )" );

        r_parameters.ValidateAndAssignDefaults(default_parameters);
        mBassetForceType = 2;

        if (r_parameters["mae_parameters"]["do_use_mae"].GetBool()){
            mBassetForceType = r_parameters["mae_parameters"]["type"].GetInt();
        }

        mQuadratureOrder = r_parameters["quadrature_order"].GetInt();
        mOldBassetTerm = ZeroVector(3);
        mOldDaitchePresentCoefficient = 0.0;
    }

    std::string BoussinesqBassetHistoryForceLaw::GetTypeOfLaw() {
        std::string type_of_law = "Boussinesq-Basset history force law";
        return type_of_law;
    }

    void BoussinesqBassetHistoryForceLaw::ComputeForce(Geometry<Node >& r_geometry,
                                                       const double reynolds_number,
                                                       double particle_radius,
                                                       double fluid_density,
                                                       double fluid_kinematic_viscosity,
                                                       array_1d<double, 3>& minus_slip_velocity,
                                                       array_1d<double, 3>& basset_force,
                                                       const ProcessInfo& r_current_process_info)
    {
        const double basset_force_coeff = 6.0 * particle_radius * particle_radius * fluid_density * std::sqrt(Globals::Pi * fluid_kinematic_viscosity);
        const double delta_time = r_current_process_info[DELTA_TIME];
        int n_steps_per_quad_step = r_current_process_info[TIME_STEPS_PER_QUADRATURE_STEP];
        const double quadrature_delta_time = n_steps_per_quad_step * delta_time;
        Node& node = r_geometry[0];

        if (r_current_process_info[TIME_STEPS] >= r_current_process_info[NUMBER_OF_INIT_BASSET_STEPS]){
            DenseVector<double>& historic_integrands = node.GetValue(BASSET_HISTORIC_INTEGRANDS);
            const double time = r_current_process_info[TIME];
            const double latest_quadrature_time_step = time + delta_time - r_current_process_info[LAST_TIME_APPENDING];
            array_1d<double, 3> fractional_derivative_of_slip_vel;
            double present_coefficient;
            const double sqrt_of_quad_h_q = std::sqrt(quadrature_delta_time);
            const double last_h_over_h = latest_quadrature_time_step / quadrature_delta_time;

            CalculateExplicitFractionalDerivative(node,
                                                  fractional_derivative_of_slip_vel,
                                                  present_coefficient,
                                                  historic_integrands,
                                                  last_h_over_h,
                                                  n_steps_per_quad_step);

            if (r_current_process_info[FRAME_OF_REFERENCE_TYPE] >= 1){
                const array_1d<double, 3>& displacement = node.FastGetSolutionStepValue(DISPLACEMENT);
                const array_1d<double, 3>& displacement_old = node.FastGetSolutionStepValue(DISPLACEMENT_OLD);
                array_1d<double, 3> aux = displacement - displacement_old;
                const array_1d<double, 3>& omega_frame = r_current_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
                array_1d<double, 3> correction;
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                noalias(fractional_derivative_of_slip_vel) += present_coefficient * correction;
            }

            if (mBassetForceType == 3){
                AddHinsbergTailContribution(node, fractional_derivative_of_slip_vel, mQuadratureOrder, n_steps_per_quad_step, time, quadrature_delta_time, last_h_over_h, historic_integrands);
            }

            if (mBassetForceType == 4){
                AddHinsbergTailContributionStrict(node, fractional_derivative_of_slip_vel, mQuadratureOrder, n_steps_per_quad_step, time, quadrature_delta_time, last_h_over_h, historic_integrands);
            }

            array_1d<double, 3> basset_term    = fractional_derivative_of_slip_vel;
            array_1d<double, 3> old_basset_term;
            const array_1d<double, 3>& vel     = node.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& old_vel = node.FastGetSolutionStepValue(VELOCITY_OLD);
            old_basset_term = mOldBassetTerm + mOldDaitchePresentCoefficient * (old_vel - vel); // the second term corresponds to the part that was treated implicitly in the last step minus a part that was added but did not correspond to the basset term
            noalias(mOldBassetTerm) = basset_term;

            if (r_current_process_info[FRAME_OF_REFERENCE_TYPE] >= 1){
                array_1d<double, 3>& displacement_old = node.FastGetSolutionStepValue(DISPLACEMENT_OLD);
                array_1d<double, 3>& displacement_old_old = node.FastGetSolutionStepValue(VELOCITY_OLD_OLD);
                array_1d<double, 3> aux = mOldDaitchePresentCoefficient * (displacement_old_old - displacement_old);
                const array_1d<double, 3>& omega_frame = r_current_process_info[ANGULAR_VELOCITY_MOVING_FRAME];
                array_1d<double, 3> correction;
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                noalias(displacement_old_old) = displacement_old;
                noalias(displacement_old) = node.FastGetSolutionStepValue(DISPLACEMENT);
                // correcting the old_basset_term for the velocities at different times being subtracted
                noalias(old_basset_term) += correction;
                noalias(aux) = delta_time * (1.5 * basset_term - 0.5 * old_basset_term);
                SWIMMING_SET_TO_CROSS_OF_FIRST_TWO_3(omega_frame, aux, correction);
                // correcting the basset force for the discretization of the time derivative
                noalias(fractional_derivative_of_slip_vel) += correction;
            }

            noalias(fractional_derivative_of_slip_vel) -= old_basset_term;

            mOldDaitchePresentCoefficient = present_coefficient;
            mLastHistoryForceAddedMass = basset_force_coeff * sqrt_of_quad_h_q * present_coefficient;

            noalias(basset_force) = basset_force_coeff * sqrt_of_quad_h_q / delta_time * fractional_derivative_of_slip_vel;
        }

        else {
            basset_force = node.FastGetSolutionStepValue(BASSET_FORCE);
            mOldDaitchePresentCoefficient = 0.0;
            mOldBassetTerm = std::sqrt(quadrature_delta_time) / basset_force_coeff * basset_force;
        }
    }

    double BoussinesqBassetHistoryForceLaw::GetDaitcheCoefficient(int order, unsigned int n, unsigned int j, const double last_h_over_h, const int n_steps_per_quad_step)
    {
        const int l = (int)(last_h_over_h * n_steps_per_quad_step + 0.5) - 1;

        if (order == 1){
            if (j < n){
                return BoussinesqBassetHistoryForceLaw::mAjs[n_steps_per_quad_step * j + l];
            }
            else {
                return BoussinesqBassetHistoryForceLaw::mBns[n_steps_per_quad_step * j + l];
            }
        }

        else if (order == 2){
            if (n > 3){
                if (j < n - 1){
                    return BoussinesqBassetHistoryForceLaw::mAjs[n_steps_per_quad_step * j + l];
                }
                else if (j == n - 1){
                    return BoussinesqBassetHistoryForceLaw::mBns[n_steps_per_quad_step * j + l];
                }
                else {
                    return BoussinesqBassetHistoryForceLaw::mCns[n_steps_per_quad_step * j + l];
                }
            }
            else {
                if (n == 1){ // use formula for the phis of first order
                    BoussinesqBassetHistoryForceLaw::GetDaitcheCoefficient(1, n, j, last_h_over_h, n_steps_per_quad_step);
                }

                else if (n == 2){

                    const double sqrt_phi_plus_1 = std::sqrt(1 + last_h_over_h);

                    if (j == 0){
                        return 4 * sqrt_phi_plus_1 * (4 * last_h_over_h - 1) / (15 * last_h_over_h);
                    }
                    else if (j == 1){
                        return 4 * SWIMMING_POW_5(sqrt_phi_plus_1) / (15 * last_h_over_h);
                    }
                    else {
                        return 2 * sqrt_phi_plus_1 * (3 - 2 * last_h_over_h) / 15;
                    }
                }
                else {
                    const double sqrt_phi_plus_1 = std::sqrt(1 + last_h_over_h);
                    const double sqrt_phi_plus_2 = std::sqrt(2 + last_h_over_h);

                    if (j == 0){
                        return 4 * sqrt_phi_plus_1 * (4 * last_h_over_h - 1) / (15 * last_h_over_h);
                    }
                    else if (j == 1){
                        return 4 * SWIMMING_POW_5(sqrt_phi_plus_1) / (15 * last_h_over_h) + 2 * (4 * SWIMMING_POW_3(sqrt_phi_plus_2) * (3 + 4 * last_h_over_h) - SWIMMING_POW_3(sqrt_phi_plus_1) * (9 + 4 * last_h_over_h)) / 15;
                    }
                    else if (j == 2){
                        return 4 / 15 * (SWIMMING_POW_3(sqrt_phi_plus_2) * (2 - 4 * last_h_over_h) + sqrt_phi_plus_1 * (4 * last_h_over_h * last_h_over_h + 7 * last_h_over_h - 2));
                        //return 2 * sqrt_phi_plus_1 * (3 - 2 * last_h_over_h) / 15 + 2 * (4 * SWIMMING_POW_3(sqrt_phi_plus_2) * (2 * last_h_over_h - 1) + sqrt_phi_plus_1 * (8 * last_h_over_h * (2 + last_h_over_h) - 7)) / 15;
                    }
                    else {
                        return 2 * (sqrt_phi_plus_1 * (1 - 3 * last_h_over_h - 4 * last_h_over_h * last_h_over_h) + sqrt_phi_plus_2 * (1 + last_h_over_h + 4 * last_h_over_h * last_h_over_h)) / 15;
                    }
                }
            }
        }

        else { // not implemented with substeping yet
            if (n > 6){
                if (j < n - 3){
                    return BoussinesqBassetHistoryForceLaw::mAjs[j];
                }
                else if (j == n - 3){
                    return BoussinesqBassetHistoryForceLaw::mBns[n];
                }
                else if (j == n - 2){
                    return BoussinesqBassetHistoryForceLaw::mCns[n];
                }
                else if (j == n - 1){
                    return BoussinesqBassetHistoryForceLaw::mDns[n];
                }
                else {
                    return BoussinesqBassetHistoryForceLaw::mEns[n];
                }
            }

            else if (n == 2){ // use formula for the betas of second order
                BoussinesqBassetHistoryForceLaw::GetDaitcheCoefficient(2, n, j, last_h_over_h, n_steps_per_quad_step);
            }

            else if (n == 3){
                long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
                if (j == 0){
                    return 68 * sqrt_3_over_105;
                }
                else if (j == 1){
                    return 90 * sqrt_3_over_105;
                }
                else if (j == 2){
                    return 36 * sqrt_3_over_105;
                }
                else {
                    return 16 * sqrt_3_over_105;
                }
            }
            else if (n == 4){
                long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
                long double OneOver315 = 1.0 / 315;

                if (j == 0){
                    return 244 * sqrt_2_over_315;
                }
                else if (j == 1){
                    return 1888 * OneOver315 - 976 * sqrt_2_over_315;
                }
                else if (j == 2){
                    return - 656 * OneOver315 + 1464 * sqrt_2_over_315;
                }
                else if (j == 3){
                    return 1632 - 976 * OneOver315;
                }
                else {
                    return - 292 * OneOver315 + 244 * sqrt_2_over_315;
                }
            }
            else if (n == 5){
                long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
                long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
                long double sqrt_5_over_63  = std::sqrt(static_cast<long double>(5)) / 63;

                if (j == 0){
                    return 244 * sqrt_2_over_315;
                }
                else if (j == 1){
                    return 362 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
                }
                else if (j == 2){
                    return 500 * sqrt_5_over_63 - 1448 * sqrt_3_over_105 + 1464 * sqrt_2_over_315;
                }
                else if (j == 3){
                    return - 870 * sqrt_5_over_63 + 2172 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
                }
                else if (j == 4){
                    return 660 * sqrt_5_over_63 - 1448 * sqrt_3_over_105 + 244 * sqrt_2_over_315;
                }
                else {
                    return - 164 * sqrt_5_over_63 + 362 * sqrt_3_over_105;
                }
            }
            else {
                long double sqrt_2_over_315 = std::sqrt(static_cast<long double>(2)) / 315;
                long double sqrt_3_over_105 = std::sqrt(static_cast<long double>(3)) / 105;
                long double sqrt_6_over_105 = std::sqrt(static_cast<long double>(6)) / 105;
                long double OneOver315 = 1.0 / 315;

                if (j == 0){
                    return 244 * sqrt_2_over_315;
                }
                else if (j == 1){
                    return 362 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
                }
                else if (j == 2){
                    return 5584 * OneOver315 - 1448 * sqrt_3_over_105 + 1464 * sqrt_2_over_315;
                }
                else if (j == 3){
                    return 1720 * sqrt_6_over_105 - 22336 * OneOver315 + 2172 * sqrt_3_over_105 - 976 * sqrt_2_over_315;
                }
                else if (j == 4){
                    return - 3564 * sqrt_6_over_105 + 33504 * OneOver315 - 1448 * sqrt_3_over_105 + 244 * sqrt_2_over_315;
                }
                else if (j == 5){
                    return 2808 * sqrt_6_over_105 - 22336 * OneOver315 + 362 * sqrt_3_over_105;
                }
                else {
                    return - 754 * sqrt_6_over_105 + 5584 * OneOver315;
                }
            }
        }
        return 0.0;
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BoussinesqBassetHistoryForceLaw::CalculateExplicitFractionalDerivative(NodeType& node, array_1d<double, 3>& fractional_derivative,
                                                                                    double& present_coefficient,
                                                                                    DenseVector<double>& historic_integrands,
                                                                                    const double last_h_over_h,
                                                                                    const int n_steps_per_quad_step)
    {
        const int N = historic_integrands.size() - 3;
        const int n = (int)N / 3;
        double fast_fractional_derivative[3] = {0.0};

        for (int j = 0; j < n + 1; ++j){
            double coefficient = GetDaitcheCoefficient(mQuadratureOrder, n + 1, j + 1, last_h_over_h, n_steps_per_quad_step);
            for (int i_comp = 0; i_comp < 3; i_comp++){
                unsigned int integrand_component_position = N - 3 * j + i_comp;
                fast_fractional_derivative[i_comp] += coefficient * historic_integrands[integrand_component_position];
            }
        }
        present_coefficient = GetDaitcheCoefficient(mQuadratureOrder, n + 1, 0, last_h_over_h, n_steps_per_quad_step);
        noalias(fractional_derivative) = present_coefficient * (node.FastGetSolutionStepValue(AUX_VEL) - node.FastGetSolutionStepValue(VELOCITY));
        SWIMMING_ADD_SECOND_TO_FIRST(fractional_derivative, fast_fractional_derivative)
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    double BoussinesqBassetHistoryForceLaw::Phi(const double x)
    {
        if (std::abs(x) < 1e-10){
            return std::expm1(x) / x;
        }
        else {
            return 1 + 0.5 * x + 1.0 / 6 * x * x;
        }
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    double BoussinesqBassetHistoryForceLaw::Ki(const double alpha, const double beta, const double time)
    {
        return alpha * std::exp(beta * time);
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************


    void BoussinesqBassetHistoryForceLaw::AddFdi(const int order, array_1d<double, 3>& F, const double t_win, const double alpha, const double beta, const double phi, const double dt, const DenseVector<double>& historic_integrands, const array_1d<double, 3>& oldest_integrand)
    {
        if (order == 1){
            const double beta_dt = beta * dt;
            const double coeff = - alpha / beta * std::exp(beta * (t_win - dt + dt * phi));
            const double coeff_N = 1 - Phi(beta_dt);
            const double coeff_N_plus_1 = std::exp(beta_dt) * (Phi(- beta_dt) - 1);

            F[0] +=  coeff * (coeff_N * historic_integrands[0] + coeff_N_plus_1 * oldest_integrand[0]);
            F[1] +=  coeff * (coeff_N * historic_integrands[1] + coeff_N_plus_1 * oldest_integrand[1]);
            F[2] +=  coeff * (coeff_N * historic_integrands[2] + coeff_N_plus_1 * oldest_integrand[2]);
        }

        else if (order == 2){
            if (false){ // unstable
                const double coeff = 0.5 * alpha / (beta * SWIMMING_POW_2(dt * beta));
                const double exp_1 = std::exp(beta * (t_win - dt + dt * phi));
                const double exp_2 = std::exp(beta * dt);
                const double f20 = oldest_integrand[0];
                const double f21 = oldest_integrand[1];
                const double f22 = oldest_integrand[2];
                const double f10 = historic_integrands[0];
                const double f11 = historic_integrands[1];
                const double f12 = historic_integrands[2];
                const double f00 = historic_integrands[3];
                const double f01 = historic_integrands[4];
                const double f02 = historic_integrands[5];
                F[0] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f00
                                    + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f10
                                    + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f20);
                F[1] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f01
                                    + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f11
                                    + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f21);
                F[2] += coeff * exp_1 * ((2 * exp_2 * (dt * beta) - dt * beta - 2)                                                * f02
                                    + (exp_2 * (4 - 4 * dt * beta) + 2 * beta * SWIMMING_POW_2(dt * beta) - 4)                 * f12
                                    + (exp_2 * (2 - 3 * dt * beta + 2 * 2 * beta * SWIMMING_POW_2(dt * beta)) + dt * beta - 2) * f22);
            }
            else {
                const double t_minus_t2 = t_win + dt * phi;
                const double t_minus_t1 = t_win + dt * phi - dt;
                const double t_minus_t0 = t_win + dt * phi - 2 * dt;
                const double Ki2 = Ki(alpha, beta, t_minus_t2);
                const double Ki1 = Ki(alpha, beta, t_minus_t1);
                const double Ki0 = Ki(alpha, beta, t_minus_t0);
                const double f20 = oldest_integrand[0] * Ki2;
                const double f21 = oldest_integrand[1] * Ki2;
                const double f22 = oldest_integrand[2] * Ki2;
                const double f10 = historic_integrands[0] * Ki1;
                const double f11 = historic_integrands[1] * Ki1;
                const double f12 = historic_integrands[2] * Ki1;
                const double f00 = historic_integrands[3] * Ki0;
                const double f01 = historic_integrands[4] * Ki0;
                const double f02 = historic_integrands[5] * Ki0;
                const double aux_coeff = dt / 12;
                F[0] += aux_coeff * (- f00 + 8 * f10 + 5 * f20);
                F[1] += aux_coeff * (- f01 + 8 * f11 + 5 * f21);
                F[2] += aux_coeff * (- f02 + 8 * f12 + 5 * f22);
            }
        }
        else {
            return;
        }
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BoussinesqBassetHistoryForceLaw::AddFre(array_1d<double, 3>& old_Fi, const double beta, const double dt)
    {
        const double exp_coeff = std::exp(beta * dt);
        old_Fi[0] = exp_coeff * old_Fi[0];
        old_Fi[1] = exp_coeff * old_Fi[1];
        old_Fi[2] = exp_coeff * old_Fi[2];
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BoussinesqBassetHistoryForceLaw::AddHinsbergTailContribution(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double quadrature_delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands)
    {
        DenseVector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
        int m = hinsberg_tail_contributions.size() / 3 - 1; // number of exponentials: the last three slots hold the components of the oldest historic integrand
        const double t_win = BoussinesqBassetHistoryForceLaw::mTimeWindow;

        if (n_steps_per_quad_step * last_h_over_h < 1.5 && 2 * n_steps_per_quad_step * (time - t_win) > quadrature_delta_time){ // first calculation after appending and already later than t_win
            array_1d<double, 3> oldest_integrand;
            oldest_integrand[0] = hinsberg_tail_contributions[3 * m];
            oldest_integrand[1] = hinsberg_tail_contributions[3 * m + 1];
            oldest_integrand[2] = hinsberg_tail_contributions[3 * m + 2];
            const std::vector<double>& Ts = BoussinesqBassetHistoryForceLaw::mTs;
            const double e = std::exp(1);
            array_1d<double, 3> Fi;

            for (int i = 0; i < m; ++i){
                const double ti = Ts[i];
                const double beta = - 0.5 / ti;
                const double alpha = std::sqrt(e / ti);
                Fi[0] = hinsberg_tail_contributions[3 * i];
                Fi[1] = hinsberg_tail_contributions[3 * i + 1];
                Fi[2] = hinsberg_tail_contributions[3 * i + 2];
                AddFre(Fi, beta, quadrature_delta_time);
                AddFdi(order, Fi, t_win, alpha, beta, 1.0, quadrature_delta_time, historic_integrands, oldest_integrand);
                hinsberg_tail_contributions[3 * i]     = Fi[0];
                hinsberg_tail_contributions[3 * i + 1] = Fi[1];
                hinsberg_tail_contributions[3 * i + 2] = Fi[2];
            }
        }

        array_1d<double, 3> F_tail = ZeroVector(3);

        for (int i = 0; i < m; ++i){
            double ai = BoussinesqBassetHistoryForceLaw::mAs[i];
            F_tail[0] += ai * hinsberg_tail_contributions[3 * i];
            F_tail[1] += ai * hinsberg_tail_contributions[3 * i + 1];
            F_tail[2] += ai * hinsberg_tail_contributions[3 * i + 2];
        }

        const double sqrt_delta_time_inv = 1.0 / std::sqrt(quadrature_delta_time); // since the multiplication by sqrt(delta_time) corresponding to the F_win part is done to F_win + F_tail outside
        noalias(fractional_derivative_of_slip_vel) += sqrt_delta_time_inv * F_tail;
    }
    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void BoussinesqBassetHistoryForceLaw::AddHinsbergTailContributionStrict(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double quadrature_delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands)
    {
          DenseVector<double>& hinsberg_tail_contributions = node.GetValue(HINSBERG_TAIL_CONTRIBUTIONS);
        int m = hinsberg_tail_contributions.size() / 3 - 1; // number of exponentials: the last three slots hold the components of the oldest historic integrand

        if (m < 1){ // trivial, 0-exponentials case (there is no tail contribution)
            return;
        }

        const double t_win = BoussinesqBassetHistoryForceLaw::mTimeWindow;
        const std::vector<double>& Ts = BoussinesqBassetHistoryForceLaw::mTs;
        const double e = std::exp(1);
        const double delta_time = quadrature_delta_time / n_steps_per_quad_step;

        if (n_steps_per_quad_step * last_h_over_h < 1.5 && 2 * n_steps_per_quad_step * (time - t_win) > quadrature_delta_time){ // calculation right after last append (but at least later than t = t_win, so there is some tail)
            array_1d<double, 3> oldest_integrand;
            oldest_integrand[0] = hinsberg_tail_contributions[3 * m];
            oldest_integrand[1] = hinsberg_tail_contributions[3 * m + 1];
            oldest_integrand[2] = hinsberg_tail_contributions[3 * m + 2];
            array_1d<double, 3> Fi;

            for (int i = 0; i < m; ++i){
                const double ti = Ts[i];
                const double beta = - 0.5 / ti;
                const double alpha = std::sqrt(e / ti);
                Fi[0] = hinsberg_tail_contributions[3 * i];
                Fi[1] = hinsberg_tail_contributions[3 * i + 1];
                Fi[2] = hinsberg_tail_contributions[3 * i + 2];
                AddFre(Fi, beta, delta_time);
                AddFdi(order, Fi, t_win, alpha, beta, last_h_over_h, quadrature_delta_time, historic_integrands, oldest_integrand);
                hinsberg_tail_contributions[3 * i]     = Fi[0];
                hinsberg_tail_contributions[3 * i + 1] = Fi[1];
                hinsberg_tail_contributions[3 * i + 2] = Fi[2];
            }
        }

        else { // intermediate step: the time assigned to the tail has not changed (only an exponential factor is needed to correct for the changing current time, which affects the approximate kernel)

            for (int i = 0; i < m; ++i){
                const double ti = Ts[i];
                const double beta = - 0.5 / ti;
                const double exp_beta_dt = std::exp(beta * delta_time);
                hinsberg_tail_contributions[3 * i]     *= exp_beta_dt;
                hinsberg_tail_contributions[3 * i + 1] *= exp_beta_dt;
                hinsberg_tail_contributions[3 * i + 2] *= exp_beta_dt;
            }
        }

        array_1d<double, 3> F_tail = ZeroVector(3);
        const std::vector<double>& As = BoussinesqBassetHistoryForceLaw::mAs;

        for (int i = 0; i < m; ++i){
            const double ai = As[i];
            F_tail[0] += ai * hinsberg_tail_contributions[3 * i];
            F_tail[1] += ai * hinsberg_tail_contributions[3 * i + 1];
            F_tail[2] += ai * hinsberg_tail_contributions[3 * i + 2];
        }

        const double sqrt_delta_time_inv = 1.0 / std::sqrt(quadrature_delta_time); // since the multiplication by sqrt(delta_time) corresponding to the F_win part is done to F_win + F_tail outside
        noalias(fractional_derivative_of_slip_vel) += sqrt_delta_time_inv * F_tail;
    }


    // Definition of static variables ( this probably needs to be moved to the .h file )
    std::vector<double> BoussinesqBassetHistoryForceLaw::mAjs;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mBns;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mCns;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mDns;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mEns;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mAs;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mTs;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mAlphas;
    std::vector<double> BoussinesqBassetHistoryForceLaw::mBetas;
    double BoussinesqBassetHistoryForceLaw::mTimeWindow;
    bool BoussinesqBassetHistoryForceLaw::mDaitcheVectorsAreFull;

} // namespace Kratos

