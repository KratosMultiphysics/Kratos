// Author: Guillermo Casas (gcasas@cimne.upc.edu)
// Date: February 2019

#if !defined(SDEM_BOUSSINESQ_BASSET_HISTORY_FORCE_LAW_H_INCLUDED)
#define SDEM_BOUSSINESQ_BASSET_HISTORY_FORCE_LAW_H_INCLUDED

#include "history_force_law.h"

namespace Kratos {

    class KRATOS_API(SWIMMING_DEM_APPLICATION) BoussinesqBassetHistoryForceLaw : public HistoryForceLaw {

    public:
        typedef Node NodeType;
        KRATOS_CLASS_POINTER_DEFINITION(BoussinesqBassetHistoryForceLaw);

        // TODO: make mDoApplyFaxenCorrections an option
        BoussinesqBassetHistoryForceLaw():
            mDoApplyFaxenCorrections(false),
            mBassetForceType(4),
            mQuadratureOrder(2),
            mOldDaitchePresentCoefficient(0.0),
            mOldBassetTerm(ZeroVector(3)){}

        BoussinesqBassetHistoryForceLaw(Parameters r_parameters);

        ~BoussinesqBassetHistoryForceLaw(){}

        HistoryForceLaw::Pointer Clone() const override;

        std::string GetTypeOfLaw() override;

        void ComputeForce(Geometry<Node >& r_geometry,
                          const double reynolds_number,
                          double particle_radius,
                          double fluid_density,
                          double fluid_kinematic_viscosity,
                          array_1d<double, 3>& minus_slip_velocity,
                          array_1d<double, 3>& history_force,
                          const ProcessInfo& r_current_process_info) override;

        // variables for Daitche's method
        static std::vector<double> mAjs;
        static std::vector<double> mBns;
        static std::vector<double> mCns;
        static std::vector<double> mDns;
        static std::vector<double> mEns;
        static bool mDaitcheVectorsAreFull;
        // variables for Hinsberg's method
        static double mTimeWindow;
        static std::vector<double> mAs;
        static std::vector<double> mTs;
        static std::vector<double> mAlphas;
        static std::vector<double> mBetas;

    private:
        bool mDoApplyFaxenCorrections;
        int mBassetForceType;
        int mQuadratureOrder;
        double mOldDaitchePresentCoefficient;
        array_1d<double, 3> mOldBassetTerm;


        double GetDaitcheCoefficient(int order, unsigned int n, unsigned int j, const double last_h_over_h, const int n_steps_per_quad_step);
        void CalculateExplicitFractionalDerivative(NodeType& node, array_1d<double, 3>& fractional_derivative, double& present_coefficient, Vector& historic_integrands, const double last_h_over_h, const int n_steps_per_quad_step);
        void AddHinsbergTailContribution(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands);
        void AddHinsbergTailContributionStrict(NodeType& node, array_1d<double, 3>& fractional_derivative_of_slip_vel, const int order, const int n_steps_per_quad_step, const double time, const double delta_time, const double last_h_over_h, DenseVector<double>& historic_integrands);
        double Phi(const double x);
        double Ki(const double alpha, const double beta, const double time);
        void AddFdi(const int order, array_1d<double, 3>& F, const double t_win, const double alpha, const double beta, const double last_h_over_h, const double delta_time, const DenseVector<double>& historic_integrands, const array_1d<double, 3>& oldest_integrand);
        void AddFre(array_1d<double, 3>& old_Fi, const double beta, const double delta_time);

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HistoryForceLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HistoryForceLaw)
        }

    }; //class BoussinesqBassetHistoryForceLaw

} // Namespace Kratos

#endif /* SDEM_BOUSSINESQ_BASSET_HISTORY_FORCE_LAW_H_INCLUDED  defined */
