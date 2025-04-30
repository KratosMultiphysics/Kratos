#if !defined(KRATOS_POISEUILLE_TORUS_FLOW_FIELD_H)
#define KRATOS_POISEUILLE_TORUS_FLOW_FIELD_H

// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "velocity_field.h"

namespace Kratos
{
    class KRATOS_API(SWIMMING_DEM_APPLICATION) PoiseuilleTorusFlowField : public VelocityField
    {

    public:
        KRATOS_CLASS_POINTER_DEFINITION(PoiseuilleTorusFlowField);

        /// Default constructor.
        PoiseuilleTorusFlowField() : VelocityField(), mMajorRadius(1.), mMinorRadius(.5), mU0(1.)
        {
            unsigned int number_of_threads = ParallelUtilities::GetNumThreads();
            ResizeVectorsForParallelism(number_of_threads);
        };

        PoiseuilleTorusFlowField(const double major_radius, const double minor_radius, const double u0) : VelocityField(), mMajorRadius(major_radius), mMinorRadius(minor_radius), mU0(u0)
        {
            unsigned int number_of_threads = ParallelUtilities::GetNumThreads();
            // std::cout << "Resizing..." << std::endl;
            ResizeVectorsForParallelism(number_of_threads);
            // std::cout << "Done." << std::endl;
        }

        /// Destructor.
        virtual ~PoiseuilleTorusFlowField() {}

        //***************************************************************************************************************
        //***************************************************************************************************************
        void ResizeVectorsForParallelism(const int n_threads) override;

        void UpdateCoordinates(const double time, const array_1d<double, 3>& coor, const int i_thread = 0) override;

        void UpdateCoordinates(const double time, const DenseVector<double>& coor, const int i_thread = 0) override;

        void LockCoordinates(const int i_thread = 0) override;

        void UnlockCoordinates(const int i_thread = 0) override;
        //***************************************************************************************************************
        //***************************************************************************************************************

        // virtual void CalculateMaterialAcceleration(const double time, const array_1d<double, 3> &coor, array_1d<double, 3> &accel, const int i_thread = 0) override;

        // Values

        double U0(const int i_thread = 0) override;
        double U1(const int i_thread = 0) override;
        double U2(const int i_thread = 0) override;

        // First-order derivatives

        double U0DT(const int i_thread = 0) override;
        double U0D0(const int i_thread = 0) override;
        double U0D1(const int i_thread = 0) override;
        double U0D2(const int i_thread = 0) override;

        double U1DT(const int i_thread = 0) override;
        double U1D0(const int i_thread = 0) override;
        double U1D1(const int i_thread = 0) override;
        double U1D2(const int i_thread = 0) override;

        double U2DT(const int i_thread = 0) override;
        double U2D0(const int i_thread = 0) override;
        double U2D1(const int i_thread = 0) override;
        double U2D2(const int i_thread = 0) override;

        // Second-order derivatives

        double U0DTDT(const int i_thread = 0) override;
        double U0DTD0(const int i_thread = 0) override;
        double U0DTD1(const int i_thread = 0) override;
        double U0DTD2(const int i_thread = 0) override;
        double U0D0D0(const int i_thread = 0) override;
        double U0D0D1(const int i_thread = 0) override;
        double U0D0D2(const int i_thread = 0) override;
        double U0D1D1(const int i_thread = 0) override;
        double U0D1D2(const int i_thread = 0) override;
        double U0D2D2(const int i_thread = 0) override;

        double U1DTDT(const int i_thread = 0) override;
        double U1DTD0(const int i_thread = 0) override;
        double U1DTD1(const int i_thread = 0) override;
        double U1DTD2(const int i_thread = 0) override;
        double U1D0D0(const int i_thread = 0) override;
        double U1D0D1(const int i_thread = 0) override;
        double U1D0D2(const int i_thread = 0) override;
        double U1D1D1(const int i_thread = 0) override;
        double U1D1D2(const int i_thread = 0) override;
        double U1D2D2(const int i_thread = 0) override;

        double U2DTDT(const int i_thread = 0) override;
        double U2DTD0(const int i_thread = 0) override;
        double U2DTD1(const int i_thread = 0) override;
        double U2DTD2(const int i_thread = 0) override;
        double U2D0D0(const int i_thread = 0) override;
        double U2D0D1(const int i_thread = 0) override;
        double U2D0D2(const int i_thread = 0) override;
        double U2D1D1(const int i_thread = 0) override;
        double U2D1D2(const int i_thread = 0) override;
        double U2D2D2(const int i_thread = 0) override;

    private:
        ///@}
        ///@name Member r_variables
        ///@{
        double mMajorRadius;
        double mMinorRadius;
        double mU0;

        std::vector<int> mCoordinatesAreUpToDate;
        std::vector<double> mXYDistance;
        std::vector<double> mRho;
        std::vector<double> mCos;
        std::vector<double> mSin;
        std::vector<double> mZ;
        std::vector<double> mCommonTerm;

        // double mX;
        // double mY;
        // double mZ;

        /// Assignment operator.
        PoiseuilleTorusFlowField & operator=(PoiseuilleTorusFlowField const& rOther);

    }; // Class PoiseuilleTorusFlowField

} // Namespace Kratos

#endif