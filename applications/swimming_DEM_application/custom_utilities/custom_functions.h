//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_CUSTOM_FUNCTIONS)
#define KRATOS_CUSTOM_FUNCTIONS

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif
// System includes

// Project includes
#include "includes/model_part.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"

//Database includes
#include "custom_utilities/discrete_particle_configure.h"
#include "discrete_particle_configure.h"
#include "includes/define.h"
#include "../../DEM_application/custom_elements/discrete_element.h"
#include "custom_elements/spheric_swimming_particle.h"
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "../../DEM_application/custom_elements/spheric_particle.h"
#include "../swimming_DEM_application.h"
#include "../../../kratos/utilities/geometry_utilities.h"

namespace Kratos
{
template <std::size_t TDim>
class CustomFunctionsCalculator
{
public:

typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
typedef ModelPart::NodesContainerType::iterator     NodeIterator;

KRATOS_CLASS_POINTER_DEFINITION(CustomFunctionsCalculator);

CustomFunctionsCalculator(): mPressuresFilled(false){}
/// Calculator

virtual ~CustomFunctionsCalculator(){}

/// Default calculator

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

void CalculatePressureGradient(ModelPart& r_model_part)
{
    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        noalias(inode->FastGetSolutionStepValue(PRESSURE_GRADIENT)) = ZeroVector(3);
    }

    array_1d <double, 3> grad = ZeroVector(3); // its dimension is always 3
    array_1d <double, TDim + 1 > elemental_pressures;
    array_1d <double, TDim + 1 > N; // shape functions vector
    boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;

    for (ModelPart::ElementIterator ielem = r_model_part.ElementsBegin(); ielem != r_model_part.ElementsEnd(); ielem++){
        // computing the shape function derivatives
        Geometry<Node<3> >& geom = ielem->GetGeometry();
        double Volume;

        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

        // getting the pressure gradients;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            elemental_pressures[i] = geom[i].FastGetSolutionStepValue(PRESSURE);
        }


        array_1d <double, TDim> grad_aux = prod(trans(DN_DX), elemental_pressures); // its dimension may be 2

        for (unsigned int i = 0; i < TDim; ++i){
            grad[i] = grad_aux[i];
        }

        double nodal_area = Volume / static_cast<double>(TDim + 1);
        grad *= nodal_area;

        for (unsigned int i = 0; i < TDim + 1; ++i){
            geom[i].FastGetSolutionStepValue(PRESSURE_GRADIENT) += grad;
        }
    }

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
        inode->FastGetSolutionStepValue(PRESSURE_GRADIENT) /= inode->FastGetSolutionStepValue(NODAL_AREA);
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

bool AssessStationarity(ModelPart& r_model_part, const double& tol)
{
    if (!mPressuresFilled){
        PerformFirstStepComputations(r_model_part);

        return(false);
    }

    else {
        double max_pressure_change_rate = 0.0; // measure of stationarity
        double mean_celerity = 0.0;            // used to adimensionalize the time step

        // filling up mPressures and calculating the mean velocities and the maximum nodal pressure change

        unsigned int i = 0;

        for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++){
            array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(VELOCITY);
            mean_celerity += sqrt(velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2]);

            double aux    = mPressures[i];
            mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
            aux           = fabs(aux - mPressures[i]);

            if (aux > max_pressure_change_rate){
                max_pressure_change_rate = aux;
            }

            ++i;
        }

        mean_celerity /= i;
        double delta_t = r_model_part.GetProcessInfo()[TIME] - mLastMeasurementTime;

        if (delta_t > 0.0){
            // calculating coefficients for adimensionalization of the pressure change rate
            double pressure_variation;
            CalculateVariationWithingVector(mPressures, pressure_variation);
            double char_length         = pow(mTotalVolume, 1/3); // characteristic length of the model. Should be improved: a hydraulic radius or such
            double time_adim_coeff     = mean_celerity / char_length;
            double pressure_adim_coeff = 0.5 * (pressure_variation + mLastPressureVariation);
            mLastPressureVariation     = pressure_variation;

            if (pressure_adim_coeff == 0.0 || time_adim_coeff == 0.0){ // unlikely
                std::cout << "Uniform problem: stationarity check being performed with dimensional values...! " << "\n";

                if (max_pressure_change_rate <= tol){ // go with the absolute value
                    return true;
                }
            }

            max_pressure_change_rate /= time_adim_coeff * delta_t * pressure_adim_coeff ;
        }

        else {
            KRATOS_THROW_ERROR(std::runtime_error,"Trying to calculate pressure variations between two coincident time steps! (null time variation since last recorded time)","");
        }

        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";
        std::cout << "The stationarity condition tolerance is " << "\n";
        KRATOS_WATCH(tol)
        std::cout << "The stationarity residual is now " << "\n";
        KRATOS_WATCH(max_pressure_change_rate)
        std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << "\n";

        if (max_pressure_change_rate <= tol){
            return true;
        }

        else {
            return false;
        }
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateDomainVolume(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    double added_volume = 0.0;

    #pragma omp parallel for reduction(+ : added_volume)
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            added_volume += CalculateElementalVolume(it->GetGeometry());
        }
    }

    return added_volume;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

void CalculateTotalHydrodynamicForceOnParticles(ModelPart& r_dem_model_part, array_1d <double, 3>& force)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_dem_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    std::vector<array_1d <double, 3> > added_force_vect;
    added_force_vect.resize(OpenMPUtils::GetNumThreads());

    for (unsigned int k = 0; k < added_force_vect.size(); ++k){
        added_force_vect[k] = ZeroVector(3);
    }

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_dem_model_part, k); it != GetElementPartitionEnd(r_dem_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            array_1d <double, 3> element_force;

            if (geom[0].SolutionStepsDataHas(HYDRODYNAMIC_FORCE)){
                element_force = geom[0].FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);
            }

            else {
                element_force = ZeroVector(3);
            }

            added_force_vect[k] += element_force;
        }
    }

    force = added_force_vect[0];

    for (unsigned int k = 1; k < added_force_vect.size(); ++k){
        force += added_force_vect[k];
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

void CalculateTotalHydrodynamicForceOnFluid(ModelPart& r_fluid_model_part, array_1d <double, 3>& force)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    std::vector<array_1d <double, 3> > added_force_vect;
    added_force_vect.resize(OpenMPUtils::GetNumThreads());

    for (unsigned int k = 0; k < added_force_vect.size(); ++k){
        added_force_vect[k] = ZeroVector(3);
    }

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            double element_volume;
            array_1d <double, 3> element_force;

            if (geom[0].SolutionStepsDataHas(HYDRODYNAMIC_REACTION) && geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_force  = CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(geom, HYDRODYNAMIC_REACTION, element_volume);
            }

            else {
                element_force = ZeroVector(3);
            }

            added_force_vect[k] += element_force;
        }
    }

    force = added_force_vect[0];

    for (unsigned int k = 1; k < added_force_vect.size(); ++k){
        force += added_force_vect[k];
    }
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************
// this function assumes linear elements are used

double CalculateGlobalFluidVolume(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_fluid_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

    double added_fluid_volume = 0.0;

    #pragma omp parallel for reduction(+ : added_fluid_volume)
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

        for (ElementIterator it = GetElementPartitionBegin(r_fluid_model_part, k); it != GetElementPartitionEnd(r_fluid_model_part, k); ++it){
            Geometry< Node<3> >& geom = it->GetGeometry();
            double element_volume;
            double element_fluid_volume;

            if (geom[0].SolutionStepsDataHas(FLUID_FRACTION)){
                element_fluid_volume = CalculateScalarIntegralOfLinearInterpolation(geom, FLUID_FRACTION, element_volume);
            }

            else {
                element_fluid_volume = CalculateElementalVolume(geom);
            }

            added_fluid_volume += element_fluid_volume;
        }
    }

    return added_fluid_volume;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

bool mPressuresFilled;
double mLastMeasurementTime;
double mLastPressureVariation;
double mTotalVolume;
std::vector<double> mPressures;

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

inline double CalculateVol(const double x0, const double y0, const double z0,
                           const double x1, const double y1, const double z1,
                           const double x2, const double y2, const double z2,
                           const double x3, const double y3, const double z3)
{
    double x10 = x1 - x0;
    double y10 = y1 - y0;
    double z10 = z1 - z0;

    double x20 = x2 - x0;
    double y20 = y2 - y0;
    double z20 = z2 - z0;

    double x30 = x3 - x0;
    double y30 = y3 - y0;
    double z30 = z3 - z0;

    double detJ = x10 * y20 * z30 - x10 * y30 * z20 +
                  y10 * z20 * x30 - y10 * x20 * z30 +
                  z10 * x20 * y30 - z10 * y20 * x30;

    return  detJ * 0.1666666666666666666666667;
}

//***************************************************************************************************************
//***************************************************************************************************************

double CalculateElementalVolume(const Geometry<Node <3> >& geom)
{
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "element with zero area found with the current geometry ", geom);
    }

    return vol;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

double CalculateScalarIntegralOfLinearInterpolation(const Geometry<Node < 3 > >& geom, const Variable<double>& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    double value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var, 0);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

array_1d <double, 3> CalculateVectorIntegralOfLinearInterpolation(const Geometry<Node < 3 > >& geom, const Variable<array_1d <double, 3> >& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    array_1d <double, 3> value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

array_1d <double, 3> CalculateVectorIntegralOfLinearInterpolationPerUnitFluidMass(const Geometry<Node < 3 > >& geom, const Variable<array_1d <double, 3> >& r_var, double& vol)
{
    array_1d<double, 4> N;
    double x0 = geom[0].X();
    double y0 = geom[0].Y();
    double z0 = geom[0].Z();
    double x1 = geom[1].X();
    double y1 = geom[1].Y();
    double z1 = geom[1].Z();
    double x2 = geom[2].X();
    double y2 = geom[2].Y();
    double z2 = geom[2].Z();
    double x3 = geom[3].X();
    double y3 = geom[3].Y();
    double z3 = geom[3].Z();

    double xc = 0.25 * (x0 + x1 + x2 + x3);
    double yc = 0.25 * (y0 + y1 + y2 + y3);
    double zc = 0.25 * (z0 + z1 + z2 + z3);

    vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    if (vol == 0.0){
        KRATOS_THROW_ERROR(std::logic_error, "Element with zero area found. Its geometry is given by", geom);
    }

    N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc);
    N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc);
    N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc);
    N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc);

    array_1d <double, 3> value_at_gauss_point = N[0] * geom[0].FastGetSolutionStepValue(r_var) * geom[0].FastGetSolutionStepValue(DENSITY) * geom[0].FastGetSolutionStepValue(FLUID_FRACTION);

    for (unsigned int i = 1; i != 4; ++i){
        value_at_gauss_point += N[i] * geom[i].FastGetSolutionStepValue(r_var) * geom[i].FastGetSolutionStepValue(DENSITY) * geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
    }

    return value_at_gauss_point;
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


void PerformFirstStepComputations(ModelPart& r_model_part)
{
    mTotalVolume = CalculateDomainVolume(r_model_part);
    mPressures.resize(r_model_part.Nodes().size());
    mLastMeasurementTime = r_model_part.GetProcessInfo()[TIME];

    unsigned int i = 0;

    for (NodeIterator inode = r_model_part.NodesBegin(); inode != r_model_part.NodesEnd(); inode++) {
        mPressures[i] = inode->FastGetSolutionStepValue(PRESSURE);
        ++i;
    }

    mPressuresFilled = true;
    CalculateVariationWithingVector(mPressures, mLastPressureVariation);
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

///@}
///@name Member r_variables
///@{
vector<unsigned int> mElementsPartition;

///@}
///@name Un accessible methods
///@{

inline void CalculateVariationWithingVector(const std::vector<double>& vector, double& variation)
{
    double min = vector[0];
    double max = vector[0];

    for (unsigned int i = 0; i != vector.size(); ++i){
        min = std::min(min, mPressures[i]);
        max = std::max(max, mPressures[i]);
    }

    variation = max - min;
}

vector<unsigned int>& GetElementPartition()
{
    return mElementsPartition;
}

ElementIterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
{
    return r_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin() + mElementsPartition[k];
}

ElementIterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
{
    return r_model_part.GetCommunicator().LocalMesh().Elements().ptr_begin() + mElementsPartition[k + 1];
}

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

}; // Class CustomFunctionsCalculator
} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


