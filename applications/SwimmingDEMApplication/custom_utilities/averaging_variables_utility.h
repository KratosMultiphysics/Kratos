//   $Author: Joaquin Gonzalez-Usua
#if !defined(KRATOS_AVERAGING_VARIABLES_UTILITY)
#define  KRATOS_AVERAGING_VARIABLES_UTILITY

// System includes

#include <limits>
#include <iostream>
#include <iomanip>
#include <list>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "DEM_application_variables.h"


/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) AveragingVariablesUtility {

public:

typedef ModelPart::NodesContainerType::iterator               NodeIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(AveragingVariablesUtility);

/// Default constructor

AveragingVariablesUtility(){}

/// Destructor

virtual ~AveragingVariablesUtility(){}

double AverageVariables(ModelPart& r_spheres_model_part,
                      const double centroid_current_layer,
                      const double layer_width,
                      const double number_of_sublayers,
                      const double plane_area,
                      Vector& layer_averaged_velocity,
                      Vector& layer_averaged_dv_dz,
                      Matrix& layer_averaged_particle_stress,
                      Matrix& layer_granular_temperature);

void CalculateVelocityVariance(ModelPart& r_spheres_model_part,
                               Matrix& layer_velocity_variance,
                               Vector sublayer_averaged_velocity,
                               double& centroid_current_sublayer);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:


/// Assignment operator
AveragingVariablesUtility & operator=(AveragingVariablesUtility const& rOther);

}; // Class AveragingVariablesUtility

} // namespace Kratos.

#endif // KRATOS_AVERAGING_VARIABLES_UTILITY