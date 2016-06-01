//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas (gcasas@cimmne.upc.edu) $
//   Date:                $Date: 2016-5-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_BASSET_FORCE_TOOLS)
#define KRATOS_BASSET_FORCE_TOOLS

// /* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

// System includes
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "processes/find_elements_neighbours_process.h"
#include "processes/find_nodal_neighbours_process.h"

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
class BassetForceTools
{
public:

typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
typedef ModelPart::NodesContainerType::iterator     NodeIterator;
typedef ModelPart::NodesContainerType               NodesArrayType;

KRATOS_CLASS_POINTER_DEFINITION(BassetForceTools);

BassetForceTools(): mFirstTimeAppending(true){}
/// Calculator

virtual ~BassetForceTools(){}

/// Default calculator

void FillDaitcheVectors(int N, int order);
void FillHinsbergVectors(ModelPart& r_model_part, int m, const double time_window);
void AppendIntegrands(ModelPart& r_model_part);
void AppendIntegrandsImplicit(ModelPart& r_model_part);
void AppendIntegrandsWindow(ModelPart& r_model_part);

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

bool mFirstTimeAppending;

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

///@}
///@name Member r_variables
///@{
vector<unsigned int> mElementsPartition;

///@}
///@name Un accessible methods
///@{

double Phi(const double x);
void AddFdi(int order, array_1d<double, 3>& F, const double t_win, const double ti, const double beta, const double delta_time, vector<double>& historic_integrands);
void AddFre(array_1d<double, 3>& old_Fi, const double beta, const double delta_time);

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

}; // Class BassetForceTools
} // namespace Kratos.

#endif // KRATOS_CREATE_AND_DESTROY  defined


