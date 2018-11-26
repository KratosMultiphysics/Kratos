//
//   Project Name:        Kratos
//   Last Modified by:    $Author: G.Casas (gcasas@cimmne.upc.edu) $
//   Date:                $Date: 2016-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_SAMPLING_TOOL)
#define KRATOS_SAMPLING_TOOL

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
#include "includes/define.h"
#include "../../../DEMApplication/custom_elements/discrete_element.h"
#include "custom_elements/spheric_swimming_particle.h"

namespace Kratos
{
class SamplingTool
{
public:

typedef ModelPart::ElementsContainerType::iterator  ElementIterator;
typedef ModelPart::NodesContainerType::iterator     NodeIterator;
typedef ModelPart::NodesContainerType               NodesArrayType;

KRATOS_CLASS_POINTER_DEFINITION(SamplingTool);

SamplingTool(){}
/// Calculator

virtual ~SamplingTool(){}

/// Default calculator

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************


//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

private:

//**************************************************************************************************************************************************
//**************************************************************************************************************************************************

///@}
///@name Member r_variables
///@{
DenseVector<unsigned int> mElementsPartition;

///@}
///@name Un accessible methods
///@{


DenseVector<unsigned int>& GetElementPartition()
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

}; // Class SamplingTool
} // namespace Kratos.

#endif // KRATOS_SAMPLING_TOOL  defined


