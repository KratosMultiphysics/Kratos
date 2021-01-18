//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

#if !defined(KRATOS_TRANSFER_SELFWEIGHT_STRESS_UTILITIES)
#define KRATOS_TRANSFER_SELFWEIGHT_STRESS_UTILITIES

// System includes
#include <cmath>

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class TransferSelfweightStressUtility
{

  public:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    TransferSelfweightStressUtility() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~TransferSelfweightStressUtility() {}

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Transforming local stress vector in global coordinates
    void Transfer(ModelPart &selfweight_model_part, ModelPart &r_model_part, const int TDim)
    {

        KRATOS_TRY;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();
        ModelPart::NodesContainerType::iterator node_begin_selfweight = selfweight_model_part.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;
            ModelPart::NodesContainerType::iterator itNodeSelf = node_begin_selfweight + i;

            const Matrix &NodalStressSelfweight = itNodeSelf->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
            Matrix &NodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);

            for (int j = 0; j < TDim; j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    NodalStress(j, k) += NodalStressSelfweight(j, k);
                }
            }
        }

        KRATOS_CATCH("");
    }

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Transforming local stress vector in global coordinates
    void TransferInitialStress(ModelPart &r_model_part, const int TDim)
    {

        KRATOS_TRY;

        const int NNodes = static_cast<int>(r_model_part.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            const Matrix &InitialNodalStress = itNode->FastGetSolutionStepValue(INITIAL_NODAL_CAUCHY_STRESS_TENSOR);
            Matrix &NodalStress = itNode->FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);

            for (int j = 0; j < TDim; j++)
            {
                for (int k = 0; k < TDim; k++)
                {
                    NodalStress(j, k) += InitialNodalStress(j, k);
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  protected:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

}; //Class

} /* namespace Kratos.*/

#endif /* KRATOS_TRANSFER_SELFWEIGHT_STRESS_UTILITIES defined */