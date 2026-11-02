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

#if !defined(KRATOS_TRANSFER_INITIAL_STRESS_UTILITY)
#define KRATOS_TRANSFER_INITIAL_STRESS_UTILITY

// Project includes
#include "includes/model_part.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class TransferInitialStressUtility
{

  public:
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Constructor
    TransferInitialStressUtility() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    ~TransferInitialStressUtility() {}

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

#endif /* KRATOS_TRANSFER_INITIAL_STRESS_UTILITY defined */
