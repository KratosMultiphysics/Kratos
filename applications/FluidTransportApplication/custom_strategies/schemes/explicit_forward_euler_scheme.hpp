//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Albert Puigferrat Perez
//                   Ignasi de Pouplana
//

#if !defined(KRATOS_EXPLICIT_FORWARD_EULER_SCHEME )
#define  KRATOS_EXPLICIT_FORWARD_EULER_SCHEME

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
// #include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/schemes/generalized_newmark_GN11_scheme.hpp"
#include "includes/convection_diffusion_settings.h"

// Application includes
#include "fluid_transport_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>

class ExplicitForwardEulerScheme : public GeneralizedNewmarkGN11Scheme<TSparseSpace,TDenseSpace>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( ExplicitForwardEulerScheme );

    typedef GeneralizedNewmarkGN11Scheme<TSparseSpace,TDenseSpace>                      BaseType;
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    ExplicitForwardEulerScheme(double theta) : GeneralizedNewmarkGN11Scheme<TSparseSpace,TDenseSpace>(theta)
    {
    }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~ExplicitForwardEulerScheme() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Predict(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        // this->UpdateVariablesDerivatives(r_model_part);
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void Update(
        ModelPart& r_model_part,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];

            //Update Phi (DOFs)
            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {
                if (itDof->IsFree())
                    itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx, itDof->EquationId());
            }
        }

        //this->UpdateVariablesDerivatives(r_model_part);

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // inline void UpdateVariablesDerivatives(ModelPart& r_model_part)
    // {
    //     KRATOS_TRY

    //     ConvectionDiffusionSettings::Pointer my_settings = r_model_part.GetProcessInfo().GetValue(CONVECTION_DIFFUSION_SETTINGS);
    //     const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    //     const int NNodes = static_cast<int>(r_model_part.Nodes().size());
    //     ModelPart::NodesContainerType::iterator node_begin = r_model_part.NodesBegin();

    //     #pragma omp parallel for
    //     for(int i = 0; i < NNodes; i++)
    //     {
    //         ModelPart::NodesContainerType::iterator itNode = node_begin + i;

    //         double& CurrentPhi = itNode->FastGetSolutionStepValue(TEMPERATURE);
    //         const double& PreviousPhi = itNode->FastGetSolutionStepValue(TEMPERATURE, 1);
    //         const double& CurrentPhiTheta = itNode->FastGetSolutionStepValue(rUnknownVar);

    //         CurrentPhi = 1.0 / mTheta * CurrentPhiTheta + (1.0 - 1.0 / mTheta) * PreviousPhi;
    //     }

    //     KRATOS_CATCH( "" )
    // }

}; // Class ExplicitForwardEulerScheme
}  // namespace Kratos

#endif // KRATOS_EXPLICIT_FORWARD_EULER_SCHEME defined
