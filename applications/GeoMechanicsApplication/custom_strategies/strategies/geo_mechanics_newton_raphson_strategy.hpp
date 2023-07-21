// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#if !defined(KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_STRATEGY)
#define KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_STRATEGY

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsNewtonRaphsonStrategy :
  public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsNewtonRaphsonStrategy);

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>              BaseType;
    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>                 TConvergenceCriteriaType;
    typedef typename BaseType::TBuilderAndSolverType                          TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType                                              TSchemeType;
    typedef typename BaseType::DofsArrayType                                          DofsArrayType;
    typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;
    using MotherType::mpScheme;
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpA; //Tangent matrix
    using MotherType::mpb; //Residual vector of iteration i
    using MotherType::mpDx; //Delta x of iteration i
    using MotherType::mMaxIterationNumber;
    using MotherType::mInitializeWasPerformed;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ///Constructor
    /**
     * @brief Default constructor 
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    GeoMechanicsNewtonRaphsonStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters& rParameters,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                         pScheme,
                                                                                         /*pNewLinearSolver,*/
                                                                                         pNewConvergenceCriteria,
                                                                                         pNewBuilderAndSolver,
                                                                                         MaxIterations,
                                                                                         CalculateReactions,
                                                                                         ReformDofSetAtEachStep,
                                                                                         MoveMeshFlag)
        {
            //only include validation with c++11 since raw_literals do not exist in c++03
            Parameters default_parameters( R"(
                {
                    "min_iteration":    2,
                    "number_cycles":    5,
                    "increase_factor":  2.0,
                    "reduction_factor": 0.5,
                    "end_time": 1.0,
                    "max_piping_iterations": 50,
                    "desired_iterations": 4,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list" : [],
                    "rebuild_level": 2
                }  )" );

            // Validate agains defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);

            mpParameters = &rParameters;

            // Set Load SubModelParts and Variable names
            if (rParameters["loads_sub_model_part_list"].size() > 0)
            {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                if ( mSubModelPartList.size() != mVariableNames.size() )
                    KRATOS_ERROR << "For each SubModelPart there must be a corresponding nodal Variable"
                                 << std::endl;

                for (unsigned int i = 0; i < mVariableNames.size(); ++i) {
                    mSubModelPartList[i] = &( model_part.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~GeoMechanicsNewtonRaphsonStrategy() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters* mpParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    int Check() override
    {
        KRATOS_TRY

        return MotherType::Check();

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateReferenceDofsNorm(DofsArrayType& rDofSet)
    {
        double ReferenceDofsNorm = 0.0;

        int NumThreads = ParallelUtilities::GetNumThreads();
        OpenMPUtils::PartitionVector DofSetPartition;
        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofSetPartition);

        #pragma omp parallel reduction(+:ReferenceDofsNorm)
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofsBegin = rDofSet.begin() + DofSetPartition[k];
            typename DofsArrayType::iterator DofsEnd = rDofSet.begin() + DofSetPartition[k+1];

            for (typename DofsArrayType::iterator itDof = DofsBegin; itDof != DofsEnd; ++itDof)
            {
                if (itDof->IsFree())
                {
                    const double& temp = itDof->GetSolutionStepValue();
                    ReferenceDofsNorm += temp*temp;
                }
            }
        }

        return sqrt(ReferenceDofsNorm);
    }

private:

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos

#endif // KRATOS_GEO_MECHANICS_NEWTON_RAPHSON_STRATEGY  defined
