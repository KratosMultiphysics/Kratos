// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Navaneeth K Narayanan
//

#if !defined(KRATOS_NEWTON_RAPHSON_WITH_HYDROSTATIC_STRATEGY)
#define KRATOS_NEWTON_RAPHSON_WITH_HYDROSTATIC_STRATEGY
// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "custom_utilities/volume_calculation_under_plane_utility.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class NewtonRaphsonWithHydrostaticLoadStrategy
 *
 * @ingroup StrucutralMechanicsApplication
 *
 * @brief inherited class from ResidualBasedNewtonRaphsonStrategy for implementing non-linear hydrostatic loading
 *
 * @details additions for formfinding: update the reference configuration for each element, initialize the elements for formfinding,
 * adaption line search for formfinding, print formfinding output (different nonlinear iterations) // TODO
 *
 * @author Navaneeth K Narayanan
 */

template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class NewtonRaphsonWithHydrostaticLoadStrategy
    : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonWithHydrostaticLoadStrategy);

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef GidIO<> IterationIOType;
    typedef IterationIOType::Pointer IterationIOPointerType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
        * Constructor.
        */

    NewtonRaphsonWithHydrostaticLoadStrategy(
        ModelPart &model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)

        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag)

    {
    }

    // constructor with Builder and Solver
    NewtonRaphsonWithHydrostaticLoadStrategy(
        ModelPart &model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part,
                                                                                       pScheme,
                                                                                       pNewLinearSolver,
                                                                                       pNewConvergenceCriteria,
                                                                                       pNewBuilderAndSolver,
                                                                                       MaxIterations,
                                                                                       CalculateReactions,
                                                                                       ReformDofSetAtEachStep,
                                                                                       MoveMeshFlag)

    {
    }

    /**
        * Destructor.
        */

    ~NewtonRaphsonWithHydrostaticLoadStrategy() override
    {
    }

    void Initialize() override
    {
        KRATOS_TRY;
        BaseType::Initialize();
        double volume;
        double radius;
        Vector centre;
        Vector w;
        double norm_w;
        VolumeCalculationUnderPlaneUtility plane_updater;

        for (ModelPart::PropertiesContainerType::iterator i_prop = BaseType::GetModelPart().PropertiesBegin(); i_prop != BaseType::GetModelPart().PropertiesEnd(); i_prop++)
        {

            if (i_prop->Has(FLUID_VOLUME) && (i_prop->Has(FREE_SURFACE_RADIUS)) && (i_prop->Has(FREE_SURFACE_CENTRE)) && (i_prop->Has(FREE_SURFACE_NORMAL)))

            {
                volume = i_prop->GetValue(FLUID_VOLUME);
                radius = i_prop->GetValue(FREE_SURFACE_RADIUS);
                centre = i_prop->GetValue(FREE_SURFACE_CENTRE);
                w = i_prop->GetValue(FREE_SURFACE_NORMAL);
                norm_w = norm_2(w);
                if (norm_w > std::numeric_limits<double>::epsilon())
                    w = w / norm_w;
                else
                    w *= 0;

                plane_updater.SetPlaneParameters(centre, radius, w);
                mVectorOfVolumes.push_back(volume);
                mVectorOfPlaneUpdaters.push_back(plane_updater);
                mListOfPropertiesId.push_back(i_prop->Id());
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators

    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access

    ///@{

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    void
    UpdateDatabase(
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b,
        const bool MoveMesh) override
    {

        BaseType::UpdateDatabase(A, Dx, b, MoveMesh);

        for (IndexType i = 0; i < mVectorOfPlaneUpdaters.size(); ++i)
        {
            std::cout << "Target Volume :: " << mVectorOfVolumes[i] << std::endl;
            std::cout << "Property Id :: " << mListOfPropertiesId[i] << std::endl;
            mVectorOfPlaneUpdaters[i].UpdatePositionOfPlaneBasedOnTargetVolume(BaseType::GetModelPart(), mVectorOfVolumes[i]);

            BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i])->GetValue(FREE_SURFACE_AREA) = mVectorOfPlaneUpdaters[i].GetIntersectedArea();

            BaseType::GetModelPart().GetMesh(0).Properties().find(mListOfPropertiesId[i])->GetValue(FREE_SURFACE_CENTRE) = mVectorOfPlaneUpdaters[i].GetPlaneCentre();
        }
    }

    ///@name Protected Operations
    ///@}
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    std::vector<VolumeCalculationUnderPlaneUtility> mVectorOfPlaneUpdaters;
    std::vector<double> mVectorOfVolumes;
    std::vector<IndexType> mListOfPropertiesId;

    ///@}
    ///@name Private Operators
    ///@{

    /**
        * Copy constructor.
        */

    //NewtonRaphsonWithHydrostaticLoadStrategy(const NewtonRaphsonWithHydrostaticLoadStrategy &Other){};

    ///@}
}; // namespace Kratos

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* KRATOS_FORMFINDING_UPDATED_REFERENCE_STRATEGY defined */
