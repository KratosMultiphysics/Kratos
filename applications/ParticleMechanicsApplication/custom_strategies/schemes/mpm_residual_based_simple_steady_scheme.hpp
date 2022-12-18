//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno
//
//


#if !defined(KRATOS_MPM_RESIDUAL_BASED_SIMPLE_STEADY_SCHEME )
#define      KRATOS_MPM_RESIDUAL_BASED_SIMPLE_STEADY_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"


namespace Kratos
{

/**
 * @class MPMResidualBasedSimpleSteadyScheme
 * @ingroup
 * @brief
 * @details
 */
template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedSimpleSteadyScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{
public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(MPMResidualBasedSimpleSteadyScheme);

  typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> BaseType;

  typedef MPMResidualBasedSimpleSteadyScheme<TSparseSpace, TDenseSpace>         ClassType;

  typedef typename BaseType::DofsArrayType                                      DofsArrayType;

  typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType                              LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType                              LocalSystemMatrixType;

  typedef Element::GeometryType                                                 GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

    /**
    * @brief Constructor. The pseudo static scheme (parameters)
    * @param ThisParameters Dummy parameters
    */

  explicit MPMResidualBasedSimpleSteadyScheme(Parameters ThisParameters)
      : BaseType()
  {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
  }

  /** Default onstructor.
   */

  explicit MPMResidualBasedSimpleSteadyScheme()
      : BaseType()
  {}

    /** Copy Constructor.
    */

   explicit MPMResidualBasedSimpleSteadyScheme(MPMResidualBasedSimpleSteadyScheme& rOther)
        :BaseType(rOther)
    {
    }

    /** Destructor.
    */

  virtual ~MPMResidualBasedSimpleSteadyScheme() override {}


    //***************************************************************************
    //***************************************************************************

  void FinalizeNonLinIteration(ModelPart& rModelPart,
                                       TSystemMatrixType& rA,
                                       TSystemVectorType& rDx,
                                       TSystemVectorType& rb) override
  {

    ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

    //if orthogonal subscales are computed
    if (CurrentProcessInfo.GetValue(STABILIZATION_TYPE) == 3) {

      KRATOS_INFO_IF("MPMResidualBasedSimpleSteadyScheme", rModelPart.GetCommunicator().MyPID() == 0)
          << "Computing OSS projections" << std::endl;

//       const int number_of_nodes = rModelPart.NumberOfNodes();
//
//       #pragma omp parallel for
//       for (int i = 0; i < number_of_nodes; i++) {
//         ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
//         noalias(it_node->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);
//         it_node->FastGetSolutionStepValue(DIVPROJ) = 0.0;
//         it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
//       }
//
//       const int number_of_elements = rModelPart.NumberOfElements();
//       array_1d<double, 3 > output;
//
//       #pragma omp parallel for private(output)
//       for (int i = 0; i < number_of_elements; i++) {
//         ModelPart::ElementIterator it_elem = rModelPart.ElementsBegin() + i;
//         it_elem->Calculate(ADVPROJ,output,CurrentProcessInfo);
//       }
//
// //       rModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
// //       rModelPart.GetCommunicator().AssembleCurrentData(DIVPROJ);
// //       rModelPart.GetCommunicator().AssembleCurrentData(ADVPROJ);
//
//       #pragma omp parallel for
//       for (int i = 0; i < number_of_nodes; i++) {
//         ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
//         if (it_node->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
//           it_node->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
//         const double area_inverse = 1.0 / it_node->FastGetSolutionStepValue(NODAL_AREA);
//         it_node->FastGetSolutionStepValue(ADVPROJ) *= area_inverse;
//         it_node->FastGetSolutionStepValue(DIVPROJ) *= area_inverse;
//       }
    }
  }





    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "static_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedIncrementalUpdateStaticScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
  ///@name Protected Operators
  ///@{

  ///@}

private:
  ///@name Member Variables
  ///@{

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

  ///@}

}; /* Class MPMResidualBasedSimpleSteadyScheme */
}  /* namespace Kratos.*/

#endif

