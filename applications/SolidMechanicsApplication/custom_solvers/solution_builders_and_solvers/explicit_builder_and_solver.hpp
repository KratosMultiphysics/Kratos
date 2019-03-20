//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:           March 2018 $
//
//

#if !defined(KRATOS_EXPLICIT_BUILDER_AND_SOLVER_H_INCLUDED)
#define  KRATOS_EXPLICIT_BUILDER_AND_SOLVER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_builders_and_solvers/solution_builder_and_solver.hpp"

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

/** @brief Explicit Solution Buider and Solver base class
 *  @details This is the explicit builder and solver
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ExplicitBuilderAndSolver : public SolutionBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
 public:


  ///@name Type Definitions
  ///@{

  /// Pointer definition of ExplicitBuilderAndSolver
  KRATOS_CLASS_POINTER_DEFINITION( ExplicitBuilderAndSolver );

  typedef SolutionBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>      BaseType;

  typedef typename BaseType::LocalFlagType                                   LocalFlagType;
  typedef typename BaseType::DofsArrayType                                   DofsArrayType;

  typedef typename BaseType::SystemMatrixType                             SystemMatrixType;
  typedef typename BaseType::SystemVectorType                             SystemVectorType;
  typedef typename BaseType::SystemMatrixPointerType               SystemMatrixPointerType;
  typedef typename BaseType::SystemVectorPointerType               SystemVectorPointerType;
  typedef typename BaseType::LocalSystemVectorType                   LocalSystemVectorType;
  typedef typename BaseType::LocalSystemMatrixType                   LocalSystemMatrixType;

  typedef typename ModelPart::NodesContainerType                        NodesContainerType;
  typedef typename ModelPart::ElementsContainerType                  ElementsContainerType;
  typedef typename ModelPart::ConditionsContainerType              ConditionsContainerType;

  typedef typename BaseType::SchemePointerType                           SchemePointerType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default Constructor.
  ExplicitBuilderAndSolver()
      : BaseType()
  {
  }

  /// Destructor.
  ~ExplicitBuilderAndSolver() override
  {
  }
  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  void BuildLHS(SchemePointerType pScheme,
                ModelPart& rModelPart,
                SystemMatrixType& rA) override
  {
    KRATOS_TRY

    //Set Nodal Mass to zero
    NodesContainerType& pNodes         = rModelPart.Nodes();
    ElementsContainerType& pElements   = rModelPart.Elements();
    ProcessInfo& rCurrentProcessInfo   = rModelPart.GetProcessInfo();

#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
#else
    int number_of_threads = 1;
#endif

    vector<unsigned int> node_partition;
    OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

    vector<unsigned int> element_partition;
    OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

#pragma omp parallel
    {
#pragma omp for

      for(int k=0; k<number_of_threads; k++)
      {
        typename NodesContainerType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
        typename NodesContainerType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

        for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
        {
          double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS);
          nodal_mass = 0.0;
        }
      }

    }

    //Calculate and assemble Mass Matrix on nodes

    bool CalculateLumpedMassMatrix = false;
    if( rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) ){
      CalculateLumpedMassMatrix = rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    }

#pragma omp parallel
    {
      int k = OpenMPUtils::ThisThread();
      typename ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
      typename ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

      for (typename ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)  //MSI: To be parallelized
      {
        Matrix MassMatrix;

        Element::GeometryType& geometry = itElem->GetGeometry();

        (itElem)->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo);

        const unsigned int dimension   = geometry.WorkingSpaceDimension();

        unsigned int index = 0;
        for (unsigned int i = 0; i <geometry.size(); i++)
        {
          index = i*dimension;

          double& mass = geometry(i)->FastGetSolutionStepValue(NODAL_MASS);

          geometry(i)->SetLock();

          if(!CalculateLumpedMassMatrix){
            for (unsigned int j = 0; j <MassMatrix.size2(); j++)
            {
              mass += MassMatrix(index,j);
            }
          }
          else{
            mass += MassMatrix(index,index);
          }

          geometry(i)->UnSetLock();
        }
      }
    }

    rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] = CalculateLumpedMassMatrix;

    KRATOS_CATCH( "" )

        }

  //**************************************************************************
  //**************************************************************************

  void BuildRHS(SchemePointerType pScheme,
                ModelPart& rModelPart,
                SystemVectorType& rb) override
  {
    KRATOS_TRY

    // Compute condition contributions to RHS.
    CalculateAndAddConditionsRHS(pScheme, rModelPart);

    // Compute element contributions to RHS.
    CalculateAndAddElementsRHS(pScheme, rModelPart);


    KRATOS_CATCH( "" )

  }


  void CalculateAndAddConditionsRHS(SchemePointerType pScheme,
                                    ModelPart& rModelPart )
  {

    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo      = rModelPart.GetProcessInfo();
    ConditionsContainerType& rConditions  = rModelPart.Conditions();

#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
#else
    int number_of_threads = 1;
#endif

    vector<unsigned int> condition_partition;
    OpenMPUtils::CreatePartition(number_of_threads, rConditions.size(), condition_partition);


#pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ConditionsContainerType::ptr_iterator it_begin=rConditions.ptr_begin()+condition_partition[k];
      typename ConditionsContainerType::ptr_iterator it_end=rConditions.ptr_begin()+condition_partition[k+1];

      for (typename ConditionsContainerType::ptr_iterator it= it_begin; it!=it_end; ++it)
      {

        LocalSystemVectorType RHS_Condition_Contribution = LocalSystemVectorType(0);

        Element::EquationIdVectorType EquationId; //Dummy

        pScheme->Condition_Calculate_RHS_Contribution(*it, RHS_Condition_Contribution, EquationId, rCurrentProcessInfo);


      }
    }

    KRATOS_CATCH("")
  }

  void CalculateAndAddElementsRHS(SchemePointerType pScheme,
                                  ModelPart& rModelPart)
  {

    KRATOS_TRY

    ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
    ElementsContainerType& pElements        = rModelPart.Elements();

#ifdef _OPENMP
    int number_of_threads = omp_get_max_threads();
#else
    int number_of_threads = 1;
#endif

    vector<unsigned int> element_partition;
    OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);

#pragma omp parallel for
    for(int k=0; k<number_of_threads; k++)
    {
      typename ElementsContainerType::ptr_iterator it_begin=pElements.ptr_begin()+element_partition[k];
      typename ElementsContainerType::ptr_iterator it_end=pElements.ptr_begin()+element_partition[k+1];
      for (typename ElementsContainerType::ptr_iterator it= it_begin; it!=it_end; ++it)
      {

        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);
        Element::EquationIdVectorType EquationId; //Dummy

        pScheme->Calculate_RHS_Contribution(*it, RHS_Contribution, EquationId, rCurrentProcessInfo);

      }
    }

    KRATOS_CATCH("")
  }


  /**
   * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed
  */
  void Clear() override
  {
    BaseType::Clear();
  }

  /**
   * This function is designed to be called once to perform all the checks needed
   * on the input provided. Checks can be "expensive" as the function is designed
   * to catch user's errors.
   * @param rModelPart
   * @return 0 all ok
   */
  int Check(ModelPart& rModelPart) override
  {
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" )
  }

  ///@}
  ///@name Access
  ///@{

  ///@}
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Friends
  ///@{

  ///@}

protected:
  ///@name Protected static Member Variables
  ///@{
  ///@}
  ///@name Protected member Variables
  ///@{
  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
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

private:

  ///@}
  ///@name Static Member Variables
  ///@{
  ///@}
  ///@name Member Variables
  ///@{
  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{
  ///@}
  ///@name Private  Access
  ///@{
  ///@}
  ///@name Private Inquiry
  ///@{
  ///@}
  ///@name Un accessible methods
  ///@{
  ///@}

}; /// Class ExplicitBuilderAndSolver

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_EXPLICIT_BUILDER_AND_SOLVER_H_INCLUDED defined


