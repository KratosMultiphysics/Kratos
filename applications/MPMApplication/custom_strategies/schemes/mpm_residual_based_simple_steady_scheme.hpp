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
#include "custom_elements/mpm_updated_lagrangian_UP_VMS.hpp"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/material_point_generator_utility.h"


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

  typedef Scheme<TSparseSpace,TDenseSpace>                                      BaseType;

  typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> StaticBaseType;

  typedef MPMResidualBasedSimpleSteadyScheme<TSparseSpace, TDenseSpace>         ClassType;

  typedef typename BaseType::DofsArrayType                                      DofsArrayType;

  typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType                              LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType                              LocalSystemMatrixType;

  typedef Element::GeometryType                                                 GeometryType;

  typedef typename BaseType::Pointer                                     BaseTypePointer;

  ///@}
  ///@name Life Cycle
  ///@{

    /**
    * @brief Constructor. The pseudo static scheme (parameters)
    * @param ThisParameters Dummy parameters
    */

    MPMResidualBasedSimpleSteadyScheme(ModelPart& rGridModelPart)
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>(),
          mGridModelPart(rGridModelPart)
    {

    }

    /** Copy Constructor.
    */

    MPMResidualBasedSimpleSteadyScheme(MPMResidualBasedSimpleSteadyScheme& rOther)
        :StaticBaseType(rOther),
         mGridModelPart(rOther.mGridModelPart)
    {}

     /**
     * @brief Clone method
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new MPMResidualBasedSimpleSteadyScheme(*this) );
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
    const int number_of_nodes = rModelPart.NumberOfNodes();
    const int nelements = static_cast<int>(rModelPart.Elements().size());
    array_1d<double, 3 > output;

    //if orthogonal subscales are computed
    if (CurrentProcessInfo.GetValue(STABILIZATION_TYPE) == 3) {

      KRATOS_INFO_IF("MPMResidualBasedSimpleSteadyScheme", rModelPart.GetCommunicator().MyPID() == 0)
          << "Computing OSS projections" << std::endl;

      // Step1 - Inizialize nodal variables:
      for (int i = 0; i < number_of_nodes; i++) {
        ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
        noalias(it_node->FastGetSolutionStepValue(RESPROJ_DISPL)) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(RESPROJ_PRESS) = 0.0;
        it_node->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
      } 

//-----------------------------------------------------------------------------------------------      
/*std::cout << ".............numberofMP" << rModelPart.Elements().size()<<"\n";
 for (auto& submodelpart : rModelPart.SubModelParts()){
        //Loop over the elements of the ModelPart
          #pragma omp parallel for private(output)
          for (ModelPart::ElementIterator it_elem = submodelpart.ElementsBegin();
                    it_elem != submodelpart.ElementsEnd(); it_elem++) {
              const Geometry< Node >& r_geometry = it_elem ->GetGeometry();
              it_elem->Calculate(RESPROJ_DISPL,output,CurrentProcessInfo);
           }

           
       } */

//-----------------------------------------------------------------------------------------------

      // Step 2 - loop over the mp for computing residuals and interpolate it to the nodes
      //std::cout << ".............numberofMP" << rModelPart.Elements().size()<<"\n";
      ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
      for (int k = 0; k < nelements; k++) {
        auto it_elem = el_begin + k;
        it_elem->Calculate(RESPROJ_DISPL,output,CurrentProcessInfo); //RESPROJ_DISPL
      } 

      // Step 3
      
      for (int i = 0; i < number_of_nodes; i++) {
        ModelPart::NodeIterator it_node = rModelPart.NodesBegin() + i;
        if (it_node->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
          it_node->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
        const double area_inverse = 1.0 / it_node->FastGetSolutionStepValue(NODAL_AREA);
        it_node->FastGetSolutionStepValue(RESPROJ_DISPL) *= area_inverse;
        it_node->FastGetSolutionStepValue(RESPROJ_PRESS) *= area_inverse; 

      }

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
  // MPM Background Grid
    ModelPart& mGridModelPart;

  ///@}

private:
  ///@name Member Variables
  ///@{

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater = TSparseSpace::CreateDofUpdater();

  ///@}

}; /* Class MPMResidualBasedSimpleSteadyScheme */
}  /* namespace Kratos.*/

#endif

