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


#if !defined(KRATOS_MPM_RESIDUAL_BASED_SIMPLE_STEADY_VELOCITY_SCHEME )
#define      KRATOS_MPM_RESIDUAL_BASED_SIMPLE_STEADY_VELOCITY_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/element.h"
#include "custom_elements/mpm_updated_lagrangian_VP_VMS.hpp"
#include "containers/array_1d.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_utilities/material_point_generator_utility.h"


namespace Kratos
{

/**
 * @class MPMResidualBasedSimpleSteadyVelocityScheme
 * @ingroup
 * @brief
 * @details
 */
template<class TSparseSpace,  class TDenseSpace >
class MPMResidualBasedSimpleSteadyVelocityScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{
public:
  ///@name Type Definitions
  ///@{

  KRATOS_CLASS_POINTER_DEFINITION(MPMResidualBasedSimpleSteadyVelocityScheme);

  typedef Scheme<TSparseSpace,TDenseSpace>                                      BaseType;

  typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace> StaticBaseType;

  typedef MPMResidualBasedSimpleSteadyVelocityScheme<TSparseSpace, TDenseSpace> ClassType;

  typedef typename BaseType::DofsArrayType                                      DofsArrayType;

  typedef typename BaseType::TSystemMatrixType                                  TSystemMatrixType;

  typedef typename BaseType::TSystemVectorType                                  TSystemVectorType;

  typedef typename BaseType::LocalSystemVectorType                              LocalSystemVectorType;

  typedef typename BaseType::LocalSystemMatrixType                              LocalSystemMatrixType;

  typedef Element::GeometryType                                                 GeometryType;

  typedef typename BaseType::Pointer                                            BaseTypePointer;

  ///@}
  ///@name Life Cycle
  ///@{

    /**
    * @brief Constructor. The pseudo static scheme (parameters)
    * @param ThisParameters Dummy parameters
    */

    MPMResidualBasedSimpleSteadyVelocityScheme(ModelPart& rGridModelPart)
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>(),
          mGridModelPart(rGridModelPart)
    {

    }

    /** Copy Constructor.
    */

    MPMResidualBasedSimpleSteadyVelocityScheme(MPMResidualBasedSimpleSteadyVelocityScheme& rOther)
        :StaticBaseType(rOther),
         mGridModelPart(rOther.mGridModelPart)
    {}

     /**
     * @brief Clone method
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new MPMResidualBasedSimpleSteadyVelocityScheme(*this) );
    }

    /** Destructor.
    */

    virtual ~MPMResidualBasedSimpleSteadyVelocityScheme() override {}

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

}; /* Class MPMResidualBasedSimpleSteadyVelocityScheme */
}  /* namespace Kratos.*/

#endif

