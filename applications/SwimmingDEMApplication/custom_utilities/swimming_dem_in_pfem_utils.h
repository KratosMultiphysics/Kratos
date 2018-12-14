#ifndef KRATOS_SWIMMING_DEM_IN_PFEM_UTILS_H
#define KRATOS_SWIMMING_DEM_IN_PFEM_UTILS_H
// /* External includes */

// System includes

// Project includes
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{
class SwimmingDemInPfemUtils
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(SwimmingDemInPfemUtils);

/// Default constructor.

 SwimmingDemInPfemUtils(){}

/// Destructor.

virtual ~SwimmingDemInPfemUtils(){}

//***************************************************************************************************************
//***************************************************************************************************************

void TransferWalls(ModelPart& r_source_model_part, ModelPart& r_destination_model_part)
{
  Properties::Pointer props;
  if(r_destination_model_part.GetMesh(0).Properties().size()) {
    int props_id = r_destination_model_part.PropertiesBegin()->Id();
    props = r_destination_model_part.pGetProperties(props_id);
  }
  else {
    props = Properties::Pointer(new Properties(0));
    props->SetValue(WALL_FRICTION, 0.5773502691896257);
    props->SetValue(WALL_COHESION, 0.0);
    props->SetValue(COMPUTE_WEAR, 0);
    props->SetValue(SEVERITY_OF_WEAR, 0.001);
    props->SetValue(IMPACT_WEAR_SEVERITY, 0.001);
    props->SetValue(BRINELL_HARDNESS, 200.0);
    props->SetValue(YOUNG_MODULUS, 1e20);
    props->SetValue(POISSON_RATIO, 0.25);
  }

  r_destination_model_part.Conditions().Sort();
  int id = 1;
  if (r_destination_model_part.Conditions().size()) {
    id = (r_destination_model_part.ConditionsEnd()-1)->Id() + 1;
  }

  ModelPart::ElementsContainerType& source_elements = r_source_model_part.Elements();

  for (unsigned int i = 0; i < source_elements.size(); i++) {
    ModelPart::ElementsContainerType::iterator it = r_source_model_part.ElementsBegin() + i;
    Geometry< Node<3> >::Pointer p_geometry =  it->pGetGeometry();
    Condition::Pointer cond = Condition::Pointer(new RigidFace3D(id, p_geometry, props));
    r_destination_model_part.AddCondition(cond);
    id++;
  }

}

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const
{
return "";
}

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const
{
}

/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const
{
}


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{


///@}
///@name Protected member r_variables
///@{ template<class T, std::size_t dim>


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


///@}

private:

///@name Static Member r_variables
///@{


///@}
///@name Member r_variables
///@{
RealField mPorosityField;
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

/// Assignment operator.
SwimmingDemInPfemUtils & operator=(SwimmingDemInPfemUtils const& rOther);


///@}

}; // Class SwimmingDemInPfemUtils

}  // namespace Python.

#endif // KRATOS_SWIMMING_DEM_IN_PFEM_UTILS_H
