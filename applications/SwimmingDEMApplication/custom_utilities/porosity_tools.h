#ifndef KRATOS_POROSITY_TOOLS_H
#define KRATOS_POROSITY_TOOLS_H
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
#include "utilities/openmp_utils.h"
#include "real_field.h"

namespace Kratos
{
class PorosityUtils
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(PorosityUtils);

/// Default constructor.

PorosityUtils(RealField& porosity_field):mPorosityField(porosity_field){}

/// Destructor.

virtual ~PorosityUtils(){}

//***************************************************************************************************************
//***************************************************************************************************************

void CalculatePorosity(ModelPart& r_model_part, const ProcessInfo& r_current_process_info)
{
    double time = r_current_process_info[TIME];
    const int nnodes = r_model_part.Nodes().size();

    #pragma omp parallel for
    for (int i = 0; i < nnodes; ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        array_1d<double, 3> coor;
        coor[0] = node_it->X();
        coor[1] = node_it->Y();
        coor[2] = node_it->Z();
        node_it->FastGetSolutionStepValue(FLUID_FRACTION) = mPorosityField.Evaluate(time, coor);
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
PorosityUtils & operator=(PorosityUtils const& rOther);


///@}

}; // Class PorosityUtils

}  // namespace Python.

#endif // KRATOS_POROSITY_TOOLS_H
