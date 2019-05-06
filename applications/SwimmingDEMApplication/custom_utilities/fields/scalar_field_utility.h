#ifndef KRATOS_FIELD_UTILITY_H
#define KRATOS_FIELD_UTILITY_H
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
#include "utilities/openmp_utils.h"
#include "real_field.h"
#include "sets/space_time_set.h"

namespace Kratos
{
class KRATOS_API(SWIMMING_DEM_APPLICATION) RealFieldUtility
{
public:

KRATOS_CLASS_POINTER_DEFINITION(RealFieldUtility);

/// Default constructor.

RealFieldUtility(SpaceTimeSet& rDomain, RealFeild& rField):
    mrDomain(rDomain), mrRealField(rField){}

/// Destructor.

virtual ~RealFieldUtility(){}

//***************************************************************************************************************
//***************************************************************************************************************

void MarkNodesInside(ModelPart& r_model_part, const ProcessInfo& r_current_process_info)
{
    const int nnodes = r_model_part.Nodes().size();
    const double time = r_current_process_info[TIME];
    mIsInArray.resize(nnodes);

    #pragma omp parallel for
    for (int i = 0; i < nnodes; ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        double coor_x = node_it->X();
        double coor_y = node_it->Y();
        double coor_z = node_it->Z();
        bool is_in = mrDomain.IsIn(time, coor_x, coor_y, coor_z);
        node_it->Set(INSIDE, is_in);
        mIsInArray[i] = is_in;
    }
}

//***************************************************************************************************************
//***************************************************************************************************************

double EvaluateFieldAtPoint(const double& time,
                            const array_1d<double, 3>& coor,
                            RealField::Pointer formula)
{
    if (mrDomain.IsIn(time, coor[0], coor[1], coor[2])){
        return(formula->Evaluate(time, coor));
    }

    return(0.0);
}

//***************************************************************************************************************
//***************************************************************************************************************

array_1d<double, 3> EvaluateFieldAtPoint(const double& time,
                                      const array_1d<double, 3>& coor,
                                      RealField::Pointer formula)
{
    if (mrDomain.IsIn(time, coor[0], coor[1], coor[2])){
        double value;
        formula->Evaluate(time, coor, value);
        return(value);
    }

    return(0.0);

}

//***************************************************************************************************************
//***************************************************************************************************************


virtual void ImposeFieldOnNodes(Variable<double>& destination_variable,
                     const double default_value,
                     RealField::Pointer formula,
                     ModelPart& r_model_part,
                     const ProcessInfo& r_current_process_info,
                     const bool recalculate_domain)
{
    const unsigned int nnodes = r_model_part.Nodes().size();
    const double time = r_current_process_info[TIME];

    if (recalculate_domain || nnodes != mIsInArray.size()){
        MarkNodesInside(r_model_part, r_current_process_info);
    }

    #pragma omp parallel for
    for (int i = 0; i < (int)nnodes; ++i){
        ModelPart::NodeIterator node_it = r_model_part.NodesBegin() + i;
        double& destination_value = node_it->FastGetSolutionStepValue(destination_variable);
        destination_value = default_value;

        if (mIsInArray[i]){
            array_1d<double, 3> coor;
            coor[0] = node_it->X();
            coor[1] = node_it->Y();
            coor[2] = node_it->Z();
            destination_value = formula->Evaluate(time, coor);
        }
    }
}


//***************************************************************************************************************
//***************************************************************************************************************

virtual void ImposeFieldOnNodes(ModelPart& r_model_part, const VariablesList& variables_to_be_imposed);

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

SpaceTimeSet& mrDomain;
RealField& mrRealField;
std::vector<bool> mIsInArray;
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
RealFieldUtility & operator=(RealFieldUtility const& rOther);


///@}

}; // Class RealFieldUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
#endif // KRATOS_FIELD_UTILITY_H
