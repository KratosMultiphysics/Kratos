/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

#ifndef KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
#define KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
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
#include "../../DEM_application/custom_conditions/RigidFace.h"
#include "../../DEM_application/DEM_application_variables.h"


namespace Kratos
{
class DemStructuresCouplingUtilities
{
public:
typedef ModelPart::NodesContainerType::ContainerType::iterator NodesIteratorType;

KRATOS_CLASS_POINTER_DEFINITION(DemStructuresCouplingUtilities);

/// Default constructor.

 DemStructuresCouplingUtilities(){}

/// Destructor.

virtual ~DemStructuresCouplingUtilities(){}

//***************************************************************************************************************
//***************************************************************************************************************

void TransferStructuresSkinToDem(ModelPart& r_source_model_part, ModelPart& r_destination_model_part, Properties::Pointer props) {
    
    std::string error = CheckProvidedProperties(props);
    
    if (error != "all_ok") KRATOS_ERROR << "The Dem Walls ModelPart has no valid Properties. Missing " << error << " . Exiting." << std::endl;

    r_destination_model_part.Conditions().Sort();
    int id = 1;
    
    if (r_destination_model_part.Conditions().size()) id = (r_destination_model_part.ConditionsEnd()-1)->Id() + 1;

    ModelPart::ConditionsContainerType& source_conditions = r_source_model_part.Conditions();

    // Adding conditions
    for (unsigned int i = 0; i < source_conditions.size(); i++) {
        ModelPart::ConditionsContainerType::iterator it = r_source_model_part.ConditionsBegin() + i;
        Geometry< Node<3> >::Pointer p_geometry =  it->pGetGeometry();
        Condition::Pointer cond = Condition::Pointer(new RigidFace3D(id, p_geometry, props));
        r_destination_model_part.AddCondition(cond);
        id++;
    }
        
    // Adding nodes
    r_destination_model_part.AddNodes(r_source_model_part.NodesBegin(), r_source_model_part.NodesEnd());
}

std::string CheckProvidedProperties(Properties::Pointer props) {
    std::vector<Variable<double> > list_of_variables_double_to_check = {WALL_FRICTION, WALL_COHESION, SEVERITY_OF_WEAR, IMPACT_WEAR_SEVERITY, BRINELL_HARDNESS, YOUNG_MODULUS, POISSON_RATIO};
    std::vector<Variable<bool> > list_of_variables_bool_to_check = {COMPUTE_WEAR};
    for (int i=0; i<(int)list_of_variables_double_to_check.size(); i++) {
        if(!props->Has(list_of_variables_double_to_check[i])) return list_of_variables_double_to_check[i].Name();
    }
    for (int i=0; i<(int)list_of_variables_bool_to_check.size(); i++) {
        if(!props->Has(list_of_variables_bool_to_check[i])) return list_of_variables_bool_to_check[i].Name();
    }
    return "all_ok";
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
DemStructuresCouplingUtilities & operator=(DemStructuresCouplingUtilities const& rOther);


///@}

}; // Class DemStructuresCouplingUtilities

}  // namespace Python.

#endif // KRATOS_STRUCTURES_DEM_COUPLING_UTILITIES_H
