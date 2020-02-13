//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Guillermo Casas
//                  Ignasi de Pouplana
//


#ifndef KRATOS_POST_PROCESS_UTILITIES
#define KRATOS_POST_PROCESS_UTILITIES

// System includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"

#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"

// Application includes
#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{

class PostProcessUtilities {

public:

KRATOS_CLASS_POINTER_DEFINITION(PostProcessUtilities);

/// Default constructor.

PostProcessUtilities(ModelPart& rModelPart);

/// Destructor.

virtual ~PostProcessUtilities();

//***************************************************************************************************************
//***************************************************************************************************************

void GetStickyStatus(pybind11::list& is_sticky_list);

void GetInitialContinuumBonds(pybind11::list& initial_continuum_bonds_list);

void GetCurrentContinuumBonds(pybind11::list& current_continuum_bonds_list);

//***************************************************************************************************************
//***************************************************************************************************************

///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a stemplate<class T, std::size_t dim> tring.

virtual std::string Info() const;

/// Print information about this object.

virtual void PrintInfo(std::ostream& rOStream) const;
/// Print object's data.

virtual void PrintData(std::ostream& rOStream) const;


///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member r_variables
///@{

    ModelPart& mrModelPart;

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

void Clear(pybind11::list& my_list);

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
PostProcessUtilities & operator=(PostProcessUtilities const& rOther){};


///@}

}; // Class PostProcessUtilities

}  // namespace Python.

#endif // KRATOS_STRESS_FAILURE_CHECK_UTILITIES
