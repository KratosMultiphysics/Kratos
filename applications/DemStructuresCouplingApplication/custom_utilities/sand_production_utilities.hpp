//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Salva Latorre
//                  Ignasi de Pouplana
//


#ifndef KRATOS_SAND_PRODUCTION_UTILITIES
#define KRATOS_SAND_PRODUCTION_UTILITIES

// System includes
#include <pybind11/pybind11.h>

// Project includes

#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/discrete_particle_configure.h"

#include "dem_structures_coupling_application_variables.h"


namespace Kratos
{

class SandProductionUtilities {

public:

static const std::size_t space_dim                  = 3; ///WARNING: generalize to 2d.
typedef DiscreteParticleConfigure<space_dim>        Configure;

KRATOS_CLASS_POINTER_DEFINITION(SandProductionUtilities);

/// Default constructor.

SandProductionUtilities();

/// Destructor.

virtual ~SandProductionUtilities();

//***************************************************************************************************************
//***************************************************************************************************************

void MarkSandProductionParticlesForErasing(ModelPart& r_model_part);

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
SandProductionUtilities & operator=(SandProductionUtilities const& rOther){};


///@}

}; // Class SandProductionUtilities

}  // namespace Python.

#endif // KRATOS_STRESS_FAILURE_CHECK_UTILITIES
