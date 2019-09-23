//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS)
#define KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS

#include "includes/model_part.h"
#include "processes/process.h"
#include <pybind11/pybind11.h>
// #include <list>
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

typedef std::size_t SizeType;

/** 
 * @class RegeneratePfemPressureConditionsProcess
 * @ingroup FemToDemApplication 
 * @brief Regenerates the pressure conditions for the PFEM coupling
 * @details when several elements are removed this methods generates the line loads
 * in order to adapt to the new geometry
 * @author Alejandro Cornejo
 */
template <SizeType TDim = 3>
class RegeneratePfemPressureConditionsProcess : public Process 
{


public:
    /// Pointer definition of RegeneratePfemPressureConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(RegeneratePfemPressureConditionsProcess);

    // Constructor
    RegeneratePfemPressureConditionsProcess(ModelPart& r_model_part);

    // Destructor
    ~RegeneratePfemPressureConditionsProcess() override = default;

    void operator()() { Execute(); }

    /**
     * @brief Regenerates the prssure load according to ne new boundary
     */
    void Execute() override;

protected:
    // Member Variables
    ModelPart& mrModelPart;

}  // namespace Kratos
#endif /* KRATOS_REGENERATE_PFEM_PRESSURE_CONDITIONS_PROCESS defined */