//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license:
// kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

#if !defined(KRATOS_UPDATE_PRESSURE_VALUE_PFEM_CONDITIONS_PROCESS)
#define KRATOS_UPDATE_PRESSURE_VALUE_PFEM_CONDITIONS_PROCESS

#include "processes/process.h"
#include "fem_to_dem_application_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"

namespace Kratos {

typedef std::size_t SizeType;

/** 
 * @class UpdatePressureValuePfemConditionsProcess
 * @ingroup FemToDemApplication 
 * @brief Assigns the pressure value according to the nodal PRESSURE at this time step
 * @author Alejandro Cornejo
 */
template <SizeType TDim = 3>
class UpdatePressureValuePfemConditionsProcess : public Process 
{

public:
    /// Pointer definition of UpdatePressureValuePfemConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(UpdatePressureValuePfemConditionsProcess);

    // Constructor
    UpdatePressureValuePfemConditionsProcess(ModelPart& r_model_part);

    // Destructor
    ~UpdatePressureValuePfemConditionsProcess() override = default;

    void operator()() { Execute(); }

    /**
     * @brief Updates the pressure values due to the PFEM fluid
     */
    void Execute() override;


protected:
    // Member Variables
    ModelPart& mrModelPart;

};  // Class

}  // namespace Kratos
#endif /* KRATOS_EXTEND_PRESSURE_PROCESS defined */