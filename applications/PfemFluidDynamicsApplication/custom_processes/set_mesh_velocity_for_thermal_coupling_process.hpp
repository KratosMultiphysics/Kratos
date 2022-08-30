//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Massimiliano Zecchetto
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined(SET_MESH_VELOCITY_FOR_THERMAL_COUPLING_PROCESS)
#define  SET_MESH_VELOCITY_FOR_THERMAL_COUPLING_PROCESS

#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

/**
 * @class SetMeshVelocityForThermalCouplingProcess
 * @ingroup PfemFluidDynamicsApplication
 * @brief This method sets the MESH_VELOCITY equal to the nodal VELOCITY
 * @author Massimiliano Zecchetto
*/
class SetMeshVelocityForThermalCouplingProcess : public Process {

public:

    KRATOS_CLASS_POINTER_DEFINITION(SetMeshVelocityForThermalCouplingProcess);

    /// Constructor
    explicit SetMeshVelocityForThermalCouplingProcess(ModelPart& model_part) : rModelPart(model_part) {}

    /// Destructor.
    ~SetMeshVelocityForThermalCouplingProcess() override {}

	void operator()() {
		Execute();
	}

    void Execute() override {
        KRATOS_TRY;

        this->SetMeshVelocity(rModelPart);

        KRATOS_CATCH("");
    }

    void ExecuteInitialize() override {}

    void ExecuteInitializeSolutionStep() override {}

protected:

    ModelPart& rModelPart;

private:

    void SetMeshVelocity(ModelPart& rModelPart) const {
        const auto& it_node_begin = rModelPart.NodesBegin();
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.Nodes().size()); i++) {
            auto it_node = it_node_begin + i;
            noalias(it_node->FastGetSolutionStepValue(MESH_VELOCITY)) = it_node->FastGetSolutionStepValue(VELOCITY);
        }
    }

}; // Class SetMeshVelocityForThermalCouplingProcess

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
	SetMeshVelocityForThermalCouplingProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
	const SetMeshVelocityForThermalCouplingProcess &rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

} // namespace Kratos.

#endif /* SET_MESH_VELOCITY_FOR_THERMAL_COUPLING_PROCESS defined */