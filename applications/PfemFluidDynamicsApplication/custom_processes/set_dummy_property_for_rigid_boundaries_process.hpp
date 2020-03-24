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

#if !defined(SET_DUMMY_PROPERTY_FOR_RIGID_ELEMENTS_PROCESS)
#define  SET_DUMMY_PROPERTY_FOR_RIGID_ELEMENTS_PROCESS

#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

class SetDummyPropertyForRigidElementsProcess : public Process {

public:

    KRATOS_CLASS_POINTER_DEFINITION(SetDummyPropertyForRigidElementsProcess);

    /// Constructor
    SetDummyPropertyForRigidElementsProcess(ModelPart &model_part,
                                            unsigned int dummy_property_id)
        : rModelPart(model_part) {
      rDummyPropertyId = dummy_property_id;
    }

    /// Destructor.
    ~SetDummyPropertyForRigidElementsProcess() override {}

	void operator()() {
		Execute();
	}

    void Execute() override {

        KRATOS_TRY;

        Properties::Pointer rDummyProperty = rModelPart.pGetProperties(rDummyPropertyId);
        const auto& it_elem_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = it_elem_begin + i;
            it_elem->SetProperties(rDummyProperty);
        }

        KRATOS_CATCH("");

    }

    void ExecuteInitialize() override {}

    void ExecuteInitializeSolutionStep() override {}

protected:

    ModelPart& rModelPart;
    unsigned int rDummyPropertyId;

private:

}; // Class SetDummyPropertyForRigidElementsProcess

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
	SetDummyPropertyForRigidElementsProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
	const SetDummyPropertyForRigidElementsProcess &rThis) {
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

} // namespace Kratos.

#endif /* SET_DUMMY_PROPERTY_FOR_RIGID_ELEMENTS_PROCESS defined */