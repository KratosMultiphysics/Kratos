// Author: Salva Latorre

#include "flex_wrapper.h"

namespace Kratos {

    FlexWrapper::FlexWrapper() {}

    /// Destructor
    FlexWrapper::~FlexWrapper() {}

    void FlexWrapper::FlexWrapperFunction() {
        
        KRATOS_TRY;
        
        KRATOS_WATCH("TAKI, EH")
        
        KRATOS_CATCH("");
        
        NvFlexSolver* g_solver;
        NvFlexParams g_params;
        NvFlexSetParams(g_solver, &g_params);
        
//        SimBuffers* g_buffers;
//        NvFlexLibrary* g_flexLib;
//        g_buffers = AllocBuffers(g_flexLib);
//        
//        NvFlexUpdateSolver(g_solver, 0.0001f, 1, false);
//        NvFlexSetVelocities(g_solver, g_buffers->velocities.buffer, NULL);
    }

    /// Turn back information as a string.
    std::string FlexWrapper::Info() const {
        
        std::stringstream buffer;
        buffer << "FlexWrapper" ;
        return buffer.str();
    }

    /// Print information about this object.
    void FlexWrapper::PrintInfo(std::ostream& rOStream) const {rOStream << "FlexWrapper";}

    /// Print object's data.
    void FlexWrapper::PrintData(std::ostream& rOStream) const {}
} // namespace Kratos
