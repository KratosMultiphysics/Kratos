//
// Last Modified by:    Salva Latorre
//

// Project includes
#include "NvidiaFlex_application.h"
#include "includes/kernel.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/sphere_3d_1.h"
#include "utilities/quaternion.h"

namespace Kratos {

    //KRATOS_CREATE_VARIABLE(double, RIGID_BODY_MASS)
    //KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RIGID_BODY_CENTER_OF_MASS)

    //FLAGS
    //KRATOS_CREATE_LOCAL_FLAG(NvidiaFlexFlags, HAS_ROTATION, 0);

    //ELEMENTS

    KratosNvidiaFlexApplication::KratosNvidiaFlexApplication() : KratosApplication("NvidiaFlexApplication"){}

    void KratosNvidiaFlexApplication::Register() {
        // Calling base class register to register Kratos components

        KratosApplication::Register();

        KRATOS_INFO("NvidiaFlex") << std::endl;
        KRATOS_INFO("NvidiaFlex") << "     KRATOS NVIDIA FLEX APPLICATION " << std::endl;
        KRATOS_INFO("NvidiaFlex") << std::endl;
        KRATOS_INFO("NvidiaFlex") << "Importing NvidiaFlexApplication... ";

        KRATOS_INFO("") << "( compiled in mode \"" << Kernel::BuildType() << "\" )";

        //KRATOS_REGISTER_VARIABLE(CONTINUUM_INI_NEIGHBOUR_ELEMENTS)
        //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINEAR_VELOCITY)
        
        // ELEMENTS
        //KRATOS_REGISTER_ELEMENT("CylinderParticle2D", mCylinderParticle2D)
        
        //KRATOS_REGISTER_CONDITION("MAPcond", mMapCon3D3N)
        
        // SERIALIZER
        //Serializer::Register("PropertiesProxy", PropertiesProxy());

        KRATOS_INFO("") << " done." << std::endl;
    }
}  // namespace Kratos
