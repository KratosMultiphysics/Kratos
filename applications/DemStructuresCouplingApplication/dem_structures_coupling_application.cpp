/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "dem_structures_coupling_application.h"

#include "geometries/triangle_3d_3.h"

namespace Kratos {

    // We define the node type
    typedef Node<3> NodeType;

    // STRUCTURAL COUPLING
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DEM_SURFACE_LOAD)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_VELOCITY)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SMOOTHED_STRUCTURAL_VELOCITY)

    KratosDemStructuresCouplingApplication::KratosDemStructuresCouplingApplication()
        : KratosApplication("DemStructuresCouplingApplication"),

        // Adding surface load conditions
        mSurfaceLoadFromDEMCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))) {}

    void KratosDemStructuresCouplingApplication::Register() {
        // Calling base class register to register Kratos components

        KratosApplication::Register();

        // STRUCTURAL COUPLING
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DEM_SURFACE_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_VELOCITY)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BACKUP_LAST_STRUCTURAL_DISPLACEMENT)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SMOOTHED_STRUCTURAL_VELOCITY)

        // Surface loads
        KRATOS_REGISTER_CONDITION("SurfaceLoadFromDEMCondition3D3N", mSurfaceLoadFromDEMCondition3D3N)

        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "     KRATOS DEM STRUCTURES COUPLING APPLICATION " << std::endl;
        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "Importing DemStructuresCouplingApplication... ";
        KRATOS_INFO("") << " done." << std::endl;
    }
}  // namespace Kratos
