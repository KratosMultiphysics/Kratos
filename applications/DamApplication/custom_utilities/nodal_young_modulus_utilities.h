//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_NODAL_YOUNG_MODULUS_UTILITIES_H_INCLUDED)
#define  KRATOS_NODAL_YOUNG_MODULUS_UTILITIES_H_INCLUDED

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/variables.h"
#include "includes/database_accessor.h"
#include "geometries/geometry.h"

namespace Kratos
{

/**
 * @brief Dam-side utilities for the historical spatially varying Young's
 * modulus field.
 * @details The nodal field is stored as nodal YOUNG_MODULUS. Accessor-aware
 * standard constitutive laws (e.g. SMA::FlexibleElasticIsotropic3D or the
 * ConstitutiveLawsApplication thermal elastic laws) retrieve it transparently
 * through the standard Kratos DatabaseAccessor:
 *
 *     Properties::GetValue(YOUNG_MODULUS, geometry, N, process_info)
 *         -> DatabaseAccessor(node_historical)
 *         -> E_gp = sum_i N_i * YOUNG_MODULUS_i
 *
 * The mechanical plane-strain/plane-stress SMA bases access YOUNG_MODULUS
 * directly from the Properties and therefore do not query Accessors; the thin
 * Dam 2D compatibility laws reuse InterpolatedYoungModulus() instead.
 */
class NodalYoungModulusUtilities
{
public:

    using SizeType     = std::size_t;
    using GeometryType = Geometry<Node>;

    /**
     * @brief Interpolated nodal Young modulus (E_gp = sum_i N_i * E_i).
     */
    static double InterpolatedYoungModulus(
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionValues)
    {
        double young_modulus = 0.0;
        const SizeType number_of_nodes = rGeometry.size();
        for (SizeType j = 0; j < number_of_nodes; ++j) {
            young_modulus += rShapeFunctionValues[j] *
                rGeometry[j].FastGetSolutionStepValue(YOUNG_MODULUS);
        }
        return young_modulus;
    }

    /**
     * @brief Installs the standard DatabaseAccessor (node_historical) that
     * exposes the nodal YOUNG_MODULUS field through Properties.
     */
    static void InstallDatabaseAccessor(ModelPart& rModelPart)
    {
        for (auto& r_properties : rModelPart.GetMesh(0).Properties()) {
            InstallDatabaseAccessor(r_properties);
        }
    }

    /**
     * @brief Installs the standard DatabaseAccessor (node_historical) on a
     * single Properties, if not already present.
     */
    static void InstallDatabaseAccessor(Properties& rProperties)
    {
        if (!rProperties.HasAccessor(YOUNG_MODULUS)) {
            rProperties.SetAccessor(
                YOUNG_MODULUS,
                Kratos::make_unique<DatabaseAccessor>("node_historical"));
        }
    }

}; // class NodalYoungModulusUtilities

} // namespace Kratos

#endif // KRATOS_NODAL_YOUNG_MODULUS_UTILITIES_H_INCLUDED defined
