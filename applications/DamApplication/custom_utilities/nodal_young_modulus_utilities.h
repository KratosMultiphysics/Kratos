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
#include "includes/node.h"
#include "geometries/geometry.h"
#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Dam-side interpolation helper for the historical NODAL_YOUNG_MODULUS.
 * @details The generic SMA/CLA elastic kernels retrieve Young's modulus through
 * Properties::GetValue(YOUNG_MODULUS, ...). Dam's nodal laws instead read the
 * legacy NODAL_YOUNG_MODULUS nodal field and interpolate it at the integration
 * point with the element shape functions:
 *
 *     E_gp = sum_i N_i * NODAL_YOUNG_MODULUS_i
 *
 * This helper is the single source of that interpolation, reused by the thin
 * nodal elastic/thermoelastic Dam adapters when overriding the E-consuming
 * generic seams (CalculatePK2Stress / CalculateElasticMatrix).
 */
class NodalYoungModulusUtilities
{
public:

    using SizeType       = std::size_t;
    using GeometryType   = Geometry<Node>;

    /**
     * @brief Interpolated nodal Young modulus at a Gauss point.
     * @param rGeometry The element geometry (reference node data for
     *                  NODAL_YOUNG_MODULUS).
     * @param rShapeFunctionValues The shape function values (vector at the
     *                             integration point);
     * @return The interpolated Young modulus.
     */
    static double InterpolatedYoungModulus(
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionValues)
    {
        double young_modulus = 0.0;
        const SizeType number_of_nodes = rGeometry.size();
        for (SizeType j = 0; j < number_of_nodes; ++j) {
            young_modulus += rShapeFunctionValues[j] *
                rGeometry[j].FastGetSolutionStepValue(NODAL_YOUNG_MODULUS);
        }
        return young_modulus;
    }

}; // class NodalYoungModulusUtilities

} // namespace Kratos

#endif // KRATOS_NODAL_YOUNG_MODULUS_UTILITIES_H_INCLUDED defined