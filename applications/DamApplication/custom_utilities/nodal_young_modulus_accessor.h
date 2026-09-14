//
//   Project Name:        $ProjectName:            $
//   Last modified by:    $Author:             $
//   Date:                $Date:               $
//   Revision:            $Revision:           $
//

#if !defined(KRATOS_NODAL_YOUNG_MODULUS_ACCESSOR_H_INCLUDED)
#define KRATOS_NODAL_YOUNG_MODULUS_ACCESSOR_H_INCLUDED

// Project includes
#include "includes/accessor.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "includes/variables.h"
#include "geometries/geometry.h"
#include "includes/serializer.h"
#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Provides Young's modulus as the shape-function interpolation of the
 * historical nodal NODAL_YOUNG_MODULUS field:
 *
 *     E_gp = sum_i N_i * NODAL_YOUNG_MODULUS_i
 *
 * @details The standard Kratos DatabaseAccessor interpolates the requested
 * variable itself, so it cannot expose NODAL_YOUNG_MODULUS as YOUNG_MODULUS.
 * This small Dam-side Accessor fills exactly that gap: once installed on a
 * Properties as the accessor of YOUNG_MODULUS, any accessor-aware standard
 * constitutive law (e.g. SMA::FlexibleElasticIsotropic3D or the
 * ConstitutiveLawsApplication thermal elastic laws) transparently retrieves the
 * spatially varying Young's modulus through
 * Properties::GetValue(YOUNG_MODULUS, geometry, N, process_info).
 *
 * Other variables fall back to the regular Properties storage, so installing
 * this accessor only affects YOUNG_MODULUS retrieval.
 */
class KRATOS_API(DAM_APPLICATION) NodalYoungModulusAccessor : public Accessor
{
public:

    ///@name Type Definitions
    ///@{

    using BaseType = Accessor;
    using GeometryType = Geometry<Node>;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(NodalYoungModulusAccessor);

    ///@}
    ///@name Life Cycle
    ///@{

    NodalYoungModulusAccessor() = default;

    ~NodalYoungModulusAccessor() override = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Interpolates NODAL_YOUNG_MODULUS at the integration point.
     */
    double GetValue(
        const Variable<double>& rVariable,
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionVector,
        const ProcessInfo& rProcessInfo
        ) const override
    {
        double young_modulus = 0.0;
        for (SizeType i = 0; i < rShapeFunctionVector.size(); ++i) {
            young_modulus += rShapeFunctionVector[i] *
                rGeometry[i].FastGetSolutionStepValue(NODAL_YOUNG_MODULUS);
        }
        return young_modulus;
    }

    Accessor::UniquePointer Clone() const override
    {
        return Kratos::make_unique<NodalYoungModulusAccessor>(*this);
    }

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "NodalYoungModulusAccessor";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NodalYoungModulusAccessor";
    }

    ///@}

private:

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
    }

    void load(Serializer& rSerializer) override
    {
    }

    ///@}
}; // class NodalYoungModulusAccessor

} // namespace Kratos

#endif // KRATOS_NODAL_YOUNG_MODULUS_ACCESSOR_H_INCLUDED defined
