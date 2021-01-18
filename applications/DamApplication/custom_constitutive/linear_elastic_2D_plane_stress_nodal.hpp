//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/linear_elastic_2D_plane_strain_nodal.hpp"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) LinearElastic2DPlaneStressNodal : public LinearElastic2DPlaneStrainNodal
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic2DPlaneStressNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStressNodal();

    // Copy Constructor
    LinearElastic2DPlaneStressNodal (const LinearElastic2DPlaneStressNodal& rOther);

    // Destructor
    virtual ~LinearElastic2DPlaneStressNodal();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, const double &rYoungModulus, const double &rPoissonCoefficient ) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic2DPlaneStrainNodal)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic2DPlaneStrainNodal)
    }

}; // Class LinearElastic2DPlaneStressNodal
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED  defined