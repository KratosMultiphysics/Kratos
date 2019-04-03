//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/thermal_linear_elastic_3D_law_nodal.hpp"

namespace Kratos
{

class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic2DPlaneStrainNodal : public ThermalLinearElastic3DLawNodal
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStrainNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStrainNodal();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStrainNodal (const ThermalLinearElastic2DPlaneStrainNodal& rOther);

    // Destructor
    virtual ~ThermalLinearElastic2DPlaneStrainNodal();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
    };

protected:

    // Member Variables

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures) override;

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, const double &rYoungModulus, const double &rPoissonCoefficient ) override;

    void CalculateThermalStrain( Vector& rThermalStrainVector, const MaterialResponseVariables & rElasticVariables, double & rTemperature, double & rNodalReferenceTemperature) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalLinearElastic3DLawNodal)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalLinearElastic3DLawNodal)
    }

}; // Class ThermalLinearElastic2DPlaneStrainNodal
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED  defined
