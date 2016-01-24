//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/linear_elastic_3D_law.hpp"

namespace Kratos
{

class LinearElastic2DPlaneStrainLaw : public LinearElastic3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LinearElastic2DPlaneStrainLaw );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStrainLaw();

    // Copy Constructor
    LinearElastic2DPlaneStrainLaw (const LinearElastic2DPlaneStrainLaw& rOther);
    
    // Destructor
    virtual ~LinearElastic2DPlaneStrainLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);

    ConstitutiveLaw::Pointer Clone() const;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,const double& rYoungModulus,const double& rPoissonCoefficient );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLaw )
    }
    
}; // Class LinearElastic2DPlaneStrainLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_LAW_H_INCLUDED  defined 
