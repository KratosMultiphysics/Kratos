//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_LAW_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/linear_elastic_3D_law.hpp"

namespace Kratos
{

class LinearElastic2DPlaneStressLaw : public LinearElastic3DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION( LinearElastic2DPlaneStressLaw );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStressLaw();

    // Copy Constructor
    LinearElastic2DPlaneStressLaw (const LinearElastic2DPlaneStressLaw& rOther);

    // Destructor
    virtual ~LinearElastic2DPlaneStressLaw();

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

}; // Class LinearElastic2DPlaneStressLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_LAW_H_INCLUDED  defined 
