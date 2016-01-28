//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED

// External includes
#include <cmath>

// Project includes
#include "includes/serializer.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class LinearElastic3DLaw : public ConstitutiveLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic3DLaw);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic3DLaw();

    // Copy Constructor
    LinearElastic3DLaw (const LinearElastic3DLaw& rOther);

    // Destructor
    virtual ~LinearElastic3DLaw();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void GetLawFeatures(Features& rFeatures);
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);
    
    ConstitutiveLaw::Pointer Clone() const;
    
    void InitializeMaterial( const Properties& rMaterialProperties,const GeometryType& rElementGeometry,const Vector& rShapeFunctionsValues );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateMaterialResponseCauchy (Parameters & rValues);
    
    void FinalizeMaterialResponseCauchy (Parameters & rValues);
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    virtual void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix,const double& rYoungModulus,const double& rPoissonCoefficient );

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class LinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 
