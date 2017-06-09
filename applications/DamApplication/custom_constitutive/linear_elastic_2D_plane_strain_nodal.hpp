//
//   Project Name:   
//   Last modified by:    $Author:   
//   Date:                $Date:  
//   Revision:            $Revision:  
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"
#include "custom_constitutive/linear_elastic_3D_law_nodal.hpp" 

namespace Kratos
{

class LinearElastic2DPlaneStrainNodal : public LinearElastic3DLawNodal
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic2DPlaneStrainNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStrainNodal();

    // Copy Constructor
    LinearElastic2DPlaneStrainNodal (const LinearElastic2DPlaneStrainNodal& rOther);

    // Destructor
    virtual ~LinearElastic2DPlaneStrainNodal();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    ConstitutiveLaw::Pointer Clone() const;
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
        return 3;
    };

protected:

    // Member Variables
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,Vector& rStrainVector );

    void CalculateLinearElasticMatrix( Matrix& rConstitutiveMatrix, const double &rYoungModulus, const double &rPoissonCoefficient );

    void GetLawFeatures(Features& rFeatures);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    // Serialization
    
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElastic3DLawNodal)
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElastic3DLawNodal)
    }

}; // Class LinearElastic2DPlaneStrainNodal
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED  defined 
