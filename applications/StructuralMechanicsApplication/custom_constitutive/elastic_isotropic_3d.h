// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined (KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED)
#define  KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ElasticIsotropic3D : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of ElasticIsotropic3D
     */

    KRATOS_CLASS_POINTER_DEFINITION( ElasticIsotropic3D );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    ElasticIsotropic3D();

    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    ElasticIsotropic3D (const ElasticIsotropic3D& rOther);


    /**
     * Destructor.
     */
    virtual ~ElasticIsotropic3D();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);
    
    SizeType GetStrainSize(){return 6;};
    
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    void CalculateMaterialResponseKirchhoff (Parameters & rValues);

    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);



protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}


private:


    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{
    void CalculateElasticMatrix(Matrix& C, const double E, const double NU)
    {
        const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        const double c2 = c1 * ( 1 - NU );
        const double c3 = c1 * NU;
        const double c4 = c1 * 0.5 * ( 1 - 2 * NU );
        
        C( 0, 0 ) = c2;
        C( 0, 1 ) = c3;
        C( 0, 2 ) = c3;
        C( 0, 3 ) = 0.0;
        C( 0, 4 ) = 0.0;
        C( 0, 5 ) = 0.0;
        C( 1, 0 ) = c3;
        C( 1, 1 ) = c2;
        C( 1, 2 ) = c3;
        C( 1, 3 ) = 0.0;
        C( 1, 4 ) = 0.0;
        C( 1, 5 ) = 0.0;
        C( 2, 0 ) = c3;
        C( 2, 1 ) = c3;
        C( 2, 2 ) = c2;
        C( 2, 3 ) = 0.0;
        C( 2, 4 ) = 0.0;
        C( 2, 5 ) = 0.0;
        C( 3, 0 ) = 0.0;
        C( 3, 1 ) = 0.0;
        C( 3, 2 ) = 0.0;
        C( 3, 3 ) = c4;
        C( 3, 4 ) = 0.0;
        C( 3, 5 ) = 0.0;
        C( 4, 0 ) = 0.0;
        C( 4, 1 ) = 0.0;
        C( 4, 2 ) = 0.0;
        C( 4, 3 ) = 0.0;
        C( 4, 4 ) = c4;
        C( 4, 5 ) = 0.0;
        C( 5, 0 ) = 0.0;
        C( 5, 1 ) = 0.0;
        C( 5, 2 ) = 0.0;
        C( 5, 3 ) = 0.0;
        C( 5, 4 ) = 0.0;
        C( 5, 5 ) = c4;
    }

    void CalculateStress(const Vector& StrainVector, Vector& StressVector, const double E, const double NU)
    {
        const double c1 = E / (( 1.00 + NU ) * ( 1 - 2 * NU ) );
        const double c2 = c1 * ( 1 - NU );
        const double c3 = c1 * NU;
        const double c4 = c1 * 0.5 * ( 1 - 2 * NU );

        StressVector[0] = c2*StrainVector[0] + c3 * (StrainVector[1] + StrainVector[2])	;
        StressVector[1] = c2*StrainVector[1] + c3 * (StrainVector[0] + StrainVector[2])	;
        StressVector[2] = c2*StrainVector[2] + c3 * (StrainVector[0] + StrainVector[1])	;
        StressVector[3] = c4*StrainVector[3];
        StressVector[4] = c4*StrainVector[4];
        StressVector[5] = c4*StrainVector[5];
    }


    ///@}
    ///@name Private Operations
    ///@{
    ///@}


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticIsotropic3D )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticIsotropic3D )
    }


}; // Class ElasticIsotropic3D
}  // namespace Kratos.
#endif // KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED  defined 
