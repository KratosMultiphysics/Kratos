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

#if !defined (KRATOS_LINEAR_PLANE_STRESS_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_PLANE_STRESS_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearPlaneStress : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearPlaneStress
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearPlaneStress );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearPlaneStress();

    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    LinearPlaneStress (const LinearPlaneStress& rOther);


    /**
     * Destructor.
     */
    virtual ~LinearPlaneStress();

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
    
    SizeType GetStrainSize(){return 3;};
    
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
        double c1 = E / (1.00 - NU*NU);
        double c2 = c1 * NU;
        double c3 = 0.5*E / (1 + NU);

        C(0,0) = c1;
        C(0,1) = c2;
        C(0,2) = 0.0;
        C(1,0) = c2;
        C(1,1) = c1;
        C(1,2) = 0.0;
        C(2,0) = 0.0;
        C(2,1) = 0.0;
        C(2,2) = c3;
    }

    void CalculateStress(const Vector& StrainVector, Vector& StressVector, const double E, const double NU)
    {
        double c1 = E / (1.00 - NU*NU);
        double c2 = c1 * NU;
        double c3 = 0.5* E / (1 + NU);


        StressVector[0] = c1*StrainVector[0] + c2 * (StrainVector[1])	;
        StressVector[1] = c1*StrainVector[1] + c2 * (StrainVector[0])	;
        StressVector[2] = c3*StrainVector[2];

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearPlaneStress )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearPlaneStress )
    }


}; // Class LinearPlaneStress
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRESS_LAW_H_INCLUDED  defined 
