//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HERSCHEL_BULKEY_LAW_3D_H_INCLUDED)
#define  KRATOS_HERSCHEL_BULKEY_LAW_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
/**
 * Defines a Herschel-Bulkey non-newtonian constitutive law
 * This material law is defined by the parameters:
 * 1) YIELD_STRESS
 * 2) REGULARIZATION_COEFFICIENT
 * 3) POWER_LAW_K
 * 4) POWER_LAW_N
 */

class HerschelBulkey3DLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HerschelBulkey3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HerschelBulkey3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HerschelBulkey3DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    HerschelBulkey3DLaw (const HerschelBulkey3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HerschelBulkey3DLaw& operator=(const HerschelBulkey3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HerschelBulkey3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */
    virtual void CalculateMaterialResponseCauchy (Parameters& rValues);

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);

    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HerschelBulkey3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HerschelBulkey3DLaw )
    }


}; // Class HerschelBulkey3DLaw
}  // namespace Kratos.
#endif // KRATOS_HERSCHEL_BULKEY_LAW_3D_H_INCLUDED  defined 
