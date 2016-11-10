//
//   Project Name:         KratosFluidDynamicsApplication $
//   Last modified by:    $Author:              RZorrilla $
//   Date:                $Date:             October 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED)
#define  KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
/**
 * Defines a bingham non-newtonian constitutive law
 * This material law is defined by the parameters:
 * 1) DYNAMIC_VISCOSITY
 */

class Newtonian2DLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of Newtonian3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(Newtonian2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    Newtonian2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    Newtonian2DLaw (const Newtonian2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //Newtonian3DLaw& operator=(const Newtonian3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~Newtonian2DLaw();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }


}; // Class Newtonian2DLaw
}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_LAW_2D_H_INCLUDED  defined 
