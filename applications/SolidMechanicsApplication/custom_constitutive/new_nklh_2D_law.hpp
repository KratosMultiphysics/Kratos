//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:               DAbadias $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_NEW_NKLH_2D_LAW_H_INCLUDED)
#define  KRATOS_NEW_NKLH_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/new_nklh_3D_law.hpp"

namespace Kratos
{
/**
 * Defines a linear isotropic constitutive law in 2D (Plane Strain)
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO 
 * As there are no further parameters the functionality is valid
 * for small and large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) NewNKLH2DLaw : public NewNKLH3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of NewNKLH3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NewNKLH2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NewNKLH2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    NewNKLH2DLaw (const NewNKLH2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElastic3DLaw& operator=(const LinearElastic3DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~NewNKLH2DLaw();

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
        return 4;
    };
    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */


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

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;
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

    virtual void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NewNKLH3DLaw )
    }

    virtual void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NewNKLH3DLaw )
    }


}; // Class LinearElastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_H_INCLUDED  defined 
