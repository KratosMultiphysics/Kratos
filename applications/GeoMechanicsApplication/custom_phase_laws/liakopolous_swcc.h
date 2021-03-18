/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2 Mar 2021 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_SWCC_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_SWCC_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"


namespace Kratos
{

/**
 * Class for the soil water characteristic curve used in Liakopolous test.
 */
class LiakopolousSWCC : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(LiakopolousSWCC);

    /**
     * Constructor.
     */
    LiakopolousSWCC()
    {}

    /**
     * Destructor.
     */
    virtual ~LiakopolousSWCC()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      LiakopolousSWCC::Pointer p_clone(new LiakopolousSWCC());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new LiakopolousSWCC());
    }

    /**
     * Operations
     */

    /// Get the value of the saturation w.r.t capillary pressure
    double GetValue(const double& capillaryPressure) const final
    {
        return 1.0-1.9722*1e-11*pow(capillaryPressure, 2.4279);
    }

    /// Get the derivative of the saturation w.r.t capillary pressure
    double GetDerivative(const double& capillaryPressure) const final
    {
        return -1.9722*2.4279*1e-11*pow(capillaryPressure, 1.4279);
    }

    /// Get the second derivative of the saturation w.r.t capillary pressure
    double GetSecondDerivative(const double& capillaryPressure) const final
    {
        return -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure, 0.4279);
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "LiakopolousSWCC";
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const final
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    }

    void load(Serializer& rSerializer) final
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }

    ///@}

}; /* Class LiakopolousSWCC */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_SWCC_H_INCLUDED  defined */

