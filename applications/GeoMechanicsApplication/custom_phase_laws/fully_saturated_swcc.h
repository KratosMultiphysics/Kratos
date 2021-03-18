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



#if !defined(KRATOS_SOIL_MECHANICS_FULLY_SATURATED_SWCC_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_FULLY_SATURATED_SWCC_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Class for the fully saturated soil water characteristic curve.
 */
class FullySaturatedSWCC : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(FullySaturatedSWCC);

    /**
     * Constructor.
     */
    FullySaturatedSWCC()
    {}

    /**
     * Destructor.
     */
    virtual ~FullySaturatedSWCC()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      FullySaturatedSWCC::Pointer p_clone(new FullySaturatedSWCC());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new FullySaturatedSWCC());
    }

    /**
     * Operations
     */

    /// Get the value of the saturation w.r.t capillary pressure
    double GetValue(const double& capillaryPressure) const final
    {
        return 1.0;
    }

    /// Get the derivative of the saturation w.r.t capillary pressure
    double GetDerivative(const double& capillaryPressure) const final
    {
        return 0.0;
    }

    /// Get the second derivative of the the saturation w.r.t capillary pressure
    double GetSecondDerivative(const double& capillaryPressure) const final
    {
        return 0.0;
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "FullySaturatedSWCC";
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

}; /* Class FullySaturatedSWCC */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_FULLY_SATURATED_SWCC_H_INCLUDED  defined */

