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



#if !defined(KRATOS_SOIL_MECHANICS_SIMPLE_RELATIVE_PERMEABILITY_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_SIMPLE_RELATIVE_PERMEABILITY_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"


namespace Kratos
{

/**
 * In this law, the relative permeability is set to be equal to saturation.
 */
class SimpleRelativePermeabilityLaw : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(SimpleRelativePermeabilityLaw);

    /**
     * Constructor.
     */
    SimpleRelativePermeabilityLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~SimpleRelativePermeabilityLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      SimpleRelativePermeabilityLaw::Pointer p_clone(new SimpleRelativePermeabilityLaw());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new SimpleRelativePermeabilityLaw());
    }

    /**
     * Operations
     */

    /// Get the value of the relative permeability w.r.t saturation
    double GetValue(const double& saturation) const final
    {
        return saturation;
    }

    /// Get the derivative of the relative permeability w.r.t saturation
    double GetDerivative(const double& saturation) const final
    {
        return 1.0;
    }

    /// Get the second derivative of the relative permeability w.r.t saturation
    double GetSecondDerivative(const double& saturation) const final
    {
        return 0.0;
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "SimpleRelativePermeabilityLaw";
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

}; /* Class SimpleRelativePermeabilityLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_SIMPLE_RELATIVE_PERMEABILITY_LAW_H_INCLUDED  defined */
