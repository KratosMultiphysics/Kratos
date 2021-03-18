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



#if !defined(KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_WATER_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_WATER_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"


namespace Kratos
{

/**
 * The relative permeability equation for water phase in Liakopolous tes.
 */
class LiakopolousRelativePermeabilityWaterLaw : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(LiakopolousRelativePermeabilityWaterLaw);

    /**
     * Constructor.
     */
    LiakopolousRelativePermeabilityWaterLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~LiakopolousRelativePermeabilityWaterLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      LiakopolousRelativePermeabilityWaterLaw::Pointer p_clone(new LiakopolousRelativePermeabilityWaterLaw());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new LiakopolousRelativePermeabilityWaterLaw());
    }

    /**
     * Operations
     */

    /// Get the value of the relative permeability w.r.t saturation
    double GetValue(const double& saturation) const final
    {
        return 1.0 - 2.207*pow(1.0-saturation, 1.0121);
    }

    /// Get the derivative of the relative permeability w.r.t saturation
    double GetDerivative(const double& saturation) const final
    {
        return 2.207*1.0121*pow(1.0-saturation, 0.0121);
    }

    /// Get the second derivative of the relative permeability w.r.t saturation
    double GetSecondDerivative(const double& saturation) const final
    {
        return -2.207*1.0121*0.0121*pow(1.0-saturation, -0.9879);
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "LiakopolousRelativePermeabilityWaterLaw";
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

}; /* Class LiakopolousRelativePermeabilityWaterLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_WATER_LAW_H_INCLUDED  defined */
