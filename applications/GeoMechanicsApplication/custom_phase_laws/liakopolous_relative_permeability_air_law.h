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



#if !defined(KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_AIR_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_AIR_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"


namespace Kratos
{

/**
 * The relative permeability equation for air phase in Liakopolous tes.
 */
class LiakopolousRelativePermeabilityAirLaw : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(LiakopolousRelativePermeabilityAirLaw);

    /**
     * Constructor.
     */
    LiakopolousRelativePermeabilityAirLaw()
    {}

    /**
     * Destructor.
     */
    virtual ~LiakopolousRelativePermeabilityAirLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      LiakopolousRelativePermeabilityAirLaw::Pointer p_clone(new LiakopolousRelativePermeabilityAirLaw());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new LiakopolousRelativePermeabilityAirLaw());
    }

    /**
     * Operations
     */

    /// Get the value of the relative permeability w.r.t saturation
    double GetValue(const double& saturation) const final
    {
        double relSat = ((1.0-saturation)-0.2)/0.8;
        double relPerm = 0.0001+pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));
        return relPerm;
    }

    /// Get the derivative of the relative permeability w.r.t saturation
    double GetDerivative(const double& saturation) const final
    {
        double relSat = ((1.0-saturation)-0.2)/0.8;
        double relPerm_pa = 2*(1.0-relSat)*(-1)/0.8*(1-pow(relSat,5.0/3.0))
            + pow((1.0-relSat),2)*(-1)*5.0/3.0*pow(relSat,2.0/3.0)/0.8;
        return -relPerm_pa;
    }

    /// Get the second derivative of the relative permeability w.r.t saturation
    double GetSecondDerivative(const double& saturation) const final
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "not implemented") // TODO
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "LiakopolousRelativePermeabilityAirLaw";
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

}; /* Class LiakopolousRelativePermeabilityAirLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_LIAKOPOLOUS_RELATIVE_PERMEABILITY_AIR_LAW_H_INCLUDED  defined */
