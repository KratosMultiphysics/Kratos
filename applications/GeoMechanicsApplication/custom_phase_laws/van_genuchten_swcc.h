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



#if !defined(KRATOS_SOIL_MECHANICS_VAN_GENUCHTEN_SWCC_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_VAN_GENUCHTEN_SWCC_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Class for the Van Genuchten soil water characteristic curve.
 */
class VanGenuchtenSWCC : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(VanGenuchtenSWCC);

    /**
     * Constructor.
     */
    VanGenuchtenSWCC()
    {}

    VanGenuchtenSWCC(const double& FirstSaturationParam, const double& SecondSaturationParam, const double& AirEntryValue)
    : mFirstSaturationParam(FirstSaturationParam)
    , mSecondSaturationParam(SecondSaturationParam)
    , mAirEntryValue(AirEntryValue)
    {}

    /**
     * Copy Constructor.
     */
    VanGenuchtenSWCC(const VanGenuchtenSWCC& rOther)
    : mFirstSaturationParam(rOther.mFirstSaturationParam)
    , mSecondSaturationParam(rOther.mSecondSaturationParam)
    , mAirEntryValue(rOther.mAirEntryValue)
    {}

    /**
     * Destructor.
     */
    virtual ~VanGenuchtenSWCC()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      VanGenuchtenSWCC::Pointer p_clone(new VanGenuchtenSWCC());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new VanGenuchtenSWCC(*this));
    }

    /**
     * Operations
     */

    bool Has( const Variable<double>& rThisVariable ) final
    {
        if (rThisVariable == FIRST_SATURATION_PARAM
         || rThisVariable == SECOND_SATURATION_PARAM
         || rThisVariable == AIR_ENTRY_VALUE)
            return true;
        return false;
    }

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final
    {
        if (rThisVariable == FIRST_SATURATION_PARAM)
            rValue = mFirstSaturationParam;
        else if (rThisVariable == SECOND_SATURATION_PARAM)
            rValue = mSecondSaturationParam;
        else if (rThisVariable == AIR_ENTRY_VALUE)
            rValue = mAirEntryValue;
        return rValue;
    }

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final
    {
        if (rThisVariable == FIRST_SATURATION_PARAM)
            mFirstSaturationParam = rValue;
        else if (rThisVariable == SECOND_SATURATION_PARAM)
            mSecondSaturationParam = rValue;
        else if (rThisVariable == AIR_ENTRY_VALUE)
            mAirEntryValue = rValue;
    }

    /// Get the value of the saturation w.r.t capillary pressure
    double GetValue(const double& capillaryPressure) const final
    {
        const double& b = mFirstSaturationParam;
        const double& c = mSecondSaturationParam;
        const double& airEntryPressure = mAirEntryValue;

        double saturation;

        if ( capillaryPressure < 0.0 )
        {
            saturation = 1.0;
        }
        else
            saturation = pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c ) );

        if(saturation < 0.0 || saturation > 1.0)
            KRATOS_THROW_ERROR(std::logic_error, "Water saturation is not in range [0.0, 1.0]. It is", saturation)

        return saturation;
    }

    /// Get the derivative of the saturation w.r.t capillary pressure
    double GetDerivative(const double& capillaryPressure) const final
    {
        const double& b = mFirstSaturationParam;
        const double& c = mSecondSaturationParam;
        const double& airEntryPressure = mAirEntryValue;

        double result = 0.0;

        if ( capillaryPressure < 0.0 )
        {
            result = 0.0;
        }
        else
            result = ( -c ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * b *
                 pow(( capillaryPressure / airEntryPressure ), ( b - 1.0 ) ) * 1.0 / airEntryPressure;

        return result;
    }

    /// Get the second derivative of the the saturation w.r.t capillary pressure
    double GetSecondDerivative(const double& capillaryPressure) const final
    {
        const double& b = mFirstSaturationParam;
        const double& c = mSecondSaturationParam;
        const double& airEntryPressure = mAirEntryValue;

        double result = 0.0;

        if ( capillaryPressure < 0.0 )
            result = 0.0;
        else
            result = ( -c ) * b / airEntryPressure * (
                     ( -c - 1.0 ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 2.0 ) )
                     * b / airEntryPressure * pow(( capillaryPressure / airEntryPressure ), ( 2.0 * ( b - 1.0 ) ) )
                     + pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * ( b - 1.0 )
                     * pow(( capillaryPressure / airEntryPressure ), ( b - 2.0 ) ) * 1.0 / airEntryPressure );

        return result;
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "VanGenuchtenSWCC";
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream& rOStream) const final
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream& rOStream) const final
    {
        rOStream << "First Saturation Param: " << mFirstSaturationParam << std::endl;
        rOStream << "Second Saturation Param: " << mSecondSaturationParam << std::endl;
        rOStream << "Air Entry Value: " << mAirEntryValue;
    }

private:

    double mFirstSaturationParam;
    double mSecondSaturationParam;
    double mAirEntryValue;

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

}; /* Class VanGenuchtenSWCC */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_VAN_GENUCHTEN_SWCC_H_INCLUDED  defined */

