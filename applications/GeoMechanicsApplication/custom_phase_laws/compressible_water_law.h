/*
==============================================================================
see soil_mechanics_appplication/LICENSE.txt
==============================================================================
*/
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 13 Mar 2021 $
//   Revision:            $Revision: 1.0 $
//
//



#if !defined(KRATOS_SOIL_MECHANICS_COMPRESSIBLE_WATER_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_COMPRESSIBLE_WATER_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Class to compute the water density using exponential rule.
 */
class CompressibleWaterLaw : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(CompressibleWaterLaw);

    /**
     * Constructor.
     */
    CompressibleWaterLaw()
    {}

    CompressibleWaterLaw(const double& DensityWater, const double& BulkWater)
    : mDensityWater(DensityWater), mBulkWater(BulkWater)
    {}

    /**
     * Copy Constructor.
     */
    CompressibleWaterLaw(const CompressibleWaterLaw& rOther)
    : mDensityWater(rOther.mDensityWater), mBulkWater(rOther.mBulkWater)
    {}

    /**
     * Destructor.
     */
    virtual ~CompressibleWaterLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      CompressibleWaterLaw::Pointer p_clone(new CompressibleWaterLaw());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new CompressibleWaterLaw(*this));
    }

    /**
     * Operations
     */

    bool Has( const Variable<double>& rThisVariable ) final
    {
        if (rThisVariable == DENSITY_WATER)
            return true;
        if (rThisVariable == BULK_WATER)
            return true;
        return false;
    }

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final
    {
        if (rThisVariable == DENSITY_WATER)
            rValue = mDensityWater;
        if (rThisVariable == BULK_WATER)
            rValue = mBulkWater;
        return rValue;
    }

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final
    {
        if (rThisVariable == DENSITY_WATER)
            mDensityWater = rValue;
        if (rThisVariable == BULK_WATER)
            mBulkWater = rValue;
    }

    /// Get the value of the water density w.r.t water pressure
    double GetValue(const double& waterPressure) const final
    {
        return mDensityWater*exp(waterPressure/mBulkWater);
    }

    /// Get the derivative of the water density w.r.t water pressure
    double GetDerivative(const double& waterPressure) const final
    {
        return mDensityWater/mBulkWater*exp(waterPressure/mBulkWater);
    }

    /// Get the second derivative of the water density w.r.t water pressure
    double GetSecondDerivative(const double& waterPressure) const final
    {
        return mDensityWater/pow(mBulkWater, 2)*exp(waterPressure/mBulkWater);
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "CompressibleWaterLaw";
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
        rOStream << "Density Water: " << mDensityWater << ", Bulk Water: " << mBulkWater;
    }

private:

    double mDensityWater, mBulkWater;

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

}; /* Class CompressibleWaterLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_COMPRESSIBLE_WATER_LAW_H_INCLUDED  defined */
