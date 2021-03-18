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



#if !defined(KRATOS_SOIL_MECHANICS_IDEAL_GAS_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_IDEAL_GAS_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "univariate_phase_law.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Class to compute the air density according to ideal gas law.
 */
class IdealGasLaw : public UnivariatePhaseLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(IdealGasLaw);

    /**
     * Constructor.
     */
    IdealGasLaw()
    {}

    IdealGasLaw(const double& DensityAir, const double& BulkAir)
    : mDensityAir(DensityAir), mBulkAir(BulkAir)
    {}

    /**
     * Copy Constructor.
     */
    IdealGasLaw(const IdealGasLaw& rOther)
    : mDensityAir(rOther.mDensityAir), mBulkAir(rOther.mBulkAir)
    {}

    /**
     * Destructor.
     */
    virtual ~IdealGasLaw()
    {}

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      IdealGasLaw::Pointer p_clone(new IdealGasLaw());
     *      return p_clone;
     */
    UnivariatePhaseLaw::Pointer Clone() const final
    {
        return UnivariatePhaseLaw::Pointer(new IdealGasLaw(*this));
    }

    /**
     * Operations
     */

    bool Has( const Variable<double>& rThisVariable ) final
    {
        if (rThisVariable == DENSITY_AIR)
            return true;
        if (rThisVariable == BULK_AIR)
            return true;
        return false;
    }

    double& GetValue( const Variable<double>& rThisVariable, double& rValue ) final
    {
        if (rThisVariable == DENSITY_AIR)
            rValue = mDensityAir;
        if (rThisVariable == BULK_AIR)
            rValue = mBulkAir;
        return rValue;
    }

    void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo ) final
    {
        if (rThisVariable == DENSITY_AIR)
            mDensityAir = rValue;
        if (rThisVariable == BULK_AIR)
            mBulkAir = rValue;
    }

    /// Get the value of the air density w.r.t air pressure
    double GetValue(const double& airPressure) const final
    {
        return mDensityAir + mBulkAir * airPressure;
    }

    /// Get the derivative of the air density w.r.t air pressure
    double GetDerivative(const double& airPressure) const final
    {
        return mBulkAir;
    }

    /// Get the second derivative of the air density w.r.t air pressure
    double GetSecondDerivative(const double& airPressure) const final
    {
        return 0.0;
    }

    /**
     * Turn back information as a string.
     */
    std::string Info() const final
    {
        return "IdealGasLaw";
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
        rOStream << "Density Air: " << mDensityAir << ", Bulk Air: " << mBulkAir;
    }

private:

    double mDensityAir, mBulkAir;

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

}; /* Class IdealGasLaw */

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_IDEAL_GAS_LAW_H_INCLUDED  defined */
