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



#if !defined(KRATOS_SOIL_MECHANICS_UNIVARIATE_PHASE_LAW_H_INCLUDED )
#define  KRATOS_SOIL_MECHANICS_UNIVARIATE_PHASE_LAW_H_INCLUDED

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/variables.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "containers/flags.h"


namespace Kratos
{

/**
 * Asbtract class for all the phase law depending on single variable.
 */
class UnivariatePhaseLaw : public Flags
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(UnivariatePhaseLaw);

    /**
     * Constructor.
     */
    UnivariatePhaseLaw();

    /**
     * Destructor.
     */
    virtual ~UnivariatePhaseLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     * NOTE: implementation scheme:
     *      UnivariatePhaseLaw::Pointer p_clone(new UnivariatePhaseLaw());
     *      return p_clone;
     */
    virtual UnivariatePhaseLaw::Pointer Clone() const;

    /**
     * Operations
     */
    virtual bool Has( const Variable<int>& rThisVariable );
    virtual bool Has( const Variable<double>& rThisVariable );
    virtual bool Has( const Variable<Vector>& rThisVariable );
    virtual bool Has( const Variable<Matrix>& rThisVariable );

    virtual int& GetValue( const Variable<int>& rThisVariable, int& rValue );
    virtual double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    virtual Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    virtual Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );

    virtual void SetValue( const Variable<int>& rThisVariable, const int& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    virtual void SetValue( const Variable<double>& rThisVariable, const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    virtual void SetValue( const Variable<array_1d<double, 3> >& rThisVariable,
                   const array_1d<double, 3>& rValue, const ProcessInfo& rCurrentProcessInfo );
    virtual void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    virtual void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );

    /// Get the value of the function
    virtual double GetValue(const double& v) const;

    /// Get the derivative of the function
    virtual double GetDerivative(const double& v) const;

    /// Get the second derivative of the function
    virtual double GetSecondDerivative(const double& v) const;

    /**
     * Turn back information as a string.
     */
    virtual std::string Info() const
    {
        return "UnivariatePhaseLaw";
    }

    /**
     * Print information about this object.
     */
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /**
     * Print object's data.
     */
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

private:

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
    }

    ///@}

}; /* Class UnivariatePhaseLaw */

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const UnivariatePhaseLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);
    rOStream << std::endl;

    return rOStream;
}

} /* namespace Kratos.*/

#endif /* KRATOS_SOIL_MECHANICS_UNIVARIATE_PHASE_LAW_H_INCLUDED  defined */
