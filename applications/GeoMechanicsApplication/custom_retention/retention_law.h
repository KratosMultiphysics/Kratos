// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

/* System includes */

/* External includes */

/* Project includes */
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/serializer.h"
#include <optional>

namespace Kratos {

/**
 * Base class of retention laws.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) RetentionLaw {
public:
    /**
     * Type definitions
     * NOTE: geometries are assumed to be of type Node for all problems
     */
    using ProcessInfoType = ProcessInfo;
    using SizeType = std::size_t;
    using GeometryType = Geometry<Node>;

    /**
     * Counted pointer of RetentionLaw
     */
    KRATOS_CLASS_POINTER_DEFINITION(RetentionLaw);

    class Parameters {
        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

        /**
         * Structure "Parameters" to be used by the element to pass the parameters into the retention law *

        * KINEMATIC PARAMETERS:

        *** NOTE: Pointers are used only to point to a certain variable,
        *   no "new" or "malloc" can be used for this Parameters ***

        * MATERIAL PROPERTIES:
        * @param mrMaterialProperties reference to the material's Properties object (input data)

        * PROCESS PROPERTIES:
        * @param mrCurrentProcessInfo reference to current ProcessInfo instance (input data)

        */

    public:
        Parameters(const Properties& rMaterialProperties,
                   const ProcessInfo& rCurrentProcessInfo)
            : mrCurrentProcessInfo(rCurrentProcessInfo),
              mrMaterialProperties(rMaterialProperties){};

        ~Parameters() = default;

        void SetFluidPressure(double FluidPressure)
        {
            mFluidPressure = FluidPressure;
        };

        double GetFluidPressure() const
        {
            KRATOS_ERROR_IF_NOT(mFluidPressure.has_value())
                << "Fluid pressure is not yet set in the retention "
                   "law when trying to retrieve it, aborting.\n";
            return mFluidPressure.value();
        }

        const ProcessInfo& GetProcessInfo() const
        {
            return mrCurrentProcessInfo;
        }

        const Properties& GetMaterialProperties() const
        {
            return mrMaterialProperties;
        }

    private:
        std::optional<double> mFluidPressure;
        const ProcessInfo& mrCurrentProcessInfo;
        const Properties& mrMaterialProperties;

    }; // class Parameters end

    RetentionLaw() = default;

    virtual ~RetentionLaw() = default;

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this retention law
     * @note implementation scheme:
     *      RetentionLaw::Pointer p_clone(new RetentionLaw());
     *      return p_clone;
     */
    virtual RetentionLaw::Pointer Clone() const = 0;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    virtual double& CalculateValue(Parameters& rParameters,
                                   const Variable<double>& rThisVariable,
                                   double& rValue) = 0;

    virtual double CalculateSaturation(Parameters& rParameters) = 0;

    virtual double CalculateEffectiveSaturation(Parameters& rParameters) = 0;

    virtual double CalculateDerivativeOfSaturation(Parameters& rParameters) = 0;

    virtual double CalculateRelativePermeability(Parameters& rParameters) = 0;

    virtual double CalculateBishopCoefficient(Parameters& rParameters) = 0;

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the retention law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rCurrentProcessInfo process info
     */
    virtual void InitializeMaterial(const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues);

    virtual void Initialize(Parameters& rParameters);

    /**
     * to be called at the beginning of each solution step
     * (e.g. from Element::InitializeSolutionStep)
     */
    virtual void InitializeSolutionStep(Parameters& rParameters);

    /**
     * to be called at the end of each solution step
     * (e.g. from Element::FinalizeSolutionStep)
     */
    virtual void FinalizeSolutionStep(Parameters& rParameters);

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    virtual void Finalize(Parameters& rParameters);

    /**
     * This can be used in order to reset all internal variables of the
     * retention law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    virtual void ResetMaterial(const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry,
                               const Vector& rShapeFunctionsValues);

    /**
     * This function is designed to be called once to perform all the checks
     * needed on the input provided. Checks can be "expensive" as the function
     * is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties& rMaterialProperties,
                      const ProcessInfo& rCurrentProcessInfo) = 0;

    /**
     * @brief This method is used to check that two Retention Laws are the same type (references)
     * @param rLHS The first argument
     * @param rRHS The second argument
     */
    inline static bool HasSameType(const RetentionLaw& rLHS, const RetentionLaw& rRHS)
    {
        return (typeid(rLHS) == typeid(rRHS));
    }

    /**
     * @brief This method is used to check that tow Retention Laws are the same type (pointers)
     * @param rLHS The first argument
     * @param rRHS The second argument
     */
    inline static bool HasSameType(const RetentionLaw* rLHS, const RetentionLaw* rRHS)
    {
        return RetentionLaw::HasSameType(*rLHS, *rRHS);
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "RetentionLaw";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "RetentionLaw has no data";
    }

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

}; /* Class RetentionLaw */

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, RetentionLaw& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RetentionLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/
