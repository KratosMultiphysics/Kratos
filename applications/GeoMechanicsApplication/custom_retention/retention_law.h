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
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/properties.h"
#include "includes/serializer.h"
#include <optional>

namespace Kratos
{

/**
 * Base class of retention laws.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) RetentionLaw
{
public:
    using GeometryType = Geometry<Node>;

    // Counted pointer of RetentionLaw
    KRATOS_CLASS_POINTER_DEFINITION(RetentionLaw);

    class KRATOS_API(GEO_MECHANICS_APPLICATION) Parameters
    {
        KRATOS_CLASS_POINTER_DEFINITION(Parameters);

        /**
         * Structure "Parameters" to be used by the element to pass the parameters into the retention law *
         */

    public:
        explicit Parameters(const Properties& rMaterialProperties);
        ~Parameters() = default;

        void SetFluidPressure(double FluidPressure);

        [[nodiscard]] double GetFluidPressure() const;

        [[nodiscard]] const Properties& GetMaterialProperties() const;

    private:
        std::optional<double> mFluidPressure;
        const Properties&     mrMaterialProperties;

    }; // class Parameters end

    virtual ~RetentionLaw() = default;

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this retention law
     * @note implementation scheme:
     *      RetentionLaw::Pointer p_clone(new RetentionLaw());
     *      return p_clone;
     */
    [[nodiscard]] virtual Pointer Clone() const = 0;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameters the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(Parameters& rParameters, const Variable<double>& rThisVariable, double& rValue) const;

    static std::vector<double> CalculateRelativePermeabilityValues(const std::vector<Pointer>& rRetentionLawVector,
                                                                   const Properties& rProperties,
                                                                   const std::vector<double>& rFluidPressures);

    virtual double CalculateSaturation(Parameters& rParameters) const = 0;

    virtual double CalculateEffectiveSaturation(Parameters& rParameters) const = 0;

    virtual double CalculateDerivativeOfSaturation(Parameters& rParameters) const = 0;

    virtual double CalculateRelativePermeability(Parameters& rParameters) const = 0;

    virtual double CalculateBishopCoefficient(Parameters& rParameters) const = 0;

    /**
     * This function is designed to be called once to perform all the checks
     * needed on the input provided. Checks can be "expensive" as the function
     * is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rCurrentProcessInfo
     * @return
     */
    virtual int Check(const Properties& rMaterialProperties, const ProcessInfo& rCurrentProcessInfo) = 0;

    static int Check(const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                     const Properties&                         rProperties,
                     const ProcessInfo&                        rCurrentProcessInfo);

    [[nodiscard]] virtual std::string Info() const;

    virtual void PrintInfo(std::ostream& rOStream) const;

    virtual void PrintData(std::ostream& rOStream) const;

private:
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

}; /* Class RetentionLaw */

// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RetentionLaw& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/