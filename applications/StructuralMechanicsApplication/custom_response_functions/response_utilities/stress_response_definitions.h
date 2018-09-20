//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//

#if !defined(KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED )
#define  KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

enum class TracedStressType
{
    FX,
    FY,
    FZ,
    MX,
    MY,
    MZ,
    FXX,
    FXY,
    FXZ,
    FYX,
    FYY,
    FYZ,
    FZX,
    FZY,
    FZZ,
    MXX,
    MXY,
    MXZ,
    MYX,
    MYY,
    MYZ,
    MZX,
    MZY,
    MZZ,
    PK2
};

enum class StressTreatment
{
    Mean,
    Node,
    GaussPoint
};

namespace StressResponseDefinitions
{

    TracedStressType ConvertStringToTracedStressType(const std::string& Str);

    StressTreatment ConvertStringToStressTreatment(const std::string& Str);

} // namespace StressResponseDefinitions.


/** \brief StressCalculation
 *
 * This class calculates a specific entry of the stress tensor of an element.
 * It calls the CalculateOnIntegrationPoints function of the element and extracts
 * the desired value from.
 */
class StressCalculation
{
public:

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    static void CalculateStressOnNode(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

    static void CalculateStressOnGP(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

private:

    static void CalculateStressOnGPLinearTruss(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

    static void CalculateStressOnGPTruss(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

    static void CalculateStressOnGPShell(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

    static void CalculateStressBeam(Element& rElement,
                            const TracedStressType rTracedStressType,
                            std::vector< array_1d<double, 3 > >& rStressVector,
                            const ProcessInfo& rCurrentProcessInfo,
                            int& rDirection);

    static void CalculateStressOnGPBeam(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

    static void CalculateStressOnNodeBeam(Element& rElement,
                            const TracedStressType rTracedStressType,
                            Vector& rOutput,
                            const ProcessInfo& rCurrentProcessInfo);

};  // class StressCalculation.

}  // namespace Kratos.

#endif // KRATOS_STRESS_RESPONSE_DEFINITIONS_H_INCLUDED  defined


