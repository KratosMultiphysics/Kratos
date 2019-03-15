//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Suneth Warnakulasuriya
//
//

#if !defined(KRATOS_SCHEME_EXTENSION_INCLUDED)
#define KRATOS_SCHEME_EXTENSION_INCLUDED

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "utilities/indirect_scalar_fwd.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * Interface extensions for elements and conditions to retrieve respective derivatives for time schemes.
 *
 * This class is an interface to get derivatives of variables of interest in an element,
 * which can be used in schemes, therefore scheme can work without knowing the variables
 * used in the element, treating them in a generic way.
 *
 * @class SchemeExtension
 * @ingroup KratosCore
 * @see ResidualBasedBossakVelocityScheme
 */
class SchemeExtension
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(SchemeExtension);

    virtual ~SchemeExtension()
    {
    }

    /**
     * @brief Get the exact double values from an element
     *
     * Retrieves $\underline{x}$ (the non-derivative) vector for a given node,
     * including all the variables of interest which needs to be updated/read
     * in the scheme
     *
     * @param[in]  NodeId                Id of the node which $\underline{x}$ is retrieved
     * @param[in]  Step                  The step which $\underline{x}$ is retrieved
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVector               Vector containing all the $\underline{x}$
     */
    virtual void GetZeroDerivativesVector(std::size_t NodeId,
                                          std::vector<IndirectScalar<double>>& rVector,
                                          std::size_t Step,
                                          ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetZeroDerivativesVector method. Please implement it "
                        "in derrived class.";
    }

    /**
     * @brief Get the exact first derivative values from an element
     *
     * Retrieves $\underline{\dot{x}}$ (the first derivative) vector for a given node,
     * including all the variables of interest which needs to be updated/read
     * in the scheme
     *
     * @param[in]  NodeId                Id of the node which $\underline{\dot{x}}$ is retrieved
     * @param[in]  Step                  The step which $\underline{\dot{x}}$ is retrieved
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVector               Vector containing all the $\underline{\dot{x}}$
     */
    virtual void GetFirstDerivativesVector(std::size_t NodeId,
                                           std::vector<IndirectScalar<double>>& rVector,
                                           std::size_t Step,
                                           ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetFirstDerivativesVector method. Please implement it "
                        "in derrived class.";
    }

    /**
     * @brief Get the exact second derivatives from an element
     *
     * Retrieves $\underline{\ddot{x}}$ (the second derivative) vector for a given node,
     * including all the variables of interest which needs to be updated/read
     * in the scheme
     *
     * @param[in]  NodeId                Id of the node which $\underline{\ddot{x}}$ is retrieved
     * @param[in]  Step                  The step which $\underline{\ddot{x}}$ is retrieved
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVector               Vector containing all the $\underline{\ddot{x}}$
     */
    virtual void GetSecondDerivativesVector(std::size_t NodeId,
                                            std::vector<IndirectScalar<double>>& rVector,
                                            std::size_t Step,
                                            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetSecondDerivativesVector method. Please implement "
                        "it in derrived class.";
    }

    /**
     * @brief Get list of variables of non-derivatives
     *
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVariables            Vector containing all the $\underline{x}$'s variables.
     */
    virtual void GetZeroDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                             ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetZeroDerivativesVariables method. Please implement "
                        "it in derrived class.";
    }

    /**
     * @brief Get list of variables of first derivatives
     *
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVariables            Vector containing all the $\underline{\dot{x}}$'s variables.
     */
    virtual void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                              ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetFirstDerivativesVariables method. Please implement "
                        "it in derrived class.";
    }

    /**
     * @brief Get list of variables of second derivatives
     *
     * @param[in]  rCurrentProcessInfo   Current process info of the model part
     * @param[out] rVariables            Vector containing all the $\underline{\ddot{x}}$'s variables.
     */
    virtual void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables,
                                               ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_ERROR << "Calling SchemeExtension base class "
                        "GetSecondDerivativesVariables method. Please "
                        "implement it in derrived class.";
    }

    virtual std::ostream& Print(std::ostream& os) const
    {
        return os;
    }

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }
};

///@} // Kratos Classes

} // namespace Kratos.

#endif // KRATOS_SCHEME_EXTENSION_INCLUDED  defined
