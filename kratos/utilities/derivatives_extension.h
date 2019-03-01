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

#if !defined(KRATOS_DERIVATIVES_EXTENSION_INCLUDED)
#define KRATOS_DERIVATIVES_EXTENSION_INCLUDED

// System includes
#include <iosfwd>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/process_info.h"
#include "utilities/indirect_scalar_fwd.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class DerivativesExtension
 * @ingroup KratosCore
 * @brief Interface extensions for elements and conditions.
 */
class DerivativesExtension
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DerivativesExtension);

    virtual ~DerivativesExtension()
    {
    }

    virtual void GetZeroDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                          std::size_t Step,
                                          ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void GetFirstDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                           std::size_t Step,
                                           ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void GetSecondDerivativesVector(std::vector<IndirectScalar<double>>& rVector,
                                            std::size_t Step,
                                            ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void GetZeroDerivativesDofsVector(std::vector<Dof<double>::Pointer>& rVector,
                                              ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void GetFirstDerivativesDofsVector(std::vector<Dof<double>::Pointer>& rVector,
                                               ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void GetSecondDerivativesDofsVector(std::vector<Dof<double>::Pointer>& rVector,
                                                ProcessInfo& rCurrentProcessInfo)
    {
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

#endif // KRATOS_DERIVATIVES_EXTENSION_INCLUDED  defined
