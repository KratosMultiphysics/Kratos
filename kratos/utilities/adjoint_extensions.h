//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//
//

#if !defined(KRATOS_ADJOINT_EXTENSIONS_INCLUDED)
#define  KRATOS_ADJOINT_EXTENSIONS_INCLUDED

// System includes
#include <iosfwd>
#include <vector>

// Project includes
#include "includes/define.h"
#include "utilities/indirect_scalar_fwd.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class AdjointExtensions
 * @ingroup KratosCore
 * @brief Interface extensions for adjoint elements and conditions.
 */
class AdjointExtensions
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointExtensions);

    virtual void GetFirstDerivativesVector(std::size_t NodeId,
                                           std::vector<IndirectScalar<double>>& rVector,
                                           std::size_t Step)
    {
    }

    virtual void GetSecondDerivativesVector(std::size_t NodeId,
                                            std::vector<IndirectScalar<double>>& rVector,
                                            std::size_t Step)
    {
    }

    virtual void GetAuxiliaryVector(std::size_t NodeId,
                                    std::vector<IndirectScalar<double>>& rVector,
                                    std::size_t Step)
    {
    }

    virtual void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const
    {
    }

    virtual void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const
    {
    }

    virtual void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const
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

}  // namespace Kratos.

#endif // KRATOS_ADJOINT_EXTENSIONS_INCLUDED  defined 
