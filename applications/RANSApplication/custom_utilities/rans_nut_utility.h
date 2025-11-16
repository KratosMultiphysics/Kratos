//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_NUT_UTILITY_H_INCLUDED)
#define KRATOS_RANS_NUT_UTILITY_H_INCLUDED

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(RANS_APPLICATION) RansNutUtility
{
public:
    ///@name Type Definitions

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using ElementType = ModelPart::ElementType;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using TLSType = std::tuple<Vector, Matrix, ShapeFunctionDerivativesArrayType>;

    /// Pointer definition of RansNutUtility
    KRATOS_CLASS_POINTER_DEFINITION(RansNutUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    RansNutUtility(
        ModelPart& rModelPart,
        const double RelativeTolerance,
        const double AbsoluteTolerance,
        const int EchoLevel);

    /**
     * Destructor
     */
    ~RansNutUtility()
    {
        mElementData.clear();
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize();

    void InitializeCalculation();

    bool CheckConvergence() const;

    void UpdateTurbulentViscosity();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "RansNutUtility";
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RansNutUtility";
    }

    ///@}

private:
    ///@name Private Members
    ///@{

    ModelPart& mrModelPart;
    std::vector<double> mElementData;

    const int mEchoLevel;
    const double mRelativeTolerance;
    const double mAbsoluteTolerance;

    ///@}
    ///@name Private Operations
    ///@{

    double CalculateTurbulentViscosity(
        ElementType& rElement,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo) const;

    ///@}

}; // Class RansNutUtility

///@}
///@name Input and output
///@{

inline std::istream& operator>>(std::istream& rIStream, RansNutUtility& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansNutUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_UTILITY_H_INCLUDED defined
