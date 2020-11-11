//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_INTEGRATION_INFO_H_INCLUDED )
#define  KRATOS_INTEGRATION_INFO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry_data.h"

#include "geometries/nurbs_shape_function_utilities/nurbs_interval.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Integration information for the creation of integration points.
/* Within this class distinct information of integration can be
 * stored and processed.
 */
class IntegrationInfo
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IntegrationPoint
    KRATOS_CLASS_POINTER_DEFINITION(IntegrationInfo);

    typedef typename Point::IndexType SizeType;
    typedef typename Point::IndexType IndexType;

    /// Integration methods implemented specified within enum.
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    ///@}
    ///@name Type Definitions
    ///@{

    enum class QuadratureMethod
    {
        Default,
        GAUSS,
        EXTENDED_GAUSS
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IntegrationInfo()
        : mIntegrationMethod(IntegrationMethod::NumberOfIntegrationMethods)
        , mQuadratureMethod(QuadratureMethod::Default)
        , mNumberOfIntegrationPointsPerSpan(0)
    {
    }

    IntegrationInfo(IntegrationMethod ThisIntegrationMethod)
        : mIntegrationMethod(ThisIntegrationMethod)
    {
        SetNumberOfIntegrationPointsPerSpan(ThisIntegrationMethod);
    }

    IntegrationInfo(SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod = QuadratureMethod::GAUSS)
        : mQuadratureMethod(ThisQuadratureMethod)
        , mNumberOfIntegrationPointsPerSpan(NumberOfIntegrationPointsPerSpan)
    {
        SetIntegrationMethod(NumberOfIntegrationPointsPerSpan, ThisQuadratureMethod);
    }

    ///@}
    ///@name integration rules
    ///@{

    void SetNumberOfIntegrationPointsPerSpan(IntegrationMethod ThisIntegrationMethod)
    {
        switch (ThisIntegrationMethod) {
        case IntegrationMethod::GI_GAUSS_1:
            mNumberOfIntegrationPointsPerSpan = 1;
            mQuadratureMethod = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_2:
            mNumberOfIntegrationPointsPerSpan = 2;
            mQuadratureMethod = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_3:
            mNumberOfIntegrationPointsPerSpan = 3;
            mQuadratureMethod = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_4:
            mNumberOfIntegrationPointsPerSpan = 4;
            mQuadratureMethod = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_5:
            mNumberOfIntegrationPointsPerSpan = 5;
            mQuadratureMethod = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_1:
            mNumberOfIntegrationPointsPerSpan = 1;
            mQuadratureMethod = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_2:
            mNumberOfIntegrationPointsPerSpan = 2;
            mQuadratureMethod = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_3:
            mNumberOfIntegrationPointsPerSpan = 3;
            mQuadratureMethod = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_4:
            mNumberOfIntegrationPointsPerSpan = 4;
            mQuadratureMethod = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_5:
            mNumberOfIntegrationPointsPerSpan = 5;
            mQuadratureMethod = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::NumberOfIntegrationMethods:
            mNumberOfIntegrationPointsPerSpan = 0;
            mQuadratureMethod = QuadratureMethod::Default;
            break;
        }
    }

    SizeType NumberOfIntegrationPointsPerSpan()
    {
        return mNumberOfIntegrationPointsPerSpan;
    }

    void SetIntegrationMethod(SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod = QuadratureMethod::GAUSS)
    {
        mQuadratureMethod = ThisQuadratureMethod;

        switch (NumberOfIntegrationPointsPerSpan) {
        case 1:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                mIntegrationMethod = IntegrationMethod::GI_GAUSS_1;
            } else {
                mIntegrationMethod = IntegrationMethod::GI_EXTENDED_GAUSS_1;
            }
            break;
        case 2:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                mIntegrationMethod = IntegrationMethod::GI_GAUSS_2;
            } else {
                mIntegrationMethod = IntegrationMethod::GI_EXTENDED_GAUSS_2;
            }
            break;
        case 3:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                mIntegrationMethod = IntegrationMethod::GI_GAUSS_3;
            } else {
                mIntegrationMethod = IntegrationMethod::GI_EXTENDED_GAUSS_3;
            }
            break;
        case 4:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                mIntegrationMethod = IntegrationMethod::GI_GAUSS_4;
            } else {
                mIntegrationMethod = IntegrationMethod::GI_EXTENDED_GAUSS_4;
            }
            break;
        case 5:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                mIntegrationMethod = IntegrationMethod::NumberOfIntegrationMethods;
            } else {
                mIntegrationMethod = IntegrationMethod::GI_EXTENDED_GAUSS_5;
            }
            break;
        case 0:
            mIntegrationMethod = IntegrationMethod::NumberOfIntegrationMethods;
            break;
        }
    }

    IntegrationMethod GetIntegrationMethod()
    {
        return mIntegrationMethod;
    }

    ///@}
    ///@name Set geometrical constraints
    ///@{

    void SetSpans(std::vector<double>& rSpans, SizeType DirectionIndex)
    {
        if (mSpansVector.size() <= DirectionIndex) {
            mSpansVector.resize(DirectionIndex + 1);
        }
        mSpansVector[DirectionIndex] = rSpans;
    }

    bool HasSpansInDirection(SizeType DirectionIndex)
    {
        if (mSpansVector.size() > DirectionIndex) {
            if (mSpansVector.size() > 0) {
                return true;
            }
        }
        return false;
    }

    std::vector<double>& GetSpans(SizeType DirectionIndex)
    {
        KRATOS_ERROR_IF(mSpansVector.size() <= DirectionIndex)
            << "Not enough Spans defined at IntegrationInfo." << std::endl;

        return mSpansVector[DirectionIndex];
    }

    static void MergeSpans(std::vector<double>& rResultSpans, const std::vector<double>& rSpans1, const std::vector<double>& rSpans2, double Tolerance = 1e-6) {
        NurbsInterval interval_1(rSpans1[0], rSpans1[rSpans1.size() - 1]);
        NurbsInterval interval_2(rSpans2[0], rSpans2[rSpans2.size() - 1]);

        for (IndexType i = 0; i < rSpans1.size(); ++i) {
            double temp = rSpans1[i];
            interval_2.IsInside(temp);
            rResultSpans.push_back(temp);
        }
        for (IndexType i = 0; i < rSpans2.size(); ++i) {
            double temp = rSpans2[i];
            interval_1.IsInside(temp);
            rResultSpans.push_back(temp);
        }

        SortUnique(rResultSpans, Tolerance);
    }

    static void SortUnique(
        std::vector<double>& rIntersectionParameters,
        const double Tolerance)
    {
        std::sort(std::begin(rIntersectionParameters), std::end(rIntersectionParameters));

        auto last = std::unique(std::begin(rIntersectionParameters), std::end(rIntersectionParameters),
            [=](double a, double b) { return b - a < Tolerance; });

        auto nb_unique = std::distance(std::begin(rIntersectionParameters), last);

        rIntersectionParameters.resize(nb_unique);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << " Integration info. Number of integration points per span: "
            << mNumberOfIntegrationPointsPerSpan;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << " Integration info. Number of integration points per span: "
            << mNumberOfIntegrationPointsPerSpan;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Integration info. Number of integration points per span: "
            << mNumberOfIntegrationPointsPerSpan;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    IntegrationMethod mIntegrationMethod;

    QuadratureMethod mQuadratureMethod;

    SizeType mNumberOfIntegrationPointsPerSpan;

    std::vector<std::vector<double>> mSpansVector;

    ///@}

}; // Class IntegrationPoint

///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::istream& operator >> (std::istream& rIStream, IntegrationInfo& rThis);

/// output stream function
template<std::size_t TDimension, class TDataType, class TWeightType>
inline std::ostream& operator << (std::ostream& rOStream, const IntegrationInfo& rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_INFO_H_INCLUDED  defined


