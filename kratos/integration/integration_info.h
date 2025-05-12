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
#include "containers/data_value_container.h"
#include "containers/flags.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Integration information for the creation of integration points.
/* Within this class distinct information of integration can be
 * stored and processed.
 */
class KRATOS_API(KRATOS_CORE) IntegrationInfo : public Flags
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
    ///@name Local Flags
    ///@{

    KRATOS_DEFINE_LOCAL_FLAG(DO_NOT_CREATE_TESSELLATION_ON_SLAVE);

    ///@}
    ///@name Type Definitions
    ///@{

    enum class QuadratureMethod
    {
        Default,
        GAUSS,
        EXTENDED_GAUSS,
        LOBATTO,
        GRID
    };

    ///@}
    ///@name Life Cycle
    ///@{

    IntegrationInfo(SizeType LocalSpaceDimension,
        IntegrationMethod ThisIntegrationMethod);

    IntegrationInfo(SizeType LocalSpaceDimension,
        SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod = QuadratureMethod::GAUSS);

    IntegrationInfo(
        const std::vector<SizeType>& NumberOfIntegrationPointsPerSpanVector,
        const std::vector<QuadratureMethod>& ThisQuadratureMethodVector);

    ///@}
    ///@name Dimension
    ///@{

    SizeType LocalSpaceDimension()
    {
        return mNumberOfIntegrationPointsPerSpanVector.size();
    }

    ///@}
    ///@name integration rules
    ///@{

    void SetIntegrationMethod(
        IndexType DimensionIndex,
        IntegrationMethod ThisIntegrationMethod);

    SizeType GetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex) const;

    void SetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex,
        SizeType NumberOfIntegrationPointsPerSpan);

    QuadratureMethod GetQuadratureMethod(IndexType DimensionIndex) const
    {
        return mQuadratureMethodVector[DimensionIndex];
    }

    void SetQuadratureMethod(IndexType DimensionIndex,
        QuadratureMethod ThisQuadratureMethod);

    /* returns the IntegrationMethod to
     * corresponding to the direction index.
     */
    IntegrationMethod GetIntegrationMethod(
        IndexType DimensionIndex) const;

    /* Evaluates the corresponding IntegrationMethod to
     * the number of points and the quadrature method.
     */
    static IntegrationMethod GetIntegrationMethod(
        SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod)
    {
        switch (NumberOfIntegrationPointsPerSpan) {
        case 1:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                return IntegrationMethod::GI_GAUSS_1;
            }
            else {
                return IntegrationMethod::GI_LOBATTO_1;
            }
            break;
        case 2:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                return IntegrationMethod::GI_GAUSS_2;
            }
            break;
        case 3:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                return IntegrationMethod::GI_GAUSS_3;
            }
            break;
        case 4:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                return IntegrationMethod::GI_GAUSS_4;
            }
            break;
        case 5:
            if (ThisQuadratureMethod == QuadratureMethod::GAUSS) {
                return IntegrationMethod::GI_GAUSS_5;
            }
            break;
        case 0:
            return IntegrationMethod::NumberOfIntegrationMethods;
            break;
        }
        KRATOS_WARNING("Evaluation of Integration Method")
            << "Chosen combination of number of points per span and quadrature method does not has a corresponding IntegrationMethod in the KRATOS core."
            << "NumberOfIntegrationPointsPerSpan: " << NumberOfIntegrationPointsPerSpan << std::endl;
        return IntegrationMethod::NumberOfIntegrationMethods;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << " Integration info with local space dimension: " << mNumberOfIntegrationPointsPerSpanVector.size()
            << " and number of integration points per spans: " << mNumberOfIntegrationPointsPerSpanVector;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::vector<SizeType> mNumberOfIntegrationPointsPerSpanVector;

    std::vector<QuadratureMethod> mQuadratureMethodVector;

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
