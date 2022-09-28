//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "integration_info.h"

namespace Kratos
{
    ///@name Local Flags
    ///@{

    KRATOS_CREATE_LOCAL_FLAG(IntegrationInfo, DO_NOT_CREATE_TESSELLATION_ON_SLAVE, 0);

    ///@}
    ///@name Constructors
    ///@{

    IntegrationInfo::IntegrationInfo(SizeType LocalSpaceDimension,
        IntegrationMethod ThisIntegrationMethod)
    {
        mNumberOfIntegrationPointsPerSpanVector = std::vector<SizeType>(LocalSpaceDimension);
        mQuadratureMethodVector = std::vector<QuadratureMethod>(LocalSpaceDimension);

        for (IndexType i = 0; i < LocalSpaceDimension; ++i) {
            SetIntegrationMethod(i, ThisIntegrationMethod);
        }
    }

    IntegrationInfo::IntegrationInfo(SizeType LocalSpaceDimension,
        SizeType NumberOfIntegrationPointsPerSpan,
        QuadratureMethod ThisQuadratureMethod)
    {
        mNumberOfIntegrationPointsPerSpanVector = std::vector<SizeType>(LocalSpaceDimension);
        mQuadratureMethodVector = std::vector<QuadratureMethod>(LocalSpaceDimension);

        for (IndexType i = 0; i < LocalSpaceDimension; ++i) {
            mNumberOfIntegrationPointsPerSpanVector[i] = NumberOfIntegrationPointsPerSpan;
            mQuadratureMethodVector[i] = ThisQuadratureMethod;
        }
    }

    IntegrationInfo::IntegrationInfo(
        const std::vector<SizeType>& NumberOfIntegrationPointsPerSpanVector,
        const std::vector<QuadratureMethod>& ThisQuadratureMethodVector)
        : mNumberOfIntegrationPointsPerSpanVector(NumberOfIntegrationPointsPerSpanVector),
          mQuadratureMethodVector(ThisQuadratureMethodVector)
    {
        KRATOS_ERROR_IF(NumberOfIntegrationPointsPerSpanVector.size() != ThisQuadratureMethodVector.size())
            << "The sizes of the NumberOfIntegrationPointsPerSpanVector: " << NumberOfIntegrationPointsPerSpanVector.size()
            << " and the ThisQuadratureMethodVector: " << ThisQuadratureMethodVector.size() << " does not coincide." << std::endl;
    }

    ///@}
    ///@name integration rules
    ///@{

    void IntegrationInfo::SetIntegrationMethod(
        IndexType DimensionIndex,
        IntegrationMethod ThisIntegrationMethod)
    {
        switch (ThisIntegrationMethod) {
        case IntegrationMethod::GI_GAUSS_1:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 1;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_2:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 2;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_3:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 3;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_4:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 4;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_GAUSS_5:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 5;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_1:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 1;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_2:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 2;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_3:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 3;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_4:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 4;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::GI_EXTENDED_GAUSS_5:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 5;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::EXTENDED_GAUSS;
            break;
        case IntegrationMethod::NumberOfIntegrationMethods:
            mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = 0;
            mQuadratureMethodVector[DimensionIndex] = QuadratureMethod::Default;
            break;
        }
    }

    IntegrationInfo::SizeType IntegrationInfo::GetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex) const
    {
        return mNumberOfIntegrationPointsPerSpanVector[DimensionIndex];
    }

    void IntegrationInfo::SetNumberOfIntegrationPointsPerSpan(IndexType DimensionIndex,
        SizeType NumberOfIntegrationPointsPerSpan)
    {
        mNumberOfIntegrationPointsPerSpanVector[DimensionIndex] = NumberOfIntegrationPointsPerSpan;
    }

    void IntegrationInfo::SetQuadratureMethod(IndexType DimensionIndex,
        QuadratureMethod ThisQuadratureMethod)
    {
        mQuadratureMethodVector[DimensionIndex] = ThisQuadratureMethod;
    }

    /* returns the IntegrationMethod to
     * corresponding to the direction index.
     */
    IntegrationInfo::IntegrationMethod IntegrationInfo::GetIntegrationMethod(
        IndexType DimensionIndex) const
    {
        return GetIntegrationMethod(
            GetNumberOfIntegrationPointsPerSpan(DimensionIndex),
            GetQuadratureMethod(DimensionIndex));
    }

    ///@}

}  // namespace Kratos.