// //    |  /           |
// //    ' /   __| _` | __|  _ \   __|
// //    . \  |   (   | |   (   |\__ `
// //   _|\_\_|  \__,_|\__|\___/ ____/
// //                   Multi-Physics
// //
// //  License:		 BSD License
// //					 Kratos default license: kratos/license.txt
// //
// //  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
// //

// // System includes
// #include <cmath>

// // External includes

// // Project includes

// // Application includes

// // Include base h
// #include "temporal_variance_method.h"

// namespace Kratos
// {
// namespace
// {
// template <typename TDataType>
// void CalculateTemporalVariance(TDataType& rMean,
//                                TDataType& rVariance,
//                                TDataType const& rNewDataPoint,
//                                double const& rDeltaTime,
//                                double const& rTotalTime)
// {
//     const double new_total_time = rTotalTime + rDeltaTime;
//     const double new_mean = (rMean * rTotalTime + rNewDataPoint) / (new_total_time);
//     rVariance = ((rVariance + std::pow(new_mean - rMean, 2)) * rTotalTime +
//                  std::pow(rNewDataPoint - new_mean, 2)) /
//                 new_total_time;
//     rMean = new_mean;
// }

// // template instantiations
// template void CalculateTemporalVariance<int>(int& rMean,
//                                              int& rVariance,
//                                              int const& rNewDataPoint,
//                                              double const& rDeltaTime,
//                                              double const& rTotalTime);

// template void CalculateTemporalVariance<double>(double& rMean,
//                                                 double& rVariance,
//                                                 double const& rNewDataPoint,
//                                                 double const& rDeltaTime,
//                                                 double const& rTotalTime);

// } // namespace

// template <typename TDataType>
// TemporalVarianceMethod<TDataType>::TemporalVarianceMethod()
// {
// }

// template <typename TDataType>
// void TemporalVarianceMethod<TDataType>::CalculateStatistics(TDataType& rMean,
//                                                             TDataType& rVariance,
//                                                             TDataType const& rNewDataPoint,
//                                                             double const& rDeltaTime,
//                                                             double const& rTotalTime)
// {
//     CalculateTemporalVariance(rMean, rVariance, rNewDataPoint, rDeltaTime, rTotalTime);
// }

// template <>
// void TemporalVarianceMethod<array_1d<double, 3>>::CalculateStatistics(
//     array_1d<double, 3>& rMean,
//     array_1d<double, 3>& rVariance,
//     array_1d<double, 3> const& rNewDataPoint,
//     double const& rDeltaTime,
//     double const& rTotalTime)
// {
//     for (int i = 0; i < 3; ++i)
//         CalculateTemporalVariance(rMean[i], rVariance[i], rNewDataPoint[i],
//                                   rDeltaTime, rTotalTime);
// }

// template <>
// void TemporalVarianceMethod<Vector>::CalculateStatistics(Vector& rMean,
//                                                          Vector& rVariance,
//                                                          Vector const& rNewDataPoint,
//                                                          double const& rDeltaTime,
//                                                          double const& rTotalTime)
// {
//     KRATOS_TRY

//     KRATOS_ERROR_IF(rMean.size() != rNewDataPoint.size())
//         << "Mean vector size not equal to data point vector size. [ "
//         << rMean.size() << " != " << rNewDataPoint.size() << " ].\n";
//     KRATOS_ERROR_IF(rVariance.size() != rNewDataPoint.size())
//         << "Variance vector size not equal to data point vector size. [ "
//         << rVariance.size() << " != " << rNewDataPoint.size() << " ].\n";

//     for (int i = 0; i < static_cast<int>(rMean.size()); ++i)
//         CalculateTemporalVariance(rMean[i], rVariance[i], rNewDataPoint[i],
//                                   rDeltaTime, rTotalTime);

//     KRATOS_CATCH("");
// }

// template <>
// void TemporalVarianceMethod<Matrix>::CalculateStatistics(Matrix& rMean,
//                                                          Matrix& rVariance,
//                                                          Matrix const& rNewDataPoint,
//                                                          double const& rDeltaTime,
//                                                          double const& rTotalTime)
// {
//     KRATOS_TRY

//     KRATOS_ERROR_IF(rMean.size1() != rNewDataPoint.size1())
//         << "Mean matrix size1 not equal to data point matrix size1. [ "
//         << rMean.size1() << " != " << rNewDataPoint.size1() << " ].\n";
//     KRATOS_ERROR_IF(rMean.size2() != rNewDataPoint.size2())
//         << "Mean matrix size2 not equal to data point matrix size2. [ "
//         << rMean.size2() << " != " << rNewDataPoint.size2() << " ].\n";

//     KRATOS_ERROR_IF(rVariance.size1() != rNewDataPoint.size1())
//         << "Variance matrix size1 not equal to data point matrix size1. [ "
//         << rVariance.size1() << " != " << rNewDataPoint.size1() << " ].\n";

//     KRATOS_ERROR_IF(rVariance.size2() != rNewDataPoint.size2())
//         << "Variance matrix size2 not equal to data point matrix size2. [ "
//         << rVariance.size2() << " != " << rNewDataPoint.size2() << " ].\n";

//     for (int i = 0; i < static_cast<int>(rMean.size1()); ++i)
//         for (int j = 0; j < static_cast<int>(rMean.size1()); ++j)
//             CalculateTemporalVariance(rMean(i, j), rVariance(i, j),
//                                       rNewDataPoint(i, j), rDeltaTime, rTotalTime);

//     KRATOS_CATCH("");
// }

// // template instantiations
// template class TemporalVarianceMethod<int>;
// template class TemporalVarianceMethod<double>;
// template class TemporalVarianceMethod<array_1d<double, 3>>;
// template class TemporalVarianceMethod<Vector>;
// template class TemporalVarianceMethod<Matrix>;

// } // namespace Kratos.
