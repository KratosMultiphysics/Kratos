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

// #if !defined(KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED)
// #define KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED

// // System includes
// #include <string>

// // External includes

// // Project includes

// // Application includes
// #include "statistics_method.h"

// namespace Kratos
// {
// ///@addtogroup RANSApplication
// ///@{

// ///@name Kratos Globals
// ///@{

// ///@}
// ///@name Type Definitions
// ///@{

// ///@}
// ///@name  Enum's
// ///@{

// ///@}
// ///@name  Functions
// ///@{

// ///@}
// ///@name Kratos Classes
// ///@{

// /**
//  * @brief Clips given scalar variable to a range
//  *
//  * This process clips a given scalar variable to a range in all nodes in the model part.
//  *
//  */
// template <typename TDataType>
// class KRATOS_API(STATISTICS_APPLICATION) TemporalVarianceMethod
//     : public StatisticsMethod<TDataType, TDataType, TDataType, double, double>
// {
// public:
//     ///@name Type Definitions
//     ///@{

//     /// Pointer definition of TemporalVarianceMethod
//     KRATOS_CLASS_POINTER_DEFINITION(TemporalVarianceMethod);

//     ///@}
//     ///@name Life Cycle
//     ///@{

//     /// Constructor

//     TemporalVarianceMethod();

//     /// Destructor.
//     ~TemporalVarianceMethod()  = default;

//     ///@}
//     ///@name Operators
//     ///@{

//     ///@}
//     ///@name Operations
//     ///@{

//     void CalculateStatistics(TDataType& rMean,
//                              TDataType& rVariance,
//                              TDataType const& rNewDataPoint,
//                              double const& rDeltaTime,
//                              double const& rTotalTime);

//     ///@}
//     ///@name Access
//     ///@{

//     ///@}
//     ///@name Inquiry
//     ///@{

//     ///@}
//     ///@name Input and output
//     ///@{

//     ///@}
//     ///@name Friends
//     ///@{

//     ///@}

// protected:
//     ///@name Protected static Member Variables
//     ///@{

//     ///@}
//     ///@name Protected member Variables
//     ///@{

//     ///@}
//     ///@name Protected Operators
//     ///@{

//     ///@}
//     ///@name Protected Operations
//     ///@{

//     ///@}
//     ///@name Protected  Access
//     ///@{

//     ///@}
//     ///@name Protected Inquiry
//     ///@{

//     ///@}
//     ///@name Protected LifeCycle
//     ///@{

//     ///@}

// private:
//     ///@name Static Member Variables
//     ///@{

//     ///@}
//     ///@name Member Variables
//     ///@{

//     ///@}
//     ///@name Private Operators
//     ///@{

//     ///@}
//     ///@name Private Operations
//     ///@{

//     ///@}
//     ///@name Private  Access
//     ///@{

//     ///@}
//     ///@name Private Inquiry
//     ///@{

//     ///@}
//     ///@name Un accessible methods
//     ///@{

//     /// Assignment operator.
//     TemporalVarianceMethod& operator=(TemporalVarianceMethod const& rOther);

//     /// Copy constructor.
//     TemporalVarianceMethod(TemporalVarianceMethod const& rOther);

//     ///@}
// };

// ///@}

// ///@name Type Definitions
// ///@{

// ///@}
// ///@name Input and output
// ///@{

// /// output stream function
// template <typename TDataType>
// inline std::ostream& operator<<(std::ostream& rOStream,
//                                 const TemporalVarianceMethod<TDataType>& rThis);

// ///@}

// ///@} addtogroup block

// } // namespace Kratos.

// #endif // KRATOS_TEMPORAL_VARIANCE_METHOD_H_INCLUDED defined
