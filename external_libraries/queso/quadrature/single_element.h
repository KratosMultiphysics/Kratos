// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef SINGLE_ELEMENT_H
#define SINGLE_ELEMENT_H

//// STL includes
#include <vector>
#include <array>
//// Project includes
#include "containers/element.hpp"
#include "quadrature/integration_points_1d/integration_points_factory_1d.h"

namespace queso {

///@name QuESo Classes
///@{

////
/**
 * @class  QuadratureSingleElement. Provides assembly opeartions for tensor-product quadrature rules of single non-trimmed element.
 * @author Manuel Messmer
 * @brief  Provides assembly for 3D quadrature rules.
 * @tparam TElementType
 * @details Available Quadrature rules:
 *          {Gauss, Gauss_Reduced1, Gauss_Reduced2}
*/
template<typename TElementType>
class QuadratureSingleElement {

public:
        ///@name Type Definitions
        ///@{

        typedef TElementType ElementType;
        typedef typename ElementType::IntegrationPointType IntegrationPointType;
        typedef typename ElementType::IntegrationPointVectorType IntegrationPointVectorType;

        ///@}
        ///@name Operations
        ///@{

        /// @brief Assemble tensor product quadrature rules.
        /// @param rElement
        /// @param rOrder Order of quadrature rule.
        /// @param Method Integration method: Default - Gauss.
        static void AssembleIPs(ElementType& rElement, const Vector3i& rOrder, IntegrationMethodType Method = IntegrationMethod::Gauss);

        /// @brief Assemble tensor product quadrature rules.
        /// @note This functions clears rIntegrationPoints.
        /// @param[out] rIntegrationPoints
        /// @param rLowerBoundParam LowerBound of element in parametric space.
        /// @param rUpperBoundParam LowerBound of element in parametric space.
        /// @param rOrder Order of quadrature rule.
        /// @param Method Integration method: Default - Gauss.
        static void AssembleIPs(IntegrationPointVectorType& rIntegrationPoints, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam,
                                const Vector3i& rOrder, IntegrationMethodType Method = IntegrationMethod::Gauss );
        ///@}

}; // End Class QuadratureSingleElement

///@}

} // End namespace queso

#endif // SINGLE_ELEMENT_H