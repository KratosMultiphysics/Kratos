// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef SINGLE_ELEMENT_INCLUDE_HPP
#define SINGLE_ELEMENT_INCLUDE_HPP

//// STL includes
#include <vector>
#include <array>
//// Project includes
#include "queso/containers/element.hpp"
#include "queso/quadrature/integration_points_1d/integration_points_factory_1d.h"

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
    static void AssembleIPs(ElementType& rElement, const Vector3i& rOrder, IntegrationMethodType Method = IntegrationMethod::Gauss) {
        auto& integration_points = rElement.GetIntegrationPoints();

        PointType lower_bound_param = rElement.GetBoundsUVW().first;
        PointType upper_bound_param = rElement.GetBoundsUVW().second;

        AssembleIPs(integration_points, lower_bound_param, upper_bound_param, rOrder, Method);
    }

    /// @brief Assemble tensor product quadrature rules.
    /// @note This functions clears rIntegrationPoints.
    /// @param[out] rIntegrationPoints
    /// @param rLowerBoundParam LowerBound of element in parametric space.
    /// @param rUpperBoundParam LowerBound of element in parametric space.
    /// @param rOrder Order of quadrature rule.
    /// @param Method Integration method: Default - Gauss.
    static void AssembleIPs(IntegrationPointVectorType& rIntegrationPoints, const PointType& rLowerBoundParam, const PointType& rUpperBoundParam,
                            const Vector3i& rOrder, IntegrationMethodType Method = IntegrationMethod::Gauss ) {
        const auto p_ip_list_u = IntegrationPointFactory1D::GetGauss(rOrder[0], Method);
        const auto p_ip_list_v = IntegrationPointFactory1D::GetGauss(rOrder[1], Method);
        const auto p_ip_list_w = IntegrationPointFactory1D::GetGauss(rOrder[2], Method);

        const SizeType n_point_u = p_ip_list_u->size();
        const SizeType n_point_v = p_ip_list_v->size();
        const SizeType n_point_w = p_ip_list_w->size();

        const SizeType n_point = (n_point_u)*(n_point_v)*(n_point_w);
        rIntegrationPoints.clear();
        rIntegrationPoints.reserve(n_point);

        const double length_u = std::abs(rUpperBoundParam[0] - rLowerBoundParam[0]);
        const double length_v = std::abs(rUpperBoundParam[1] - rLowerBoundParam[1]);
        const double length_w = std::abs(rUpperBoundParam[2] - rLowerBoundParam[2]);

        for (SizeType u = 0; u < p_ip_list_u->size(); ++u) {
            for (SizeType v = 0; v < p_ip_list_v->size(); ++v) {
                for( SizeType w = 0; w < p_ip_list_w->size(); ++w) {
                    rIntegrationPoints.push_back( IntegrationPointType( rLowerBoundParam[0] + length_u * (*p_ip_list_u)[u][0],
                                                                        rLowerBoundParam[1] + length_v * (*p_ip_list_v)[v][0],
                                                                        rLowerBoundParam[2] + length_w * (*p_ip_list_w)[w][0],
                                                                        (*p_ip_list_u)[u][1] * length_u *
                                                                        (*p_ip_list_v)[v][1] * length_v *
                                                                        (*p_ip_list_w)[w][1] * length_w ) );
                }
            }
        }
    }
    ///@}

}; // End Class QuadratureSingleElement

///@}

} // End namespace queso

#endif // SINGLE_ELEMENT_INCLUDE_HPP