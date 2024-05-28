// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef BOUNDARY_INTEGRATION_POINT_INCLUDE_H
#define BOUNDARY_INTEGRATION_POINT_INCLUDE_H

//// Project includes
#include "containers/integration_point.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  BoundaryIntegrationPoint
 * @author Manuel Messmer
 * @brief Derives from IntegrationPoint. Contains also Normal Vector.
*/
class BoundaryIntegrationPoint : public IntegrationPoint
{
public:

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    BoundaryIntegrationPoint() : IntegrationPoint(), mNormal{0.0, 0.0, 0.0}
    {}

    /// 3D Constructor
    BoundaryIntegrationPoint(double X, double Y, double Z, double Weigth, const Vector3d& rNormal) :
        IntegrationPoint(X,Y,Z, Weigth), mNormal(rNormal)
    {
    }

    /// Destructor
    ~BoundaryIntegrationPoint() = default;

    /// Copy Constructor
    BoundaryIntegrationPoint(const BoundaryIntegrationPoint& rOther)
        : IntegrationPoint(rOther), mNormal(rOther.mNormal)
    {
    }

    /// Assignement operator
    BoundaryIntegrationPoint& operator=(const BoundaryIntegrationPoint& rOther)
    {
        IntegrationPoint::operator=(rOther);
        mNormal[0] = rOther.mNormal[0];
        mNormal[1] = rOther.mNormal[1];
        mNormal[2] = rOther.mNormal[2];

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// @brief Get Normal vector
    /// @return
    const Vector3d& Normal() const{
        return mNormal;
    }

    ///@}

private:
    ///@name Private Member Variables
    ///@{
    Vector3d mNormal;
    ///@}
}; // End class BoundaryIntegrationPoint
///@} // End QuESo classes

} // End namespace queso

#endif // BOUNDARY_INTEGRATION_POINT_INCLUDE_H