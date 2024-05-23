// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef INTEGRATION_POINT_INCLUDE_H
#define INTEGRATION_POINT_INCLUDE_H

//// Project includes
#include "includes/define.hpp"

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  IntegrationPoint
 * @author Manuel Messmer
*/
class IntegrationPoint : public PointType
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor
    IntegrationPoint() : PointType(), mWeight(0.0)
    {
    }

    /// 2D Constructor
    IntegrationPoint(double x, double y, double weigth_) :
        PointType(x,y,0.0), mWeight(weigth_)
    {
    }

    /// 3D Constructor
    IntegrationPoint(double x, double y, double z, double weigth_) :
        PointType(x,y,z), mWeight(weigth_)
    {
    }

    /// 3D Constructor
    IntegrationPoint(const PointType& rPoint, double weigth_) :
        PointType(rPoint), mWeight(weigth_)
    {
    }

    /// Destructor
    virtual ~IntegrationPoint() = default;

    /// Copy Constructor
    IntegrationPoint(const IntegrationPoint& rOther)
        : PointType(rOther), mWeight(rOther.mWeight)
    {
    }

    /// Assignement operator
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {

        PointType::operator=(rOther);
        mWeight = rOther.mWeight;

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Get integration weight
    double Weight() const{
        return mWeight;
    }

    /// Set integration weight
    void SetWeight(double weigth_){
        mWeight = weigth_;
    }

    ///@}
private:
    ///@name Private member variables
    ///@{
    double mWeight;
    ///@}
}; // End class IntegrationPoint
///@} End QuESo classes

} // End namespace queso

#endif // INTEGRATION_POINT_INCLUDE_H