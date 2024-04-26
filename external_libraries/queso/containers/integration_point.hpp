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
    IntegrationPoint() : PointType()
    {}

    /// 2D Constructor
    IntegrationPoint(double x, double y, double weigth_) :
        PointType(x,y,0.0), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    /// 3D Constructor
    IntegrationPoint(double x, double y, double z, double weigth_) :
        PointType(x,y,z), mWeight(weigth_)
    {
        mActiveFlag = true;
    }
    /// 3D Constructor
    IntegrationPoint(const PointType& rPoint, double weigth_) :
        PointType(rPoint), mWeight(weigth_)
    {
        mActiveFlag = true;
    }

    /// Destructor
    ~IntegrationPoint() = default;

    /// Copy Constructor
    IntegrationPoint(const IntegrationPoint& rOther)
        : PointType(rOther), mWeight(rOther.mWeight), mActiveFlag(rOther.mActiveFlag)
    {
    }

    /// Assignement operator
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {
        PointType::operator=(rOther);
        if( this != &rOther) {
            mWeight = rOther.mWeight;
            mActiveFlag = rOther.mActiveFlag;

        }
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    double GetWeight() const{
        return mWeight;
    }

    void SetWeight(double weigth_){
        mWeight = weigth_;
    }

    void SetActive(bool value){
        mActiveFlag = value;
    }

    bool IsActive(){
        return mActiveFlag;
    }

    ///@}
private:
    ///@name Private member variables
    ///@{
    double mWeight;
    bool mActiveFlag;
    ///@}
}; // End class IntegrationPoint
///@} End QuESo classes

} // End namespace queso

#endif // INTEGRATION_POINT_INCLUDE_H