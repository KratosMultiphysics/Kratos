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
class IntegrationPoint
{
public:

    ///@name Type defintion
    ///@{
    typedef std::array<double,3> BaseType;

    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor
    IntegrationPoint() : mPoint{}, mWeight(0.0)
    {
    }

    /// 2D Constructor
    IntegrationPoint(double x, double y, double weigth_) :
        mPoint{x,y,0.0}, mWeight(weigth_)
    {
    }

    /// 3D Constructor
    IntegrationPoint(double x, double y, double z, double weigth_) :
        mPoint{x,y,z}, mWeight(weigth_)
    {
    }

    /// 3D Constructor
    IntegrationPoint(const BaseType& rPoint, double weigth_) :
        mPoint(rPoint), mWeight(weigth_)
    {
    }

    /// Destructor
    virtual ~IntegrationPoint() = default;

    /// Copy Constructor
    IntegrationPoint(const IntegrationPoint& rOther)
        : mPoint(rOther.data()), mWeight(rOther.mWeight)
    {
    }

    /// Assignement operator
    IntegrationPoint& operator=(const IntegrationPoint& rOther)
    {
        mPoint = rOther.mPoint;
        mWeight = rOther.mWeight;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Returns underlying point data.
    const BaseType& data() const {
        return mPoint;
    }

    /// Access first element
    double& X(){
        return mPoint[0];
    }

    /// Access second element
    double& Y(){
        return mPoint[1];
    }

    /// Access third element
    double& Z(){
        return mPoint[2];
    }

    /// Access first element const
    double X() const{
        return mPoint[0];
    }

    /// Access second element const
    double Y() const{
        return mPoint[1];
    }

    /// Access third element const
    double Z() const{
        return mPoint[2];
    }

    /// Access elements by index
    double& operator [] (std::size_t i){
        return mPoint[i];
    }

    /// Access elements by index (const)
    double operator [] (std::size_t i) const{
        return mPoint[i];
    }

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
    BaseType mPoint;
    double mWeight;
    ///@}
}; // End class IntegrationPoint
///@} End QuESo classes

} // End namespace queso

#endif // INTEGRATION_POINT_INCLUDE_H