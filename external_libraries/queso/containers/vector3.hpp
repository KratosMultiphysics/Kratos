// Author: Manuel Meßmer
// Email: manuel.messmer@tum.de

#ifndef VECTOR3_INCLUDE_H
#define VECTOR3_INCLUDE_H

//// STL includes
#include <array>
#include <cmath>
#include <iostream>

namespace queso {

///@name QuESo Classes
///@{

/**
 * @class  Vector3
 * @brief Aggregates std::array
 * @author Manuel Messmer
*/
template<typename type>
class Vector3
{
public:
    ///@name Type defintion
    ///@{
    typedef std::array<type,3> BaseType;

    ///@}
    ///@name Life cycle
    ///@{

    /// Default constructor
    Vector3() : mData() {
    };

    /// Constructor
    Vector3(type x, type y, type z) {
        mData[0] = x;
        mData[1] = y;
        mData[2] = z;
    }

    /// Destructor
    virtual ~Vector3() = default;

    /// Copy Constructor from BaseType
    Vector3(const BaseType& rOther) {
        mData[0] = rOther[0];
        mData[1] = rOther[1];
        mData[2] = rOther[2];
    }

    /// Copy Constructor from Vector3
    Vector3(const Vector3& rOther) {
        mData[0] = rOther[0];
        mData[1] = rOther[1];
        mData[2] = rOther[2];
    }

    /// Copy Assignement from Vector3
    Vector3& operator=(const Vector3& rOther) {
        mData[0] = rOther[0];
        mData[1] = rOther[1];
        mData[2] = rOther[2];

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Returns member std::array
    const BaseType& data() const {
        return mData;
    }

    /// Returns member std::array
    const BaseType& Coordinates() const {
        return mData;
    }

    /// Returns const_iterator to std::array
    typename BaseType::const_iterator begin() const {
        return mData.begin();
    }

    /// Returns const_iterator to std::array
    typename BaseType::const_iterator end() const {
        return mData.end();
    }

    /// Access first element
    type& X(){
        return mData[0];
    }

    /// Access second element
    type& Y(){
        return mData[1];
    }

    /// Access third element
    type& Z(){
        return mData[2];
    }

    /// Access first element const
    type X() const{
        return mData[0];
    }

    /// Access second element const
    type Y() const{
        return mData[1];
    }

    /// Access third element const
    type Z() const{
        return mData[2];
    }

    /// Access elements by index
    type& operator [] (std::size_t i){
        return mData[i];
    }

    /// Access elements by index (const)
    type operator [] (std::size_t i) const{
        return mData[i];
    }

    ///@}
    ///@name Math operators
    ///@{

    Vector3 operator+ (const type& rValue) const {
        return Vector3(mData[0] + rValue, mData[1] + rValue, mData[2] + rValue);
    }

    Vector3 operator+ (const Vector3& rOther) const {
        return Vector3(mData[0] + rOther[0], mData[1] + rOther[1], mData[2] + rOther[2]);
    }

    Vector3 operator- (const type& rValue) const {
        return Vector3(mData[0] - rValue, mData[1] - rValue, mData[2] - rValue);
    }

    Vector3 operator- (const Vector3& rOther) const {
        return Vector3(mData[0] - rOther[0], mData[1] - rOther[1], mData[2] - rOther[2]);
    }

    Vector3 operator* (const type& rValue) const {
        return Vector3(mData[0] * rValue, mData[1] * rValue, mData[2] * rValue);
    }

    Vector3 operator* (const Vector3& rOther) const {
        return Vector3(mData[0] * rOther[0], mData[1] * rOther[1], mData[2] * rOther[2]);
    }

    Vector3 operator/ (const type& rValue) const {
        return Vector3(mData[0] / rValue, mData[1] / rValue, mData[2] / rValue);
    }

    Vector3 operator/ (const Vector3& rOther) const {
        return Vector3(mData[0] / rOther[0], mData[1] / rOther[1], mData[2] / rOther[2]);
    }

    Vector3 operator+= (const Vector3& rOther)  {
        mData[0] += rOther[0];
        mData[1] += rOther[1];
        mData[2] += rOther[2];
        return *this;
    }

    Vector3 operator-= (const Vector3& rOther)  {
        mData[0] -= rOther[0];
        mData[1] -= rOther[1];
        mData[2] -= rOther[2];
        return *this;
    }

    Vector3 operator*= (const type& rValue)  {
        mData[0] *= rValue;
        mData[1] *= rValue;
        mData[2] *= rValue;
        return *this;
    }

    Vector3 operator/= (const type& rValue)  {
        mData[0] /= rValue;
        mData[1] /= rValue;
        mData[2] /= rValue;
        return *this;
    }

    double Norm() const {
        return std::sqrt( mData[0]*mData[0]
                         +mData[1]*mData[1]
                         +mData[2]*mData[2] );
    }

    ///@}
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << '(' << mData[0] << ", " << mData[1] << ", " << mData[2] << ')';
    }
private:
    std::array<type,3> mData;
    ///@}
}; // End Vector3 class
///@} // End QuESo classes

template<typename type>
std::ostream& operator<<(std::ostream& rOStream, const Vector3<type>& rThis)  {
    rThis.PrintInfo(rOStream);
    return rOStream;
}

} // End namespace queso

#endif // VECTOR3_INCLUDE_H