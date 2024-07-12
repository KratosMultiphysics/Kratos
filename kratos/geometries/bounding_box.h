//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Representing a bounding box by storing the min and max points
 * @details It stores the min and max points and have constructor for it construction with any container of points.
 *  TPointType should provide access operator [] to its coordinate and deep copy operator=
 * @tparam TPointType The type of point considered
 * @author Pooyan Dadvand
 */
template <typename TPointType>
class BoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(BoundingBox);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    BoundingBox()
    {
        std::fill(GetMinPoint().begin(), GetMinPoint().end(), 0.0);
        std::fill(GetMaxPoint().begin(), GetMaxPoint().end(), 0.0);
    };

    /// Constructor with min and max points
    BoundingBox(TPointType const& MinPoint, TPointType const& MaxPoint)
    {
        noalias(GetMinPoint().Coordinates()) = MinPoint.Coordinates();
        noalias(GetMaxPoint().Coordinates()) = MaxPoint.Coordinates();
    }

    /// Copy constructor
    BoundingBox( const BoundingBox &Other) :
		mMinMaxPoints(Other.mMinMaxPoints) {}


    /// Construction with container of points.
    template<typename TIteratorType>
    BoundingBox(TIteratorType const& itPointsBegin, TIteratorType const& itPointsEnd) 
    {
        Set(itPointsBegin, itPointsEnd);
    }

    /// Destructor.
    virtual ~BoundingBox(){}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BoundingBox& operator=(BoundingBox const& rOther)
    {
        GetMinPoint() = rOther.GetMinPoint();
        GetMaxPoint() = rOther.GetMaxPoint();

        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Sets the minimum and maximum points based on a range of input points.
     * @details This function sets the minimum and maximum points of the object based on a range
     * of input points specified by the iterators `itPointsBegin` and `itPointsEnd`. If the
     * range is empty (itPointsBegin == itPointsEnd), it initializes both the minimum and
     * maximum points to zero vectors.
     * @tparam TIteratorType The iterator type for the input points.
     * @param itPointsBegin The iterator pointing to the beginning of the input point range.
     * @param itPointsEnd The iterator pointing to the end of the input point range.
     */
    template<typename TIteratorType>
    void Set(TIteratorType const& itPointsBegin, TIteratorType const& itPointsEnd)
    {
        if (itPointsBegin == itPointsEnd) {
            std::fill(GetMinPoint().begin(), GetMinPoint().end(), 0.0);
            std::fill(GetMaxPoint().begin(), GetMaxPoint().end(), 0.0);
            return;
        }

        // Initialize the min and max points to the first point
        auto& r_min_point = GetMinPoint();
        auto& r_max_point = GetMaxPoint();
        const auto& r_coordinates = itPointsBegin->Coordinates();
        for (unsigned int i = 0; i < Dimension; i++) {
            r_min_point[i] = r_coordinates[i];
            r_max_point[i] = r_coordinates[i];
        }

        Extend(itPointsBegin, itPointsEnd);
    }

    /**
     * @brief Extends the bounding box to include a range of input points.
     * @details This function extends the bounding box defined by the minimum and maximum points
     * to include a range of input points specified by the iterators `itPointsBegin` and
     * `itPointsEnd`. It adjusts the minimum and maximum points as necessary to encompass
     * all input points.
     * @tparam TIteratorType The iterator type for the input points.
     * @param itPointsBegin The iterator pointing to the beginning of the input point range.
     * @param itPointsEnd The iterator pointing to the end of the input point range.
     */
    template<typename TIteratorType>
    void Extend(TIteratorType const& itPointsBegin, TIteratorType const& itPointsEnd)
    {
        // Extend the min and max points
        auto& r_min_point = GetMinPoint();
        auto& r_max_point = GetMaxPoint();
        for (TIteratorType it_point = itPointsBegin; it_point != itPointsEnd; it_point++){
            for (unsigned int i = 0; i < Dimension; i++) {
                if ((*it_point)[i] < r_min_point[i]) r_min_point[i] = (*it_point)[i];
                if ((*it_point)[i] > r_max_point[i]) r_max_point[i] = (*it_point)[i];
            }
        }
    }

    /**
     * @brief Extends the bounding box by adding a margin to its dimensions.
     * @details This function extends the bounding box by adding a specified margin to each dimension
     * of both the minimum and maximum points. It effectively enlarges the bounding box
     * in all directions.
     * @param Margin The margin value to be added to each dimension.
     */
    void Extend(const double Margin)
    {
        // Extend the min and max points
        auto& r_min_point = GetMinPoint();
        auto& r_max_point = GetMaxPoint();
        for (unsigned int i = 0; i < Dimension; i++){
            r_min_point[i] -= Margin;
            r_max_point[i] += Margin;
        }
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Gets a reference to the minimum point.
     * @details This function returns a reference to the minimum point stored in the object.
     * @return A reference to the minimum point.
     */
    TPointType& GetMinPoint()
    {
        return mMinMaxPoints[0];
    }

    /**
     * @brief Gets a constant reference to the minimum point (read-only).
     * @details This function returns a constant reference to the minimum point stored in the object. It allows you to access the minimum point without modifying it.
     * @return A constant reference to the minimum point.
     */
    TPointType const& GetMinPoint() const
    {
        return mMinMaxPoints[0];
    }

    /**
     * @brief Gets a reference to the maximum point.
     * @details This function returns a reference to the maximum point stored in the object.
     * @return A reference to the maximum point.
     */
    TPointType& GetMaxPoint()
    {
        return mMinMaxPoints[1];
    }

    /**
     * @brief Gets a constant reference to the maximum point (read-only).
     * @details This function returns a constant reference to the maximum point stored in the object. It allows you to access the maximum point without modifying it.
     * @return A constant reference to the maximum point.
     */
    TPointType const& GetMaxPoint() const 
    {
        return mMinMaxPoints[1];
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "BoundingBox" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "BoundingBox";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        rOStream << "   MinPoint : [" << GetMinPoint()[0] << ","  << GetMinPoint()[1] << ","  << GetMinPoint()[2] << "]" << std::endl;
        rOStream << "   MaxPoint : [" << GetMaxPoint()[0] << ","  << GetMaxPoint()[1] << ","  << GetMaxPoint()[2] << "]" << std::endl;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    static constexpr unsigned int Dimension = 3;

    ///@}
    ///@name Member Variables
    ///@{

    std::array<TPointType, 2> mMinMaxPoints;  /// The min and max points 

    ///@}

}; // Class BoundingBox

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <typename TPointType>
inline std::istream& operator >> (std::istream& rIStream,
                BoundingBox<TPointType>& rThis){
                    return rIStream;
                }

/// output stream function
template <typename TPointType>
inline std::ostream& operator << (std::ostream& rOStream,
                const BoundingBox<TPointType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
