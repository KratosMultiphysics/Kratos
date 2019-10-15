// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(POINT_ITEM_DEFINED )
#define  POINT_ITEM_DEFINED

// System includes

// External includes

// Project includes
#include "geometries/point.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class PointItem
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Custom Point container to be used by the search
 * @details It stores the pointer of a certain entity
 * @tparam TEntity The entity which pointer will be stored
 * @author Vicente Mataix Ferrandiz
 */
template<class TEntity>
class PointItem
    : public Point
{
public:

    ///@name Type Definitions
    ///@{

    typedef Point BaseType;

    /// Counted pointer of PointItem
    KRATOS_CLASS_POINTER_DEFINITION( PointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointItem():
        BaseType(),
        mpOriginEntity(nullptr)
    {}

    PointItem(const array_1d<double, 3>& rCoordinates)
        :BaseType(rCoordinates),
         mpOriginEntity(nullptr)
    {}

    PointItem(typename TEntity::Pointer pEntity):
        mpOriginEntity(pEntity)
    {
        UpdatePoint();
    }

    PointItem(
        const array_1d<double, 3>& rCoordinates,
        typename TEntity::Pointer pEntity
    ):
        BaseType(rCoordinates),
        mpOriginEntity(pEntity)
    {}

    ///Copy constructor  (not really required)
    PointItem(const PointItem& rRHS):
        BaseType(rRHS),
        mpOriginEntity(rRHS.mpOriginEntity)
    {
    }

    /// Destructor.
    ~PointItem() override= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the point
     * @return The point
     */
    BaseType GetPoint()
    {
        BaseType Point(this->Coordinates());

        return Point;
    }

    /**
     * @brief Set the point
     * @param rPoint The point
     */
    void SetPoint(const BaseType& rPoint)
    {
        this->Coordinates() = rPoint.Coordinates();
    }

    /**
     * @brief Sets the condition associated to the point
     * @param pEntity The pointer to the condition
     */

    void SetEntity(typename TEntity::Pointer pEntity)
    {
        mpOriginEntity = pEntity;
    }

    /**
     * @brief Returns the condition associated to the point
     * @return mpOriginEntity The pointer to the condition associated to the point
     */

    typename TEntity::Pointer GetEntity()
    {
    #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(mpOriginEntity.get() == nullptr) << "TEntity no initialized in the PointItem class" << std::endl;
    #endif
        return mpOriginEntity;
    }

    /**
     * This method checks everything is right
     */

    void Check()
    {
        KRATOS_TRY;

        auto aux_coord = std::make_shared<array_1d<double, 3>>(this->Coordinates());
        KRATOS_ERROR_IF(!aux_coord) << "Coordinates no initialized in the PointItem class" << std::endl;
        KRATOS_ERROR_IF(mpOriginEntity.get() == nullptr) << "TEntity no initialized in the PointItem class" << std::endl;

        KRATOS_CATCH("Error checking the PointItem");
    }

    /**
     * @brief This function updates the database, using as base for the coordinates the condition center
     */
    void UpdatePoint()
    {
        noalias(this->Coordinates()) = mpOriginEntity->GetGeometry().Center().Coordinates();
    }

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    typename TEntity::Pointer mpOriginEntity; // Entity pointer

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PointItem

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // POINT_ITEM_DEFINED  defined
