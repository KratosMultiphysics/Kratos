// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_GAUSS_POINT_ITEM )
#define  KRATOS_GAUSS_POINT_ITEM

// System includes

// External includes

// Project includes
#include "geometries/point.h"
#include "includes/constitutive_law.h"

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
 * @class GaussPointItem
 * @ingroup MeshingApplication
 * @brief Custom Gauss Point container to be used by the search
 * @author Vicente Mataix Ferrandiz
 */
class GaussPointItem 
    : public Point
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of GaussPointItem
    KRATOS_CLASS_POINTER_DEFINITION( GaussPointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief  Default constructor
     * @details It just computes the default constructor of the point 
     */
    GaussPointItem():
        Point()
    {
    }

    /**
     * @brief Constructor with coordinates
     * @details It just computes the coordinates constructor of the point 
     */
    GaussPointItem(const array_1d<double, 3>& Coordinates):
        Point(Coordinates)
    {
    }

    /**
     * @brief Complete constructor
     * @details Computes the point constructor + Considers the CL pointer and the integration weight 
     */
    GaussPointItem(
        const array_1d<double, 3>& Coordinates,
        ConstitutiveLaw::Pointer pConstitutiveLaw,
        const double Weight
        ):Point(Coordinates),
          mpConstitutiveLaw(std::move(pConstitutiveLaw)),
          mWeight(Weight)
    {
    }

    ///Copy constructor  (not really required)
    GaussPointItem(const GaussPointItem& GP)= default;

    /// Destructor.
    ~GaussPointItem() override= default;

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

    Point GetPoint()
    {
        Point Point(this->Coordinates());

        return Point;
    }

    /**
     * @brief Set the point
     * @param Point The point
     */

    void SetPoint(const Point& Point)
    {
        this->Coordinates() = Point.Coordinates();
    }

    /**
     * @brief Sets the Constitutive Law associated to the point
     * @param pConstitutiveLaw The pointer to the Constitutive Law
     */

    void SetConstitutiveLaw(ConstitutiveLaw::Pointer pConstitutiveLaw)
    {
        mpConstitutiveLaw = pConstitutiveLaw;
    }

    /**
     * @brief Returns the Constitutive Law associated to the point
     * @return mpConstitutiveLaw: The pointer to the Constitutive Law associated to the point
     */

    ConstitutiveLaw::Pointer GetConstitutiveLaw()
    {
        return mpConstitutiveLaw;
    }

    /**
     * @brief Returns the integration weigth associated to the point
     * @return mWeight: The pointer to the Constitutive Law associated to the point
     */

    double GetWeight() const
    {
        return mWeight;
    }

    /**
     * @brief Sets the integration weigth associated to the point
     * @param Weight The integration weight
     */

    void SetWeight(const double Weight)
    {
        mWeight = Weight;
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

    ConstitutiveLaw::Pointer mpConstitutiveLaw; // The constitutive law pointer
    double mWeight;                             // The integration weight of the GP

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
}; // Class GaussPointItem 
///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_GAUSS_POINT_ITEM  defined 
