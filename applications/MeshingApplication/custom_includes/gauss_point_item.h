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

    /// The type used to idenify index and key
    typedef std::size_t IndexType;

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

    /**
    * It checks if an ID exists in the map
    * @param VariableKey The condition ID to remove
    * @return If the ID already exists or not
    */
    bool Has(const IndexType VariableKey)
    {
        auto double_set = mMapDoubleVariables.find(VariableKey);
        if(double_set != mMapDoubleVariables.end()) {
            return true;
        }
        auto array_set = mMapArrayVariables.find(VariableKey);
        if(array_set != mMapArrayVariables.end()) {
            return true;
        }
        auto vector_set = mMapVectorVariables.find(VariableKey);
        if(vector_set != mMapVectorVariables.end()) {
            return true;
        }
        auto matrix_set = mMapMatrixVariables.find(VariableKey);
        if(matrix_set != mMapMatrixVariables.end()) {
            return true;
        }

        return false;
    }

    /**
    * It adds a new value to the map (double)
    * @param VariableKey The variable ID to set
    * @param rValue The value to assign
    */
    void SetValue(
        const IndexType VariableKey,
        const double rValue
        )
    {
        auto double_set = mMapDoubleVariables.find(VariableKey);
        if(double_set != mMapDoubleVariables.end()) {
            mMapDoubleVariables[VariableKey] = rValue;
            return void();
        }
        mMapDoubleVariables.insert({VariableKey, rValue});
    }

    /**
    * It adds a new value to the map (array_1d<double, 3>)
    * @param VariableKey The variable ID to set
    * @param rValue The value to assign
    */
    void SetValue(
        const IndexType VariableKey,
        const array_1d<double, 3>& rValue
        )
    {
        auto array_set = mMapArrayVariables.find(VariableKey);
        if(array_set != mMapArrayVariables.end()) {
            mMapArrayVariables[VariableKey] = rValue;
            return void();
        }
        mMapArrayVariables.insert({VariableKey, rValue});
    }

    /**
    * It adds a new value to the map (Vector)
    * @param VariableKey The variable ID to set
    * @param rValue The value to assign
    */
    void SetValue(
        const IndexType VariableKey,
        const Vector& rValue
        )
    {
        auto vector_set = mMapVectorVariables.find(VariableKey);
        if(vector_set != mMapVectorVariables.end()) {
            mMapVectorVariables[VariableKey] = rValue;
            return void();
        }
        mMapVectorVariables.insert({VariableKey, rValue});
    }

    /**
    * It adds a new value to the map (Matrix)
    * @param VariableKey The variable ID to set
    * @param rValue The value to assign
    */
    void SetValue(
        const IndexType VariableKey,
        const Matrix& rValue
        )
    {
        auto matrix_set = mMapMatrixVariables.find(VariableKey);
        if(matrix_set != mMapMatrixVariables.end()) {
            mMapMatrixVariables[VariableKey] = rValue;
            return void();
        }
        mMapMatrixVariables.insert({VariableKey, rValue});
    }

    /**
    * It return a value from the map (double)
    * @param VariableKey The variable ID to set
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    double GetValue(
        const IndexType VariableKey,
        double& rValue
        )
    {

        rValue = mMapDoubleVariables[VariableKey];
        return rValue;
    }

    /**
    * It return a value from the map (array_1d<double, 3>)
    * @param VariableKey The variable ID to set
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    array_1d<double, 3> GetValue(
        const IndexType VariableKey,
        array_1d<double, 3>& rValue
        )
    {
        rValue = mMapArrayVariables[VariableKey];
        return rValue;
    }

    /**
    * It return a value from the map (Vector)
    * @param VariableKey The variable ID to set
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    Vector GetValue(
        const IndexType VariableKey,
        Vector& rValue
        )
    {
        rValue = mMapVectorVariables[VariableKey];
        return rValue;
    }

    /**
    * It return a value from the map (Matrix)
    * @param VariableKey The variable ID to set
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    Matrix GetValue(
        const IndexType VariableKey,
        Matrix& rValue
        )
    {
        rValue = mMapMatrixVariables[VariableKey];
        return rValue;
    }

    /**
    * It removes one particular pair from the map
    * @param VariableKey The variable ID to remove
    */
    void RemoveId(const IndexType VariableKey)
    {
        auto double_set = mMapDoubleVariables.find(VariableKey);
        if(double_set != mMapDoubleVariables.end()) {
            mMapDoubleVariables.erase(double_set);
            return void();
        }
        auto array_set = mMapArrayVariables.find(VariableKey);
        if(array_set != mMapArrayVariables.end()) {
            mMapArrayVariables.erase(array_set);
            return void();
        }
        auto vector_set = mMapVectorVariables.find(VariableKey);
        if(vector_set != mMapVectorVariables.end()) {
            mMapVectorVariables.erase(vector_set);
            return void();
        }
        auto matrix_set = mMapMatrixVariables.find(VariableKey);
        if(matrix_set != mMapMatrixVariables.end()) {
            mMapMatrixVariables.erase(matrix_set);
            return void();
        }
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

    ConstitutiveLaw::Pointer mpConstitutiveLaw; /// The constitutive law pointer
    double mWeight;                             /// The integration weight of the GP

    /* For values not available on the constitutive law */                                        /// This is the position on the list of GP inside the element
    std::unordered_map<IndexType,double> mMapDoubleVariables;             /// This maps stores auxiliar doubles to interpolate later
    std::unordered_map<IndexType,array_1d<double, 3>> mMapArrayVariables; /// This maps stores auxiliar arrays to interpolate later
    std::unordered_map<IndexType,Vector> mMapVectorVariables;             /// This maps stores auxiliar vectors to interpolate later
    std::unordered_map<IndexType,Matrix> mMapMatrixVariables;             /// This maps stores auxiliar matrixes to interpolate later

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
