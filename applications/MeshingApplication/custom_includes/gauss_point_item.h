// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
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
    * @brief It checks if an ID exists in the map
    * @param rVariable The  variable to be check
    * @return If the ID already exists or not
    */
    bool Has(const Variable<double>& rVariable)
    {
        const IndexType variable_key = rVariable.Key();
        auto double_set = mMapDoubleVariables.find(variable_key);
        if(double_set != mMapDoubleVariables.end()) {
            return true;
        }

        return false;
    }

    /**
    * @brief It checks if an ID exists in the map
    * @param rVariable The  variable to be check
    * @return If the ID already exists or not
    */
    bool Has(const Variable<array_1d<double, 3>>& rVariable)
    {
        const IndexType variable_key = rVariable.Key();
        auto array_set = mMapArrayVariables.find(variable_key);
        if(array_set != mMapArrayVariables.end()) {
            return true;
        }

        return false;
    }

    /**
    * @brief It checks if an ID exists in the map
    * @param rVariable The  variable to be check
    * @return If the ID already exists or not
    */
    bool Has(const Variable<Vector>& rVariable)
    {
        const IndexType variable_key = rVariable.Key();
        auto vector_set = mMapVectorVariables.find(variable_key);
        if(vector_set != mMapVectorVariables.end()) {
            return true;
        }

        return false;
    }

    /**
    * @brief It checks if an ID exists in the map
    * @param rVariable The  variable to be check
    * @return If the ID already exists or not
    */
    bool Has(const Variable<Matrix>& rVariable)
    {
        const IndexType variable_key = rVariable.Key();
        auto matrix_set = mMapMatrixVariables.find(variable_key);
        if(matrix_set != mMapMatrixVariables.end()) {
            return true;
        }

        return false;
    }

    /**
    * @brief It adds a new value to the map (double)
    * @param rVariable The variable being set
    * @param rValue The value to assign
    */
    void SetValue(
        const Variable<double>& rVariable,
        const double rValue
        )
    {
        const IndexType variable_key = rVariable.Key();
        auto double_set = mMapDoubleVariables.find(variable_key);
        if(double_set != mMapDoubleVariables.end()) {
            mMapDoubleVariables[variable_key] = rValue;
            return void();
        }
        mMapDoubleVariables.insert({variable_key, rValue});
    }

    /**
    * @brief It adds a new value to the map (array_1d<double, 3>)
    * @param rVariable The variable being set
    * @param rValue The value to assign
    */
    void SetValue(
        const Variable<array_1d<double, 3>>& rVariable,
        const array_1d<double, 3>& rValue
        )
    {
        const IndexType variable_key = rVariable.Key();
        auto array_set = mMapArrayVariables.find(variable_key);
        if(array_set != mMapArrayVariables.end()) {
            mMapArrayVariables[variable_key] = rValue;
            return void();
        }
        mMapArrayVariables.insert({variable_key, rValue});
    }

    /**
    * @brief It adds a new value to the map (Vector)
    * @param rVariable The variable being set
    * @param rValue The value to assign
    */
    void SetValue(
        const Variable<Vector>& rVariable,
        const Vector& rValue
        )
    {
        const IndexType variable_key = rVariable.Key();
        auto vector_set = mMapVectorVariables.find(variable_key);
        if(vector_set != mMapVectorVariables.end()) {
            mMapVectorVariables[variable_key] = rValue;
            return void();
        }
        mMapVectorVariables.insert({variable_key, rValue});
    }

    /**
    * @brief It adds a new value to the map (Matrix)
    * @param rVariable The variable being set
    * @param rValue The value to assign
    */
    void SetValue(
        const Variable<Matrix>& rVariable,
        const Matrix& rValue
        )
    {
        const IndexType variable_key = rVariable.Key();
        auto matrix_set = mMapMatrixVariables.find(variable_key);
        if(matrix_set != mMapMatrixVariables.end()) {
            mMapMatrixVariables[variable_key] = rValue;
            return void();
        }
        mMapMatrixVariables.insert({variable_key, rValue});
    }

    /**
    * @brief It return a value from the map (double)
    * @param rVariable The variable being recovered
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    double GetValue(
        const Variable<double>& rVariable,
        double& rValue
        )
    {
        rValue = mMapDoubleVariables[rVariable.Key()];
        return rValue;
    }

    /**
    * @brief It return a value from the map (array_1d<double, 3>)
    * @param rVariable The variable being recovered
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    array_1d<double, 3> GetValue(
        const Variable<array_1d<double, 3>>& rVariable,
        array_1d<double, 3>& rValue
        )
    {
        rValue = mMapArrayVariables[rVariable.Key()];
        return rValue;
    }

    /**
    * @brief It return a value from the map (Vector)
    * @param rVariable The variable being recovered
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    Vector GetValue(
        const Variable<Vector>& rVariable,
        Vector& rValue
        )
    {
        rValue = mMapVectorVariables[rVariable.Key()];
        return rValue;
    }

    /**
    * @brief It return a value from the map (Matrix)
    * @param rVariable The variable being recovered
    * @param rValue The value to recover
    * @return rValue The value to recover
    */
    Matrix GetValue(
        const Variable<Matrix>& rVariable,
        Matrix& rValue
        )
    {
        rValue = mMapMatrixVariables[rVariable.Key()];
        return rValue;
    }

    /**
    * @brief It removes one particular pair from the map
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
