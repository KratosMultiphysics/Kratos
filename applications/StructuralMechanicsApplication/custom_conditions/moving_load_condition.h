// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once 

// System includes

// System includes

// External includes

// Project includes
#include "custom_conditions/base_load_condition.h"

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
 * @class MovingLoadCondition
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the responsible to add the contributions of the RHS and LHS of the moving loads of the structure
 * @details Bending moment and reaction forces are calculated on the nodes of the condition element, following a load on an
 * arbitrary position within the element
 * @tparam TDim The dimension of the condition
 * @tparam TNumNodes The number of nodes of the condition
 * @author Aron Noordam
 */
template<std::size_t TDim, std::size_t TNumNodes>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  MovingLoadCondition
    : public BaseLoadCondition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PointLoadCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(MovingLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MovingLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    MovingLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~MovingLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& ThisNodes
        ) const override;


    /**
     * @brief Check if Rotational Dof existant
     * @return Trues if exists, false otherwise
     */
    bool HasRotDof() const override;


    /**
     * \brief Initializes solution step. It determines wether the moving load reactions are to be calculated
     * \param rCurrentProcessInfo current process info
     */
    void InitializeSolutionStep(const ProcessInfo & rCurrentProcessInfo) override;

    /**
     * \brief Initializes non linear iteration. It calculates the displacement and the rotation at the location of the moving load
     * \param rCurrentProcessInfo current process info
     */
    void InitializeNonLinearIteration(const ProcessInfo & rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Point load Condition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Point load Condition #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

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

    /**
     * This functions calculates both the RHS and the LHS following a moving load
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * \brief Calculates exact shape functions for a local load in normal direction
     * \param rShapeFunctionsVector  vector of exact shape functions
     * \param LocalXCoord local x coordinate within condition element
     */
    void CalculateExactNormalShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const;

    /**
     * \brief Calculates exact shape functions for a local load in perpendicular direction
     * \param rShapeFunctionsVector vector of exact shape functions
     * \param LocalXCoord local x coordinate within condition element
     */
    void CalculateExactShearShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const;

    /**
     * \brief Calculates exact shape functions for a local moment around z-axis.
     * \param rShapeFunctionsVector vector of exact shape functions
     * \param LocalXCoord local x coordinate within condition element
     */
    void CalculateExactRotationalShapeFunctions(VectorType& rShapeFunctionsVector, const double LocalXCoord) const;

    /**
    * \brief Calculates derivatives of exact shape functions for a local load in perpendicular direction
    * \param rShapeFunctionsVector vector of exact shape functions
    * \param LocalXCoord local x coordinate within condition element
    */
    void CalculateExactShearShapeFunctionsDerivatives(VectorType& rShapeFunctionsVector, const double LocalXCoord) const;

    /**
     * \brief Calculates exact shape functions for a local moment around z-axis.
     * \param rShapeFunctionsVector vector of exact shape functions
     * \param LocalXCoord local x coordinate within condition element
     */
    void CalculateExactRotationalShapeFunctionsDerivatives(VectorType& rShapeFunctionsVector, const double LocalXCoord) const;

    /**
     * \brief Calculates rotation matrix 
     * \param rRotationMatrix rotation matrix for current condition element
     * \param rGeom condition element
     */
    void CalculateRotationMatrix(BoundedMatrix<double, TDim, TDim>& rRotationMatrix, const GeometryType& rGeom);


    /**
     * \brief Calculates the global bending moment matrix
     * \param rRotationalShapeFunctionVector shape functions vector for rotation
     * \param rLocalMovingLoad array for the value if the local moving load
     * \return global bending moment matrix
     */
    Matrix CalculateGlobalMomentMatrix(const VectorType& rRotationalShapeFunctionVector, const array_1d<double, TDim>& rLocalMovingLoad) const;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    MovingLoadCondition() {};

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    bool mIsMovingLoad = false;


    ///@}
    ///@name Private Operators
    ///@{

    /**
     * \brief Gets the nodal rotation vector
     * \param rRotationsVector nodal rotation vector
     * \param Step step from which the rotations needs to be retrieved
     */
    void GetRotationsVector(Vector& rRotationsVector, const int Step) const;

    ///@}
    ///@name Private Operations
    ///@{


     /**
     * \brief Calculates displacement vector at the location of the moving load
     */
    Vector CalculateLoadPointDisplacementVector();

    /**
    * \brief Calculates rotation vector at the location of the moving load
    */
    Vector CalculateLoadPointRotationVector();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseLoadCondition );
        rSerializer.save("mIsMovingLoad", mIsMovingLoad);
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseLoadCondition );
        rSerializer.load("mIsMovingLoad", mIsMovingLoad);
    }

    ///@}
    ///@name Un accessible methods
    ///@{
    
    ///@}

}; // Class MovingLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<std::size_t TDim, std::size_t TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
        MovingLoadCondition<TDim, TNumNodes>& rThis);

/// output stream function
template<std::size_t TDim, std::size_t TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
        const MovingLoadCondition<TDim, TNumNodes>& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

}  // namespace Kratos.
