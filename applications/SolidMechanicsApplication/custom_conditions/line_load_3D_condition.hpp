//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_LOAD_3D_CONDITION_H_INCLUDED )
#define  KRATOS_LINE_LOAD_3D_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/force_load_condition.hpp"

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

/// Force Load Condition for 3D and 2D geometries. (base class)

/**
 * Implements a Force Load definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D (base class)
 */
class KRATOS_API(SOLID_MECHANICS_APPLICATION) LineLoad3DCondition
    : public ForceLoadCondition
{
public:

    ///@name Type Definitions
    ///@{
    // Counted pointer of LineLoad3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( LineLoad3DCondition );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LineLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry );

    LineLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties );

    /// Copy constructor
    LineLoad3DCondition( LineLoad3DCondition const& rOther);

    /// Destructor
    virtual ~LineLoad3DCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * creates a new condition pointer
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(IndexType NewId,
			      NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties ) const;


    /**
     * clones the selected condition variables, creating a new one
     * @param NewId: the ID of the new condition
     * @param ThisNodes: the nodes of the new condition
     * @param pProperties: the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone(IndexType NewId, 
			     NodesArrayType const& ThisNodes) const;



    //************* COMPUTING  METHODS


    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo );

    ///@}
    ///@name Access
    ///@{
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    ///@name Friends
    ///@{
    ///@}

	//calculate load direction for Follower Force

	void CalculateFollowerForceDirection(Vector& rVectorForce);
	bool mIsFirstStep = true;
	double mCosAngle;

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    LineLoad3DCondition() {};
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Initialize System Matrices
     */
    virtual void InitializeGeneralVariables(GeneralVariables& rVariables, 
					    const ProcessInfo& rCurrentProcessInfo);

    /**
     * Calculate Condition Kinematics
     */
    virtual void CalculateKinematics(GeneralVariables& rVariables, 
				     const double& rPointNumber);

    /**
     * Calculation of the Position Increment
     */
    virtual Matrix& CalculateDeltaPosition(Matrix & rDeltaPosition);

    /**
     * Calculation of the Vector Force of the Condition
     */
    virtual Vector& CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables);


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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

	//follower force variables
	Vector mQuaternionVEC_A, mQuaternionVEC_B;
	double mQuaternionSCA_A, mQuaternionSCA_B;
	Vector mDirectionVectorXOriginal, mDirectionVectorYOriginal,
		mDirectionVectorZOriginal;

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


}; // class LineLoad3DCondition.


class OrientationGeometry
{
public:
	OrientationGeometry(array_1d<double, 3>& v1, const double theta = 0.00);
	OrientationGeometry(array_1d<double, 3>& v1, array_1d<double, 3>& v2);


	void CalculateRotationMatrix(bounded_matrix<double, 3, 3>& R);
	void CalculateBasisVectors(array_1d<double, 3>& v1,
		array_1d<double, 3>& v2,
		array_1d<double, 3>& v3);

	Quaternion<double>& GetQuaternion() { return morientation; }

private:
	Quaternion<double> morientation;
};

} // namespace Kratos.

#endif // KRATOS_LINE_LOAD_3D_CONDITION_H_INCLUDED defined 
