// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

#ifndef FACE_ANGLE_RESPONSE_FUNCTION_UTILITY_H
#define FACE_ANGLE_RESPONSE_FUNCTION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/geometry_utilities.h"

// ==============================================================================

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

/// Short class definition.
/** Detail class definition.

 */

class FaceAngleResponseFunctionUtility
{
public:
	///@name Type Definitions
	///@{

	typedef array_1d<double, 3> array_3d;

	/// Pointer definition of FaceAngleResponseFunctionUtility
	KRATOS_CLASS_POINTER_DEFINITION(FaceAngleResponseFunctionUtility);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	FaceAngleResponseFunctionUtility(ModelPart& rModelPart, Parameters ResponseSettings)
	: mrModelPart(rModelPart)
	{
		const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
		KRATOS_ERROR_IF(domain_size != 3) << "FaceAngleResponseFunctionUtility can only be used on 3D geometries!" << std::endl;

		mMainDirection = ResponseSettings["main_direction"].GetVector();
		const double direction_norm = norm_2(mMainDirection);
		KRATOS_ERROR_IF(domain_size != 3) << "FaceAngleResponseFunctionUtility: 'main_direction' vector norm is 0!" << std::endl;
		mMainDirection /= direction_norm;

		double min_angle = ResponseSettings["min_angle"].GetDouble();
		mSinMinAngle = std::sin(min_angle * M_PI / 180);

		std::string gradient_mode = ResponseSettings["gradient_mode"].GetString();
		if (gradient_mode.compare("finite_differencing") == 0)
		{
			double delta = ResponseSettings["step_size"].GetDouble();
			mDelta = delta;
		}
		else
			KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: finite_differencing" << std::endl;
	}


	/// Destructor.
	virtual ~FaceAngleResponseFunctionUtility()
	{
	}

	///@}
	///@name Operators
	///@{

	///@}
	///@name Operations
	///@{

	void Initialize()
	{}

	double CalculateValue()
	{
		KRATOS_TRY;

		double value = 0.0;

		for (auto& cond_i : mrModelPart.Conditions()){
			const double g_i = CalculateConditionValue(cond_i);
			if (g_i <= 0) {
				continue;
			}
			value += g_i*g_i;
		}

		mValue = sqrt(value);

		return mValue;

		KRATOS_CATCH("");
	}

	void CalculateGradient()
	{
		KRATOS_TRY;
		// First gradients are initialized
		VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());

		for (auto& cond_i : mrModelPart.Conditions()){

			const double g_i = CalculateConditionValue(cond_i);
			if (g_i <= 0) {
				continue;
			}

			// Compute sensitivities using finite differencing in the three spatial direction
			array_3d gradient(3, 0.0);

			for (auto& node_i : cond_i.GetGeometry()){

				// Apply pertubation in X-direction
				double g_i_after_fd = 0.0;
				node_i.X() += mDelta;
				node_i.X0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i);
				gradient[0] = (g_i_after_fd - g_i) / mDelta;
				node_i.X() -= mDelta;
				node_i.X0() -= mDelta;

				// Apply pertubation in Y-direction
				g_i_after_fd = 0.0;
				node_i.Y() += mDelta;
				node_i.Y0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i);
				gradient[1] = (g_i_after_fd - g_i) / mDelta;
				node_i.Y() -= mDelta;
				node_i.Y0() -= mDelta;

				// Apply pertubation in Z-direction
				g_i_after_fd = 0.0;
				node_i.Z() += mDelta;
				node_i.Z0() += mDelta;
				g_i_after_fd = CalculateConditionValue(cond_i);
				gradient[2] = (g_i_after_fd - g_i) / mDelta;
				node_i.Z() -= mDelta;
				node_i.Z0() -= mDelta;

				// Add to aggregated sensitivities
				noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += 1/mValue * g_i * gradient;
			}
		}

		KRATOS_CATCH("");
	}

	// ==============================================================================

	///@name Input and output
	///@{

	/// Turn back information as a string.
	std::string Info() const
	{
		return "FaceAngleResponseFunctionUtility";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream &rOStream) const
	{
		rOStream << "FaceAngleResponseFunctionUtility";
	}

	/// Print object's data.
	virtual void PrintData(std::ostream &rOStream) const
	{
	}

	///@}

protected:

	// ==============================================================================

private:

	///@name Operations
	///@{

	double CalculateConditionValue(Condition& rFace)
	{
		// face normal
		array_3d face_normal;
		array_3d v1;
		array_3d v2;

		//calculate the normal on the given condition
		if (rFace.GetGeometry().PointsNumber() == 3)
			GeometryUtilities::CalculateNormal3DTriangle(rFace,face_normal,v1,v2);
		else if (rFace.GetGeometry().PointsNumber() == 4)
			GeometryUtilities::CalculateNormal3DQuad(rFace,face_normal,v1,v2);
		else
			KRATOS_ERROR << "Calculation of surface normal not implemented for the given surface condition!";

		face_normal /= norm_2(face_normal);

		// positive inner product results in negative nodal value to conform to g_i < 0
		return - ((inner_prod(mMainDirection, face_normal)) - mSinMinAngle);
	}

	///@}
	///@name Static Member Variables
	///@{

	///@}
	///@name Member Variables
	///@{

	ModelPart &mrModelPart;
	double mDelta;
	array_3d mMainDirection;
	double mSinMinAngle;
	double mValue;

	///@}

}; // Class FaceAngleResponseFunctionUtility

///@}

} // namespace Kratos.

#endif // FACE_ANGLE_RESPONSE_FUNCTION_UTILITY_H
