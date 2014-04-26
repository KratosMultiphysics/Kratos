#include "werner_wengle_wall_condition.h"

namespace Kratos {

///@name Specialized implementation for functions that depend on TDim
///@{

/**
 * @see WernerWengleWallCondition::EquationIdVector
 */
template<>
void WernerWengleWallCondition<2, 2>::EquationIdVector(
		EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
	{
		const unsigned int NumNodes = 2;
		const unsigned int LocalSize = 4;
		unsigned int LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
		{
			rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(
					VELOCITY_X).EquationId();
			rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(
					VELOCITY_Y).EquationId();
		}
	}
	else if (this->Is(INTERFACE) && rCurrentProcessInfo[FRACTIONAL_STEP] == 5)
	{
		const SizeType NumNodes = 2;
		const SizeType LocalSize = 2;
		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rResult[LocalIndex++] =
					this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
		}
	}
	else
	{
		rResult.resize(0, false);
	}
}

/**
 * @see WernerWengleWallCondition::EquationIdVector
 */
template<>
void WernerWengleWallCondition<3, 3>::EquationIdVector(
		EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 9;
		unsigned int LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
		{
			rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(
					VELOCITY_X).EquationId();
			rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(
					VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(
					VELOCITY_Z).EquationId();
		}
	}
	else if (this->Is(INTERFACE) && rCurrentProcessInfo[FRACTIONAL_STEP] == 5)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 3;
		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rResult[LocalIndex++] =
					this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
		}
	}
	else
	{
		rResult.resize(0, false);
	}
}

/**
 * @see WernerWengleWallCondition::GetDofList
 */
template<>
void WernerWengleWallCondition<2, 2>::GetDofList(
		DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
	{
		const SizeType NumNodes = 2;
		const SizeType LocalSize = 4;

		if (rConditionDofList.size() != LocalSize)
			rConditionDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
		}
	}
	else if (this->Is(INTERFACE) && rCurrentProcessInfo[FRACTIONAL_STEP] == 5)
	{
		const SizeType NumNodes = 2;
		const SizeType LocalSize = 2;

		if (rConditionDofList.size() != LocalSize)
			rConditionDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(PRESSURE);
		}
	}
	else
	{
		rConditionDofList.resize(0);
	}
}

/**
 * @see WernerWengleWallCondition::GetDofList
 */
template<>
void WernerWengleWallCondition<3, 3>::GetDofList(
		DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 9;

		if (rConditionDofList.size() != LocalSize)
			rConditionDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
		}
	}
	else if (this->Is(INTERFACE) && rCurrentProcessInfo[FRACTIONAL_STEP] == 5)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 3;

		if (rConditionDofList.size() != LocalSize)
			rConditionDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType iNode = 0; iNode < NumNodes; ++iNode)
		{
			rConditionDofList[LocalIndex++] =
					this->GetGeometry()[iNode].pGetDof(PRESSURE);
		}
	}
	else
	{
		rConditionDofList.resize(0);
	}
}

/**
 * @see WernerWengleWallCondition::CalculateWallParameters
 */
template<>
void WernerWengleWallCondition<2, 2>::CalculateWallParameters(
		double& rWallHeight, array_1d<double, 3>& rWallVel, double& rArea)
{
	KRATOS_TRY;
	const double Small = 1.0e-12;
	double DetM, s, w1, Proj;
	array_1d<double, 3> Rhs;
	MatrixType M(3, 3), InvM(3, 3);
	const array_1d<double, 3>& Normal = this->GetValue(NORMAL);
	GeometryType& rElemGeom = pGetElement()->GetGeometry();
	const GeometriesArrayType& rElemFaces = rElemGeom.Faces();
	const array_1d<double, 3>& Center = this->GetGeometry().Center();
	const array_1d<double, 3> t1 = rElemGeom[1].Coordinates()
			- rElemGeom[0].Coordinates();
	const array_1d<double, 3> t2 = rElemGeom[2].Coordinates()
			- rElemGeom[0].Coordinates();

	// generate a dummy base vector orthogonal to 2D domain
	M(0, 1) = t1[1] * t2[2] - t1[2] * t2[1];
	M(1, 1) = t1[2] * t2[0] - t1[0] * t2[2];
	M(2, 1) = t1[0] * t2[1] - t1[1] * t2[0];

	rWallHeight = 0.0;
	rArea = norm_2(Normal);
	for (SizeType i = 0; i < rElemFaces.size(); i++)
	{
		const GeometryType& rFace = rElemFaces[i];

		// rFace[0] + w1*(rFace[1] - rFace[0]) + w2*(Dummy Vector) = Center - s*Normal
		M(0, 0) = rFace[1].X() - rFace[0].X();
		M(1, 0) = rFace[1].Y() - rFace[0].Y();
		M(2, 0) = rFace[1].Z() - rFace[0].Z();
		M(0, 2) = Normal[0];
		M(1, 2) = Normal[1];
		M(2, 2) = Normal[2];

		if (fabs(MathUtils<double>::Det3(M)) < Small * pow(mMinEdgeLength, 4))
		{
			continue;
		}

		Rhs = Center - rFace[0].Coordinates();

		MathUtils<double>::InvertMatrix3(M, InvM, DetM);
		w1 = InvM(0, 0) * Rhs[0] + InvM(0, 1) * Rhs[1] + InvM(0, 2) * Rhs[2];
		s = InvM(2, 0) * Rhs[0] + InvM(2, 1) * Rhs[1] + InvM(2, 2) * Rhs[2];
		if (w1 >= -Small && w1 <= 1.0 + Small) // check if normal intersects this face
		{
			// rWallHeight = ||s*Normal|| = |s| * ||Normal|| = |s| * rArea
			rWallHeight = 2.0 * fabs(s) * rArea;
			if (rWallHeight > Small * mMinEdgeLength) // don't count condition's face
			{
				const array_1d<double, 3> v0 =
						rFace[0].FastGetSolutionStepValue(VELOCITY, 1)
								- rFace[0].FastGetSolutionStepValue(
										MESH_VELOCITY, 1);
				const array_1d<double, 3> v1 =
						rFace[1].FastGetSolutionStepValue(VELOCITY, 1)
								- rFace[1].FastGetSolutionStepValue(
										MESH_VELOCITY, 1);

				rWallVel[0] = w1 * v1[0] + (1.0 - w1) * v0[0];
				rWallVel[1] = w1 * v1[1] + (1.0 - w1) * v0[1];
				rWallVel[2] = w1 * v1[2] + (1.0 - w1) * v0[2];

				// make velocity tangent
				Proj = (rWallVel[0] * Normal[0] + rWallVel[1] * Normal[1]
						+ rWallVel[2] * Normal[2]) / (rArea * rArea);
				rWallVel[0] -= Proj * Normal[0];
				rWallVel[1] -= Proj * Normal[1];
				rWallVel[2] -= Proj * Normal[2];

				break;
			}
		}
	}

	KRATOS_CATCH("");
}

/**
 * @see WernerWengleWallCondition::CalculateWallParameters
 */
template<>
void WernerWengleWallCondition<3, 3>::CalculateWallParameters(
		double& rWallHeight, array_1d<double, 3>& rWallVel, double& rArea)
{
	KRATOS_TRY;

	const double Small = 1.0e-12;
	double DetM, s, w1, w2, Proj;
	array_1d<double, 3> Rhs;
	MatrixType M(3, 3), InvM(3, 3);
	ElementPointerType pElem = pGetElement();
	const array_1d<double, 3>& Normal = this->GetValue(NORMAL);
	const GeometriesArrayType& rElemFaces = pElem->GetGeometry().Faces();
	const array_1d<double, 3>& center = this->GetGeometry().Center();

	rWallHeight = 0.0;
	rArea = norm_2(Normal);
	for (SizeType i = 0; i < rElemFaces.size(); i++)
	{
		const GeometryType& rFace = rElemFaces[i];

		// rFace[0] + w1*(rFace[1] - rFace[0]) + w2*(rFace[2] - rFace[0]) = center - s*Normal
		M(0, 0) = rFace[1].X() - rFace[0].X();
		M(1, 0) = rFace[1].Y() - rFace[0].Y();
		M(2, 0) = rFace[1].Z() - rFace[0].Z();
		M(0, 1) = rFace[2].X() - rFace[0].X();
		M(1, 1) = rFace[2].Y() - rFace[0].Y();
		M(2, 1) = rFace[2].Z() - rFace[0].Z();
		M(0, 2) = Normal[0];
		M(1, 2) = Normal[1];
		M(2, 2) = Normal[2];

		if (fabs(MathUtils<double>::Det3(M)) < Small * pow(mMinEdgeLength, 4))
			continue;

		Rhs = center - rFace[0].Coordinates();

		MathUtils<double>::InvertMatrix3(M, InvM, DetM);
		w1 = InvM(0, 0) * Rhs[0] + InvM(0, 1) * Rhs[1] + InvM(0, 2) * Rhs[2];
		w2 = InvM(1, 0) * Rhs[0] + InvM(1, 1) * Rhs[1] + InvM(1, 2) * Rhs[2];
		s = InvM(2, 0) * Rhs[0] + InvM(2, 1) * Rhs[1] + InvM(2, 2) * Rhs[2];
		if (w1 >= -Small && w2 >= -Small && (w1 + w2) <= 1. + Small) // check if normal intersects this face
		{
			// rWallHeight = ||s*Normal|| = |s| * ||Normal|| = |s| * rArea
			rWallHeight = 2.0 * fabs(s) * rArea;
			if (rWallHeight > Small * mMinEdgeLength) // don't count condition's face
			{
				const array_1d<double, 3> v0 =
						rFace[0].FastGetSolutionStepValue(VELOCITY, 1)
								- rFace[0].FastGetSolutionStepValue(
										MESH_VELOCITY, 1);
				const array_1d<double, 3> v1 =
						rFace[1].FastGetSolutionStepValue(VELOCITY, 1)
								- rFace[1].FastGetSolutionStepValue(
										MESH_VELOCITY, 1);
				const array_1d<double, 3> v2 =
						rFace[2].FastGetSolutionStepValue(VELOCITY, 1)
								- rFace[2].FastGetSolutionStepValue(
										MESH_VELOCITY, 1);

				rWallVel[0] = w1 * v1[0] + w2 * v2[0] + (1.0 - w1 - w2) * v0[0];
				rWallVel[1] = w1 * v1[1] + w2 * v2[1] + (1.0 - w1 - w2) * v0[1];
				rWallVel[2] = w1 * v1[2] + w2 * v2[2] + (1.0 - w1 - w2) * v0[2];

				// make velocity tangent
				Proj = (rWallVel[0] * Normal[0] + rWallVel[1] * Normal[1]
						+ rWallVel[2] * Normal[2]) / (rArea * rArea);
				rWallVel[0] -= Proj * Normal[0];
				rWallVel[1] -= Proj * Normal[1];
				rWallVel[2] -= Proj * Normal[2];

				break;
			}
		}
	}

	KRATOS_CATCH("");
}

} // namespace Kratos
