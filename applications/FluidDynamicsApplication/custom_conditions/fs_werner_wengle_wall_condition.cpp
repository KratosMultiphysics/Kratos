#include "fs_werner_wengle_wall_condition.h"

namespace Kratos {

///@name Specialized implementation for functions that depend on TDim
///@{

/**
 * @see FSWernerWengleWallCondition::EquationIdVector
 */
template<>
void FSWernerWengleWallCondition<2, 2>::EquationIdVector(
		EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const 
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
 * @see FSWernerWengleWallCondition::EquationIdVector
 */
template<>
void FSWernerWengleWallCondition<3, 3>::EquationIdVector(
		EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const 
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
 * @see FSWernerWengleWallCondition::GetDofList
 */
template<>
void FSWernerWengleWallCondition<2, 2>::GetDofList(
		DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const 
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
 * @see FSWernerWengleWallCondition::GetDofList
 */
template<>
void FSWernerWengleWallCondition<3, 3>::GetDofList(
		DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const 
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



template<unsigned int TDim, unsigned int TNumNodes>
void FSWernerWengleWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<array_1d<double,3> > &rVariable,
        std::vector<array_1d<double,3> > &rValues,
        const ProcessInfo &rCurrentProcessInfo)
{
    rValues.resize(1);
    /* The cast is done to avoid modification of the element's data. Data modification
     * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
     * with associated value of 0.0). This is catastrophic if the variable referenced
     * goes out of scope.
     */
    const FSWernerWengleWallCondition* const_this = static_cast< const FSWernerWengleWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void FSWernerWengleWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    /*
       The cast is done to avoid modification of the element's data. Data modification
       would happen if rVariable is not stored now (would initialize a pointer to &rVariable
       with associated value of 0.0). This is catastrophic if the variable referenced
       goes out of scope.
       */
    const FSWernerWengleWallCondition* const_this = static_cast< const FSWernerWengleWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void FSWernerWengleWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<array_1d<double, 6 > >& rVariable,
        std::vector<array_1d<double, 6 > >& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const FSWernerWengleWallCondition* const_this = static_cast< const FSWernerWengleWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void FSWernerWengleWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const FSWernerWengleWallCondition* const_this = static_cast< const FSWernerWengleWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void FSWernerWengleWallCondition<TDim,TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const FSWernerWengleWallCondition* const_this = static_cast< const FSWernerWengleWallCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}




/**
 * @see FSWernerWengleWallCondition::CalculateWallParameters
 */
template<>
void FSWernerWengleWallCondition<2, 2>::CalculateWallParameters(
		double& rWallHeight, array_1d<double, 3>& rWallVel, double& rArea)
{
	KRATOS_TRY;

	const double Small = 1.0e-12;
        double DetM, s, w1, Proj;
        array_1d<double, 3> Rhs;
        MatrixType M(2, 2), InvM(2, 2);
        ElementPointerType pElem = pGetElement();
        const array_1d<double, 3>& Normal = this->GetValue(NORMAL);
        GeometryType& rElemGeom = pElem->GetGeometry();
        const GeometriesArrayType& edges = rElemGeom.GenerateEdges();
        const array_1d<double, 3>& center = this->GetGeometry().Center();

        rWallHeight = 0.0;
        rArea = norm_2(Normal);
	for (SizeType i = 0; i < edges.size(); i++)
	  {
	    const GeometryType& rEdge = edges[i];

	    // rEdge[0] + w1*(rEdge[1] - rEdge[0]) = center - s*Normal
	    M(0, 0) = rEdge[1].X() - rEdge[0].X();
	    M(1, 0) = rEdge[1].Y() - rEdge[0].Y();
	    M(0, 1) = Normal[0];
	    M(1, 1) = Normal[1];

	    if (fabs(MathUtils<double>::Det2(M)) < Small * pow(mMinEdgeLength, 2))
	      {
		continue;
	      }

	    Rhs = center - rEdge[0].Coordinates();

	    MathUtils<double>::InvertMatrix2(M, InvM, DetM);
	    w1 = InvM(0, 0) * Rhs[0] + InvM(0, 1) * Rhs[1];
	    s  = InvM(1, 0) * Rhs[0] + InvM(1, 1) * Rhs[1];
	    if (w1 >= -Small && w1 <= 1.0 + Small) // check if normal intersects this edge
	      {
		// rWallHeight = ||s*Normal|| = |s| * ||Normal|| = |s| * rArea
		rWallHeight = fabs(s) * rArea;
		if (rWallHeight > Small * mMinEdgeLength) // don't count condition's face
		  {
		    const array_1d<double, 3> v0 =
		      rEdge[0].FastGetSolutionStepValue(VELOCITY, 1)
		      - rEdge[0].FastGetSolutionStepValue(MESH_VELOCITY, 1);
		    const array_1d<double, 3> v1 =
		      rEdge[1].FastGetSolutionStepValue(VELOCITY, 1)
		      - rEdge[1].FastGetSolutionStepValue(MESH_VELOCITY, 1);

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
 * @see FSWernerWengleWallCondition::CalculateWallParameters
 */
template<>
void FSWernerWengleWallCondition<3, 3>::CalculateWallParameters(
		double& rWallHeight, array_1d<double, 3>& rWallVel, double& rArea)
{
	KRATOS_TRY;

	const double Small = 1.0e-12;
	double DetM, s, w1, w2, Proj;
	array_1d<double, 3> Rhs;
	MatrixType M(3, 3), InvM(3, 3);
	ElementPointerType pElem = pGetElement();
	const array_1d<double, 3>& Normal = this->GetValue(NORMAL);
	const GeometriesArrayType& rElemFaces = pElem->GetGeometry().GenerateFaces();
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

template class FSWernerWengleWallCondition<2,2>;
template class FSWernerWengleWallCondition<3,3>;

} // namespace Kratos
