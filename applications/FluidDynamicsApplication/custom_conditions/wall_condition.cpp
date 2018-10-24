#include "wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see WallCondition::EquationIdVector
 */
template <>
void WallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
	int step = r_process_info[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const unsigned int NumNodes = 2;
        const unsigned int LocalSize = 4;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        }
    }
    else
    {
        		if(this->Is(INTERFACE) && step==5 )
        {
                //add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const SizeType NumNodes = 2;

                if (rResult.size() != NumNodes)
					rResult.resize(NumNodes, false);

                unsigned int LocalIndex = 0;

                for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
					rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
        else
        {
			rResult.resize(0,false);
		}
    }
}

/**
 * @see WallCondition::EquationIdVector
 */
template <>
void WallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
	int step = r_process_info[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        }
    }
    else
    {
		if(this->Is(INTERFACE) && step==5 )
        {
                //add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const SizeType NumNodes = 3;

                if (rResult.size() != NumNodes)
					rResult.resize(NumNodes, false);

                unsigned int LocalIndex = 0;

                for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
					rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
        else
        {
			rResult.resize(0,false);
		}
    }
}

/**
 * @see WallCondition::GetDofList
 */
template <>
void WallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                              ProcessInfo& rCurrentProcessInfo)
{
 	const ProcessInfo& r_process_info = rCurrentProcessInfo;
	int step = r_process_info[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 2;
        const SizeType LocalSize = 4;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        }
    }
    else
    {
                if(this->Is(INTERFACE) && step==5 )
        {
                //add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const SizeType NumNodes = 2;

                if (rElementalDofList.size() != NumNodes)
                    rElementalDofList.resize(NumNodes);

                unsigned int LocalIndex = 0;

                for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
        else
        {
            rElementalDofList.resize(0);
        }
    }
}

/**
 * @see WallCondition::GetDofList
 */
template <>
void WallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                    ProcessInfo& rCurrentProcessInfo)
{
    const ProcessInfo& r_process_info = rCurrentProcessInfo;
	int step = r_process_info[FRACTIONAL_STEP];
    if ( step == 1 )
    {
        const SizeType NumNodes = 3;
        const SizeType LocalSize = 9;

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        }
    }
    else
    {
        if(this->Is(INTERFACE) && step==5 )
        {
                //add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
                const SizeType NumNodes = 3;

                if (rElementalDofList.size() != NumNodes)
                    rElementalDofList.resize(NumNodes);

                unsigned int LocalIndex = 0;

                for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
                    rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
        else
        {
            rElementalDofList.resize(0);
        }
    }
}

template <>
 void WallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
    {
        Geometry<Node<3> >& pGeometry = this->GetGeometry();

        An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
        An[1] = - (pGeometry[1].X() - pGeometry[0].X());
        An[2] =    0.00;

    }

template <>
void WallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
    {
        Geometry<Node<3> >& pGeometry = this->GetGeometry();

        array_1d<double,3> v1,v2;
        v1[0] = pGeometry[1].X() - pGeometry[0].X();
        v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
        v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

        v2[0] = pGeometry[2].X() - pGeometry[0].X();
        v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
        v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

        MathUtils<double>::CrossProduct(An,v1,v2);
        An *= 0.5;
    }

template<unsigned int TDim, unsigned int TNumNodes>
void WallCondition<TDim,TNumNodes>::ApplyNeumannCondition(MatrixType &rLocalMatrix, VectorType &rLocalVector)
{
    const WallCondition<TDim,TNumNodes>& rConstThis = *this;
    if (rConstThis.GetValue(IS_STRUCTURE) == 0.0)
    {
        const unsigned int LocalSize = TDim;
        const GeometryType& rGeom = this->GetGeometry();
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        const unsigned int NumGauss = IntegrationPoints.size();
        Vector GaussWeights = ZeroVector(NumGauss);

        MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        array_1d<double,3> Normal;
        this->CalculateNormal(Normal); //this already contains the area
        double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
        Normal /= A;

        for (unsigned int g = 0; g < NumGauss; g++)
            GaussWeights[g] = 2.0*A * IntegrationPoints[g].Weight();

        for (unsigned int g = 0; g < NumGauss; g++)
        {
            Vector N = row(NContainer,g);
            double Weight = GaussWeights[g];

            // Neumann boundary condition
            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                //unsigned int row = i*LocalSize;
                const NodeType& rConstNode = this->GetGeometry()[i];
                if ( rConstNode.IsFixed(PRESSURE)==true)
                {
                    const double pext = rConstNode.FastGetSolutionStepValue(PRESSURE);
                    for (unsigned int j = 0; j < TNumNodes; j++)
                    {
                        unsigned int row = j*LocalSize;
                        for (unsigned int d = 0; d < TDim;d++)
                            rLocalVector[row+d] -= Weight*N[j]*N[i]*pext*Normal[d];
                    }
                }
            }

            // Velocity inflow correction
            array_1d<double,3> Vel = ZeroVector(3);
            double Density = 0.0;

            for (unsigned int i = 0; i < TNumNodes; i++)
            {
                const NodeType& rConstNode = this->GetGeometry()[i];
                Vel += N[i]*rConstNode.FastGetSolutionStepValue(VELOCITY);
                Density += N[i]*rConstNode.FastGetSolutionStepValue(DENSITY);
            }

            double Proj = Vel[0]*Normal[0] + Vel[1]*Normal[1] + Vel[2]*Normal[2];

            if (Proj < 0)
            {
                const double W = Weight*Density*Proj;
                for (unsigned int i = 0; i < TNumNodes; i++)
                {
                    double row = i*LocalSize;
                    for (unsigned int j = 0; j < TNumNodes; j++)
                    {
                        double col = j*LocalSize;
                        const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                        for (unsigned int d = 0; d < TDim; d++)
                        {
                            double Tij = W*N[i]*N[j];
                            rLocalMatrix(row+d,col+d) -= Tij;
                            rLocalVector[row+d] += Tij*rVel[d];
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class WallCondition<2,2>;
template class WallCondition<3,3>;

} // namespace Kratos
