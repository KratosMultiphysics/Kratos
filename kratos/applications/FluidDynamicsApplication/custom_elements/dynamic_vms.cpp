#include "dynamic_vms.h"
#include "utilities/math_utils.h"

namespace Kratos
{
    ///@name Specialized implementation of DynamicVMS for functions that depend on TDim
    ///@{

    /**
     * @see DynamicVMS::GetValueOnIntegrationPoints
     */
    template <>
    void DynamicVMS<2>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                                     std::vector<array_1d<double,3> >& rOutput,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int Dim(2),NumNodes(3);
        if (rVariable == VORTICITY)
        {
            // Set output vector (for a single integration point)
            rOutput.resize(1);
            array_1d<double,3> & rVorticity = rOutput[0];
            rVorticity[0] = 0.0;
            rVorticity[1] = 0.0;
            rVorticity[2] = 0.0;

            double Area;
            array_1d<double, NumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
                rVorticity[2] += N[iNode] * ( DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0] );
            }
            rVorticity *= 0.5; // vorticity = 1/2 (nabla x velocity)
        }
        else if (rVariable == SUBSCALE)
        {
            rOutput.resize(1);
            array_1d<double,3> & rSubscale = rOutput[0];
            rSubscale[0] = mSubscaleVel[0];
            rSubscale[1] = mSubscaleVel[1];
            rSubscale[2] = 0.0;
        }
        else // Default behaviour (returns elemental data)
        {
            rOutput.resize(1);
            /*
             The cast is done to avoid modification of the element's data. Data modification
             would happen if rVariable is not stored now (would initialize a pointer to &rVariable
             with associated value of 0.0). This is catastrophic if the variable referenced
             goes out of scope.
             */
            const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
            rOutput[0] = const_this->GetValue(rVariable);
        }
    }

    /**
     * @see DynamicVMS::GetValueOnIntegrationPoints
     */
    template <>
    void DynamicVMS<3>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                                     std::vector<array_1d<double,3> >& rOutput,
                                                     const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int Dim(3),NumNodes(4);
        if (rVariable == VORTICITY)
        {
            // Set output vector (for a single integration point)
            rOutput.resize(1);
            array_1d<double, 3 > & rVorticity = rOutput[0];
            rVorticity[0] = 0.0;
            rVorticity[1] = 0.0;
            rVorticity[2] = 0.0;

            double Area;
            array_1d<double, NumNodes> N;
            boost::numeric::ublas::bounded_matrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
                rVorticity[0] += N[iNode] * ( DN_DX(iNode,1)*rVelocity[2] - DN_DX(iNode,2)*rVelocity[1] );
                rVorticity[1] += N[iNode] * ( DN_DX(iNode,2)*rVelocity[0] - DN_DX(iNode,0)*rVelocity[2] );
                rVorticity[2] += N[iNode] * ( DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0] );
            }
            rVorticity *= 0.5; // vorticity = 1/2 (nabla x velocity)
        }
        else if (rVariable == SUBSCALE)
        {
            rOutput.resize(1);
            array_1d<double,3> & rSubscale = rOutput[0];
            rSubscale[0] = mSubscaleVel[0];
            rSubscale[1] = mSubscaleVel[1];
            rSubscale[2] = mSubscaleVel[2];
        }
        else // Default behaviour (returns elemental data)
        {
            rOutput.resize(1);
            /*
             The cast is done to avoid modification of the element's data. Data modification
             would happen if rVariable is not stored now (would initialize a pointer to &rVariable
             with associated value of 0.0). This is catastrophic if the variable referenced
             goes out of scope.
             */
            const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
            rOutput[0] = const_this->GetValue(rVariable);
        }
    }


    template <>
    void DynamicVMS<2,3>::DenseSystemSolve(const MatrixType& rA,
                                           VectorType& rx,
                                           const VectorType& rb)
    {
        MatrixType Inverse(2,2);
        double Det;
        MathUtils<double>::InvertMatrix2(rA,Inverse,Det);

        noalias(rx) = prod(Inverse,rb); // rx = Inverse * rb
    }

    template <>
    void DynamicVMS<3,4>::DenseSystemSolve(const MatrixType& rA,
                                           VectorType& rx,
                                           const VectorType& rb)
    {
        MatrixType Inverse(3,3);
        double Det;
        MathUtils<double>::InvertMatrix3(rA,Inverse,Det);

        noalias(rx) = prod(Inverse,rb); // rx = Inverse * rb
    }

    template <>
    void DynamicVMS<2,3>::save(Serializer& rSerializer) const
    {
        rSerializer.save("Name","DynamicVMS2D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElementType );
    }

    template <>
    void DynamicVMS<3,4>::save(Serializer& rSerializer) const
    {
        rSerializer.save("Name","DynamicVMS3D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElementType );
    }

    ///@} // Specialized implementations
}
