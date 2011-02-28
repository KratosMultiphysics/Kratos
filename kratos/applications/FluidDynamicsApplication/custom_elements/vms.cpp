#include "vms.h"

namespace Kratos
{
    ///@name Specialized implementation of VMS for functions that depend on TDim
    ///@{

    /**
     @see VMS::EquationIdVector
     */
    template <>
    void VMS<2>::EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int NumNodes(3),LocalSize(9);
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
    }

    /**
     @see VMS::EquationIdVector
     */
    template <>
    void VMS<3>::EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int NumNodes(4),LocalSize(16);
        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
        }
    }

    /**
     @see VMS::GetDofList
     */
    template <>
    void VMS<2>::GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int NumNodes(3),LocalSize(9);
        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
    }

    /**
     @see VMS::GetDofList
     */
    template <>
    void VMS<3>::GetDofList(DofsVectorType& rElementalDofList,
                            ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int NumNodes(4),LocalSize(16);
        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
    }

    /**
     @see VMS::GetFirstDerivativesVector
     */
    template <>
    void VMS<2>::GetFirstDerivativesVector(Vector& Values, int Step)
    {
        const unsigned int NumNodes(3),LocalSize(9);
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
            Values[LocalIndex++] = rVelocity[0];
            Values[LocalIndex++] = rVelocity[1];
            Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
        }
    }

    /**
     @see VMS::GetFirstDerivativesVector
     */
    template <>
    void VMS<3>::GetFirstDerivativesVector(Vector& Values, int Step)
    {
        const unsigned int NumNodes(4),LocalSize(16);
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
            Values[LocalIndex++] = rVelocity[0];
            Values[LocalIndex++] = rVelocity[1];
            Values[LocalIndex++] = rVelocity[2];
            Values[LocalIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE, Step);
        }
    }

    /**
     @see VMS::GetSecondDerivativesVector
     */
    template <>
    void VMS<2>::GetSecondDerivativesVector(Vector& Values, int Step)
    {
        const unsigned int NumNodes(3),LocalSize(9);
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
            Values[LocalIndex++] = rAcceleration[0];
            Values[LocalIndex++] = rAcceleration[1];
            Values[LocalIndex++] = 0.0; // Pressure Dof
        }
    }

    /**
     @see VMS::GetSecondDerivativesVector
     */
    template <>
    void VMS<3>::GetSecondDerivativesVector(Vector& Values, int Step)
    {
        const unsigned int NumNodes(4),LocalSize(16);
        unsigned int LocalIndex = 0;

        if (Values.size() != LocalSize)
            Values.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            array_1d<double,3>& rAcceleration = this->GetGeometry()[iNode].FastGetSolutionStepValue(ACCELERATION, Step);
            Values[LocalIndex++] = rAcceleration[0];
            Values[LocalIndex++] = rAcceleration[1];
            Values[LocalIndex++] = rAcceleration[2];
            Values[LocalIndex++] = 0.0; // Pressure Dof
        }
    }

    /**
     @see VMS::GetValueOnIntegrationPoints
     */
    template <>
    void VMS<2>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                              std::vector<array_1d<double,3> >& rOutput,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int Dim(2),NumNodes(3);
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
                rVorticity[2] += N[iNode] * ( DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0] );
            }
            rVorticity *= 0.5; // vorticity = 1/2 (nabla x velocity)
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
     @see VMS::GetValueOnIntegrationPoints
     */
    template <>
    void VMS<3>::GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
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
     The size of the 2D element is estimated as the diameter of a circle of the same area.
     Area = Pi * (h/2)^2
     @see VMS::ElementSize
     */
    template <>
    double VMS<2>::ElementSize(const double Area)
    {
        return 1.128379167 * sqrt(Area); //Diameter of circumference of given Area
    }

    /**
     The size of the 3D element is estimated as the diameter of the sphere
     circumscribed to a regular tetrahedron with the same volume.
     @see VMS::ElementSize
     */
    template <>
    double VMS<3>::ElementSize(const double Volume)
    {
        return 0.60046878 * pow(Volume,0.333333333333333333333);
    }

    /**
     Returns the squared element size, estimated as h^2 = 2*Area
     @see VMS::FilterWidth
     */
    template <>
    double VMS<2>::FilterWidth()
    {
        double FilterWidth = GeometryUtils::CalculateVolume2D(this->GetGeometry());
        return 2.0 * FilterWidth;
    }

    /**
     Returns the squared element size, estimated from the assumption V = (1/6) * h^3
     @see VMS::FilterWidth
     */
    template <>
    double VMS<3>::FilterWidth()
    {
        double FilterWidth = GeometryUtils::CalculateVolume3D(this->GetGeometry());
        FilterWidth *= 6.0;
        return pow(FilterWidth, 2.0/3.0);
    }

    template <>
    void VMS<2>::save(Serializer& rSerializer) const
    {
        rSerializer.save("Name","VMS2D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    template <>
    void VMS<3>::save(Serializer& rSerializer) const
    {
        rSerializer.save("Name","VMS3D");
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    }

    ///@} // Specialized implementations
}
