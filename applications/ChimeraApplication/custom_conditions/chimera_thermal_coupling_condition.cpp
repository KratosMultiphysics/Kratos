//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Sonja Schneider
//                   Jordi Cotela
//

// From S. Schneider Implementation of a Chimera Technique TUM Master Thesis 2015

#include "chimera_thermal_coupling_condition.h"
#include "chimera_application_variables.h"

namespace Kratos
{

template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::ChimeraThermalCouplingCondition(IndexType NewId):
    Condition(NewId)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    mGPFlux = std::vector< double >(NumGauss,0.0);
}


template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::ChimeraThermalCouplingCondition(IndexType NewId, const NodesArrayType &ThisNodes):
    Condition(NewId, ThisNodes)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    mGPFlux = std::vector< double >(NumGauss,0.0);
}


template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::ChimeraThermalCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry):
    Condition(NewId, pGeometry)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    mGPFlux = std::vector< double >(NumGauss,0.0);
}


template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::ChimeraThermalCouplingCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Condition(NewId, pGeometry, pProperties)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    mGPFlux = std::vector< double >(NumGauss,0.0);
}


template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::ChimeraThermalCouplingCondition(const ChimeraThermalCouplingCondition<TDim> &rOther):
    Condition(rOther)
{
    this->mGPFlux = rOther.mGPFlux;
}


template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>::~ChimeraThermalCouplingCondition()
{}

///////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
ChimeraThermalCouplingCondition<TDim>& ChimeraThermalCouplingCondition<TDim>::operator=(const ChimeraThermalCouplingCondition<TDim> &rOther)
{
    Condition::operator=(rOther);
    this->mGPFlux = rOther.mGPFlux;
    return *this;
}


template< unsigned int TDim >
Condition::Pointer ChimeraThermalCouplingCondition<TDim>::Create(IndexType NewId, const NodesArrayType &ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ChimeraThermalCouplingCondition(NewId,this->GetGeometry().Create(ThisNodes), pProperties));
}


///////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::Initialize()
{
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                                                                 VectorType &rRightHandSideVector,
                                                                 ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Initialize local contributions
    const SizeType LocalSize = this->GetGeometry().PointsNumber();

    if (rLeftHandSideMatrix.size1() != LocalSize)
    {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize);
    }

    if (rRightHandSideVector.size() != LocalSize)
    {
        rRightHandSideVector.resize(LocalSize);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    this->ImposeFlux(rLeftHandSideMatrix, rRightHandSideVector);

    KRATOS_CATCH("");
}

template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix, ProcessInfo &rCurrentProcessInfo)
{
    VectorType TempRHS;
    this->CalculateLocalSystem(rLeftHandSideMatrix,TempRHS,rCurrentProcessInfo);
}

template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::CalculateRightHandSide(VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    MatrixType TempLHS;
    this->CalculateLocalSystem(TempLHS,rRightHandSideVector,rCurrentProcessInfo);
}




template< unsigned int TDim >
int ChimeraThermalCouplingCondition<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area >= 0

    if (Check != 0)
    {
        return Check;
    }
    else
    {
        // Check that all required variables have been registered
        if(TEMPERATURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"TEMPERATURE Key is 0. Check if the application was correctly registered.","");
        if(FLUX.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"FLUX Key is 0. Check if the application was correctly registered.","");
        if(NORMAL.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"NORMAL Key is 0. Check if the application was correctly registered.","");

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(TEMPERATURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing TEMPERATURE variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(FLUX) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing FLUX variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(NORMAL) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data for node ",this->GetGeometry()[i].Id());

            if(this->GetGeometry()[i].HasDofFor(TEMPERATURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing TEMPERATURE component degree of freedom on node ",this->GetGeometry()[i].Id());
        }

        return Check;
    }

    KRATOS_CATCH("");
}


///////////////////////////////////////////////////////////////////////////////

template <const unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::EquationIdVector(EquationIdVectorType& rResult,
                                                             ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int LocalSize = this->GetGeometry().PointsNumber();
    unsigned int LocalIndex = 0;

    unsigned int tpos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < LocalSize; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(TEMPERATURE,tpos).EquationId();
    }
}


///////////////////////////////////////////////////////////////////////////////


template <const unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetDofList(DofsVectorType& rElementalDofList,
                                                       ProcessInfo& rCurrentProcessInfo)
{
    const SizeType LocalSize = this->GetGeometry().PointsNumber();

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < LocalSize; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(TEMPERATURE);
    }
}

///////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
GeometryData::IntegrationMethod ChimeraThermalCouplingCondition<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////



template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::GetValuesVector(Vector &rValues, int Step)
{
    const SizeType LocalSize = this->GetGeometry().PointsNumber();
    unsigned int LocalIndex = 0;

    if (rValues.size() != LocalSize)
    {
        rValues.resize(LocalSize, false);
    }

    for (unsigned int i = 0; i < LocalSize; ++i)
    {
        double Temp = this->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE, Step);
        rValues[LocalIndex++] = Temp;
    }
}

///////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::Calculate(const Variable< array_1d<double,3> >& rVariable,
                                                      array_1d<double,3>& Output,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == NORMAL)
        this->CalculateNormal(Output);
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::SetValueOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable,
                                                                        std::vector<array_1d<double,3> > rValues,
                                                                        const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);

    if (rValues.size() != NumGauss)
    {
        std::stringstream Msg;
        Msg << "ChimeraThermalCouplingCondition error: Wrong number of Gauss point values" << std::endl;
        Msg << rValues.size() << " values recieved for a condition with " << NumGauss << " integration points." << std::endl;

        KRATOS_THROW_ERROR(std::logic_error,Msg.str(),"");
    }
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::SetValueOnIntegrationPoints(const Variable< double > &rVariable,
                                                                        std::vector< double >& rValues,
                                                                        const ProcessInfo &rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);

    if (rValues.size() != NumGauss)
    {
        std::stringstream Msg;
        Msg << "ChimeraThermalCouplingCondition error: Wrong number of Gauss point values" << std::endl;
        Msg << rValues.size() << " values recieved for a condition with " << NumGauss << " integration points." << std::endl;

        KRATOS_THROW_ERROR(std::logic_error,Msg.str(),"");
    }

    if (rVariable == FLUX)
    {
        mGPFlux = rValues;
    }
}


///////////////////////////////////////////////////////////////////////////////


template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable,
                                                                        std::vector<array_1d<double,3> > &rValues,
                                                                        const ProcessInfo &rCurrentProcessInfo
                                                                        )
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    rValues.resize(NumGauss);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
        for (unsigned int g = 1; g < NumGauss; g++)
            rValues[g] = rValues[0];
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ChimeraThermalCouplingCondition* const_this = static_cast< const ChimeraThermalCouplingCondition* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }
}


///////////////////////////////////////////////////////////////////////////////

template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                                        std::vector<double>& rValues,
                                                                        const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int NumGauss = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
    rValues.resize(NumGauss);
    if (rVariable == FLUX)
    {
        for (unsigned int g = 0; g < NumGauss; g++)
            rValues[g] = mGPFlux[g];
    }
    else
    {
        /*
         * The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const ChimeraThermalCouplingCondition* const_this = static_cast< const ChimeraThermalCouplingCondition* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }
}


template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                                      std::vector<array_1d<double, 6 > >& rValues,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const ChimeraThermalCouplingCondition* const_this = static_cast< const ChimeraThermalCouplingCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                      std::vector<Vector>& rValues,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const ChimeraThermalCouplingCondition* const_this = static_cast< const ChimeraThermalCouplingCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                                      std::vector<Matrix>& rValues,
                                                                      const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const ChimeraThermalCouplingCondition* const_this = static_cast< const ChimeraThermalCouplingCondition* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


///////////////////////////////////////////////////////////////////////////////



template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                                                         std::vector<double>& rOutput,
                                                                         const ProcessInfo& rCurrentProcessInfo)
{
    this->GetValueOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
}

template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                                                         std::vector< array_1d<double, 3 > >& Output,
                                                                         const ProcessInfo& rCurrentProcessInfo)
{
    this->GetValueOnIntegrationPoints(rVariable,Output,rCurrentProcessInfo);
}


template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
                                                                       std::vector< Vector >& Output,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{

}

template<unsigned int TDim>
void ChimeraThermalCouplingCondition<TDim>::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                                                                       std::vector< Matrix >& Output,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{

}

///////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
std::string ChimeraThermalCouplingCondition<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "ChimeraThermalCouplingCondition" << TDim << "D";
    return buffer.str();
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "ChimeraThermalCouplingCondition";
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::PrintData(std::ostream &rOStream) const
{}



///////////////////////////////////////////////////////////////////////////////
// protected methods
///////////////////////////////////////////////////////////////////////////////



template <>
void ChimeraThermalCouplingCondition<2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

template <>
void ChimeraThermalCouplingCondition<3>::CalculateNormal(array_1d<double,3>& An )
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


///////////////////////////////////////////////////////////////////////////////

template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::ImposeFlux(MatrixType &rLocalMatrix,
                                                       VectorType &rLocalVector)
{
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();

    MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    array_1d<double,3> Normal;
    this->CalculateNormal(Normal); //this already contains the area
    double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
    //Normal /= A;

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    double J = (TDim == 2) ? 0.5*A : 2.0*A;

    // Traction term is the sum on gauss points of weight_gauss * N_i * N_j * traction
    // (i,j nodal indices)
    for (unsigned int g = 0; g < NumGauss; g++)
    {
        // Shape functions
        Vector N = row(NContainer,g);
        // Integration point weight
        double Weight = J * IntegrationPoints[g].Weight();

        /*
        // Interpolate the traction at the integration point from nodal values
        array_1d<double,3> TGauss = N[0] * rGeom[0].FastGetSolutionStepValue(TRACTION);
        for (unsigned int j = 1; j < NumNodes; j++)
            TGauss += N[j]*rGeom[j].FastGetSolutionStepValue(TRACTION);
        */

        // Read the traction at integration point
        const double& FluxGauss = mGPFlux[g];


//        KRATOS_WATCH(g)
//                KRATOS_WATCH(J)
//                KRATOS_WATCH(Weight)
//        KRATOS_WATCH(mGPFlux[g])

        for (unsigned int i = 0; i < NumNodes; i++)
            rLocalVector[i] += Weight * N[i] * FluxGauss;
    }
}


///////////////////////////////////////////////////////////////////////////////
// private methods
///////////////////////////////////////////////////////////////////////////////


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    rSerializer.save("mGPFlux",mGPFlux);
}


template< unsigned int TDim >
void ChimeraThermalCouplingCondition<TDim>::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    rSerializer.load("mGPFlux",mGPFlux);
}



///////////////////////////////////////////////////////////////////////////////
// template class instantiation
///////////////////////////////////////////////////////////////////////////////

template class ChimeraThermalCouplingCondition<2>;
template class ChimeraThermalCouplingCondition<3>;

} // namespace Kratos
