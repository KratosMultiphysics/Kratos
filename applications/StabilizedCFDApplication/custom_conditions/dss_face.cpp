//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "custom_conditions/dss_face.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::DSSFace(DSSFace::IndexType NewId):
    Condition(NewId)
{}

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::DSSFace(DSSFace::IndexType NewId,
                 const DSSFace::NodesArrayType &ThisNodes):
    Condition(NewId,ThisNodes)
{}

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::DSSFace(DSSFace::IndexType NewId,
                                 GeometryType::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::DSSFace(DSSFace::IndexType NewId,
                                 GeometryType::Pointer pGeometry,
                                 PropertiesType::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::DSSFace(const DSSFace &rOther):
    Condition(rOther)
{
}

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>::~DSSFace()
{
}

///////////////////////////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNumNodes>
DSSFace<TDim,TNumNodes>& DSSFace<TDim,TNumNodes>::operator=(const Kratos::DSSFace<TDim,TNumNodes> &rOther)
{
    Condition::operator=(rOther);

    return *this;
}

///////////////////////////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer DSSFace<TDim,TNumNodes>::Create(IndexType NewId,
                                                   NodesArrayType const& ThisNodes,
                                                   PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DSSFace(NewId, GetGeometry().Create(ThisNodes), pProperties));
}



template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                    ProcessInfo& rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalSystem(rLeftHandSideMatrix,RHS,rCurrentProcessInfo);
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                   VectorType& rRightHandSideVector,
                                                   ProcessInfo& rCurrentProcessInfo)
{
    // Initialize local contributions
    const SizeType LocalSize = (TDim+1) * TNumNodes;

    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
}

template< unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
                                                                 VectorType& rRightHandSideVector,
                                                                 ProcessInfo& rCurrentProcessInfo)
{
    // Initialize local contributions
    const SizeType LocalSize = (TDim+1) * TNumNodes;

    if (rDampMatrix.size1() != LocalSize)
        rDampMatrix.resize(LocalSize,LocalSize,false);
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    this->AddBoundaryTerms(rDampMatrix,rRightHandSideVector,rCurrentProcessInfo);
}


template<unsigned int TDim, unsigned int TNumNodes>
int DSSFace<TDim,TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    if (Check != 0)
    {
        return Check;
    }
    else
    {
        // Check that all required variables have been registered
        if(VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check if the application was correctly registered.","");
        if(MESH_VELOCITY.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check if the application was correctly registered.","");
        if(PRESSURE.Key() == 0)
            KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check if the application was correctly registered.","");

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {

            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        }

        return Check;
    }

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////
// Instanced functions, depending on dimension and node number
///////////////////////////////////////////////////////////////////////////////

template <>
void DSSFace<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int LocalSize = 6;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < 2; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template <>
void DSSFace<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int LocalSize = 12;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < 3; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}

template <>
void DSSFace<3,4>::EquationIdVector(EquationIdVectorType& rResult,
                                    ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int LocalSize = 16;

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (unsigned int i = 0; i < 4; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE,ppos).EquationId();
    }
}


template <>
void DSSFace<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                              ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}


template <>
void DSSFace<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                              ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = 3;
    const unsigned int LocalSize = 4*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}


template <>
void DSSFace<3,4>::GetDofList(DofsVectorType& rElementalDofList,
                              ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = 4;
    const unsigned int LocalSize = 4*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}

///////////////////////////////////////////////////////////////////////////////


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValuesVector(Vector& Values,
                                              int Step)
{
    const SizeType LocalSize = TDim * TNumNodes;
    unsigned int LocalIndex = 0;

    if (Values.size() != LocalSize)
        Values.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        array_1d<double,3>& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step);
        for (unsigned int d = 0; d < TDim; ++d)
            Values[LocalIndex++] = rVelocity[d];
    }
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> > &rVariable,
                                                          std::vector<array_1d<double,3> > &rValues,
                                                          const ProcessInfo &rCurrentProcessInfo)
{
    rValues.resize(1);
    if (rVariable == NORMAL)
    {
        this->CalculateNormal(rValues[0]);
    }
    else
    {
        /* The cast is done to avoid modification of the element's data. Data modification
         * would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         * with associated value of 0.0). This is catastrophic if the variable referenced
         * goes out of scope.
         */
        const DSSFace* const_this = static_cast< const DSSFace* >(this);
        rValues[0] = const_this->GetValue(rVariable);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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
    const DSSFace* const_this = static_cast< const DSSFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                          std::vector<array_1d<double, 6 > >& rValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const DSSFace* const_this = static_cast< const DSSFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                          std::vector<Vector>& rValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const DSSFace* const_this = static_cast< const DSSFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                          std::vector<Matrix>& rValues,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    rValues.resize(1);
    const DSSFace* const_this = static_cast< const DSSFace* >(this);
    rValues[0] = const_this->GetValue(rVariable);
}

template<unsigned int TDim, unsigned int TNumNodes>
std::string DSSFace<TDim,TNumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "DSSFace" << TDim << "D";
    return buffer.str();
}

template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << "DSSFace";
}


template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::PrintData(std::ostream &rOStream) const
{
    this->PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////
// protected functions
///////////////////////////////////////////////////////////////////////////////


template <>
void DSSFace<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.00;

}

template <>
void DSSFace<3,3>::CalculateNormal(array_1d<double,3>& An )
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

template <>
void DSSFace<3,4>::CalculateNormal(array_1d<double,3>& An )
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
    //An *= 0.5;
}

template <>
double DSSFace<2,2>::CalculateJacobian(double Area) const
{
    return 0.5*Area;
}

template <>
double DSSFace<3,3>::CalculateJacobian(double Area) const
{
    return 2.0*Area;
}

template <>
double DSSFace<3,4>::CalculateJacobian(double Area) const
{
    return 0.25*Area;
}



template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::AddBoundaryTerms(MatrixType& rLocalMatrix,
                                               VectorType& rLocalVector,
                                               const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    const unsigned int BlockSize = TDim+1;

    MatrixType NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

    array_1d<double,3> Normal;
    this->CalculateNormal(Normal); //this already contains the area
    double A = std::sqrt(Normal[0]*Normal[0]+Normal[1]*Normal[1]+Normal[2]*Normal[2]);
    Normal /= A;

    // CAUTION: "Jacobian" is 2.0*A for triangles but 0.5*A for lines
    //double J = (TDim == 2) ? 0.5*A : 2.0*A;
    double J = this->CalculateJacobian(A);

    for (unsigned int g = 0; g < NumGauss; g++)
    {
        double GaussWeight = J * IntegrationPoints[g].Weight();
        Vector N = row(NContainer,g);

//        // Neumann boundary condition
//        for (unsigned int i = 0; i < TNumNodes; i++)
//        {
//            //unsigned int row = i*LocalSize;
//            const NodeType& rConstNode = this->GetGeometry()[i];
//            if ( rConstNode.IsFixed(PRESSURE)==true)
//            {
//                const double pext = rConstNode.FastGetSolutionStepValue(PRESSURE);
//                for (unsigned int j = 0; j < TNumNodes; j++)
//                {
//                    unsigned int row = j*BlockSize;
//                    for (unsigned int d = 0; d < TDim;d++)
//                        rLocalVector[row+d] -= GaussWeight*N[j]*N[i]*pext*Normal[d];
//                }
//            }
//        }

        // Velocity boundary term (for skew-symmetric formulation)
        array_1d<double,3> ConvVel(3,0.0);
        this->ResolvedConvectiveVelocity(ConvVel,N);
        double Density = 0.0;
        this->EvaluateInPoint(Density,DENSITY,N);

        double Proj = ConvVel[0]*Normal[0] + ConvVel[1]*Normal[1] + ConvVel[2]*Normal[2];

        const double W = 0.5*GaussWeight*Density*Proj; // 0.5 is for skew-symmetric formulation
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            double row = i*BlockSize;
            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                double col = j*BlockSize;
                const array_1d<double,3>& rVel = this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                double Tij = W*N[i]*N[j];
                for (unsigned int d = 0; d < TDim; d++)
                {
                    rLocalMatrix(row+d,col+d) += Tij;
                    rLocalVector[row+d] -= Tij*rVel[d];
                }
            }
        }
    }

    /*
    Vector Values;
    this->GetValuesVector(Values);
    noalias(rLocalVector) -= prod(rLocalMatrix,rLocalVector);*/
}

template< unsigned int TDim, unsigned int TNumNodes >
void DSSFace<TDim,TNumNodes>::ResolvedConvectiveVelocity(array_1d<double,3> &rConvVel, const ShapeFunctionsType &rN)
{
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    array_1d<double,3> NodeVel = rGeom[0].FastGetSolutionStepValue(VELOCITY);
    NodeVel -= rGeom[0].FastGetSolutionStepValue(MESH_VELOCITY);
    rConvVel = rN[0] * NodeVel;

    for (unsigned int i = 1; i < NumNodes; i++)
    {
        NodeVel = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        NodeVel -= rGeom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        rConvVel += rN[i] * NodeVel;
    }
}

///////////////////////////////////////////////////////////////////////////////

template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
}

template<unsigned int TDim, unsigned int TNumNodes>
void DSSFace<TDim,TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
}



///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class DSSFace<2,2>;
template class DSSFace<3,3>;
template class DSSFace<3,4>;

}
