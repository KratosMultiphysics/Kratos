//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#include "d_vms.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
//#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId):
    FluidElement<TElementData>(NewId)
{}

template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, const NodesArrayType& ThisNodes):
    FluidElement<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, GeometryType::Pointer pGeometry):
    FluidElement<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
DVMS<TElementData>::DVMS(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
    FluidElement<TElementData>(NewId,pGeometry,pProperties)
{}


template< class TElementData >
DVMS<TElementData>::~DVMS()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer DVMS<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<DVMS>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TElementData >
Element::Pointer DVMS<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<DVMS>(NewId, pGeom, pProperties);
}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<double>& rVariable,
    double& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMS<TElementData>::Calculate(
    const Variable<array_1d<double, 3>>& rVariable,
    array_1d<double, 3>& rOutput, const ProcessInfo& rCurrentProcessInfo) {
    // Lumped projection terms
    if (rVariable == ADVPROJ) {
        this->CalculateProjections(rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<Vector>& rVariable,
    Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void DVMS<TElementData>::Calculate(const Variable<Matrix>& rVariable,
    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int DVMS<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int out = FluidElement<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Extra variables
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    // Output variables (for Calculate() functions)
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_PRESSURE);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMS<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_VELOCITY) {

    }
    else {
        FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}


template< class TElementData >
void DVMS<TElementData>::GetValueOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    if (rVariable == SUBSCALE_PRESSURE) {
    }
    else {
        FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

template <class TElementData>
void DVMS<TElementData>::GetValueOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::GetValueOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

template <class TElementData>
void DVMS<TElementData>::GetValueOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    FluidElement<TElementData>::GetValueOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template< class TElementData >
std::string DVMS<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "DVMS #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void DVMS<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "DVMS" << Dim << "D";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions

template< class TElementData >
void DVMS<TElementData>::ASGSMomentumResidual(
    TElementData& rData,
    array_1d<double,3> &rMomentumRes)
{

}


template< class TElementData >
void DVMS<TElementData>::ASGSMassResidual(
    TElementData& rData,
    double &rMomentumRes)
{

}


template< class TElementData >
void DVMS<TElementData>::OSSMomentumResidual(
    TElementData& rData,
    array_1d<double,3> &rMomentumRes)
{

}


template< class TElementData >
void DVMS<TElementData>::OSSMassResidual(
    TElementData& rData,
    double &rMassRes)
{

}


template< class TElementData >
void DVMS<TElementData>::MomentumProjTerm(
    TElementData& rData,
    array_1d<double,3> &rMomentumRHS)
{

}


template< class TElementData >
void DVMS<TElementData>::MassProjTerm(
    TElementData& rData,
    double &rMassRHS)
{

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Evaluation of system terms on Gauss Points

template <class TElementData>
void DVMS<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {

    // Call specialized implementation (it is on a helper class to avoid partial template specialization problems)
    Internals::SpecializedAddTimeIntegratedSystem<TElementData,
        TElementData::ElementManagesTimeIntegration>::AddSystem(this, rData,
        rLHS, rRHS);
}

template <class TElementData>
void DVMS<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {

    KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl;
}

template <class TElementData>
void DVMS<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {

    KRATOS_ERROR << "AddTimeIntegratedRHS is not implemented." << std::endl;
}

template< class TElementData >
void DVMS<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType &rLocalLHS,
    VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMS<TElementData>::AddMassLHS(
    TElementData& rData,
    MatrixType &rMassMatrix)
{
    const double density = this->GetAtCoordinate(rData.Density,rData.N);

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;
            const double Mij = rData.Weight * density * rData.N[i] * rData.N[j];
            for (unsigned int d = 0; d < Dim; d++)
                rMassMatrix(row+d,col+d) += Mij;
        }
    }

    /* Note on OSS and full projection: Riccardo says that adding the terms provided by
     * AddMassStabilization (and incluiding their corresponding terms in the projeciton)
     * could help reduce the non-linearity of the coupling between projection and u,p
     * However, leaving them on gives a lot of trouble whith the Bossak scheme:
     * think that we solve F - (1-alpha)*M*u^(n+1) - alpha*M*u^(n) - K(u^(n+1)) = 0
     * so the projection of the dynamic terms should be Pi( (1-alpha)*u^(n+1) - alpha*u^(n) )
     */
    if ( rData.UseOSS != 1.0 )
        this->AddMassStabilization(rData,rMassMatrix);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMS<TElementData>::AddMassStabilization(
    TElementData& rData,
    MatrixType &rMassMatrix)
{

}

template <class TElementData>
void DVMS<TElementData>::AddBoundaryIntegral(TElementData& rData,
    const Vector& rUnitNormal, MatrixType& rLHS, VectorType& rRHS) {

    boost::numeric::ublas::bounded_matrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;
    
    boost::numeric::ublas::bounded_matrix<double,StrainSize,LocalSize> shear_stress_matrix = boost::numeric::ublas::prod(constitutive_matrix,strain_matrix);

    boost::numeric::ublas::bounded_matrix<double,Dim,StrainSize> normal_projection = ZeroMatrix(Dim,StrainSize);
    FluidElementUtilities<NumNodes>::VoigtTransformForProduct(rUnitNormal,normal_projection);
    
    // Contribution to boundary stress from 2*mu*symmetric_gradient(velocity)*n
    boost::numeric::ublas::bounded_matrix<double,Dim,LocalSize> normal_stress_operator = boost::numeric::ublas::prod(normal_projection,shear_stress_matrix);

    // Contribution to boundary stress from p*n
    for (unsigned int i = 0; i < NumNodes; i++) {
        const double ni = rData.N[i];
        for (unsigned int d = 0; d < Dim; d++) {
            const std::size_t pressure_column = i*BlockSize + Dim;
            normal_stress_operator(d,pressure_column) = -rUnitNormal[d]*ni;
        }
    }

    // RHS: stress computed using current solution
    array_1d<double,Dim> shear_stress = boost::numeric::ublas::prod(normal_projection,rData.ShearStress);
    const double p_gauss = this->GetAtCoordinate(rData.Pressure,rData.N);

    // Add -Ni*normal_stress_operator to the LHS, Ni*current_stress to the RHS
    for (unsigned int i = 0; i < NumNodes; i++) {
        const double wni = rData.Weight*rData.N[i];
        for (unsigned int d = 0; d < Dim; d++) {
            const unsigned int row = i*BlockSize + d;
            for (unsigned int col = 0; col < LocalSize; col++) {
                rLHS(row,col) -= wni*normal_stress_operator(d,col);
            }
            rRHS[row] += wni*(shear_stress[d]-p_gauss*rUnitNormal[d]);
        }
    }
}

template <class TElementData>
void DVMS<TElementData>::AddViscousTerm(
    const TElementData& rData,
    boost::numeric::ublas::bounded_matrix<double,LocalSize,LocalSize>& rLHS,
    VectorType& rRHS) {

    boost::numeric::ublas::bounded_matrix<double,StrainSize,LocalSize> strain_matrix = ZeroMatrix(StrainSize,LocalSize);
    FluidElementUtilities<NumNodes>::GetStrainMatrix(rData.DN_DX,strain_matrix);

    const auto& constitutive_matrix = rData.C;
    boost::numeric::ublas::bounded_matrix<double,StrainSize,LocalSize> shear_stress_matrix = boost::numeric::ublas::prod(constitutive_matrix,strain_matrix);

    // Multiply times integration point weight (I do this here to avoid a temporal in LHS += weight * Bt * C * B)
    strain_matrix *= rData.Weight;

    noalias(rLHS) += boost::numeric::ublas::prod(boost::numeric::ublas::trans(strain_matrix),shear_stress_matrix);
    noalias(rRHS) -= boost::numeric::ublas::prod(boost::numeric::ublas::trans(strain_matrix),rData.ShearStress);
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void DVMS<TElementData>::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{

}

template< class TElementData >
void DVMS<TElementData>::SubscaleVelocity(
    TElementData& rData,
    const ProcessInfo &rProcessInfo,
    array_1d<double,3> &rVelocitySubscale)
{

}

template< class TElementData >
void DVMS<TElementData>::SubscalePressure(
        TElementData& rData,
        const ProcessInfo& rProcessInfo,
        double &rPressureSubscale)
{

}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void DVMS<TElementData>::save(Serializer& rSerializer) const
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void DVMS<TElementData>::load(Serializer& rSerializer)
{
    typedef FluidElement<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Internals
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace Internals {

///////////////////////////////////////////////////////////////////////////////////////////////////
// For Standard data: Time integration is not available
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, false>::AddSystem(
    DVMS<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
    KRATOS_TRY;
    KRATOS_ERROR << "Trying to use time-integrated element functions with a "
                    "data type that does not know previous time step data"
                 << std::endl;
    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Specialized time integration
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
void SpecializedAddTimeIntegratedSystem<TElementData, true>::AddSystem(
    DVMS<TElementData>* pElement, TElementData& rData, Matrix& rLHS,
    Vector& rRHS) {
        Matrix mass_matrix = ZeroMatrix(rLHS.size1(),rLHS.size2());
        Matrix velocity_lhs = ZeroMatrix(rLHS.size1(),rLHS.size2());

        pElement->AddVelocitySystem(rData,velocity_lhs,rRHS);
        pElement->AddMassLHS(rData,mass_matrix);

        noalias(rLHS) += rData.bdf0*mass_matrix + velocity_lhs;
        
        Vector acceleration = ZeroVector(rRHS.size());

        int LocalIndex = 0;
        const auto& r_velocities = rData.Velocity;
        const auto& r_velocities_step1 = rData.Velocity_OldStep1;
        const auto& r_velocities_step2 = rData.Velocity_OldStep2;

        for (unsigned int i = 0; i < TElementData::NumNodes; ++i) {
            for (unsigned int d = 0; d < TElementData::Dim; ++d)  {
                // Velocity Dofs
                acceleration[LocalIndex] = rData.bdf0*r_velocities(i,d);
                acceleration[LocalIndex] += rData.bdf1*r_velocities_step1(i,d);
                acceleration[LocalIndex] += rData.bdf2*r_velocities_step2(i,d);
                ++LocalIndex;
            }
            ++LocalIndex;
        }

        noalias(rRHS) -= prod(mass_matrix,acceleration);
}

} // namespace Internals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class DVMS< QSVMSData<2,3> >;

} // namespace Kratos
