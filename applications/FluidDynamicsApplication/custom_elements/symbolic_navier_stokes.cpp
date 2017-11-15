#include "symbolic_navier_stokes.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class SymbolicNavierStokes< SymbolicNavierStokesData<2,3> >;
template class SymbolicNavierStokes< SymbolicNavierStokesData<3,4> >;

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template <class TElementData>
SymbolicNavierStokes<TElementData>::SymbolicNavierStokes(IndexType NewId)
    : FluidElement<TElementData>(NewId) {}

template <class TElementData>
SymbolicNavierStokes<TElementData>::SymbolicNavierStokes(
    IndexType NewId, const NodesArrayType& ThisNodes)
    : FluidElement<TElementData>(NewId, ThisNodes) {}

template <class TElementData>
SymbolicNavierStokes<TElementData>::SymbolicNavierStokes(
    IndexType NewId, GeometryType::Pointer pGeometry)
    : FluidElement<TElementData>(NewId, pGeometry) {}

template <class TElementData>
SymbolicNavierStokes<TElementData>::SymbolicNavierStokes(IndexType NewId,
    GeometryType::Pointer pGeometry, Properties::Pointer pProperties)
    : FluidElement<TElementData>(NewId, pGeometry, pProperties) {}

template <class TElementData>
SymbolicNavierStokes<TElementData>::~SymbolicNavierStokes() {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TElementData >
Element::Pointer SymbolicNavierStokes<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Element::Pointer(new SymbolicNavierStokes(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}


template< class TElementData >
Element::Pointer SymbolicNavierStokes<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Element::Pointer(new SymbolicNavierStokes(NewId, pGeom, pProperties));
}

template <class TElementData>
void SymbolicNavierStokes<TElementData>::Initialize() {
    KRATOS_TRY

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();
        if (r_properties.Has(CONSTITUTIVE_LAW)) {
            mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();
            mpConstitutiveLaw->InitializeMaterial(r_properties,
            this->GetGeometry(),
            row(this->GetGeometry().ShapeFunctionsValues(), 0));
        }
        else {
            KRATOS_ERROR << "No constitutive law provided in the element's properties for " << this->Info();
        }
    }

    KRATOS_CATCH("")
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template <class TElementData>
int SymbolicNavierStokes<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) {
    
    KRATOS_TRY;
    FluidElement<TElementData>::Check(rCurrentProcessInfo);

    // Check constitutive law
    KRATOS_ERROR_IF( mpConstitutiveLaw == nullptr ) << "No consitutive law defined. Please call the element's Initialize() method before attempting the check." << std::endl;

    mpConstitutiveLaw->Check(this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo);

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template <class TElementData>
std::string SymbolicNavierStokes<TElementData>::Info() const {
    std::stringstream buffer;
    buffer << "SymbolicNavierStokes" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void SymbolicNavierStokes<TElementData>::PrintInfo(
    std::ostream& rOStream) const {
    rOStream << this->Info() << std::endl;

    if (this->mpConstitutiveLaw != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->mpConstitutiveLaw->PrintInfo(rOStream);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected operations

template <class TElementData>
void SymbolicNavierStokes<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) {}

template <class TElementData>
void SymbolicNavierStokes<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) {}

template <class TElementData>
void SymbolicNavierStokes<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< class TElementData >
void SymbolicNavierStokes<TElementData>::save(Serializer& rSerializer) const
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType );
    rSerializer.save("mpConstitutiveLaw",mpConstitutiveLaw);
}


template< class TElementData >
void SymbolicNavierStokes<TElementData>::load(Serializer& rSerializer)
{
    using BaseType = FluidElement<TElementData>;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
    rSerializer.load("mpConstitutiveLaw",mpConstitutiveLaw);
}



}