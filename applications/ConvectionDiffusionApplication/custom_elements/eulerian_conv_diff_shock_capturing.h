// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Aniol Sala, Guillermo Casas
//

#if !defined(KRATOS_EULERIAN_CONVECTION_DIFFUSION_SHOCK_CAPTURING_ELEMENT_INCLUDED )
#define  KRATOS_EULERIAN_CONVECTION_DIFFUSION_SHOCK_CAPTURING_ELEMENT_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"



namespace Kratos
{

///formulation described in https://docs.google.com/document/d/13a_zGLj6xORDuLgoOG5LwHI6BwShvfO166opZ815zLY/edit?usp=sharing
template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) EulerianConvectionDiffusionShockCapturingElement
    : public Element
{
public:
    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(EulerianConvectionDiffusionShockCapturingElement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor.

    EulerianConvectionDiffusionShockCapturingElement() : Element()
    {
    }

    EulerianConvectionDiffusionShockCapturingElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EulerianConvectionDiffusionShockCapturingElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EulerianConvectionDiffusionShockCapturingElement() {};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<EulerianConvectionDiffusionShockCapturingElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<EulerianConvectionDiffusionShockCapturingElement>(NewId, pGeom, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::string Info() const override
    {
        return "EulerianConvectionDiffusionShockCapturingElement #";
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << Id();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    struct ElementVariables
    {
        double theta;
        double dyn_st_beta;
        double dt_inv;
        double lumping_factor;
        double conductivity;
        double specific_heat;
        double density;
        double beta;
        double div_v;
        double discontinuity_capturing_constant;
        bool stationary;
        bool use_anisotropic_disc_capturing;

        array_1d<double,TNumNodes> phi;
        array_1d<double,TNumNodes> phi_old;
        array_1d<double, TDim> gradient_of_phi;
        array_1d<double,TNumNodes> volumetric_source;
        array_1d< array_1d<double,3 >, TNumNodes> v;
        array_1d< array_1d<double,3 >, TNumNodes> vold;
        BoundedMatrix<double, TDim, TDim> nonlinear_diffusion_matrix;
    };

    void InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume);

    double ComputeH(BoundedMatrix<double,TNumNodes,TDim>& rDN_DX);

    void GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const;

    double CalculateTau(const ElementVariables& rVariables, double norm_vel, double h);

    void CalculateNonlinearDiffusionMatrix(const ElementVariables& rVariables,
                                           const double h,
                                           const array_1d<double, TDim >& velocity,
                                           BoundedMatrix<double, TDim, TDim>& diffusion_matrix);

    // Member Variables


//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

        // Serialization


    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }


}; // Class EulerianConvectionDiffusionShockCapturingElement

} // namespace Kratos.

#endif // KRATOS_EULERIAN_CONVECTION_DIFFUSION_SHOCK_CAPTURING_ELEMENT_INCLUDED  defined
