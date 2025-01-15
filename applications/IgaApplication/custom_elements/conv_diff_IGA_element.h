// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

#if !defined(KRATOS_CONV_DIFF_IGA_ELEMENT_INCLUDED )
#define  KRATOS_CONV_DIFF_IGA_ELEMENT_ELEMENT_INCLUDED


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

class KRATOS_API(IGA_APPLICATION) ConvDiffIGAElement
    : public Element
{
public:
    /// Counted pointer of
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConvDiffIGAElement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor.

    ConvDiffIGAElement() : Element()
    {
    }

    ConvDiffIGAElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    ConvDiffIGAElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ConvDiffIGAElement() {};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes, 
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<ConvDiffIGAElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }
    
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override
    {
        return Kratos::make_intrusive<ConvDiffIGAElement>(NewId, pGeom, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::string Info() const override
    {
        return "ConvDiffIGAElement #";
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

        array_1d<double,1> phi;
        array_1d<double,1> phi_old;
        array_1d<double,1> volumetric_source;
        array_1d< array_1d<double,3 >, 1> v;
        array_1d< array_1d<double,3 >, 1> vold;
    };

    void InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(Matrix& rDN_DX, double& rVolume);

    // double ComputeH(Matrix& rDN_DX);

    void GetNodalValues(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo) const;

    double CalculateTau(const ElementVariables& rVariables, double norm_vel, double h);

    double CalculateTauHigherOrder(const ElementVariables& rVariables, double norm_vel, double h, int p);

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

    Parameters ReadParamatersFile(const std::string& rDataFileName) const;


}; // Class ConvDiffIGAElement

} // namespace Kratos.

#endif // KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED  defined
