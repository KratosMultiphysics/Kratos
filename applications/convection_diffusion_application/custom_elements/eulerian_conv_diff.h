/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */


#if !defined(KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED )
#define  KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED


// System includes


// External includes
#include "boost/smart_ptr.hpp"


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
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) EulerianConvectionDiffusionElement
    : public Element
{
public:
    /// Counted pointer of
    KRATOS_CLASS_POINTER_DEFINITION(EulerianConvectionDiffusionElement);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor.

    EulerianConvectionDiffusionElement() : Element()
    {
    }

    EulerianConvectionDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {}

    EulerianConvectionDiffusionElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~EulerianConvectionDiffusionElement() {};

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Element::Pointer(new EulerianConvectionDiffusionElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

	void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    std::string Info() const override
    {
        return "LevelSetConvectionElementSimplex #";
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

        array_1d<double,TNumNodes> phi;
        array_1d<double,TNumNodes> phi_old;
        array_1d<double,TNumNodes> volumetric_source;
        array_1d< array_1d<double,3 >, TNumNodes> v;
        array_1d< array_1d<double,3 >, TNumNodes> vold;
    };

    void InitializeEulerianElement(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometry(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim>& rDN_DX, double& rVolume);

    double ComputeH(boost::numeric::ublas::bounded_matrix<double,TNumNodes,TDim>& rDN_DX);

    void GetNodalValues(ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo);

    double CalculateTau(const ElementVariables& rVariables, double norm_vel, double h);

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


}; // Class EulerianConvectionDiffusionElement

} // namespace Kratos.

#endif // KRATOS_EULERIAN_CONVECTION_DIFFUSION_ELEMENT_INCLUDED  defined
