// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  
//

#if !defined(KRATOS_CONSERVATIVE_LEVELSET_ELEMENT_INCLUDED )
#define  KRATOS_CONSERVATIVE_LEVELSET_ELEMENT_INCLUDED 


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
#include "geometries/geometry_data.h"

// Application includes
#include "custom_elements/eulerian_conv_diff.h"

namespace Kratos
{
template< unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(CONVECTION_DIFFUSION_APPLICATION) ConservativeLevelsetElement
    : public EulerianConvectionDiffusionElement<TDim, TNumNodes>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConservativeLevelsetElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ConservativeLevelsetElement);

    /// Base convection-diffusion element type
    using BaseType = EulerianConvectionDiffusionElement<TDim, TNumNodes>;

    /// Geometry type
    using GeometryType = typename BaseType::GeometryType;

    /// Properties type
    using PropertiesType = typename BaseType::PropertiesType;

    /// Index type
    using IndexType = typename BaseType::IndexType;

    /// Size type
    using SizeType = typename BaseType::SizeType;

    /// Vector type
    using VectorType = typename BaseType::VectorType;

    /// Matrix type
    using MatrixType = typename BaseType::MatrixType;

    /// Nodes array type
    using NodesArrayType = typename BaseType::NodesArrayType;

    /// Shape functions gradient container type
    using ShapeFunctionsGradientsType = typename GeometryType::ShapeFunctionsGradientsType;

    ///@}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default constructor.

    ConservativeLevelsetElement() : BaseType()
    {
    }

    ConservativeLevelsetElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {}

    ConservativeLevelsetElement(IndexType NewId, typename GeometryType::Pointer pGeometry,
      typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {}

    /// Destructor.
    virtual ~ConservativeLevelsetElement(){};

    ///@name Operations
    ///@{

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ConservativeLevelsetElement<TDim, TNumNodes>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(
        IndexType NewId,
        typename GeometryType::Pointer pGeom,
        typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ConservativeLevelsetElement<TDim, TNumNodes>>(NewId, pGeom, pProperties);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

};

}// namespace Kratos.

#endif // KRATOS_CONSERVATIVE_LEVELSET_ELEMENT_INCLUDED