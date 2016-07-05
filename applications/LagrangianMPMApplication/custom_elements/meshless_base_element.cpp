///    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//




// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/meshless_base_element.h"

#include "utilities/math_utils.h"

#include "geometries/geometry.h"
//#include "custom_geometries/meshless_geometry.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessBaseElement::MeshlessBaseElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
MeshlessBaseElement::MeshlessBaseElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer MeshlessBaseElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new MeshlessBaseElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

MeshlessBaseElement::~MeshlessBaseElement()
{
}

//************************************************************************************
//************************************************************************************
void MeshlessBaseElement::GetGeometryData(double& integration_weight,
        Vector& N,
        Matrix& DN_Dx
                                         )
{
    integration_weight = this->GetValue(GAUSS_AREA);
    N = this->GetValue(SHAPE_FUNCTIONS);
    DN_Dx = this->GetValue(SHAPE_FUNCTIONS_DERIVATIVES);
}


//************************************************************************************
//************************************************************************************
} // Namespace Kratos


