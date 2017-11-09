/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"

#include "custom_python/add_spatial_containers_to_python.h"

#include "spatial_containers/bins_object_static_OCL.h"

#include "spatial_containers/spatial_containers.h"
#include "utilities/spatial_containers_configure.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_2d_3.h"
#include "includes/model_part.h"

namespace Kratos
{

namespace Python
{
void  AddSpatialContainersToPython()
{

    static const std::size_t Dim = 3;

    typedef Kratos::Point 		PointType;
    typedef Kratos::Tetrahedra3D4<PointType>ObjectType;

    typedef PointType*			PtrPointType;
    typedef ObjectType*			PtrObjectType;

    typedef PtrPointType*			PointVector;
    typedef PtrPointType*			PointIterator;

    typedef PtrObjectType*			ObjectVector;
    typedef PtrObjectType*			ObjectIterator;

    typedef double*				DistanceVector;
    typedef double*				DistanceIterator;

    typedef Kratos::SearchUtils::SquaredDistanceFunction<Dim,PointType> DistanceFunction;

    typedef Kratos::BinsObjectStaticOCL< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, DistanceFunction > StaticBinsOCL;

    using namespace boost::python;

    class_< StaticBinsOCL, boost::noncopyable > ("BinsObjectStaticOCL", init<Kratos::ModelPart * , Kratos::ModelPart * ,OpenCL::DeviceGroup& > ())
    .def("GenerateBins", 	   		&StaticBinsOCL::GenerateBins)
    .def("LoadSample", 	   			&StaticBinsOCL::LoadSample)
    .def("AllocateOCLBuffers", 		&StaticBinsOCL::AllocateOCLBuffers)
    .def("InitializeBuffers",  		&StaticBinsOCL::InitializeBuffers)
    .def("SearchParticles",	   		&StaticBinsOCL::SearchParticles)
    .def("TransferStaticMeshToGPU",	&StaticBinsOCL::TransferStaticMeshToGPU)
    .def("TransferStaticMeshToCPU",	&StaticBinsOCL::TransferStaticMeshToCPU)
    .def("CopyStaticmeshData",		&StaticBinsOCL::CopyStaticmeshData)
    .def("TransferParticMeshToGPU", &StaticBinsOCL::TransferParticMeshToGPU)
    ;
}

}  // namespace Python.

} // Namespace Kratos

