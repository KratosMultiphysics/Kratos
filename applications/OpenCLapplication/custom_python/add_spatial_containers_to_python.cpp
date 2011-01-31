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

#include "custom_utilities/opencl_interface.h"
#include "spatial_containers/bins_static_OCL.h"

template< std::size_t dim_type>
class Point {

  public:

      double       coord[dim_type];
      std::size_t  id;
      std::size_t  tag;
      //int id;

      double& operator[](std::size_t i) {return coord[i];}

      double const & operator[](std::size_t i) const {return coord[i];}
};

template< class T, std::size_t dim >
class PointDistance{
   public:
      double operator()( T const& p1, T const& p2 ){
         double dist = 0.0;
         for( std::size_t i = 0 ; i < dim ; i++){
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
         }
         return sqrt(dist);
      }
};

template< class T, std::size_t dim >
class PointDistance2{
   public:
      double operator()( T const& p1, T const& p2 ){
         double dist = 0.0;
         for( std::size_t i = 0 ; i < dim ; i++){
            double tmp = p1[i] - p2[i];
            dist += tmp*tmp;
         }
         return dist;
      }
};

template< std::size_t dim_type >
std::ostream & operator<<( std::ostream& rOut, Point<dim_type> & rPoint){
	 rOut << "(" << rPoint.id << ") ";
   for(std::size_t i = 0 ; i < dim_type ; i++)
      rOut << rPoint[i] << " "; 
   return rOut; 
};

template< std::size_t dim_type >
std::istream & operator>>( std::istream& rIn, Point<dim_type> & rPoint){
   for(std::size_t i = 0 ; i < dim_type ; i++)
      rIn >> rPoint[i]; 
   return rIn; 
};

namespace Kratos
{
    
namespace Python
{
    void  AddSpatialContainersToPython()
    {

        static const std::size_t Dim = 3;
	
        typedef Point<Dim>              PointType;

	typedef PointType*        	PtrPointType;
	typedef PtrPointType*     	PointVector;
	typedef PtrPointType*           PointIterator;
	
	typedef double*                 DistanceVector;
        typedef double*                 DistanceIterator;
	
	typedef Kratos::BinsOCL< Dim, PointType, PointVector, PtrPointType, PointIterator, DistanceIterator, PointDistance2<PointType,Dim> > StaticBinsOCL;
      
	using namespace boost::python;
	
	class_< OpenCL::DeviceGroup, boost::noncopyable > ("OpenCLDeviceGroup", init<cl_device_type, bool> ())
		.def("AddCLSearchPath", &OpenCL::DeviceGroup::AddCLSearchPath)
		;
		
	class_< StaticBinsOCL, boost::noncopyable > ("BinsOCL", init<PointVector,PointVector,cl_int4 *,int,int,OpenCL::DeviceGroup&,size_t> ())
		.def("GenerateBins", 	   &StaticBinsOCL::generateBins)
		.def("allocateOCLBuffers", &StaticBinsOCL::allocateOCLBuffers)
		.def("searchInRadiusOCL",  &StaticBinsOCL::searchInRadiusOCL)
		.def("searchNearestOCL",   &StaticBinsOCL::searchNearestOCL)
		;
		
	enum_<cl_device_type>("cl_device_type")
                          .value("CL_DEVICE_TYPE_CPU", CL_DEVICE_TYPE_CPU)
                          .value("CL_DEVICE_TYPE_GPU", CL_DEVICE_TYPE_GPU)
                          .value("CL_DEVICE_TYPE_ALL", CL_DEVICE_TYPE_ALL)
                        ;
  }
	
}  // namespace Python.

} // Namespace Kratos

