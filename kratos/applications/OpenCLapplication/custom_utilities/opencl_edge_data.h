/*
==============================================================================
KratosOpenCLApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
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

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: mossaiby $
//   Date:                $Date: 2010-09-30 18:26:51 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_OPENCL_EDGE_DATA_H_INCLUDED)
#define KRATOS_OPENCL_EDGE_DATA_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
//#include "incompressible_fluid_application.h"

// Some useful macros, will be renamed if not consistent
#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	a.s0
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	a.s1
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	a.s2
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	a.s3
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	a.s4
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	a.s5
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	a.s6
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	a.s7
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	a.s8

#define KRATOS_OCL_MASS(a)				a.s9

#define KRATOS_OCL_NI_DNJ_0(a)			a.sa
#define KRATOS_OCL_NI_DNJ_1(a)			a.sb
#define KRATOS_OCL_NI_DNJ_2(a)			a.sc

#define KRATOS_OCL_DNI_NJ_0(a)			a.sd
#define KRATOS_OCL_DNI_NJ_1(a)			a.se
#define KRATOS_OCL_DNI_NJ_2(a)			a.sf

#define KRATOS_OCL_COMP(a, n)			a.s[n]

#define KRATOS_OCL_COMP_0(a)			a.x
#define KRATOS_OCL_COMP_1(a)			a.y
#define KRATOS_OCL_COMP_2(a)			a.z

namespace Kratos
{

	//
	// AllocateArray
	//
	// Helper function to allocate an array
	// Pass the address of the array variable

	template <typename _Type> void AllocateArray(_Type **_Array, unsigned int _Size1, unsigned int _Size2 = 1)
	{
		// Allocate memory
		*_Array = new _Type[_Size1 * _Size2];
	}

	//
	// FreeArray
	//
	// Helper function to free an array
	// Pass the address of the array variable

	template <typename _Type> void FreeArray(_Type **_Array)
	{
		// Free memory
		delete [] *_Array;
	}

	//
	// Array
	//
	// Helper function to access a 2D array using its 1D memory with given _Dim2

	template <typename _Type, unsigned int _Dim2> inline _Type &Array(_Type *_Array, unsigned int _i, unsigned int _j)
	{
		return *(_Array + _i * _Dim2 + _j);
	}

	//
	// Array3
	//
	// Helper function to access a 2D array using its 1D memory with _Dim2 = 3

	template <typename _Type> inline _Type &Array3(_Type *_Array, unsigned int _i, unsigned int _j)
	{
		return *(_Array + _i * 3 + _j);
	}

	//
	// Sqr
	//
	// Helper function to return square of a number

	template <typename _Type> inline _Type Sqr(_Type _X)
	{
		return _X * _X;
	}

	//
	// Norm2
	//
	// Helper function to calculate norm 2 of a vector

	template <typename _Type, unsigned int _Dim2> inline _Type Norm2(_Type *_Array, unsigned int _i)
	{
		_Type Sum();

		for (unsigned int Component = 0; Component < _Dim2; Component++)
		{
			Sum += Sqr(Array <_Type, _Dim2> (_Array, _i, Component));
		}

		return Sqrt(Sum);
	}

	//
	// Norm2_3
	//
	// Helper function to calculate norm 2 of a vector with _Dim2 = 3

	template <typename _Type> inline _Type Norm2_3(_Type *_Array, unsigned int _i)
	{
		return sqrt(Sqr(Array3(_Array, _i, 0)) + Sqr(Array3(_Array, _i, 1)) + Sqr(Array3(_Array, _i, 2)));
	}

/*
	void FillOpenCLEdge(cl_double16& a, const EdgesStructureType& b)
	{
		KRATOS_OCL_LAPLACIANIJ_0_0(a) = b.LaplacianIJ(0, 0);
		KRATOS_OCL_LAPLACIANIJ_0_1(a) = b.LaplacianIJ(0, 1);
		KRATOS_OCL_LAPLACIANIJ_0_2(a) = b.LaplacianIJ(0, 2);
		KRATOS_OCL_LAPLACIANIJ_1_0(a) = b.LaplacianIJ(1, 0);
		KRATOS_OCL_LAPLACIANIJ_1_1(a) = b.LaplacianIJ(1, 1);
		KRATOS_OCL_LAPLACIANIJ_1_2(a) = b.LaplacianIJ(1, 2);
		KRATOS_OCL_LAPLACIANIJ_2_0(a) = b.LaplacianIJ(2, 0);
		KRATOS_OCL_LAPLACIANIJ_2_1(a) = b.LaplacianIJ(2, 1);
		KRATOS_OCL_LAPLACIANIJ_2_2(a) = b.LaplacianIJ(2, 2);

		KRATOS_OCL_MASS(a) = b.Mass;

		KRATOS_OCL_NI_DNJ_0(a) = b.Ni_DNj[0];
		KRATOS_OCL_NI_DNJ_1(a) = b.Ni_DNj[1];
		KRATOS_OCL_NI_DNJ_2(a) = b.Ni_DNj[2];

		KRATOS_OCL_DNI_NJ_0(a) = b.DNi_Nj[0];
		KRATOS_OCL_DNI_NJ_1(a) = b.DNi_Nj[1];
		KRATOS_OCL_DNI_NJ_2(a) = b.DNi_Nj[2];
	}
*/
	//
	// OpenCLMatrixContainer
	//
	// A helper class for edge based mathods

	class OpenCLMatrixContainer
	{
		public:

			//
			// Used types

			typedef std::vector<cl_double16> EdgesVectorType;
			typedef unsigned int *IndicesVectorType;
			typedef double *CalcVectorType;
			typedef double *ValuesVectorType;

			//
			// OpenCLMatrixContainer
			//
			// Constructor

			OpenCLMatrixContainer(/*MatrixContainer _matrix_container*/)/*: matrix_container(_matrix_container)*/
			{
				// Nothing to do!
			}

			//
			// ~OpenCLMatrixContainer
			//
			// Destructor

			~OpenCLMatrixContainer()
			{
				// Nothing to do!
			}

			//
			// Functions to return private values

			// TODO: Remove unneeded ones

			inline unsigned int GetNumberEdges()
			{
				return mNumberEdges;
			}

			inline EdgesVectorType GetEdgeValues()
			{
				return mNonzeroEdgeValues;
			}

			inline IndicesVectorType GetColumnIndex()
			{
				return mColumnIndex;
			}

			inline IndicesVectorType GetRowStartIndex()
			{
				return mRowStartIndex;
			}

			inline ValuesVectorType GetLumpedMass()
			{
				return mLumpedMassMatrix;
			}

			inline ValuesVectorType GetInvertedMass()
			{
				return mInvertedMassMatrix;
			}

			inline CalcVectorType GetDiagGradient()
			{
				return mDiagGradientMatrix;
			}

			inline ValuesVectorType GetHmin()
			{
				return mHmin;
			}

			//
			// FillCoordinatesFromDatabase
			//
			// Function to access database

			void FillCoordinatesFromDatabase(CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Save value in the destination vector
					for (unsigned int component = 0; component < 3; component++)
					{
						Array3(rDestination, i_node, component) = (*node_it)[component];
					}
				}

				KRATOS_CATCH("");
			}

			//
			// FillVectorFromDatabase
			//
			// Function to access database

			void FillVectorFromDatabase(Variable <array_1d<double, 3> > &rVariable, CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested value in vector form
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable);

					// Save value in the destination vector
					for (unsigned int component = 0; component < 3; component++)
					{
						Array3(rDestination, i_node, component) = vector[component];
					}
				}

				KRATOS_CATCH("");
			}

			//
			// FillOldVectorFromDatabase
			//
			// Function to access database

			void FillOldVectorFromDatabase(Variable <array_1d<double, 3> > &rVariable, CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over alle nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested value in vector form
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable, 1);

					// Save value in the destination vector
					for(unsigned int component = 0; component < 3; component++)
					{
						Array3(rDestination, i_node, component) = vector[component];
					}
				}

				KRATOS_CATCH("");
			}

			//
			// FillScalarFromDatabase
			//
			// Function to access database

			void FillScalarFromDatabase(Variable <double> &rVariable, ValuesVectorType rDestination, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested scalar value
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable);

					// Save value in the destination vector
					rDestination[i_node] = scalar;
				}

				KRATOS_CATCH("");
			}

			//
			// FillOldScalarFromDatabase
			//
			// Function to access database

			void FillOldScalarFromDatabase(Variable <double> &rVariable, ValuesVectorType rDestination, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested scalar value
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable, 1);

					// Save value in the destination vector
					rDestination[i_node] = scalar;
				}

				KRATOS_CATCH("");
			}

			//
			// WriteVectorToDatabase
			//
			// Function to access database

			void WriteVectorToDatabase(Variable <array_1d<double, 3> > &rVariable, CalcVectorType rOrigin, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get reference of destination
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable);

					// Save vector in database
					for(unsigned int component = 0; component < 3; component++)
					{
						vector[component] = Array3(rOrigin, i_node, component);
					}
				}

				KRATOS_CATCH("");
			}

			//
			// WriteScalarToDatabase
			//
			// Function to access database

			void WriteScalarToDatabase(Variable <double> &rVariable, ValuesVectorType rOrigin, ModelPart::NodesContainerType &rNodes)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get reference of destination
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable);

					// Save scalar in database
					scalar = rOrigin[i_node];
				}

				KRATOS_CATCH("");
			}

		private:

			//MatrixContainer matrix_container;

			// Number of edges
			unsigned int mNumberEdges;

			// CSR data vector for storage of the G, L and consistent M components of edge ij
			EdgesVectorType mNonzeroEdgeValues;

			// Vector to store column indices of nonzero matrix elements for each row
			IndicesVectorType mColumnIndex;

			// Index vector to access the start of matrix row i in the column vector
			IndicesVectorType mRowStartIndex;

			// Inverse of the mass matrix ... for parallel calculation each subdomain should contain this correctly calculated (including contributions of the neighbours)
			ValuesVectorType mInvertedMassMatrix;

			// Minimum height around one node
			ValuesVectorType mHmin;

			// Lumped mass matrix (separately stored due to lack of diagonal elements of the consistent mass matrix)
			ValuesVectorType mLumpedMassMatrix;

			// Diagonal of the gradient matrix (separately stored due to special calculations)
			CalcVectorType mDiagGradientMatrix;

	};

} //namespace Kratos

#endif //KRATOS_OPENCL_EDGE_DATA_H_INCLUDED defined


