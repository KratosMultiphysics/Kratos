/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
/* *********************************************************   
*          
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2008-09-03 
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_COMPOSE_MATERIAL_H_INCLUDED)
#define KRATOS_COMPOSE_MATERIAL_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "includes/properties.h"
#include "geometries/geometry.h"
#include "includes/process_info.h"
#include "containers/data_value_container.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

	class ComposeMaterial : public ConstitutiveLaw<Node<3> >
	{
		public:


			  ///@name Type Definitions

			  typedef Matrix MatrixType;
		
                          typedef Vector VectorType;

			  typedef ConstitutiveLaw<Node<3> > BaseType;
			  
			  typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor;
			  
			  typedef boost::numeric::ublas::vector<Matrix> Third_Order_Tensor;
			  
			  typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;

			  typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; 
			  			  
			  typedef unsigned int SizeType;
			  
			  typedef unsigned int NumMaterial;
			  
			  typedef std::vector<ConstitutiveLaw<Node<3> >::Pointer> Materials;
                          

			  /**
			  * Counted pointer of ComposeMaterial
			  */
			  typedef boost::shared_ptr<ComposeMaterial> Pointer;

			  /**
			  * Life Cycle 
			  */
			  /**
			  * Default constructor.
			  */
                          ComposeMaterial();
			  ComposeMaterial(const Materials& Mat);


			  virtual boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const
			  {
			  boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new ComposeMaterial(mMaterials));
			  return p_clone;
			  }

			  /**
			  * Destructor.
			  */
			  virtual ~ComposeMaterial();



      bool Has( const Variable<double>& rThisVariable );

      bool Has( const Variable<Vector>& rThisVariable );

      bool Has( const Variable<Matrix>& rThisVariable );

      double GetValue( const Variable<double>& rThisVariable );

      Vector GetValue( const Variable<Vector>& rThisVariable );

      Matrix GetValue( const Variable<Matrix>& rThisVariable );

      void SetValue( const Variable<double>& rThisVariable, const double& rValue, 
      const ProcessInfo& rCurrentProcessInfo );

      void SetValue( const Variable<Vector>& rThisVariable, const Vector& rValue, 
      const ProcessInfo& rCurrentProcessInfo );

      void SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, 
      const ProcessInfo& rCurrentProcessInfo );

      void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
      const ProcessInfo& rCurrentProcessInfo);

      void InitializeMaterial( const Properties& props,
      const GeometryType& geom,
      const Vector& ShapeFunctionsValues );

      void InitializeSolutionStep( const Properties& props,
      const GeometryType& geom,
      const Vector& ShapeFunctionsValues ,
      const ProcessInfo& CurrentProcessInfo);

      void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& ConstitutiveMatrix);
     
      void  FinalizeSolutionStep( const Properties& props,
      const GeometryType& geom, 
      const Vector& ShapeFunctionsValues ,
      const ProcessInfo& CurrentProcessInfo);

      void CalculateStress( const Vector& StrainVector, Vector& StressVector);

      void CalculateCauchyStresses(
      Vector& rCauchy_StressVector,
      const Matrix& rF,
      const Vector& rPK2_StressVector,
      const Vector& rGreenLagrangeStrainVector);

      void CalculateStressAndTangentMatrix(Vector& StressVector,
      const Vector& StrainVector,
      Matrix& algorithmicTangent);


                  
		private:
		// atributos
		NumMaterial mNumMat; 
                Materials mMaterials;

      }; // Class Compose_Material
 }  // namespace Kratos.
#endif 
