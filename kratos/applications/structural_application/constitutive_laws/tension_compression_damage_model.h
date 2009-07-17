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

#if !defined(KRATOS_TENSION_COMPRESSION_DAMAGE_MODEL_H_INCLUDED)
#define  KRATOS_TENSION_COMPRESSION_DAMAGE_MODEL_H_INCLUDED

// System includes 

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

	class Tension_Compression_Damage_Model: public ConstitutiveLaw<Node<3> >
	{
		public:


			  ///@name Type Definitions
			  typedef ConstitutiveLaw<Node<3> > BaseType;
			  
			  typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector
			  
			  typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
			  typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
			  typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.

                          //typedef FluencyCriteria<Node<3> >::Pointer FluencyCriteriaType; 
		      
		      
			  /**
			  * Counted pointer of Tension_Compression_Damage_Model
			  */
			  typedef boost::shared_ptr<Tension_Compression_Damage_Model> Pointer;

			  /**
			  * Life Cycle 
			  */
			  /**
			  * Default constructor.
			  */
			  Tension_Compression_Damage_Model();

			  virtual boost::shared_ptr<ConstitutiveLaw<Node<3> > > Clone() const
			  {
			  boost::shared_ptr<ConstitutiveLaw<Node<3> > > p_clone(new Tension_Compression_Damage_Model());
			  return p_clone;
			  }

			  /**
			  * Destructor.
			  */
			  virtual ~Tension_Compression_Damage_Model();

			  /**
			  * Operators 
			  */
			  /**
			  * Operations
			  */
			  bool Has( const Variable<double>& rThisVariable );
			  bool Has( const Variable<Vector>& rThisVariable );
			  bool Has( const Variable<Matrix>& rThisVariable );

			  double GetValue( const Variable<double>& rThisVariable );
			  Vector GetValue( const Variable<Vector>& rThisVariable );
			  Matrix GetValue( const Variable<Matrix>& rThisVariable );

			  void SetValue( const Variable<double>& rVariable, 
			  const double& Value, 
			  const ProcessInfo& rCurrentProcessInfo );
			  void SetValue( const Variable<Vector>& rThisVariable, 
			  const Vector& rValue, 
			  const ProcessInfo& rCurrentProcessInfo );
			  void SetValue( const Variable<Matrix>& rThisVariable, 
			  const Matrix& rValue, 
			  const ProcessInfo& rCurrentProcessInfo );


			  void InitializeMaterial( const Properties& props,
			  const GeometryType& geom,
			  const Vector& ShapeFunctionsValues );

			  void CalculateConstitutiveMatrix(const Vector& StrainVector, Matrix& rResult);

			  void CalculateStress( const Vector& StrainVector, 
			  Vector& StressVector);

			  void CalculateCauchyStresses( Vector& Cauchy_StressVector,
			  const Matrix& F,
			  const Vector& PK2_StressVector,
			  const Vector& GreenLagrangeStrainVector);


			  void Calculate(const Variable<Matrix >& rVariable, Matrix& rResult, 
			  const ProcessInfo& rCurrentProcessInfo);

			  void InitializeSolutionStep( const Properties& props,
			  const GeometryType& geom,
			  const Vector& ShapeFunctionsValues ,
			  const ProcessInfo& CurrentProcessInfo);

			  void FinalizeSolutionStep( const Properties& props,
			  const GeometryType& geom, 
			  const Vector& ShapeFunctionsValues ,
			  const ProcessInfo& CurrentProcessInfo);

			  void CalculateStressAndTangentMatrix( Vector& StressVector,
			  const Vector& StrainVector,
			  Matrix& algorithmicTangent);


		private:
			  double md_pos;
			  double md_neg;
			  double mr_pos_old;
			  double mr_pos_new;
			  double mr_neg_old;
			  double mr_neg_new;
			  double mEc;
			  double mEt; 
			  double mFc;
			  double mFt;
			  double mGE;
			  double mNU;
                          double ml;	
                          double malfa;	

			  Vector msplit_negative_stress;
			  Vector msplit_positive_stress;
                          Vector mstressvector;
			  Vector mCurrentStress; 
                          Matrix mCtangent; 
                          Vector mInSituStress; 
			  Matrix_Second_Tensor mp;
			  Fourth_Order_Tensor mP;

                          //FluencyCriteriaType mFluencyCriteria;

						
			  // Miembros Privados
			  void CalculateNoDamageElasticMatrix(Matrix& C, const double E, const double NU);
			  void CalculateNoDamageStress(const Vector& StrainVector, Vector& rResult);
			  void CalculateDamage(const Vector& Stress);
                          void Compute_Principal_Stress(const Vector& StressVector, Vector& Result);
                          void Verify_Matrix(Matrix& Result);
			  void Compute_pij_Tensor(const Second_Order_Tensor& A, const Second_Order_Tensor& B);
                          void Compute_P_Fourth_Order_Tensor(Vector& EigenValues);






      }; // Class Tension_Compression_Damage_Model 
 }  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_2D_H_INCLUDED  defined 
