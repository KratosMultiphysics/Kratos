
//   Project Name:        Kratos       
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-10 11:10:16 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_BEAM_ELEMENT_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "structural_application.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

    class BeamElement
    :public Element
	{
               
	      typedef GeometryData::IntegrationMethod IntegrationMethod;

      
     private:
     ///@name Static Member Variables  

				

	      double mArea;                            // Area de la seccion tranversal de la viga. 
	      double mInertia_x;                       // Momento de Inercia alredor del eje Ix local.
	      double mInertia_y;                       // Momento de Inercia alrededor del eje Iy local.
	      double mInertia_Polar;                   // Momento Polar de Inercia
	      double mlength;                          // Longitud del Elemento.  


	      void CalculateSectionProperties();

	      void CalculateLocalMatrix(Matrix& LocalMatrix);

	      void CalculateTransformationMatrix(Matrix& Rotation);

	      void CalculateBodyForce(Matrix& Rotation,  Vector& LocalBody, Vector& GlobalBody);
	      
	      void CalculateLocalNodalStress(Vector& Stress); 
	      
	      double CalculateInternalForces(const double& Mo, const double& Po, const double& Load, const double& X);
	      
	      void CalculateDistrubuitedBodyForce(const int Direction, Vector& Load);
	      

	      public:
	      
	      

	      KRATOS_CLASS_POINTER_DEFINITION(BeamElement); 

	      /// Default constructor.
	      BeamElement(IndexType NewId, GeometryType::Pointer pGeometry);
	      BeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

	      /// Destructor.
	      virtual ~BeamElement();


	      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

	      void Initialize();

	      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

	      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

	      void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

	      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

	      void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

	      void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);

	      void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo);

	      void GetValuesVector(Vector& values, int Step);

	      void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

	      void CalculateRHS(Vector& rRightHandSideVector);

	      void CalculateLHS(Matrix& rLeftHandSideMatrix);

	      void CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
	      ProcessInfo& rCurrentProcessInfo,
	      bool CalculateStiffnessMatrixFlag,
	      bool CalculateResidualVectorFlag);

	      void  GetFirstDerivativesVector(Vector& values, int Step);
	      void  GetSecondDerivativesVector(Vector& values, int Step);
	      				     
              void CalculateOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                            std::vector< array_1d<double,3> >& Output, 
                                            const ProcessInfo& rCurrentProcessInfo);
	      
	      void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                           std::vector<array_1d<double,3> >& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo);

             IntegrationMethod GetIntegrationMethod();
             

		
	      
	      

}; // Class BeamElement

} // Namespace Kratos.
			   
				
#endif 

