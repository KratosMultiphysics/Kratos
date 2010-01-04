
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
               
      
     private:
     ///@name Static Member Variables  

				
				Matrix mGlobalMatrix;					             
				//Vector mUniformLoads;					               // Valores de carga uniformenete distribuida. Axial, X or Y, Gravedad, Momento Distribuido. 
				Vector mLoads;						                 // Vector cuyo contenido  son las cargas equivalentes de load: Fuerzas y Momentos. 
				Vector mCurrentDisplacement;			          // Este Matriz incluye la rotacion y los desplazamientos preescritos

				double mArea;                            // Area de la seccion tranversal de la viga. 
				double mInertia_x;                       // Momento de Inercia alredor del eje Ix local.
				double mInertia_y;                       // Momento de Inercia alrededor del eje Iy local.
				double mInertia_Polar;                   // Momento Polar de Inercia
				double mPoisson;                         // razon de Poisson
				double mYoungs;                          // Modulo de Young
				double mWeight;                          // Peso Especifico  
 			double mElasticidad_Cortante;            // Modulo de elasticidad al Cortante.
				double mReference_length;                // Longitud del Elemento.  
	

				double Unitarios(double a, double b);

				void CalculateSectionProperties();

				void CalculateLocalMatrix(Matrix& LocalMatrix);

				void CalculateTransformationMatrix(Matrix& Rotation);

    void CalculateLoads(Matrix Rotation, Vector& mLoads); 


				public:
      
				typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef ConstitutiveLaw<Node<3> > ConstitutiveLawType;

    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
     

    KRATOS_CLASS_POINTER_DEFINITION(BeamElement); 
      
      /// Default constructor.
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

      /// Destructor.
    virtual ~BeamElement();
      
			   // IntegrationMethod GetIntegrationMethod();

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


		

}; // Class BeamElement

} // Namespace Kratos.
			   
				
#endif 

