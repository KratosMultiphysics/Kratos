
/* *********************************************************   
*          
*   Last Modified by:    $Author: nelson $
*   Date:                $Date: 2009-01-21 09:56:09 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes 

//custom_elements
// External includes 

// Project includes 

#include "includes/define.h"
#include "structural_application.h"
#include "custom_elements/beam_element.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"


namespace Kratos
{	
	//*****************************************************************************
	//*****************************************************************************
  
	BeamElement::BeamElement(IndexType NewId,GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	  {		
		      //DO NOT ADD DOFS HERE!!!
		      //THIS IS THE DEFAULT CONSTRUCTOR
	  }

	BeamElement::BeamElement(IndexType NewId,GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
	{
	KRATOS_TRY

	unsigned int dimension = GetGeometry().WorkingSpaceDimension();    // Dimension de trabajo: en 2D o 3D
	unsigned int Nodos = GetGeometry().size();                         // Cantidad de Nodos en el elemento


	if (dimension != 3.0)
	{
	std::cout<<"This element works only with a 2 node line and 3D dimension"<<std::endl;
	return;
	}
	for (unsigned int i=0; i < Nodos; i++)
	{
	GetGeometry()[i].pAddDof(DISPLACEMENT_X, REACTION_X);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	GetGeometry()[i].pAddDof(DISPLACEMENT_Y, REACTION_Y);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	GetGeometry()[i].pAddDof(DISPLACEMENT_Z, REACTION_Z);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	GetGeometry()[i].pAddDof(ROTATION_X,     MOMENTUM_X);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	GetGeometry()[i].pAddDof(ROTATION_Y,     MOMENTUM_Y);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	GetGeometry()[i].pAddDof(ROTATION_Z,     MOMENTUM_Z);      //	GRADOS DE LIBERTAD DEL ELEMENTO.
	}

	KRATOS_CATCH("")

	}


Element::Pointer BeamElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new BeamElement(NewId, GetGeometry().Create(ThisNodes), 
                                pProperties));
    }


	BeamElement::~BeamElement()
    {
    }

	//************************************************************************************
	//THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
	//************************************************************************************


void BeamElement::Initialize()
	{
	  KRATOS_TRY		
	  CalculateSectionProperties();
	  if(mCurrentDisplacement.size()!=12)	
	  {		
	  mCurrentDisplacement.resize(12,false);
	  mLoads.resize(12,false);
	  mGlobalMatrix.resize(12,12,false);
	  }
	  KRATOS_CATCH("")

	}

//************************************************************************************
//************************************************************************************
 
void BeamElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {

    }
//************************************************************************************
//************************************************************************************
 
void BeamElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {

    }


//************************************************************************************
//************************************************************************************


void BeamElement::CalculateAll(MatrixType& rLeftHandSideMatrix, 
    VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo,
    bool CalculateStiffnessMatrixFlag,bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
      //update the current displacement
      unsigned	int dimension = GetGeometry().WorkingSpaceDimension(); // Dimension de trabajo: en 2D o 3D
      unsigned int Nodos = GetGeometry().size();                      // Cantidad de Nodos en el elemento
      unsigned int Dof=Nodos*6;	

      if (dimension != 3)
      {
      std::cout<<"this element works only with a 2 node line and 3D dimension"<<std::endl;
      return;
      }


      Matrix LocalMatrix;
      Matrix Rotation;
      Matrix aux_matrix;  // Matriz auxiliar para relaizar el triple producto matricial. (R K Rt)
      //Matrix rLeftHandSideMatrix;

      LocalMatrix.resize(12,12, false);
      Rotation.resize(12,12, false);
      aux_matrix.resize(12,12, false);
      //rLeftHandSideMatrix.resize(12,12,false);


      CalculateLocalMatrix(LocalMatrix);
      CalculateTransformationMatrix(Rotation);
      noalias(aux_matrix) = prod(Rotation, LocalMatrix); 
      noalias(mGlobalMatrix)= prod(aux_matrix,trans(Rotation));
      CalculateLoads(Rotation, mLoads);

      mCurrentDisplacement(0)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);	
      mCurrentDisplacement(1)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
      mCurrentDisplacement(2)		=   GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
      mCurrentDisplacement(3)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_X);
      mCurrentDisplacement(4)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Y);
      mCurrentDisplacement(5)		=   GetGeometry()[0].GetSolutionStepValue(ROTATION_Z);
      mCurrentDisplacement(6)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);	
      mCurrentDisplacement(7)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
      mCurrentDisplacement(8)		=   GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);
      mCurrentDisplacement(9)		=   GetGeometry()[1].GetSolutionStepValue(ROTATION_X);
      mCurrentDisplacement(10)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Y);
      mCurrentDisplacement(11)	        =   GetGeometry()[1].GetSolutionStepValue(ROTATION_Z);


      //resizing and calculate the LHS contribution if needed
      if (CalculateStiffnessMatrixFlag == true) 
      {
      if(rLeftHandSideMatrix.size1() != Dof)
      rLeftHandSideMatrix.resize(Dof,Dof,false);
      noalias(rLeftHandSideMatrix) = ZeroMatrix(Dof,Dof); 
      CalculateLHS(rLeftHandSideMatrix);
      }
      //resizing and calculate the RHS contribution if needed
      if (CalculateResidualVectorFlag == true) 
      {
      if(rRightHandSideVector.size() != Dof)
      rRightHandSideVector.resize(Dof,false);
      noalias(rRightHandSideVector) = ZeroVector(Dof); 
      CalculateRHS(rRightHandSideVector);
      }

      //KRATOS_WATCH(rRightHandSideVector);
        KRATOS_CATCH("")
	}
	

//************************************************************************************
//************************************************************************************	

void  BeamElement::CalculateRightHandSide(VectorType& rRightHandSideVector, 
             ProcessInfo& rCurrentProcessInfo)
      {
	   //calculation flags
           bool CalculateStiffnessMatrixFlag = false;
           bool CalculateResidualVectorFlag = true;
           MatrixType temp = Matrix();
		
           CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag);
       }	

//************************************************************************************
//************************************************************************************	

 
void BeamElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
                VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
	   //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;
            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, 
            CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);
	
        }

//************************************************************************************
//************************************************************************************

void BeamElement::EquationIdVector(EquationIdVectorType& rResult, 
                ProcessInfo& CurrentProcessInfo)
 {
            if(rResult.size() != 12)
                    rResult.resize(12,false);

	    rResult[0]	= GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
	    rResult[1]	= GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[2]	= GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
	    rResult[3]	= GetGeometry()[0].GetDof(ROTATION_X).EquationId();
	    rResult[4]	= GetGeometry()[0].GetDof(ROTATION_Y).EquationId();
	    rResult[5]	= GetGeometry()[0].GetDof(ROTATION_Z).EquationId();
	    rResult[6]	= GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
	    rResult[7]	= GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
	    rResult[8]	= GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
	    rResult[9]	= GetGeometry()[1].GetDof(ROTATION_X).EquationId();
	    rResult[10] = GetGeometry()[1].GetDof(ROTATION_Y).EquationId();
	    rResult[11] = GetGeometry()[1].GetDof(ROTATION_Z).EquationId();
			
}


	//************************************************************************************
	//************************************************************************************
        
void BeamElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& 
                CurrentProcessInfo)
        {
	  ElementalDofList.resize(0);

	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_X));
	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Y));
	  ElementalDofList.push_back(GetGeometry()[0].pGetDof(ROTATION_Z));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_X));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Y));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_X));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_Y));
	  ElementalDofList.push_back(GetGeometry()[1].pGetDof(ROTATION_Z));
 
         }

	//************************************************************************************
	//************************************************************************************

void BeamElement::GetValuesVector(Vector& values, int Step)
        {
                if(values.size() != 12)
                    values.resize(12,false);

		  values(0)	= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X,Step);
		  values(1)	= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y,Step); 
		  values(2)	= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z,Step); 
		  values(3)	= GetGeometry()[0].GetSolutionStepValue(ROTATION_X,Step);
		  values(4)	= GetGeometry()[0].GetSolutionStepValue(ROTATION_Y,Step);
		  values(5)	= GetGeometry()[0].GetSolutionStepValue(ROTATION_Z,Step); 
		  values(6)	= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X,Step); 
		  values(7)	= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y,Step);
		  values(8)	= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z,Step);
		  values(9)	= GetGeometry()[1].GetSolutionStepValue(ROTATION_X,Step); 
		  values(10)	= GetGeometry()[1].GetSolutionStepValue(ROTATION_Y,Step); 
		  values(11)	= GetGeometry()[1].GetSolutionStepValue(ROTATION_Z,Step);
//	          values(12)	= GetGeometry()[0].GetSolutionStepValue(MOMENTUM_X,Step);
// 		  values(13)	= GetGeometry()[0].GetSolutionStepValue(MOMENTUM_Y,Step);
// 		  values(14)	= GetGeometry()[0].GetSolutionStepValue(FORCE_Z,Step);
// 		  values(15)	= GetGeometry()[1].GetSolutionStepValue(MOMENTUM_X,Step);
// 		  values(16)	= GetGeometry()[1].GetSolutionStepValue(MOMENTUM_Y,Step);
// 		  values(17)	= GetGeometry()[1].GetSolutionStepValue(FORCE_Z,Step);

              
          }               
				
//************************************************************************************
//************************************************************************************

void BeamElement::CalculateLHS(Matrix& rLeftHandSideMatrix)
    {
    noalias(rLeftHandSideMatrix)= mGlobalMatrix;
    return;
    }

   
//************************************************************************************
//************************************************************************************
void BeamElement::CalculateRHS(Vector& rRightHandSideVector)

{
    noalias(rRightHandSideVector) = mLoads;
    noalias(rRightHandSideVector) -= prod(mGlobalMatrix,mCurrentDisplacement); 
    return;
}




//************************************************************************************
//************************************************************************************

	void BeamElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
		{
		
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = false;
		Vector temp = Vector();
		CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
		
		}

//************************************************************************************
//************************************************************************************
		
double BeamElement::Unitarios(double a, double b) // Calculo de los vectores ortonormales.
		{		
				double normal;
				normal = a/b;
				return normal;
		}


//************************************************************************************
//************************************************************************************

void BeamElement::CalculateSectionProperties()

{
KRATOS_TRY

Vector x_zero(6);        // Vector que contiene coordenadas de los nodos.
Vector Vector_zero(3);   // Vector que contiene la direccion de la barra. 
double minimo,maximo,B;
double b,h;
//double const gravity = 9.80665;
b        = GetProperties()[BASE];                               
h        = GetProperties()[HEIGHT];
mPoisson = GetProperties()[POISSON_RATIO];                     
mYoungs  = GetProperties()[YOUNG_MODULUS];
 
mWeight        =  GetProperties()[BODY_FORCE](1);  
mElasticidad_Cortante	= mYoungs /(2*(1+mPoisson));                  
mInertia_x     = b*h*h*h/12.0;                                      
mInertia_y     = b*b*b*h/12.0;                                      
minimo         = std::min(b,h);
maximo	       = std::max(b,h);                                              
B	       = (1.00-0.63*(minimo/maximo)*(1-(pow(minimo,4)/(12*pow(maximo,4)))))/3;	// constante torsional. Solo para secciones rectangulares.
mInertia_Polar = B*minimo*minimo*minimo*maximo;											                            
mArea	       = b*h;														                                                 

x_zero(0)= GetGeometry()[0].X0();  x_zero(1)= GetGeometry()[0].Y0();  x_zero(2)= GetGeometry()[0].Z0();
x_zero(3)= GetGeometry()[1].X0();  x_zero(4)= GetGeometry()[1].Y0();  x_zero(5)= GetGeometry()[1].Z0();

for (unsigned int i=0; i<3; i++)
{
Vector_zero[i] = x_zero[i+3] - x_zero[i];
}

mReference_length = norm_2(Vector_zero); 


KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************

void BeamElement::CalculateLocalMatrix(Matrix& LocalMatrix)
{
 KRATOS_TRY

	int dimension = GetGeometry().WorkingSpaceDimension();     // Dimension de trabajo: en 2D o 3D
	//int Nodos = GetGeometry().size();                          // Cantidad de Nodos en el elemento
	//int Dof=Nodos*6;										                                 // Cantidad de grados de libertad en la barra.
	if (dimension != 3.0)



		if(LocalMatrix.size1()!=12 || LocalMatrix.size2()!=12)   // Matriz local de rigidez de la Estructura. 
		LocalMatrix.resize(12,12,false);
		noalias(LocalMatrix)   = zero_matrix<double>(12,12);

		//Inicializando matriz local de la estructura

		 double L	  =	         mReference_length;
		 double LL  =	         mReference_length* mReference_length;
		 double LLL	=	         mReference_length* mReference_length* mReference_length;

	         LocalMatrix(0,0)	=   (mArea * mYoungs)/(L);
		 LocalMatrix(6,0)	=   -(mArea * mYoungs)/(L);

		 LocalMatrix(1,1)	=   (12*mInertia_x*mYoungs)/(LLL);
		 LocalMatrix(5,1)	=   (6*mInertia_x* mYoungs)/(LL); 
 		LocalMatrix(7,1)	=   -(12*mInertia_x*mYoungs)/(LLL);
		 LocalMatrix(11,1)	=   (6*mInertia_x* mYoungs)/(LL); 

		 LocalMatrix(2,2)	=   (12*mInertia_y*mYoungs)/(LLL);
		 LocalMatrix(4,2)	=   -(6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(8,2)	=   -(12*mInertia_y*mYoungs)/(LLL);
		 LocalMatrix(10,2)	=   -(6*mInertia_y* mYoungs)/(LL);

		 LocalMatrix(3,3)	=   (mInertia_Polar* mElasticidad_Cortante)/L;
		 LocalMatrix(9,3)	=   -(mInertia_Polar* mElasticidad_Cortante)/L; 
	
		 LocalMatrix(2,4)	=   -(6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(4,4)	=   (4*mInertia_y*mYoungs)/L;
		 LocalMatrix(8,4)	=    (6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(10,4)	=    (2*mInertia_y*mYoungs)/L;

		 LocalMatrix(1,5)	=   (6*mInertia_x* mYoungs)/(LL); 
		 LocalMatrix(5,5)	=   (4*mInertia_x*mYoungs)/L;
		 LocalMatrix(7,5)	=   -(6*mInertia_x* mYoungs)/(LL); 
		 LocalMatrix(11,5)	=   (2*mInertia_x*mYoungs)/L;

		 LocalMatrix(0,6)	=    -(mArea * mYoungs)/( L);
		 LocalMatrix(6,6)	=   (mArea * mYoungs)/( L);
	
		 LocalMatrix(1,7)	=   -(12*mInertia_x*mYoungs)/(LLL);
		 LocalMatrix(5,7)	=   -(6*mInertia_x* mYoungs)/(LL); 
		 LocalMatrix(7,7)	=    (12*mInertia_x*mYoungs)/(LLL);
		 LocalMatrix(11,7)	=	-(6*mInertia_x* mYoungs)/(LL); 
		
		 LocalMatrix(2,8)	=	 -(12*mInertia_y*mYoungs)/(LLL);
		 LocalMatrix(4,8)	=    (6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(8,8)	=    (12*mInertia_y*mYoungs)/(LLL);
		 LocalMatrix(10,8)	=    (6*mInertia_y* mYoungs)/(LL);
		
		 LocalMatrix(3,9)	=   -(mInertia_Polar* mElasticidad_Cortante)/L;
		 LocalMatrix(9,9)	=    (mInertia_Polar* mElasticidad_Cortante)/L;
		
		 LocalMatrix(2,10)	=   -(6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(4,10)	=  (2*mInertia_y*mYoungs)/L;
		 LocalMatrix(8,10)	=  (6*mInertia_y* mYoungs)/(LL);
		 LocalMatrix(10,10)	=  (4*mInertia_y*mYoungs)/L;
		
		 LocalMatrix(1,11)	=   (6*mInertia_x* mYoungs)/(LL); 
		 LocalMatrix(5,11)	=   (2*mInertia_x*mYoungs)/L;
		 LocalMatrix(7,11)	=  -(6*mInertia_x* mYoungs)/(LL); 
		 LocalMatrix(11,11)	=  (4*mInertia_x*mYoungs)/L;

 KRATOS_CATCH("")

}


//*****************************************************************************
//*****************************************************************************

void BeamElement::CalculateTransformationMatrix(Matrix& Rotation)

{

	KRATOS_TRY

 Vector Normal_zero(9); // vector que contiene los cosenos directores.
	Vector x_zero(6);
 Vector Vector_zero(3);
 noalias(Normal_zero)	=	zero_vector<double>(9);
 noalias(x_zero)      =	zero_vector<double>(6);
 noalias(Vector_zero) =	zero_vector<double>(3);
 noalias(Rotation)    =	zero_matrix<double> (12,12);
	
 double nx, ny, nz,teta, phi;

	x_zero(0)= GetGeometry()[0].X0();  x_zero(1)= GetGeometry()[0].Y0();  x_zero(2)= GetGeometry()[0].Z0();
	x_zero(3)= GetGeometry()[1].X0();  x_zero(4)= GetGeometry()[1].Y0();  x_zero(5)= GetGeometry()[1].Z0();
		
	for (unsigned int i=0; i<3; i++)
	  {
	    Vector_zero[i] = x_zero[i+3] - x_zero[i];
		 }

	for (unsigned int i=0; i<3; i++)
	{
		Normal_zero[i] = Unitarios(Vector_zero[i], mReference_length);
	}

	//double a = 1.0;  // Variable axiliar para calcular valor de pi. 
	//pi= 4*atan(a);
	nx = Normal_zero[0];
 ny = Normal_zero[1];
 nz = Normal_zero[2];
 //KRATOS_WATCH(nx)
 //KRATOS_WATCH(ny)
 //KRATOS_WATCH(nz)	

	if (nx ==0.0)
		{
		 teta = PI/2;
			 if (ny == 0.0)
			 {
			 teta = 0.0;
			 phi  = PI/2;
			 }
			else
			phi = atan(nz/sqrt(nx*nx+ny*ny));
	}
    else
	{
        teta = atan(ny/nx);
        phi  = atan(nz/sqrt(nx*nx+ny*ny));
	}
      
      if(nx < 0.0) 
		teta = teta + PI;

	    Normal_zero[3] = -sin(teta);
	    Normal_zero[4] =  cos(teta);
	    Normal_zero[5]=  0.0;
	    Normal_zero[6]= -nz*cos(teta);
	    Normal_zero[7]= -nz*sin(teta);
	    Normal_zero[8]=  nx*cos(teta) + ny*sin(teta);


 // Creacion de la matriz de transformacion.
	for (unsigned int kk=0; kk < 12; kk += 3)
	{
		for (unsigned int i=0; i<3; i++)
		{
			for(unsigned int j=0; j<3; j++)
			{
				Rotation(i+kk,j+kk)=Normal_zero(3*j+i);
			}
		}
	}

KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************
 void BeamElement::CalculateLoads(Matrix Rotation, Vector& mLoads)

{
		KRATOS_TRY
  //Creacion de los vectores de cargas externas.
  // Fuerzas externas Uniformente distriduida. Por lo general es un dato suministrado por el usuario.
  // Cambiaro de posicion una vez terminado el programa.
		
			double alpha=0.00;
			double signo=1.00;
			Vector Cargas;
			Vector Load;
			Vector Normal_Loads;
   Vector x_zero;
   Vector Vector_zero;
		 Cargas.resize(12,false);
   Load.resize(2,false);
   Normal_Loads.resize(3,false);
   x_zero.resize(6,false);
   Vector_zero.resize(3,false);


	  x_zero(0)= GetGeometry()[0].X0();  x_zero(1)= GetGeometry()[0].Y0();  x_zero(2)= GetGeometry()[0].Z0();
		 x_zero(3)= GetGeometry()[1].X0();  x_zero(4)= GetGeometry()[1].Y0();  x_zero(5)= GetGeometry()[1].Z0();
		
		  for (unsigned int i=0; i<3; i++)
	      	{
	        	Vector_zero[i] = x_zero[i+3] - x_zero[i];
		      }

			Normal_Loads[0]	= Vector_zero[0] ;
			Normal_Loads[1]	= 0.00 ;
			Normal_Loads[2]	= Vector_zero[2];

			if (Vector_zero[1]<0) {signo =-1.00;}
			if( norm_2(Normal_Loads)==0 || norm_2( Vector_zero)==0  )
			{
				alpha = signo*PI/2;
			}
			else
			{
				alpha = inner_prod(Normal_Loads,Vector_zero)/(norm_2(Vector_zero)*norm_2( Normal_Loads));
				alpha	= signo*acos(alpha);
			}

			// las fuerzas consideradas son las de peso propio.
			Load[0]= mArea*mWeight*sin(alpha);         // Carga Axialmente Distribuida.
 		Load[1]= mArea*mWeight*cos(alpha);         // Carga en la Direccion gravedad
 			
	/*	if(Loads.size()!=12)   // Matriz local de rigidez de la Estructura. 
		Loads.resize(12,false);
		noalias(Loads)	= zero_vector<double> (12); 
	*/
		
		Cargas[0]=   -Load[0]*mReference_length/2.00;							                        // Fuerza en X;
		Cargas[1]=   -(Load[1]*mReference_length)/2.00;						                       // Fuerza en Y; graveded
		Cargas[2]=   0.00;														                                             // Fuerza en Z
		Cargas[3]=   0.00;														                                             // Momento Tersor X;
		Cargas[4]=   0.00;														                                             // Momento Y
		Cargas[5]=  -(Load[1])*mReference_length*mReference_length/12.00;;										// Momento Z
		Cargas[6]=   -Load[0]*mReference_length/2.00;
		Cargas[7]=   -(Load[1])*mReference_length/2.00;
		Cargas[8]=	   0.00;
		Cargas[9]=    0.00;
		Cargas[10]=   0.00;
		Cargas[11]=   (Load[1])*mReference_length*mReference_length/12.00;
		 
	 noalias(mLoads)= prod(Rotation,Cargas);		// Cargas externas en coordenadas globales. 
  //KRATOS_WATCH(mLoads)
  KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
 
 void BeamElement::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
            {
		
	      KRATOS_TRY
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		rMassMatrix = ZeroMatrix(MatSize,MatSize);

		double TotalMass = mArea*mReference_length*GetProperties()[DENSITY];
		
		Vector LumpFact;
		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		for(unsigned int i=0; i<NumberOfNodes; i++)
		{
			double temp = LumpFact[i]*TotalMass;
			for(unsigned int j=0; j<dimension; j++)
			{
				unsigned int index = i*dimension + j;
                                rMassMatrix(index,index) = temp;
				if (index==3 or index==4 or index==5)
				      rMassMatrix(index,index) = 0.00;
			        
			}
		}
		KRATOS_CATCH("")
	    }



//************************************************************************************
	//************************************************************************************
	  void BeamElement::GetFirstDerivativesVector(Vector& values, int Step)
	{
		 KRATOS_TRY
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = 2*number_of_nodes*dim;
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index =  i*number_of_nodes*dim;


			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
			values[index + 3] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_X,Step);
			values[index + 4] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_Y,Step);
			values[index + 5] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_VELOCITY_Z,Step);
		}

	 
          KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	  void BeamElement::GetSecondDerivativesVector(Vector& values, int Step)
	{
		 KRATOS_TRY
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = 2*number_of_nodes*dim;
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*number_of_nodes*dim;
	                //KRATOS_WATCH(index)
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
		        values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
			values[index + 3] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_X,Step);
			values[index + 4] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Y,Step);
			values[index + 5] = 0.00; // GetGeometry()[i].GetSolutionStepValue(ANGULAR_ACCELERATION_Z,Step);

		}
              
               KRATOS_CATCH("")

	}
	//************************************************************************************
	//************************************************************************************


} // Namespace Kratos


