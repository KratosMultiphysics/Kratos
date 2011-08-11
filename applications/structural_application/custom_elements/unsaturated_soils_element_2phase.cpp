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
*   Last Modified by:    $Author: janosch $
*   Date:                $Date: 2009-01-14 09:24:26 $
*   Revision:            $Revision: 1.9 $
*
* ***********************************************************/

// System includes 


// External includes 

// Project includes 
// #include "includes/define.h"
#include "structural_application.h"
#include "custom_elements/unsaturated_soils_element_2phase.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"

namespace Kratos
{

    UnsaturatedSoilsElement_2phase::UnsaturatedSoilsElement_2phase(IndexType NewId, 
            GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {		
		//DO NOT ADD DOFS HERE!!!
    }

	//************************************************************************************
	//************************************************************************************
    UnsaturatedSoilsElement_2phase::UnsaturatedSoilsElement_2phase(IndexType NewId, 
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
		//setting up the nodal degrees of freedom
        //with DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z
        //DOFs at the end of time step
        //All calculations are made on the general midpoint alpha
        // Variables DOF_ALPHA are updated in the scheme
        if(GetGeometry().size()== 27 ||GetGeometry().size()== 20 || GetGeometry().size()== 10  || GetGeometry().size()== 15)
        {
			if(GetGeometry().size()== 27 )
                {
                        mNodesPressMin=1;
                        mNodesPressMax=8;
                        mNodesDispMin=1;
                        mNodesDispMax=27;
                        mpPressureGeometry = Geometry< Node<3> >::Pointer(new Hexahedra3D8 <Node<3> >(
                                        GetGeometry()(0),GetGeometry()(1),GetGeometry()(2),GetGeometry()(3),
                                        GetGeometry()(4),GetGeometry()(5),GetGeometry()(6),GetGeometry()(7)));
                        mThisIntegrationMethod= GeometryData::GI_GAUSS_3;
                }
                if(GetGeometry().size()== 20 )
                {
                        mNodesPressMin=1;
                        mNodesPressMax=8;
                        mNodesDispMin=1;
                        mNodesDispMax=20;
                        mpPressureGeometry = Geometry< Node<3> >::Pointer(new Hexahedra3D8 <Node<3> >(
                                        GetGeometry()(0),GetGeometry()(1),GetGeometry()(2),GetGeometry()(3),
                                        GetGeometry()(4),GetGeometry()(5),GetGeometry()(6),GetGeometry()(7)));
                        mThisIntegrationMethod= GeometryData::GI_GAUSS_3;
                }
                if(GetGeometry().size()== 10 )
                {
                        mNodesPressMin=1;
                        mNodesPressMax=4;
                        mNodesDispMin=1;
                        mNodesDispMax=10;
                        mpPressureGeometry = Geometry< Node<3> >::Pointer(new Tetrahedra3D4 <Node<3> >(
                                        GetGeometry()(0),GetGeometry()(1),GetGeometry()(2),GetGeometry()(3)));
                        mThisIntegrationMethod= GeometryData::GI_GAUSS_5;
                }

                if(GetGeometry().size()== 15 )
                {
                        mNodesPressMin=1;
                        mNodesPressMax=6;
                        mNodesDispMin=1;
                        mNodesDispMax=15;
                        mpPressureGeometry = Geometry< Node<3> >::Pointer(new Prism3D6 <Node<3> >(
                                        GetGeometry()(0),GetGeometry()(1),GetGeometry()(2),GetGeometry()(3),GetGeometry()(4),GetGeometry()(5)));
                        mThisIntegrationMethod= GeometryData::GI_GAUSS_3;
                }
        }     
        else
            KRATOS_ERROR(std::logic_error, "This element matches only with a quadratic hexaeder (20 or 27), tetraeder (10) or prism (15) geometry" , *this);

    }

    Element::Pointer UnsaturatedSoilsElement_2phase::Create(IndexType NewId, 
            NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new UnsaturatedSoilsElement_2phase(NewId, GetGeometry().Create(ThisNodes), 
                                pProperties));
    }

    UnsaturatedSoilsElement_2phase::~UnsaturatedSoilsElement_2phase()
    {
    }
	//************************************************************************************
	//************************************************************************************
    void UnsaturatedSoilsElement_2phase::Initialize()
    {
       KRATOS_TRY

		unsigned int dim = GetGeometry().WorkingSpaceDimension();

       	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		mInvJ0.resize(integration_points.size());

     	for(unsigned int i=0; i< integration_points.size(); i++)
      	{
  			mInvJ0[i].resize(dim,dim,false);
    		noalias(mInvJ0[i])= ZeroMatrix(dim,dim);
    	}

		mDetJ0.resize(integration_points.size());
        noalias(mDetJ0)= ZeroVector(integration_points.size());

        GeometryType::JacobiansType J0(integration_points.size());

        J0 = GetGeometry().Jacobian(J0, mThisIntegrationMethod);

		//calculating the inverse J0
        for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
        {
			//calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);
        }

		//Constitutive Law initialisation
        if(mConstitutiveLawVector.size() != integration_points.size())
        {
            mConstitutiveLawVector.resize(integration_points.size());
            InitializeMaterial();
        }

       for(unsigned int i = (mNodesDispMin-1) ; i < mNodesDispMax ; i++)
        {  
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_NULL)=
			(GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT);
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_EINS)=
			(GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT);
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_DT)=ZeroVector(3);
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_NULL_DT)=ZeroVector(3);
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_EINS_DT)=ZeroVector(3);
   	    	(GetGeometry()[i]).GetSolutionStepValue(ACCELERATION)=ZeroVector(3);
   	    	(GetGeometry()[i]).GetSolutionStepValue(ACCELERATION_NULL)=ZeroVector(3);
   	    	(GetGeometry()[i]).GetSolutionStepValue(ACCELERATION_EINS)=ZeroVector(3);
            (GetGeometry()[i]).GetSolutionStepValue(DISPLACEMENT_OLD)=ZeroVector(3);
        }
        for(unsigned int i = (mNodesPressMin-1) ; i < mNodesPressMax ; i++)
        { 
            (GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_NULL)=
			(GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE);
            (GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_EINS)=
			(GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE);
            (GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_DT)=0;
            (GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_NULL_DT)=0;
            (GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_EINS_DT)=0;
 	    	(GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_ACCELERATION)=0;
 	    	(GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)=0;
 	    	(GetGeometry()[i]).GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)=0;
//             (GetGeometry()[i]).GetSolutionStepValue(DARCY_VELO_WATER)=ZeroVector(3);
	       (GetGeometry()[i]).GetSolutionStepValue(REACTION_WATER_PRESSURE)=0;
	       (GetGeometry()[i]).GetSolutionStepValue(REACTION_AIR_PRESSURE)=0;
	           (GetGeometry()[i]).GetSolutionStepValue(PRESSURE)=0;
        }

        KRATOS_CATCH("")
    }

	//************************************************************************************
	//************************************************************************************
    void UnsaturatedSoilsElement_2phase::CalculateAll(MatrixType& rLeftHandSideMatrix, 
            VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag,bool CalculateResidualVectorFlag)
    {
         KRATOS_TRY
         unsigned int number_of_nodes_disp = (mNodesDispMax-mNodesDispMin+1);
         unsigned int number_of_nodes_press = (mNodesPressMax-mNodesPressMin+1);
//          unsigned int number_of_nodes = number_of_nodes_disp+number_of_nodes_press;
         unsigned int dim = GetGeometry().WorkingSpaceDimension();

         //resizing as needed the LHS
         unsigned int MatSize1=(number_of_nodes_disp*dim+number_of_nodes_press);
         unsigned int MatSizeU= number_of_nodes_disp*dim;
         unsigned int MatSizeP= number_of_nodes_press;      

         if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
         {
                 if(rLeftHandSideMatrix.size1() != MatSize1)
                 rLeftHandSideMatrix.resize(MatSize1,MatSize1,false);
                 noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize1,MatSize1); //resetting LHS
         }

         //resizing as needed the RHS
         if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
         {
                 if(rRightHandSideVector.size() != MatSize1)
                       rRightHandSideVector.resize(MatSize1,false);
                 noalias(rRightHandSideVector) = ZeroVector(MatSize1); //resetting RHS    
         }

         //reading integration points and local gradients
         const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =   
                GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
        const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =   
                mpPressureGeometry->ShapeFunctionsLocalGradients(mThisIntegrationMethod);

        const Matrix& Ncontainer_Displacement = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

        const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues(mThisIntegrationMethod);

        Matrix Help_K_UU(MatSizeU,MatSizeU);
        Matrix Help_K_UW(MatSizeU,MatSizeP);
                 
        Matrix Help_K_WU(MatSizeP,MatSizeU);
        Matrix Help_K_WW(MatSizeP,MatSizeP);
        
        Vector Help_R_U(MatSizeU);
        Vector Help_R_W(MatSizeP);
 
//         std::vector<std::vector<Matrix> > tanC_U(dim,dim);
        array_1d<double,81> tanC_U;
        Matrix tanC_W(dim,dim);
                
        Matrix DN_DX_DISP(number_of_nodes_disp, dim);
                
        Matrix DN_DX_PRESS(number_of_nodes_press, dim);
                
        double Weight;

        Vector N_DISP(number_of_nodes_disp);
                
        Vector N_PRESS(number_of_nodes_press);
                
        Matrix StrainTensor;
                
        Matrix StressTensor;
    
        double capillaryPressure;
                        
        double waterPressure;
                        
        double porosity;
                
        double density;

        Matrix du_dx(dim,dim); 

        double DetJ=0.0;

               if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
               {
                        noalias(Help_K_UU)= ZeroMatrix(MatSizeU,MatSizeU);
                        noalias(Help_K_UW)=ZeroMatrix(MatSizeU,MatSizeP);
                        
                        noalias(Help_K_WU)= ZeroMatrix(MatSizeP,MatSizeU);
                        noalias(Help_K_WW)=ZeroMatrix(MatSizeP,MatSizeP);

                  }
                  if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
                  {
                         noalias(Help_R_U)= ZeroVector(MatSizeU);
                         noalias(Help_R_W)= ZeroVector(MatSizeP);
                   }

                  /////////////////////////////////////////////////////////////////////////
                  //// Integration in space sum_(beta=0)^(number of quadrature points)
                  /////////////////////////////////////////////////////////////////////////
                  for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
                  {  
               			if(DN_DX_PRESS.size1()!= (number_of_nodes_press) || DN_DX_PRESS.size2()!= dim)
                     		DN_DX_PRESS.resize(number_of_nodes_press,dim);
						noalias(DN_DX_PRESS)= ZeroMatrix(number_of_nodes_press, dim);

						noalias(DN_DX_PRESS)= prod(DN_De_Pressure[PointNumber],mInvJ0[PointNumber]);

               			if(DN_DX_DISP.size1()!= (number_of_nodes_disp) || DN_DX_DISP.size2()!= dim)
                     		DN_DX_DISP.resize(number_of_nodes_disp,dim);
						noalias(DN_DX_DISP)= ZeroMatrix(number_of_nodes_disp, dim);

						noalias(DN_DX_DISP)= prod(DN_De_Displacement[PointNumber],mInvJ0[PointNumber]);

						double DetDef= Determinant_DeformationTensor(DN_DX_DISP);
                        //Initializing 4dim tangential constitutive Matrix
//                           for(unsigned int i=0; i<dim; i++)
//                           {
//                                    for(unsigned int j=0; j<dim; j++)
//                                     {
//                                            if(tanC_U[i][j].size1()!= dim || 
//                                               tanC_U[i][j].size2()!= dim)
//                                                 tanC_U[i][j].resize(dim,dim);
//                                            noalias(tanC_U[i][j]) = ZeroMatrix(dim,dim);
//                                      }
//                            }

                           for( int i=0; i<81; i++ )
                               tanC_U[i]=0.0;
                           if(tanC_W.size1()!= dim || tanC_W.size2()!= dim)
                                    tanC_W.resize(dim,dim,false);
                           noalias(tanC_W)= ZeroMatrix(dim,dim); 

                           Weight= integration_points[PointNumber].Weight();

                           DetJ= mDetJ0[PointNumber];

                        	// Shape Functions on current spatial quadrature point
                           if(N_PRESS.size()!=number_of_nodes_press)
                                      N_PRESS.resize(number_of_nodes_press);
                           noalias(N_PRESS)= row(Ncontainer_Pressure,PointNumber);
                           if(N_DISP.size()!=number_of_nodes_disp)
                                      N_DISP.resize(number_of_nodes_disp);
                           noalias(N_DISP)= row(Ncontainer_Displacement, PointNumber);

                           if(StrainTensor.size1()!=dim || StrainTensor.size2()!=dim)
                                       StrainTensor.resize(dim,dim);
                           if(StressTensor.size1()!=dim || StressTensor.size2()!=dim)
                                       StressTensor.resize(dim,dim);

                           noalias(StressTensor) = ZeroMatrix(dim,dim);

                           GetPressures(N_PRESS, capillaryPressure, waterPressure);

                           porosity= GetPorosity(DN_DX_DISP);

							if(porosity >1.0 || porosity < 0.0)
							{
// 								rCurrentProcessInfo[ABORTED]= 1;
								return;
							}

                           density= GetAveragedDensity(capillaryPressure, porosity);

                           CalculateStressAndTangentialStiffnessUnsaturatedSoils(StressTensor, tanC_U, tanC_W, StrainTensor, DN_DX_DISP, waterPressure,  PointNumber, rCurrentProcessInfo);

                           noalias(du_dx)= CalculateDisplacementGradient(DN_DX_DISP);

                           if (CalculateStiffnessMatrixFlag == true)
                                    {
                                        //Calculation of spatial Stiffnes and Mass Matrix
                                          CalculateStiffnesMatrixUU(Help_K_UU, tanC_U, DN_DX_DISP, N_DISP, density, du_dx, capillaryPressure, Weight, DetJ, DetDef);

                                          CalculateStiffnesMatrixUW(Help_K_UW, tanC_W, DN_DX_DISP,N_DISP, N_PRESS, capillaryPressure, Weight, DetJ, DetDef);

                                          CalculateStiffnesMatrixWU(Help_K_WU, DN_DX_DISP, 
                                               DN_DX_PRESS, 
                                                  N_PRESS,du_dx,capillaryPressure,  Weight, DetJ, DetDef);

                                          CalculateStiffnesMatrixWW(Help_K_WW, DN_DX_DISP, 
                                                  DN_DX_PRESS,N_PRESS, capillaryPressure, Weight, DetJ, DetDef);
                                      }

                 	if (CalculateResidualVectorFlag == true) 
                 	{
                                        //Calculation of spatial Loadvector
                                       AddBodyForcesToRHSVectorU(Help_R_U, N_DISP, density, 
                                               Weight, DetJ , DetDef);
                                       AddInternalForcesToRHSU(Help_R_U, DN_DX_DISP, 
                                               StressTensor, Weight, DetJ, DetDef);
                                       AddInternalForcesToRHSW(Help_R_W, DN_DX_DISP, 
                                               DN_DX_PRESS,N_PRESS, capillaryPressure, Weight, 
                                               DetJ, DetDef);
                   	}

///////////////////////////////////////////////////////////////////////
// END Integration in space sum_(beta=0)^(number of quadrature points)
///////////////////////////////////////////////////////////////////////
                    }
                     if (CalculateStiffnessMatrixFlag == true)
                     { 
                         
                           AssembleTimeSpaceStiffnessFromStiffSubMatrices(rLeftHandSideMatrix,
                                      Help_K_UU, Help_K_UW, 
                                      Help_K_WU, Help_K_WW);
                     }

                     if (CalculateResidualVectorFlag == true)
                     {
                             AssembleTimeSpaceRHSFromSubVectors(rRightHandSideVector, Help_R_U,Help_R_W);
                     }
                       // // 
                KRATOS_CATCH("")
    }

	//************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
 
     void UnsaturatedSoilsElement_2phase::CalculateRightHandSide(VectorType& rRightHandSideVector, 
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
 
        void UnsaturatedSoilsElement_2phase::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
                VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
	   //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;
            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, 
                     CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);
        }
        
	////************************************************************************************
	////************************************************************************************
        
        void UnsaturatedSoilsElement_2phase::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY
            unsigned int number_of_nodes_disp = (mNodesDispMax-mNodesDispMin+1);
            unsigned int number_of_nodes_press = (mNodesPressMax-mNodesPressMin+1);
//             unsigned int number_of_nodes = number_of_nodes_disp+number_of_nodes_press;
            unsigned int dim = GetGeometry().WorkingSpaceDimension();

         //resizing as needed the LHS
            unsigned int MatSize1=(number_of_nodes_disp*dim+number_of_nodes_press);
            unsigned int MatSizeU= number_of_nodes_disp*dim;
            unsigned int MatSizeP= number_of_nodes_press;      


            if(rDampMatrix.size1() != MatSize1)
                rDampMatrix.resize(MatSize1,MatSize1);
            noalias(rDampMatrix) = ZeroMatrix(MatSize1,MatSize1); //resetting LHS

         //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =   
                    GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

//             const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =   
//                     mPressureGeometry.ShapeFunctionsLocalGradients(mThisIntegrationMethod);

//             const Matrix& Ncontainer_Displacement = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

            const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues(mThisIntegrationMethod);

        //auxiliary terms
//         InitializeGalerkinScheme(rCurrentProcessInfo);

            Matrix Help_D_UU(MatSizeU,MatSizeU);
            Matrix Help_D_UW(MatSizeU,MatSizeP);
                 
            Matrix Help_D_WU(MatSizeP,MatSizeU);
            Matrix Help_D_WW(MatSizeP,MatSizeP);

            Matrix DN_DX_DISP(number_of_nodes_disp, dim);

            double Weight;
 
            Vector N_PRESS(number_of_nodes_press);
                
            double capillaryPressure;

            double waterPressure;

            noalias(Help_D_UU)= ZeroMatrix(MatSizeU,MatSizeU);
            noalias(Help_D_UW)=ZeroMatrix(MatSizeU,MatSizeP);
                        
            noalias(Help_D_WU)= ZeroMatrix(MatSizeP,MatSizeU);
            noalias(Help_D_WW)=ZeroMatrix(MatSizeP,MatSizeP);

            double DetJ = 0.0;
            /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
            /////////////////////////////////////////////////////////////////////////
            for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
            {  
               	if(DN_DX_DISP.size1()!= (number_of_nodes_disp) || DN_DX_DISP.size2()!= dim)
                     DN_DX_DISP.resize(number_of_nodes_disp,dim);
				noalias(DN_DX_DISP)= ZeroMatrix(number_of_nodes_press, dim);

				noalias(DN_DX_DISP)= prod(DN_De_Displacement[PointNumber],mInvJ0[PointNumber]);

				double DetDef= Determinant_DeformationTensor(DN_DX_DISP);

                Weight=integration_points[PointNumber].Weight();
                        // Jacobian on current quadrature point

                DetJ= mDetJ0[PointNumber];
                // Shape Functions on current spatial quadrature point
 				if(N_PRESS.size()!=number_of_nodes_press)
 					N_PRESS.resize(number_of_nodes_press);
                noalias(N_PRESS)= row(Ncontainer_Pressure,PointNumber);

                GetPressures(N_PRESS, capillaryPressure, waterPressure);
//                            KRATOS_WATCH(waterPressure);
                           
//                         Calculation of spatial Stiffnes and Mass Matrix
                CalculateDampingMatrixWU(Help_D_WU, DN_DX_DISP, N_PRESS, 
                                              capillaryPressure, Weight, DetJ, DetDef);
                CalculateDampingMatrixWW( Help_D_WW, DN_DX_DISP, N_PRESS, 
                                               capillaryPressure, Weight, DetJ, DetDef);
//                          CalculateMassMatrix(HelpMassMatrix, N, Weight,DetJ,density);

                ///////////////////////////////////////////////////////////////////////
// END Integration in space sum_(beta=0)^(number of quadrature points)
                ///////////////////////////////////////////////////////////////////////
            }


            AssembleTimeSpaceStiffnessFromDampSubMatrices(rDampMatrix,
                        Help_D_UU, Help_D_UW, 
                        Help_D_WU, Help_D_WW);

            KRATOS_CATCH("")
        }
        
	////************************************************************************************
	////************************************************************************************
 
        void UnsaturatedSoilsElement_2phase::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
        {
            unsigned int number_of_nodes_disp = (mNodesDispMax-mNodesDispMin+1);
            unsigned int number_of_nodes_press = (mNodesPressMax-mNodesPressMin+1);
//             unsigned int number_of_nodes = number_of_nodes_disp+ number_of_nodes_press;
            unsigned int dim = GetGeometry().WorkingSpaceDimension();

         //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            Matrix DN_DX_DISP(number_of_nodes_disp, dim);

            const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =   
                    GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

//             std::vector<std::vector<Matrix> > tanC_U(dim,dim);
            array_1d<double,81> tanC_U;


            for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
            {
         
				//calculating and storing inverse of the jacobian and the parameters needed
				if(DN_DX_DISP.size1() != number_of_nodes_disp || DN_DX_DISP.size2() != dim)
                    DN_DX_DISP.resize(number_of_nodes_disp, dim);
                noalias(DN_DX_DISP) = ZeroMatrix(number_of_nodes_disp, dim);
                          
				noalias(DN_DX_DISP) = prod(DN_De_Displacement[PointNumber],mInvJ0[PointNumber]);

				Matrix StrainTensor(dim,dim);

				Matrix StressTensor(dim,dim);

                noalias(StrainTensor)= CalculateElasticNonlinearStrainTensorTrial(DN_DX_DISP, PointNumber);
                
                Vector StrainVector = SD_MathUtils<double>::TensorToStrainVector( StrainTensor );
                Vector StressVector = ZeroVector(6);
                SD_MathUtils<double>::TensorToVector( StressTensor, StressVector );
                Matrix tanC_U_matrix = ZeroMatrix(6,6);
                SD_MathUtils<double>::TensorToMatrix( tanC_U, tanC_U_matrix );

//              CalculateStressAndTangentialStiffness(StressTensor, tanC_U, StrainTensor); 
//                 mConstitutiveLawVector[PointNumber]->CalculateStressAndTangentMatrix( 
//                            StressTensor, StrainTensor, 
//                            tanC_U);	
                mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                    StrainVector,
                    ZeroMatrix(1),
                    StressVector,
                    tanC_U_matrix,
                    CurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod),PointNumber),
                    true,
                    1,
                    true );
                
                StrainTensor = SD_MathUtils<double>::StrainVectorToTensor(StrainVector);
                SD_MathUtils<double>::VectorToTensor(StressVector,StressTensor);
                SD_MathUtils<double>::MatrixToTensor( tanC_U_matrix, tanC_U );
            }
			//Calculate Darcy Velo at Nodes
	        Matrix DN_DX_PRESS_De(number_of_nodes_press, dim);
	        Matrix DN_DX_PRESS(number_of_nodes_press, dim);
			Matrix nodesLocalCoords(number_of_nodes_press, dim);
			
            nodesLocalCoords(0,0)=-1.0;
            nodesLocalCoords(0,1)=-1.0;
            nodesLocalCoords(0,2)=-1.0;

            nodesLocalCoords(1,0)=1.0;
            nodesLocalCoords(1,1)=-1.0;
            nodesLocalCoords(1,2)=-1.0;

            nodesLocalCoords(2,0)=1.0;
            nodesLocalCoords(2,1)=1.0;
            nodesLocalCoords(2,2)=-1.0;

            nodesLocalCoords(3,0)=-1.0;
            nodesLocalCoords(3,1)=1.0;
            nodesLocalCoords(3,2)=-1.0;

            nodesLocalCoords(4,0)=-1.0;
            nodesLocalCoords(4,1)=-1.0;
            nodesLocalCoords(4,2)=1.0;

            nodesLocalCoords(5,0)=1.0;
            nodesLocalCoords(5,1)=-1.0;
            nodesLocalCoords(5,2)=1.0;

            nodesLocalCoords(6,0)=1.0;
            nodesLocalCoords(6,1)=1.0;
            nodesLocalCoords(6,2)=1.0;

            nodesLocalCoords(7,0)=-1.0;
            nodesLocalCoords(7,1)=1.0;
            nodesLocalCoords(7,2)=1.0;

			for(unsigned int node=0; node< number_of_nodes_press; node++)
			{
				Vector local_coords(3);
				local_coords(0)=nodesLocalCoords(node,0);
				local_coords(1)=nodesLocalCoords(node,1);
				local_coords(2)=nodesLocalCoords(node,2);

				noalias(DN_DX_PRESS_De)= mpPressureGeometry->ShapeFunctionsLocalGradients(DN_DX_PRESS_De, local_coords);
                noalias(DN_DX_PRESS)= prod(DN_DX_PRESS_De,mInvJ0[0]);
				(GetGeometry()[node]).GetSolutionStepValue(REACTION_WATER_PRESSURE)=
										GetFlowWater(DN_DX_PRESS, (GetGeometry()[node]).GetSolutionStepValue(AIR_PRESSURE)-(GetGeometry()[node]).GetSolutionStepValue(WATER_PRESSURE))(2);
				(GetGeometry()[node]).GetSolutionStepValue(REACTION_AIR_PRESSURE)=
										GetSaturation(-(GetGeometry()[node]).GetSolutionStepValue(WATER_PRESSURE));
				//calculating and storing inverse of the jacobian and the parameters needed
				if(DN_DX_DISP.size1() != number_of_nodes_disp || DN_DX_DISP.size2() != dim)
                    DN_DX_DISP.resize(number_of_nodes_disp, dim);
                noalias(DN_DX_DISP) = ZeroMatrix(number_of_nodes_disp, dim);                         
				noalias(DN_DX_DISP) = prod(DN_De_Displacement[0],mInvJ0[0]);
				Matrix StrainTensor(dim,dim);
				Matrix StressTensor(dim,dim);
                noalias(StrainTensor)= CalculateElasticNonlinearStrainTensorTrial
                           (DN_DX_DISP, 0);
                
                Vector StrainVector = SD_MathUtils<double>::TensorToStrainVector( StrainTensor );
                Vector StressVector = ZeroVector(6);
                SD_MathUtils<double>::TensorToVector( StressTensor, StressVector );
                Matrix tanC_U_matrix = ZeroMatrix(6,6);
                SD_MathUtils<double>::TensorToMatrix( tanC_U, tanC_U_matrix );
                
//              CalculateStressAndTangentialStiffness(StressTensor, tanC_U, StrainTensor); 
//                 mConstitutiveLawVector[0]->CalculateStressAndTangentMatrix( 
//                            StressTensor, StrainTensor, 
//                            tanC_U);	
                mConstitutiveLawVector[0]->CalculateMaterialResponse(
                    StrainVector,
                    ZeroMatrix(1),
                    StressVector,
                    tanC_U_matrix,
                    CurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod),0),
                    true,
                    1,
                    true );
                
                StrainTensor = SD_MathUtils<double>::StrainVectorToTensor(StrainVector);
                SD_MathUtils<double>::VectorToTensor(StressVector,StressTensor);
                SD_MathUtils<double>::MatrixToTensor( tanC_U_matrix, tanC_U );
                
				(GetGeometry()[node]).GetSolutionStepValue(PRESSURE)=StressTensor(2,2);
			}


            for(unsigned int Point=0; Point< integration_points.size(); Point++)
            {
                mConstitutiveLawVector[Point]->FinalizeSolutionStep(GetProperties(), GetGeometry(), ZeroVector(0), CurrentProcessInfo);
            }
        }

	//************************************************************************************
	//************************************************************************************
	//************************************************************************************
        void UnsaturatedSoilsElement_2phase::CalculateOnIntegrationPoints(const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY

            unsigned int number_of_nodes_press = (mNodesPressMax-mNodesPressMin+1);

         //reading integration points and local gradients
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

            if(Output.size() != integration_points.size())
                Output.resize(integration_points.size());

            const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues(mThisIntegrationMethod);

            Vector N_PRESS(number_of_nodes_press);

            double capillaryPressure;
                        
            double waterPressure;

            double saturation;
            /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
            /////////////////////////////////////////////////////////////////////////
            for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
            {  
                // Shape Functions on current spatial quadrature point
                if(N_PRESS.size()!=number_of_nodes_press)
                   	N_PRESS.resize(number_of_nodes_press);
              	noalias(N_PRESS)= row(Ncontainer_Pressure,PointNumber);

                GetPressures(N_PRESS, capillaryPressure, waterPressure);

                saturation= GetSaturation(capillaryPressure);

                if(rVariable==SATURATION)
                {
                    Output[PointNumber]= saturation;
                }
				if(rVariable==WATER_PRESSURE)
				{
					Output[PointNumber]= waterPressure;
				}
				if(rVariable==AIR_PRESSURE)
				{
					Output[PointNumber]= 0.0;
				}
            }

            KRATOS_CATCH("")            
        }
        
        //************************************************************************************
        //************************************************************************************
        void UnsaturatedSoilsElement_2phase::CalculateOnIntegrationPoints(const Variable<Vector >& rVariable, std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo)
        {
            if(Output.size() != GetGeometry().IntegrationPoints(mThisIntegrationMethod).size())
                Output.resize(GetGeometry().IntegrationPoints(mThisIntegrationMethod).size());
            for(unsigned int ii = 0; ii<mConstitutiveLawVector.size(); ii++)
                Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
        }

        //************************************************************************************
	//************************************************************************************

        inline void UnsaturatedSoilsElement_2phase::CalculateAndAddExtForceContribution(const Vector& N, 
                const ProcessInfo& CurrentProcessInfo, Vector& BodyForce, VectorType& rRightHandSideVector, 
                double weight)
        {
            KRATOS_TRY

            KRATOS_CATCH("")
        }

        //************************************************************************************
	//************************************************************************************

        void UnsaturatedSoilsElement_2phase::EquationIdVector(EquationIdVectorType& rResult, 
                ProcessInfo& CurrentProcessInfo)
        {
            unsigned int dim_press = 1;//two pressure dofs 
            unsigned int dim_disp = (GetGeometry().WorkingSpaceDimension());//3 displacement dofs 
            unsigned int MatSize = 
                    (mNodesPressMax-mNodesPressMin+1)*dim_press+
                    (mNodesDispMax-mNodesDispMin+1)*dim_disp;
            if(rResult.size() != MatSize)
                    rResult.resize(MatSize);
            unsigned int adddisp= (mNodesPressMax-mNodesPressMin+1)*dim_press;

            for (unsigned int i=(mNodesPressMin-1);i<mNodesPressMax;i++)
            {
                    int index = (i-mNodesPressMin+1)*dim_press;
                    rResult[index] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
            }
            for (unsigned int i=(mNodesDispMin-1);i<mNodesDispMax;i++)
            {
                    int index = adddisp+(i-mNodesDispMin+1)*dim_disp;
                    rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
                    rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
                    rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
            }

        }

	//************************************************************************************
	//************************************************************************************
        
        void UnsaturatedSoilsElement_2phase::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& 
                CurrentProcessInfo)
        {
                ElementalDofList.resize(0);
                for (unsigned int i=(mNodesPressMin-1);i<mNodesPressMax;i++)
                {      
                    ElementalDofList.push_back(GetGeometry()[i].pGetDof(WATER_PRESSURE));    
                }
                for (unsigned int i=(mNodesDispMin-1);i<mNodesDispMax;i++)
                {      
                    ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                    ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
                    ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
                }
         }

	//************************************************************************************
	//************************************************************************************
        void UnsaturatedSoilsElement_2phase::GetValuesVector(Vector& values, int Step)
        {
                unsigned int dim_press = 1;//two pressure dofs 
                unsigned int dim_disp = (GetGeometry().WorkingSpaceDimension());//3 displacement dofs 
                unsigned int MatSize = 
                    (mNodesPressMax-mNodesPressMin+1)*dim_press+(mNodesDispMax-
                        mNodesDispMin+1)*dim_disp;
                if(values.size() != MatSize)
                    values.resize(MatSize);
                unsigned int adddisp= (mNodesPressMax-mNodesPressMin+1)*dim_press;

                for (unsigned int i=(mNodesPressMin-1);i<mNodesPressMax;i++)
                {
                        int index = (i-mNodesPressMin+1)*dim_press;
                        values(index) = 
                                GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE,Step); 
                }
             
                for (unsigned int i=(mNodesDispMin-1);i<mNodesDispMax;i++)
                {
                        int index =adddisp+(i-mNodesDispMin+1)*dim_disp;
                        
                        values(index) = 
                                GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X,Step);
                        values(index+1) = 
                                GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y,Step); 
                        values(index+2) = 
                                GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z,Step); 
                }

          }

	//************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
	//************************************************************************************  
         void UnsaturatedSoilsElement_2phase::AssembleTimeSpaceStiffnessFromDampSubMatrices
                 (MatrixType& rLeftHandSideMatrix, const Matrix& D_UU,const Matrix& D_UW,
                 const Matrix&  D_WU,const Matrix&  D_WW)
         {

                unsigned int dimension = GetGeometry().WorkingSpaceDimension();
                unsigned int number_of_nodes_disp= mNodesDispMax-mNodesDispMin+1;
                unsigned int number_of_nodes_press= mNodesPressMax-mNodesPressMin+1;
                unsigned int dim_disp= dimension;
                unsigned int dim_press= 1;
                unsigned int index_time_prim;
                unsigned int index_time_sec;
                unsigned int index_space_prim;
                unsigned int index_space_sec;
                unsigned int addIndex_disp= number_of_nodes_press*dim_press;

                for(unsigned int prim=0; prim<number_of_nodes_disp; prim++)
                {
                       for(unsigned int i=0; i< dim_disp; i++)
                       {
                              index_space_prim=prim*dim_disp+i;
                              
                              index_time_prim=addIndex_disp+prim*dim_disp+i;
                              
                              for(unsigned int sec=0; sec<number_of_nodes_disp; sec++)
                              {
                                    for(unsigned int j=0; j< dim_disp; j++)
                                    {
                                            index_space_sec=sec*dim_disp+j;
                                            
                                            index_time_sec=addIndex_disp+sec*dim_disp+j;
                                        
                                             rLeftHandSideMatrix(index_time_prim,
                                                     index_time_sec)
                                                  +=
                                                  (-1)*
                                                  D_UU(index_space_prim,index_space_sec);
                                     }
                                }
                         }
                   }

                   for(unsigned int prim=0; prim<number_of_nodes_disp; prim++)
                   {
                          for(unsigned int i=0; i< dim_disp; i++)
                          {
                                  index_space_prim=prim*dim_disp+i;

                                  index_time_prim=addIndex_disp+prim*dim_disp+i;
                                  
                                  for(unsigned int sec=0; sec<number_of_nodes_press; sec++)
                                  {
                                         index_space_sec=sec;

                                         index_time_sec=sec*dim_press;

                                         rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                            +=
                                                 (-1)* 
                                                 D_UW(index_space_prim,index_space_sec);
                                    }
                            }
                    }

                    for(unsigned int prim=0; prim<number_of_nodes_press; prim++)
                    {

                         index_space_prim=prim;
                         
                         index_time_prim= prim*dim_press;
                         
                         for(unsigned int sec=0; sec<number_of_nodes_disp; sec++)
                         {
                               for(unsigned int j=0; j< dim_disp; j++)
                               {
                                      index_space_sec=sec*dim_disp+j;

                                      index_time_sec=addIndex_disp+sec*dim_disp+j;                     

                                      rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                            +=
                                                (-1)* 
                                                D_WU(index_space_prim,index_space_sec);
                                 }
                            }
                    }

                    for(unsigned int prim=0; prim<number_of_nodes_press; prim++)
                    {
                            index_space_prim=prim;
                            
                            index_time_prim=prim*dim_press;

                            for(unsigned int sec=0; sec<number_of_nodes_press; sec++)
                            {
                                 index_space_sec=sec;

                                 index_time_sec=sec*dim_press;

                                 rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                    +=
                                        (-1)* 
                                        D_WW(index_space_prim,index_space_sec);
                              }
                     }
         }
        
         //************************************************************************************
	//************************************************************************************  
        
        void UnsaturatedSoilsElement_2phase::AssembleTimeSpaceStiffnessFromStiffSubMatrices
                (MatrixType& rLeftHandSideMatrix,const Matrix& K_UU,const Matrix& K_UW,
                const Matrix& K_WU,const Matrix& K_WW)
        {
             KRATOS_TRY
            
             unsigned int dimension = GetGeometry().WorkingSpaceDimension();
             unsigned int number_of_nodes_disp= mNodesDispMax-mNodesDispMin+1;
             unsigned int number_of_nodes_press= mNodesPressMax-mNodesPressMin+1;
             unsigned int dim_disp= dimension;
             unsigned int dim_press= 1;
             unsigned int index_time_prim;
             unsigned int index_time_sec;
             unsigned int index_space_prim;
             unsigned int index_space_sec;
             unsigned int addIndex_disp= number_of_nodes_press*dim_press;

             for(unsigned int prim=0; prim<number_of_nodes_disp; prim++)
             {
                 for(unsigned int i=0; i< dim_disp; i++)
                 {
                     index_space_prim=prim*dim_disp+i;
                              
                     index_time_prim=addIndex_disp+prim*dim_disp+i;
                              
                     for(unsigned int sec=0; sec<number_of_nodes_disp; sec++)
                     {
                         for(unsigned int j=0; j< dim_disp; j++)
                         {
                             index_space_sec=sec*dim_disp+j;
                                            
                             index_time_sec=addIndex_disp+sec*dim_disp+j;
                                        
                             rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                     +=
                                     (-1)*
                                     K_UU(index_space_prim,index_space_sec);
                         }
                     }
                 }
             }

             for(unsigned int prim=0; prim<number_of_nodes_disp; prim++)
             {
                 for(unsigned int i=0; i< dim_disp; i++)
                 {
                     index_space_prim=prim*dim_disp+i;

                     index_time_prim=addIndex_disp+prim*dim_disp+i;
                                  
                     for(unsigned int sec=0; sec<number_of_nodes_press; sec++)
                     {
                         index_space_sec=sec;

                         index_time_sec=sec*dim_press;

                         rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                 +=
                                 (-1)* 
                                 K_UW(index_space_prim,index_space_sec);
                     }
                 }
             }

             for(unsigned int prim=0; prim<number_of_nodes_press; prim++)
             {

                 index_space_prim=prim;
                         
                 index_time_prim= prim*dim_press;
                         
                 for(unsigned int sec=0; sec<number_of_nodes_disp; sec++)
                 {
                     for(unsigned int j=0; j< dim_disp; j++)
                     {
                         index_space_sec=sec*dim_disp+j;

                         index_time_sec=addIndex_disp+sec*dim_disp+j;                     

                         rLeftHandSideMatrix(index_time_prim,index_time_sec)
                                 +=
                                 (-1)* 
                                 K_WU(index_space_prim,index_space_sec);
                     }
                 }
             }

             for(unsigned int prim=0; prim<number_of_nodes_press; prim++)
             {
                 index_space_prim=prim;
                            
                 index_time_prim=prim*dim_press;

                 for(unsigned int sec=0; sec<number_of_nodes_press; sec++)
                 {
                     index_space_sec=sec;

                     index_time_sec=sec*dim_press;

                     rLeftHandSideMatrix(index_time_prim,index_time_sec)
                             +=
                             (-1)* 
                             K_WW(index_space_prim,index_space_sec);
                 }
             }
        KRATOS_CATCH("")
       }
        
        //************************************************************************************
	//************************************************************************************  
             
       void UnsaturatedSoilsElement_2phase::AssembleTimeSpaceRHSFromSubVectors(VectorType& 
               rRightHandSideVector,const Vector& R_U,const Vector& R_W)
       {
               KRATOS_TRY
            
            
               unsigned int dimension = GetGeometry().WorkingSpaceDimension();
               unsigned int number_of_nodes_disp= (mNodesDispMax-mNodesDispMin+1);
               unsigned int number_of_nodes_press= (mNodesPressMax-mNodesPressMin+1);
               unsigned int dim_disp= dimension;
               unsigned int dim_press= 1;
               unsigned int index_time;
               unsigned int addIndex_disp= number_of_nodes_press*dim_press;
            
//             noalias(rRightHandSideVector)= ZeroVector(number_of_nodes_press*dim_1+(number_of_nodes_disp-number_of_nodes_press)*dim_2);

               for(unsigned int prim=0; prim < number_of_nodes_disp; prim++)
               {
                      for(unsigned int i=0; i< dim_disp; i++)
                      {
                               index_time= addIndex_disp+prim*dim_disp+i;

                               rRightHandSideVector(index_time) +=
                                    R_U(prim*dim_disp+i);
                       }
               }
               for(unsigned int prim=0; prim < number_of_nodes_press; prim++)
               {
                                index_time= prim*dim_press;
 
                                rRightHandSideVector(index_time) +=
                                    R_W(prim);
                }

                KRATOS_CATCH("")

      }

        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
             
        //CALCULATE FORCEVECTORS DISPLACEMENT
             
       void UnsaturatedSoilsElement_2phase::AddBodyForcesToRHSVectorU(Vector& R, Vector& N_DISP, 
               double density, double Weight, double detJ, double DetDef)
       {
               KRATOS_TRY
            
               unsigned int dim= GetGeometry().WorkingSpaceDimension();

               Vector gravity(dim);
            
               noalias(gravity)=GetProperties()[GRAVITY];

               for(unsigned int prim=0; prim< (mNodesDispMax-mNodesDispMin+1); prim++)
               {
                      for(unsigned int i=0; i< dim; i++)
                      {
                               R(prim*dim+i) +=
                                    N_DISP(prim)*density*gravity(i)
									*Weight*DetDef*detJ;
                                    
                      }
                }

                KRATOS_CATCH("")
        }    
           
        //************************************************************************************
	//************************************************************************************  
             
        void UnsaturatedSoilsElement_2phase::AddInternalForcesToRHSU
                (Vector& R, const Matrix& DN_DX_DISP, 
                Matrix& StressTensor, double Weight, double detJ, double DetDef)
        {
                 KRATOS_TRY
             
                 unsigned int dim= GetGeometry().WorkingSpaceDimension();
            
                 for(unsigned int prim=0; prim< (mNodesDispMax-mNodesDispMin+1); prim++)
                 {
                        for(unsigned int i=0; i< dim; i++)
                        {
                               for(unsigned int gamma=0; gamma<dim; gamma++)
                               {
                                     R(prim*dim+i) +=
                                        (-1)*(DN_DX_DISP(prim,gamma)*StressTensor(i,gamma)*
                                        detJ*Weight);     
                                }
                         }
                    }
             
                    KRATOS_CATCH("")
          }

        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
             
        //CALCULATE STIFFNESS MATRICES DISPLACEMENT
         
        void UnsaturatedSoilsElement_2phase::CalculateStiffnesMatrixUU(Matrix& K,const 
                std::vector<std::vector<Matrix> >& C,const Matrix& DN_DX_DISP, Vector& N_DISP, 
               	double density, Matrix& du_dx, double capillaryPressure, double Weight, double detJ, double DetDef)
        {
			KRATOS_TRY
			unsigned int dim= GetGeometry().WorkingSpaceDimension();
            unsigned int number_of_nodes_disp= mNodesDispMax-mNodesDispMin+1;
            Vector Gravity(dim);
            noalias(Gravity)= GetProperties()[GRAVITY];

            double density_soil= GetProperties()[DENSITY];
            double density_air= GetProperties()[DENSITY_AIR];
       		double density_water= GetProperties()[DENSITY_WATER];
            double saturation=GetSaturation(capillaryPressure);

	       	double porosity_divu= GetDerivativeDPorosityDDivU(DN_DX_DISP);

        	double DrhoDdivU= 
            	porosity_divu*(-density_soil+(1-saturation)*density_air+saturation*density_water);

			vector<unsigned int> help;
			help.resize(2);

			for(unsigned int prim=0; prim< number_of_nodes_disp; prim++)
            {
				for(unsigned int i=0; i< dim; i++)
                {
           			for(unsigned int sec=0; sec< number_of_nodes_disp; sec++)
                    {
						for(unsigned int j=0; j< dim; j++)
                        {

                          	K(prim*dim+i,sec*dim+j) +=
                                N_DISP(prim)*DrhoDdivU*Gravity(i)
                                *DN_DX_DISP(sec,j)
								*Weight*DetDef*detJ;

							//BEGIN:Part ascociated with DefDet
							switch(j)
							{
								case 0:
									help[0]= 1;
									help[1]= 2;
									break;
								case 1:
									help[0]= 2;
									help[1]= 0;
									break;
								case 2:
									help[0]= 0;
									help[1]= 1;
									break;
								default:
									help[0]= 0;
									help[1]= 0;
								break;
							}
							double DefDet_U= 
								DN_DX_DISP(sec,j)*(1.0+du_dx(help[0],help[0]))*(1.0+du_dx(help[1],help[1]))
								+
								DN_DX_DISP(sec,help[0])*du_dx(help[0],help[1])*du_dx(help[1],j)
								+
								DN_DX_DISP(sec,help[1])*du_dx(help[0],j)*du_dx(help[1],help[0])
								-
								DN_DX_DISP(sec,j)*du_dx(help[0],help[1])*du_dx(help[1],help[0])
								-
								DN_DX_DISP(sec,help[0])*du_dx(help[1],help[1])*du_dx(help[0],j)
								-
								DN_DX_DISP(sec,help[1])*du_dx(help[0],help[0])*du_dx(help[1],j);

							 K(prim*dim+i,sec*dim+j) +=
									N_DISP(prim)*density*Gravity(i)*DefDet_U
									*Weight*detJ;
							//END:Part ascociated with DefDet
                            for(unsigned int alpha=0; alpha<dim; alpha++)
                            {
                            	for(unsigned int beta=0; beta<dim; beta++)
                                {
									K(prim*dim+i,sec*dim+j) +=
                                    	(-1)*DN_DX_DISP(prim,alpha)
                                        *(C[i][alpha](beta,j)+C[i][alpha](j,beta))
                                        *DN_DX_DISP(sec,beta)*detJ*Weight;
                                    for(unsigned int gamma=0; gamma<dim; gamma++)
                                    {
                            			K(prim*dim+i,sec*dim+j) +=
                                        	(-1)*DN_DX_DISP(prim,alpha)
                                          	*(C[i][alpha](beta,j)
                                           	+C[i][alpha](j,beta))
                                            *du_dx(beta,gamma)
                                            *DN_DX_DISP(sec,gamma)*detJ*Weight;
                              		}
                       			}
                       		}
     					}
                   	}
               	}
			}
			KRATOS_CATCH("")
         }
         
         void UnsaturatedSoilsElement_2phase::CalculateStiffnesMatrixUU(Matrix& K,const 
                 array_1d<double,81>& C,const Matrix& DN_DX_DISP, Vector& N_DISP, 
                 double density, Matrix& du_dx, double capillaryPressure, double Weight, double detJ, double DetDef)
                 {
                     KRATOS_TRY
                             unsigned int dim= GetGeometry().WorkingSpaceDimension();
                     unsigned int number_of_nodes_disp= mNodesDispMax-mNodesDispMin+1;
                     Vector Gravity(dim);
                     noalias(Gravity)= GetProperties()[GRAVITY];

                     double density_soil= GetProperties()[DENSITY];
                     double density_air= GetProperties()[DENSITY_AIR];
                     double density_water= GetProperties()[DENSITY_WATER];
                     double saturation=GetSaturation(capillaryPressure);

                     double porosity_divu= GetDerivativeDPorosityDDivU(DN_DX_DISP);

                     double DrhoDdivU= 
                                 porosity_divu*(-density_soil+(1-saturation)*density_air+saturation*density_water);

                     vector<unsigned int> help;
                     help.resize(2);

                     for(unsigned int prim=0; prim< number_of_nodes_disp; prim++)
                     {
                         for(unsigned int i=0; i< dim; i++)
                         {
                             for(unsigned int sec=0; sec< number_of_nodes_disp; sec++)
                             {
                                 for(unsigned int j=0; j< dim; j++)
                                 {

                                     K(prim*dim+i,sec*dim+j) +=
                                             N_DISP(prim)*DrhoDdivU*Gravity(i)
                                             *DN_DX_DISP(sec,j)
                                             *Weight*DetDef*detJ;

							//BEGIN:Part ascociated with DefDet
                                     switch(j)
                                     {
                                         case 0:
                                             help[0]= 1;
                                             help[1]= 2;
                                             break;
                                         case 1:
                                             help[0]= 2;
                                             help[1]= 0;
                                             break;
                                         case 2:
                                             help[0]= 0;
                                             help[1]= 1;
                                             break;
                                         default:
                                             help[0]= 0;
                                             help[1]= 0;
                                             break;
                                     }
                                     double DefDet_U= 
                                                 DN_DX_DISP(sec,j)*(1.0+du_dx(help[0],help[0]))*(1.0+du_dx(help[1],help[1]))
                                                 +
                                                 DN_DX_DISP(sec,help[0])*du_dx(help[0],help[1])*du_dx(help[1],j)
                                                 +
                                                 DN_DX_DISP(sec,help[1])*du_dx(help[0],j)*du_dx(help[1],help[0])
                                                 -
                                                 DN_DX_DISP(sec,j)*du_dx(help[0],help[1])*du_dx(help[1],help[0])
                                                 -
                                                 DN_DX_DISP(sec,help[0])*du_dx(help[1],help[1])*du_dx(help[0],j)
                                                 -
                                                 DN_DX_DISP(sec,help[1])*du_dx(help[0],help[0])*du_dx(help[1],j);

                                     K(prim*dim+i,sec*dim+j) +=
                                             N_DISP(prim)*density*Gravity(i)*DefDet_U
                                             *Weight*detJ;
							//END:Part ascociated with DefDet
                                     for(unsigned int alpha=0; alpha<dim; alpha++)
                                     {
                                         for(unsigned int beta=0; beta<dim; beta++)
                                         {
                                             K(prim*dim+i,sec*dim+j) +=
                                                     (-1)*DN_DX_DISP(prim,alpha)
                                                     *(C[27*i+9*alpha+3*beta+j]+C[27*i+9*alpha+3*j+beta])
                                                     *DN_DX_DISP(sec,beta)*detJ*Weight;
                                             for(unsigned int gamma=0; gamma<dim; gamma++)
                                             {
                                                 K(prim*dim+i,sec*dim+j) +=
                                                         (-1)*DN_DX_DISP(prim,alpha)
                                                         *(C[27*i+9*alpha+3*beta+j]
                                                                 +C[27*i+9*alpha+3*j+beta])
                                                         *du_dx(beta,gamma)
                                                         *DN_DX_DISP(sec,gamma)*detJ*Weight;
                                             }
                                         }
                                     }
                                 }
                             }
                         }
                     }
                     KRATOS_CATCH("")
                 }
        //************************************************************************************
	//************************************************************************************ 

        void UnsaturatedSoilsElement_2phase::CalculateStiffnesMatrixUW(Matrix& 
                Help_K_UW,Matrix& 
                tanC_W,const Matrix& DN_DX_DISP,Vector& N_DISP,Vector& N_PRESS, double 
                        capillaryPressure, 
                double Weight,double DetJ, double DetDef)
        {
                unsigned int dim= GetGeometry().WorkingSpaceDimension();
             
                unsigned int pressure_size= mNodesPressMax-mNodesPressMin+1;
             
                unsigned int displacement_size= mNodesDispMax-mNodesDispMin+1;
             
                Vector Gravity(dim);
             
                noalias(Gravity)= GetProperties()[GRAVITY];

                double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure);
             
                double porosity= GetPorosity(DN_DX_DISP);

                double density_air= GetProperties()[DENSITY_AIR];
                
                double DrhoDp_w= 
                            porosity*(density_air
                            -GetProperties()[DENSITY_WATER])*DSDpc;

                for(unsigned int prim=0; prim< displacement_size; prim++)
                {
                       for(unsigned int i=0; i<dim; i++)
                       {
                               for(unsigned int sec=0; sec< pressure_size; sec++)
                               {
                                     Help_K_UW(prim*dim+i,sec) +=
                                           N_DISP(prim)*DrhoDp_w*Gravity(i)*N_PRESS(sec)
                                             *Weight*DetDef*DetJ;
                                     for(unsigned int gamma=0; gamma<dim; gamma++)
                                     {
                                           Help_K_UW(prim*dim+i,sec) += 
                                               (-1)*(DN_DX_DISP(prim,gamma)*tanC_W(i,gamma)
                                                       *N_PRESS(sec))
                                                   *Weight*DetJ;
                                     }

                                 }
                         }
                  }
 
         }
         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
             
        //CALCULATE FORCEVECTORS WATER
         
        void UnsaturatedSoilsElement_2phase::AddInternalForcesToRHSW(Vector& Help_R_W,const 
                Matrix& DN_DX_DISP,const Matrix& DN_DX_PRESS, Vector& N_PRESS, double 
                        capillaryPressure, double 
                        Weight,double  DetJ, double DetDef)
       {
                unsigned int pressure_size=mNodesPressMax-mNodesPressMin+1;
             
                unsigned int dim= GetGeometry().WorkingSpaceDimension();
             
                double porosity= GetPorosity(DN_DX_DISP);

                double DS_Dpc= GetDerivativeDSaturationDpc(capillaryPressure); 

                double saturation= GetSaturation(capillaryPressure);

                double Dpc_Dt= GetDerivativeDCapillaryPressureDt(N_PRESS);

                double div_Dt= GetDerivativeDDivUDt(DN_DX_DISP);

                Vector flow_water(dim);
                noalias(flow_water)= GetFlowWater(DN_DX_PRESS,capillaryPressure);

                for(unsigned int prim=0; prim< pressure_size; prim++)
                {
                       Help_R_W(prim) += 
                            N_PRESS(prim)*porosity*DS_Dpc*Dpc_Dt
							*Weight*DetDef*DetJ*GetProperties()[SCALE];
                       Help_R_W(prim) += 
                           	N_PRESS(prim)*saturation*div_Dt
							*Weight*DetDef*DetJ*GetProperties()[SCALE];
                       for(unsigned int gamma=0; gamma<dim; gamma++)
                       {
                           Help_R_W(prim) += 
                                   (-1)*(DN_DX_PRESS(prim,gamma)*flow_water(gamma))
                                   *Weight*DetDef*DetJ*GetProperties()[SCALE];
                       }
                  }
         }

         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************   
         
          
        //CALCULATE STIFFNESS MATRICES WATER

        void UnsaturatedSoilsElement_2phase::CalculateStiffnesMatrixWU
                (Matrix& Help_K_WU,const Matrix& DN_DX_DISP,const Matrix& DN_DX_PRESS,Vector& N_PRESS, Matrix& du_dx,
                 double capillaryPressure,double Weight,double DetJ, double DetDef)
        {
			unsigned int dim= GetGeometry().WorkingSpaceDimension();
             
			unsigned int pressure_size= mNodesPressMax-mNodesPressMin+1;
             
			unsigned int displacement_size= mNodesDispMax-mNodesDispMin+1;

            double porosity= GetPorosity(DN_DX_DISP);

			double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure);

			double Dpc_Dt= GetDerivativeDCapillaryPressureDt(N_PRESS);
             
			double DnDdivU= GetDerivativeDPorosityDDivU(DN_DX_DISP);

            double saturation= GetSaturation(capillaryPressure);

            double div_Dt= GetDerivativeDDivUDt(DN_DX_DISP);

            Vector flow_water(dim);
            noalias(flow_water)= GetFlowWater(DN_DX_PRESS,capillaryPressure);

			vector<unsigned int> help;
			help.resize(2);

  			for(unsigned int prim=0; prim< pressure_size; prim++)
   			{
             	for(unsigned int sec=0; sec< displacement_size; sec++)
                {
              		for(unsigned int j=0; j<dim; j++)
                 	{
						//BEGIN:Part ascociated with DefDet
						switch(j)
						{
							case 0:
								help[0]= 1;
								help[1]= 2;
								break;
							case 1:
								help[0]= 2;
								help[1]= 0;
								break;
							case 2:
								help[0]= 0;
								help[1]= 1;
								break;
							default:
								help[0]= 0;
								help[1]= 0;
								break;
						}
						double DefDet_U= 
							DN_DX_DISP(sec,j)*(1.0+du_dx(help[0],help[0]))*(1.0+du_dx(help[1],help[1]))
							+
							DN_DX_DISP(sec,help[0])*du_dx(help[0],help[1])*du_dx(help[1],j)
							+
							DN_DX_DISP(sec,help[1])*du_dx(help[0],j)*du_dx(help[1],help[0])
							-
							DN_DX_DISP(sec,j)*du_dx(help[0],help[1])*du_dx(help[1],help[0])
							-
							DN_DX_DISP(sec,help[0])*du_dx(help[1],help[1])*du_dx(help[0],j)
							-
							DN_DX_DISP(sec,help[1])*du_dx(help[0],help[0])*du_dx(help[1],j);

						double HelpForceVecEntry_W= 0.0;

                       	HelpForceVecEntry_W += 
                            N_PRESS(prim)*porosity*DSDpc*Dpc_Dt
							*Weight*DetDef*DetJ*GetProperties()[SCALE];
                       	HelpForceVecEntry_W += 
                           	N_PRESS(prim)*saturation*div_Dt
							*Weight*DetDef*DetJ*GetProperties()[SCALE];

                       	for(unsigned int gamma=0; gamma<dim; gamma++)
                       	{
                           HelpForceVecEntry_W += 
                                   (-1)*(DN_DX_PRESS(prim,gamma)*flow_water(gamma))
                                   *Weight*DetDef*DetJ*GetProperties()[SCALE];
                       	}

                        Help_K_WU(prim,sec*dim+j) += HelpForceVecEntry_W*DefDet_U;
						//END:Part ascociated with DefDet

        				Help_K_WU(prim,sec*dim+j) +=
                  			N_PRESS(prim)*DnDdivU*DSDpc*Dpc_Dt*DN_DX_DISP(sec,j)
							*Weight*DetDef*DetJ*GetProperties()[SCALE];
                   	}
               	}
			}
        }
         
        //************************************************************************************
	//************************************************************************************
        void UnsaturatedSoilsElement_2phase::CalculateStiffnesMatrixWW
                (Matrix& Help_K_WW,const Matrix& DN_DX_DISP,const Matrix& DN_DX_PRESS,
                 Vector& N_PRESS, double capillaryPressure,double 
                        Weight,double DetJ, double DetDef)
        {
             unsigned int dim= GetGeometry().WorkingSpaceDimension();
             
             unsigned int pressure_size=mNodesPressMax-mNodesPressMin+1;
             
             double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure);
                     
             double porosity= GetPorosity(DN_DX_DISP);
             
             double D2S_Dpc2= GetSecondDerivativeD2SaturationDpc2(capillaryPressure);
             
             double Dpc_Dt= GetDerivativeDCapillaryPressureDt(N_PRESS);
             
             Vector Dflow_waterDpw(dim);
             
             noalias(Dflow_waterDpw)= GetDerivativeDWaterFlowDpw(DN_DX_PRESS, capillaryPressure);

             double Dflow_waterDgradpw= GetDerivativeDWaterFlowDGradpw(capillaryPressure);

             double Ddiv_Dt= GetDerivativeDDivUDt(DN_DX_DISP);
             
             for(unsigned int prim=0; prim< pressure_size; prim++)
             {
                    for(unsigned int sec=0; sec< pressure_size; sec++)
                    {
                            Help_K_WW(prim,sec) +=
                                 	(-1)*N_PRESS(prim)*porosity*D2S_Dpc2*Dpc_Dt*N_PRESS(sec)
                                    *Weight*DetDef*DetJ*GetProperties()[SCALE];
                         
                            Help_K_WW(prim,sec) +=
                                 	(-1)*N_PRESS(prim)*DSDpc*Ddiv_Dt*N_PRESS(sec)
									*Weight*DetDef*DetJ*GetProperties()[SCALE];

                            for(unsigned int gamma=0; gamma<dim; gamma++)
                            {
                                  Help_K_WW(prim,sec) +=
                                          	(-1)*DN_DX_PRESS(prim,gamma)*Dflow_waterDpw(gamma)
                                          	*N_PRESS(sec)
                                          	*Weight*DetDef*DetJ*GetProperties()[SCALE];
                                  Help_K_WW(prim,sec) +=
                                          	(-1)*DN_DX_PRESS(prim,gamma)*Dflow_waterDgradpw
                                          	*DN_DX_PRESS(sec,gamma)
											*Weight*DetDef*DetJ*GetProperties()[SCALE];
                             }

                      }
              }

        }
         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
              
        //CALCULATE DAMPING MATRICES WATER
        void UnsaturatedSoilsElement_2phase::CalculateDampingMatrixWU
                (Matrix& Help_D_WU,const Matrix& 
                DN_DX_DISP,Vector& N_PRESS, double capillaryPressure,double Weight,double DetJ, double DetDef)
        {
                unsigned int dim= GetGeometry().WorkingSpaceDimension();

                unsigned int pressure_size=mNodesPressMax-mNodesPressMin+1;
                
                unsigned int displacement_size=mNodesDispMax-mNodesDispMin+1;
                 
                double saturation= GetSaturation(capillaryPressure);

                for(unsigned int prim=0; prim< pressure_size; prim++)
                {
                      for(unsigned int sec=0; sec< displacement_size; sec++)
                      {
                           for(unsigned int j=0; j<dim; j++)
                           {
                                   Help_D_WU(prim,sec*dim+j) +=
                                           N_PRESS(prim)*saturation*DN_DX_DISP(sec,j)
                                           *Weight*DetDef*DetJ*GetProperties()[SCALE];
                           }
                       }
                }
        }
         
        //************************************************************************************
	//************************************************************************************
         
        void UnsaturatedSoilsElement_2phase::CalculateDampingMatrixWW(Matrix& Help_D_WW,const Matrix& 
                DN_DX_DISP,Vector& N_PRESS,double capillaryPressure,double Weight,double DetJ, double DetDef)
       {
               unsigned int pressure_size=mNodesPressMax-mNodesPressMin+1;
             
              double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure);
              double porosity= GetPorosity(DN_DX_DISP);
             
              for(unsigned int prim=0; prim< pressure_size; prim++)
              {
                     for(unsigned int sec=0; sec< pressure_size; sec++)
                     {
                             Help_D_WW(prim,sec) +=
                                     (-1)*N_PRESS(prim)*porosity*DSDpc*N_PRESS(sec)
                                     *Weight*DetDef*DetJ*GetProperties()[SCALE];
                     }
               }
        }

        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************

              
        //PRIMARY VARIABLES AND THEIR DERIVATIVES
         
         Matrix UnsaturatedSoilsElement_2phase::CalculateDisplacementGradient(const Matrix& DN_DX_DISP)
         { 
                  unsigned int dim= GetGeometry().WorkingSpaceDimension(); 
            
                  Matrix result(dim,dim);
                  noalias(result)= ZeroMatrix(dim,dim);
                    
                  Vector disp_alpha(dim);
                     
                  for(unsigned int point=(mNodesDispMin-1); point<mNodesDispMax; point++)
                  {
                      noalias(disp_alpha)= GetGeometry()[point].GetSolutionStepValue(DISPLACEMENT);
            
                           for(unsigned int j=0; j<3; j++)
                           {
                                  for(unsigned int k=0; k<dim; k++)
                                  {
                                      result(j,k) += disp_alpha(j)
                                              *DN_DX_DISP((point-mNodesDispMin+1),k);
                                   }
                            }
                   }

                   return result;
           }

         //************************************************************************************
	//************************************************************************************     
              
         Vector UnsaturatedSoilsElement_2phase::GetGradientWater(const Matrix& DN_DX_PRESS)
         {
                  unsigned int dim = GetGeometry().WorkingSpaceDimension();

                  Vector result(dim);
                  noalias(result)=ZeroVector(dim);
             
                  double presW_alpha;
                         
                  for(unsigned int i = mNodesPressMin-1 ; i < mNodesPressMax ;i++)
                  {

//                    //nodal Displacements
                      presW_alpha= GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE); 
                                
                      //thus strain in vector format

                      for(unsigned int k=0; k<3; k++)
                      {
                              result(k) +=
                                      presW_alpha
                                        *DN_DX_PRESS((i-mNodesPressMin+1),k); 
                      }
                   }
            
                   return result;
           }

        //************************************************************************************
	//************************************************************************************   
              
            void UnsaturatedSoilsElement_2phase::GetPressures(const Vector& N_PRESS, 
                    double& capillaryPressure, double& waterPressure)
            {
                    capillaryPressure= 0.0;
             
                    waterPressure= 0.0;

                    double presW_alpha;
                    
            //Calculating Strain on space points
                    for(unsigned int i = (mNodesPressMin-1) ; i < mNodesPressMax ;i++)
                    {

                          presW_alpha= GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE); 

                                //thus strain in vector format
                          capillaryPressure += 
                                  (-presW_alpha)
                                 *N_PRESS(i-mNodesPressMin+1);

                          waterPressure += 
                                  presW_alpha
                                  *N_PRESS(i-mNodesPressMin+1);

                      }
         }
         
        //************************************************************************************
	//************************************************************************************
         
         void UnsaturatedSoilsElement_2phase::GetDerivativeDPressuresDt(const Vector& N_PRESS, double& capillaryPressure_Dt, double& waterPressure_Dt)
        {
                      capillaryPressure_Dt= 0.0;
             
                      waterPressure_Dt= 0.0;

                      double presW_alpha_Dt;

            //Calculating Strain on space points
                      for(unsigned int i = mNodesPressMin-1 ; i < mNodesPressMax ;i++)
                      {
                            presW_alpha_Dt=
                                    GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE_DT);
                            
                            capillaryPressure_Dt += 
                                    (-presW_alpha_Dt)                      
                                    *N_PRESS(i-mNodesPressMin+1);

                            waterPressure_Dt += 
                                    presW_alpha_Dt
                                    *N_PRESS(i-mNodesPressMin+1);
                       }
              }
         
        //************************************************************************************
	//************************************************************************************
         
             double UnsaturatedSoilsElement_2phase::GetDerivativeDCapillaryPressureDt(const Vector& N_PRESS)
             {
                 double capillaryPressure_Dt= 0.0;

                 double presW_alpha_Dt;

            //Calculating Strain on space points
                 for(unsigned int i = mNodesPressMin-1 ; i < mNodesPressMax ;i++)
                 {
                     presW_alpha_Dt= GetGeometry()[i].GetSolutionStepValue(WATER_PRESSURE_DT);

         		capillaryPressure_Dt += 
                             (-presW_alpha_Dt)                      
                             *N_PRESS(i-mNodesPressMin+1);

                 }
                 
                 return capillaryPressure_Dt;
        }

        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
              
        //POROSITY AND ITS DERIVATIVES
         
       double UnsaturatedSoilsElement_2phase::GetPorosity(const Matrix& DN_DX_DISP)
       {
             double initialPorosity= GetProperties()[POROSITY];
             double div= GetDivU(DN_DX_DISP);

             double porosity= 1-(1-initialPorosity)*exp(-div);
             
//              if(porosity < 0)
//              {
// //             		KRATOS_ERROR(std::logic_error, "Porosity is less than zero" , *this);
//              }
//              if(porosity > 1)
//              {
// //             		KRATOS_ERROR(std::logic_error, "Porosity is bigger than one" , *this);
//              }

             return porosity;
        }

       double UnsaturatedSoilsElement_2phase::GetDerivativeDPorosityDDivU(const Matrix& DN_DX_DISP)
       {
             double initialPorosity= GetProperties()[POROSITY];
             double div= GetDivU(DN_DX_DISP);

             double porosity_divu= (1-initialPorosity)*exp(-div);
             

             return porosity_divu;
        }
         

         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
        
        Vector UnsaturatedSoilsElement_2phase::GetGravity()
        {
            unsigned int dim = GetGeometry().WorkingSpaceDimension();
            Vector gravity(dim);
            noalias(gravity)= GetProperties()[GRAVITY];
            return gravity;
        }
        
        //************************************************************************************
	//************************************************************************************  
        double UnsaturatedSoilsElement_2phase::GetDivU(const Matrix& DN_DX_DISP)
        {
                   unsigned int dim = GetGeometry().WorkingSpaceDimension();

                   double div= 0.0;
             
                   Vector u_alpha(dim);

                   for(unsigned int i = mNodesDispMin-1 ; i < mNodesDispMax ;i++)
                   {
//                     //nodal Displacements
                       noalias(u_alpha)= 
                               GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT);
                       //thus strain in vector format

                         for(unsigned int k=0; k<dim; k++)
                         {
                               div += (u_alpha(k))
                                        *DN_DX_DISP(i-mNodesDispMin+1,k); 
                         }
                    }       
                   
                    return div;
         }
        //************************************************************************************
	//************************************************************************************
        double UnsaturatedSoilsElement_2phase::GetDerivativeDDivUDt(const Matrix& DN_DX_DISP)
        {
                   unsigned int dim = GetGeometry().WorkingSpaceDimension();

                   double div= 0.0;
                   
                   Vector u_alpha_Dt (dim);
                   for(unsigned int i = mNodesDispMin-1 ; i < mNodesDispMax ;i++)
                   {
//                                //nodal Displacements
                        noalias(u_alpha_Dt)= 
                                GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_DT);
                                
                                //thus strain in vector format

                        for(unsigned int k=0; k<dim; k++)
                        {
                              div += u_alpha_Dt(k)*DN_DX_DISP(i-mNodesDispMin+1,k); 
                        }
                    }           
                    return div;
          }
        //AVERAGED DENSITY
        //************************************************************************************
	//************************************************************************************
          
         double UnsaturatedSoilsElement_2phase::GetAveragedDensity(double capillaryPressure,
                 double porosity)
         {
                 double result =0.0;

                 double density_soil= GetProperties()[DENSITY];
                 double density_air= GetProperties()[DENSITY_AIR];
                 double density_water= GetProperties()[DENSITY_WATER]; 
                 double saturation= GetSaturation(capillaryPressure);

                 result= (1-porosity)*density_soil+ 
                         porosity*(saturation*density_water+(1-saturation)*density_air);
             
                 return result;
          }
         
         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
              
        //SATURATION AND ITS DERIVATIVES
         
         double UnsaturatedSoilsElement_2phase::GetSaturation(double capillaryPressure)
         {
		double airEntryPressure= GetProperties()[AIR_ENTRY_VALUE];

		if(airEntryPressure<=0)
			airEntryPressure= 1.0;

		double b= GetProperties()[FIRST_SATURATION_PARAM];

		double c= GetProperties()[SECOND_SATURATION_PARAM];
             
                double saturation = 0.0;
//              
                if(capillaryPressure< 0.0)
                       capillaryPressure=0.0;

		saturation= pow((1.0+pow((capillaryPressure/airEntryPressure),b)),(-c));

// For Liakopolous Benchmark
// saturation =  1.0-1.9722*1e-11*pow(capillaryPressure,2.4279);
             

                return saturation;
          }
         
        //************************************************************************************
	//************************************************************************************
                  
          double UnsaturatedSoilsElement_2phase::GetDerivativeDSaturationDpc(double capillaryPressure)
          {
				double airEntryPressure= GetProperties()[AIR_ENTRY_VALUE];
				if(airEntryPressure<=0)
					airEntryPressure= 1.0;

				double b= GetProperties()[FIRST_SATURATION_PARAM];

				double c= GetProperties()[SECOND_SATURATION_PARAM];

                double result = 0.0;
//              
                if(capillaryPressure< 0.0)
                         capillaryPressure=0.0;

				result= (-c)*pow((1.0+pow((capillaryPressure/airEntryPressure),b)),(-c-1.0))*b* 
					pow((capillaryPressure/airEntryPressure),(b-1))*1.0/airEntryPressure;   
 
// For Liakopolous Benchmark 
// result =  -1.9722*2.4279*1e-11*pow(capillaryPressure,1.4279);

                return result;
           } 
         
        //************************************************************************************
	//************************************************************************************
         
            double UnsaturatedSoilsElement_2phase::GetSecondDerivativeD2SaturationDpc2(double capillaryPressure)
            {
				double airEntryPressure= GetProperties()[AIR_ENTRY_VALUE];
				if(airEntryPressure<=0)
					airEntryPressure= 1.0;

				double b= GetProperties()[FIRST_SATURATION_PARAM];

				double c= GetProperties()[SECOND_SATURATION_PARAM];
             
                double result = 0.0;

                if(capillaryPressure< 0.0)
                         capillaryPressure=0.0;
             
				result= (-c)*b/airEntryPressure*(
					(-c-1.0)*pow((1+pow((capillaryPressure/airEntryPressure),b)),(-c-2.0))
					*b/airEntryPressure*pow((capillaryPressure/airEntryPressure),(2.0*(b-1.0)))
					+pow((1.0+pow((capillaryPressure/airEntryPressure),b)),(-c-1.0))*(b-1.0)
					*pow((capillaryPressure/airEntryPressure),(b-2.0))*1.0/airEntryPressure);

// For Liakopolous Benschmark
// result =  -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure,0.4279);

                 return result;
              } 
         
        //************************************************************************************
	//************************************************************************************
        //************************************************************************************
	//************************************************************************************
              
        //WATER FLOW AND ITS DERIVATIVES

              Vector UnsaturatedSoilsElement_2phase::GetFlowWater(const Matrix& DN_DX_PRESS, double capillaryPressure)
              {
                   	unsigned int dim = GetGeometry().WorkingSpaceDimension();
             
					//Calculation of relative Permeability after Mualem
					double relPerm= GetSaturation(capillaryPressure);
					if(relPerm<= 0.01)
						relPerm= 0.01;

// For Liakopolous Benchmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

                       Vector gravity(dim);
                       noalias(gravity)= GetProperties()[GRAVITY];
             
                       Vector result(dim);
                       noalias(result)=ZeroVector(dim);
             
                       Vector grad_water(dim);
                       noalias(grad_water)= GetGradientWater(DN_DX_PRESS);
                       
                       for(unsigned int i=0; i< dim; i++)
                       {
                               result(i) = -relPerm*GetProperties()[PERMEABILITY_WATER]/
                                       (GetProperties()[DENSITY_WATER]*9.81)
                                       *(grad_water(i)-GetProperties()[DENSITY_WATER]
                                       *gravity(i));
                       }
             
                       return result;
              }

        //************************************************************************************
	//************************************************************************************
         
             Vector UnsaturatedSoilsElement_2phase::GetDerivativeDWaterFlowDpw(const Matrix& DN_DX_PRESS, double capillaryPressure)
             {
          		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		     	double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure); 
// 		    	Calculation of Derivative relative Permeability after Mualem
		     	double relPerm= GetSaturation(capillaryPressure);
		     	double relPerm_pw= (-1)*DSDpc;
				if(relPerm<= 0.0001)
				{
					relPerm= 0.0001;
					relPerm_pw= 0.0;
				}

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm_pw= -2.207*1.0121*pow((1-saturation),0.0121)*(-DSDpc)*(-1);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

                Vector result(dim);
                noalias(result)= ZeroVector(dim);

                Vector gravity(dim);
                noalias(gravity)= GetProperties()[GRAVITY];

                Vector grad_water(dim);
                noalias(grad_water)= GetGradientWater(DN_DX_PRESS);
             
                for(unsigned int i=0; i< dim; i++)
                {
                    result(i) = -relPerm_pw*GetProperties()[PERMEABILITY_WATER]/
               			(GetProperties()[DENSITY_WATER]*9.81)
                 		*(grad_water(i)-GetProperties()[DENSITY_WATER]
                        *gravity(i));
                }
              
                return result;
             }
         
        //************************************************************************************
	//************************************************************************************
         
             double UnsaturatedSoilsElement_2phase::GetDerivativeDWaterFlowDGradpw(double capillaryPressure)
             {
                     unsigned int dim = GetGeometry().WorkingSpaceDimension();

             		//Calculation of Derivative relative Permeability after Mualem
					double relPerm= GetSaturation(capillaryPressure);

	  				if(relPerm<= 0.0001)
					relPerm= 0.0001;
// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

                    double result(dim);

                    result = 
						-relPerm*GetProperties()[PERMEABILITY_WATER]/(GetProperties()[DENSITY_WATER]*9.81);

                    return result;
              }

        //************************************************************************************
	//************************************************************************************
        //STRESSES, STRAINS AND CONSTITUTIVE MODELL (UNSATURATED CASE)
        //************************************************************************************
        //************************************************************************************
         
         void UnsaturatedSoilsElement_2phase::CalculateEffectiveStress(Matrix& 
                 StressTensor, Matrix& tanC_W, const 
                 double waterPressure)
         {
             double capillaryPressure= -waterPressure;

             double saturation= GetSaturation(capillaryPressure);
             
             double DSDpc= GetDerivativeDSaturationDpc(capillaryPressure);
             
             Matrix Unity(3,3);
             noalias(Unity)=ZeroMatrix(3,3);
 
             for(unsigned int i=0; i<3; i++)
                 Unity(i,i)=1.0;
             
             for(unsigned int i=0; i<3; i++)
             {
                 for(unsigned int j=0; j<3; j++)
                 {
                     StressTensor(i,j)= 
                             StressTensor(i,j)-
                             (saturation*waterPressure)*Unity(i,j);
                 }
             }
             for(unsigned int i=0; i<3; i++)
             {
                 for(unsigned int j=0; j<3; j++)
                 {
                        tanC_W(i,j)=(DSDpc*waterPressure
                                                       -saturation)*Unity(i,j);

                 }
             }
             
          }
          
        //************************************************************************************
	//************************************************************************************
          /**
           * removed due to redefinition of material tensor to array_1d<double, 81>)
           */
/*          void UnsaturatedSoilsElement_2phase::CalculateStressAndTangentialStiffnessUnsaturatedSoils
                  (Matrix& StressTensor,std::vector<std::vector<Matrix> >& tanC_U,
                   Matrix& tanC_W, Matrix& StrainTensor,
                   const Matrix& DN_DX_DISP, double waterPressure, int PointNumber)
          {
                   KRATOS_TRY

                   unsigned int dim = GetGeometry().WorkingSpaceDimension();
                   if(tanC_W.size1()!=dim || tanC_W.size2()!=dim)
                           tanC_W.resize(dim,dim);
		   			noalias(tanC_W)= ZeroMatrix(dim,dim);

		   			if(StrainTensor.size1()!=dim || StrainTensor.size2()!=dim)
						StrainTensor.resize(dim,dim);
		  			if(StressTensor.size1()!=dim || StressTensor.size2()!=dim)
						StressTensor.resize(dim,dim);

                   	noalias(StrainTensor)= CalculateElasticNonlinearStrainTensorTrial
                           (DN_DX_DISP, PointNumber);

//              CalculateStressAndTangentialStiffness(StressTensor, tanC_U, StrainTensor); 

                   mConstitutiveLawVector[PointNumber]->CalculateStressAndTangentMatrix( 
                           StressTensor, StrainTensor, 
                           tanC_U);

                   CalculateEffectiveStress(StressTensor, tanC_W, waterPressure);

                   KRATOS_CATCH("")             
           }*/
           
           void UnsaturatedSoilsElement_2phase::CalculateStressAndTangentialStiffnessUnsaturatedSoils
                   (Matrix& StressTensor, array_1d<double,81>& tanC_U,
                    Matrix& tanC_W, Matrix& StrainTensor,
                    const Matrix& DN_DX_DISP, double waterPressure, int PointNumber, ProcessInfo& rCurrentProcessInfo)
                   {
                       KRATOS_TRY

                               unsigned int dim = GetGeometry().WorkingSpaceDimension();
                       if(tanC_W.size1()!=dim || tanC_W.size2()!=dim)
                           tanC_W.resize(dim,dim);
                       noalias(tanC_W)= ZeroMatrix(dim,dim);

                       if(StrainTensor.size1()!=dim || StrainTensor.size2()!=dim)
                           StrainTensor.resize(dim,dim);
                       if(StressTensor.size1()!=dim || StressTensor.size2()!=dim)
                           StressTensor.resize(dim,dim);

                       noalias(StrainTensor)= CalculateElasticNonlinearStrainTensorTrial
                               (DN_DX_DISP, PointNumber);

//              CalculateStressAndTangentialStiffness(StressTensor, tanC_U, StrainTensor); 
                       
                       Vector StrainVector = SD_MathUtils<double>::TensorToStrainVector( StrainTensor );
                       Vector StressVector = ZeroVector(6);
                       SD_MathUtils<double>::TensorToVector( StressTensor, StressVector );
                       Matrix tanC_U_matrix = ZeroMatrix(6,6);
                       SD_MathUtils<double>::TensorToMatrix( tanC_U, tanC_U_matrix );
                       
                       mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                               StrainVector,
                               ZeroMatrix(1),
                               StressVector,
                               tanC_U_matrix,
                               rCurrentProcessInfo,
                               GetProperties(),
                               GetGeometry(),
                               row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod),PointNumber),
                               true,
                               1,
                               true );
                       
                       StrainTensor = SD_MathUtils<double>::StrainVectorToTensor(StrainVector);
                       SD_MathUtils<double>::VectorToTensor(StressVector,StressTensor);
                       SD_MathUtils<double>::MatrixToTensor( tanC_U_matrix, tanC_U );

//                        mConstitutiveLawVector[PointNumber]->CalculateStressAndTangentMatrix( 
//                                StressTensor, StrainTensor, 
//                                tanC_U);

                       CalculateEffectiveStress(StressTensor, tanC_W, waterPressure);

                       KRATOS_CATCH("")             
                   }


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
               
           Matrix UnsaturatedSoilsElement_2phase::CalculateElasticNonlinearStrainTensorTrial
                   (const Matrix& DN_DX_DISP, int PointNumber)
         {
               KRATOS_TRY
                       
               unsigned int dim= GetGeometry().WorkingSpaceDimension();
                       
               Matrix LeftCauchyGreenTensorOld(dim,dim);
               noalias(LeftCauchyGreenTensorOld)=
                       mConstitutiveLawVector[PointNumber]->GetValue(ELASTIC_LEFT_CAUCHY_GREEN_OLD, LeftCauchyGreenTensorOld);

              //B=F*F_T
              //Calculate relative Deformation Gradient
              Matrix kronecker(3,3);
              noalias(kronecker)=ZeroMatrix(3,3);
                                   
               for(unsigned int i=0; i<3;i++)
               { 
                    kronecker(i,i)=1.0;
               }

               Matrix DefGrad(3,3);
               noalias(DefGrad)=ZeroMatrix(3,3);
                       
               Vector disp_alpha(3);
               Vector disp_old(3);

               for(unsigned int point=(mNodesDispMin-1); point<mNodesDispMax; point++)
               {
                   noalias(disp_alpha)= 
                           GetGeometry()[point].GetSolutionStepValue(DISPLACEMENT);

                   noalias(disp_old)= 
                           GetGeometry()[point].GetSolutionStepValue(DISPLACEMENT_OLD);

                   for(unsigned int j=0; j<3; j++)
                   {
                       for(unsigned int k=0; k<3; k++)
                       {
                           DefGrad(j,k) += (disp_alpha(j)-disp_old(j))*DN_DX_DISP((point-mNodesDispMin+1),k);
                       }
                   }
               }
               for(unsigned int i=0; i<dim; i++)
               {
                   DefGrad(i,i) += kronecker(i,i);
               }

              Matrix LeftCauchyGreenTensor(dim, dim);
              noalias(LeftCauchyGreenTensor)= ZeroMatrix(dim, dim);
              
              Matrix LeftCauchyGreenTensor_help(dim, dim);
              noalias(LeftCauchyGreenTensor_help)= ZeroMatrix(dim, dim);

              for(unsigned int i=0; i<dim; i++)
              {
                  for(unsigned int j=0; j<dim; j++)
                  {
                      for(unsigned int k=0; k<dim; k++)
                      {
                          LeftCauchyGreenTensor_help(i,j) += 
                                  DefGrad(i,k)*LeftCauchyGreenTensorOld(k,j);
                      }
                  }
              }
              
              for(unsigned int i=0; i<dim; i++)
              {
                  for(unsigned int j=0; j<dim; j++)
                  {
                      for(unsigned int k=0; k<dim; k++)
                      {
                          LeftCauchyGreenTensor(i,j) += 
                                  LeftCauchyGreenTensor_help(i,k)*DefGrad(j,k);
                      }
                  }
              }

              return LeftCauchyGreenTensor;
             
              KRATOS_CATCH("")
        }
        

        //************************************************************************************
	//************************************************************************************
        void UnsaturatedSoilsElement_2phase::InitializeMaterial
        ()
        {
                 KRATOS_TRY

                 for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
                 {
                     mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                     mConstitutiveLawVector[i]->InitializeMaterial( 
                             GetProperties(), GetGeometry(), ZeroVector(0));
                 }

                 KRATOS_CATCH("")
         } 
         
         void UnsaturatedSoilsElement_2phase::ResetConstitutiveLaw()
        {
            KRATOS_TRY
            for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
                    mConstitutiveLawVector[i]->ResetMaterial(GetProperties(), GetGeometry(), row(GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod), i) );
            KRATOS_CATCH("")
        }


 		void UnsaturatedSoilsElement_2phase::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
		{
			if(rValues.size() != mConstitutiveLawVector.size())
				rValues.resize(mConstitutiveLawVector.size());

            if(rVariable==ELASTIC_LEFT_CAUCHY_GREEN_OLD)
            {
				for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
                {
					if(rValues[i].size1() != 3 || rValues[i].size2() != 3 )
						rValues[i].resize(3,3);
					noalias(rValues[i])=mConstitutiveLawVector[i]->GetValue(ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i]);
				}
			}
		}

         void UnsaturatedSoilsElement_2phase::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
//              std::cout << mConstitutiveLawVector[0] << std::endl;
                if (rVariable==INSITU_STRESS)
                {
                    for( unsigned int PointNumber = 0; PointNumber<GetGeometry().IntegrationPoints(mThisIntegrationMethod).size(); PointNumber++ )
                        {
                                mConstitutiveLawVector[PointNumber]->SetValue(INSITU_STRESS, rValues[PointNumber],
                                                rCurrentProcessInfo );
                        }
                }
                if (rVariable==MATERIAL_PARAMETERS)
                {
                    for( unsigned int PointNumber = 0; PointNumber<GetGeometry().IntegrationPoints(mThisIntegrationMethod).size(); PointNumber++ )
                    {
                        mConstitutiveLawVector[PointNumber]->SetValue( MATERIAL_PARAMETERS,
                                rValues[PointNumber], rCurrentProcessInfo );
                    }
                }
        }

	void UnsaturatedSoilsElement_2phase::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() != mConstitutiveLawVector.size())
		{
			std::cout<<"wrong size: "<< rValues.size()<<"!="<<mConstitutiveLawVector.size()<<std::endl;
			return;
		}

		if(rVariable==ELASTIC_LEFT_CAUCHY_GREEN_OLD)
                {

			for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
                        {
						mConstitutiveLawVector[i]->SetValue(ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i], rCurrentProcessInfo);
			}		
		}
	}

	void UnsaturatedSoilsElement_2phase::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if(rValues.size() != mConstitutiveLawVector.size())
		{
			std::cout<<"wrong Size!: "<< rValues.size()<<"!="<<mConstitutiveLawVector.size()<<std::endl;
			return;
		}

		if(rVariable==INSITU_STRESS)
                {

			for( unsigned int i=0; i<mConstitutiveLawVector.size(); i++ )
                        {
						mConstitutiveLawVector[i]->SetValue(INSITU_STRESS, rValues[i], rCurrentProcessInfo);
			}		
		}
	}

	UnsaturatedSoilsElement_2phase::IntegrationMethod UnsaturatedSoilsElement_2phase::GetIntegrationMethod()
	{
			return mThisIntegrationMethod;
	}

	double UnsaturatedSoilsElement_2phase::Determinant_DeformationTensor(const Matrix& DN_DX_DISP)
	{
 		unsigned int dim= GetGeometry().WorkingSpaceDimension();

		//B=F*F_T
		//Calculate relative Deformation Gradient
		Matrix kronecker(3,3);
		noalias(kronecker)=ZeroMatrix(3,3);

                for(unsigned int i=0; i<3;i++)
                { 
 			kronecker(i,i)=1.0;
                }
                                   
		Matrix DefGrad(3,3);
		noalias(DefGrad)=ZeroMatrix(3,3);
                       
		Vector disp_alpha(3);

		for(unsigned int point=(mNodesDispMin-1); point<mNodesDispMax; point++)
		{
			noalias(disp_alpha)= 
				GetGeometry()[point].GetSolutionStepValue(DISPLACEMENT);

			for(unsigned int j=0; j<3; j++)
			{
				for(unsigned int k=0; k<3; k++)
				{
					DefGrad(j,k) += disp_alpha(j)*DN_DX_DISP((point-mNodesDispMin+1),k);
				}
			}
		}

		for(unsigned int i=0; i<dim; i++)
		{
			DefGrad(i,i) += kronecker(i,i);
		}

		double result=  DefGrad(0,0)*DefGrad(1,1)*DefGrad(2,2)
					+DefGrad(0,1)*DefGrad(1,2)*DefGrad(2,0)
					+DefGrad(0,2)*DefGrad(1,0)*DefGrad(2,1)
					-DefGrad(0,2)*DefGrad(1,1)*DefGrad(2,0)
					-DefGrad(0,0)*DefGrad(1,2)*DefGrad(2,1)
					-DefGrad(0,1)*DefGrad(1,0)*DefGrad(2,2);

		return result;
	}
} // Namespace Kratos


