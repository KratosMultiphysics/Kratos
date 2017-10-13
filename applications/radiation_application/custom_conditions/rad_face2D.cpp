//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


// Project includes 
#include "includes/define.h"
#include "custom_conditions/rad_face2D.h"
#include "utilities/math_utils.h"
#include "radiation_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  RadFace2D::RadFace2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
  {		
    //DO NOT ADD DOFS HERE!!!
  }
  
  //************************************************************************************
  //************************************************************************************
  RadFace2D::RadFace2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
  {
  }
  
  Condition::Pointer RadFace2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new RadFace2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }
  
  RadFace2D::~RadFace2D()
  {
  }
  
  
  //************************************************************************************
  //************************************************************************************
  void RadFace2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
    
    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
  }
  
  //************************************************************************************
  //************************************************************************************
  void RadFace2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
  }
  
  //************************************************************************************
  //************************************************************************************
  void RadFace2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag,bool CalculateResidualVectorFlag)
  {
    KRATOS_TRY
      
      
      unsigned int number_of_nodes = GetGeometry().size();
    //resizing as needed the LHS
    unsigned int MatSize=number_of_nodes;
    
    //calculate lenght
    double x21 = GetGeometry()[1].X() - GetGeometry()[0].X();
    double y21 = GetGeometry()[1].Y() - GetGeometry()[0].Y();
    
    double lenght = x21*x21 + y21*y21;
    lenght = sqrt(lenght);
    
    const Kratos::Condition::PropertiesType pproperties = GetProperties();
    //const double& 
    
    double StefenBoltzmann = 5.67e-8;
    double emissivity = pproperties[EMISSIVITY];
    emissivity = 1.0;
    
    const double& T0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
    const double& T1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);
    
    const double& I0 = GetGeometry()[0].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
    const double& I1 = GetGeometry()[1].FastGetSolutionStepValue(INCIDENT_RADIATION_FUNCTION);
    
    emissivity=1.0;
    
    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
      {
	if(rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize(MatSize,MatSize,false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize);	
	
	rLeftHandSideMatrix(0,0) = 1.0 * ( emissivity/(4.0-2.0*emissivity))* 0.5 * lenght; 
	rLeftHandSideMatrix(1,1) = 1.0 * ( emissivity/(4.0-2.0*emissivity))* 0.5 * lenght;
      }
    
    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
      {
	if(rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize(MatSize,false);
	
	rRightHandSideVector[0] = 1.0 * (emissivity*StefenBoltzmann*4.0 * pow(T0,4)) /(4.0-2.0*emissivity)- ( emissivity * I0/(4.0-2.0*emissivity)); 

	rRightHandSideVector[1] = 1.0 * (emissivity*StefenBoltzmann*4.0 * pow(T1,4)) /(4.0-2.0*emissivity)- ( emissivity * I1/(4.0-2.0*emissivity)); 
	
	rRightHandSideVector *= 0.5*lenght;
	//
      }

    KRATOS_CATCH("")
      }
  
  //************************************************************************************
  //************************************************************************************
  void RadFace2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
  {
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
      rResult.resize(number_of_nodes,false);
    for (unsigned int i=0;i<number_of_nodes;i++)
      {
	rResult[i] = (GetGeometry()[i].GetDof(INCIDENT_RADIATION_FUNCTION)).EquationId();
      }
  }
  
  //************************************************************************************
  //************************************************************************************
  void RadFace2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
  {
    ConditionalDofList.resize(GetGeometry().size());
    for (unsigned int i=0;i<GetGeometry().size();i++)
      {
		    ConditionalDofList[i] = (GetGeometry()[i].pGetDof(INCIDENT_RADIATION_FUNCTION));
      }
  }
  
} // Namespace Kratos


