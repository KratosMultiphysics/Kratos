//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Kazem Kamran
//                   Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/environment_contact.h"
#include "utilities/math_utils.h"
#include "thermo_mechanical_application.h"
#include "includes/convection_diffusion_settings.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
EnvironmentContact::EnvironmentContact(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
EnvironmentContact::EnvironmentContact(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{

}

Condition::Pointer EnvironmentContact::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new EnvironmentContact(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

EnvironmentContact::~EnvironmentContact()
{
}


//************************************************************************************
//************************************************************************************
void EnvironmentContact::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    MatrixType temp = Matrix();
// 		CalculateLocalSystem(temp, rRightHandSideVector,  rCurrentProcessInfo);

}

//************************************************************************************
//************************************************************************************
void EnvironmentContact::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

// 	  int nodes_number = GetGeometry().size();

    if(rLeftHandSideMatrix.size1() != 1)
    {
        rLeftHandSideMatrix.resize(1,1,false);
        rRightHandSideVector.resize(1,false);

    }

    rRightHandSideVector = ZeroVector(1);
//     array_1d<double,3> length_normal;
//     length_normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
//
// //     double length = norm_2(length_normal);
    double length =  GetGeometry()[0].FastGetSolutionStepValue(NODAL_PAUX);

    //take thermal properties
    ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rTransferCoefficientVar = my_settings->GetTransferCoefficientVariable();

    //         double HTC_Alpha = GetProperties()[rDiffusionVar];
    double HTC_Alpha = GetGeometry()[0].FastGetSolutionStepValue(rTransferCoefficientVar);
    HTC_Alpha *= length;

    if(length == 0.0)
        KRATOS_WATCH("NORMAL is ZERO")
    rLeftHandSideMatrix(0,0) = HTC_Alpha;

    //Residual
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
    array_1d<double, 1 > unknown_vec;
    unknown_vec[0] =  GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);

    double amb_T = rCurrentProcessInfo[AMBIENT_TEMPERATURE];
    unknown_vec[0] -= amb_T;

    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, unknown_vec);

}

//************************************************************************************
//************************************************************************************

void EnvironmentContact::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

// 	int nodes_number = 4;
// 	int dim = 2;
// 	unsigned int matsize = nodes_number*(dim);
//
// 	if(rDampingMatrix.size1() != matsize)
// 			rDampingMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!




    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void EnvironmentContact::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo,
                                      bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag)
{
    KRATOS_TRY


    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************
void EnvironmentContact::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void EnvironmentContact::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
    const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    }
    KRATOS_CATCH("");
}




} // Namespace Kratos


