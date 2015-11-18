/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-02-14 09:41:09 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/fluid3d_neumann.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "includes/deprecated_variables.h"
#include "includes/cfd_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************

Fluid3DNeumann::Fluid3DNeumann(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

Fluid3DNeumann::Fluid3DNeumann(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer Fluid3DNeumann::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new Fluid3DNeumann(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

Fluid3DNeumann::~Fluid3DNeumann()
{
}


//************************************************************************************
//************************************************************************************

void Fluid3DNeumann::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    const unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

    if (FractionalStepNumber < 4)
    {
        if (rRightHandSideVector.size() != 3)
            rRightHandSideVector.resize(3, false);

        //calculate normal to element
        array_1d<double, 3 > An, v1, v2;
        v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
        v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
        v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();

        v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
        v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
        v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

        MathUtils<double>::CrossProduct(An, v1, v2);
        An *= 0.5;

        unsigned int is_structure = this->GetValue(IS_STRUCTURE);

        if (is_structure != 1.0)
        {
            unsigned int component = FractionalStepNumber - 1;
            for (unsigned int i = 0; i < GetGeometry().size(); i++)
                rRightHandSideVector[i] = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE) * An[component] * 0.3333333333333333;
        }
        else
        {
            unsigned int component = FractionalStepNumber - 1;
            const double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            const double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[0] = -(2.0 * p0 + p1 + p2) * An[component] / 12.0;
            rRightHandSideVector[1] = -(p0 + 2.0 * p1 + p2) * An[component] / 12.0;
            rRightHandSideVector[2] = -(p0 + p1 + 2.0 * p2) * An[component] / 12.0;

        }

    }
    if (FractionalStepNumber == 4) //pressure step
    {
        if (rRightHandSideVector.size() != 3)
            rRightHandSideVector.resize(3, false);

        noalias(rRightHandSideVector) = ZeroVector(3);

    }
    else
    {
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }
}

//************************************************************************************
//************************************************************************************

void Fluid3DNeumann::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const unsigned int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];


    //calculate normal to element
    array_1d<double, 3 > An, v1, v2;
    v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
    v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
    v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();

    v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
    v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
    v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

    MathUtils<double>::CrossProduct(An, v1, v2);
    An *= 0.5;

    if (FractionalStepNumber < 4)
    {
        if (rLeftHandSideMatrix.size1() != 3)
        {
            rLeftHandSideMatrix.resize(3, 3, false);
            rRightHandSideVector.resize(3, false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(3, 3);
        }

        unsigned int is_structure = this->GetValue(IS_STRUCTURE);

        if (is_structure == 1.0)
        {
            unsigned int component = FractionalStepNumber - 1;
            const double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            const double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            rRightHandSideVector[0] = -(2.0 * p0 + p1 + p2) * An[component] / 12.0;
            rRightHandSideVector[1] = -(p0 + 2.0 * p1 + p2) * An[component] / 12.0;
            rRightHandSideVector[2] = -(p0 + p1 + 2.0 * p2) * An[component] / 12.0;
        }
        else
        {
            unsigned int component = FractionalStepNumber - 1;
            for (unsigned int i = 0; i < GetGeometry().size(); i++)
                rRightHandSideVector[i] = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE) * An[component] * 0.3333333333333333;
        }
    }
    else if (FractionalStepNumber == 4)   //pressure step
    {
        if (rLeftHandSideMatrix.size1() != 3)
        {
            rLeftHandSideMatrix.resize(3, 3, false);
            rRightHandSideVector.resize(3, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(3, 3);
        }

        unsigned int is_structure = this->GetValue(IS_STRUCTURE);

        if (is_structure == 1)
        {
            double vn0 = inner_prod(An, GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL));
            double vn1 = inner_prod(An, GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL));
            double vn2 = inner_prod(An, GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL));
            double v0 = 2.0 * vn0 + vn1 + vn2;
            double v1 = vn0 + 2.0 * vn1 + vn2;
            double v2 = vn0 + vn1 + 2.0 * vn2;

            rRightHandSideVector[0] = v0 / 12.0;
            rRightHandSideVector[1] = v1 / 12.0;
            rRightHandSideVector[2] = v2 / 12.0;

        }
        else
            noalias(rRightHandSideVector) = ZeroVector(3);
    }
    else
    {
        if (rLeftHandSideMatrix.size1() != 0)
        {
            rLeftHandSideMatrix.resize(0, 0, false);
            rRightHandSideVector.resize(0, false);
        }
    }
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

void Fluid3DNeumann::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    const unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

    if (FractionalStepNumber < 4)
    {
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
        else if (FractionalStepNumber == 3) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Z).EquationId();
    }
    else if (FractionalStepNumber == 4)
    {
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }
    else if (rResult.size() != 0)
        rResult.resize(0, false);

}

//************************************************************************************
//************************************************************************************

void Fluid3DNeumann::GetDofList(DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo)
{
    unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];


    if (FractionalStepNumber < 4)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (ConditionalDofList.size() != number_of_nodes)
            ConditionalDofList.resize(number_of_nodes);

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ConditionalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ConditionalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
        else if (FractionalStepNumber == 3) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ConditionalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Z);
    }
    if (FractionalStepNumber == 4) //pressure step
    {
        if (ConditionalDofList.size() != GetGeometry().size())
            ConditionalDofList.resize(GetGeometry().size());
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            ConditionalDofList[i] = (GetGeometry()[i].pGetDof(PRESSURE));
        }
    }
    else if (ConditionalDofList.size() != 0)
        ConditionalDofList.resize(0);
}


//************************************************************************************
//************************************************************************************

void Fluid3DNeumann::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int is_structure = this->GetValue(IS_STRUCTURE);

    if (is_structure == 1.0)
    {

        int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
        if (FractionalStepNumber == 6) //calculation of stabilization terms
        {

            //getting data for the given geometry
            //calculate normal to element
            array_1d<double, 3 > An, v1, v2;
            v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
            v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
            v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();

            v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
            v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
            v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

            MathUtils<double>::CrossProduct(An, v1, v2);
            An *= 0.5;


            double dp0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE) - GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            double dp1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE) - GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            double dp2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE) - GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);

            double p0 = 2.0 * dp0 + dp1 + dp2;
            double p1 = dp0 + 2.0 * dp1 + dp2;
            double p2 = dp0 + dp1 + 2.0 * dp2;

            p0 /= 12.0;
            p1 /= 12.0;
            p2 /= 12.0;

            GetGeometry()[0].SetLock();
            noalias(GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL)) += -p0 * An;
            GetGeometry()[0].UnSetLock();

            GetGeometry()[1].SetLock();
            noalias(GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL)) += -p1 * An;
            GetGeometry()[1].UnSetLock();

            GetGeometry()[2].SetLock();
            noalias(GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL)) += -p2 * An;
            GetGeometry()[2].UnSetLock();
        }
    }

    KRATOS_CATCH("");
}

} // Namespace Kratos


