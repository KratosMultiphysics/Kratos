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

//
//   Project Name:        Kratos
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/irriducible_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mixedelement_application.h"
#include "utilities/geometry_utilities.h" 
//#include <omp.h>

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    IrriducibleElement::IrriducibleElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    IrriducibleElement::IrriducibleElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer IrriducibleElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new IrriducibleElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    IrriducibleElement::~IrriducibleElement()
    {
    }

    std::string IrriducibleElement::Info() const
    {
        std::stringstream buffer;
        buffer << "Linear Element" << std::endl;
        return buffer.str();
    }


    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::Initialize()
    {
        KRATOS_TRY

                unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int nintegration_points = GetGeometry().size(); //we will use a nodal integration rule
//        KRATOS_WATCH("ln114")
        if (dim == 2)
        {
            nintegration_points = 1;

            boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
            array_1d<double, 3 > N;

            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, mArea0);
            mDN_DX.resize(3, 2, false);
            noalias(mDN_DX) = DN_DX;
        } else
        {
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "");
        }


        //compute and save original area and shape functions

        //Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != nintegration_points)
        {
            mConstitutiveLawVector.resize(nintegration_points);
        }
        InitializeMaterial();
//        KRATOS_WATCH(mConstitutiveLawVector)
        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateAll(MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag,
            bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        const unsigned int nnodes = GetGeometry().size();

        unsigned int StrainSize;
        if (dim == 2)
            StrainSize = 3;
        else
            StrainSize = 6;

        //resizing LHS and RHS to the correct size
        const unsigned int var_block_size = dim;
        const unsigned int MatSize = var_block_size*nnodes;
        if(rLeftHandSideMatrix.size1() != MatSize)
            rLeftHandSideMatrix.resize(MatSize,MatSize);

         if(rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize);

        Matrix F;

        //compute B (it is constant over the whole element)
        Matrix B(StrainSize,dim*nnodes);
        CalculateB(B, mDN_DX, StrainSize);

        Vector discontinuous_strain(StrainSize);
        Vector stress_discontinuous(StrainSize);
        Vector u_list(dim*nnodes);
        GetValuesVector(u_list);
        noalias(discontinuous_strain) = prod(B,u_list);


        //nodal integration weight is used
        Matrix C(StrainSize,StrainSize);
        Vector N(nnodes,0.3333333333333);

        mConstitutiveLawVector[0]->CalculateMaterialResponse(discontinuous_strain,F,stress_discontinuous,C,rCurrentProcessInfo,GetProperties(),GetGeometry(),N,true,1,false);
        
        //compute block 22 (note that Cavg is used here)
        C *= mArea0;
        Matrix tmp = prod(C,B);
        noalias(rLeftHandSideMatrix) = prod(trans(B), tmp);

        stress_discontinuous *= mArea0;

        noalias(rRightHandSideVector) = -prod(trans(B),stress_discontinuous);

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp;

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo & rCurrentProcessInfo)
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);


    }



    ////************************************************************************************
    ////************************************************************************************

    void IrriducibleElement::InitializeSolutionStep(ProcessInfo & CurrentProcessInfo)
    {
//        double factor = 1.0 / static_cast<double> (geom.size());
        Vector N(GetGeometry().size());

        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            //nodal integration rule is used
            for (unsigned int j = 0; j < N.size(); j++)
            {
                 N[j] = 1.0/static_cast<double>(GetGeometry().size());
            }
            mConstitutiveLawVector[i]->InitializeSolutionStep(GetProperties(),
                    GetGeometry(), N,
                    CurrentProcessInfo);
        }
    }

    ////************************************************************************************
    ////************************************************************************************

    void IrriducibleElement::FinalizeSolutionStep(ProcessInfo & CurrentProcessInfo)
    {
        Vector N(GetGeometry().size());
        for (unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
        {
            //nodal integration rule is used
            for (unsigned int j = 0; j < N.size(); j++)
            {
                 N[j] = 1.0/static_cast<double>(GetGeometry().size());
            }

            mConstitutiveLawVector[i]->FinalizeSolutionStep(GetProperties(),
                    GetGeometry(),
                    N,
                    CurrentProcessInfo);
        }
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::InitializeMaterial()
    {
        KRATOS_TRY

        Vector N(GetGeometry().size());

        if (GetProperties()[CONSTITUTIVE_LAW] != NULL)
        {
            
                for (unsigned int j = 0; j < N.size(); j++)
                {
                     N[j] = 1.0/static_cast<double>(GetGeometry().size());
                }

                mConstitutiveLawVector[0] = GetProperties()[CONSTITUTIVE_LAW]->Clone();

                mConstitutiveLawVector[0]->InitializeMaterial(GetProperties(), GetGeometry(), N);
            
        } else
            KRATOS_ERROR(std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id())

            KRATOS_CATCH("")
        }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateB(
            Matrix& B,
            Matrix& DN_DX,
            unsigned int StrainSize)
    {
        KRATOS_TRY
                const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();


        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            unsigned int index = dimension * i;

            if (dimension == 2)
            {
                B(0, index + 0) = DN_DX(i, 0);
                B(0, index + 1) = 0.0;
                B(1, index + 0) = 0.0;
                B(1, index + 1) = DN_DX(i, 1);
                B(2, index + 0) = DN_DX(i, 1);
                B(2, index + 1) = DN_DX(i, 0);
            } else
            {
                KRATOS_ERROR(std::logic_error, "3d not yet implemented", "");
            }
        }

        KRATOS_CATCH("")
    }



    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo & CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (dim == 2)
        {
            unsigned int block_size = 2;
            unsigned int MatSize = number_of_nodes * block_size;
            if (rResult.size() != MatSize) rResult.resize(MatSize, false);

            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                unsigned int index = i * block_size;

                rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            }
        } else
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "");

    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo & CurrentProcessInfo)
    {
        ElementalDofList.resize(0);
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (dim == 2)
        {
            for (unsigned int i = 0; i < GetGeometry().size(); i++)
            {
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

            }
        } else
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "");
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::MassMatrix(MatrixType& rMassMatrix, ProcessInfo & rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::DampMatrix(MatrixType& rDampMatrix, ProcessInfo & rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& Output, const ProcessInfo & rCurrentProcessInfo)
    {
        if (Output.size() != mConstitutiveLawVector.size())
            Output.resize(mConstitutiveLawVector.size(), false);

        for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
            Output[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, Output[ii]);
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& Output, const ProcessInfo & rCurrentProcessInfo)
    {
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int StrainSize;

        if (dim == 2)
            StrainSize = 3;
        else
            StrainSize = 6;

        Vector StrainVector(StrainSize);

        if (Output.size() != mConstitutiveLawVector.size())
            Output.resize(mConstitutiveLawVector.size());

        for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
            Output[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, Output[ii]);

    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo & rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo & rCurrentProcessInfo)
    {
        for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            mConstitutiveLawVector[PointNumber]->SetValue(rVariable,
                    rValues[PointNumber], rCurrentProcessInfo);
        }

    }


    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo & rCurrentProcessInfo)
    {
        for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++)
        {
            mConstitutiveLawVector[PointNumber]->SetValue(rVariable,
                    rValues[PointNumber], rCurrentProcessInfo);
        }

    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo & rCurrentProcessInfo)
    {
        if (rValues.size() != mConstitutiveLawVector.size())
            rValues.resize(mConstitutiveLawVector.size(), false);

        for (unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++)
            rValues[ii] = mConstitutiveLawVector[ii]->GetValue(rVariable, rValues[ii]);
    }


    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo & rCurrentProcessInfo)
    {

        for (unsigned int PointNumber = 0;
             PointNumber < mConstitutiveLawVector.size();
             PointNumber++)
        {
            rValues[PointNumber] =
                    mConstitutiveLawVector[PointNumber]->GetValue(rVariable, rValues[PointNumber]);
        }

    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
            std::vector<Matrix>& rValues, const ProcessInfo & rCurrentProcessInfo)
    {
        CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetValuesVector(Vector& values, int Step)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (dim == 2)
        {
            unsigned int block_size = 2;
            unsigned int MatSize = number_of_nodes * block_size;
            if (values.size() != MatSize) values.resize(MatSize, false);

            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                unsigned int index = i * block_size;
                values[index ] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X, Step);
                values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y, Step);
            }
        } else
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "");
    }


    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetFirstDerivativesVector(Vector& values, int Step)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (dim == 2)
        {
            unsigned int block_size = 2;
            unsigned int MatSize = number_of_nodes * block_size;
            if (values.size() != MatSize) values.resize(MatSize, false);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                unsigned int index = i * block_size;
                values[index ] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X, Step);
                values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y, Step);
            }
        } else
        {
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "")
        }


    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::GetSecondDerivativesVector(Vector& values, int Step)
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if (dim == 2)
        {
            unsigned int block_size = 2;
            unsigned int MatSize = number_of_nodes * block_size;
            if (values.size() != MatSize) values.resize(MatSize, false);
            for (unsigned int i = 0; i < number_of_nodes; i++)
            {
                unsigned int index = i * block_size;
                values[index ] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_X, Step);
                values[index + 1] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Y, Step);
            }
        } else
        {
            KRATOS_ERROR(std::logic_error, "3d not yet implemented", "")
        }
    }

    //************************************************************************************
    //************************************************************************************

    void IrriducibleElement::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo & rCurrentProcessInfo)
    {
    }



} // Namespace Kratos


