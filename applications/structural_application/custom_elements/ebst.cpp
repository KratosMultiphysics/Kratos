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

/* **************************************************************************************
 *
 *   Last Modified by:    $Author: rrossi $
 *   Date:                $Date: 2008-10-13 07:00:53 $
 *   Revision:            $Revision: 1.12 $
 *
 * ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/ebst.h"
#include "includes/constitutive_law.h"
#include "structural_application.h"


namespace Kratos {

    namespace EbstAuxiliaries {
        Matrix msB(0, 0);
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > msQ = ZeroMatrix(3, 3);
        Matrix msD = ZeroMatrix(3, 3);
        Vector msStrainVector = ZeroVector(3);
        Vector msStressVector = ZeroVector(3);
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > msC = ZeroMatrix(2, 2);
        Matrix msDN_DX(0, 0);
    }

    using namespace EbstAuxiliaries;


    //***********************************************************************************
    //***********************************************************************************
    // -------- //
    //  PUBLIC  //
    // -------- //

    // Constructor

    Ebst::Ebst() {
    }

    // Constructor

    Ebst::Ebst(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
    }

    // Constructor

    Ebst::Ebst(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
        //initializing static variables
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = number_of_nodes * 3;
        msB.resize(3, dim);
        msDN_DX.resize(number_of_nodes, 3);
    }

    //***********************************************************************************
    //***********************************************************************************

    Element::Pointer Ebst::Create(
            IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new Ebst(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    //***********************************************************************************
    //***********************************************************************************
    // Destructor

    Ebst::~Ebst() {
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::EquationIdVector(
            EquationIdVectorType& rResult,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);

        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);
        unsigned int dim = number_of_nodes * 3;

        if (rResult.size() != dim)
            rResult.resize(dim, false);

        //nodes of the central element
        for (int i = 0; i < 3; i++) {
            int index = i*3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
        }

//        KRATOS_WATCH(Id());
//
//        for(unsigned int k = 0; k<rResult.size(); k++)
//            std::cout << rResult[k] << " ";
//        std::cout << std::endl;

        //adding the ids ofthe neighbouring nodes
        int index = 9;
        for (unsigned int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
                rResult[index] = neigb[i].GetDof(DISPLACEMENT_X).EquationId();
                rResult[index + 1] = neigb[i].GetDof(DISPLACEMENT_Y).EquationId();
                rResult[index + 2] = neigb[i].GetDof(DISPLACEMENT_Z).EquationId();
                index += 3;
            }
        }
        
//        for(unsigned int k = 0; k<rResult.size(); k++)
//            std::cout << rResult[k] << " ";
//        std::cout << std::endl;

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void Ebst::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        ElementalDofList.resize(0);

        //nodes of the central element
        for (unsigned int i = 0; i < GetGeometry().size(); i++)
        {
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        //adding the dofs ofthe neighbouring nodes
        for (unsigned int i = 0; i < 3; i++)
        {
            if (HasNeighbour(i,neigb[i]))
            {
                ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_X));
                ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_Y));
                  ElementalDofList.push_back(neigb[i].pGetDof(DISPLACEMENT_Z));
            }
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::GetValuesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            const array_1d<double, 3 > & disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            unsigned int index = i * 3;
            values[index] = disp[0];
            values[index + 1] = disp[1];
            values[index + 2] = disp[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
                const array_1d<double, 3 > & disp = neigb[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
                values[index] = disp[0];
                values[index + 1] = disp[1];
                values[index + 2] = disp[2];
                index += 3;
            }
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::GetFirstDerivativesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            const array_1d<double, 3 > & vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
            unsigned int index = i * 3;
            values[index] = vel[0];
            values[index + 1] = vel[1];
            values[index + 2] = vel[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
                const array_1d<double, 3 > & vel = neigb[i].FastGetSolutionStepValue(VELOCITY, Step);
                values[index] = vel[0];
                values[index + 1] = vel[1];
                values[index + 2] = vel[2];
                index += 3;
            }
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::GetSecondDerivativesVector(
            Vector& values,
            int Step) {
        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);
        unsigned int number_of_nodes = GetGeometry().size() + NumberOfActiveNeighbours(neigb);

        const unsigned int MatSize = number_of_nodes * 3;
        if (values.size() != MatSize)
            values.resize(MatSize);

        //nodes of the central element
        for (unsigned int i = 0; i < number_of_nodes; i++) {
            const array_1d<double, 3 > & acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
            unsigned int index = i * 3;
            values[index] = acc[0];
            values[index + 1] = acc[1];
            values[index + 2] = acc[2];
        }

        //neighbour nodes
        int index = 9;
        for (int i = 0; i < 3; i++) {
            if (HasNeighbour(i,neigb[i])) {
                const array_1d<double, 3 > & acc = neigb[i].FastGetSolutionStepValue(ACCELERATION, Step);
                values[index] = acc[0];
                values[index + 1] = acc[1];
                values[index + 2] = acc[2];
                index += 3;
            }
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo) {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;

        CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateOnIntegrationPoints(
            const Variable<Matrix>& rVariable,
            std::vector<Matrix>& Output,
            const ProcessInfo& rCurrentProcessInfo) {

    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::MassMatrix(
            MatrixType& rMassMatrix,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

        WeakPointerVector< Node<3> >& neigb = this->GetValue(NEIGHBOUR_NODES);

        //rMassMatrix.resize(0,0);
        // LUMPED MASS MATRIX
        unsigned int number_of_nodes = 3 + NumberOfActiveNeighbours(neigb);
        unsigned int MatSize = number_of_nodes * 3;
        if (rMassMatrix.size1() != MatSize)
            rMassMatrix.resize(MatSize, MatSize, false);
        rMassMatrix = ZeroMatrix(MatSize, MatSize);

        double Area = GetGeometry().Area();
        double TotalMass = Area * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
        Vector LumpFact;
        LumpFact = GetGeometry().LumpingFactors(LumpFact);

        for (unsigned int i = 0; i < number_of_nodes; i++) {
            double temp = LumpFact[i] * TotalMass;
            for (unsigned int j = 0; j < 3; j++) {
                unsigned int index = i * 3 + j;
                rMassMatrix(index, index) = temp;
            }
        }

        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::DampMatrix(
            MatrixType& rDampMatrix,
            ProcessInfo& rCurrentProcessInfo) {
        KRATOS_TRY

        if (rDampMatrix.size1() != 0)
            rDampMatrix.resize(0, 0, false);

        KRATOS_CATCH("")
    }

    //***********************************************************************************
    //***********************************************************************************

    void Ebst::FinalizeSolutionStep(
            ProcessInfo& rCurrentProcessInfo) {
    }





    //***********************************************************************************
    //***********************************************************************************

    void Ebst::CalculateAll(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag,
            bool CalculateResidualVectorFlag) {
        KRATOS_TRY

                const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int MatSize = number_of_nodes * 3;

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != MatSize)
                rLeftHandSideMatrix.resize(MatSize, MatSize);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize, MatSize); //resetting LHS
        }

        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != MatSize)
                rRightHandSideVector.resize(MatSize);
            rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
        }


        // LEFT HAND SIDE MATRIX
        if (CalculateStiffnessMatrixFlag == true) {
        }

        // RIGHT HAND SIDE VECTOR
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
        }

        KRATOS_CATCH("");
    }


    //***********************************************************************************
    //***********************************************************************************
    bool Ebst::HasNeighbour(unsigned int index, const Node < 3 > & neighb)
    {
        if (neighb.Id() == GetGeometry()[index].Id())
            return false;
        else
            return true;
    }

    //***********************************************************************************
    //***********************************************************************************

    unsigned int Ebst::NumberOfActiveNeighbours(WeakPointerVector< Node < 3 > >& neighbs)
    {
        unsigned int active_neighbours = 0;
        for (unsigned int i = 0; i < neighbs.size(); i++)
            if (HasNeighbour(i,neighbs[i])) active_neighbours++;
        return active_neighbours;
    }

    //***********************************************************************************
    //***********************************************************************************
    void Ebst::Initialize()
    {
        KRATOS_TRY
        //find the "nodal neighbours" given the elemental neighbours
        WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
        if(elem_neigb.size() == 0) KRATOS_ERROR(std::logic_error,"the neighbour elements are not calculated","")
        WeakPointerVector< Node<3> >& nodal_neigb = this->GetValue(NEIGHBOUR_NODES);
        nodal_neigb.resize(3);
        Geometry< Node<3> >& center_geom = GetGeometry();

        std::cout << "I am elem" << Id() <<std::endl;
        std::cout << "neighbours =" << elem_neigb[0].Id() << " " << elem_neigb[1].Id() << " " << elem_neigb[2].Id() << " " << std::endl;
        for (unsigned int i = 0; i < center_geom.size(); i++)
        {
            if(elem_neigb[i].Id() != Id() ) //if the elemental neighbour exists
            {
                Geometry< Node<3> >& geom = elem_neigb[i].GetGeometry();
                for (unsigned int j = 0; j < geom.size(); j++)
                {
                    bool aux = false;
                    for (unsigned int k = 0; k < center_geom.size(); k++)
                    {
                        if(geom[j].Id()==center_geom[k].Id())
                            aux = true;
                    }

                    if(aux == false) nodal_neigb(i) = Node<3>::WeakPointer( geom(j) );
                }
            }
            else //the elemenetal neighbour does not exist
	      nodal_neigb(i) = Node<3>::WeakPointer( center_geom(i) );
        }

        std::cout << "node1" << GetGeometry()[0].Id() << "opposite node =" << nodal_neigb[0].Id() << std::endl;
        std::cout << "node2" << GetGeometry()[1].Id() << "opposite node =" << nodal_neigb[1].Id() << std::endl;
        std::cout << "node3" << GetGeometry()[2].Id() << "opposite node =" << nodal_neigb[2].Id() << std::endl;

        KRATOS_CATCH("");
    }


} // Namespace Kratos.
