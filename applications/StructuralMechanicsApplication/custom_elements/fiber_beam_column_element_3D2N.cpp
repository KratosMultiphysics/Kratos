// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Long Chen
//
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/fiber_beam_column_element_3D2N.hpp"
#include "includes/define.h"
#include "structural_machanics_application_variables.h"


namespace Kratos {

FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId,
                                                        GeometryType::Pointer pGeometry)
                                                        :Element(NewId, pGeometry) {}

FiberBeamColumnElement3D2N::FiberBeamColumnElement3D2N(IndexType NewId,
                                                        GeometryType::Pointer pGeometry,
                                                        PropertiesType::Pointer pProperties)
                                                        :Element(NewId, pGeometry, pProperties) {}

Element::Pointer FiberBeamColumnElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
PropertiesType::Pointer pProperties) const {
    return Kratos::make_shared<FiberBeamColumnElement3D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}


Element::Pointer FiberBeamColumnElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const {
    return Kratos::make_shared<FiberBeamColumnElement3D2N>(NewId, pGeom, pProperties);
}

FiberBeamColumnElement3D2N::~FiberBeamColumnElement3D2N() {}

void FiberBeamColumnElement3D2N::EquationIdVector(EquationIdVector &rResult,
ProcessInfo &rCurrentProcessInfo) {
    if(rResult.size() != msElementSize)
    rResult.resize(msElementSize);

    for (int i = 0; i < msNumberOfNodes; ++i)
    {
        int index = i * msLocalSize;
        rResult[index] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index + 2] = this->Getgeometry()[i].GetDof(DISPLACEMENT_Z).EquaitonId();
        rResult[index + 3] = this->GetGeometry()[i].GetDof(ROTATION_X).EquationId();
        rResult[index + 4] = this->GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
    }
}

void FiberBeamColumnElement3D2N::GetDofList(DofsVectorType &rElementalDofList,
ProcessInfo &rCurrentProcessInfo)
{
    if(rElementalDofList.size() != msElementSize) {
        rElementalDofList.resize(msElementSiez);
    }

    for(int i = 0; i < msNumberOfNodes; ++i){
        int index = i * msLocalSize;
        rElementalDofList[index] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index + 3] = this->GetGeometry()[i].pGetDof(ROTATION_X);
        rElementalDofList[index + 4] = this->GetGeometry()[i].pGetDof(ROTATION_Y);
    }
}

void FiberBeamColumnElement3D2N::Initialize() {
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void FiberBeamColumnElement3D2N::getDeformationValuesVector(Vector &rValues, int Step){
    KRATOS_TRY;
    if (rValue.size() != msElementSize)
        rValue.resize(msElementSize, false);

    for(int i = 0; i < msNumberOfNodes; ++i)
    {
        int index = i * msLocalSize;
        const auto &disp = 
        this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const auto &rot =
        this->GetGeometry()[i].FastGetSolutionStepValue(ROTATION, Step)

        rValues[index] = disp[0];
        rValues[index + 1] = disp[1];
        rValues[index + 2] = disp[2];
        rValues[index + 3] = rot[0];
        rValues[index + 4] = rot[1];
    }
    KRATOS_CATCH("");
}

BoundedVector<double, FiberBeamColumnElement3D2N::msElementSize,
FiberBeamColumnElement3D2N::msElementSize>
FiberBeamColumnElement3D2N::CalculateInitialLocalCS(){
    
    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    array_1d<double, msDimension> direction_vector_x = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_y = ZeroVector(msDimension);
    array_1d<double, msDimension> direction_vector_z = ZeroVector(msDimension);
    array_1d<double, msLocalSize> reference_coordinates = ZeroVector(msLocalSize);

    reference_coordinates[0] = this->GetGeometry()[0].X0();
    reference_coordinates[1] = this->GetGeometry()[0].Y0();
    reference_coordinates[2] = this->GetGeometry()[0].Z0();
    reference_coordinates[3] = this->GetGeometry()[1].X0();
    reference_coordinates[4] = this->GetGeometry()[1].Y0();
    reference_coordinates[5] = this->GetGeometry()[1].Z0();

    for (unsigned int i = 0; i < msDimension; ++i)
    {
        direction_vector_x[i] = (reference_coordinates[i + msDimension] - reference_coordinates[i]);
    }
    Matrix temp_matrix = ZeroMatrix(msDimension);

    // TODO
}

    /**
     * Initialize before element state determination
     * 
     */
    void InitializeBeforeElementStateDetermination()
    {
        mIterationCount = 0;
    }


    /**
     * Element state determination process consists of the following functions
     * 
     */ 

    /**
     * (4) Compute the element deformation increments
     * 
     * using the compatibility matrix L_ele, the element deformation increments ddq_i
     * is computed from the structure displacement increments ddp_i
     */
    void ComputeChangeInElementDeformationIncr()
    {
        KRATOS_TRY;
        // Structure displacement increments ddp_i calculated from the difference of displacement of 
        // current N-R iteration and previous N-R iteration
        Vector ddp_i = ZeroMatrix(10);
        Vector dp_previous = ZeroMatrix(10);;
        Vecotr dp_current = ZeroMatrix(10);;
        //ddp_i = this->getDeformationValuesVector(dp_current,0) - this->getDeformationValuesVector(dp_current,0);

        // construct the compatibility matrix L_ele
        Matrix L_ele = ZeroMatrix(5,10);
        mChangeInElementDeformationIncr = L_ele * ddp_i;
        // MathUtils<double>::dot_product(L_ele, ddp_i);

        KRATOS_CATCH("");
    }

    void ComputeChangeInElementDeformationIncr()
    {
        KRATOS_TRY;
        mElementDeformationIncr += 

        KRATOS_CATCH("");
    }



} // end: namespace Kratos