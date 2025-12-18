//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Gianmarco Boscolo
//

#include "fluid_topology_optimization_element.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "utilities/element_size_calculator.h"

// include Fluid Topology Optimization Data
#include "custom_utilities/fluid_topology_optimization_element_data.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TElementData >
FluidTopologyOptimizationElement<TElementData>::FluidTopologyOptimizationElement(IndexType NewId):
    Element(NewId)
{}

template< class TElementData >
FluidTopologyOptimizationElement<TElementData>::FluidTopologyOptimizationElement(IndexType NewId, const NodesArrayType& ThisNodes):
    Element(NewId,ThisNodes)
{}


template< class TElementData >
FluidTopologyOptimizationElement<TElementData>::FluidTopologyOptimizationElement(IndexType NewId, GeometryType::Pointer pGeometry):
    Element(NewId,pGeometry)
{}


template< class TElementData >
FluidTopologyOptimizationElement<TElementData>::FluidTopologyOptimizationElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Element(NewId,pGeometry,pProperties)
{}


template< class TElementData >
FluidTopologyOptimizationElement<TElementData>::~FluidTopologyOptimizationElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <class TElementData>
Element::Pointer FluidTopologyOptimizationElement<TElementData>::Create(
    IndexType NewId,
    NodesArrayType const &ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FluidTopologyOptimizationElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template <class TElementData>
Element::Pointer FluidTopologyOptimizationElement<TElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FluidTopologyOptimizationElement>(NewId, pGeom, pProperties);
}

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo) {
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();
        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);
        mpConstitutiveLaw->InitializeMaterial(r_properties,r_geometry,row(r_shape_functions,0));
    }

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {

                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedSystem(
                    data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) {

                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedSystemAdjoint(
                    data, rLeftHandSideMatrix, rRightHandSideVector);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
            }
        }
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                       const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedLHS(data, rLeftHandSideMatrix);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            // PRIMAL NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedLHSAdjoint(data, rLeftHandSideMatrix);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
        
            }
        }
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if constexpr (TElementData::ElementManagesTimeIntegration) 
    {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = data.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedRHS(data, rRightHandSideVector);
            }
        }
        else if (problem_physics == 2)
        {
            // ADJOINT NAVIER-STOKES
            // Iterate over integration points to evaluate local contribution
            for (unsigned int g = 0; g < number_of_gauss_points; g++) 
            {
                this->UpdateIntegrationPointData(
                    data, g, gauss_weights[g],
                    row(shape_functions, g),shape_derivatives[g]);

                this->AddTimeIntegratedRHSAdjoint(data, rRightHandSideVector);
            }
        }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
        
            }
        }
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateLocalVelocityContribution(
    MatrixType& rDampMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if( rDampMatrix.size1() != LocalSize )
        rDampMatrix.resize(LocalSize,LocalSize,false);

    if( rRightHandSideVector.size() != LocalSize )
        rRightHandSideVector.resize(LocalSize,false);

    noalias(rDampMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            const auto& r_dndx = shape_derivatives[g];
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),r_dndx);

            this->AddVelocitySystem(data, rDampMatrix, rRightHandSideVector);
        }
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix,
                                                     const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    noalias(rMassMatrix) = ZeroMatrix(LocalSize, LocalSize);

    if (!TElementData::ElementManagesTimeIntegration) {
        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(
            gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Iterate over integration points to evaluate local contribution
        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            this->UpdateIntegrationPointData(
                data, g, gauss_weights[g],
                row(shape_functions, g),shape_derivatives[g]);

            this->AddMassLHS(data, rMassMatrix);
        }
    }
}

template< class TElementData >
void FluidTopologyOptimizationElement< TElementData >::EquationIdVector(EquationIdVectorType &rResult, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
    unsigned int problem_physics = rCurrentProcessInfo[FLUID_TOP_OPT_PROBLEM_STAGE];
    // problem_physics = 1: NS equations
    // problem_physics = 2: ADJOINT NS equations
    if (problem_physics == 1)
    {
        // PRIMAL NAVIER-STOKES
        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
        const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_X,xpos).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
            if (Dim == 3) rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(PRESSURE,ppos).EquationId();
        }
    }
    else if (problem_physics == 2)
    {
        // ADJOINT NAVIER-STOKES
        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_ADJ_X);
        const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE_ADJ);

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_ADJ_X,xpos).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_ADJ_Y,xpos+1).EquationId();
            if (Dim == 3) rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_ADJ_Z,xpos+2).EquationId();
            rResult[LocalIndex++] = r_geometry[i].GetDof(PRESSURE_ADJ,ppos).EquationId();
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
    
        }
    }
}


template< class TElementData >
void FluidTopologyOptimizationElement< TElementData >::GetDofList(DofsVectorType &rElementalDofList, const ProcessInfo &rCurrentProcessInfo) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rElementalDofList.size() != LocalSize)
         rElementalDofList.resize(LocalSize);

    // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
    unsigned int problem_physics = rCurrentProcessInfo[FLUID_TOP_OPT_PROBLEM_STAGE];
    // problem_physics = 1: NS equations
    // problem_physics = 2: ADJOINT NS equations
    if (problem_physics == 1)
    {
        // PRIMAL NAVIER-STOKES
        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
        const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_X,xpos);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Y,xpos+1);
            if (Dim == 3) rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Z,xpos+2);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PRESSURE,ppos);
        }
    }
    else if (problem_physics == 2)
    {
        // ADJOINT NAVIER-STOKES
        const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_ADJ_X);
        const unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE_ADJ);

        unsigned int LocalIndex = 0;
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_ADJ_X,xpos);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_ADJ_Y,xpos+1);
            if (Dim == 3) rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_ADJ_Z,xpos+2);
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PRESSURE_ADJ,ppos);
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
    
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::GetFirstDerivativesVector(Vector &rValues, int Step) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    unsigned int Index = 0;

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        const array_1d<double,3>& rVel = r_geometry[i].FastGetSolutionStepValue(VELOCITY,Step);
        for (unsigned int d = 0; d < Dim; d++)
            rValues[Index++] = rVel[d];
        rValues[Index++] = r_geometry[i].FastGetSolutionStepValue(PRESSURE,Step);
    }
}


template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::GetSecondDerivativesVector(Vector &rValues, int Step) const
{
    const GeometryType& r_geometry = this->GetGeometry();

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);

    unsigned int index = 0;

    for (unsigned int i = 0; i < NumNodes; i++) {
        const array_1d<double,3>& r_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION,Step);
        for (unsigned int d = 0; d < Dim; d++)
            rValues[index++] = r_acceleration[d];
        rValues[index++] = 0.0; // skip pressure Dof
    }
}


template< class TElementData >
GeometryData::IntegrationMethod FluidTopologyOptimizationElement<TElementData>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int FluidTopologyOptimizationElement<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0) {
        return out;
    }

    // Check variables used by TElementData
    out = TElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    const GeometryType& r_geometry = this->GetGeometry();

    unsigned int problem_physics = rCurrentProcessInfo[FLUID_TOP_OPT_PROBLEM_STAGE];
    if (problem_physics == 1)
    {
        for(unsigned int i=0; i<NumNodes; ++i)
        {
            const Node& rNode = r_geometry[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
            
            // Check that required dofs exist for PRIMAL NS
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X,rNode);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y,rNode);
            if (Dim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z,rNode);
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
        }
    }
    else if (problem_physics == 2)
    {
        for(unsigned int i=0; i<NumNodes; ++i)
        {
            const Node& rNode = r_geometry[i];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION_ADJ,rNode);
            
            // Check that required dofs exist for ADJOINT NS
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_ADJ_X,rNode);
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_ADJ_Y,rNode);
            if (Dim == 3) KRATOS_CHECK_DOF_IN_NODE(VELOCITY_ADJ_Z,rNode);
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE_ADJ,rNode);
        }
    }
    else
    {
        if (problem_physics != 0)
        {
            KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
    
        }    
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if ( Dim == 2)
    {
        for (unsigned int i=0; i<NumNodes; ++i) {
            if (r_geometry[i].Z() != 0.0)
                KRATOS_ERROR << "Node " << r_geometry[i].Id() << "has non-zero Z coordinate." << std::endl;
        }
    }

    // Check the constitutive law
    KRATOS_ERROR_IF(mpConstitutiveLaw == nullptr) << "Constitutive Law not initialized for Element " << this->Info() << std::endl;

    constexpr auto dimension = Dim;  // I need to set this here otherwise it gives me a linking error when attempting to '<<' Dim.

    KRATOS_ERROR_IF(mpConstitutiveLaw->WorkingSpaceDimension() != Dim)
        << "Wrong dimension: The " << mpConstitutiveLaw->WorkingSpaceDimension()
        << "D constitutive law " << mpConstitutiveLaw->Info()
        << " is not compatible with " << dimension << "D element " << this->Info()
        << "." << std::endl;

    out = mpConstitutiveLaw->Check(this->GetProperties(),r_geometry,rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT( out == 0) << "The Constitutive Law provided for Element " << this->Info() << " is not correct." << std::endl;

    return out;
}


///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 3 > > const& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    // if (rVariable == VORTICITY)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateVorticityVector(this->GetGeometry(),shape_function_gradients,rValues);
    // }
}


template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<double> const& rVariable,
    std::vector<double>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{
    // if (rVariable == Q_VALUE)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateQValue(this->GetGeometry(),shape_function_gradients,rValues);
    // }
    // else if (rVariable == VORTICITY_MAGNITUDE)
    // {
    //     // Get Shape function data
    //     Vector gauss_weights;
    //     Matrix shape_functions;
    //     ShapeFunctionDerivativesArrayType shape_function_gradients;
    //     this->CalculateGeometryData(gauss_weights,shape_functions,shape_function_gradients);

    //     VorticityUtilities<Dim>::CalculateVorticityMagnitude(this->GetGeometry(),shape_function_gradients,rValues);
    // }
    // else if (rVariable == UPDATE_STATISTICS)
    // {
    //     KRATOS_DEBUG_ERROR_IF(!rCurrentProcessInfo.Has(STATISTICS_CONTAINER)) << "Trying to compute turbulent statistics, but ProcessInfo does not have STATISTICS_CONTAINER defined." << std::endl;
    //     rCurrentProcessInfo.GetValue(STATISTICS_CONTAINER)->UpdateStatistics(this);
    // }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<array_1d<double, 6>> const& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Vector> const& rVariable,
    std::vector<Vector>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateOnIntegrationPoints(
    Variable<Matrix> const& rVariable,
    std::vector<Matrix>& rValues,
    ProcessInfo const& rCurrentProcessInfo)
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TElementData>
std::string FluidTopologyOptimizationElement<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidTopologyOptimizationElement" << Dim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (this->GetConstitutiveLaw() != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        this->GetConstitutiveLaw()->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <class TElementData>
double FluidTopologyOptimizationElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalScalarData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    double result = 0.0;

    for (size_t i = 0; i < NumNodes; i++) {
        result += rN[i] * rValues[i];
    }

    return result;
}

template <class TElementData>
array_1d<double, 3> FluidTopologyOptimizationElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalVectorData& rValues,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    array_1d<double, 3> result = ZeroVector(3);

    for (size_t i = 0; i < NumNodes; i++) {
        for (size_t j = 0; j < Dim; j++) {
            result[j] += rN[i] * rValues(i, j);
        }
    }

    return result;
}

template <class TElementData>
BoundedMatrix<double, TElementData::Dim, TElementData::Dim> FluidTopologyOptimizationElement<TElementData>::GetAtCoordinate(
    const typename TElementData::NodalTensorData &rValues,
    const typename TElementData::ShapeFunctionsType &rN) const
{
    BoundedMatrix<double,Dim,Dim> result = ZeroMatrix(Dim,Dim);

    for (size_t i = 0; i < NumNodes; i++) {
        noalias(result) += rN[i] * rValues[i];
    }

    return result;
}

template <class TElementData>
double FluidTopologyOptimizationElement<TElementData>::GetAtCoordinate(
    const double Value,
    const typename TElementData::ShapeFunctionsType& rN) const
{
    return Value;
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::UpdateIntegrationPointData(
    TElementData& rData,
    unsigned int IntegrationPointIndex,
    double Weight,
    const typename TElementData::MatrixRowType& rN,
    const typename TElementData::ShapeDerivativesType& rDN_DX) const
{
    rData.UpdateGeometryValues(IntegrationPointIndex, Weight, rN, rDN_DX);

    this->CalculateMaterialResponse(rData);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateMaterialResponse(TElementData& rData) const
{
    this->CalculateStrainRate(rData);

    auto& Values = rData.ConstitutiveLawValues;

    const Vector& shape_functions_vector = rData.N;
    const Matrix& shape_functions_derivative_matrix = rData.DN_DX;
    Values.SetShapeFunctionsValues(shape_functions_vector);
    Values.SetShapeFunctionsDerivatives(shape_functions_derivative_matrix);

    //ATTENTION: here we assume that only one constitutive law is employed for all of the gauss points in the element.
    //this is ok under the hypothesis that no history dependent behavior is employed
    mpConstitutiveLaw->CalculateMaterialResponseCauchy(Values);

    mpConstitutiveLaw->CalculateValue(Values,EFFECTIVE_VISCOSITY,rData.EffectiveViscosity);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::CalculateStrainRate(TElementData& rData) const
{
    // Now the distinction of what to do for the primal or adjoint problem must be done calling the relative Kratos property
        unsigned int problem_physics = rData.TopOptProblemStage;
        // problem_physics = 1: NS equations
        // problem_physics = 2: ADJOINT NS equations
        if (problem_physics == 1)
        {
            // PRIMAL NAVIER-STOKES
            Internals::StrainRateSpecialization<TElementData,Dim>::Calculate(
            rData.StrainRate,
            rData.Velocity,
            rData.DN_DX);
        }
        else if (problem_physics == 2)
        {
            // ADJOINT NAVIER-STOKES
            Internals::StrainRateSpecialization<TElementData,Dim>::Calculate(
            rData.StrainRate,
            rData.Velocity_adj,
            rData.DN_DX);
            }
        else
        {
            if (problem_physics != 0)
            {
                KRATOS_ERROR << "\nInvalid value for the variable FLUID_TOP_OPT_PROBLEM_STAGE for the Fluid Topology Optimization aplication. |\t TopOptProblemStage = " << problem_physics << " |\t Accepted values are 1->NS, 2->ADJ_NS.\n";
            }
        }
    
}

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::CalculateGeometryData(Vector &rGaussWeights,
                                      Matrix &rNContainer,
                                      ShapeFunctionDerivativesArrayType &rDN_DX) const
{
    const GeometryData::IntegrationMethod integration_method = this->GetIntegrationMethod();
    const GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_gauss_points = r_geometry.IntegrationPointsNumber(integration_method);

    Vector DetJ;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,integration_method);

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != NumNodes) {
        rNContainer.resize(number_of_gauss_points,NumNodes,false);
    }
    rNContainer = r_geometry.ShapeFunctionsValues(integration_method);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = r_geometry.IntegrationPoints(integration_method);

    if (rGaussWeights.size() != number_of_gauss_points) {
        rGaussWeights.resize(number_of_gauss_points,false);
    }

    for (unsigned int g = 0; g < number_of_gauss_points; g++)
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<double> &rVariable,
    double &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<Vector >& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo )
{
    noalias( rOutput ) = ZeroVector( StrainSize );

    if (rVariable == FLUID_STRESS) {

        // creating a new data container that goes out of scope after the function is left
        TElementData data_local;

        // transferring the velocity (among other variables)
        data_local.Initialize(*this, rCurrentProcessInfo);

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;

        // computing DN_DX values for the strain rate
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int number_of_gauss_points = gauss_weights.size();

        double sum_of_gauss_weights = 0.0;

        for (unsigned int g = 0; g < number_of_gauss_points; g++){

            this->UpdateIntegrationPointData(data_local, g, gauss_weights[g], row(shape_functions, g), shape_derivatives[g]);

            const Vector gauss_point_contribution = data_local.ShearStress;

            noalias( rOutput ) += gauss_point_contribution * gauss_weights[g];
            sum_of_gauss_weights += gauss_weights[g];
        }

        for (unsigned int i = 0; i < StrainSize; i++){
            rOutput[i] = ( 1.0 / sum_of_gauss_weights ) * rOutput[i];
        }

    } else {
        Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    Element::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::ConvectionOperator(Vector &rResult,
                                   const array_1d<double,3> &rConvVel,
                                   const ShapeFunctionDerivativesType &DN_DX) const
{
    if(rResult.size() != NumNodes) rResult.resize(NumNodes,false);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        rResult[i] = rConvVel[0]*DN_DX(i,0);
        for(unsigned int k = 1; k < Dim; k++)
            rResult[i] += rConvVel[k]*DN_DX(i,k);
    }
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedSystem(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedLHS(
    TElementData& rData, MatrixType& rLHS) 
{
    this->ComputeGaussPointLHSContribution(rData, rLHS);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedRHS(
    TElementData& rData, VectorType& rRHS) 
{
    this->ComputeGaussPointRHSContribution(rData, rRHS);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedSystemAdjoint(
    TElementData& rData, MatrixType& rLHS, VectorType& rRHS) 
{
    this->ComputeGaussPointLHSContributionAdjoint(rData, rLHS);
    this->ComputeGaussPointRHSContributionAdjoint(rData, rRHS);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedLHSAdjoint(
    TElementData& rData, MatrixType& rLHS) 
{
    this->ComputeGaussPointLHSContributionAdjoint(rData, rLHS);
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddTimeIntegratedRHSAdjoint(
    TElementData& rData, VectorType& rRHS) 
{
    this->ComputeGaussPointRHSContributionAdjoint(rData, rRHS);
}

template <>
void FluidTopologyOptimizationElement< FluidTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointLHSContribution(
    FluidTopologyOptimizationElementData<2,3,true> & rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,3>& c = rData.SoundVelocity;

    const array_1d<double,3> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;
    
    // Get constitutive matrix
    const BoundedMatrix<double,3,3>& C = rData.C;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> vconv = rData.Velocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = pow(DN(0,0), 2);
const double crLHS4 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crLHS5 = crLHS4*stab_c3;
const double crLHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS8 = rho*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2));
const double crLHS9 = h*(crLHS5*h + crLHS8)/stab_c1 + mu;
const double crLHS10 = pow(N[0], 2);
const double crLHS11 = DN(0,0)*crLHS6;
const double crLHS12 = DN(0,1)*crLHS7;
const double crLHS13 = crLHS11 + crLHS12;
const double crLHS14 = N[0]*rho;
const double crLHS15 = bdf0*rho;
const double crLHS16 = bdf0*crLHS14;
const double crLHS17 = N[0]*crLHS4;
const double crLHS18 = crLHS11*rho;
const double crLHS19 = crLHS12*rho;
const double crLHS20 = crLHS16 + crLHS17 + crLHS18 + crLHS19;
const double crLHS21 = 1.0/(crLHS5 + crLHS8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crLHS22 = 1.0*crLHS21;
const double crLHS23 = crLHS22*rho;
const double crLHS24 = crLHS13*crLHS23;
const double crLHS25 = 1.0*crLHS17;
const double crLHS26 = crLHS21*crLHS25;
const double crLHS27 = crLHS20*crLHS22;
const double crLHS28 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS29 = crLHS14*crLHS28;
const double crLHS30 = crLHS10*crLHS15 + crLHS10*crLHS4 + crLHS13*crLHS14 + crLHS20*crLHS24 - crLHS20*crLHS26 + crLHS27*crLHS29;
const double crLHS31 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS32 = C(1,2)*DN(0,1);
const double crLHS33 = C(2,2)*DN(0,0) + crLHS32;
const double crLHS34 = DN(0,0)*crLHS9;
const double crLHS35 = DN(0,1)*crLHS34;
const double crLHS36 = crLHS22*crLHS28;
const double crLHS37 = gauss_weight*(-N[0] + crLHS14*crLHS36 + crLHS24 - crLHS26);
const double crLHS38 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS39 = C(0,2)*DN(1,0);
const double crLHS40 = C(2,2)*DN(1,1) + crLHS39;
const double crLHS41 = DN(0,0)*DN(1,0);
const double crLHS42 = N[1]*crLHS16 + N[1]*crLHS17;
const double crLHS43 = crLHS41*crLHS9 + crLHS42;
const double crLHS44 = DN(1,0)*crLHS6;
const double crLHS45 = DN(1,1)*crLHS7;
const double crLHS46 = crLHS44 + crLHS45;
const double crLHS47 = N[1]*crLHS15;
const double crLHS48 = N[1]*crLHS4;
const double crLHS49 = crLHS44*rho;
const double crLHS50 = crLHS45*rho;
const double crLHS51 = crLHS47 + crLHS48 + crLHS49 + crLHS50;
const double crLHS52 = crLHS22*crLHS51;
const double crLHS53 = crLHS14*crLHS46 + crLHS24*crLHS51 - crLHS26*crLHS51 + crLHS29*crLHS52;
const double crLHS54 = C(0,1)*DN(1,1) + crLHS39;
const double crLHS55 = C(1,2)*DN(1,1);
const double crLHS56 = C(2,2)*DN(1,0) + crLHS55;
const double crLHS57 = DN(1,1)*crLHS34;
const double crLHS58 = DN(0,0)*N[1];
const double crLHS59 = DN(1,0)*N[0];
const double crLHS60 = crLHS23*crLHS28;
const double crLHS61 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS62 = C(0,2)*DN(2,0);
const double crLHS63 = C(2,2)*DN(2,1) + crLHS62;
const double crLHS64 = DN(0,0)*DN(2,0);
const double crLHS65 = N[2]*crLHS16 + N[2]*crLHS17;
const double crLHS66 = crLHS64*crLHS9 + crLHS65;
const double crLHS67 = DN(2,0)*crLHS6;
const double crLHS68 = DN(2,1)*crLHS7;
const double crLHS69 = crLHS67 + crLHS68;
const double crLHS70 = N[2]*crLHS15;
const double crLHS71 = N[2]*crLHS4;
const double crLHS72 = crLHS67*rho;
const double crLHS73 = crLHS68*rho;
const double crLHS74 = crLHS70 + crLHS71 + crLHS72 + crLHS73;
const double crLHS75 = crLHS22*crLHS74;
const double crLHS76 = crLHS14*crLHS69 + crLHS24*crLHS74 - crLHS26*crLHS74 + crLHS29*crLHS75;
const double crLHS77 = C(0,1)*DN(2,1) + crLHS62;
const double crLHS78 = C(1,2)*DN(2,1);
const double crLHS79 = C(2,2)*DN(2,0) + crLHS78;
const double crLHS80 = DN(2,1)*crLHS34;
const double crLHS81 = DN(0,0)*N[2];
const double crLHS82 = DN(2,0)*N[0];
const double crLHS83 = C(0,1)*DN(0,0) + crLHS32;
const double crLHS84 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS85 = pow(DN(0,1), 2);
const double crLHS86 = C(0,1)*DN(1,0) + crLHS55;
const double crLHS87 = DN(0,1)*crLHS9;
const double crLHS88 = DN(1,0)*crLHS87;
const double crLHS89 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS90 = DN(0,1)*DN(1,1);
const double crLHS91 = crLHS42 + crLHS9*crLHS90;
const double crLHS92 = DN(0,1)*N[1];
const double crLHS93 = DN(1,1)*N[0];
const double crLHS94 = C(0,1)*DN(2,0) + crLHS78;
const double crLHS95 = DN(2,0)*crLHS87;
const double crLHS96 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS97 = DN(0,1)*DN(2,1);
const double crLHS98 = crLHS65 + crLHS9*crLHS97;
const double crLHS99 = DN(0,1)*N[2];
const double crLHS100 = DN(2,1)*N[0];
const double crLHS101 = gauss_weight*(N[0] + crLHS21*(1.0*crLHS16 + 1.0*crLHS18 + 1.0*crLHS19 + crLHS25));
const double crLHS102 = crLHS22*gauss_weight;
const double crLHS103 = crLHS102*(crLHS41 + crLHS90);
const double crLHS104 = crLHS102*(crLHS64 + crLHS97);
const double crLHS105 = N[1]*rho;
const double crLHS106 = crLHS23*crLHS46;
const double crLHS107 = 1.0*crLHS48;
const double crLHS108 = crLHS107*crLHS21;
const double crLHS109 = crLHS105*crLHS28;
const double crLHS110 = crLHS105*crLHS13 + crLHS106*crLHS20 - crLHS108*crLHS20 + crLHS109*crLHS27;
const double crLHS111 = pow(DN(1,0), 2);
const double crLHS112 = pow(N[1], 2);
const double crLHS113 = crLHS105*crLHS46 + crLHS106*crLHS51 - crLHS108*crLHS51 + crLHS109*crLHS52 + crLHS112*crLHS15 + crLHS112*crLHS4;
const double crLHS114 = DN(1,0)*crLHS9;
const double crLHS115 = DN(1,1)*crLHS114;
const double crLHS116 = gauss_weight*(-N[1] + crLHS105*crLHS36 + crLHS106 - crLHS108);
const double crLHS117 = DN(1,0)*DN(2,0);
const double crLHS118 = N[2]*crLHS47 + N[2]*crLHS48;
const double crLHS119 = crLHS117*crLHS9 + crLHS118;
const double crLHS120 = crLHS105*crLHS69 + crLHS106*crLHS74 - crLHS108*crLHS74 + crLHS109*crLHS75;
const double crLHS121 = DN(2,1)*crLHS114;
const double crLHS122 = DN(1,0)*N[2];
const double crLHS123 = DN(2,0)*N[1];
const double crLHS124 = pow(DN(1,1), 2);
const double crLHS125 = DN(2,0)*crLHS9;
const double crLHS126 = DN(1,1)*crLHS125;
const double crLHS127 = DN(1,1)*DN(2,1);
const double crLHS128 = crLHS118 + crLHS127*crLHS9;
const double crLHS129 = DN(1,1)*N[2];
const double crLHS130 = DN(2,1)*N[1];
const double crLHS131 = gauss_weight*(N[1] + crLHS21*(crLHS107 + 1.0*crLHS47 + 1.0*crLHS49 + 1.0*crLHS50));
const double crLHS132 = crLHS102*(crLHS117 + crLHS127);
const double crLHS133 = N[2]*rho;
const double crLHS134 = crLHS23*crLHS69;
const double crLHS135 = 1.0*crLHS71;
const double crLHS136 = crLHS135*crLHS21;
const double crLHS137 = crLHS133*crLHS28;
const double crLHS138 = crLHS13*crLHS133 + crLHS134*crLHS20 - crLHS136*crLHS20 + crLHS137*crLHS27;
const double crLHS139 = crLHS133*crLHS46 + crLHS134*crLHS51 - crLHS136*crLHS51 + crLHS137*crLHS52;
const double crLHS140 = pow(DN(2,0), 2);
const double crLHS141 = pow(N[2], 2);
const double crLHS142 = crLHS133*crLHS69 + crLHS134*crLHS74 - crLHS136*crLHS74 + crLHS137*crLHS75 + crLHS141*crLHS15 + crLHS141*crLHS4;
const double crLHS143 = DN(2,1)*crLHS125;
const double crLHS144 = gauss_weight*(-N[2] + crLHS133*crLHS36 + crLHS134 - crLHS136);
const double crLHS145 = pow(DN(2,1), 2);
const double crLHS146 = gauss_weight*(N[2] + crLHS21*(crLHS135 + 1.0*crLHS70 + 1.0*crLHS72 + 1.0*crLHS73));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS3*crLHS9 + crLHS30);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS31 + DN(0,1)*crLHS33 + crLHS35);
rLHS(0,2)+=DN(0,0)*crLHS37;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS40 + crLHS43 + crLHS53);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS56 + crLHS57);
rLHS(0,5)+=-gauss_weight*(-DN(1,0)*crLHS24 + DN(1,0)*crLHS26 + crLHS58 - crLHS59*crLHS60);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + crLHS66 + crLHS76);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS79 + crLHS80);
rLHS(0,8)+=-gauss_weight*(-DN(2,0)*crLHS24 + DN(2,0)*crLHS26 - crLHS60*crLHS82 + crLHS81);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS83 + crLHS35);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS33 + DN(0,1)*crLHS84 + crLHS30 + crLHS85*crLHS9);
rLHS(1,2)+=DN(0,1)*crLHS37;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS86 + crLHS88);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS56 + DN(0,1)*crLHS89 + crLHS53 + crLHS91);
rLHS(1,5)+=-gauss_weight*(-DN(1,1)*crLHS24 + DN(1,1)*crLHS26 - crLHS60*crLHS93 + crLHS92);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS94 + crLHS95);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS96 + crLHS76 + crLHS98);
rLHS(1,8)+=-gauss_weight*(-DN(2,1)*crLHS24 + DN(2,1)*crLHS26 - crLHS100*crLHS60 + crLHS99);
rLHS(2,0)+=DN(0,0)*crLHS101;
rLHS(2,1)+=DN(0,1)*crLHS101;
rLHS(2,2)+=crLHS102*(crLHS3 + crLHS85);
rLHS(2,3)+=gauss_weight*(DN(0,0)*crLHS52 + crLHS59);
rLHS(2,4)+=gauss_weight*(DN(0,1)*crLHS52 + crLHS93);
rLHS(2,5)+=crLHS103;
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS75 + crLHS82);
rLHS(2,7)+=gauss_weight*(DN(0,1)*crLHS75 + crLHS100);
rLHS(2,8)+=crLHS104;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS110 + crLHS43);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS31 + DN(1,1)*crLHS33 + crLHS88);
rLHS(3,2)+=gauss_weight*(DN(0,0)*crLHS106 - DN(0,0)*crLHS108 + crLHS58*crLHS60 - crLHS59);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS40 + crLHS111*crLHS9 + crLHS113);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS56 + crLHS115);
rLHS(3,5)+=DN(1,0)*crLHS116;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + crLHS119 + crLHS120);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS79 + crLHS121);
rLHS(3,8)+=-gauss_weight*(-DN(2,0)*crLHS106 + DN(2,0)*crLHS108 + crLHS122 - crLHS123*crLHS60);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS83 + crLHS57);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS33 + DN(1,1)*crLHS84 + crLHS110 + crLHS91);
rLHS(4,2)+=gauss_weight*(DN(0,1)*crLHS106 - DN(0,1)*crLHS108 + crLHS60*crLHS92 - crLHS93);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS86 + crLHS115);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS56 + DN(1,1)*crLHS89 + crLHS113 + crLHS124*crLHS9);
rLHS(4,5)+=DN(1,1)*crLHS116;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS94 + crLHS126);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS96 + crLHS120 + crLHS128);
rLHS(4,8)+=-gauss_weight*(-DN(2,1)*crLHS106 + DN(2,1)*crLHS108 + crLHS129 - crLHS130*crLHS60);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS27 + crLHS58);
rLHS(5,1)+=gauss_weight*(DN(1,1)*crLHS27 + crLHS92);
rLHS(5,2)+=crLHS103;
rLHS(5,3)+=DN(1,0)*crLHS131;
rLHS(5,4)+=DN(1,1)*crLHS131;
rLHS(5,5)+=crLHS102*(crLHS111 + crLHS124);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS75 + crLHS123);
rLHS(5,7)+=gauss_weight*(DN(1,1)*crLHS75 + crLHS130);
rLHS(5,8)+=crLHS132;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS138 + crLHS66);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS31 + DN(2,1)*crLHS33 + crLHS95);
rLHS(6,2)+=gauss_weight*(DN(0,0)*crLHS134 - DN(0,0)*crLHS136 + crLHS60*crLHS81 - crLHS82);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS40 + crLHS119 + crLHS139);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS56 + crLHS126);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS134 - DN(1,0)*crLHS136 + crLHS122*crLHS60 - crLHS123);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + crLHS140*crLHS9 + crLHS142);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS79 + crLHS143);
rLHS(6,8)+=DN(2,0)*crLHS144;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS83 + crLHS80);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS33 + DN(2,1)*crLHS84 + crLHS138 + crLHS98);
rLHS(7,2)+=gauss_weight*(DN(0,1)*crLHS134 - DN(0,1)*crLHS136 - crLHS100 + crLHS60*crLHS99);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS86 + crLHS121);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS56 + DN(2,1)*crLHS89 + crLHS128 + crLHS139);
rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS134 - DN(1,1)*crLHS136 + crLHS129*crLHS60 - crLHS130);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS94 + crLHS143);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS96 + crLHS142 + crLHS145*crLHS9);
rLHS(7,8)+=DN(2,1)*crLHS144;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS27 + crLHS81);
rLHS(8,1)+=gauss_weight*(DN(2,1)*crLHS27 + crLHS99);
rLHS(8,2)+=crLHS104;
rLHS(8,3)+=gauss_weight*(DN(2,0)*crLHS52 + crLHS122);
rLHS(8,4)+=gauss_weight*(DN(2,1)*crLHS52 + crLHS129);
rLHS(8,5)+=crLHS132;
rLHS(8,6)+=DN(2,0)*crLHS146;
rLHS(8,7)+=DN(2,1)*crLHS146;
rLHS(8,8)+=crLHS102*(crLHS140 + crLHS145);
    
}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContribution(
    FluidTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,4>& c = rData.SoundVelocity;

    const array_1d<double,4> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;
    
    // Get constitutive matrix
    const BoundedMatrix<double,6,6>& C = rData.C;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> vconv = rData.Velocity - rData.MeshVelocity;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS1 = C(0,3)*DN(0,0);
const double crLHS2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS1;
const double crLHS3 = C(0,5)*DN(0,0);
const double crLHS4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS3;
const double crLHS5 = pow(DN(0,0), 2);
const double crLHS6 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crLHS7 = crLHS6*stab_c3;
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS11 = rho*stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2));
const double crLHS12 = h*(crLHS11 + crLHS7*h)/stab_c1 + mu;
const double crLHS13 = pow(N[0], 2);
const double crLHS14 = DN(0,0)*crLHS8;
const double crLHS15 = DN(0,1)*crLHS9;
const double crLHS16 = DN(0,2)*crLHS10;
const double crLHS17 = crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = N[0]*rho;
const double crLHS19 = bdf0*rho;
const double crLHS20 = bdf0*crLHS18;
const double crLHS21 = N[0]*crLHS6;
const double crLHS22 = crLHS14*rho;
const double crLHS23 = crLHS15*rho;
const double crLHS24 = crLHS16*rho;
const double crLHS25 = crLHS20 + crLHS21 + crLHS22 + crLHS23 + crLHS24;
const double crLHS26 = 1.0/(crLHS11/h + crLHS7 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crLHS27 = 1.0*crLHS26;
const double crLHS28 = crLHS27*rho;
const double crLHS29 = crLHS17*crLHS28;
const double crLHS30 = 1.0*crLHS21;
const double crLHS31 = crLHS26*crLHS30;
const double crLHS32 = crLHS25*crLHS27;
const double crLHS33 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS34 = crLHS18*crLHS33;
const double crLHS35 = crLHS13*crLHS19 + crLHS13*crLHS6 + crLHS17*crLHS18 + crLHS25*crLHS29 - crLHS25*crLHS31 + crLHS32*crLHS34;
const double crLHS36 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS37 = C(1,3)*DN(0,1);
const double crLHS38 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS37;
const double crLHS39 = C(3,5)*DN(0,0);
const double crLHS40 = C(4,5)*DN(0,2);
const double crLHS41 = C(1,5)*DN(0,1) + crLHS39 + crLHS40;
const double crLHS42 = DN(0,0)*crLHS12;
const double crLHS43 = DN(0,1)*crLHS42;
const double crLHS44 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS45 = C(3,4)*DN(0,1);
const double crLHS46 = C(2,3)*DN(0,2) + crLHS39 + crLHS45;
const double crLHS47 = C(2,5)*DN(0,2);
const double crLHS48 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS47;
const double crLHS49 = DN(0,2)*crLHS42;
const double crLHS50 = crLHS27*crLHS33;
const double crLHS51 = gauss_weight*(-N[0] + crLHS18*crLHS50 + crLHS29 - crLHS31);
const double crLHS52 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS53 = C(0,3)*DN(1,0);
const double crLHS54 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS53;
const double crLHS55 = C(0,5)*DN(1,0);
const double crLHS56 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS55;
const double crLHS57 = DN(0,0)*DN(1,0);
const double crLHS58 = N[1]*crLHS20 + N[1]*crLHS21;
const double crLHS59 = crLHS12*crLHS57 + crLHS58;
const double crLHS60 = DN(1,0)*crLHS8;
const double crLHS61 = DN(1,1)*crLHS9;
const double crLHS62 = DN(1,2)*crLHS10;
const double crLHS63 = crLHS60 + crLHS61 + crLHS62;
const double crLHS64 = N[1]*crLHS19;
const double crLHS65 = N[1]*crLHS6;
const double crLHS66 = crLHS60*rho;
const double crLHS67 = crLHS61*rho;
const double crLHS68 = crLHS62*rho;
const double crLHS69 = crLHS64 + crLHS65 + crLHS66 + crLHS67 + crLHS68;
const double crLHS70 = crLHS27*crLHS69;
const double crLHS71 = crLHS18*crLHS63 + crLHS29*crLHS69 - crLHS31*crLHS69 + crLHS34*crLHS70;
const double crLHS72 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS53;
const double crLHS73 = C(1,3)*DN(1,1);
const double crLHS74 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS73;
const double crLHS75 = C(3,5)*DN(1,0);
const double crLHS76 = C(4,5)*DN(1,2);
const double crLHS77 = C(1,5)*DN(1,1) + crLHS75 + crLHS76;
const double crLHS78 = DN(1,1)*crLHS42;
const double crLHS79 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS55;
const double crLHS80 = C(3,4)*DN(1,1);
const double crLHS81 = C(2,3)*DN(1,2) + crLHS75 + crLHS80;
const double crLHS82 = C(2,5)*DN(1,2);
const double crLHS83 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS82;
const double crLHS84 = DN(1,2)*crLHS42;
const double crLHS85 = DN(0,0)*N[1];
const double crLHS86 = DN(1,0)*N[0];
const double crLHS87 = crLHS28*crLHS33;
const double crLHS88 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS89 = C(0,3)*DN(2,0);
const double crLHS90 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS89;
const double crLHS91 = C(0,5)*DN(2,0);
const double crLHS92 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS91;
const double crLHS93 = DN(0,0)*DN(2,0);
const double crLHS94 = N[2]*crLHS20 + N[2]*crLHS21;
const double crLHS95 = crLHS12*crLHS93 + crLHS94;
const double crLHS96 = DN(2,0)*crLHS8;
const double crLHS97 = DN(2,1)*crLHS9;
const double crLHS98 = DN(2,2)*crLHS10;
const double crLHS99 = crLHS96 + crLHS97 + crLHS98;
const double crLHS100 = N[2]*crLHS19;
const double crLHS101 = N[2]*crLHS6;
const double crLHS102 = crLHS96*rho;
const double crLHS103 = crLHS97*rho;
const double crLHS104 = crLHS98*rho;
const double crLHS105 = crLHS100 + crLHS101 + crLHS102 + crLHS103 + crLHS104;
const double crLHS106 = crLHS105*crLHS27;
const double crLHS107 = crLHS105*crLHS29 - crLHS105*crLHS31 + crLHS106*crLHS34 + crLHS18*crLHS99;
const double crLHS108 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS89;
const double crLHS109 = C(1,3)*DN(2,1);
const double crLHS110 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS109;
const double crLHS111 = C(3,5)*DN(2,0);
const double crLHS112 = C(4,5)*DN(2,2);
const double crLHS113 = C(1,5)*DN(2,1) + crLHS111 + crLHS112;
const double crLHS114 = DN(2,1)*crLHS42;
const double crLHS115 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS91;
const double crLHS116 = C(3,4)*DN(2,1);
const double crLHS117 = C(2,3)*DN(2,2) + crLHS111 + crLHS116;
const double crLHS118 = C(2,5)*DN(2,2);
const double crLHS119 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS118;
const double crLHS120 = DN(2,2)*crLHS42;
const double crLHS121 = DN(0,0)*N[2];
const double crLHS122 = DN(2,0)*N[0];
const double crLHS123 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS124 = C(0,3)*DN(3,0);
const double crLHS125 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS124;
const double crLHS126 = C(0,5)*DN(3,0);
const double crLHS127 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS126;
const double crLHS128 = DN(0,0)*DN(3,0);
const double crLHS129 = N[3]*crLHS20 + N[3]*crLHS21;
const double crLHS130 = crLHS12*crLHS128 + crLHS129;
const double crLHS131 = DN(3,0)*crLHS8;
const double crLHS132 = DN(3,1)*crLHS9;
const double crLHS133 = DN(3,2)*crLHS10;
const double crLHS134 = crLHS131 + crLHS132 + crLHS133;
const double crLHS135 = N[3]*crLHS19;
const double crLHS136 = N[3]*crLHS6;
const double crLHS137 = crLHS131*rho;
const double crLHS138 = crLHS132*rho;
const double crLHS139 = crLHS133*rho;
const double crLHS140 = crLHS135 + crLHS136 + crLHS137 + crLHS138 + crLHS139;
const double crLHS141 = crLHS140*crLHS27;
const double crLHS142 = crLHS134*crLHS18 + crLHS140*crLHS29 - crLHS140*crLHS31 + crLHS141*crLHS34;
const double crLHS143 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS124;
const double crLHS144 = C(1,3)*DN(3,1);
const double crLHS145 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS144;
const double crLHS146 = C(3,5)*DN(3,0);
const double crLHS147 = C(4,5)*DN(3,2);
const double crLHS148 = C(1,5)*DN(3,1) + crLHS146 + crLHS147;
const double crLHS149 = DN(3,1)*crLHS42;
const double crLHS150 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS126;
const double crLHS151 = C(3,4)*DN(3,1);
const double crLHS152 = C(2,3)*DN(3,2) + crLHS146 + crLHS151;
const double crLHS153 = C(2,5)*DN(3,2);
const double crLHS154 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS153;
const double crLHS155 = DN(3,2)*crLHS42;
const double crLHS156 = DN(0,0)*N[3];
const double crLHS157 = DN(3,0)*N[0];
const double crLHS158 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS37;
const double crLHS159 = C(0,4)*DN(0,0) + crLHS40 + crLHS45;
const double crLHS160 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS161 = C(1,4)*DN(0,1);
const double crLHS162 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS161;
const double crLHS163 = pow(DN(0,1), 2);
const double crLHS164 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS161;
const double crLHS165 = C(2,4)*DN(0,2);
const double crLHS166 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS165;
const double crLHS167 = DN(0,1)*crLHS12;
const double crLHS168 = DN(0,2)*crLHS167;
const double crLHS169 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS73;
const double crLHS170 = C(0,4)*DN(1,0) + crLHS76 + crLHS80;
const double crLHS171 = DN(1,0)*crLHS167;
const double crLHS172 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS173 = C(1,4)*DN(1,1);
const double crLHS174 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS173;
const double crLHS175 = DN(0,1)*DN(1,1);
const double crLHS176 = crLHS12*crLHS175;
const double crLHS177 = crLHS58 + crLHS71;
const double crLHS178 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS173;
const double crLHS179 = C(2,4)*DN(1,2);
const double crLHS180 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS179;
const double crLHS181 = DN(1,2)*crLHS167;
const double crLHS182 = DN(0,1)*N[1];
const double crLHS183 = DN(1,1)*N[0];
const double crLHS184 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS109;
const double crLHS185 = C(0,4)*DN(2,0) + crLHS112 + crLHS116;
const double crLHS186 = DN(2,0)*crLHS167;
const double crLHS187 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS188 = C(1,4)*DN(2,1);
const double crLHS189 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS188;
const double crLHS190 = DN(0,1)*DN(2,1);
const double crLHS191 = crLHS12*crLHS190;
const double crLHS192 = crLHS107 + crLHS94;
const double crLHS193 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS188;
const double crLHS194 = C(2,4)*DN(2,2);
const double crLHS195 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS194;
const double crLHS196 = DN(2,2)*crLHS167;
const double crLHS197 = DN(0,1)*N[2];
const double crLHS198 = DN(2,1)*N[0];
const double crLHS199 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS144;
const double crLHS200 = C(0,4)*DN(3,0) + crLHS147 + crLHS151;
const double crLHS201 = DN(3,0)*crLHS167;
const double crLHS202 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS203 = C(1,4)*DN(3,1);
const double crLHS204 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS203;
const double crLHS205 = DN(0,1)*DN(3,1);
const double crLHS206 = crLHS12*crLHS205;
const double crLHS207 = crLHS129 + crLHS142;
const double crLHS208 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS203;
const double crLHS209 = C(2,4)*DN(3,2);
const double crLHS210 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS209;
const double crLHS211 = DN(3,2)*crLHS167;
const double crLHS212 = DN(0,1)*N[3];
const double crLHS213 = DN(3,1)*N[0];
const double crLHS214 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS47;
const double crLHS215 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS165;
const double crLHS216 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS217 = pow(DN(0,2), 2);
const double crLHS218 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS82;
const double crLHS219 = DN(0,2)*crLHS12;
const double crLHS220 = DN(1,0)*crLHS219;
const double crLHS221 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS179;
const double crLHS222 = DN(1,1)*crLHS219;
const double crLHS223 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS224 = DN(0,2)*DN(1,2);
const double crLHS225 = crLHS12*crLHS224;
const double crLHS226 = DN(0,2)*N[1];
const double crLHS227 = DN(1,2)*N[0];
const double crLHS228 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS118;
const double crLHS229 = DN(2,0)*crLHS219;
const double crLHS230 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS194;
const double crLHS231 = DN(2,1)*crLHS219;
const double crLHS232 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS233 = DN(0,2)*DN(2,2);
const double crLHS234 = crLHS12*crLHS233;
const double crLHS235 = DN(0,2)*N[2];
const double crLHS236 = DN(2,2)*N[0];
const double crLHS237 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS153;
const double crLHS238 = DN(3,0)*crLHS219;
const double crLHS239 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS209;
const double crLHS240 = DN(3,1)*crLHS219;
const double crLHS241 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS242 = DN(0,2)*DN(3,2);
const double crLHS243 = crLHS12*crLHS242;
const double crLHS244 = DN(0,2)*N[3];
const double crLHS245 = DN(3,2)*N[0];
const double crLHS246 = gauss_weight*(N[0] + crLHS26*(1.0*crLHS20 + 1.0*crLHS22 + 1.0*crLHS23 + 1.0*crLHS24 + crLHS30));
const double crLHS247 = crLHS27*gauss_weight;
const double crLHS248 = crLHS247*(crLHS175 + crLHS224 + crLHS57);
const double crLHS249 = crLHS247*(crLHS190 + crLHS233 + crLHS93);
const double crLHS250 = crLHS247*(crLHS128 + crLHS205 + crLHS242);
const double crLHS251 = N[1]*rho;
const double crLHS252 = crLHS28*crLHS63;
const double crLHS253 = 1.0*crLHS65;
const double crLHS254 = crLHS253*crLHS26;
const double crLHS255 = crLHS251*crLHS33;
const double crLHS256 = crLHS17*crLHS251 + crLHS25*crLHS252 - crLHS25*crLHS254 + crLHS255*crLHS32;
const double crLHS257 = pow(DN(1,0), 2);
const double crLHS258 = pow(N[1], 2);
const double crLHS259 = crLHS19*crLHS258 + crLHS251*crLHS63 + crLHS252*crLHS69 - crLHS254*crLHS69 + crLHS255*crLHS70 + crLHS258*crLHS6;
const double crLHS260 = DN(1,0)*crLHS12;
const double crLHS261 = DN(1,1)*crLHS260;
const double crLHS262 = DN(1,2)*crLHS260;
const double crLHS263 = gauss_weight*(-N[1] + crLHS251*crLHS50 + crLHS252 - crLHS254);
const double crLHS264 = DN(1,0)*DN(2,0);
const double crLHS265 = N[2]*crLHS64 + N[2]*crLHS65;
const double crLHS266 = crLHS12*crLHS264 + crLHS265;
const double crLHS267 = crLHS105*crLHS252 - crLHS105*crLHS254 + crLHS106*crLHS255 + crLHS251*crLHS99;
const double crLHS268 = DN(2,1)*crLHS260;
const double crLHS269 = DN(2,2)*crLHS260;
const double crLHS270 = DN(1,0)*N[2];
const double crLHS271 = DN(2,0)*N[1];
const double crLHS272 = DN(1,0)*DN(3,0);
const double crLHS273 = N[3]*crLHS64 + N[3]*crLHS65;
const double crLHS274 = crLHS12*crLHS272 + crLHS273;
const double crLHS275 = crLHS134*crLHS251 + crLHS140*crLHS252 - crLHS140*crLHS254 + crLHS141*crLHS255;
const double crLHS276 = DN(3,1)*crLHS260;
const double crLHS277 = DN(3,2)*crLHS260;
const double crLHS278 = DN(1,0)*N[3];
const double crLHS279 = DN(3,0)*N[1];
const double crLHS280 = crLHS256 + crLHS58;
const double crLHS281 = pow(DN(1,1), 2);
const double crLHS282 = DN(1,1)*crLHS12;
const double crLHS283 = DN(1,2)*crLHS282;
const double crLHS284 = DN(2,0)*crLHS282;
const double crLHS285 = DN(1,1)*DN(2,1);
const double crLHS286 = crLHS12*crLHS285;
const double crLHS287 = crLHS265 + crLHS267;
const double crLHS288 = DN(2,2)*crLHS282;
const double crLHS289 = DN(1,1)*N[2];
const double crLHS290 = DN(2,1)*N[1];
const double crLHS291 = DN(3,0)*crLHS282;
const double crLHS292 = DN(1,1)*DN(3,1);
const double crLHS293 = crLHS12*crLHS292;
const double crLHS294 = crLHS273 + crLHS275;
const double crLHS295 = DN(3,2)*crLHS282;
const double crLHS296 = DN(1,1)*N[3];
const double crLHS297 = DN(3,1)*N[1];
const double crLHS298 = pow(DN(1,2), 2);
const double crLHS299 = DN(1,2)*crLHS12;
const double crLHS300 = DN(2,0)*crLHS299;
const double crLHS301 = DN(2,1)*crLHS299;
const double crLHS302 = DN(1,2)*DN(2,2);
const double crLHS303 = crLHS12*crLHS302;
const double crLHS304 = DN(1,2)*N[2];
const double crLHS305 = DN(2,2)*N[1];
const double crLHS306 = DN(3,0)*crLHS299;
const double crLHS307 = DN(3,1)*crLHS299;
const double crLHS308 = DN(1,2)*DN(3,2);
const double crLHS309 = crLHS12*crLHS308;
const double crLHS310 = DN(1,2)*N[3];
const double crLHS311 = DN(3,2)*N[1];
const double crLHS312 = gauss_weight*(N[1] + crLHS26*(crLHS253 + 1.0*crLHS64 + 1.0*crLHS66 + 1.0*crLHS67 + 1.0*crLHS68));
const double crLHS313 = crLHS247*(crLHS264 + crLHS285 + crLHS302);
const double crLHS314 = crLHS247*(crLHS272 + crLHS292 + crLHS308);
const double crLHS315 = N[2]*rho;
const double crLHS316 = crLHS28*crLHS99;
const double crLHS317 = 1.0*crLHS101;
const double crLHS318 = crLHS26*crLHS317;
const double crLHS319 = crLHS315*crLHS33;
const double crLHS320 = crLHS17*crLHS315 + crLHS25*crLHS316 - crLHS25*crLHS318 + crLHS319*crLHS32;
const double crLHS321 = crLHS315*crLHS63 + crLHS316*crLHS69 - crLHS318*crLHS69 + crLHS319*crLHS70;
const double crLHS322 = pow(DN(2,0), 2);
const double crLHS323 = pow(N[2], 2);
const double crLHS324 = crLHS105*crLHS316 - crLHS105*crLHS318 + crLHS106*crLHS319 + crLHS19*crLHS323 + crLHS315*crLHS99 + crLHS323*crLHS6;
const double crLHS325 = DN(2,0)*crLHS12;
const double crLHS326 = DN(2,1)*crLHS325;
const double crLHS327 = DN(2,2)*crLHS325;
const double crLHS328 = gauss_weight*(-N[2] + crLHS315*crLHS50 + crLHS316 - crLHS318);
const double crLHS329 = DN(2,0)*DN(3,0);
const double crLHS330 = N[3]*crLHS100 + N[3]*crLHS101;
const double crLHS331 = crLHS12*crLHS329 + crLHS330;
const double crLHS332 = crLHS134*crLHS315 + crLHS140*crLHS316 - crLHS140*crLHS318 + crLHS141*crLHS319;
const double crLHS333 = DN(3,1)*crLHS325;
const double crLHS334 = DN(3,2)*crLHS325;
const double crLHS335 = DN(2,0)*N[3];
const double crLHS336 = DN(3,0)*N[2];
const double crLHS337 = crLHS320 + crLHS94;
const double crLHS338 = crLHS265 + crLHS321;
const double crLHS339 = pow(DN(2,1), 2);
const double crLHS340 = DN(2,1)*crLHS12;
const double crLHS341 = DN(2,2)*crLHS340;
const double crLHS342 = DN(3,0)*crLHS340;
const double crLHS343 = DN(2,1)*DN(3,1);
const double crLHS344 = crLHS12*crLHS343;
const double crLHS345 = crLHS330 + crLHS332;
const double crLHS346 = DN(3,2)*crLHS340;
const double crLHS347 = DN(2,1)*N[3];
const double crLHS348 = DN(3,1)*N[2];
const double crLHS349 = pow(DN(2,2), 2);
const double crLHS350 = DN(2,2)*crLHS12;
const double crLHS351 = DN(3,0)*crLHS350;
const double crLHS352 = DN(3,1)*crLHS350;
const double crLHS353 = DN(2,2)*DN(3,2);
const double crLHS354 = crLHS12*crLHS353;
const double crLHS355 = DN(2,2)*N[3];
const double crLHS356 = DN(3,2)*N[2];
const double crLHS357 = gauss_weight*(N[2] + crLHS26*(1.0*crLHS100 + 1.0*crLHS102 + 1.0*crLHS103 + 1.0*crLHS104 + crLHS317));
const double crLHS358 = crLHS247*(crLHS329 + crLHS343 + crLHS353);
const double crLHS359 = N[3]*rho;
const double crLHS360 = crLHS134*crLHS28;
const double crLHS361 = 1.0*crLHS136;
const double crLHS362 = crLHS26*crLHS361;
const double crLHS363 = crLHS33*crLHS359;
const double crLHS364 = crLHS17*crLHS359 + crLHS25*crLHS360 - crLHS25*crLHS362 + crLHS32*crLHS363;
const double crLHS365 = crLHS359*crLHS63 + crLHS360*crLHS69 - crLHS362*crLHS69 + crLHS363*crLHS70;
const double crLHS366 = crLHS105*crLHS360 - crLHS105*crLHS362 + crLHS106*crLHS363 + crLHS359*crLHS99;
const double crLHS367 = pow(DN(3,0), 2);
const double crLHS368 = pow(N[3], 2);
const double crLHS369 = crLHS134*crLHS359 + crLHS140*crLHS360 - crLHS140*crLHS362 + crLHS141*crLHS363 + crLHS19*crLHS368 + crLHS368*crLHS6;
const double crLHS370 = DN(3,0)*crLHS12;
const double crLHS371 = DN(3,1)*crLHS370;
const double crLHS372 = DN(3,2)*crLHS370;
const double crLHS373 = gauss_weight*(-N[3] + crLHS359*crLHS50 + crLHS360 - crLHS362);
const double crLHS374 = crLHS129 + crLHS364;
const double crLHS375 = crLHS273 + crLHS365;
const double crLHS376 = crLHS330 + crLHS366;
const double crLHS377 = pow(DN(3,1), 2);
const double crLHS378 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS379 = pow(DN(3,2), 2);
const double crLHS380 = gauss_weight*(N[3] + crLHS26*(1.0*crLHS135 + 1.0*crLHS137 + 1.0*crLHS138 + 1.0*crLHS139 + crLHS361));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS35);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS36 + DN(0,1)*crLHS38 + DN(0,2)*crLHS41 + crLHS43);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS44 + DN(0,1)*crLHS46 + DN(0,2)*crLHS48 + crLHS49);
rLHS(0,3)+=DN(0,0)*crLHS51;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS52 + DN(0,1)*crLHS54 + DN(0,2)*crLHS56 + crLHS59 + crLHS71);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS72 + DN(0,1)*crLHS74 + DN(0,2)*crLHS77 + crLHS78);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS81 + DN(0,2)*crLHS83 + crLHS84);
rLHS(0,7)+=-gauss_weight*(-DN(1,0)*crLHS29 + DN(1,0)*crLHS31 + crLHS85 - crLHS86*crLHS87);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS88 + DN(0,1)*crLHS90 + DN(0,2)*crLHS92 + crLHS107 + crLHS95);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS108 + DN(0,1)*crLHS110 + DN(0,2)*crLHS113 + crLHS114);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS115 + DN(0,1)*crLHS117 + DN(0,2)*crLHS119 + crLHS120);
rLHS(0,11)+=-gauss_weight*(-DN(2,0)*crLHS29 + DN(2,0)*crLHS31 + crLHS121 - crLHS122*crLHS87);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS123 + DN(0,1)*crLHS125 + DN(0,2)*crLHS127 + crLHS130 + crLHS142);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS143 + DN(0,1)*crLHS145 + DN(0,2)*crLHS148 + crLHS149);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS150 + DN(0,1)*crLHS152 + DN(0,2)*crLHS154 + crLHS155);
rLHS(0,15)+=-gauss_weight*(-DN(3,0)*crLHS29 + DN(3,0)*crLHS31 + crLHS156 - crLHS157*crLHS87);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS158 + DN(0,2)*crLHS159 + crLHS43);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS160 + DN(0,2)*crLHS162 + crLHS12*crLHS163 + crLHS35);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS46 + DN(0,1)*crLHS164 + DN(0,2)*crLHS166 + crLHS168);
rLHS(1,3)+=DN(0,1)*crLHS51;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS169 + DN(0,2)*crLHS170 + crLHS171);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS172 + DN(0,2)*crLHS174 + crLHS176 + crLHS177);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS81 + DN(0,1)*crLHS178 + DN(0,2)*crLHS180 + crLHS181);
rLHS(1,7)+=-gauss_weight*(-DN(1,1)*crLHS29 + DN(1,1)*crLHS31 + crLHS182 - crLHS183*crLHS87);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS90 + DN(0,1)*crLHS184 + DN(0,2)*crLHS185 + crLHS186);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS187 + DN(0,2)*crLHS189 + crLHS191 + crLHS192);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS117 + DN(0,1)*crLHS193 + DN(0,2)*crLHS195 + crLHS196);
rLHS(1,11)+=-gauss_weight*(-DN(2,1)*crLHS29 + DN(2,1)*crLHS31 + crLHS197 - crLHS198*crLHS87);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS125 + DN(0,1)*crLHS199 + DN(0,2)*crLHS200 + crLHS201);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS145 + DN(0,1)*crLHS202 + DN(0,2)*crLHS204 + crLHS206 + crLHS207);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS152 + DN(0,1)*crLHS208 + DN(0,2)*crLHS210 + crLHS211);
rLHS(1,15)+=-gauss_weight*(-DN(3,1)*crLHS29 + DN(3,1)*crLHS31 + crLHS212 - crLHS213*crLHS87);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS159 + DN(0,2)*crLHS214 + crLHS49);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS41 + DN(0,1)*crLHS162 + DN(0,2)*crLHS215 + crLHS168);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS48 + DN(0,1)*crLHS166 + DN(0,2)*crLHS216 + crLHS12*crLHS217 + crLHS35);
rLHS(2,3)+=DN(0,2)*crLHS51;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS56 + DN(0,1)*crLHS170 + DN(0,2)*crLHS218 + crLHS220);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS174 + DN(0,2)*crLHS221 + crLHS222);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS83 + DN(0,1)*crLHS180 + DN(0,2)*crLHS223 + crLHS177 + crLHS225);
rLHS(2,7)+=-gauss_weight*(-DN(1,2)*crLHS29 + DN(1,2)*crLHS31 + crLHS226 - crLHS227*crLHS87);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS92 + DN(0,1)*crLHS185 + DN(0,2)*crLHS228 + crLHS229);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS113 + DN(0,1)*crLHS189 + DN(0,2)*crLHS230 + crLHS231);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS119 + DN(0,1)*crLHS195 + DN(0,2)*crLHS232 + crLHS192 + crLHS234);
rLHS(2,11)+=-gauss_weight*(-DN(2,2)*crLHS29 + DN(2,2)*crLHS31 + crLHS235 - crLHS236*crLHS87);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS200 + DN(0,2)*crLHS237 + crLHS238);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS148 + DN(0,1)*crLHS204 + DN(0,2)*crLHS239 + crLHS240);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS154 + DN(0,1)*crLHS210 + DN(0,2)*crLHS241 + crLHS207 + crLHS243);
rLHS(2,15)+=-gauss_weight*(-DN(3,2)*crLHS29 + DN(3,2)*crLHS31 + crLHS244 - crLHS245*crLHS87);
rLHS(3,0)+=DN(0,0)*crLHS246;
rLHS(3,1)+=DN(0,1)*crLHS246;
rLHS(3,2)+=DN(0,2)*crLHS246;
rLHS(3,3)+=crLHS247*(crLHS163 + crLHS217 + crLHS5);
rLHS(3,4)+=gauss_weight*(DN(0,0)*crLHS70 + crLHS86);
rLHS(3,5)+=gauss_weight*(DN(0,1)*crLHS70 + crLHS183);
rLHS(3,6)+=gauss_weight*(DN(0,2)*crLHS70 + crLHS227);
rLHS(3,7)+=crLHS248;
rLHS(3,8)+=gauss_weight*(DN(0,0)*crLHS106 + crLHS122);
rLHS(3,9)+=gauss_weight*(DN(0,1)*crLHS106 + crLHS198);
rLHS(3,10)+=gauss_weight*(DN(0,2)*crLHS106 + crLHS236);
rLHS(3,11)+=crLHS249;
rLHS(3,12)+=gauss_weight*(DN(0,0)*crLHS141 + crLHS157);
rLHS(3,13)+=gauss_weight*(DN(0,1)*crLHS141 + crLHS213);
rLHS(3,14)+=gauss_weight*(DN(0,2)*crLHS141 + crLHS245);
rLHS(3,15)+=crLHS250;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS256 + crLHS59);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS36 + DN(1,1)*crLHS38 + DN(1,2)*crLHS41 + crLHS171);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS44 + DN(1,1)*crLHS46 + DN(1,2)*crLHS48 + crLHS220);
rLHS(4,3)+=gauss_weight*(DN(0,0)*crLHS252 - DN(0,0)*crLHS254 + crLHS85*crLHS87 - crLHS86);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS52 + DN(1,1)*crLHS54 + DN(1,2)*crLHS56 + crLHS12*crLHS257 + crLHS259);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS72 + DN(1,1)*crLHS74 + DN(1,2)*crLHS77 + crLHS261);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS81 + DN(1,2)*crLHS83 + crLHS262);
rLHS(4,7)+=DN(1,0)*crLHS263;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS88 + DN(1,1)*crLHS90 + DN(1,2)*crLHS92 + crLHS266 + crLHS267);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS108 + DN(1,1)*crLHS110 + DN(1,2)*crLHS113 + crLHS268);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS115 + DN(1,1)*crLHS117 + DN(1,2)*crLHS119 + crLHS269);
rLHS(4,11)+=-gauss_weight*(-DN(2,0)*crLHS252 + DN(2,0)*crLHS254 + crLHS270 - crLHS271*crLHS87);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS123 + DN(1,1)*crLHS125 + DN(1,2)*crLHS127 + crLHS274 + crLHS275);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS143 + DN(1,1)*crLHS145 + DN(1,2)*crLHS148 + crLHS276);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS150 + DN(1,1)*crLHS152 + DN(1,2)*crLHS154 + crLHS277);
rLHS(4,15)+=-gauss_weight*(-DN(3,0)*crLHS252 + DN(3,0)*crLHS254 + crLHS278 - crLHS279*crLHS87);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS158 + DN(1,2)*crLHS159 + crLHS78);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS160 + DN(1,2)*crLHS162 + crLHS176 + crLHS280);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS46 + DN(1,1)*crLHS164 + DN(1,2)*crLHS166 + crLHS222);
rLHS(5,3)+=gauss_weight*(DN(0,1)*crLHS252 - DN(0,1)*crLHS254 + crLHS182*crLHS87 - crLHS183);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS169 + DN(1,2)*crLHS170 + crLHS261);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS172 + DN(1,2)*crLHS174 + crLHS12*crLHS281 + crLHS259);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS81 + DN(1,1)*crLHS178 + DN(1,2)*crLHS180 + crLHS283);
rLHS(5,7)+=DN(1,1)*crLHS263;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS90 + DN(1,1)*crLHS184 + DN(1,2)*crLHS185 + crLHS284);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS187 + DN(1,2)*crLHS189 + crLHS286 + crLHS287);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS117 + DN(1,1)*crLHS193 + DN(1,2)*crLHS195 + crLHS288);
rLHS(5,11)+=-gauss_weight*(-DN(2,1)*crLHS252 + DN(2,1)*crLHS254 + crLHS289 - crLHS290*crLHS87);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS125 + DN(1,1)*crLHS199 + DN(1,2)*crLHS200 + crLHS291);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS145 + DN(1,1)*crLHS202 + DN(1,2)*crLHS204 + crLHS293 + crLHS294);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS152 + DN(1,1)*crLHS208 + DN(1,2)*crLHS210 + crLHS295);
rLHS(5,15)+=-gauss_weight*(-DN(3,1)*crLHS252 + DN(3,1)*crLHS254 + crLHS296 - crLHS297*crLHS87);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS159 + DN(1,2)*crLHS214 + crLHS84);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS41 + DN(1,1)*crLHS162 + DN(1,2)*crLHS215 + crLHS181);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS48 + DN(1,1)*crLHS166 + DN(1,2)*crLHS216 + crLHS225 + crLHS280);
rLHS(6,3)+=gauss_weight*(DN(0,2)*crLHS252 - DN(0,2)*crLHS254 + crLHS226*crLHS87 - crLHS227);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS56 + DN(1,1)*crLHS170 + DN(1,2)*crLHS218 + crLHS262);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS174 + DN(1,2)*crLHS221 + crLHS283);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS83 + DN(1,1)*crLHS180 + DN(1,2)*crLHS223 + crLHS12*crLHS298 + crLHS259);
rLHS(6,7)+=DN(1,2)*crLHS263;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS92 + DN(1,1)*crLHS185 + DN(1,2)*crLHS228 + crLHS300);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS113 + DN(1,1)*crLHS189 + DN(1,2)*crLHS230 + crLHS301);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS119 + DN(1,1)*crLHS195 + DN(1,2)*crLHS232 + crLHS287 + crLHS303);
rLHS(6,11)+=-gauss_weight*(-DN(2,2)*crLHS252 + DN(2,2)*crLHS254 + crLHS304 - crLHS305*crLHS87);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS200 + DN(1,2)*crLHS237 + crLHS306);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS148 + DN(1,1)*crLHS204 + DN(1,2)*crLHS239 + crLHS307);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS154 + DN(1,1)*crLHS210 + DN(1,2)*crLHS241 + crLHS294 + crLHS309);
rLHS(6,15)+=-gauss_weight*(-DN(3,2)*crLHS252 + DN(3,2)*crLHS254 + crLHS310 - crLHS311*crLHS87);
rLHS(7,0)+=gauss_weight*(DN(1,0)*crLHS32 + crLHS85);
rLHS(7,1)+=gauss_weight*(DN(1,1)*crLHS32 + crLHS182);
rLHS(7,2)+=gauss_weight*(DN(1,2)*crLHS32 + crLHS226);
rLHS(7,3)+=crLHS248;
rLHS(7,4)+=DN(1,0)*crLHS312;
rLHS(7,5)+=DN(1,1)*crLHS312;
rLHS(7,6)+=DN(1,2)*crLHS312;
rLHS(7,7)+=crLHS247*(crLHS257 + crLHS281 + crLHS298);
rLHS(7,8)+=gauss_weight*(DN(1,0)*crLHS106 + crLHS271);
rLHS(7,9)+=gauss_weight*(DN(1,1)*crLHS106 + crLHS290);
rLHS(7,10)+=gauss_weight*(DN(1,2)*crLHS106 + crLHS305);
rLHS(7,11)+=crLHS313;
rLHS(7,12)+=gauss_weight*(DN(1,0)*crLHS141 + crLHS279);
rLHS(7,13)+=gauss_weight*(DN(1,1)*crLHS141 + crLHS297);
rLHS(7,14)+=gauss_weight*(DN(1,2)*crLHS141 + crLHS311);
rLHS(7,15)+=crLHS314;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS320 + crLHS95);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS36 + DN(2,1)*crLHS38 + DN(2,2)*crLHS41 + crLHS186);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS44 + DN(2,1)*crLHS46 + DN(2,2)*crLHS48 + crLHS229);
rLHS(8,3)+=gauss_weight*(DN(0,0)*crLHS316 - DN(0,0)*crLHS318 + crLHS121*crLHS87 - crLHS122);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS52 + DN(2,1)*crLHS54 + DN(2,2)*crLHS56 + crLHS266 + crLHS321);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS72 + DN(2,1)*crLHS74 + DN(2,2)*crLHS77 + crLHS284);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS81 + DN(2,2)*crLHS83 + crLHS300);
rLHS(8,7)+=gauss_weight*(DN(1,0)*crLHS316 - DN(1,0)*crLHS318 + crLHS270*crLHS87 - crLHS271);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS88 + DN(2,1)*crLHS90 + DN(2,2)*crLHS92 + crLHS12*crLHS322 + crLHS324);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS108 + DN(2,1)*crLHS110 + DN(2,2)*crLHS113 + crLHS326);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS115 + DN(2,1)*crLHS117 + DN(2,2)*crLHS119 + crLHS327);
rLHS(8,11)+=DN(2,0)*crLHS328;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS123 + DN(2,1)*crLHS125 + DN(2,2)*crLHS127 + crLHS331 + crLHS332);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS143 + DN(2,1)*crLHS145 + DN(2,2)*crLHS148 + crLHS333);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS150 + DN(2,1)*crLHS152 + DN(2,2)*crLHS154 + crLHS334);
rLHS(8,15)+=-gauss_weight*(-DN(3,0)*crLHS316 + DN(3,0)*crLHS318 + crLHS335 - crLHS336*crLHS87);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS158 + DN(2,2)*crLHS159 + crLHS114);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS160 + DN(2,2)*crLHS162 + crLHS191 + crLHS337);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS46 + DN(2,1)*crLHS164 + DN(2,2)*crLHS166 + crLHS231);
rLHS(9,3)+=gauss_weight*(DN(0,1)*crLHS316 - DN(0,1)*crLHS318 + crLHS197*crLHS87 - crLHS198);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS169 + DN(2,2)*crLHS170 + crLHS268);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS172 + DN(2,2)*crLHS174 + crLHS286 + crLHS338);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS81 + DN(2,1)*crLHS178 + DN(2,2)*crLHS180 + crLHS301);
rLHS(9,7)+=gauss_weight*(DN(1,1)*crLHS316 - DN(1,1)*crLHS318 + crLHS289*crLHS87 - crLHS290);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS90 + DN(2,1)*crLHS184 + DN(2,2)*crLHS185 + crLHS326);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS187 + DN(2,2)*crLHS189 + crLHS12*crLHS339 + crLHS324);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS117 + DN(2,1)*crLHS193 + DN(2,2)*crLHS195 + crLHS341);
rLHS(9,11)+=DN(2,1)*crLHS328;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS125 + DN(2,1)*crLHS199 + DN(2,2)*crLHS200 + crLHS342);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS145 + DN(2,1)*crLHS202 + DN(2,2)*crLHS204 + crLHS344 + crLHS345);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS152 + DN(2,1)*crLHS208 + DN(2,2)*crLHS210 + crLHS346);
rLHS(9,15)+=-gauss_weight*(-DN(3,1)*crLHS316 + DN(3,1)*crLHS318 + crLHS347 - crLHS348*crLHS87);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS159 + DN(2,2)*crLHS214 + crLHS120);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS41 + DN(2,1)*crLHS162 + DN(2,2)*crLHS215 + crLHS196);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS48 + DN(2,1)*crLHS166 + DN(2,2)*crLHS216 + crLHS234 + crLHS337);
rLHS(10,3)+=gauss_weight*(DN(0,2)*crLHS316 - DN(0,2)*crLHS318 + crLHS235*crLHS87 - crLHS236);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS56 + DN(2,1)*crLHS170 + DN(2,2)*crLHS218 + crLHS269);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS174 + DN(2,2)*crLHS221 + crLHS288);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS83 + DN(2,1)*crLHS180 + DN(2,2)*crLHS223 + crLHS303 + crLHS338);
rLHS(10,7)+=gauss_weight*(DN(1,2)*crLHS316 - DN(1,2)*crLHS318 + crLHS304*crLHS87 - crLHS305);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS92 + DN(2,1)*crLHS185 + DN(2,2)*crLHS228 + crLHS327);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS113 + DN(2,1)*crLHS189 + DN(2,2)*crLHS230 + crLHS341);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS119 + DN(2,1)*crLHS195 + DN(2,2)*crLHS232 + crLHS12*crLHS349 + crLHS324);
rLHS(10,11)+=DN(2,2)*crLHS328;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS200 + DN(2,2)*crLHS237 + crLHS351);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS148 + DN(2,1)*crLHS204 + DN(2,2)*crLHS239 + crLHS352);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS154 + DN(2,1)*crLHS210 + DN(2,2)*crLHS241 + crLHS345 + crLHS354);
rLHS(10,15)+=-gauss_weight*(-DN(3,2)*crLHS316 + DN(3,2)*crLHS318 + crLHS355 - crLHS356*crLHS87);
rLHS(11,0)+=gauss_weight*(DN(2,0)*crLHS32 + crLHS121);
rLHS(11,1)+=gauss_weight*(DN(2,1)*crLHS32 + crLHS197);
rLHS(11,2)+=gauss_weight*(DN(2,2)*crLHS32 + crLHS235);
rLHS(11,3)+=crLHS249;
rLHS(11,4)+=gauss_weight*(DN(2,0)*crLHS70 + crLHS270);
rLHS(11,5)+=gauss_weight*(DN(2,1)*crLHS70 + crLHS289);
rLHS(11,6)+=gauss_weight*(DN(2,2)*crLHS70 + crLHS304);
rLHS(11,7)+=crLHS313;
rLHS(11,8)+=DN(2,0)*crLHS357;
rLHS(11,9)+=DN(2,1)*crLHS357;
rLHS(11,10)+=DN(2,2)*crLHS357;
rLHS(11,11)+=crLHS247*(crLHS322 + crLHS339 + crLHS349);
rLHS(11,12)+=gauss_weight*(DN(2,0)*crLHS141 + crLHS336);
rLHS(11,13)+=gauss_weight*(DN(2,1)*crLHS141 + crLHS348);
rLHS(11,14)+=gauss_weight*(DN(2,2)*crLHS141 + crLHS356);
rLHS(11,15)+=crLHS358;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS130 + crLHS364);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS36 + DN(3,1)*crLHS38 + DN(3,2)*crLHS41 + crLHS201);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS44 + DN(3,1)*crLHS46 + DN(3,2)*crLHS48 + crLHS238);
rLHS(12,3)+=gauss_weight*(DN(0,0)*crLHS360 - DN(0,0)*crLHS362 + crLHS156*crLHS87 - crLHS157);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS52 + DN(3,1)*crLHS54 + DN(3,2)*crLHS56 + crLHS274 + crLHS365);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS72 + DN(3,1)*crLHS74 + DN(3,2)*crLHS77 + crLHS291);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS79 + DN(3,1)*crLHS81 + DN(3,2)*crLHS83 + crLHS306);
rLHS(12,7)+=gauss_weight*(DN(1,0)*crLHS360 - DN(1,0)*crLHS362 + crLHS278*crLHS87 - crLHS279);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS88 + DN(3,1)*crLHS90 + DN(3,2)*crLHS92 + crLHS331 + crLHS366);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS108 + DN(3,1)*crLHS110 + DN(3,2)*crLHS113 + crLHS342);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS115 + DN(3,1)*crLHS117 + DN(3,2)*crLHS119 + crLHS351);
rLHS(12,11)+=gauss_weight*(DN(2,0)*crLHS360 - DN(2,0)*crLHS362 + crLHS335*crLHS87 - crLHS336);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS123 + DN(3,1)*crLHS125 + DN(3,2)*crLHS127 + crLHS12*crLHS367 + crLHS369);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS143 + DN(3,1)*crLHS145 + DN(3,2)*crLHS148 + crLHS371);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS150 + DN(3,1)*crLHS152 + DN(3,2)*crLHS154 + crLHS372);
rLHS(12,15)+=DN(3,0)*crLHS373;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS158 + DN(3,2)*crLHS159 + crLHS149);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS38 + DN(3,1)*crLHS160 + DN(3,2)*crLHS162 + crLHS206 + crLHS374);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS46 + DN(3,1)*crLHS164 + DN(3,2)*crLHS166 + crLHS240);
rLHS(13,3)+=gauss_weight*(DN(0,1)*crLHS360 - DN(0,1)*crLHS362 + crLHS212*crLHS87 - crLHS213);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS54 + DN(3,1)*crLHS169 + DN(3,2)*crLHS170 + crLHS276);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS74 + DN(3,1)*crLHS172 + DN(3,2)*crLHS174 + crLHS293 + crLHS375);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS81 + DN(3,1)*crLHS178 + DN(3,2)*crLHS180 + crLHS307);
rLHS(13,7)+=gauss_weight*(DN(1,1)*crLHS360 - DN(1,1)*crLHS362 + crLHS296*crLHS87 - crLHS297);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS90 + DN(3,1)*crLHS184 + DN(3,2)*crLHS185 + crLHS333);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS110 + DN(3,1)*crLHS187 + DN(3,2)*crLHS189 + crLHS344 + crLHS376);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS117 + DN(3,1)*crLHS193 + DN(3,2)*crLHS195 + crLHS352);
rLHS(13,11)+=gauss_weight*(DN(2,1)*crLHS360 - DN(2,1)*crLHS362 + crLHS347*crLHS87 - crLHS348);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS125 + DN(3,1)*crLHS199 + DN(3,2)*crLHS200 + crLHS371);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS145 + DN(3,1)*crLHS202 + DN(3,2)*crLHS204 + crLHS12*crLHS377 + crLHS369);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS152 + DN(3,1)*crLHS208 + DN(3,2)*crLHS210 + crLHS378);
rLHS(13,15)+=DN(3,1)*crLHS373;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS159 + DN(3,2)*crLHS214 + crLHS155);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS41 + DN(3,1)*crLHS162 + DN(3,2)*crLHS215 + crLHS211);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS48 + DN(3,1)*crLHS166 + DN(3,2)*crLHS216 + crLHS243 + crLHS374);
rLHS(14,3)+=gauss_weight*(DN(0,2)*crLHS360 - DN(0,2)*crLHS362 + crLHS244*crLHS87 - crLHS245);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS56 + DN(3,1)*crLHS170 + DN(3,2)*crLHS218 + crLHS277);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS174 + DN(3,2)*crLHS221 + crLHS295);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS83 + DN(3,1)*crLHS180 + DN(3,2)*crLHS223 + crLHS309 + crLHS375);
rLHS(14,7)+=gauss_weight*(DN(1,2)*crLHS360 - DN(1,2)*crLHS362 + crLHS310*crLHS87 - crLHS311);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS92 + DN(3,1)*crLHS185 + DN(3,2)*crLHS228 + crLHS334);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS113 + DN(3,1)*crLHS189 + DN(3,2)*crLHS230 + crLHS346);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS119 + DN(3,1)*crLHS195 + DN(3,2)*crLHS232 + crLHS354 + crLHS376);
rLHS(14,11)+=gauss_weight*(DN(2,2)*crLHS360 - DN(2,2)*crLHS362 + crLHS355*crLHS87 - crLHS356);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS200 + DN(3,2)*crLHS237 + crLHS372);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS148 + DN(3,1)*crLHS204 + DN(3,2)*crLHS239 + crLHS378);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS154 + DN(3,1)*crLHS210 + DN(3,2)*crLHS241 + crLHS12*crLHS379 + crLHS369);
rLHS(14,15)+=DN(3,2)*crLHS373;
rLHS(15,0)+=gauss_weight*(DN(3,0)*crLHS32 + crLHS156);
rLHS(15,1)+=gauss_weight*(DN(3,1)*crLHS32 + crLHS212);
rLHS(15,2)+=gauss_weight*(DN(3,2)*crLHS32 + crLHS244);
rLHS(15,3)+=crLHS250;
rLHS(15,4)+=gauss_weight*(DN(3,0)*crLHS70 + crLHS278);
rLHS(15,5)+=gauss_weight*(DN(3,1)*crLHS70 + crLHS296);
rLHS(15,6)+=gauss_weight*(DN(3,2)*crLHS70 + crLHS310);
rLHS(15,7)+=crLHS314;
rLHS(15,8)+=gauss_weight*(DN(3,0)*crLHS106 + crLHS335);
rLHS(15,9)+=gauss_weight*(DN(3,1)*crLHS106 + crLHS347);
rLHS(15,10)+=gauss_weight*(DN(3,2)*crLHS106 + crLHS355);
rLHS(15,11)+=crLHS358;
rLHS(15,12)+=DN(3,0)*crLHS380;
rLHS(15,13)+=DN(3,1)*crLHS380;
rLHS(15,14)+=DN(3,2)*crLHS380;
rLHS(15,15)+=crLHS247*(crLHS367 + crLHS377 + crLHS379);

}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContribution(
    FluidTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,3>& c = rData.SoundVelocity;

    const array_1d<double,3> alpha = rData.Resistance;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3>& v = rData.Velocity;
    const BoundedMatrix<double,2,3>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,2,3>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,2,3> vconv = rData.Velocity - rData.MeshVelocity;
    const BoundedMatrix<double,2,3>& f = rData.BodyForce;
    const array_1d<double,3>& p = rData.Pressure;
    // const array_1d<double,3>& pn = rData.Pressure_OldStep1;
    // const array_1d<double,3>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,3>& stress = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // NAVIER-STOKES ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crRHS1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crRHS2 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crRHS3 = crRHS2*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
const double crRHS4 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crRHS5 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS7 = crRHS5*crRHS6;
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS9 = crRHS8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
const double crRHS10 = crRHS7 + crRHS9;
const double crRHS11 = N[0]*rho;
const double crRHS12 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS13 = crRHS12 + crRHS5;
const double crRHS14 = crRHS2*stab_c3;
const double crRHS15 = rho*stab_c2*sqrt(pow(crRHS6, 2) + pow(crRHS8, 2));
const double crRHS16 = crRHS13*(h*(crRHS14*h + crRHS15)/stab_c1 + mu);
const double crRHS17 = 1.0/(crRHS14 + crRHS15/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS18 = crRHS17*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS1 + crRHS3 + crRHS4 + crRHS7*rho + crRHS9*rho);
const double crRHS19 = N[0]*crRHS2;
const double crRHS20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crRHS21 = crRHS11*crRHS20;
const double crRHS22 = rho*(DN(0,0)*crRHS6 + DN(0,1)*crRHS8);
const double crRHS23 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS24 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crRHS25 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crRHS26 = crRHS6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
const double crRHS27 = crRHS12*crRHS8;
const double crRHS28 = crRHS26 + crRHS27;
const double crRHS29 = crRHS17*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS23 + crRHS24 + crRHS25 + crRHS26*rho + crRHS27*rho);
const double crRHS30 = N[1]*rho;
const double crRHS31 = N[1]*crRHS2;
const double crRHS32 = crRHS20*crRHS30;
const double crRHS33 = rho*(DN(1,0)*crRHS6 + DN(1,1)*crRHS8);
const double crRHS34 = N[2]*rho;
const double crRHS35 = N[2]*crRHS2;
const double crRHS36 = crRHS20*crRHS34;
const double crRHS37 = rho*(DN(2,0)*crRHS6 + DN(2,1)*crRHS8);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS16 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS1 + N[0]*crRHS3 + N[0]*crRHS4 + crRHS10*crRHS11 - crRHS18*crRHS19 + crRHS18*crRHS21 + crRHS18*crRHS22);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS16 + DN(0,1)*stress[1] - N[0]*crRHS23 + N[0]*crRHS24 + N[0]*crRHS25 + crRHS11*crRHS28 - crRHS19*crRHS29 + crRHS21*crRHS29 + crRHS22*crRHS29);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS18 + DN(0,1)*crRHS29 + N[0]*crRHS13);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS16 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS1 + N[1]*crRHS3 + N[1]*crRHS4 + crRHS10*crRHS30 - crRHS18*crRHS31 + crRHS18*crRHS32 + crRHS18*crRHS33);
rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS16 + DN(1,1)*stress[1] - N[1]*crRHS23 + N[1]*crRHS24 + N[1]*crRHS25 + crRHS28*crRHS30 - crRHS29*crRHS31 + crRHS29*crRHS32 + crRHS29*crRHS33);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS18 + DN(1,1)*crRHS29 + N[1]*crRHS13);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS16 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS1 + N[2]*crRHS3 + N[2]*crRHS4 + crRHS10*crRHS34 - crRHS18*crRHS35 + crRHS18*crRHS36 + crRHS18*crRHS37);
rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS16 + DN(2,1)*stress[1] - N[2]*crRHS23 + N[2]*crRHS24 + N[2]*crRHS25 + crRHS28*crRHS34 - crRHS29*crRHS35 + crRHS29*crRHS36 + crRHS29*crRHS37);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS18 + DN(2,1)*crRHS29 + N[2]*crRHS13);

}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContribution(
    FluidTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,4>& c = rData.SoundVelocity;

    const array_1d<double,4> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4>& v = rData.Velocity;
    const BoundedMatrix<double,3,4>& vn = rData.Velocity_OldStep1;
    const BoundedMatrix<double,3,4>& vnn = rData.Velocity_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh = rData.MeshVelocity;
    const BoundedMatrix<double,3,4> vconv = rData.Velocity - rData.MeshVelocity;
    const BoundedMatrix<double,3,4>& f = rData.BodyForce;
    const array_1d<double,4>& p = rData.Pressure;
    // const array_1d<double,4>& pn = rData.Pressure_OldStep1;
    // const array_1d<double,4>& pnn = rData.Pressure_OldStep2;
    const array_1d<double,6>& stress = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // NAVIER-STOKES ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crRHS1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS2 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crRHS3 = crRHS2*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crRHS4 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crRHS5 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS7 = crRHS5*crRHS6;
const double crRHS8 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS9 = crRHS8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0));
const double crRHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS11 = crRHS10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
const double crRHS12 = crRHS11 + crRHS7 + crRHS9;
const double crRHS13 = N[0]*rho;
const double crRHS14 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS15 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS16 = crRHS14 + crRHS15 + crRHS5;
const double crRHS17 = crRHS2*stab_c3;
const double crRHS18 = rho*stab_c2*sqrt(pow(crRHS10, 2) + pow(crRHS6, 2) + pow(crRHS8, 2));
const double crRHS19 = crRHS16*(h*(crRHS17*h + crRHS18)/stab_c1 + mu);
const double crRHS20 = 1.0/(crRHS17 + crRHS18/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS21 = crRHS20*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crRHS1 + crRHS11*rho + crRHS3 + crRHS4 + crRHS7*rho + crRHS9*rho);
const double crRHS22 = N[0]*crRHS2;
const double crRHS23 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crRHS24 = crRHS13*crRHS23;
const double crRHS25 = rho*(DN(0,0)*crRHS6 + DN(0,1)*crRHS8 + DN(0,2)*crRHS10);
const double crRHS26 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS27 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crRHS28 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crRHS29 = crRHS6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1));
const double crRHS30 = crRHS14*crRHS8;
const double crRHS31 = crRHS10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
const double crRHS32 = crRHS29 + crRHS30 + crRHS31;
const double crRHS33 = crRHS20*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS26 + crRHS27 + crRHS28 + crRHS29*rho + crRHS30*rho + crRHS31*rho);
const double crRHS34 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS35 = crRHS2*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crRHS36 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crRHS37 = crRHS6*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2));
const double crRHS38 = crRHS8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
const double crRHS39 = crRHS10*crRHS15;
const double crRHS40 = crRHS37 + crRHS38 + crRHS39;
const double crRHS41 = crRHS20*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS34 + crRHS35 + crRHS36 + crRHS37*rho + crRHS38*rho + crRHS39*rho);
const double crRHS42 = N[1]*rho;
const double crRHS43 = N[1]*crRHS2;
const double crRHS44 = crRHS23*crRHS42;
const double crRHS45 = rho*(DN(1,0)*crRHS6 + DN(1,1)*crRHS8 + DN(1,2)*crRHS10);
const double crRHS46 = N[2]*rho;
const double crRHS47 = N[2]*crRHS2;
const double crRHS48 = crRHS23*crRHS46;
const double crRHS49 = rho*(DN(2,0)*crRHS6 + DN(2,1)*crRHS8 + DN(2,2)*crRHS10);
const double crRHS50 = N[3]*rho;
const double crRHS51 = N[3]*crRHS2;
const double crRHS52 = crRHS23*crRHS50;
const double crRHS53 = rho*(DN(3,0)*crRHS6 + DN(3,1)*crRHS8 + DN(3,2)*crRHS10);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS19 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] - N[0]*crRHS1 + N[0]*crRHS3 + N[0]*crRHS4 + crRHS12*crRHS13 - crRHS21*crRHS22 + crRHS21*crRHS24 + crRHS21*crRHS25);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS19 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS26 + N[0]*crRHS27 + N[0]*crRHS28 + crRHS13*crRHS32 - crRHS22*crRHS33 + crRHS24*crRHS33 + crRHS25*crRHS33);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS19 + DN(0,2)*stress[2] - N[0]*crRHS34 + N[0]*crRHS35 + N[0]*crRHS36 + crRHS13*crRHS40 - crRHS22*crRHS41 + crRHS24*crRHS41 + crRHS25*crRHS41);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS21 + DN(0,1)*crRHS33 + DN(0,2)*crRHS41 + N[0]*crRHS16);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS19 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] - N[1]*crRHS1 + N[1]*crRHS3 + N[1]*crRHS4 + crRHS12*crRHS42 - crRHS21*crRHS43 + crRHS21*crRHS44 + crRHS21*crRHS45);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS19 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS26 + N[1]*crRHS27 + N[1]*crRHS28 + crRHS32*crRHS42 - crRHS33*crRHS43 + crRHS33*crRHS44 + crRHS33*crRHS45);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS19 + DN(1,2)*stress[2] - N[1]*crRHS34 + N[1]*crRHS35 + N[1]*crRHS36 + crRHS40*crRHS42 - crRHS41*crRHS43 + crRHS41*crRHS44 + crRHS41*crRHS45);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS21 + DN(1,1)*crRHS33 + DN(1,2)*crRHS41 + N[1]*crRHS16);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS19 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] - N[2]*crRHS1 + N[2]*crRHS3 + N[2]*crRHS4 + crRHS12*crRHS46 - crRHS21*crRHS47 + crRHS21*crRHS48 + crRHS21*crRHS49);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS19 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS26 + N[2]*crRHS27 + N[2]*crRHS28 + crRHS32*crRHS46 - crRHS33*crRHS47 + crRHS33*crRHS48 + crRHS33*crRHS49);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS19 + DN(2,2)*stress[2] - N[2]*crRHS34 + N[2]*crRHS35 + N[2]*crRHS36 + crRHS40*crRHS46 - crRHS41*crRHS47 + crRHS41*crRHS48 + crRHS41*crRHS49);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS21 + DN(2,1)*crRHS33 + DN(2,2)*crRHS41 + N[2]*crRHS16);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS19 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] - N[3]*crRHS1 + N[3]*crRHS3 + N[3]*crRHS4 + crRHS12*crRHS50 - crRHS21*crRHS51 + crRHS21*crRHS52 + crRHS21*crRHS53);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS19 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS26 + N[3]*crRHS27 + N[3]*crRHS28 + crRHS32*crRHS50 - crRHS33*crRHS51 + crRHS33*crRHS52 + crRHS33*crRHS53);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS19 + DN(3,2)*stress[2] - N[3]*crRHS34 + N[3]*crRHS35 + N[3]*crRHS36 + crRHS40*crRHS50 - crRHS41*crRHS51 + crRHS41*crRHS52 + crRHS41*crRHS53);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS21 + DN(3,1)*crRHS33 + DN(3,2)*crRHS41 + N[3]*crRHS16);

}

template <>
void FluidTopologyOptimizationElement< FluidTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointLHSContributionAdjoint(
    FluidTopologyOptimizationElementData<2,3,true>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,3>& c = rData.SoundVelocity;

    const array_1d<double,3> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;
    
    // Get constitutive matrix
    const BoundedMatrix<double,3,3>& C = rData.C;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,2,3> v_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = pow(DN(0,0), 2);
const double crLHS4 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crLHS5 = crLHS4*stab_c3;
const double crLHS6 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0);
const double crLHS7 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1);
const double crLHS8 = rho*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2));
const double crLHS9 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crLHS10 = DN(0,0)*v_ns(0,1);
const double crLHS11 = DN(1,0)*v_ns(1,1);
const double crLHS12 = DN(2,0)*v_ns(2,1);
const double crLHS13 = crLHS10 + crLHS11 + crLHS12;
const double crLHS14 = DN(0,1)*v_ns(0,0);
const double crLHS15 = DN(1,1)*v_ns(1,0);
const double crLHS16 = DN(2,1)*v_ns(2,0);
const double crLHS17 = crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crLHS19 = rho*stab_c3*sqrt(pow(crLHS13, 2) + pow(crLHS17, 2) + pow(crLHS18, 2) + pow(crLHS9, 2));
const double crLHS20 = h*(crLHS19*h + crLHS5*h + crLHS8)/stab_c1 + mu;
const double crLHS21 = pow(N[0], 2);
const double crLHS22 = crLHS21*rho;
const double crLHS23 = N[0]*crLHS4;
const double crLHS24 = N[0]*rho;
const double crLHS25 = crLHS24*crLHS9;
const double crLHS26 = DN(0,0)*crLHS6;
const double crLHS27 = crLHS26*rho;
const double crLHS28 = DN(0,1)*crLHS7;
const double crLHS29 = crLHS28*rho;
const double crLHS30 = -crLHS23 + crLHS27 + crLHS29;
const double crLHS31 = -crLHS25 + crLHS30;
const double crLHS32 = 1.0/(crLHS19 + crLHS5 + crLHS8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crLHS33 = 1.0*crLHS32;
const double crLHS34 = crLHS31*crLHS33;
const double crLHS35 = crLHS26 + crLHS28;
const double crLHS36 = crLHS35*rho;
const double crLHS37 = crLHS13*crLHS17;
const double crLHS38 = crLHS24*crLHS37;
const double crLHS39 = -crLHS31*crLHS9 + crLHS38;
const double crLHS40 = crLHS24*crLHS33;
const double crLHS41 = crLHS21*crLHS4;
const double crLHS42 = bdf0*crLHS22 + crLHS24*crLHS35 + crLHS41;
const double crLHS43 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS44 = C(1,2)*DN(0,1);
const double crLHS45 = C(2,2)*DN(0,0) + crLHS44;
const double crLHS46 = DN(0,0)*crLHS20;
const double crLHS47 = DN(0,1)*crLHS46;
const double crLHS48 = crLHS13*crLHS33;
const double crLHS49 = crLHS41*rho;
const double crLHS50 = pow(rho, 2);
const double crLHS51 = crLHS35*crLHS50;
const double crLHS52 = crLHS48*crLHS51;
const double crLHS53 = crLHS18*crLHS24;
const double crLHS54 = crLHS23 + crLHS25 - crLHS27 - crLHS29 + crLHS53;
const double crLHS55 = 1.0*crLHS10 + 1.0*crLHS11 + 1.0*crLHS12;
const double crLHS56 = crLHS24*crLHS32;
const double crLHS57 = crLHS55*crLHS56;
const double crLHS58 = DN(0,0)*N[0];
const double crLHS59 = DN(0,0)*crLHS33;
const double crLHS60 = DN(0,0)*crLHS9 + DN(0,1)*crLHS13;
const double crLHS61 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS62 = C(0,2)*DN(1,0);
const double crLHS63 = C(2,2)*DN(1,1) + crLHS62;
const double crLHS64 = N[1]*rho;
const double crLHS65 = crLHS64*crLHS9;
const double crLHS66 = N[1]*crLHS4;
const double crLHS67 = DN(1,0)*crLHS6;
const double crLHS68 = crLHS67*rho;
const double crLHS69 = DN(1,1)*crLHS7;
const double crLHS70 = crLHS69*rho;
const double crLHS71 = -crLHS66 + crLHS68 + crLHS70;
const double crLHS72 = -crLHS65 + crLHS71;
const double crLHS73 = crLHS33*crLHS72;
const double crLHS74 = crLHS37*crLHS64;
const double crLHS75 = -crLHS72*crLHS9 + crLHS74;
const double crLHS76 = N[1]*crLHS23;
const double crLHS77 = N[1]*crLHS24;
const double crLHS78 = bdf0*crLHS77 + crLHS76;
const double crLHS79 = crLHS35*crLHS64 + crLHS78;
const double crLHS80 = DN(0,0)*DN(1,0);
const double crLHS81 = N[1]*crLHS25 + crLHS20*crLHS80;
const double crLHS82 = C(0,1)*DN(1,1) + crLHS62;
const double crLHS83 = C(1,2)*DN(1,1);
const double crLHS84 = C(2,2)*DN(1,0) + crLHS83;
const double crLHS85 = DN(1,1)*crLHS46;
const double crLHS86 = crLHS18*crLHS64;
const double crLHS87 = crLHS65 + crLHS66 - crLHS68 - crLHS70 + crLHS86;
const double crLHS88 = crLHS48*rho;
const double crLHS89 = crLHS13*crLHS77 - crLHS76*crLHS88;
const double crLHS90 = DN(0,0)*N[1];
const double crLHS91 = DN(1,0)*crLHS33;
const double crLHS92 = DN(1,0)*crLHS9 + DN(1,1)*crLHS13;
const double crLHS93 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS94 = C(0,2)*DN(2,0);
const double crLHS95 = C(2,2)*DN(2,1) + crLHS94;
const double crLHS96 = N[2]*rho;
const double crLHS97 = crLHS9*crLHS96;
const double crLHS98 = N[2]*crLHS4;
const double crLHS99 = DN(2,0)*crLHS6;
const double crLHS100 = crLHS99*rho;
const double crLHS101 = DN(2,1)*crLHS7;
const double crLHS102 = crLHS101*rho;
const double crLHS103 = crLHS100 + crLHS102 - crLHS98;
const double crLHS104 = crLHS103 - crLHS97;
const double crLHS105 = crLHS104*crLHS33;
const double crLHS106 = crLHS37*crLHS96;
const double crLHS107 = -crLHS104*crLHS9 + crLHS106;
const double crLHS108 = N[2]*crLHS23;
const double crLHS109 = N[2]*crLHS24;
const double crLHS110 = bdf0*crLHS109 + crLHS108;
const double crLHS111 = crLHS110 + crLHS35*crLHS96;
const double crLHS112 = DN(0,0)*DN(2,0);
const double crLHS113 = N[2]*crLHS25 + crLHS112*crLHS20;
const double crLHS114 = C(0,1)*DN(2,1) + crLHS94;
const double crLHS115 = C(1,2)*DN(2,1);
const double crLHS116 = C(2,2)*DN(2,0) + crLHS115;
const double crLHS117 = DN(2,1)*crLHS46;
const double crLHS118 = crLHS18*crLHS96;
const double crLHS119 = -crLHS100 - crLHS102 + crLHS118 + crLHS97 + crLHS98;
const double crLHS120 = -crLHS108*crLHS88 + crLHS109*crLHS13;
const double crLHS121 = DN(0,0)*N[2];
const double crLHS122 = DN(2,0)*crLHS33;
const double crLHS123 = DN(2,0)*crLHS9 + DN(2,1)*crLHS13;
const double crLHS124 = C(0,1)*DN(0,0) + crLHS44;
const double crLHS125 = crLHS17*crLHS33;
const double crLHS126 = crLHS125*crLHS51;
const double crLHS127 = 1.0*crLHS14 + 1.0*crLHS15 + 1.0*crLHS16;
const double crLHS128 = crLHS127*crLHS56;
const double crLHS129 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS130 = pow(DN(0,1), 2);
const double crLHS131 = crLHS30 - crLHS53;
const double crLHS132 = crLHS131*crLHS33;
const double crLHS133 = -crLHS131*crLHS18 + crLHS38;
const double crLHS134 = DN(0,1)*N[0];
const double crLHS135 = DN(0,1)*crLHS33;
const double crLHS136 = DN(0,0)*crLHS17 + DN(0,1)*crLHS18;
const double crLHS137 = C(0,1)*DN(1,0) + crLHS83;
const double crLHS138 = DN(0,1)*crLHS20;
const double crLHS139 = DN(1,0)*crLHS138;
const double crLHS140 = crLHS125*rho;
const double crLHS141 = -crLHS140*crLHS76 + crLHS17*crLHS77;
const double crLHS142 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS143 = crLHS71 - crLHS86;
const double crLHS144 = crLHS143*crLHS33;
const double crLHS145 = -crLHS143*crLHS18 + crLHS74;
const double crLHS146 = DN(0,1)*DN(1,1);
const double crLHS147 = N[1]*crLHS53 + crLHS146*crLHS20;
const double crLHS148 = DN(0,1)*N[1];
const double crLHS149 = DN(1,1)*crLHS33;
const double crLHS150 = DN(1,0)*crLHS17 + DN(1,1)*crLHS18;
const double crLHS151 = C(0,1)*DN(2,0) + crLHS115;
const double crLHS152 = DN(2,0)*crLHS138;
const double crLHS153 = -crLHS108*crLHS140 + crLHS109*crLHS17;
const double crLHS154 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS155 = crLHS103 - crLHS118;
const double crLHS156 = crLHS155*crLHS33;
const double crLHS157 = crLHS106 - crLHS155*crLHS18;
const double crLHS158 = DN(0,1)*DN(2,1);
const double crLHS159 = N[2]*crLHS53 + crLHS158*crLHS20;
const double crLHS160 = DN(0,1)*N[2];
const double crLHS161 = DN(2,1)*crLHS33;
const double crLHS162 = DN(2,0)*crLHS17 + DN(2,1)*crLHS18;
const double crLHS163 = crLHS33*gauss_weight;
const double crLHS164 = DN(1,0)*N[0];
const double crLHS165 = DN(1,1)*N[0];
const double crLHS166 = crLHS163*(crLHS146 + crLHS80);
const double crLHS167 = DN(2,0)*N[0];
const double crLHS168 = DN(2,1)*N[0];
const double crLHS169 = crLHS163*(crLHS112 + crLHS158);
const double crLHS170 = crLHS67 + crLHS69;
const double crLHS171 = crLHS170*rho;
const double crLHS172 = crLHS33*crLHS64;
const double crLHS173 = crLHS170*crLHS24 + crLHS78;
const double crLHS174 = crLHS32*crLHS64;
const double crLHS175 = crLHS174*crLHS55;
const double crLHS176 = crLHS170*crLHS50;
const double crLHS177 = crLHS176*crLHS48;
const double crLHS178 = pow(DN(1,0), 2);
const double crLHS179 = pow(N[1], 2);
const double crLHS180 = crLHS179*rho;
const double crLHS181 = crLHS179*crLHS4;
const double crLHS182 = bdf0*crLHS180 + crLHS170*crLHS64 + crLHS181;
const double crLHS183 = DN(1,0)*crLHS20;
const double crLHS184 = DN(1,1)*crLHS183;
const double crLHS185 = DN(1,0)*N[1];
const double crLHS186 = N[2]*crLHS64;
const double crLHS187 = N[2]*crLHS66 + bdf0*crLHS186;
const double crLHS188 = crLHS170*crLHS96 + crLHS187;
const double crLHS189 = DN(1,0)*DN(2,0);
const double crLHS190 = N[2]*crLHS65 + crLHS189*crLHS20;
const double crLHS191 = DN(2,1)*crLHS183;
const double crLHS192 = crLHS33*crLHS96;
const double crLHS193 = crLHS192*crLHS66;
const double crLHS194 = crLHS13*crLHS186 - crLHS13*crLHS193;
const double crLHS195 = DN(1,0)*N[2];
const double crLHS196 = crLHS127*crLHS174;
const double crLHS197 = crLHS125*crLHS176;
const double crLHS198 = pow(DN(1,1), 2);
const double crLHS199 = DN(1,1)*N[1];
const double crLHS200 = DN(2,0)*crLHS20;
const double crLHS201 = DN(1,1)*crLHS200;
const double crLHS202 = crLHS17*crLHS186 - crLHS17*crLHS193;
const double crLHS203 = DN(1,1)*DN(2,1);
const double crLHS204 = N[2]*crLHS86 + crLHS20*crLHS203;
const double crLHS205 = DN(1,1)*N[2];
const double crLHS206 = DN(2,0)*N[1];
const double crLHS207 = DN(2,1)*N[1];
const double crLHS208 = crLHS163*(crLHS189 + crLHS203);
const double crLHS209 = crLHS101 + crLHS99;
const double crLHS210 = crLHS209*rho;
const double crLHS211 = crLHS110 + crLHS209*crLHS24;
const double crLHS212 = crLHS32*crLHS96;
const double crLHS213 = crLHS212*crLHS55;
const double crLHS214 = crLHS209*crLHS50;
const double crLHS215 = crLHS214*crLHS48;
const double crLHS216 = crLHS187 + crLHS209*crLHS64;
const double crLHS217 = pow(DN(2,0), 2);
const double crLHS218 = pow(N[2], 2);
const double crLHS219 = crLHS218*rho;
const double crLHS220 = crLHS218*crLHS4;
const double crLHS221 = bdf0*crLHS219 + crLHS209*crLHS96 + crLHS220;
const double crLHS222 = DN(2,1)*crLHS200;
const double crLHS223 = DN(2,0)*N[2];
const double crLHS224 = crLHS127*crLHS212;
const double crLHS225 = crLHS125*crLHS214;
const double crLHS226 = pow(DN(2,1), 2);
const double crLHS227 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS20*crLHS3 + crLHS22*crLHS9 + crLHS23*crLHS34 + crLHS34*crLHS36 - crLHS39*crLHS40 + crLHS42);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS45 - N[0]*crLHS52 + crLHS13*crLHS22 + crLHS47 - crLHS48*crLHS49 - crLHS54*crLHS57);
rLHS(0,2)+=-gauss_weight*(crLHS23*crLHS59 + crLHS36*crLHS59 + crLHS40*crLHS60 + crLHS58);
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + crLHS23*crLHS73 + crLHS36*crLHS73 - crLHS40*crLHS75 + crLHS79 + crLHS81);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS82 + DN(0,1)*crLHS84 - N[1]*crLHS52 - crLHS57*crLHS87 + crLHS85 + crLHS89);
rLHS(0,5)+=-gauss_weight*(crLHS23*crLHS91 + crLHS36*crLHS91 + crLHS40*crLHS92 + crLHS90);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS93 + DN(0,1)*crLHS95 + crLHS105*crLHS23 + crLHS105*crLHS36 - crLHS107*crLHS40 + crLHS111 + crLHS113);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS116 - N[2]*crLHS52 + crLHS117 - crLHS119*crLHS57 + crLHS120);
rLHS(0,8)+=-gauss_weight*(crLHS121 + crLHS122*crLHS23 + crLHS122*crLHS36 + crLHS123*crLHS40);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS124 - N[0]*crLHS126 - crLHS125*crLHS49 - crLHS128*crLHS54 + crLHS17*crLHS22 + crLHS47);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS129 + crLHS130*crLHS20 + crLHS132*crLHS23 + crLHS132*crLHS36 - crLHS133*crLHS40 + crLHS18*crLHS22 + crLHS42);
rLHS(1,2)+=-gauss_weight*(crLHS134 + crLHS135*crLHS23 + crLHS135*crLHS36 + crLHS136*crLHS40);
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS137 - N[1]*crLHS126 - crLHS128*crLHS87 + crLHS139 + crLHS141);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS84 + DN(0,1)*crLHS142 + crLHS144*crLHS23 + crLHS144*crLHS36 - crLHS145*crLHS40 + crLHS147 + crLHS79);
rLHS(1,5)+=-gauss_weight*(crLHS148 + crLHS149*crLHS23 + crLHS149*crLHS36 + crLHS150*crLHS40);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS95 + DN(0,1)*crLHS151 - N[2]*crLHS126 - crLHS119*crLHS128 + crLHS152 + crLHS153);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS116 + DN(0,1)*crLHS154 + crLHS111 + crLHS156*crLHS23 + crLHS156*crLHS36 - crLHS157*crLHS40 + crLHS159);
rLHS(1,8)+=-gauss_weight*(crLHS160 + crLHS161*crLHS23 + crLHS161*crLHS36 + crLHS162*crLHS40);
rLHS(2,0)+=gauss_weight*(crLHS134*crLHS140 - crLHS31*crLHS59 + crLHS58);
rLHS(2,1)+=gauss_weight*(-crLHS131*crLHS135 + crLHS134 + crLHS58*crLHS88);
rLHS(2,2)+=crLHS163*(crLHS130 + crLHS3);
rLHS(2,3)+=gauss_weight*(crLHS140*crLHS148 + crLHS164 - crLHS59*crLHS72);
rLHS(2,4)+=gauss_weight*(-crLHS135*crLHS143 + crLHS165 + crLHS88*crLHS90);
rLHS(2,5)+=crLHS166;
rLHS(2,6)+=gauss_weight*(-crLHS104*crLHS59 + crLHS140*crLHS160 + crLHS167);
rLHS(2,7)+=gauss_weight*(crLHS121*crLHS88 - crLHS135*crLHS155 + crLHS168);
rLHS(2,8)+=crLHS169;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS171*crLHS34 - crLHS172*crLHS39 + crLHS173 + crLHS34*crLHS66 + crLHS81);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS45 - N[0]*crLHS177 + crLHS139 - crLHS175*crLHS54 + crLHS89);
rLHS(3,2)+=-gauss_weight*(crLHS164 + crLHS171*crLHS59 + crLHS172*crLHS60 + crLHS59*crLHS66);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + crLHS171*crLHS73 - crLHS172*crLHS75 + crLHS178*crLHS20 + crLHS180*crLHS9 + crLHS182 + crLHS66*crLHS73);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS82 + DN(1,1)*crLHS84 - N[1]*crLHS177 + crLHS13*crLHS180 - crLHS175*crLHS87 - crLHS181*crLHS88 + crLHS184);
rLHS(3,5)+=-gauss_weight*(crLHS171*crLHS91 + crLHS172*crLHS92 + crLHS185 + crLHS66*crLHS91);
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS93 + DN(1,1)*crLHS95 + crLHS105*crLHS171 + crLHS105*crLHS66 - crLHS107*crLHS172 + crLHS188 + crLHS190);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS116 - N[2]*crLHS177 - crLHS119*crLHS175 + crLHS191 + crLHS194);
rLHS(3,8)+=-gauss_weight*(crLHS122*crLHS171 + crLHS122*crLHS66 + crLHS123*crLHS172 + crLHS195);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS124 - N[0]*crLHS197 + crLHS141 - crLHS196*crLHS54 + crLHS85);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS129 + crLHS132*crLHS171 + crLHS132*crLHS66 - crLHS133*crLHS172 + crLHS147 + crLHS173);
rLHS(4,2)+=-gauss_weight*(crLHS135*crLHS171 + crLHS135*crLHS66 + crLHS136*crLHS172 + crLHS165);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS137 - N[1]*crLHS197 - crLHS140*crLHS181 + crLHS17*crLHS180 + crLHS184 - crLHS196*crLHS87);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS84 + DN(1,1)*crLHS142 + crLHS144*crLHS171 + crLHS144*crLHS66 - crLHS145*crLHS172 + crLHS18*crLHS180 + crLHS182 + crLHS198*crLHS20);
rLHS(4,5)+=-gauss_weight*(crLHS149*crLHS171 + crLHS149*crLHS66 + crLHS150*crLHS172 + crLHS199);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS95 + DN(1,1)*crLHS151 - N[2]*crLHS197 - crLHS119*crLHS196 + crLHS201 + crLHS202);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS116 + DN(1,1)*crLHS154 + crLHS156*crLHS171 + crLHS156*crLHS66 - crLHS157*crLHS172 + crLHS188 + crLHS204);
rLHS(4,8)+=-gauss_weight*(crLHS161*crLHS171 + crLHS161*crLHS66 + crLHS162*crLHS172 + crLHS205);
rLHS(5,0)+=gauss_weight*(crLHS140*crLHS165 - crLHS31*crLHS91 + crLHS90);
rLHS(5,1)+=gauss_weight*(-crLHS131*crLHS149 + crLHS148 + crLHS164*crLHS88);
rLHS(5,2)+=crLHS166;
rLHS(5,3)+=gauss_weight*(crLHS140*crLHS199 + crLHS185 - crLHS72*crLHS91);
rLHS(5,4)+=gauss_weight*(-crLHS143*crLHS149 + crLHS185*crLHS88 + crLHS199);
rLHS(5,5)+=crLHS163*(crLHS178 + crLHS198);
rLHS(5,6)+=gauss_weight*(-crLHS104*crLHS91 + crLHS140*crLHS205 + crLHS206);
rLHS(5,7)+=gauss_weight*(-crLHS149*crLHS155 + crLHS195*crLHS88 + crLHS207);
rLHS(5,8)+=crLHS208;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS113 - crLHS192*crLHS39 + crLHS210*crLHS34 + crLHS211 + crLHS34*crLHS98);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS45 - N[0]*crLHS215 + crLHS120 + crLHS152 - crLHS213*crLHS54);
rLHS(6,2)+=-gauss_weight*(crLHS167 + crLHS192*crLHS60 + crLHS210*crLHS59 + crLHS59*crLHS98);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + crLHS190 - crLHS192*crLHS75 + crLHS210*crLHS73 + crLHS216 + crLHS73*crLHS98);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS82 + DN(2,1)*crLHS84 - N[1]*crLHS215 + crLHS194 + crLHS201 - crLHS213*crLHS87);
rLHS(6,5)+=-gauss_weight*(crLHS192*crLHS92 + crLHS206 + crLHS210*crLHS91 + crLHS91*crLHS98);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS93 + DN(2,1)*crLHS95 + crLHS105*crLHS210 + crLHS105*crLHS98 - crLHS107*crLHS192 + crLHS20*crLHS217 + crLHS219*crLHS9 + crLHS221);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS116 - N[2]*crLHS215 - crLHS119*crLHS213 + crLHS13*crLHS219 - crLHS220*crLHS88 + crLHS222);
rLHS(6,8)+=-gauss_weight*(crLHS122*crLHS210 + crLHS122*crLHS98 + crLHS123*crLHS192 + crLHS223);
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS124 - N[0]*crLHS225 + crLHS117 + crLHS153 - crLHS224*crLHS54);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS129 + crLHS132*crLHS210 + crLHS132*crLHS98 - crLHS133*crLHS192 + crLHS159 + crLHS211);
rLHS(7,2)+=-gauss_weight*(crLHS135*crLHS210 + crLHS135*crLHS98 + crLHS136*crLHS192 + crLHS168);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS137 - N[1]*crLHS225 + crLHS191 + crLHS202 - crLHS224*crLHS87);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS84 + DN(2,1)*crLHS142 + crLHS144*crLHS210 + crLHS144*crLHS98 - crLHS145*crLHS192 + crLHS204 + crLHS216);
rLHS(7,5)+=-gauss_weight*(crLHS149*crLHS210 + crLHS149*crLHS98 + crLHS150*crLHS192 + crLHS207);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS95 + DN(2,1)*crLHS151 - N[2]*crLHS225 - crLHS119*crLHS224 - crLHS140*crLHS220 + crLHS17*crLHS219 + crLHS222);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS116 + DN(2,1)*crLHS154 + crLHS156*crLHS210 + crLHS156*crLHS98 - crLHS157*crLHS192 + crLHS18*crLHS219 + crLHS20*crLHS226 + crLHS221);
rLHS(7,8)+=-gauss_weight*(crLHS161*crLHS210 + crLHS161*crLHS98 + crLHS162*crLHS192 + crLHS227);
rLHS(8,0)+=gauss_weight*(crLHS121 - crLHS122*crLHS31 + crLHS140*crLHS168);
rLHS(8,1)+=gauss_weight*(-crLHS131*crLHS161 + crLHS160 + crLHS167*crLHS88);
rLHS(8,2)+=crLHS169;
rLHS(8,3)+=gauss_weight*(-crLHS122*crLHS72 + crLHS140*crLHS207 + crLHS195);
rLHS(8,4)+=gauss_weight*(-crLHS143*crLHS161 + crLHS205 + crLHS206*crLHS88);
rLHS(8,5)+=crLHS208;
rLHS(8,6)+=gauss_weight*(-crLHS104*crLHS122 + crLHS140*crLHS227 + crLHS223);
rLHS(8,7)+=gauss_weight*(-crLHS155*crLHS161 + crLHS223*crLHS88 + crLHS227);
rLHS(8,8)+=crLHS163*(crLHS217 + crLHS226);
    
}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointLHSContributionAdjoint(
    FluidTopologyOptimizationElementData<3,4,true>& rData,
    MatrixType& rLHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,4>& c = rData.SoundVelocity;

    const array_1d<double,4> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;

    const double dyn_tau = rData.DynamicTau;
    
    // Get constitutive matrix
    const BoundedMatrix<double,6,6>& C = rData.C;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> v_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS1 = C(0,3)*DN(0,0);
const double crLHS2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS1;
const double crLHS3 = C(0,5)*DN(0,0);
const double crLHS4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS3;
const double crLHS5 = pow(DN(0,0), 2);
const double crLHS6 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crLHS7 = crLHS6*stab_c3;
const double crLHS8 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0) + N[3]*v_ns(3,0);
const double crLHS9 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crLHS10 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crLHS11 = rho*stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2));
const double crLHS12 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crLHS13 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crLHS14 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crLHS15 = DN(0,1)*v_ns(0,0) + DN(1,1)*v_ns(1,0) + DN(2,1)*v_ns(2,0) + DN(3,1)*v_ns(3,0);
const double crLHS16 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crLHS17 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crLHS18 = DN(0,2)*v_ns(0,0) + DN(1,2)*v_ns(1,0) + DN(2,2)*v_ns(2,0) + DN(3,2)*v_ns(3,0);
const double crLHS19 = DN(0,2)*v_ns(0,1) + DN(1,2)*v_ns(1,1) + DN(2,2)*v_ns(2,1) + DN(3,2)*v_ns(3,1);
const double crLHS20 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crLHS21 = rho*stab_c3*sqrt(pow(crLHS12, 2) + pow(crLHS13, 2) + pow(crLHS14, 2) + pow(crLHS15, 2) + pow(crLHS16, 2) + pow(crLHS17, 2) + pow(crLHS18, 2) + pow(crLHS19, 2) + pow(crLHS20, 2));
const double crLHS22 = h*(crLHS11 + crLHS21*h + crLHS7*h)/stab_c1 + mu;
const double crLHS23 = pow(N[0], 2);
const double crLHS24 = crLHS23*rho;
const double crLHS25 = N[0]*crLHS6;
const double crLHS26 = N[0]*rho;
const double crLHS27 = crLHS12*crLHS26;
const double crLHS28 = DN(0,0)*crLHS8;
const double crLHS29 = DN(0,1)*crLHS9;
const double crLHS30 = DN(0,2)*crLHS10;
const double crLHS31 = -crLHS25 + crLHS28*rho + crLHS29*rho + crLHS30*rho;
const double crLHS32 = -crLHS27 + crLHS31;
const double crLHS33 = 1.0/(crLHS11/h + crLHS21 + crLHS7 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crLHS34 = crLHS32*crLHS33;
const double crLHS35 = crLHS28 + crLHS29 + crLHS30;
const double crLHS36 = crLHS35*rho;
const double crLHS37 = crLHS13*crLHS26;
const double crLHS38 = crLHS15*crLHS37;
const double crLHS39 = crLHS14*crLHS26;
const double crLHS40 = crLHS18*crLHS39;
const double crLHS41 = -crLHS12*crLHS32 + crLHS38 + crLHS40;
const double crLHS42 = crLHS26*crLHS33;
const double crLHS43 = crLHS23*crLHS6;
const double crLHS44 = bdf0*crLHS24 + crLHS26*crLHS35 + crLHS43;
const double crLHS45 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS46 = C(1,3)*DN(0,1);
const double crLHS47 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS46;
const double crLHS48 = C(3,5)*DN(0,0);
const double crLHS49 = C(4,5)*DN(0,2);
const double crLHS50 = C(1,5)*DN(0,1) + crLHS48 + crLHS49;
const double crLHS51 = DN(0,0)*crLHS22;
const double crLHS52 = DN(0,1)*crLHS51;
const double crLHS53 = crLHS13*crLHS33;
const double crLHS54 = crLHS43*rho;
const double crLHS55 = pow(rho, 2);
const double crLHS56 = crLHS35*crLHS55;
const double crLHS57 = N[0]*crLHS56;
const double crLHS58 = crLHS16*crLHS26;
const double crLHS59 = crLHS31 - crLHS58;
const double crLHS60 = crLHS13*crLHS27 - crLHS13*crLHS59 + crLHS19*crLHS39;
const double crLHS61 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS62 = C(3,4)*DN(0,1);
const double crLHS63 = C(2,3)*DN(0,2) + crLHS48 + crLHS62;
const double crLHS64 = C(2,5)*DN(0,2);
const double crLHS65 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS64;
const double crLHS66 = DN(0,2)*crLHS51;
const double crLHS67 = crLHS14*crLHS33;
const double crLHS68 = crLHS20*crLHS26;
const double crLHS69 = crLHS31 - crLHS68;
const double crLHS70 = crLHS14*crLHS27 - crLHS14*crLHS69 + crLHS17*crLHS37;
const double crLHS71 = DN(0,0)*N[0];
const double crLHS72 = DN(0,0)*crLHS33;
const double crLHS73 = DN(0,0)*crLHS12 + DN(0,1)*crLHS13 + DN(0,2)*crLHS14;
const double crLHS74 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS75 = C(0,3)*DN(1,0);
const double crLHS76 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS75;
const double crLHS77 = C(0,5)*DN(1,0);
const double crLHS78 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS77;
const double crLHS79 = N[1]*rho;
const double crLHS80 = crLHS12*crLHS79;
const double crLHS81 = N[1]*crLHS6;
const double crLHS82 = DN(1,0)*crLHS8;
const double crLHS83 = DN(1,1)*crLHS9;
const double crLHS84 = DN(1,2)*crLHS10;
const double crLHS85 = -crLHS81 + crLHS82*rho + crLHS83*rho + crLHS84*rho;
const double crLHS86 = -crLHS80 + crLHS85;
const double crLHS87 = crLHS33*crLHS86;
const double crLHS88 = crLHS13*crLHS79;
const double crLHS89 = crLHS15*crLHS88;
const double crLHS90 = crLHS14*crLHS79;
const double crLHS91 = crLHS18*crLHS90;
const double crLHS92 = -crLHS12*crLHS86 + crLHS89 + crLHS91;
const double crLHS93 = N[1]*crLHS25;
const double crLHS94 = bdf0*crLHS26;
const double crLHS95 = N[1]*crLHS94 + crLHS93;
const double crLHS96 = crLHS35*crLHS79 + crLHS95;
const double crLHS97 = DN(0,0)*DN(1,0);
const double crLHS98 = N[1]*crLHS27 + crLHS22*crLHS97;
const double crLHS99 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS75;
const double crLHS100 = C(1,3)*DN(1,1);
const double crLHS101 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS100;
const double crLHS102 = C(3,5)*DN(1,0);
const double crLHS103 = C(4,5)*DN(1,2);
const double crLHS104 = C(1,5)*DN(1,1) + crLHS102 + crLHS103;
const double crLHS105 = DN(1,1)*crLHS51;
const double crLHS106 = crLHS16*crLHS79;
const double crLHS107 = -crLHS106 + crLHS85;
const double crLHS108 = -crLHS107*crLHS13 + crLHS13*crLHS80 + crLHS19*crLHS90;
const double crLHS109 = N[1]*crLHS56;
const double crLHS110 = crLHS93*rho;
const double crLHS111 = N[1]*crLHS37 - crLHS110*crLHS53;
const double crLHS112 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS77;
const double crLHS113 = C(3,4)*DN(1,1);
const double crLHS114 = C(2,3)*DN(1,2) + crLHS102 + crLHS113;
const double crLHS115 = C(2,5)*DN(1,2);
const double crLHS116 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS115;
const double crLHS117 = DN(1,2)*crLHS51;
const double crLHS118 = crLHS20*crLHS79;
const double crLHS119 = -crLHS118 + crLHS85;
const double crLHS120 = -crLHS119*crLHS14 + crLHS14*crLHS80 + crLHS17*crLHS88;
const double crLHS121 = N[1]*crLHS39 - crLHS110*crLHS67;
const double crLHS122 = DN(0,0)*N[1];
const double crLHS123 = DN(1,0)*crLHS33;
const double crLHS124 = DN(1,0)*crLHS12 + DN(1,1)*crLHS13 + DN(1,2)*crLHS14;
const double crLHS125 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS126 = C(0,3)*DN(2,0);
const double crLHS127 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS126;
const double crLHS128 = C(0,5)*DN(2,0);
const double crLHS129 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS128;
const double crLHS130 = N[2]*rho;
const double crLHS131 = crLHS12*crLHS130;
const double crLHS132 = N[2]*crLHS6;
const double crLHS133 = DN(2,0)*crLHS8;
const double crLHS134 = DN(2,1)*crLHS9;
const double crLHS135 = DN(2,2)*crLHS10;
const double crLHS136 = -crLHS132 + crLHS133*rho + crLHS134*rho + crLHS135*rho;
const double crLHS137 = -crLHS131 + crLHS136;
const double crLHS138 = crLHS137*crLHS33;
const double crLHS139 = crLHS13*crLHS130;
const double crLHS140 = crLHS139*crLHS15;
const double crLHS141 = crLHS130*crLHS14;
const double crLHS142 = crLHS141*crLHS18;
const double crLHS143 = -crLHS12*crLHS137 + crLHS140 + crLHS142;
const double crLHS144 = N[2]*crLHS25;
const double crLHS145 = N[2]*crLHS94 + crLHS144;
const double crLHS146 = crLHS130*crLHS35 + crLHS145;
const double crLHS147 = DN(0,0)*DN(2,0);
const double crLHS148 = N[2]*crLHS27 + crLHS147*crLHS22;
const double crLHS149 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS126;
const double crLHS150 = C(1,3)*DN(2,1);
const double crLHS151 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS150;
const double crLHS152 = C(3,5)*DN(2,0);
const double crLHS153 = C(4,5)*DN(2,2);
const double crLHS154 = C(1,5)*DN(2,1) + crLHS152 + crLHS153;
const double crLHS155 = DN(2,1)*crLHS51;
const double crLHS156 = crLHS130*crLHS16;
const double crLHS157 = crLHS136 - crLHS156;
const double crLHS158 = crLHS13*crLHS131 - crLHS13*crLHS157 + crLHS141*crLHS19;
const double crLHS159 = N[2]*crLHS56;
const double crLHS160 = crLHS144*rho;
const double crLHS161 = N[2]*crLHS37 - crLHS160*crLHS53;
const double crLHS162 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS128;
const double crLHS163 = C(3,4)*DN(2,1);
const double crLHS164 = C(2,3)*DN(2,2) + crLHS152 + crLHS163;
const double crLHS165 = C(2,5)*DN(2,2);
const double crLHS166 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS165;
const double crLHS167 = DN(2,2)*crLHS51;
const double crLHS168 = crLHS130*crLHS20;
const double crLHS169 = crLHS136 - crLHS168;
const double crLHS170 = crLHS131*crLHS14 + crLHS139*crLHS17 - crLHS14*crLHS169;
const double crLHS171 = N[2]*crLHS39 - crLHS160*crLHS67;
const double crLHS172 = DN(0,0)*N[2];
const double crLHS173 = DN(2,0)*crLHS33;
const double crLHS174 = DN(2,0)*crLHS12 + DN(2,1)*crLHS13 + DN(2,2)*crLHS14;
const double crLHS175 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS176 = C(0,3)*DN(3,0);
const double crLHS177 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS176;
const double crLHS178 = C(0,5)*DN(3,0);
const double crLHS179 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS178;
const double crLHS180 = N[3]*rho;
const double crLHS181 = crLHS12*crLHS180;
const double crLHS182 = N[3]*crLHS6;
const double crLHS183 = DN(3,0)*crLHS8;
const double crLHS184 = DN(3,1)*crLHS9;
const double crLHS185 = DN(3,2)*crLHS10;
const double crLHS186 = -crLHS182 + crLHS183*rho + crLHS184*rho + crLHS185*rho;
const double crLHS187 = -crLHS181 + crLHS186;
const double crLHS188 = crLHS187*crLHS33;
const double crLHS189 = crLHS13*crLHS180;
const double crLHS190 = crLHS15*crLHS189;
const double crLHS191 = crLHS14*crLHS180;
const double crLHS192 = crLHS18*crLHS191;
const double crLHS193 = -crLHS12*crLHS187 + crLHS190 + crLHS192;
const double crLHS194 = N[3]*crLHS25;
const double crLHS195 = N[3]*crLHS94 + crLHS194;
const double crLHS196 = crLHS180*crLHS35 + crLHS195;
const double crLHS197 = DN(0,0)*DN(3,0);
const double crLHS198 = N[3]*crLHS27 + crLHS197*crLHS22;
const double crLHS199 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS176;
const double crLHS200 = C(1,3)*DN(3,1);
const double crLHS201 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS200;
const double crLHS202 = C(3,5)*DN(3,0);
const double crLHS203 = C(4,5)*DN(3,2);
const double crLHS204 = C(1,5)*DN(3,1) + crLHS202 + crLHS203;
const double crLHS205 = DN(3,1)*crLHS51;
const double crLHS206 = crLHS16*crLHS180;
const double crLHS207 = crLHS186 - crLHS206;
const double crLHS208 = crLHS13*crLHS181 - crLHS13*crLHS207 + crLHS19*crLHS191;
const double crLHS209 = N[3]*crLHS56;
const double crLHS210 = crLHS194*rho;
const double crLHS211 = N[3]*crLHS37 - crLHS210*crLHS53;
const double crLHS212 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS178;
const double crLHS213 = C(3,4)*DN(3,1);
const double crLHS214 = C(2,3)*DN(3,2) + crLHS202 + crLHS213;
const double crLHS215 = C(2,5)*DN(3,2);
const double crLHS216 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS215;
const double crLHS217 = DN(3,2)*crLHS51;
const double crLHS218 = crLHS180*crLHS20;
const double crLHS219 = crLHS186 - crLHS218;
const double crLHS220 = crLHS14*crLHS181 - crLHS14*crLHS219 + crLHS17*crLHS189;
const double crLHS221 = N[3]*crLHS39 - crLHS210*crLHS67;
const double crLHS222 = DN(0,0)*N[3];
const double crLHS223 = DN(3,0)*crLHS33;
const double crLHS224 = DN(3,0)*crLHS12 + DN(3,1)*crLHS13 + DN(3,2)*crLHS14;
const double crLHS225 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS46;
const double crLHS226 = C(0,4)*DN(0,0) + crLHS49 + crLHS62;
const double crLHS227 = crLHS15*crLHS33;
const double crLHS228 = crLHS17*crLHS26;
const double crLHS229 = -crLHS15*crLHS32 + crLHS15*crLHS58 + crLHS18*crLHS228;
const double crLHS230 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS231 = C(1,4)*DN(0,1);
const double crLHS232 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS231;
const double crLHS233 = pow(DN(0,1), 2);
const double crLHS234 = crLHS33*crLHS59;
const double crLHS235 = crLHS19*crLHS228;
const double crLHS236 = -crLHS16*crLHS59 + crLHS235 + crLHS38;
const double crLHS237 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS231;
const double crLHS238 = C(2,4)*DN(0,2);
const double crLHS239 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS238;
const double crLHS240 = DN(0,1)*crLHS22;
const double crLHS241 = DN(0,2)*crLHS240;
const double crLHS242 = crLHS17*crLHS33;
const double crLHS243 = crLHS15*crLHS39 + crLHS17*crLHS58 - crLHS17*crLHS69;
const double crLHS244 = DN(0,1)*N[0];
const double crLHS245 = DN(0,1)*crLHS33;
const double crLHS246 = DN(0,0)*crLHS15 + DN(0,1)*crLHS16 + DN(0,2)*crLHS17;
const double crLHS247 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS100;
const double crLHS248 = C(0,4)*DN(1,0) + crLHS103 + crLHS113;
const double crLHS249 = DN(1,0)*crLHS240;
const double crLHS250 = crLHS17*crLHS79;
const double crLHS251 = crLHS106*crLHS15 - crLHS15*crLHS86 + crLHS18*crLHS250;
const double crLHS252 = crLHS15*crLHS26;
const double crLHS253 = N[1]*crLHS252 - crLHS110*crLHS227;
const double crLHS254 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS255 = C(1,4)*DN(1,1);
const double crLHS256 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS255;
const double crLHS257 = crLHS107*crLHS33;
const double crLHS258 = crLHS19*crLHS250;
const double crLHS259 = -crLHS107*crLHS16 + crLHS258 + crLHS89;
const double crLHS260 = DN(0,1)*DN(1,1);
const double crLHS261 = N[1]*crLHS58 + crLHS22*crLHS260;
const double crLHS262 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS255;
const double crLHS263 = C(2,4)*DN(1,2);
const double crLHS264 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS263;
const double crLHS265 = DN(1,2)*crLHS240;
const double crLHS266 = crLHS106*crLHS17 - crLHS119*crLHS17 + crLHS15*crLHS90;
const double crLHS267 = N[1]*crLHS228 - crLHS110*crLHS242;
const double crLHS268 = DN(0,1)*N[1];
const double crLHS269 = DN(1,1)*crLHS33;
const double crLHS270 = DN(1,0)*crLHS15 + DN(1,1)*crLHS16 + DN(1,2)*crLHS17;
const double crLHS271 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS150;
const double crLHS272 = C(0,4)*DN(2,0) + crLHS153 + crLHS163;
const double crLHS273 = DN(2,0)*crLHS240;
const double crLHS274 = crLHS130*crLHS17;
const double crLHS275 = -crLHS137*crLHS15 + crLHS15*crLHS156 + crLHS18*crLHS274;
const double crLHS276 = N[2]*crLHS252 - crLHS160*crLHS227;
const double crLHS277 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS278 = C(1,4)*DN(2,1);
const double crLHS279 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS278;
const double crLHS280 = crLHS157*crLHS33;
const double crLHS281 = crLHS19*crLHS274;
const double crLHS282 = crLHS140 - crLHS157*crLHS16 + crLHS281;
const double crLHS283 = DN(0,1)*DN(2,1);
const double crLHS284 = N[2]*crLHS58 + crLHS22*crLHS283;
const double crLHS285 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS278;
const double crLHS286 = C(2,4)*DN(2,2);
const double crLHS287 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS286;
const double crLHS288 = DN(2,2)*crLHS240;
const double crLHS289 = crLHS141*crLHS15 + crLHS156*crLHS17 - crLHS169*crLHS17;
const double crLHS290 = N[2]*crLHS228 - crLHS160*crLHS242;
const double crLHS291 = DN(0,1)*N[2];
const double crLHS292 = DN(2,1)*crLHS33;
const double crLHS293 = DN(2,0)*crLHS15 + DN(2,1)*crLHS16 + DN(2,2)*crLHS17;
const double crLHS294 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS200;
const double crLHS295 = C(0,4)*DN(3,0) + crLHS203 + crLHS213;
const double crLHS296 = DN(3,0)*crLHS240;
const double crLHS297 = crLHS17*crLHS180;
const double crLHS298 = -crLHS15*crLHS187 + crLHS15*crLHS206 + crLHS18*crLHS297;
const double crLHS299 = N[3]*crLHS252 - crLHS210*crLHS227;
const double crLHS300 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS301 = C(1,4)*DN(3,1);
const double crLHS302 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS301;
const double crLHS303 = crLHS207*crLHS33;
const double crLHS304 = crLHS19*crLHS297;
const double crLHS305 = -crLHS16*crLHS207 + crLHS190 + crLHS304;
const double crLHS306 = DN(0,1)*DN(3,1);
const double crLHS307 = N[3]*crLHS58 + crLHS22*crLHS306;
const double crLHS308 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS301;
const double crLHS309 = C(2,4)*DN(3,2);
const double crLHS310 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS309;
const double crLHS311 = DN(3,2)*crLHS240;
const double crLHS312 = crLHS15*crLHS191 + crLHS17*crLHS206 - crLHS17*crLHS219;
const double crLHS313 = N[3]*crLHS228 - crLHS210*crLHS242;
const double crLHS314 = DN(0,1)*N[3];
const double crLHS315 = DN(3,1)*crLHS33;
const double crLHS316 = DN(3,0)*crLHS15 + DN(3,1)*crLHS16 + DN(3,2)*crLHS17;
const double crLHS317 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS64;
const double crLHS318 = crLHS18*crLHS33;
const double crLHS319 = -crLHS18*crLHS32 + crLHS18*crLHS68 + crLHS19*crLHS252;
const double crLHS320 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS238;
const double crLHS321 = crLHS19*crLHS33;
const double crLHS322 = crLHS18*crLHS37 - crLHS19*crLHS59 + crLHS19*crLHS68;
const double crLHS323 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS324 = pow(DN(0,2), 2);
const double crLHS325 = crLHS33*crLHS69;
const double crLHS326 = -crLHS20*crLHS69 + crLHS235 + crLHS40;
const double crLHS327 = DN(0,2)*N[0];
const double crLHS328 = DN(0,2)*crLHS33;
const double crLHS329 = DN(0,0)*crLHS18 + DN(0,1)*crLHS19 + DN(0,2)*crLHS20;
const double crLHS330 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS115;
const double crLHS331 = DN(0,2)*crLHS22;
const double crLHS332 = DN(1,0)*crLHS331;
const double crLHS333 = crLHS15*crLHS19;
const double crLHS334 = crLHS118*crLHS18 - crLHS18*crLHS86 + crLHS333*crLHS79;
const double crLHS335 = N[1]*crLHS26;
const double crLHS336 = -crLHS110*crLHS318 + crLHS18*crLHS335;
const double crLHS337 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS263;
const double crLHS338 = DN(1,1)*crLHS331;
const double crLHS339 = -crLHS107*crLHS19 + crLHS118*crLHS19 + crLHS18*crLHS88;
const double crLHS340 = -crLHS110*crLHS321 + crLHS19*crLHS335;
const double crLHS341 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS342 = crLHS119*crLHS33;
const double crLHS343 = -crLHS119*crLHS20 + crLHS258 + crLHS91;
const double crLHS344 = DN(0,2)*DN(1,2);
const double crLHS345 = N[1]*crLHS68 + crLHS22*crLHS344;
const double crLHS346 = DN(0,2)*N[1];
const double crLHS347 = DN(1,2)*crLHS33;
const double crLHS348 = DN(1,0)*crLHS18 + DN(1,1)*crLHS19 + DN(1,2)*crLHS20;
const double crLHS349 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS165;
const double crLHS350 = DN(2,0)*crLHS331;
const double crLHS351 = crLHS130*crLHS333 - crLHS137*crLHS18 + crLHS168*crLHS18;
const double crLHS352 = N[2]*crLHS26;
const double crLHS353 = -crLHS160*crLHS318 + crLHS18*crLHS352;
const double crLHS354 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS286;
const double crLHS355 = DN(2,1)*crLHS331;
const double crLHS356 = crLHS139*crLHS18 - crLHS157*crLHS19 + crLHS168*crLHS19;
const double crLHS357 = -crLHS160*crLHS321 + crLHS19*crLHS352;
const double crLHS358 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS359 = crLHS169*crLHS33;
const double crLHS360 = crLHS142 - crLHS169*crLHS20 + crLHS281;
const double crLHS361 = DN(0,2)*DN(2,2);
const double crLHS362 = N[2]*crLHS68 + crLHS22*crLHS361;
const double crLHS363 = DN(0,2)*N[2];
const double crLHS364 = DN(2,2)*crLHS33;
const double crLHS365 = DN(2,0)*crLHS18 + DN(2,1)*crLHS19 + DN(2,2)*crLHS20;
const double crLHS366 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS215;
const double crLHS367 = DN(3,0)*crLHS331;
const double crLHS368 = -crLHS18*crLHS187 + crLHS18*crLHS218 + crLHS180*crLHS333;
const double crLHS369 = N[3]*crLHS26;
const double crLHS370 = crLHS18*crLHS369 - crLHS210*crLHS318;
const double crLHS371 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS309;
const double crLHS372 = DN(3,1)*crLHS331;
const double crLHS373 = crLHS18*crLHS189 - crLHS19*crLHS207 + crLHS19*crLHS218;
const double crLHS374 = crLHS19*crLHS369 - crLHS210*crLHS321;
const double crLHS375 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS376 = crLHS219*crLHS33;
const double crLHS377 = crLHS192 - crLHS20*crLHS219 + crLHS304;
const double crLHS378 = DN(0,2)*DN(3,2);
const double crLHS379 = N[3]*crLHS68 + crLHS22*crLHS378;
const double crLHS380 = DN(0,2)*N[3];
const double crLHS381 = DN(3,2)*crLHS33;
const double crLHS382 = DN(3,0)*crLHS18 + DN(3,1)*crLHS19 + DN(3,2)*crLHS20;
const double crLHS383 = crLHS244*rho;
const double crLHS384 = crLHS327*rho;
const double crLHS385 = crLHS71*rho;
const double crLHS386 = crLHS33*gauss_weight;
const double crLHS387 = DN(1,0)*N[0];
const double crLHS388 = crLHS268*rho;
const double crLHS389 = crLHS346*rho;
const double crLHS390 = DN(1,1)*N[0];
const double crLHS391 = crLHS122*rho;
const double crLHS392 = DN(1,2)*N[0];
const double crLHS393 = crLHS386*(crLHS260 + crLHS344 + crLHS97);
const double crLHS394 = DN(2,0)*N[0];
const double crLHS395 = crLHS291*rho;
const double crLHS396 = crLHS363*rho;
const double crLHS397 = DN(2,1)*N[0];
const double crLHS398 = crLHS172*rho;
const double crLHS399 = DN(2,2)*N[0];
const double crLHS400 = crLHS386*(crLHS147 + crLHS283 + crLHS361);
const double crLHS401 = DN(3,0)*N[0];
const double crLHS402 = crLHS314*rho;
const double crLHS403 = crLHS380*rho;
const double crLHS404 = DN(3,1)*N[0];
const double crLHS405 = crLHS222*rho;
const double crLHS406 = DN(3,2)*N[0];
const double crLHS407 = crLHS386*(crLHS197 + crLHS306 + crLHS378);
const double crLHS408 = crLHS82 + crLHS83 + crLHS84;
const double crLHS409 = crLHS408*rho;
const double crLHS410 = crLHS33*crLHS79;
const double crLHS411 = crLHS26*crLHS408 + crLHS95;
const double crLHS412 = crLHS408*crLHS55;
const double crLHS413 = N[0]*crLHS412;
const double crLHS414 = pow(DN(1,0), 2);
const double crLHS415 = pow(N[1], 2);
const double crLHS416 = crLHS415*rho;
const double crLHS417 = crLHS415*crLHS6;
const double crLHS418 = bdf0*crLHS416 + crLHS408*crLHS79 + crLHS417;
const double crLHS419 = DN(1,0)*crLHS22;
const double crLHS420 = DN(1,1)*crLHS419;
const double crLHS421 = crLHS417*rho;
const double crLHS422 = N[1]*crLHS412;
const double crLHS423 = DN(1,2)*crLHS419;
const double crLHS424 = DN(1,0)*N[1];
const double crLHS425 = bdf0*crLHS79;
const double crLHS426 = N[2]*crLHS425 + N[2]*crLHS81;
const double crLHS427 = crLHS130*crLHS408 + crLHS426;
const double crLHS428 = DN(1,0)*DN(2,0);
const double crLHS429 = N[2]*crLHS80 + crLHS22*crLHS428;
const double crLHS430 = DN(2,1)*crLHS419;
const double crLHS431 = N[2]*crLHS412;
const double crLHS432 = crLHS33*crLHS81;
const double crLHS433 = N[2]*crLHS88 - crLHS139*crLHS432;
const double crLHS434 = DN(2,2)*crLHS419;
const double crLHS435 = N[2]*crLHS90 - crLHS141*crLHS432;
const double crLHS436 = DN(1,0)*N[2];
const double crLHS437 = N[3]*crLHS425 + N[3]*crLHS81;
const double crLHS438 = crLHS180*crLHS408 + crLHS437;
const double crLHS439 = DN(1,0)*DN(3,0);
const double crLHS440 = N[3]*crLHS80 + crLHS22*crLHS439;
const double crLHS441 = DN(3,1)*crLHS419;
const double crLHS442 = N[3]*crLHS412;
const double crLHS443 = N[3]*crLHS88 - crLHS189*crLHS432;
const double crLHS444 = DN(3,2)*crLHS419;
const double crLHS445 = N[3]*crLHS90 - crLHS191*crLHS432;
const double crLHS446 = DN(1,0)*N[3];
const double crLHS447 = pow(DN(1,1), 2);
const double crLHS448 = DN(1,1)*crLHS22;
const double crLHS449 = DN(1,2)*crLHS448;
const double crLHS450 = DN(1,1)*N[1];
const double crLHS451 = DN(2,0)*crLHS448;
const double crLHS452 = crLHS15*crLHS79;
const double crLHS453 = crLHS130*crLHS33;
const double crLHS454 = crLHS15*crLHS81;
const double crLHS455 = N[2]*crLHS452 - crLHS453*crLHS454;
const double crLHS456 = DN(1,1)*DN(2,1);
const double crLHS457 = N[2]*crLHS106 + crLHS22*crLHS456;
const double crLHS458 = DN(2,2)*crLHS448;
const double crLHS459 = N[2]*crLHS250 - crLHS274*crLHS432;
const double crLHS460 = DN(1,1)*N[2];
const double crLHS461 = DN(3,0)*crLHS448;
const double crLHS462 = crLHS180*crLHS33;
const double crLHS463 = N[3]*crLHS452 - crLHS454*crLHS462;
const double crLHS464 = DN(1,1)*DN(3,1);
const double crLHS465 = N[3]*crLHS106 + crLHS22*crLHS464;
const double crLHS466 = DN(3,2)*crLHS448;
const double crLHS467 = N[3]*crLHS250 - crLHS297*crLHS432;
const double crLHS468 = DN(1,1)*N[3];
const double crLHS469 = pow(DN(1,2), 2);
const double crLHS470 = DN(1,2)*N[1];
const double crLHS471 = DN(1,2)*crLHS22;
const double crLHS472 = DN(2,0)*crLHS471;
const double crLHS473 = N[2]*crLHS79;
const double crLHS474 = crLHS453*crLHS81;
const double crLHS475 = crLHS18*crLHS473 - crLHS18*crLHS474;
const double crLHS476 = DN(2,1)*crLHS471;
const double crLHS477 = crLHS19*crLHS473 - crLHS19*crLHS474;
const double crLHS478 = DN(1,2)*DN(2,2);
const double crLHS479 = N[2]*crLHS118 + crLHS22*crLHS478;
const double crLHS480 = DN(1,2)*N[2];
const double crLHS481 = DN(3,0)*crLHS471;
const double crLHS482 = N[3]*crLHS79;
const double crLHS483 = crLHS462*crLHS81;
const double crLHS484 = crLHS18*crLHS482 - crLHS18*crLHS483;
const double crLHS485 = DN(3,1)*crLHS471;
const double crLHS486 = crLHS19*crLHS482 - crLHS19*crLHS483;
const double crLHS487 = DN(1,2)*DN(3,2);
const double crLHS488 = N[3]*crLHS118 + crLHS22*crLHS487;
const double crLHS489 = DN(1,2)*N[3];
const double crLHS490 = crLHS390*rho;
const double crLHS491 = crLHS392*rho;
const double crLHS492 = crLHS387*rho;
const double crLHS493 = crLHS450*rho;
const double crLHS494 = crLHS470*rho;
const double crLHS495 = crLHS424*rho;
const double crLHS496 = DN(2,0)*N[1];
const double crLHS497 = crLHS460*rho;
const double crLHS498 = crLHS480*rho;
const double crLHS499 = DN(2,1)*N[1];
const double crLHS500 = crLHS436*rho;
const double crLHS501 = DN(2,2)*N[1];
const double crLHS502 = crLHS386*(crLHS428 + crLHS456 + crLHS478);
const double crLHS503 = DN(3,0)*N[1];
const double crLHS504 = crLHS468*rho;
const double crLHS505 = crLHS489*rho;
const double crLHS506 = DN(3,1)*N[1];
const double crLHS507 = crLHS446*rho;
const double crLHS508 = DN(3,2)*N[1];
const double crLHS509 = crLHS386*(crLHS439 + crLHS464 + crLHS487);
const double crLHS510 = crLHS133 + crLHS134 + crLHS135;
const double crLHS511 = crLHS510*rho;
const double crLHS512 = crLHS145 + crLHS26*crLHS510;
const double crLHS513 = crLHS510*crLHS55;
const double crLHS514 = N[0]*crLHS513;
const double crLHS515 = crLHS426 + crLHS510*crLHS79;
const double crLHS516 = N[1]*crLHS513;
const double crLHS517 = pow(DN(2,0), 2);
const double crLHS518 = pow(N[2], 2);
const double crLHS519 = crLHS518*rho;
const double crLHS520 = crLHS518*crLHS6;
const double crLHS521 = bdf0*crLHS519 + crLHS130*crLHS510 + crLHS520;
const double crLHS522 = DN(2,0)*crLHS22;
const double crLHS523 = DN(2,1)*crLHS522;
const double crLHS524 = crLHS520*rho;
const double crLHS525 = N[2]*crLHS513;
const double crLHS526 = DN(2,2)*crLHS522;
const double crLHS527 = DN(2,0)*N[2];
const double crLHS528 = N[3]*crLHS130;
const double crLHS529 = N[3]*crLHS132 + bdf0*crLHS528;
const double crLHS530 = crLHS180*crLHS510 + crLHS529;
const double crLHS531 = DN(2,0)*DN(3,0);
const double crLHS532 = N[3]*crLHS131 + crLHS22*crLHS531;
const double crLHS533 = DN(3,1)*crLHS522;
const double crLHS534 = N[3]*crLHS513;
const double crLHS535 = crLHS132*crLHS33;
const double crLHS536 = N[3]*crLHS139 - crLHS189*crLHS535;
const double crLHS537 = DN(3,2)*crLHS522;
const double crLHS538 = N[3]*crLHS141 - crLHS191*crLHS535;
const double crLHS539 = DN(2,0)*N[3];
const double crLHS540 = pow(DN(2,1), 2);
const double crLHS541 = DN(2,1)*crLHS22;
const double crLHS542 = DN(2,2)*crLHS541;
const double crLHS543 = DN(2,1)*N[2];
const double crLHS544 = DN(3,0)*crLHS541;
const double crLHS545 = crLHS132*crLHS462;
const double crLHS546 = crLHS15*crLHS528 - crLHS15*crLHS545;
const double crLHS547 = DN(2,1)*DN(3,1);
const double crLHS548 = N[3]*crLHS156 + crLHS22*crLHS547;
const double crLHS549 = DN(3,2)*crLHS541;
const double crLHS550 = N[3]*crLHS274 - crLHS297*crLHS535;
const double crLHS551 = DN(2,1)*N[3];
const double crLHS552 = pow(DN(2,2), 2);
const double crLHS553 = DN(2,2)*N[2];
const double crLHS554 = DN(2,2)*crLHS22;
const double crLHS555 = DN(3,0)*crLHS554;
const double crLHS556 = crLHS18*crLHS528 - crLHS18*crLHS545;
const double crLHS557 = DN(3,1)*crLHS554;
const double crLHS558 = crLHS19*crLHS528 - crLHS19*crLHS545;
const double crLHS559 = DN(2,2)*DN(3,2);
const double crLHS560 = N[3]*crLHS168 + crLHS22*crLHS559;
const double crLHS561 = DN(2,2)*N[3];
const double crLHS562 = crLHS397*rho;
const double crLHS563 = crLHS399*rho;
const double crLHS564 = crLHS394*rho;
const double crLHS565 = crLHS499*rho;
const double crLHS566 = crLHS501*rho;
const double crLHS567 = crLHS496*rho;
const double crLHS568 = crLHS543*rho;
const double crLHS569 = crLHS553*rho;
const double crLHS570 = crLHS527*rho;
const double crLHS571 = DN(3,0)*N[2];
const double crLHS572 = crLHS551*rho;
const double crLHS573 = crLHS561*rho;
const double crLHS574 = DN(3,1)*N[2];
const double crLHS575 = crLHS539*rho;
const double crLHS576 = DN(3,2)*N[2];
const double crLHS577 = crLHS386*(crLHS531 + crLHS547 + crLHS559);
const double crLHS578 = crLHS183 + crLHS184 + crLHS185;
const double crLHS579 = crLHS578*rho;
const double crLHS580 = crLHS195 + crLHS26*crLHS578;
const double crLHS581 = crLHS55*crLHS578;
const double crLHS582 = N[0]*crLHS581;
const double crLHS583 = crLHS437 + crLHS578*crLHS79;
const double crLHS584 = N[1]*crLHS581;
const double crLHS585 = crLHS130*crLHS578 + crLHS529;
const double crLHS586 = N[2]*crLHS581;
const double crLHS587 = pow(DN(3,0), 2);
const double crLHS588 = pow(N[3], 2);
const double crLHS589 = crLHS588*rho;
const double crLHS590 = crLHS588*crLHS6;
const double crLHS591 = bdf0*crLHS589 + crLHS180*crLHS578 + crLHS590;
const double crLHS592 = DN(3,0)*crLHS22;
const double crLHS593 = DN(3,1)*crLHS592;
const double crLHS594 = crLHS590*rho;
const double crLHS595 = N[3]*crLHS581;
const double crLHS596 = DN(3,2)*crLHS592;
const double crLHS597 = DN(3,0)*N[3];
const double crLHS598 = pow(DN(3,1), 2);
const double crLHS599 = DN(3,1)*DN(3,2)*crLHS22;
const double crLHS600 = DN(3,1)*N[3];
const double crLHS601 = pow(DN(3,2), 2);
const double crLHS602 = DN(3,2)*N[3];
const double crLHS603 = crLHS404*rho;
const double crLHS604 = crLHS406*rho;
const double crLHS605 = crLHS401*rho;
const double crLHS606 = crLHS506*rho;
const double crLHS607 = crLHS508*rho;
const double crLHS608 = crLHS503*rho;
const double crLHS609 = crLHS574*rho;
const double crLHS610 = crLHS576*rho;
const double crLHS611 = crLHS571*rho;
const double crLHS612 = crLHS600*rho;
const double crLHS613 = crLHS602*rho;
const double crLHS614 = crLHS597*rho;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS24 + crLHS22*crLHS5 + crLHS25*crLHS34 + crLHS34*crLHS36 - crLHS41*crLHS42 + crLHS44);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS47 + DN(0,2)*crLHS50 + crLHS13*crLHS24 - crLHS42*crLHS60 + crLHS52 - crLHS53*crLHS54 - crLHS53*crLHS57);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + DN(0,2)*crLHS65 + crLHS14*crLHS24 - crLHS42*crLHS70 - crLHS54*crLHS67 - crLHS57*crLHS67 + crLHS66);
rLHS(0,3)+=-gauss_weight*(crLHS25*crLHS72 + crLHS36*crLHS72 + crLHS42*crLHS73 + crLHS71);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS76 + DN(0,2)*crLHS78 + crLHS25*crLHS87 + crLHS36*crLHS87 - crLHS42*crLHS92 + crLHS96 + crLHS98);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS99 + DN(0,1)*crLHS101 + DN(0,2)*crLHS104 + crLHS105 - crLHS108*crLHS42 - crLHS109*crLHS53 + crLHS111);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS114 + DN(0,2)*crLHS116 - crLHS109*crLHS67 + crLHS117 - crLHS120*crLHS42 + crLHS121);
rLHS(0,7)+=-gauss_weight*(crLHS122 + crLHS123*crLHS25 + crLHS123*crLHS36 + crLHS124*crLHS42);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS125 + DN(0,1)*crLHS127 + DN(0,2)*crLHS129 + crLHS138*crLHS25 + crLHS138*crLHS36 - crLHS143*crLHS42 + crLHS146 + crLHS148);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS149 + DN(0,1)*crLHS151 + DN(0,2)*crLHS154 + crLHS155 - crLHS158*crLHS42 - crLHS159*crLHS53 + crLHS161);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS162 + DN(0,1)*crLHS164 + DN(0,2)*crLHS166 - crLHS159*crLHS67 + crLHS167 - crLHS170*crLHS42 + crLHS171);
rLHS(0,11)+=-gauss_weight*(crLHS172 + crLHS173*crLHS25 + crLHS173*crLHS36 + crLHS174*crLHS42);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS175 + DN(0,1)*crLHS177 + DN(0,2)*crLHS179 + crLHS188*crLHS25 + crLHS188*crLHS36 - crLHS193*crLHS42 + crLHS196 + crLHS198);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS199 + DN(0,1)*crLHS201 + DN(0,2)*crLHS204 + crLHS205 - crLHS208*crLHS42 - crLHS209*crLHS53 + crLHS211);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS212 + DN(0,1)*crLHS214 + DN(0,2)*crLHS216 - crLHS209*crLHS67 + crLHS217 - crLHS220*crLHS42 + crLHS221);
rLHS(0,15)+=-gauss_weight*(crLHS222 + crLHS223*crLHS25 + crLHS223*crLHS36 + crLHS224*crLHS42);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS225 + DN(0,2)*crLHS226 + crLHS15*crLHS24 - crLHS227*crLHS54 - crLHS227*crLHS57 - crLHS229*crLHS42 + crLHS52);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS47 + DN(0,1)*crLHS230 + DN(0,2)*crLHS232 + crLHS16*crLHS24 + crLHS22*crLHS233 + crLHS234*crLHS25 + crLHS234*crLHS36 - crLHS236*crLHS42 + crLHS44);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS237 + DN(0,2)*crLHS239 + crLHS17*crLHS24 + crLHS241 - crLHS242*crLHS54 - crLHS242*crLHS57 - crLHS243*crLHS42);
rLHS(1,3)+=-gauss_weight*(crLHS244 + crLHS245*crLHS25 + crLHS245*crLHS36 + crLHS246*crLHS42);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS76 + DN(0,1)*crLHS247 + DN(0,2)*crLHS248 - crLHS109*crLHS227 + crLHS249 - crLHS251*crLHS42 + crLHS253);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS101 + DN(0,1)*crLHS254 + DN(0,2)*crLHS256 + crLHS25*crLHS257 + crLHS257*crLHS36 - crLHS259*crLHS42 + crLHS261 + crLHS96);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS262 + DN(0,2)*crLHS264 - crLHS109*crLHS242 + crLHS265 - crLHS266*crLHS42 + crLHS267);
rLHS(1,7)+=-gauss_weight*(crLHS25*crLHS269 + crLHS268 + crLHS269*crLHS36 + crLHS270*crLHS42);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS271 + DN(0,2)*crLHS272 - crLHS159*crLHS227 + crLHS273 - crLHS275*crLHS42 + crLHS276);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS151 + DN(0,1)*crLHS277 + DN(0,2)*crLHS279 + crLHS146 + crLHS25*crLHS280 + crLHS280*crLHS36 - crLHS282*crLHS42 + crLHS284);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS164 + DN(0,1)*crLHS285 + DN(0,2)*crLHS287 - crLHS159*crLHS242 + crLHS288 - crLHS289*crLHS42 + crLHS290);
rLHS(1,11)+=-gauss_weight*(crLHS25*crLHS292 + crLHS291 + crLHS292*crLHS36 + crLHS293*crLHS42);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS177 + DN(0,1)*crLHS294 + DN(0,2)*crLHS295 - crLHS209*crLHS227 + crLHS296 - crLHS298*crLHS42 + crLHS299);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS201 + DN(0,1)*crLHS300 + DN(0,2)*crLHS302 + crLHS196 + crLHS25*crLHS303 + crLHS303*crLHS36 - crLHS305*crLHS42 + crLHS307);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS214 + DN(0,1)*crLHS308 + DN(0,2)*crLHS310 - crLHS209*crLHS242 + crLHS311 - crLHS312*crLHS42 + crLHS313);
rLHS(1,15)+=-gauss_weight*(crLHS25*crLHS315 + crLHS314 + crLHS315*crLHS36 + crLHS316*crLHS42);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS226 + DN(0,2)*crLHS317 + crLHS18*crLHS24 - crLHS318*crLHS54 - crLHS318*crLHS57 - crLHS319*crLHS42 + crLHS66);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS50 + DN(0,1)*crLHS232 + DN(0,2)*crLHS320 + crLHS19*crLHS24 + crLHS241 - crLHS321*crLHS54 - crLHS321*crLHS57 - crLHS322*crLHS42);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS65 + DN(0,1)*crLHS239 + DN(0,2)*crLHS323 + crLHS20*crLHS24 + crLHS22*crLHS324 + crLHS25*crLHS325 + crLHS325*crLHS36 - crLHS326*crLHS42 + crLHS44);
rLHS(2,3)+=-gauss_weight*(crLHS25*crLHS328 + crLHS327 + crLHS328*crLHS36 + crLHS329*crLHS42);
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS78 + DN(0,1)*crLHS248 + DN(0,2)*crLHS330 - crLHS109*crLHS318 + crLHS332 - crLHS334*crLHS42 + crLHS336);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS104 + DN(0,1)*crLHS256 + DN(0,2)*crLHS337 - crLHS109*crLHS321 + crLHS338 - crLHS339*crLHS42 + crLHS340);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS116 + DN(0,1)*crLHS264 + DN(0,2)*crLHS341 + crLHS25*crLHS342 + crLHS342*crLHS36 - crLHS343*crLHS42 + crLHS345 + crLHS96);
rLHS(2,7)+=-gauss_weight*(crLHS25*crLHS347 + crLHS346 + crLHS347*crLHS36 + crLHS348*crLHS42);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS129 + DN(0,1)*crLHS272 + DN(0,2)*crLHS349 - crLHS159*crLHS318 + crLHS350 - crLHS351*crLHS42 + crLHS353);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS154 + DN(0,1)*crLHS279 + DN(0,2)*crLHS354 - crLHS159*crLHS321 + crLHS355 - crLHS356*crLHS42 + crLHS357);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS166 + DN(0,1)*crLHS287 + DN(0,2)*crLHS358 + crLHS146 + crLHS25*crLHS359 + crLHS359*crLHS36 - crLHS360*crLHS42 + crLHS362);
rLHS(2,11)+=-gauss_weight*(crLHS25*crLHS364 + crLHS36*crLHS364 + crLHS363 + crLHS365*crLHS42);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS179 + DN(0,1)*crLHS295 + DN(0,2)*crLHS366 - crLHS209*crLHS318 + crLHS367 - crLHS368*crLHS42 + crLHS370);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS204 + DN(0,1)*crLHS302 + DN(0,2)*crLHS371 - crLHS209*crLHS321 + crLHS372 - crLHS373*crLHS42 + crLHS374);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS216 + DN(0,1)*crLHS310 + DN(0,2)*crLHS375 + crLHS196 + crLHS25*crLHS376 + crLHS36*crLHS376 - crLHS377*crLHS42 + crLHS379);
rLHS(2,15)+=-gauss_weight*(crLHS25*crLHS381 + crLHS36*crLHS381 + crLHS380 + crLHS382*crLHS42);
rLHS(3,0)+=gauss_weight*(crLHS227*crLHS383 + crLHS318*crLHS384 - crLHS32*crLHS72 + crLHS71);
rLHS(3,1)+=gauss_weight*(crLHS244 - crLHS245*crLHS59 + crLHS321*crLHS384 + crLHS385*crLHS53);
rLHS(3,2)+=gauss_weight*(crLHS242*crLHS383 + crLHS327 - crLHS328*crLHS69 + crLHS385*crLHS67);
rLHS(3,3)+=crLHS386*(crLHS233 + crLHS324 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS227*crLHS388 + crLHS318*crLHS389 + crLHS387 - crLHS72*crLHS86);
rLHS(3,5)+=gauss_weight*(-crLHS107*crLHS245 + crLHS321*crLHS389 + crLHS390 + crLHS391*crLHS53);
rLHS(3,6)+=gauss_weight*(-crLHS119*crLHS328 + crLHS242*crLHS388 + crLHS391*crLHS67 + crLHS392);
rLHS(3,7)+=crLHS393;
rLHS(3,8)+=gauss_weight*(-crLHS137*crLHS72 + crLHS227*crLHS395 + crLHS318*crLHS396 + crLHS394);
rLHS(3,9)+=gauss_weight*(-crLHS157*crLHS245 + crLHS321*crLHS396 + crLHS397 + crLHS398*crLHS53);
rLHS(3,10)+=gauss_weight*(-crLHS169*crLHS328 + crLHS242*crLHS395 + crLHS398*crLHS67 + crLHS399);
rLHS(3,11)+=crLHS400;
rLHS(3,12)+=gauss_weight*(-crLHS187*crLHS72 + crLHS227*crLHS402 + crLHS318*crLHS403 + crLHS401);
rLHS(3,13)+=gauss_weight*(-crLHS207*crLHS245 + crLHS321*crLHS403 + crLHS404 + crLHS405*crLHS53);
rLHS(3,14)+=gauss_weight*(-crLHS219*crLHS328 + crLHS242*crLHS402 + crLHS405*crLHS67 + crLHS406);
rLHS(3,15)+=crLHS407;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS34*crLHS409 + crLHS34*crLHS81 - crLHS41*crLHS410 + crLHS411 + crLHS98);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS47 + DN(1,2)*crLHS50 + crLHS111 + crLHS249 - crLHS410*crLHS60 - crLHS413*crLHS53);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + DN(1,2)*crLHS65 + crLHS121 + crLHS332 - crLHS410*crLHS70 - crLHS413*crLHS67);
rLHS(4,3)+=-gauss_weight*(crLHS387 + crLHS409*crLHS72 + crLHS410*crLHS73 + crLHS72*crLHS81);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS76 + DN(1,2)*crLHS78 + crLHS12*crLHS416 + crLHS22*crLHS414 + crLHS409*crLHS87 - crLHS410*crLHS92 + crLHS418 + crLHS81*crLHS87);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS99 + DN(1,1)*crLHS101 + DN(1,2)*crLHS104 - crLHS108*crLHS410 + crLHS13*crLHS416 + crLHS420 - crLHS421*crLHS53 - crLHS422*crLHS53);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS114 + DN(1,2)*crLHS116 - crLHS120*crLHS410 + crLHS14*crLHS416 - crLHS421*crLHS67 - crLHS422*crLHS67 + crLHS423);
rLHS(4,7)+=-gauss_weight*(crLHS123*crLHS409 + crLHS123*crLHS81 + crLHS124*crLHS410 + crLHS424);
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS125 + DN(1,1)*crLHS127 + DN(1,2)*crLHS129 + crLHS138*crLHS409 + crLHS138*crLHS81 - crLHS143*crLHS410 + crLHS427 + crLHS429);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS149 + DN(1,1)*crLHS151 + DN(1,2)*crLHS154 - crLHS158*crLHS410 + crLHS430 - crLHS431*crLHS53 + crLHS433);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS162 + DN(1,1)*crLHS164 + DN(1,2)*crLHS166 - crLHS170*crLHS410 - crLHS431*crLHS67 + crLHS434 + crLHS435);
rLHS(4,11)+=-gauss_weight*(crLHS173*crLHS409 + crLHS173*crLHS81 + crLHS174*crLHS410 + crLHS436);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS175 + DN(1,1)*crLHS177 + DN(1,2)*crLHS179 + crLHS188*crLHS409 + crLHS188*crLHS81 - crLHS193*crLHS410 + crLHS438 + crLHS440);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS199 + DN(1,1)*crLHS201 + DN(1,2)*crLHS204 - crLHS208*crLHS410 + crLHS441 - crLHS442*crLHS53 + crLHS443);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS212 + DN(1,1)*crLHS214 + DN(1,2)*crLHS216 - crLHS220*crLHS410 - crLHS442*crLHS67 + crLHS444 + crLHS445);
rLHS(4,15)+=-gauss_weight*(crLHS223*crLHS409 + crLHS223*crLHS81 + crLHS224*crLHS410 + crLHS446);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS225 + DN(1,2)*crLHS226 + crLHS105 - crLHS227*crLHS413 - crLHS229*crLHS410 + crLHS253);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS47 + DN(1,1)*crLHS230 + DN(1,2)*crLHS232 + crLHS234*crLHS409 + crLHS234*crLHS81 - crLHS236*crLHS410 + crLHS261 + crLHS411);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS237 + DN(1,2)*crLHS239 - crLHS242*crLHS413 - crLHS243*crLHS410 + crLHS267 + crLHS338);
rLHS(5,3)+=-gauss_weight*(crLHS245*crLHS409 + crLHS245*crLHS81 + crLHS246*crLHS410 + crLHS390);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS76 + DN(1,1)*crLHS247 + DN(1,2)*crLHS248 + crLHS15*crLHS416 - crLHS227*crLHS421 - crLHS227*crLHS422 - crLHS251*crLHS410 + crLHS420);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS101 + DN(1,1)*crLHS254 + DN(1,2)*crLHS256 + crLHS16*crLHS416 + crLHS22*crLHS447 + crLHS257*crLHS409 + crLHS257*crLHS81 - crLHS259*crLHS410 + crLHS418);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS262 + DN(1,2)*crLHS264 + crLHS17*crLHS416 - crLHS242*crLHS421 - crLHS242*crLHS422 - crLHS266*crLHS410 + crLHS449);
rLHS(5,7)+=-gauss_weight*(crLHS269*crLHS409 + crLHS269*crLHS81 + crLHS270*crLHS410 + crLHS450);
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS271 + DN(1,2)*crLHS272 - crLHS227*crLHS431 - crLHS275*crLHS410 + crLHS451 + crLHS455);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS151 + DN(1,1)*crLHS277 + DN(1,2)*crLHS279 + crLHS280*crLHS409 + crLHS280*crLHS81 - crLHS282*crLHS410 + crLHS427 + crLHS457);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS164 + DN(1,1)*crLHS285 + DN(1,2)*crLHS287 - crLHS242*crLHS431 - crLHS289*crLHS410 + crLHS458 + crLHS459);
rLHS(5,11)+=-gauss_weight*(crLHS292*crLHS409 + crLHS292*crLHS81 + crLHS293*crLHS410 + crLHS460);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS177 + DN(1,1)*crLHS294 + DN(1,2)*crLHS295 - crLHS227*crLHS442 - crLHS298*crLHS410 + crLHS461 + crLHS463);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS201 + DN(1,1)*crLHS300 + DN(1,2)*crLHS302 + crLHS303*crLHS409 + crLHS303*crLHS81 - crLHS305*crLHS410 + crLHS438 + crLHS465);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS214 + DN(1,1)*crLHS308 + DN(1,2)*crLHS310 - crLHS242*crLHS442 - crLHS312*crLHS410 + crLHS466 + crLHS467);
rLHS(5,15)+=-gauss_weight*(crLHS315*crLHS409 + crLHS315*crLHS81 + crLHS316*crLHS410 + crLHS468);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS226 + DN(1,2)*crLHS317 + crLHS117 - crLHS318*crLHS413 - crLHS319*crLHS410 + crLHS336);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS50 + DN(1,1)*crLHS232 + DN(1,2)*crLHS320 + crLHS265 - crLHS321*crLHS413 - crLHS322*crLHS410 + crLHS340);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS65 + DN(1,1)*crLHS239 + DN(1,2)*crLHS323 + crLHS325*crLHS409 + crLHS325*crLHS81 - crLHS326*crLHS410 + crLHS345 + crLHS411);
rLHS(6,3)+=-gauss_weight*(crLHS328*crLHS409 + crLHS328*crLHS81 + crLHS329*crLHS410 + crLHS392);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS78 + DN(1,1)*crLHS248 + DN(1,2)*crLHS330 + crLHS18*crLHS416 - crLHS318*crLHS421 - crLHS318*crLHS422 - crLHS334*crLHS410 + crLHS423);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS104 + DN(1,1)*crLHS256 + DN(1,2)*crLHS337 + crLHS19*crLHS416 - crLHS321*crLHS421 - crLHS321*crLHS422 - crLHS339*crLHS410 + crLHS449);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS116 + DN(1,1)*crLHS264 + DN(1,2)*crLHS341 + crLHS20*crLHS416 + crLHS22*crLHS469 + crLHS342*crLHS409 + crLHS342*crLHS81 - crLHS343*crLHS410 + crLHS418);
rLHS(6,7)+=-gauss_weight*(crLHS347*crLHS409 + crLHS347*crLHS81 + crLHS348*crLHS410 + crLHS470);
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS129 + DN(1,1)*crLHS272 + DN(1,2)*crLHS349 - crLHS318*crLHS431 - crLHS351*crLHS410 + crLHS472 + crLHS475);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS154 + DN(1,1)*crLHS279 + DN(1,2)*crLHS354 - crLHS321*crLHS431 - crLHS356*crLHS410 + crLHS476 + crLHS477);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS166 + DN(1,1)*crLHS287 + DN(1,2)*crLHS358 + crLHS359*crLHS409 + crLHS359*crLHS81 - crLHS360*crLHS410 + crLHS427 + crLHS479);
rLHS(6,11)+=-gauss_weight*(crLHS364*crLHS409 + crLHS364*crLHS81 + crLHS365*crLHS410 + crLHS480);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS179 + DN(1,1)*crLHS295 + DN(1,2)*crLHS366 - crLHS318*crLHS442 - crLHS368*crLHS410 + crLHS481 + crLHS484);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS204 + DN(1,1)*crLHS302 + DN(1,2)*crLHS371 - crLHS321*crLHS442 - crLHS373*crLHS410 + crLHS485 + crLHS486);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS216 + DN(1,1)*crLHS310 + DN(1,2)*crLHS375 + crLHS376*crLHS409 + crLHS376*crLHS81 - crLHS377*crLHS410 + crLHS438 + crLHS488);
rLHS(6,15)+=-gauss_weight*(crLHS381*crLHS409 + crLHS381*crLHS81 + crLHS382*crLHS410 + crLHS489);
rLHS(7,0)+=gauss_weight*(crLHS122 - crLHS123*crLHS32 + crLHS227*crLHS490 + crLHS318*crLHS491);
rLHS(7,1)+=gauss_weight*(crLHS268 - crLHS269*crLHS59 + crLHS321*crLHS491 + crLHS492*crLHS53);
rLHS(7,2)+=gauss_weight*(crLHS242*crLHS490 + crLHS346 - crLHS347*crLHS69 + crLHS492*crLHS67);
rLHS(7,3)+=crLHS393;
rLHS(7,4)+=gauss_weight*(-crLHS123*crLHS86 + crLHS227*crLHS493 + crLHS318*crLHS494 + crLHS424);
rLHS(7,5)+=gauss_weight*(-crLHS107*crLHS269 + crLHS321*crLHS494 + crLHS450 + crLHS495*crLHS53);
rLHS(7,6)+=gauss_weight*(-crLHS119*crLHS347 + crLHS242*crLHS493 + crLHS470 + crLHS495*crLHS67);
rLHS(7,7)+=crLHS386*(crLHS414 + crLHS447 + crLHS469);
rLHS(7,8)+=gauss_weight*(-crLHS123*crLHS137 + crLHS227*crLHS497 + crLHS318*crLHS498 + crLHS496);
rLHS(7,9)+=gauss_weight*(-crLHS157*crLHS269 + crLHS321*crLHS498 + crLHS499 + crLHS500*crLHS53);
rLHS(7,10)+=gauss_weight*(-crLHS169*crLHS347 + crLHS242*crLHS497 + crLHS500*crLHS67 + crLHS501);
rLHS(7,11)+=crLHS502;
rLHS(7,12)+=gauss_weight*(-crLHS123*crLHS187 + crLHS227*crLHS504 + crLHS318*crLHS505 + crLHS503);
rLHS(7,13)+=gauss_weight*(-crLHS207*crLHS269 + crLHS321*crLHS505 + crLHS506 + crLHS507*crLHS53);
rLHS(7,14)+=gauss_weight*(-crLHS219*crLHS347 + crLHS242*crLHS504 + crLHS507*crLHS67 + crLHS508);
rLHS(7,15)+=crLHS509;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS132*crLHS34 + crLHS148 + crLHS34*crLHS511 - crLHS41*crLHS453 + crLHS512);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS47 + DN(2,2)*crLHS50 + crLHS161 + crLHS273 - crLHS453*crLHS60 - crLHS514*crLHS53);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + DN(2,2)*crLHS65 + crLHS171 + crLHS350 - crLHS453*crLHS70 - crLHS514*crLHS67);
rLHS(8,3)+=-gauss_weight*(crLHS132*crLHS72 + crLHS394 + crLHS453*crLHS73 + crLHS511*crLHS72);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS76 + DN(2,2)*crLHS78 + crLHS132*crLHS87 + crLHS429 - crLHS453*crLHS92 + crLHS511*crLHS87 + crLHS515);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS99 + DN(2,1)*crLHS101 + DN(2,2)*crLHS104 - crLHS108*crLHS453 + crLHS433 + crLHS451 - crLHS516*crLHS53);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS114 + DN(2,2)*crLHS116 - crLHS120*crLHS453 + crLHS435 + crLHS472 - crLHS516*crLHS67);
rLHS(8,7)+=-gauss_weight*(crLHS123*crLHS132 + crLHS123*crLHS511 + crLHS124*crLHS453 + crLHS496);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS125 + DN(2,1)*crLHS127 + DN(2,2)*crLHS129 + crLHS12*crLHS519 + crLHS132*crLHS138 + crLHS138*crLHS511 - crLHS143*crLHS453 + crLHS22*crLHS517 + crLHS521);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS149 + DN(2,1)*crLHS151 + DN(2,2)*crLHS154 + crLHS13*crLHS519 - crLHS158*crLHS453 + crLHS523 - crLHS524*crLHS53 - crLHS525*crLHS53);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS162 + DN(2,1)*crLHS164 + DN(2,2)*crLHS166 + crLHS14*crLHS519 - crLHS170*crLHS453 - crLHS524*crLHS67 - crLHS525*crLHS67 + crLHS526);
rLHS(8,11)+=-gauss_weight*(crLHS132*crLHS173 + crLHS173*crLHS511 + crLHS174*crLHS453 + crLHS527);
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS175 + DN(2,1)*crLHS177 + DN(2,2)*crLHS179 + crLHS132*crLHS188 + crLHS188*crLHS511 - crLHS193*crLHS453 + crLHS530 + crLHS532);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS199 + DN(2,1)*crLHS201 + DN(2,2)*crLHS204 - crLHS208*crLHS453 - crLHS53*crLHS534 + crLHS533 + crLHS536);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS212 + DN(2,1)*crLHS214 + DN(2,2)*crLHS216 - crLHS220*crLHS453 - crLHS534*crLHS67 + crLHS537 + crLHS538);
rLHS(8,15)+=-gauss_weight*(crLHS132*crLHS223 + crLHS223*crLHS511 + crLHS224*crLHS453 + crLHS539);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS225 + DN(2,2)*crLHS226 + crLHS155 - crLHS227*crLHS514 - crLHS229*crLHS453 + crLHS276);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS47 + DN(2,1)*crLHS230 + DN(2,2)*crLHS232 + crLHS132*crLHS234 + crLHS234*crLHS511 - crLHS236*crLHS453 + crLHS284 + crLHS512);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS237 + DN(2,2)*crLHS239 - crLHS242*crLHS514 - crLHS243*crLHS453 + crLHS290 + crLHS355);
rLHS(9,3)+=-gauss_weight*(crLHS132*crLHS245 + crLHS245*crLHS511 + crLHS246*crLHS453 + crLHS397);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS76 + DN(2,1)*crLHS247 + DN(2,2)*crLHS248 - crLHS227*crLHS516 - crLHS251*crLHS453 + crLHS430 + crLHS455);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS101 + DN(2,1)*crLHS254 + DN(2,2)*crLHS256 + crLHS132*crLHS257 + crLHS257*crLHS511 - crLHS259*crLHS453 + crLHS457 + crLHS515);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS262 + DN(2,2)*crLHS264 - crLHS242*crLHS516 - crLHS266*crLHS453 + crLHS459 + crLHS476);
rLHS(9,7)+=-gauss_weight*(crLHS132*crLHS269 + crLHS269*crLHS511 + crLHS270*crLHS453 + crLHS499);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS271 + DN(2,2)*crLHS272 + crLHS15*crLHS519 - crLHS227*crLHS524 - crLHS227*crLHS525 - crLHS275*crLHS453 + crLHS523);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS151 + DN(2,1)*crLHS277 + DN(2,2)*crLHS279 + crLHS132*crLHS280 + crLHS16*crLHS519 + crLHS22*crLHS540 + crLHS280*crLHS511 - crLHS282*crLHS453 + crLHS521);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS164 + DN(2,1)*crLHS285 + DN(2,2)*crLHS287 + crLHS17*crLHS519 - crLHS242*crLHS524 - crLHS242*crLHS525 - crLHS289*crLHS453 + crLHS542);
rLHS(9,11)+=-gauss_weight*(crLHS132*crLHS292 + crLHS292*crLHS511 + crLHS293*crLHS453 + crLHS543);
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS177 + DN(2,1)*crLHS294 + DN(2,2)*crLHS295 - crLHS227*crLHS534 - crLHS298*crLHS453 + crLHS544 + crLHS546);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS201 + DN(2,1)*crLHS300 + DN(2,2)*crLHS302 + crLHS132*crLHS303 + crLHS303*crLHS511 - crLHS305*crLHS453 + crLHS530 + crLHS548);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS214 + DN(2,1)*crLHS308 + DN(2,2)*crLHS310 - crLHS242*crLHS534 - crLHS312*crLHS453 + crLHS549 + crLHS550);
rLHS(9,15)+=-gauss_weight*(crLHS132*crLHS315 + crLHS315*crLHS511 + crLHS316*crLHS453 + crLHS551);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS226 + DN(2,2)*crLHS317 + crLHS167 - crLHS318*crLHS514 - crLHS319*crLHS453 + crLHS353);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS50 + DN(2,1)*crLHS232 + DN(2,2)*crLHS320 + crLHS288 - crLHS321*crLHS514 - crLHS322*crLHS453 + crLHS357);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS65 + DN(2,1)*crLHS239 + DN(2,2)*crLHS323 + crLHS132*crLHS325 + crLHS325*crLHS511 - crLHS326*crLHS453 + crLHS362 + crLHS512);
rLHS(10,3)+=-gauss_weight*(crLHS132*crLHS328 + crLHS328*crLHS511 + crLHS329*crLHS453 + crLHS399);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS78 + DN(2,1)*crLHS248 + DN(2,2)*crLHS330 - crLHS318*crLHS516 - crLHS334*crLHS453 + crLHS434 + crLHS475);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS104 + DN(2,1)*crLHS256 + DN(2,2)*crLHS337 - crLHS321*crLHS516 - crLHS339*crLHS453 + crLHS458 + crLHS477);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS116 + DN(2,1)*crLHS264 + DN(2,2)*crLHS341 + crLHS132*crLHS342 + crLHS342*crLHS511 - crLHS343*crLHS453 + crLHS479 + crLHS515);
rLHS(10,7)+=-gauss_weight*(crLHS132*crLHS347 + crLHS347*crLHS511 + crLHS348*crLHS453 + crLHS501);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS129 + DN(2,1)*crLHS272 + DN(2,2)*crLHS349 + crLHS18*crLHS519 - crLHS318*crLHS524 - crLHS318*crLHS525 - crLHS351*crLHS453 + crLHS526);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS154 + DN(2,1)*crLHS279 + DN(2,2)*crLHS354 + crLHS19*crLHS519 - crLHS321*crLHS524 - crLHS321*crLHS525 - crLHS356*crLHS453 + crLHS542);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS166 + DN(2,1)*crLHS287 + DN(2,2)*crLHS358 + crLHS132*crLHS359 + crLHS20*crLHS519 + crLHS22*crLHS552 + crLHS359*crLHS511 - crLHS360*crLHS453 + crLHS521);
rLHS(10,11)+=-gauss_weight*(crLHS132*crLHS364 + crLHS364*crLHS511 + crLHS365*crLHS453 + crLHS553);
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS179 + DN(2,1)*crLHS295 + DN(2,2)*crLHS366 - crLHS318*crLHS534 - crLHS368*crLHS453 + crLHS555 + crLHS556);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS204 + DN(2,1)*crLHS302 + DN(2,2)*crLHS371 - crLHS321*crLHS534 - crLHS373*crLHS453 + crLHS557 + crLHS558);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS216 + DN(2,1)*crLHS310 + DN(2,2)*crLHS375 + crLHS132*crLHS376 + crLHS376*crLHS511 - crLHS377*crLHS453 + crLHS530 + crLHS560);
rLHS(10,15)+=-gauss_weight*(crLHS132*crLHS381 + crLHS381*crLHS511 + crLHS382*crLHS453 + crLHS561);
rLHS(11,0)+=gauss_weight*(crLHS172 - crLHS173*crLHS32 + crLHS227*crLHS562 + crLHS318*crLHS563);
rLHS(11,1)+=gauss_weight*(crLHS291 - crLHS292*crLHS59 + crLHS321*crLHS563 + crLHS53*crLHS564);
rLHS(11,2)+=gauss_weight*(crLHS242*crLHS562 + crLHS363 - crLHS364*crLHS69 + crLHS564*crLHS67);
rLHS(11,3)+=crLHS400;
rLHS(11,4)+=gauss_weight*(-crLHS173*crLHS86 + crLHS227*crLHS565 + crLHS318*crLHS566 + crLHS436);
rLHS(11,5)+=gauss_weight*(-crLHS107*crLHS292 + crLHS321*crLHS566 + crLHS460 + crLHS53*crLHS567);
rLHS(11,6)+=gauss_weight*(-crLHS119*crLHS364 + crLHS242*crLHS565 + crLHS480 + crLHS567*crLHS67);
rLHS(11,7)+=crLHS502;
rLHS(11,8)+=gauss_weight*(-crLHS137*crLHS173 + crLHS227*crLHS568 + crLHS318*crLHS569 + crLHS527);
rLHS(11,9)+=gauss_weight*(-crLHS157*crLHS292 + crLHS321*crLHS569 + crLHS53*crLHS570 + crLHS543);
rLHS(11,10)+=gauss_weight*(-crLHS169*crLHS364 + crLHS242*crLHS568 + crLHS553 + crLHS570*crLHS67);
rLHS(11,11)+=crLHS386*(crLHS517 + crLHS540 + crLHS552);
rLHS(11,12)+=gauss_weight*(-crLHS173*crLHS187 + crLHS227*crLHS572 + crLHS318*crLHS573 + crLHS571);
rLHS(11,13)+=gauss_weight*(-crLHS207*crLHS292 + crLHS321*crLHS573 + crLHS53*crLHS575 + crLHS574);
rLHS(11,14)+=gauss_weight*(-crLHS219*crLHS364 + crLHS242*crLHS572 + crLHS575*crLHS67 + crLHS576);
rLHS(11,15)+=crLHS577;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS182*crLHS34 + crLHS198 + crLHS34*crLHS579 - crLHS41*crLHS462 + crLHS580);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS45 + DN(3,1)*crLHS47 + DN(3,2)*crLHS50 + crLHS211 + crLHS296 - crLHS462*crLHS60 - crLHS53*crLHS582);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS61 + DN(3,1)*crLHS63 + DN(3,2)*crLHS65 + crLHS221 + crLHS367 - crLHS462*crLHS70 - crLHS582*crLHS67);
rLHS(12,3)+=-gauss_weight*(crLHS182*crLHS72 + crLHS401 + crLHS462*crLHS73 + crLHS579*crLHS72);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS74 + DN(3,1)*crLHS76 + DN(3,2)*crLHS78 + crLHS182*crLHS87 + crLHS440 - crLHS462*crLHS92 + crLHS579*crLHS87 + crLHS583);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS99 + DN(3,1)*crLHS101 + DN(3,2)*crLHS104 - crLHS108*crLHS462 + crLHS443 + crLHS461 - crLHS53*crLHS584);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS112 + DN(3,1)*crLHS114 + DN(3,2)*crLHS116 - crLHS120*crLHS462 + crLHS445 + crLHS481 - crLHS584*crLHS67);
rLHS(12,7)+=-gauss_weight*(crLHS123*crLHS182 + crLHS123*crLHS579 + crLHS124*crLHS462 + crLHS503);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS125 + DN(3,1)*crLHS127 + DN(3,2)*crLHS129 + crLHS138*crLHS182 + crLHS138*crLHS579 - crLHS143*crLHS462 + crLHS532 + crLHS585);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS149 + DN(3,1)*crLHS151 + DN(3,2)*crLHS154 - crLHS158*crLHS462 - crLHS53*crLHS586 + crLHS536 + crLHS544);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS162 + DN(3,1)*crLHS164 + DN(3,2)*crLHS166 - crLHS170*crLHS462 + crLHS538 + crLHS555 - crLHS586*crLHS67);
rLHS(12,11)+=-gauss_weight*(crLHS173*crLHS182 + crLHS173*crLHS579 + crLHS174*crLHS462 + crLHS571);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS175 + DN(3,1)*crLHS177 + DN(3,2)*crLHS179 + crLHS12*crLHS589 + crLHS182*crLHS188 + crLHS188*crLHS579 - crLHS193*crLHS462 + crLHS22*crLHS587 + crLHS591);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS199 + DN(3,1)*crLHS201 + DN(3,2)*crLHS204 + crLHS13*crLHS589 - crLHS208*crLHS462 - crLHS53*crLHS594 - crLHS53*crLHS595 + crLHS593);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS212 + DN(3,1)*crLHS214 + DN(3,2)*crLHS216 + crLHS14*crLHS589 - crLHS220*crLHS462 - crLHS594*crLHS67 - crLHS595*crLHS67 + crLHS596);
rLHS(12,15)+=-gauss_weight*(crLHS182*crLHS223 + crLHS223*crLHS579 + crLHS224*crLHS462 + crLHS597);
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS225 + DN(3,2)*crLHS226 + crLHS205 - crLHS227*crLHS582 - crLHS229*crLHS462 + crLHS299);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS47 + DN(3,1)*crLHS230 + DN(3,2)*crLHS232 + crLHS182*crLHS234 + crLHS234*crLHS579 - crLHS236*crLHS462 + crLHS307 + crLHS580);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS63 + DN(3,1)*crLHS237 + DN(3,2)*crLHS239 - crLHS242*crLHS582 - crLHS243*crLHS462 + crLHS313 + crLHS372);
rLHS(13,3)+=-gauss_weight*(crLHS182*crLHS245 + crLHS245*crLHS579 + crLHS246*crLHS462 + crLHS404);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS76 + DN(3,1)*crLHS247 + DN(3,2)*crLHS248 - crLHS227*crLHS584 - crLHS251*crLHS462 + crLHS441 + crLHS463);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS101 + DN(3,1)*crLHS254 + DN(3,2)*crLHS256 + crLHS182*crLHS257 + crLHS257*crLHS579 - crLHS259*crLHS462 + crLHS465 + crLHS583);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS114 + DN(3,1)*crLHS262 + DN(3,2)*crLHS264 - crLHS242*crLHS584 - crLHS266*crLHS462 + crLHS467 + crLHS485);
rLHS(13,7)+=-gauss_weight*(crLHS182*crLHS269 + crLHS269*crLHS579 + crLHS270*crLHS462 + crLHS506);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS271 + DN(3,2)*crLHS272 - crLHS227*crLHS586 - crLHS275*crLHS462 + crLHS533 + crLHS546);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS151 + DN(3,1)*crLHS277 + DN(3,2)*crLHS279 + crLHS182*crLHS280 + crLHS280*crLHS579 - crLHS282*crLHS462 + crLHS548 + crLHS585);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS164 + DN(3,1)*crLHS285 + DN(3,2)*crLHS287 - crLHS242*crLHS586 - crLHS289*crLHS462 + crLHS550 + crLHS557);
rLHS(13,11)+=-gauss_weight*(crLHS182*crLHS292 + crLHS292*crLHS579 + crLHS293*crLHS462 + crLHS574);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS177 + DN(3,1)*crLHS294 + DN(3,2)*crLHS295 + crLHS15*crLHS589 - crLHS227*crLHS594 - crLHS227*crLHS595 - crLHS298*crLHS462 + crLHS593);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS201 + DN(3,1)*crLHS300 + DN(3,2)*crLHS302 + crLHS16*crLHS589 + crLHS182*crLHS303 + crLHS22*crLHS598 + crLHS303*crLHS579 - crLHS305*crLHS462 + crLHS591);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS214 + DN(3,1)*crLHS308 + DN(3,2)*crLHS310 + crLHS17*crLHS589 - crLHS242*crLHS594 - crLHS242*crLHS595 - crLHS312*crLHS462 + crLHS599);
rLHS(13,15)+=-gauss_weight*(crLHS182*crLHS315 + crLHS315*crLHS579 + crLHS316*crLHS462 + crLHS600);
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS226 + DN(3,2)*crLHS317 + crLHS217 - crLHS318*crLHS582 - crLHS319*crLHS462 + crLHS370);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS50 + DN(3,1)*crLHS232 + DN(3,2)*crLHS320 + crLHS311 - crLHS321*crLHS582 - crLHS322*crLHS462 + crLHS374);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS65 + DN(3,1)*crLHS239 + DN(3,2)*crLHS323 + crLHS182*crLHS325 + crLHS325*crLHS579 - crLHS326*crLHS462 + crLHS379 + crLHS580);
rLHS(14,3)+=-gauss_weight*(crLHS182*crLHS328 + crLHS328*crLHS579 + crLHS329*crLHS462 + crLHS406);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS78 + DN(3,1)*crLHS248 + DN(3,2)*crLHS330 - crLHS318*crLHS584 - crLHS334*crLHS462 + crLHS444 + crLHS484);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS104 + DN(3,1)*crLHS256 + DN(3,2)*crLHS337 - crLHS321*crLHS584 - crLHS339*crLHS462 + crLHS466 + crLHS486);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS116 + DN(3,1)*crLHS264 + DN(3,2)*crLHS341 + crLHS182*crLHS342 + crLHS342*crLHS579 - crLHS343*crLHS462 + crLHS488 + crLHS583);
rLHS(14,7)+=-gauss_weight*(crLHS182*crLHS347 + crLHS347*crLHS579 + crLHS348*crLHS462 + crLHS508);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS129 + DN(3,1)*crLHS272 + DN(3,2)*crLHS349 - crLHS318*crLHS586 - crLHS351*crLHS462 + crLHS537 + crLHS556);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS154 + DN(3,1)*crLHS279 + DN(3,2)*crLHS354 - crLHS321*crLHS586 - crLHS356*crLHS462 + crLHS549 + crLHS558);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS166 + DN(3,1)*crLHS287 + DN(3,2)*crLHS358 + crLHS182*crLHS359 + crLHS359*crLHS579 - crLHS360*crLHS462 + crLHS560 + crLHS585);
rLHS(14,11)+=-gauss_weight*(crLHS182*crLHS364 + crLHS364*crLHS579 + crLHS365*crLHS462 + crLHS576);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS179 + DN(3,1)*crLHS295 + DN(3,2)*crLHS366 + crLHS18*crLHS589 - crLHS318*crLHS594 - crLHS318*crLHS595 - crLHS368*crLHS462 + crLHS596);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS204 + DN(3,1)*crLHS302 + DN(3,2)*crLHS371 + crLHS19*crLHS589 - crLHS321*crLHS594 - crLHS321*crLHS595 - crLHS373*crLHS462 + crLHS599);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS216 + DN(3,1)*crLHS310 + DN(3,2)*crLHS375 + crLHS182*crLHS376 + crLHS20*crLHS589 + crLHS22*crLHS601 + crLHS376*crLHS579 - crLHS377*crLHS462 + crLHS591);
rLHS(14,15)+=-gauss_weight*(crLHS182*crLHS381 + crLHS381*crLHS579 + crLHS382*crLHS462 + crLHS602);
rLHS(15,0)+=gauss_weight*(crLHS222 - crLHS223*crLHS32 + crLHS227*crLHS603 + crLHS318*crLHS604);
rLHS(15,1)+=gauss_weight*(crLHS314 - crLHS315*crLHS59 + crLHS321*crLHS604 + crLHS53*crLHS605);
rLHS(15,2)+=gauss_weight*(crLHS242*crLHS603 + crLHS380 - crLHS381*crLHS69 + crLHS605*crLHS67);
rLHS(15,3)+=crLHS407;
rLHS(15,4)+=gauss_weight*(-crLHS223*crLHS86 + crLHS227*crLHS606 + crLHS318*crLHS607 + crLHS446);
rLHS(15,5)+=gauss_weight*(-crLHS107*crLHS315 + crLHS321*crLHS607 + crLHS468 + crLHS53*crLHS608);
rLHS(15,6)+=gauss_weight*(-crLHS119*crLHS381 + crLHS242*crLHS606 + crLHS489 + crLHS608*crLHS67);
rLHS(15,7)+=crLHS509;
rLHS(15,8)+=gauss_weight*(-crLHS137*crLHS223 + crLHS227*crLHS609 + crLHS318*crLHS610 + crLHS539);
rLHS(15,9)+=gauss_weight*(-crLHS157*crLHS315 + crLHS321*crLHS610 + crLHS53*crLHS611 + crLHS551);
rLHS(15,10)+=gauss_weight*(-crLHS169*crLHS381 + crLHS242*crLHS609 + crLHS561 + crLHS611*crLHS67);
rLHS(15,11)+=crLHS577;
rLHS(15,12)+=gauss_weight*(-crLHS187*crLHS223 + crLHS227*crLHS612 + crLHS318*crLHS613 + crLHS597);
rLHS(15,13)+=gauss_weight*(-crLHS207*crLHS315 + crLHS321*crLHS613 + crLHS53*crLHS614 + crLHS600);
rLHS(15,14)+=gauss_weight*(-crLHS219*crLHS381 + crLHS242*crLHS612 + crLHS602 + crLHS614*crLHS67);
rLHS(15,15)+=crLHS386*(crLHS587 + crLHS598 + crLHS601);

}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<2,3,true>>::ComputeGaussPointRHSContributionAdjoint(
    FluidTopologyOptimizationElementData<2,3,true>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,3>& c = rData.SoundVelocity;

    const array_1d<double,3> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,3>& N = rData.N;
    const BoundedMatrix<double,3,2>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const BoundedMatrix<double,2,3>& v_adj = rData.Velocity_adj;
    const BoundedMatrix<double,2,3>& vn_adj = rData.Velocity_adj_OldStep1;
    const BoundedMatrix<double,2,3>& vnn_adj = rData.Velocity_adj_OldStep2;
    const BoundedMatrix<double,2,3>& vmesh_adj = rData.MeshVelocity_adj;
    const BoundedMatrix<double,2,3> v_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
    const BoundedMatrix<double,2,3>& f_adj = rData.BodyForce_adj;
    const array_1d<double,3>& p_adj = rData.Pressure_adj;
    const array_1d<double,3>& t = rData.Temperature;
    const array_1d<double,3>& t_adj = rData.Temperature_adj;
    const array_1d<double,3> t_ConvCoeff = rData.TransportCouplingConvectionCoefficient;
    // const array_1d<double,3>& pn_adj = rData.Pressure_adj_OldStep1;
    // const array_1d<double,3>& pnn_adj = rData.Pressure_adj_OldStep2;
    const array_1d<double,3>& stress_adj = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*p_adj[0] + N[1]*p_adj[1] + N[2]*p_adj[2];
const double crRHS1 = rho*(N[0]*f_adj(0,0) + N[1]*f_adj(1,0) + N[2]*f_adj(2,0));
const double crRHS2 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crRHS3 = N[0]*v_adj(0,0) + N[1]*v_adj(1,0) + N[2]*v_adj(2,0);
const double crRHS4 = crRHS2*crRHS3;
const double crRHS5 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0);
const double crRHS6 = 2.0*crRHS2*functional_weights[0];
const double crRHS7 = crRHS5*crRHS6;
const double crRHS8 = DN(0,1)*v_ns(0,0);
const double crRHS9 = DN(1,1)*v_ns(1,0);
const double crRHS10 = DN(2,1)*v_ns(2,0);
const double crRHS11 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1);
const double crRHS12 = 2.0*functional_weights[2]*mu*(-crRHS10 + crRHS11 - crRHS8 - crRHS9);
const double crRHS13 = (N[0]*t_ConvCoeff[0] + N[1]*t_ConvCoeff[1] + N[2]*t_ConvCoeff[2])*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2]);
const double crRHS14 = crRHS13*(DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2]);
const double crRHS15 = N[0]*crRHS14;
const double crRHS16 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1);
const double crRHS17 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS16);
const double crRHS18 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crRHS19 = 1.0*crRHS18;
const double crRHS20 = crRHS10 + crRHS8 + crRHS9;
const double crRHS21 = 0.5*crRHS11 + 0.5*crRHS20;
const double crRHS22 = 4.0*functional_weights[1]*mu;
const double crRHS23 = N[0]*(bdf0*v_adj(0,0) + bdf1*vn_adj(0,0) + bdf2*vnn_adj(0,0)) + N[1]*(bdf0*v_adj(1,0) + bdf1*vn_adj(1,0) + bdf2*vnn_adj(1,0)) + N[2]*(bdf0*v_adj(2,0) + bdf1*vn_adj(2,0) + bdf2*vnn_adj(2,0));
const double crRHS24 = N[0]*rho;
const double crRHS25 = crRHS18*crRHS3;
const double crRHS26 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1);
const double crRHS27 = crRHS11*crRHS26;
const double crRHS28 = crRHS25 + crRHS27;
const double crRHS29 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0);
const double crRHS30 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1);
const double crRHS31 = crRHS29 + crRHS30;
const double crRHS32 = crRHS2*stab_c3;
const double crRHS33 = rho*stab_c2*sqrt(pow(crRHS16, 2) + pow(crRHS5, 2));
const double crRHS34 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crRHS35 = rho*stab_c3*sqrt(pow(crRHS11, 2) + pow(crRHS18, 2) + pow(crRHS20, 2) + pow(crRHS34, 2));
const double crRHS36 = crRHS31*(h*(crRHS32*h + crRHS33 + crRHS35*h)/stab_c1 + mu);
const double crRHS37 = crRHS5*rho;
const double crRHS38 = crRHS16*rho;
const double crRHS39 = crRHS14*functional_weights[6];
const double crRHS40 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS14 + crRHS25*rho + crRHS27*rho - crRHS29*crRHS37 - crRHS38*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0)) + crRHS39 + crRHS4 + crRHS7;
const double crRHS41 = 1.0/(crRHS32 + crRHS33/h + crRHS35 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS42 = crRHS40*crRHS41;
const double crRHS43 = N[0]*crRHS2;
const double crRHS44 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS45 = crRHS2*crRHS26;
const double crRHS46 = crRHS20*crRHS3;
const double crRHS47 = crRHS26*crRHS34;
const double crRHS48 = crRHS16*crRHS6;
const double crRHS49 = crRHS13*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2]);
const double crRHS50 = crRHS49*functional_weights[6];
const double crRHS51 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] - crRHS30*crRHS38 - crRHS37*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS44 + crRHS45 + crRHS46*rho + crRHS47*rho + crRHS48 + crRHS49 + crRHS50;
const double crRHS52 = crRHS11*crRHS51 + crRHS18*crRHS40;
const double crRHS53 = crRHS24*crRHS41;
const double crRHS54 = 1.0*crRHS34;
const double crRHS55 = N[0]*(bdf0*v_adj(0,1) + bdf1*vn_adj(0,1) + bdf2*vnn_adj(0,1)) + N[1]*(bdf0*v_adj(1,1) + bdf1*vn_adj(1,1) + bdf2*vnn_adj(1,1)) + N[2]*(bdf0*v_adj(2,1) + bdf1*vn_adj(2,1) + bdf2*vnn_adj(2,1));
const double crRHS56 = crRHS46 + crRHS47;
const double crRHS57 = crRHS41*crRHS51;
const double crRHS58 = crRHS20*crRHS40 + crRHS34*crRHS51;
const double crRHS59 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS16);
const double crRHS60 = N[1]*rho;
const double crRHS61 = N[1]*crRHS2;
const double crRHS62 = crRHS41*crRHS60;
const double crRHS63 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS16);
const double crRHS64 = N[2]*rho;
const double crRHS65 = N[2]*crRHS2;
const double crRHS66 = crRHS41*crRHS64;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS36 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS12 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS15*functional_weights[6] + crRHS15 + crRHS17*crRHS3 - crRHS17*crRHS42 + crRHS22*(DN(0,0)*crRHS19 + DN(0,1)*crRHS21) + crRHS23*crRHS24 + crRHS24*crRHS28 - crRHS42*crRHS43 - crRHS52*crRHS53);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS12 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS36 + DN(0,1)*stress_adj[1] - N[0]*crRHS44 + N[0]*crRHS45 + N[0]*crRHS48 + N[0]*crRHS49 + N[0]*crRHS50 + crRHS17*crRHS26 - crRHS17*crRHS57 + crRHS22*(DN(0,0)*crRHS21 + DN(0,1)*crRHS54) + crRHS24*crRHS55 + crRHS24*crRHS56 - crRHS43*crRHS57 - crRHS53*crRHS58);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS42 + DN(0,1)*crRHS57 + N[0]*crRHS31);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS36 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS12 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS14 + N[1]*crRHS39 + N[1]*crRHS4 + N[1]*crRHS7 + crRHS22*(DN(1,0)*crRHS19 + DN(1,1)*crRHS21) + crRHS23*crRHS60 + crRHS28*crRHS60 + crRHS3*crRHS59 - crRHS42*crRHS59 - crRHS42*crRHS61 - crRHS52*crRHS62);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS12 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS36 + DN(1,1)*stress_adj[1] - N[1]*crRHS44 + N[1]*crRHS45 + N[1]*crRHS48 + N[1]*crRHS49 + N[1]*crRHS50 + crRHS22*(DN(1,0)*crRHS21 + DN(1,1)*crRHS54) + crRHS26*crRHS59 + crRHS55*crRHS60 + crRHS56*crRHS60 - crRHS57*crRHS59 - crRHS57*crRHS61 - crRHS58*crRHS62);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS42 + DN(1,1)*crRHS57 + N[1]*crRHS31);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS36 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS12 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS14 + N[2]*crRHS39 + N[2]*crRHS4 + N[2]*crRHS7 + crRHS22*(DN(2,0)*crRHS19 + DN(2,1)*crRHS21) + crRHS23*crRHS64 + crRHS28*crRHS64 + crRHS3*crRHS63 - crRHS42*crRHS63 - crRHS42*crRHS65 - crRHS52*crRHS66);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS12 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS36 + DN(2,1)*stress_adj[1] - N[2]*crRHS44 + N[2]*crRHS45 + N[2]*crRHS48 + N[2]*crRHS49 + N[2]*crRHS50 + crRHS22*(DN(2,0)*crRHS21 + DN(2,1)*crRHS54) + crRHS26*crRHS63 + crRHS55*crRHS64 + crRHS56*crRHS64 - crRHS57*crRHS63 - crRHS57*crRHS65 - crRHS58*crRHS66);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS42 + DN(2,1)*crRHS57 + N[2]*crRHS31);

}

template <>
void FluidTopologyOptimizationElement<FluidTopologyOptimizationElementData<3,4,true>>::ComputeGaussPointRHSContributionAdjoint(
    FluidTopologyOptimizationElementData<3,4,true>& rData,
    VectorType& rRHS)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;
    // const array_1d<double,4>& c = rData.SoundVelocity;

    const array_1d<double,4> alpha = rData.Resistance;

    const double h = rData.ElementSize;
    
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.bdf0;
    const double bdf1 = rData.bdf1;
    const double bdf2 = rData.bdf2;

    const double dyn_tau = rData.DynamicTau;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    Vector functional_weights = rData.Functional_Weights; //  functional terms weights

    const BoundedMatrix<double,3,4>& v_adj = rData.Velocity_adj;
    const BoundedMatrix<double,3,4>& vn_adj = rData.Velocity_adj_OldStep1;
    const BoundedMatrix<double,3,4>& vnn_adj = rData.Velocity_adj_OldStep2;
    const BoundedMatrix<double,3,4>& vmesh_adj = rData.MeshVelocity_adj;
    const BoundedMatrix<double,3,4> v_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
    const BoundedMatrix<double,3,4>& f_adj = rData.BodyForce_adj;
    const array_1d<double,4>& p_adj = rData.Pressure_adj;
    // const array_1d<double,4>& pn_adj = rData.Pressure_adj_OldStep1;
    // const array_1d<double,4>& pnn_adj = rData.Pressure_adj_OldStep2;
    const array_1d<double,4>& t = rData.Temperature;
    const array_1d<double,4>& t_adj = rData.Temperature_adj;
    const array_1d<double,4> t_ConvCoeff = rData.TransportCouplingConvectionCoefficient;
    const array_1d<double,6>& stress_adj = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*p_adj[0] + N[1]*p_adj[1] + N[2]*p_adj[2] + N[3]*p_adj[3];
const double crRHS1 = rho*(N[0]*f_adj(0,0) + N[1]*f_adj(1,0) + N[2]*f_adj(2,0) + N[3]*f_adj(3,0));
const double crRHS2 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crRHS3 = N[0]*v_adj(0,0) + N[1]*v_adj(1,0) + N[2]*v_adj(2,0) + N[3]*v_adj(3,0);
const double crRHS4 = crRHS2*crRHS3;
const double crRHS5 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0) + N[3]*v_ns(3,0);
const double crRHS6 = 2.0*crRHS2*functional_weights[0];
const double crRHS7 = crRHS5*crRHS6;
const double crRHS8 = (N[0]*t_ConvCoeff[0] + N[1]*t_ConvCoeff[1] + N[2]*t_ConvCoeff[2] + N[3]*t_ConvCoeff[3])*(N[0]*t_adj[0] + N[1]*t_adj[1] + N[2]*t_adj[2] + N[3]*t_adj[3]);
const double crRHS9 = crRHS8*(DN(0,0)*t[0] + DN(1,0)*t[1] + DN(2,0)*t[2] + DN(3,0)*t[3]);
const double crRHS10 = N[0]*crRHS9;
const double crRHS11 = N[0]*(bdf0*v_adj(0,0) + bdf1*vn_adj(0,0) + bdf2*vnn_adj(0,0)) + N[1]*(bdf0*v_adj(1,0) + bdf1*vn_adj(1,0) + bdf2*vnn_adj(1,0)) + N[2]*(bdf0*v_adj(2,0) + bdf1*vn_adj(2,0) + bdf2*vnn_adj(2,0)) + N[3]*(bdf0*v_adj(3,0) + bdf1*vn_adj(3,0) + bdf2*vnn_adj(3,0));
const double crRHS12 = N[0]*rho;
const double crRHS13 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crRHS14 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crRHS15 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS13 + DN(0,2)*crRHS14);
const double crRHS16 = DN(0,1)*v_ns(0,0);
const double crRHS17 = DN(1,1)*v_ns(1,0);
const double crRHS18 = DN(2,1)*v_ns(2,0);
const double crRHS19 = DN(3,1)*v_ns(3,0);
const double crRHS20 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crRHS21 = -crRHS16 - crRHS17 - crRHS18 - crRHS19 + crRHS20;
const double crRHS22 = DN(0,2)*v_ns(0,0);
const double crRHS23 = DN(1,2)*v_ns(1,0);
const double crRHS24 = DN(2,2)*v_ns(2,0);
const double crRHS25 = DN(3,2)*v_ns(3,0);
const double crRHS26 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crRHS27 = -crRHS22 - crRHS23 - crRHS24 - crRHS25 + crRHS26;
const double crRHS28 = 2.0*functional_weights[2]*mu;
const double crRHS29 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crRHS30 = 1.0*crRHS29;
const double crRHS31 = crRHS16 + crRHS17 + crRHS18 + crRHS19;
const double crRHS32 = 0.5*crRHS20 + 0.5*crRHS31;
const double crRHS33 = crRHS22 + crRHS23 + crRHS24 + crRHS25;
const double crRHS34 = crRHS26 + crRHS33;
const double crRHS35 = 0.5*DN(0,2);
const double crRHS36 = 4.0*functional_weights[1]*mu;
const double crRHS37 = crRHS29*crRHS3;
const double crRHS38 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1) + N[3]*v_adj(3,1);
const double crRHS39 = crRHS20*crRHS38;
const double crRHS40 = N[0]*v_adj(0,2) + N[1]*v_adj(1,2) + N[2]*v_adj(2,2) + N[3]*v_adj(3,2);
const double crRHS41 = crRHS26*crRHS40;
const double crRHS42 = crRHS37 + crRHS39 + crRHS41;
const double crRHS43 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0) + DN(3,0)*v_adj(3,0);
const double crRHS44 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1) + DN(3,1)*v_adj(3,1);
const double crRHS45 = DN(0,2)*v_adj(0,2) + DN(1,2)*v_adj(1,2) + DN(2,2)*v_adj(2,2) + DN(3,2)*v_adj(3,2);
const double crRHS46 = crRHS43 + crRHS44 + crRHS45;
const double crRHS47 = crRHS2*stab_c3;
const double crRHS48 = rho*stab_c2*sqrt(pow(crRHS13, 2) + pow(crRHS14, 2) + pow(crRHS5, 2));
const double crRHS49 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crRHS50 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crRHS51 = DN(0,2)*v_ns(0,1);
const double crRHS52 = DN(1,2)*v_ns(1,1);
const double crRHS53 = DN(2,2)*v_ns(2,1);
const double crRHS54 = DN(3,2)*v_ns(3,1);
const double crRHS55 = crRHS51 + crRHS52 + crRHS53 + crRHS54;
const double crRHS56 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crRHS57 = rho*stab_c3*sqrt(pow(crRHS20, 2) + pow(crRHS26, 2) + pow(crRHS29, 2) + pow(crRHS31, 2) + pow(crRHS33, 2) + pow(crRHS49, 2) + pow(crRHS50, 2) + pow(crRHS55, 2) + pow(crRHS56, 2));
const double crRHS58 = crRHS46*(h*(crRHS47*h + crRHS48 + crRHS57*h)/stab_c1 + mu);
const double crRHS59 = crRHS5*rho;
const double crRHS60 = crRHS13*rho;
const double crRHS61 = crRHS14*rho;
const double crRHS62 = crRHS9*functional_weights[6];
const double crRHS63 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 + crRHS37*rho + crRHS39*rho + crRHS4 + crRHS41*rho - crRHS43*crRHS59 - crRHS60*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS61*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0)) + crRHS62 + crRHS7 + crRHS9;
const double crRHS64 = 1.0/(crRHS47 + crRHS48/h + crRHS57 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS65 = crRHS63*crRHS64;
const double crRHS66 = N[0]*crRHS2;
const double crRHS67 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS68 = crRHS2*crRHS38;
const double crRHS69 = crRHS3*crRHS31;
const double crRHS70 = crRHS38*crRHS49;
const double crRHS71 = crRHS40*crRHS50;
const double crRHS72 = crRHS13*crRHS6;
const double crRHS73 = crRHS8*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3]);
const double crRHS74 = crRHS73*functional_weights[6];
const double crRHS75 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] - crRHS44*crRHS60 - crRHS59*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS61*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS67 + crRHS68 + crRHS69*rho + crRHS70*rho + crRHS71*rho + crRHS72 + crRHS73 + crRHS74;
const double crRHS76 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS77 = crRHS2*crRHS40;
const double crRHS78 = crRHS3*crRHS33;
const double crRHS79 = crRHS38*crRHS55;
const double crRHS80 = crRHS40*crRHS56;
const double crRHS81 = crRHS14*crRHS6;
const double crRHS82 = crRHS8*(DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3]);
const double crRHS83 = crRHS82*functional_weights[6];
const double crRHS84 = DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] - crRHS45*crRHS61 - crRHS59*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS60*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS76 + crRHS77 + crRHS78*rho + crRHS79*rho + crRHS80*rho + crRHS81 + crRHS82 + crRHS83;
const double crRHS85 = crRHS20*crRHS75 + crRHS26*crRHS84 + crRHS29*crRHS63;
const double crRHS86 = crRHS12*crRHS64;
const double crRHS87 = N[0]*(bdf0*v_adj(0,1) + bdf1*vn_adj(0,1) + bdf2*vnn_adj(0,1)) + N[1]*(bdf0*v_adj(1,1) + bdf1*vn_adj(1,1) + bdf2*vnn_adj(1,1)) + N[2]*(bdf0*v_adj(2,1) + bdf1*vn_adj(2,1) + bdf2*vnn_adj(2,1)) + N[3]*(bdf0*v_adj(3,1) + bdf1*vn_adj(3,1) + bdf2*vnn_adj(3,1));
const double crRHS88 = crRHS50 - crRHS51 - crRHS52 - crRHS53 - crRHS54;
const double crRHS89 = 1.0*crRHS49;
const double crRHS90 = crRHS50 + crRHS55;
const double crRHS91 = crRHS69 + crRHS70 + crRHS71;
const double crRHS92 = crRHS64*crRHS75;
const double crRHS93 = crRHS31*crRHS63 + crRHS49*crRHS75 + crRHS50*crRHS84;
const double crRHS94 = N[0]*(bdf0*v_adj(0,2) + bdf1*vn_adj(0,2) + bdf2*vnn_adj(0,2)) + N[1]*(bdf0*v_adj(1,2) + bdf1*vn_adj(1,2) + bdf2*vnn_adj(1,2)) + N[2]*(bdf0*v_adj(2,2) + bdf1*vn_adj(2,2) + bdf2*vnn_adj(2,2)) + N[3]*(bdf0*v_adj(3,2) + bdf1*vn_adj(3,2) + bdf2*vnn_adj(3,2));
const double crRHS95 = 1.0*crRHS56;
const double crRHS96 = 0.5*crRHS34;
const double crRHS97 = 0.5*crRHS90;
const double crRHS98 = crRHS78 + crRHS79 + crRHS80;
const double crRHS99 = crRHS64*crRHS84;
const double crRHS100 = crRHS33*crRHS63 + crRHS55*crRHS75 + crRHS56*crRHS84;
const double crRHS101 = N[1]*rho;
const double crRHS102 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS13 + DN(1,2)*crRHS14);
const double crRHS103 = N[1]*crRHS2;
const double crRHS104 = crRHS101*crRHS64;
const double crRHS105 = N[2]*rho;
const double crRHS106 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS13 + DN(2,2)*crRHS14);
const double crRHS107 = N[2]*crRHS2;
const double crRHS108 = crRHS105*crRHS64;
const double crRHS109 = N[3]*rho;
const double crRHS110 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS13 + DN(3,2)*crRHS14);
const double crRHS111 = N[3]*crRHS2;
const double crRHS112 = crRHS109*crRHS64;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS58 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS10*functional_weights[6] + crRHS10 + crRHS11*crRHS12 + crRHS12*crRHS42 + crRHS15*crRHS3 - crRHS15*crRHS65 - crRHS28*(DN(0,1)*crRHS21 + DN(0,2)*crRHS27) + crRHS36*(DN(0,0)*crRHS30 + DN(0,1)*crRHS32 + crRHS34*crRHS35) - crRHS65*crRHS66 - crRHS85*crRHS86);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS58 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS67 + N[0]*crRHS68 + N[0]*crRHS72 + N[0]*crRHS73 + N[0]*crRHS74 + crRHS12*crRHS87 + crRHS12*crRHS91 + crRHS15*crRHS38 - crRHS15*crRHS92 + crRHS28*(DN(0,0)*crRHS21 - DN(0,2)*crRHS88) + crRHS36*(DN(0,0)*crRHS32 + DN(0,1)*crRHS89 + crRHS35*crRHS90) - crRHS66*crRHS92 - crRHS86*crRHS93);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS58 + DN(0,2)*stress_adj[2] - N[0]*crRHS76 + N[0]*crRHS77 + N[0]*crRHS81 + N[0]*crRHS82 + N[0]*crRHS83 - crRHS100*crRHS86 + crRHS12*crRHS94 + crRHS12*crRHS98 + crRHS15*crRHS40 - crRHS15*crRHS99 + crRHS28*(DN(0,0)*crRHS27 + DN(0,1)*crRHS88) + crRHS36*(DN(0,0)*crRHS96 + DN(0,1)*crRHS97 + DN(0,2)*crRHS95) - crRHS66*crRHS99);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS65 + DN(0,1)*crRHS92 + DN(0,2)*crRHS99 + N[0]*crRHS46);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS58 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS4 + N[1]*crRHS62 + N[1]*crRHS7 + N[1]*crRHS9 + crRHS101*crRHS11 + crRHS101*crRHS42 + crRHS102*crRHS3 - crRHS102*crRHS65 - crRHS103*crRHS65 - crRHS104*crRHS85 - crRHS28*(DN(1,1)*crRHS21 + DN(1,2)*crRHS27) + crRHS36*(DN(1,0)*crRHS30 + DN(1,1)*crRHS32 + DN(1,2)*crRHS96));
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS58 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS67 + N[1]*crRHS68 + N[1]*crRHS72 + N[1]*crRHS73 + N[1]*crRHS74 + crRHS101*crRHS87 + crRHS101*crRHS91 + crRHS102*crRHS38 - crRHS102*crRHS92 - crRHS103*crRHS92 - crRHS104*crRHS93 + crRHS28*(DN(1,0)*crRHS21 - DN(1,2)*crRHS88) + crRHS36*(DN(1,0)*crRHS32 + DN(1,1)*crRHS89 + DN(1,2)*crRHS97));
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS58 + DN(1,2)*stress_adj[2] - N[1]*crRHS76 + N[1]*crRHS77 + N[1]*crRHS81 + N[1]*crRHS82 + N[1]*crRHS83 - crRHS100*crRHS104 + crRHS101*crRHS94 + crRHS101*crRHS98 + crRHS102*crRHS40 - crRHS102*crRHS99 - crRHS103*crRHS99 + crRHS28*(DN(1,0)*crRHS27 + DN(1,1)*crRHS88) + crRHS36*(DN(1,0)*crRHS96 + DN(1,1)*crRHS97 + DN(1,2)*crRHS95));
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS65 + DN(1,1)*crRHS92 + DN(1,2)*crRHS99 + N[1]*crRHS46);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS58 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS4 + N[2]*crRHS62 + N[2]*crRHS7 + N[2]*crRHS9 + crRHS105*crRHS11 + crRHS105*crRHS42 + crRHS106*crRHS3 - crRHS106*crRHS65 - crRHS107*crRHS65 - crRHS108*crRHS85 - crRHS28*(DN(2,1)*crRHS21 + DN(2,2)*crRHS27) + crRHS36*(DN(2,0)*crRHS30 + DN(2,1)*crRHS32 + DN(2,2)*crRHS96));
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS58 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS67 + N[2]*crRHS68 + N[2]*crRHS72 + N[2]*crRHS73 + N[2]*crRHS74 + crRHS105*crRHS87 + crRHS105*crRHS91 + crRHS106*crRHS38 - crRHS106*crRHS92 - crRHS107*crRHS92 - crRHS108*crRHS93 + crRHS28*(DN(2,0)*crRHS21 - DN(2,2)*crRHS88) + crRHS36*(DN(2,0)*crRHS32 + DN(2,1)*crRHS89 + DN(2,2)*crRHS97));
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS58 + DN(2,2)*stress_adj[2] - N[2]*crRHS76 + N[2]*crRHS77 + N[2]*crRHS81 + N[2]*crRHS82 + N[2]*crRHS83 - crRHS100*crRHS108 + crRHS105*crRHS94 + crRHS105*crRHS98 + crRHS106*crRHS40 - crRHS106*crRHS99 - crRHS107*crRHS99 + crRHS28*(DN(2,0)*crRHS27 + DN(2,1)*crRHS88) + crRHS36*(DN(2,0)*crRHS96 + DN(2,1)*crRHS97 + DN(2,2)*crRHS95));
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS65 + DN(2,1)*crRHS92 + DN(2,2)*crRHS99 + N[2]*crRHS46);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS58 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS4 + N[3]*crRHS62 + N[3]*crRHS7 + N[3]*crRHS9 + crRHS109*crRHS11 + crRHS109*crRHS42 + crRHS110*crRHS3 - crRHS110*crRHS65 - crRHS111*crRHS65 - crRHS112*crRHS85 - crRHS28*(DN(3,1)*crRHS21 + DN(3,2)*crRHS27) + crRHS36*(DN(3,0)*crRHS30 + DN(3,1)*crRHS32 + DN(3,2)*crRHS96));
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS58 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS67 + N[3]*crRHS68 + N[3]*crRHS72 + N[3]*crRHS73 + N[3]*crRHS74 + crRHS109*crRHS87 + crRHS109*crRHS91 + crRHS110*crRHS38 - crRHS110*crRHS92 - crRHS111*crRHS92 - crRHS112*crRHS93 + crRHS28*(DN(3,0)*crRHS21 - DN(3,2)*crRHS88) + crRHS36*(DN(3,0)*crRHS32 + DN(3,1)*crRHS89 + DN(3,2)*crRHS97));
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS58 + DN(3,2)*stress_adj[2] - N[3]*crRHS76 + N[3]*crRHS77 + N[3]*crRHS81 + N[3]*crRHS82 + N[3]*crRHS83 - crRHS100*crRHS112 + crRHS109*crRHS94 + crRHS109*crRHS98 + crRHS110*crRHS40 - crRHS110*crRHS99 - crRHS111*crRHS99 + crRHS28*(DN(3,0)*crRHS27 + DN(3,1)*crRHS88) + crRHS36*(DN(3,0)*crRHS96 + DN(3,1)*crRHS97 + DN(3,2)*crRHS95));
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS65 + DN(3,1)*crRHS92 + DN(3,2)*crRHS99 + N[3]*crRHS46);

}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddVelocitySystem(
    TElementData& rData,
    MatrixType& rLocalLHS,
    VectorType& rLocalRHS) 
    {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidTopologyOptimizationElement::AddVelocitySystem "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::AddMassLHS(
    TElementData& rData, MatrixType& rMassMatrix) {
    KRATOS_TRY;

    KRATOS_ERROR << "Calling base FluidTopologyOptimizationElement::AddMassLHS "
                    "implementation. This method is not supported by your "
                    "element."
                 << std::endl;

    KRATOS_CATCH("");
}

template <class TElementData>
void FluidTopologyOptimizationElement<TElementData>::GetCurrentValuesVector(
    const TElementData& rData,
    array_1d<double,LocalSize>& rValues) const {

    int local_index = 0;

    const auto& r_velocities = rData.Velocity;
    const auto& r_pressures = rData.Pressure;

    for (unsigned int i = 0; i < NumNodes; ++i) {
        for (unsigned int d = 0; d < Dim; ++d)  // Velocity Dofs
            rValues[local_index++] = r_velocities(i, d);
        rValues[local_index++] = r_pressures[i];  // Pressure Dof
    }
}

template< class TElementData >
ConstitutiveLaw::Pointer FluidTopologyOptimizationElement<TElementData>::GetConstitutiveLaw() {
    return this->mpConstitutiveLaw;
}

template< class TElementData >
const ConstitutiveLaw::Pointer FluidTopologyOptimizationElement<TElementData>::GetConstitutiveLaw() const {
    return this->mpConstitutiveLaw;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
    rSerializer.save("mpConstitutiveLaw",this->mpConstitutiveLaw);
}


template< class TElementData >
void FluidTopologyOptimizationElement<TElementData>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    rSerializer.load("mpConstitutiveLaw",this->mpConstitutiveLaw);
}

///////////////////////////////////////////////////////////////////////////////////////////////////

namespace Internals {

template< class TElementData >
void StrainRateSpecialization<TElementData,2>::Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX) {

    noalias(rStrainRate) = ZeroVector(3);
    for (unsigned int i = 0; i < TElementData::NumNodes; i++) {
        rStrainRate[0] += rDNDX(i,0)*rVelocities(i,0);
        rStrainRate[1] += rDNDX(i,1)*rVelocities(i,1);
        rStrainRate[2] += rDNDX(i,0)*rVelocities(i,1) + rDNDX(i,1)*rVelocities(i,0);
    }
}

template< class TElementData >
void StrainRateSpecialization<TElementData,3>::Calculate(
    Vector& rStrainRate,
    const typename TElementData::NodalVectorData& rVelocities,
    const typename TElementData::ShapeDerivativesType& rDNDX) {

    noalias(rStrainRate) = ZeroVector(6);
    for (unsigned int i = 0; i < TElementData::NumNodes; i++) {
        rStrainRate[0] += rDNDX(i,0)*rVelocities(i,0);
        rStrainRate[1] += rDNDX(i,1)*rVelocities(i,1);
        rStrainRate[2] += rDNDX(i,2)*rVelocities(i,2);
        rStrainRate[3] += rDNDX(i,0)*rVelocities(i,1) + rDNDX(i,1)*rVelocities(i,0);
        rStrainRate[4] += rDNDX(i,1)*rVelocities(i,2) + rDNDX(i,2)*rVelocities(i,1);
        rStrainRate[5] += rDNDX(i,0)*rVelocities(i,2) + rDNDX(i,2)*rVelocities(i,0);
    }
}

} // namespace Internals

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class FluidTopologyOptimizationElement< FluidTopologyOptimizationElementData<2,3,true> >;
template class FluidTopologyOptimizationElement< FluidTopologyOptimizationElementData<3,4,true> >;

///////////////////////////////////////////////////////////////////////////////////////////////////
}