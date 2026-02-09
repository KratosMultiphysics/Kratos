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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
const double crLHS15 = rho*time_coeff;
const double crLHS16 = bdf0*crLHS15;
const double crLHS17 = bdf0*time_coeff;
const double crLHS18 = crLHS14*crLHS17;
const double crLHS19 = N[0]*crLHS4;
const double crLHS20 = crLHS11*rho;
const double crLHS21 = crLHS12*rho;
const double crLHS22 = crLHS18 + crLHS19 + crLHS20 + crLHS21;
const double crLHS23 = 1.0/(crLHS15*dyn_tau/dt + crLHS5 + crLHS8/h + mu*stab_c1/pow(h, 2));
const double crLHS24 = 1.0*crLHS23;
const double crLHS25 = crLHS24*rho;
const double crLHS26 = crLHS13*crLHS25;
const double crLHS27 = 1.0*crLHS19;
const double crLHS28 = crLHS23*crLHS27;
const double crLHS29 = crLHS22*crLHS24;
const double crLHS30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS31 = crLHS14*crLHS30;
const double crLHS32 = crLHS10*crLHS16 + crLHS10*crLHS4 + crLHS13*crLHS14 + crLHS22*crLHS26 - crLHS22*crLHS28 + crLHS29*crLHS31;
const double crLHS33 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS34 = C(1,2)*DN(0,1);
const double crLHS35 = C(2,2)*DN(0,0) + crLHS34;
const double crLHS36 = DN(0,0)*crLHS9;
const double crLHS37 = DN(0,1)*crLHS36;
const double crLHS38 = crLHS24*crLHS30;
const double crLHS39 = gauss_weight*(-N[0] + crLHS14*crLHS38 + crLHS26 - crLHS28);
const double crLHS40 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS41 = C(0,2)*DN(1,0);
const double crLHS42 = C(2,2)*DN(1,1) + crLHS41;
const double crLHS43 = DN(0,0)*DN(1,0);
const double crLHS44 = N[1]*crLHS18 + N[1]*crLHS19;
const double crLHS45 = crLHS43*crLHS9 + crLHS44;
const double crLHS46 = DN(1,0)*crLHS6;
const double crLHS47 = DN(1,1)*crLHS7;
const double crLHS48 = crLHS46 + crLHS47;
const double crLHS49 = N[1]*rho;
const double crLHS50 = crLHS17*crLHS49;
const double crLHS51 = N[1]*crLHS4;
const double crLHS52 = crLHS46*rho;
const double crLHS53 = crLHS47*rho;
const double crLHS54 = crLHS50 + crLHS51 + crLHS52 + crLHS53;
const double crLHS55 = crLHS24*crLHS54;
const double crLHS56 = crLHS14*crLHS48 + crLHS26*crLHS54 - crLHS28*crLHS54 + crLHS31*crLHS55;
const double crLHS57 = C(0,1)*DN(1,1) + crLHS41;
const double crLHS58 = C(1,2)*DN(1,1);
const double crLHS59 = C(2,2)*DN(1,0) + crLHS58;
const double crLHS60 = DN(1,1)*crLHS36;
const double crLHS61 = DN(0,0)*N[1];
const double crLHS62 = DN(1,0)*N[0];
const double crLHS63 = crLHS25*crLHS30;
const double crLHS64 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS65 = C(0,2)*DN(2,0);
const double crLHS66 = C(2,2)*DN(2,1) + crLHS65;
const double crLHS67 = DN(0,0)*DN(2,0);
const double crLHS68 = N[2]*crLHS18 + N[2]*crLHS19;
const double crLHS69 = crLHS67*crLHS9 + crLHS68;
const double crLHS70 = DN(2,0)*crLHS6;
const double crLHS71 = DN(2,1)*crLHS7;
const double crLHS72 = crLHS70 + crLHS71;
const double crLHS73 = N[2]*rho;
const double crLHS74 = crLHS17*crLHS73;
const double crLHS75 = N[2]*crLHS4;
const double crLHS76 = crLHS70*rho;
const double crLHS77 = crLHS71*rho;
const double crLHS78 = crLHS74 + crLHS75 + crLHS76 + crLHS77;
const double crLHS79 = crLHS24*crLHS78;
const double crLHS80 = crLHS14*crLHS72 + crLHS26*crLHS78 - crLHS28*crLHS78 + crLHS31*crLHS79;
const double crLHS81 = C(0,1)*DN(2,1) + crLHS65;
const double crLHS82 = C(1,2)*DN(2,1);
const double crLHS83 = C(2,2)*DN(2,0) + crLHS82;
const double crLHS84 = DN(2,1)*crLHS36;
const double crLHS85 = DN(0,0)*N[2];
const double crLHS86 = DN(2,0)*N[0];
const double crLHS87 = C(0,1)*DN(0,0) + crLHS34;
const double crLHS88 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS89 = pow(DN(0,1), 2);
const double crLHS90 = C(0,1)*DN(1,0) + crLHS58;
const double crLHS91 = DN(0,1)*crLHS9;
const double crLHS92 = DN(1,0)*crLHS91;
const double crLHS93 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS94 = DN(0,1)*DN(1,1);
const double crLHS95 = crLHS44 + crLHS9*crLHS94;
const double crLHS96 = DN(0,1)*N[1];
const double crLHS97 = DN(1,1)*N[0];
const double crLHS98 = C(0,1)*DN(2,0) + crLHS82;
const double crLHS99 = DN(2,0)*crLHS91;
const double crLHS100 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS101 = DN(0,1)*DN(2,1);
const double crLHS102 = crLHS101*crLHS9 + crLHS68;
const double crLHS103 = DN(0,1)*N[2];
const double crLHS104 = DN(2,1)*N[0];
const double crLHS105 = gauss_weight*(N[0] + crLHS23*(1.0*crLHS18 + 1.0*crLHS20 + 1.0*crLHS21 + crLHS27));
const double crLHS106 = crLHS24*gauss_weight;
const double crLHS107 = crLHS106*(crLHS43 + crLHS94);
const double crLHS108 = crLHS106*(crLHS101 + crLHS67);
const double crLHS109 = crLHS25*crLHS48;
const double crLHS110 = 1.0*crLHS51;
const double crLHS111 = crLHS110*crLHS23;
const double crLHS112 = crLHS30*crLHS49;
const double crLHS113 = crLHS109*crLHS22 - crLHS111*crLHS22 + crLHS112*crLHS29 + crLHS13*crLHS49;
const double crLHS114 = pow(DN(1,0), 2);
const double crLHS115 = pow(N[1], 2);
const double crLHS116 = crLHS109*crLHS54 - crLHS111*crLHS54 + crLHS112*crLHS55 + crLHS115*crLHS16 + crLHS115*crLHS4 + crLHS48*crLHS49;
const double crLHS117 = DN(1,0)*crLHS9;
const double crLHS118 = DN(1,1)*crLHS117;
const double crLHS119 = gauss_weight*(-N[1] + crLHS109 - crLHS111 + crLHS38*crLHS49);
const double crLHS120 = DN(1,0)*DN(2,0);
const double crLHS121 = N[2]*crLHS50 + N[2]*crLHS51;
const double crLHS122 = crLHS120*crLHS9 + crLHS121;
const double crLHS123 = crLHS109*crLHS78 - crLHS111*crLHS78 + crLHS112*crLHS79 + crLHS49*crLHS72;
const double crLHS124 = DN(2,1)*crLHS117;
const double crLHS125 = DN(1,0)*N[2];
const double crLHS126 = DN(2,0)*N[1];
const double crLHS127 = pow(DN(1,1), 2);
const double crLHS128 = DN(2,0)*crLHS9;
const double crLHS129 = DN(1,1)*crLHS128;
const double crLHS130 = DN(1,1)*DN(2,1);
const double crLHS131 = crLHS121 + crLHS130*crLHS9;
const double crLHS132 = DN(1,1)*N[2];
const double crLHS133 = DN(2,1)*N[1];
const double crLHS134 = gauss_weight*(N[1] + crLHS23*(crLHS110 + 1.0*crLHS50 + 1.0*crLHS52 + 1.0*crLHS53));
const double crLHS135 = crLHS106*(crLHS120 + crLHS130);
const double crLHS136 = crLHS25*crLHS72;
const double crLHS137 = 1.0*crLHS75;
const double crLHS138 = crLHS137*crLHS23;
const double crLHS139 = crLHS30*crLHS73;
const double crLHS140 = crLHS13*crLHS73 + crLHS136*crLHS22 - crLHS138*crLHS22 + crLHS139*crLHS29;
const double crLHS141 = crLHS136*crLHS54 - crLHS138*crLHS54 + crLHS139*crLHS55 + crLHS48*crLHS73;
const double crLHS142 = pow(DN(2,0), 2);
const double crLHS143 = pow(N[2], 2);
const double crLHS144 = crLHS136*crLHS78 - crLHS138*crLHS78 + crLHS139*crLHS79 + crLHS143*crLHS16 + crLHS143*crLHS4 + crLHS72*crLHS73;
const double crLHS145 = DN(2,1)*crLHS128;
const double crLHS146 = gauss_weight*(-N[2] + crLHS136 - crLHS138 + crLHS38*crLHS73);
const double crLHS147 = pow(DN(2,1), 2);
const double crLHS148 = gauss_weight*(N[2] + crLHS23*(crLHS137 + 1.0*crLHS74 + 1.0*crLHS76 + 1.0*crLHS77));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS3*crLHS9 + crLHS32);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS33 + DN(0,1)*crLHS35 + crLHS37);
rLHS(0,2)+=DN(0,0)*crLHS39;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS42 + crLHS45 + crLHS56);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS57 + DN(0,1)*crLHS59 + crLHS60);
rLHS(0,5)+=-gauss_weight*(-DN(1,0)*crLHS26 + DN(1,0)*crLHS28 + crLHS61 - crLHS62*crLHS63);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS64 + DN(0,1)*crLHS66 + crLHS69 + crLHS80);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS81 + DN(0,1)*crLHS83 + crLHS84);
rLHS(0,8)+=-gauss_weight*(-DN(2,0)*crLHS26 + DN(2,0)*crLHS28 - crLHS63*crLHS86 + crLHS85);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS87 + crLHS37);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS35 + DN(0,1)*crLHS88 + crLHS32 + crLHS89*crLHS9);
rLHS(1,2)+=DN(0,1)*crLHS39;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS42 + DN(0,1)*crLHS90 + crLHS92);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS59 + DN(0,1)*crLHS93 + crLHS56 + crLHS95);
rLHS(1,5)+=-gauss_weight*(-DN(1,1)*crLHS26 + DN(1,1)*crLHS28 - crLHS63*crLHS97 + crLHS96);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS66 + DN(0,1)*crLHS98 + crLHS99);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS83 + DN(0,1)*crLHS100 + crLHS102 + crLHS80);
rLHS(1,8)+=-gauss_weight*(-DN(2,1)*crLHS26 + DN(2,1)*crLHS28 + crLHS103 - crLHS104*crLHS63);
rLHS(2,0)+=DN(0,0)*crLHS105;
rLHS(2,1)+=DN(0,1)*crLHS105;
rLHS(2,2)+=crLHS106*(crLHS3 + crLHS89);
rLHS(2,3)+=gauss_weight*(DN(0,0)*crLHS55 + crLHS62);
rLHS(2,4)+=gauss_weight*(DN(0,1)*crLHS55 + crLHS97);
rLHS(2,5)+=crLHS107;
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS79 + crLHS86);
rLHS(2,7)+=gauss_weight*(DN(0,1)*crLHS79 + crLHS104);
rLHS(2,8)+=crLHS108;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS113 + crLHS45);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS33 + DN(1,1)*crLHS35 + crLHS92);
rLHS(3,2)+=gauss_weight*(DN(0,0)*crLHS109 - DN(0,0)*crLHS111 + crLHS61*crLHS63 - crLHS62);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS42 + crLHS114*crLHS9 + crLHS116);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS57 + DN(1,1)*crLHS59 + crLHS118);
rLHS(3,5)+=DN(1,0)*crLHS119;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS64 + DN(1,1)*crLHS66 + crLHS122 + crLHS123);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS81 + DN(1,1)*crLHS83 + crLHS124);
rLHS(3,8)+=-gauss_weight*(-DN(2,0)*crLHS109 + DN(2,0)*crLHS111 + crLHS125 - crLHS126*crLHS63);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS87 + crLHS60);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS35 + DN(1,1)*crLHS88 + crLHS113 + crLHS95);
rLHS(4,2)+=gauss_weight*(DN(0,1)*crLHS109 - DN(0,1)*crLHS111 + crLHS63*crLHS96 - crLHS97);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS42 + DN(1,1)*crLHS90 + crLHS118);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS59 + DN(1,1)*crLHS93 + crLHS116 + crLHS127*crLHS9);
rLHS(4,5)+=DN(1,1)*crLHS119;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS66 + DN(1,1)*crLHS98 + crLHS129);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS83 + DN(1,1)*crLHS100 + crLHS123 + crLHS131);
rLHS(4,8)+=-gauss_weight*(-DN(2,1)*crLHS109 + DN(2,1)*crLHS111 + crLHS132 - crLHS133*crLHS63);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS29 + crLHS61);
rLHS(5,1)+=gauss_weight*(DN(1,1)*crLHS29 + crLHS96);
rLHS(5,2)+=crLHS107;
rLHS(5,3)+=DN(1,0)*crLHS134;
rLHS(5,4)+=DN(1,1)*crLHS134;
rLHS(5,5)+=crLHS106*(crLHS114 + crLHS127);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS79 + crLHS126);
rLHS(5,7)+=gauss_weight*(DN(1,1)*crLHS79 + crLHS133);
rLHS(5,8)+=crLHS135;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS140 + crLHS69);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS33 + DN(2,1)*crLHS35 + crLHS99);
rLHS(6,2)+=gauss_weight*(DN(0,0)*crLHS136 - DN(0,0)*crLHS138 + crLHS63*crLHS85 - crLHS86);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS42 + crLHS122 + crLHS141);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS57 + DN(2,1)*crLHS59 + crLHS129);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS136 - DN(1,0)*crLHS138 + crLHS125*crLHS63 - crLHS126);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS64 + DN(2,1)*crLHS66 + crLHS142*crLHS9 + crLHS144);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS81 + DN(2,1)*crLHS83 + crLHS145);
rLHS(6,8)+=DN(2,0)*crLHS146;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS87 + crLHS84);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS35 + DN(2,1)*crLHS88 + crLHS102 + crLHS140);
rLHS(7,2)+=gauss_weight*(DN(0,1)*crLHS136 - DN(0,1)*crLHS138 + crLHS103*crLHS63 - crLHS104);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS42 + DN(2,1)*crLHS90 + crLHS124);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS59 + DN(2,1)*crLHS93 + crLHS131 + crLHS141);
rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS136 - DN(1,1)*crLHS138 + crLHS132*crLHS63 - crLHS133);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS66 + DN(2,1)*crLHS98 + crLHS145);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS83 + DN(2,1)*crLHS100 + crLHS144 + crLHS147*crLHS9);
rLHS(7,8)+=DN(2,1)*crLHS146;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS29 + crLHS85);
rLHS(8,1)+=gauss_weight*(DN(2,1)*crLHS29 + crLHS103);
rLHS(8,2)+=crLHS108;
rLHS(8,3)+=gauss_weight*(DN(2,0)*crLHS55 + crLHS125);
rLHS(8,4)+=gauss_weight*(DN(2,1)*crLHS55 + crLHS132);
rLHS(8,5)+=crLHS135;
rLHS(8,6)+=DN(2,0)*crLHS148;
rLHS(8,7)+=DN(2,1)*crLHS148;
rLHS(8,8)+=crLHS106*(crLHS142 + crLHS147);
    
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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
const double crLHS19 = rho*time_coeff;
const double crLHS20 = bdf0*crLHS19;
const double crLHS21 = bdf0*time_coeff;
const double crLHS22 = crLHS18*crLHS21;
const double crLHS23 = N[0]*crLHS6;
const double crLHS24 = crLHS14*rho;
const double crLHS25 = crLHS15*rho;
const double crLHS26 = crLHS16*rho;
const double crLHS27 = crLHS22 + crLHS23 + crLHS24 + crLHS25 + crLHS26;
const double crLHS28 = 1.0/(crLHS11/h + crLHS19*dyn_tau/dt + crLHS7 + mu*stab_c1/pow(h, 2));
const double crLHS29 = 1.0*crLHS28;
const double crLHS30 = crLHS29*rho;
const double crLHS31 = crLHS17*crLHS30;
const double crLHS32 = 1.0*crLHS23;
const double crLHS33 = crLHS28*crLHS32;
const double crLHS34 = crLHS27*crLHS29;
const double crLHS35 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS36 = crLHS18*crLHS35;
const double crLHS37 = crLHS13*crLHS20 + crLHS13*crLHS6 + crLHS17*crLHS18 + crLHS27*crLHS31 - crLHS27*crLHS33 + crLHS34*crLHS36;
const double crLHS38 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS39 = C(1,3)*DN(0,1);
const double crLHS40 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS39;
const double crLHS41 = C(3,5)*DN(0,0);
const double crLHS42 = C(4,5)*DN(0,2);
const double crLHS43 = C(1,5)*DN(0,1) + crLHS41 + crLHS42;
const double crLHS44 = DN(0,0)*crLHS12;
const double crLHS45 = DN(0,1)*crLHS44;
const double crLHS46 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS47 = C(3,4)*DN(0,1);
const double crLHS48 = C(2,3)*DN(0,2) + crLHS41 + crLHS47;
const double crLHS49 = C(2,5)*DN(0,2);
const double crLHS50 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS49;
const double crLHS51 = DN(0,2)*crLHS44;
const double crLHS52 = crLHS29*crLHS35;
const double crLHS53 = gauss_weight*(-N[0] + crLHS18*crLHS52 + crLHS31 - crLHS33);
const double crLHS54 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS55 = C(0,3)*DN(1,0);
const double crLHS56 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS55;
const double crLHS57 = C(0,5)*DN(1,0);
const double crLHS58 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS57;
const double crLHS59 = DN(0,0)*DN(1,0);
const double crLHS60 = N[1]*crLHS22 + N[1]*crLHS23;
const double crLHS61 = crLHS12*crLHS59 + crLHS60;
const double crLHS62 = DN(1,0)*crLHS8;
const double crLHS63 = DN(1,1)*crLHS9;
const double crLHS64 = DN(1,2)*crLHS10;
const double crLHS65 = crLHS62 + crLHS63 + crLHS64;
const double crLHS66 = N[1]*rho;
const double crLHS67 = crLHS21*crLHS66;
const double crLHS68 = N[1]*crLHS6;
const double crLHS69 = crLHS62*rho;
const double crLHS70 = crLHS63*rho;
const double crLHS71 = crLHS64*rho;
const double crLHS72 = crLHS67 + crLHS68 + crLHS69 + crLHS70 + crLHS71;
const double crLHS73 = crLHS29*crLHS72;
const double crLHS74 = crLHS18*crLHS65 + crLHS31*crLHS72 - crLHS33*crLHS72 + crLHS36*crLHS73;
const double crLHS75 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS55;
const double crLHS76 = C(1,3)*DN(1,1);
const double crLHS77 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS76;
const double crLHS78 = C(3,5)*DN(1,0);
const double crLHS79 = C(4,5)*DN(1,2);
const double crLHS80 = C(1,5)*DN(1,1) + crLHS78 + crLHS79;
const double crLHS81 = DN(1,1)*crLHS44;
const double crLHS82 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS57;
const double crLHS83 = C(3,4)*DN(1,1);
const double crLHS84 = C(2,3)*DN(1,2) + crLHS78 + crLHS83;
const double crLHS85 = C(2,5)*DN(1,2);
const double crLHS86 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS85;
const double crLHS87 = DN(1,2)*crLHS44;
const double crLHS88 = DN(0,0)*N[1];
const double crLHS89 = DN(1,0)*N[0];
const double crLHS90 = crLHS30*crLHS35;
const double crLHS91 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS92 = C(0,3)*DN(2,0);
const double crLHS93 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS92;
const double crLHS94 = C(0,5)*DN(2,0);
const double crLHS95 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS94;
const double crLHS96 = DN(0,0)*DN(2,0);
const double crLHS97 = N[2]*crLHS22 + N[2]*crLHS23;
const double crLHS98 = crLHS12*crLHS96 + crLHS97;
const double crLHS99 = DN(2,0)*crLHS8;
const double crLHS100 = DN(2,1)*crLHS9;
const double crLHS101 = DN(2,2)*crLHS10;
const double crLHS102 = crLHS100 + crLHS101 + crLHS99;
const double crLHS103 = N[2]*rho;
const double crLHS104 = crLHS103*crLHS21;
const double crLHS105 = N[2]*crLHS6;
const double crLHS106 = crLHS99*rho;
const double crLHS107 = crLHS100*rho;
const double crLHS108 = crLHS101*rho;
const double crLHS109 = crLHS104 + crLHS105 + crLHS106 + crLHS107 + crLHS108;
const double crLHS110 = crLHS109*crLHS29;
const double crLHS111 = crLHS102*crLHS18 + crLHS109*crLHS31 - crLHS109*crLHS33 + crLHS110*crLHS36;
const double crLHS112 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS92;
const double crLHS113 = C(1,3)*DN(2,1);
const double crLHS114 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS113;
const double crLHS115 = C(3,5)*DN(2,0);
const double crLHS116 = C(4,5)*DN(2,2);
const double crLHS117 = C(1,5)*DN(2,1) + crLHS115 + crLHS116;
const double crLHS118 = DN(2,1)*crLHS44;
const double crLHS119 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS94;
const double crLHS120 = C(3,4)*DN(2,1);
const double crLHS121 = C(2,3)*DN(2,2) + crLHS115 + crLHS120;
const double crLHS122 = C(2,5)*DN(2,2);
const double crLHS123 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS122;
const double crLHS124 = DN(2,2)*crLHS44;
const double crLHS125 = DN(0,0)*N[2];
const double crLHS126 = DN(2,0)*N[0];
const double crLHS127 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS128 = C(0,3)*DN(3,0);
const double crLHS129 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS128;
const double crLHS130 = C(0,5)*DN(3,0);
const double crLHS131 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS130;
const double crLHS132 = DN(0,0)*DN(3,0);
const double crLHS133 = N[3]*crLHS22 + N[3]*crLHS23;
const double crLHS134 = crLHS12*crLHS132 + crLHS133;
const double crLHS135 = DN(3,0)*crLHS8;
const double crLHS136 = DN(3,1)*crLHS9;
const double crLHS137 = DN(3,2)*crLHS10;
const double crLHS138 = crLHS135 + crLHS136 + crLHS137;
const double crLHS139 = N[3]*rho;
const double crLHS140 = crLHS139*crLHS21;
const double crLHS141 = N[3]*crLHS6;
const double crLHS142 = crLHS135*rho;
const double crLHS143 = crLHS136*rho;
const double crLHS144 = crLHS137*rho;
const double crLHS145 = crLHS140 + crLHS141 + crLHS142 + crLHS143 + crLHS144;
const double crLHS146 = crLHS145*crLHS29;
const double crLHS147 = crLHS138*crLHS18 + crLHS145*crLHS31 - crLHS145*crLHS33 + crLHS146*crLHS36;
const double crLHS148 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS128;
const double crLHS149 = C(1,3)*DN(3,1);
const double crLHS150 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS149;
const double crLHS151 = C(3,5)*DN(3,0);
const double crLHS152 = C(4,5)*DN(3,2);
const double crLHS153 = C(1,5)*DN(3,1) + crLHS151 + crLHS152;
const double crLHS154 = DN(3,1)*crLHS44;
const double crLHS155 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS130;
const double crLHS156 = C(3,4)*DN(3,1);
const double crLHS157 = C(2,3)*DN(3,2) + crLHS151 + crLHS156;
const double crLHS158 = C(2,5)*DN(3,2);
const double crLHS159 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS158;
const double crLHS160 = DN(3,2)*crLHS44;
const double crLHS161 = DN(0,0)*N[3];
const double crLHS162 = DN(3,0)*N[0];
const double crLHS163 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS39;
const double crLHS164 = C(0,4)*DN(0,0) + crLHS42 + crLHS47;
const double crLHS165 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS166 = C(1,4)*DN(0,1);
const double crLHS167 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS166;
const double crLHS168 = pow(DN(0,1), 2);
const double crLHS169 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS166;
const double crLHS170 = C(2,4)*DN(0,2);
const double crLHS171 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS170;
const double crLHS172 = DN(0,1)*crLHS12;
const double crLHS173 = DN(0,2)*crLHS172;
const double crLHS174 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS76;
const double crLHS175 = C(0,4)*DN(1,0) + crLHS79 + crLHS83;
const double crLHS176 = DN(1,0)*crLHS172;
const double crLHS177 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS178 = C(1,4)*DN(1,1);
const double crLHS179 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS178;
const double crLHS180 = DN(0,1)*DN(1,1);
const double crLHS181 = crLHS12*crLHS180;
const double crLHS182 = crLHS60 + crLHS74;
const double crLHS183 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS178;
const double crLHS184 = C(2,4)*DN(1,2);
const double crLHS185 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS184;
const double crLHS186 = DN(1,2)*crLHS172;
const double crLHS187 = DN(0,1)*N[1];
const double crLHS188 = DN(1,1)*N[0];
const double crLHS189 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS113;
const double crLHS190 = C(0,4)*DN(2,0) + crLHS116 + crLHS120;
const double crLHS191 = DN(2,0)*crLHS172;
const double crLHS192 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS193 = C(1,4)*DN(2,1);
const double crLHS194 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS193;
const double crLHS195 = DN(0,1)*DN(2,1);
const double crLHS196 = crLHS12*crLHS195;
const double crLHS197 = crLHS111 + crLHS97;
const double crLHS198 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS193;
const double crLHS199 = C(2,4)*DN(2,2);
const double crLHS200 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS199;
const double crLHS201 = DN(2,2)*crLHS172;
const double crLHS202 = DN(0,1)*N[2];
const double crLHS203 = DN(2,1)*N[0];
const double crLHS204 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS149;
const double crLHS205 = C(0,4)*DN(3,0) + crLHS152 + crLHS156;
const double crLHS206 = DN(3,0)*crLHS172;
const double crLHS207 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS208 = C(1,4)*DN(3,1);
const double crLHS209 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS208;
const double crLHS210 = DN(0,1)*DN(3,1);
const double crLHS211 = crLHS12*crLHS210;
const double crLHS212 = crLHS133 + crLHS147;
const double crLHS213 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS208;
const double crLHS214 = C(2,4)*DN(3,2);
const double crLHS215 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS214;
const double crLHS216 = DN(3,2)*crLHS172;
const double crLHS217 = DN(0,1)*N[3];
const double crLHS218 = DN(3,1)*N[0];
const double crLHS219 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS49;
const double crLHS220 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS170;
const double crLHS221 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS222 = pow(DN(0,2), 2);
const double crLHS223 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS85;
const double crLHS224 = DN(0,2)*crLHS12;
const double crLHS225 = DN(1,0)*crLHS224;
const double crLHS226 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS184;
const double crLHS227 = DN(1,1)*crLHS224;
const double crLHS228 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS229 = DN(0,2)*DN(1,2);
const double crLHS230 = crLHS12*crLHS229;
const double crLHS231 = DN(0,2)*N[1];
const double crLHS232 = DN(1,2)*N[0];
const double crLHS233 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS122;
const double crLHS234 = DN(2,0)*crLHS224;
const double crLHS235 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS199;
const double crLHS236 = DN(2,1)*crLHS224;
const double crLHS237 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS238 = DN(0,2)*DN(2,2);
const double crLHS239 = crLHS12*crLHS238;
const double crLHS240 = DN(0,2)*N[2];
const double crLHS241 = DN(2,2)*N[0];
const double crLHS242 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS158;
const double crLHS243 = DN(3,0)*crLHS224;
const double crLHS244 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS214;
const double crLHS245 = DN(3,1)*crLHS224;
const double crLHS246 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS247 = DN(0,2)*DN(3,2);
const double crLHS248 = crLHS12*crLHS247;
const double crLHS249 = DN(0,2)*N[3];
const double crLHS250 = DN(3,2)*N[0];
const double crLHS251 = gauss_weight*(N[0] + crLHS28*(1.0*crLHS22 + 1.0*crLHS24 + 1.0*crLHS25 + 1.0*crLHS26 + crLHS32));
const double crLHS252 = crLHS29*gauss_weight;
const double crLHS253 = crLHS252*(crLHS180 + crLHS229 + crLHS59);
const double crLHS254 = crLHS252*(crLHS195 + crLHS238 + crLHS96);
const double crLHS255 = crLHS252*(crLHS132 + crLHS210 + crLHS247);
const double crLHS256 = crLHS30*crLHS65;
const double crLHS257 = 1.0*crLHS68;
const double crLHS258 = crLHS257*crLHS28;
const double crLHS259 = crLHS35*crLHS66;
const double crLHS260 = crLHS17*crLHS66 + crLHS256*crLHS27 - crLHS258*crLHS27 + crLHS259*crLHS34;
const double crLHS261 = pow(DN(1,0), 2);
const double crLHS262 = pow(N[1], 2);
const double crLHS263 = crLHS20*crLHS262 + crLHS256*crLHS72 - crLHS258*crLHS72 + crLHS259*crLHS73 + crLHS262*crLHS6 + crLHS65*crLHS66;
const double crLHS264 = DN(1,0)*crLHS12;
const double crLHS265 = DN(1,1)*crLHS264;
const double crLHS266 = DN(1,2)*crLHS264;
const double crLHS267 = gauss_weight*(-N[1] + crLHS256 - crLHS258 + crLHS52*crLHS66);
const double crLHS268 = DN(1,0)*DN(2,0);
const double crLHS269 = N[2]*crLHS67 + N[2]*crLHS68;
const double crLHS270 = crLHS12*crLHS268 + crLHS269;
const double crLHS271 = crLHS102*crLHS66 + crLHS109*crLHS256 - crLHS109*crLHS258 + crLHS110*crLHS259;
const double crLHS272 = DN(2,1)*crLHS264;
const double crLHS273 = DN(2,2)*crLHS264;
const double crLHS274 = DN(1,0)*N[2];
const double crLHS275 = DN(2,0)*N[1];
const double crLHS276 = DN(1,0)*DN(3,0);
const double crLHS277 = N[3]*crLHS67 + N[3]*crLHS68;
const double crLHS278 = crLHS12*crLHS276 + crLHS277;
const double crLHS279 = crLHS138*crLHS66 + crLHS145*crLHS256 - crLHS145*crLHS258 + crLHS146*crLHS259;
const double crLHS280 = DN(3,1)*crLHS264;
const double crLHS281 = DN(3,2)*crLHS264;
const double crLHS282 = DN(1,0)*N[3];
const double crLHS283 = DN(3,0)*N[1];
const double crLHS284 = crLHS260 + crLHS60;
const double crLHS285 = pow(DN(1,1), 2);
const double crLHS286 = DN(1,1)*crLHS12;
const double crLHS287 = DN(1,2)*crLHS286;
const double crLHS288 = DN(2,0)*crLHS286;
const double crLHS289 = DN(1,1)*DN(2,1);
const double crLHS290 = crLHS12*crLHS289;
const double crLHS291 = crLHS269 + crLHS271;
const double crLHS292 = DN(2,2)*crLHS286;
const double crLHS293 = DN(1,1)*N[2];
const double crLHS294 = DN(2,1)*N[1];
const double crLHS295 = DN(3,0)*crLHS286;
const double crLHS296 = DN(1,1)*DN(3,1);
const double crLHS297 = crLHS12*crLHS296;
const double crLHS298 = crLHS277 + crLHS279;
const double crLHS299 = DN(3,2)*crLHS286;
const double crLHS300 = DN(1,1)*N[3];
const double crLHS301 = DN(3,1)*N[1];
const double crLHS302 = pow(DN(1,2), 2);
const double crLHS303 = DN(1,2)*crLHS12;
const double crLHS304 = DN(2,0)*crLHS303;
const double crLHS305 = DN(2,1)*crLHS303;
const double crLHS306 = DN(1,2)*DN(2,2);
const double crLHS307 = crLHS12*crLHS306;
const double crLHS308 = DN(1,2)*N[2];
const double crLHS309 = DN(2,2)*N[1];
const double crLHS310 = DN(3,0)*crLHS303;
const double crLHS311 = DN(3,1)*crLHS303;
const double crLHS312 = DN(1,2)*DN(3,2);
const double crLHS313 = crLHS12*crLHS312;
const double crLHS314 = DN(1,2)*N[3];
const double crLHS315 = DN(3,2)*N[1];
const double crLHS316 = gauss_weight*(N[1] + crLHS28*(crLHS257 + 1.0*crLHS67 + 1.0*crLHS69 + 1.0*crLHS70 + 1.0*crLHS71));
const double crLHS317 = crLHS252*(crLHS268 + crLHS289 + crLHS306);
const double crLHS318 = crLHS252*(crLHS276 + crLHS296 + crLHS312);
const double crLHS319 = crLHS102*crLHS30;
const double crLHS320 = 1.0*crLHS105;
const double crLHS321 = crLHS28*crLHS320;
const double crLHS322 = crLHS103*crLHS35;
const double crLHS323 = crLHS103*crLHS17 + crLHS27*crLHS319 - crLHS27*crLHS321 + crLHS322*crLHS34;
const double crLHS324 = crLHS103*crLHS65 + crLHS319*crLHS72 - crLHS321*crLHS72 + crLHS322*crLHS73;
const double crLHS325 = pow(DN(2,0), 2);
const double crLHS326 = pow(N[2], 2);
const double crLHS327 = crLHS102*crLHS103 + crLHS109*crLHS319 - crLHS109*crLHS321 + crLHS110*crLHS322 + crLHS20*crLHS326 + crLHS326*crLHS6;
const double crLHS328 = DN(2,0)*crLHS12;
const double crLHS329 = DN(2,1)*crLHS328;
const double crLHS330 = DN(2,2)*crLHS328;
const double crLHS331 = gauss_weight*(-N[2] + crLHS103*crLHS52 + crLHS319 - crLHS321);
const double crLHS332 = DN(2,0)*DN(3,0);
const double crLHS333 = N[3]*crLHS104 + N[3]*crLHS105;
const double crLHS334 = crLHS12*crLHS332 + crLHS333;
const double crLHS335 = crLHS103*crLHS138 + crLHS145*crLHS319 - crLHS145*crLHS321 + crLHS146*crLHS322;
const double crLHS336 = DN(3,1)*crLHS328;
const double crLHS337 = DN(3,2)*crLHS328;
const double crLHS338 = DN(2,0)*N[3];
const double crLHS339 = DN(3,0)*N[2];
const double crLHS340 = crLHS323 + crLHS97;
const double crLHS341 = crLHS269 + crLHS324;
const double crLHS342 = pow(DN(2,1), 2);
const double crLHS343 = DN(2,1)*crLHS12;
const double crLHS344 = DN(2,2)*crLHS343;
const double crLHS345 = DN(3,0)*crLHS343;
const double crLHS346 = DN(2,1)*DN(3,1);
const double crLHS347 = crLHS12*crLHS346;
const double crLHS348 = crLHS333 + crLHS335;
const double crLHS349 = DN(3,2)*crLHS343;
const double crLHS350 = DN(2,1)*N[3];
const double crLHS351 = DN(3,1)*N[2];
const double crLHS352 = pow(DN(2,2), 2);
const double crLHS353 = DN(2,2)*crLHS12;
const double crLHS354 = DN(3,0)*crLHS353;
const double crLHS355 = DN(3,1)*crLHS353;
const double crLHS356 = DN(2,2)*DN(3,2);
const double crLHS357 = crLHS12*crLHS356;
const double crLHS358 = DN(2,2)*N[3];
const double crLHS359 = DN(3,2)*N[2];
const double crLHS360 = gauss_weight*(N[2] + crLHS28*(1.0*crLHS104 + 1.0*crLHS106 + 1.0*crLHS107 + 1.0*crLHS108 + crLHS320));
const double crLHS361 = crLHS252*(crLHS332 + crLHS346 + crLHS356);
const double crLHS362 = crLHS138*crLHS30;
const double crLHS363 = 1.0*crLHS141;
const double crLHS364 = crLHS28*crLHS363;
const double crLHS365 = crLHS139*crLHS35;
const double crLHS366 = crLHS139*crLHS17 + crLHS27*crLHS362 - crLHS27*crLHS364 + crLHS34*crLHS365;
const double crLHS367 = crLHS139*crLHS65 + crLHS362*crLHS72 - crLHS364*crLHS72 + crLHS365*crLHS73;
const double crLHS368 = crLHS102*crLHS139 + crLHS109*crLHS362 - crLHS109*crLHS364 + crLHS110*crLHS365;
const double crLHS369 = pow(DN(3,0), 2);
const double crLHS370 = pow(N[3], 2);
const double crLHS371 = crLHS138*crLHS139 + crLHS145*crLHS362 - crLHS145*crLHS364 + crLHS146*crLHS365 + crLHS20*crLHS370 + crLHS370*crLHS6;
const double crLHS372 = DN(3,0)*crLHS12;
const double crLHS373 = DN(3,1)*crLHS372;
const double crLHS374 = DN(3,2)*crLHS372;
const double crLHS375 = gauss_weight*(-N[3] + crLHS139*crLHS52 + crLHS362 - crLHS364);
const double crLHS376 = crLHS133 + crLHS366;
const double crLHS377 = crLHS277 + crLHS367;
const double crLHS378 = crLHS333 + crLHS368;
const double crLHS379 = pow(DN(3,1), 2);
const double crLHS380 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS381 = pow(DN(3,2), 2);
const double crLHS382 = gauss_weight*(N[3] + crLHS28*(1.0*crLHS140 + 1.0*crLHS142 + 1.0*crLHS143 + 1.0*crLHS144 + crLHS363));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS37);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS40 + DN(0,2)*crLHS43 + crLHS45);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS46 + DN(0,1)*crLHS48 + DN(0,2)*crLHS50 + crLHS51);
rLHS(0,3)+=DN(0,0)*crLHS53;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS56 + DN(0,2)*crLHS58 + crLHS61 + crLHS74);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS77 + DN(0,2)*crLHS80 + crLHS81);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS82 + DN(0,1)*crLHS84 + DN(0,2)*crLHS86 + crLHS87);
rLHS(0,7)+=-gauss_weight*(-DN(1,0)*crLHS31 + DN(1,0)*crLHS33 + crLHS88 - crLHS89*crLHS90);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS91 + DN(0,1)*crLHS93 + DN(0,2)*crLHS95 + crLHS111 + crLHS98);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS114 + DN(0,2)*crLHS117 + crLHS118);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS119 + DN(0,1)*crLHS121 + DN(0,2)*crLHS123 + crLHS124);
rLHS(0,11)+=-gauss_weight*(-DN(2,0)*crLHS31 + DN(2,0)*crLHS33 + crLHS125 - crLHS126*crLHS90);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS129 + DN(0,2)*crLHS131 + crLHS134 + crLHS147);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS148 + DN(0,1)*crLHS150 + DN(0,2)*crLHS153 + crLHS154);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS155 + DN(0,1)*crLHS157 + DN(0,2)*crLHS159 + crLHS160);
rLHS(0,15)+=-gauss_weight*(-DN(3,0)*crLHS31 + DN(3,0)*crLHS33 + crLHS161 - crLHS162*crLHS90);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS163 + DN(0,2)*crLHS164 + crLHS45);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS165 + DN(0,2)*crLHS167 + crLHS12*crLHS168 + crLHS37);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS48 + DN(0,1)*crLHS169 + DN(0,2)*crLHS171 + crLHS173);
rLHS(1,3)+=DN(0,1)*crLHS53;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS56 + DN(0,1)*crLHS174 + DN(0,2)*crLHS175 + crLHS176);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS177 + DN(0,2)*crLHS179 + crLHS181 + crLHS182);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS84 + DN(0,1)*crLHS183 + DN(0,2)*crLHS185 + crLHS186);
rLHS(1,7)+=-gauss_weight*(-DN(1,1)*crLHS31 + DN(1,1)*crLHS33 + crLHS187 - crLHS188*crLHS90);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS93 + DN(0,1)*crLHS189 + DN(0,2)*crLHS190 + crLHS191);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS192 + DN(0,2)*crLHS194 + crLHS196 + crLHS197);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS121 + DN(0,1)*crLHS198 + DN(0,2)*crLHS200 + crLHS201);
rLHS(1,11)+=-gauss_weight*(-DN(2,1)*crLHS31 + DN(2,1)*crLHS33 + crLHS202 - crLHS203*crLHS90);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS129 + DN(0,1)*crLHS204 + DN(0,2)*crLHS205 + crLHS206);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS150 + DN(0,1)*crLHS207 + DN(0,2)*crLHS209 + crLHS211 + crLHS212);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS157 + DN(0,1)*crLHS213 + DN(0,2)*crLHS215 + crLHS216);
rLHS(1,15)+=-gauss_weight*(-DN(3,1)*crLHS31 + DN(3,1)*crLHS33 + crLHS217 - crLHS218*crLHS90);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS164 + DN(0,2)*crLHS219 + crLHS51);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS167 + DN(0,2)*crLHS220 + crLHS173);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS50 + DN(0,1)*crLHS171 + DN(0,2)*crLHS221 + crLHS12*crLHS222 + crLHS37);
rLHS(2,3)+=DN(0,2)*crLHS53;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS58 + DN(0,1)*crLHS175 + DN(0,2)*crLHS223 + crLHS225);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS80 + DN(0,1)*crLHS179 + DN(0,2)*crLHS226 + crLHS227);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS86 + DN(0,1)*crLHS185 + DN(0,2)*crLHS228 + crLHS182 + crLHS230);
rLHS(2,7)+=-gauss_weight*(-DN(1,2)*crLHS31 + DN(1,2)*crLHS33 + crLHS231 - crLHS232*crLHS90);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS95 + DN(0,1)*crLHS190 + DN(0,2)*crLHS233 + crLHS234);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS117 + DN(0,1)*crLHS194 + DN(0,2)*crLHS235 + crLHS236);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS123 + DN(0,1)*crLHS200 + DN(0,2)*crLHS237 + crLHS197 + crLHS239);
rLHS(2,11)+=-gauss_weight*(-DN(2,2)*crLHS31 + DN(2,2)*crLHS33 + crLHS240 - crLHS241*crLHS90);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS131 + DN(0,1)*crLHS205 + DN(0,2)*crLHS242 + crLHS243);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS153 + DN(0,1)*crLHS209 + DN(0,2)*crLHS244 + crLHS245);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS159 + DN(0,1)*crLHS215 + DN(0,2)*crLHS246 + crLHS212 + crLHS248);
rLHS(2,15)+=-gauss_weight*(-DN(3,2)*crLHS31 + DN(3,2)*crLHS33 + crLHS249 - crLHS250*crLHS90);
rLHS(3,0)+=DN(0,0)*crLHS251;
rLHS(3,1)+=DN(0,1)*crLHS251;
rLHS(3,2)+=DN(0,2)*crLHS251;
rLHS(3,3)+=crLHS252*(crLHS168 + crLHS222 + crLHS5);
rLHS(3,4)+=gauss_weight*(DN(0,0)*crLHS73 + crLHS89);
rLHS(3,5)+=gauss_weight*(DN(0,1)*crLHS73 + crLHS188);
rLHS(3,6)+=gauss_weight*(DN(0,2)*crLHS73 + crLHS232);
rLHS(3,7)+=crLHS253;
rLHS(3,8)+=gauss_weight*(DN(0,0)*crLHS110 + crLHS126);
rLHS(3,9)+=gauss_weight*(DN(0,1)*crLHS110 + crLHS203);
rLHS(3,10)+=gauss_weight*(DN(0,2)*crLHS110 + crLHS241);
rLHS(3,11)+=crLHS254;
rLHS(3,12)+=gauss_weight*(DN(0,0)*crLHS146 + crLHS162);
rLHS(3,13)+=gauss_weight*(DN(0,1)*crLHS146 + crLHS218);
rLHS(3,14)+=gauss_weight*(DN(0,2)*crLHS146 + crLHS250);
rLHS(3,15)+=crLHS255;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS260 + crLHS61);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS40 + DN(1,2)*crLHS43 + crLHS176);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS46 + DN(1,1)*crLHS48 + DN(1,2)*crLHS50 + crLHS225);
rLHS(4,3)+=gauss_weight*(DN(0,0)*crLHS256 - DN(0,0)*crLHS258 + crLHS88*crLHS90 - crLHS89);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS56 + DN(1,2)*crLHS58 + crLHS12*crLHS261 + crLHS263);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS77 + DN(1,2)*crLHS80 + crLHS265);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS82 + DN(1,1)*crLHS84 + DN(1,2)*crLHS86 + crLHS266);
rLHS(4,7)+=DN(1,0)*crLHS267;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS91 + DN(1,1)*crLHS93 + DN(1,2)*crLHS95 + crLHS270 + crLHS271);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS114 + DN(1,2)*crLHS117 + crLHS272);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS119 + DN(1,1)*crLHS121 + DN(1,2)*crLHS123 + crLHS273);
rLHS(4,11)+=-gauss_weight*(-DN(2,0)*crLHS256 + DN(2,0)*crLHS258 + crLHS274 - crLHS275*crLHS90);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS129 + DN(1,2)*crLHS131 + crLHS278 + crLHS279);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS148 + DN(1,1)*crLHS150 + DN(1,2)*crLHS153 + crLHS280);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS155 + DN(1,1)*crLHS157 + DN(1,2)*crLHS159 + crLHS281);
rLHS(4,15)+=-gauss_weight*(-DN(3,0)*crLHS256 + DN(3,0)*crLHS258 + crLHS282 - crLHS283*crLHS90);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS163 + DN(1,2)*crLHS164 + crLHS81);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS165 + DN(1,2)*crLHS167 + crLHS181 + crLHS284);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS48 + DN(1,1)*crLHS169 + DN(1,2)*crLHS171 + crLHS227);
rLHS(5,3)+=gauss_weight*(DN(0,1)*crLHS256 - DN(0,1)*crLHS258 + crLHS187*crLHS90 - crLHS188);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS56 + DN(1,1)*crLHS174 + DN(1,2)*crLHS175 + crLHS265);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS177 + DN(1,2)*crLHS179 + crLHS12*crLHS285 + crLHS263);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS84 + DN(1,1)*crLHS183 + DN(1,2)*crLHS185 + crLHS287);
rLHS(5,7)+=DN(1,1)*crLHS267;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS93 + DN(1,1)*crLHS189 + DN(1,2)*crLHS190 + crLHS288);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS192 + DN(1,2)*crLHS194 + crLHS290 + crLHS291);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS121 + DN(1,1)*crLHS198 + DN(1,2)*crLHS200 + crLHS292);
rLHS(5,11)+=-gauss_weight*(-DN(2,1)*crLHS256 + DN(2,1)*crLHS258 + crLHS293 - crLHS294*crLHS90);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS129 + DN(1,1)*crLHS204 + DN(1,2)*crLHS205 + crLHS295);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS150 + DN(1,1)*crLHS207 + DN(1,2)*crLHS209 + crLHS297 + crLHS298);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS157 + DN(1,1)*crLHS213 + DN(1,2)*crLHS215 + crLHS299);
rLHS(5,15)+=-gauss_weight*(-DN(3,1)*crLHS256 + DN(3,1)*crLHS258 + crLHS300 - crLHS301*crLHS90);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS164 + DN(1,2)*crLHS219 + crLHS87);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS167 + DN(1,2)*crLHS220 + crLHS186);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS50 + DN(1,1)*crLHS171 + DN(1,2)*crLHS221 + crLHS230 + crLHS284);
rLHS(6,3)+=gauss_weight*(DN(0,2)*crLHS256 - DN(0,2)*crLHS258 + crLHS231*crLHS90 - crLHS232);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS58 + DN(1,1)*crLHS175 + DN(1,2)*crLHS223 + crLHS266);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS80 + DN(1,1)*crLHS179 + DN(1,2)*crLHS226 + crLHS287);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS86 + DN(1,1)*crLHS185 + DN(1,2)*crLHS228 + crLHS12*crLHS302 + crLHS263);
rLHS(6,7)+=DN(1,2)*crLHS267;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS95 + DN(1,1)*crLHS190 + DN(1,2)*crLHS233 + crLHS304);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS117 + DN(1,1)*crLHS194 + DN(1,2)*crLHS235 + crLHS305);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS123 + DN(1,1)*crLHS200 + DN(1,2)*crLHS237 + crLHS291 + crLHS307);
rLHS(6,11)+=-gauss_weight*(-DN(2,2)*crLHS256 + DN(2,2)*crLHS258 + crLHS308 - crLHS309*crLHS90);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS131 + DN(1,1)*crLHS205 + DN(1,2)*crLHS242 + crLHS310);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS153 + DN(1,1)*crLHS209 + DN(1,2)*crLHS244 + crLHS311);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS159 + DN(1,1)*crLHS215 + DN(1,2)*crLHS246 + crLHS298 + crLHS313);
rLHS(6,15)+=-gauss_weight*(-DN(3,2)*crLHS256 + DN(3,2)*crLHS258 + crLHS314 - crLHS315*crLHS90);
rLHS(7,0)+=gauss_weight*(DN(1,0)*crLHS34 + crLHS88);
rLHS(7,1)+=gauss_weight*(DN(1,1)*crLHS34 + crLHS187);
rLHS(7,2)+=gauss_weight*(DN(1,2)*crLHS34 + crLHS231);
rLHS(7,3)+=crLHS253;
rLHS(7,4)+=DN(1,0)*crLHS316;
rLHS(7,5)+=DN(1,1)*crLHS316;
rLHS(7,6)+=DN(1,2)*crLHS316;
rLHS(7,7)+=crLHS252*(crLHS261 + crLHS285 + crLHS302);
rLHS(7,8)+=gauss_weight*(DN(1,0)*crLHS110 + crLHS275);
rLHS(7,9)+=gauss_weight*(DN(1,1)*crLHS110 + crLHS294);
rLHS(7,10)+=gauss_weight*(DN(1,2)*crLHS110 + crLHS309);
rLHS(7,11)+=crLHS317;
rLHS(7,12)+=gauss_weight*(DN(1,0)*crLHS146 + crLHS283);
rLHS(7,13)+=gauss_weight*(DN(1,1)*crLHS146 + crLHS301);
rLHS(7,14)+=gauss_weight*(DN(1,2)*crLHS146 + crLHS315);
rLHS(7,15)+=crLHS318;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS323 + crLHS98);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS40 + DN(2,2)*crLHS43 + crLHS191);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS46 + DN(2,1)*crLHS48 + DN(2,2)*crLHS50 + crLHS234);
rLHS(8,3)+=gauss_weight*(DN(0,0)*crLHS319 - DN(0,0)*crLHS321 + crLHS125*crLHS90 - crLHS126);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS56 + DN(2,2)*crLHS58 + crLHS270 + crLHS324);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS77 + DN(2,2)*crLHS80 + crLHS288);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS82 + DN(2,1)*crLHS84 + DN(2,2)*crLHS86 + crLHS304);
rLHS(8,7)+=gauss_weight*(DN(1,0)*crLHS319 - DN(1,0)*crLHS321 + crLHS274*crLHS90 - crLHS275);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS91 + DN(2,1)*crLHS93 + DN(2,2)*crLHS95 + crLHS12*crLHS325 + crLHS327);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS114 + DN(2,2)*crLHS117 + crLHS329);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS119 + DN(2,1)*crLHS121 + DN(2,2)*crLHS123 + crLHS330);
rLHS(8,11)+=DN(2,0)*crLHS331;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS129 + DN(2,2)*crLHS131 + crLHS334 + crLHS335);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS148 + DN(2,1)*crLHS150 + DN(2,2)*crLHS153 + crLHS336);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS155 + DN(2,1)*crLHS157 + DN(2,2)*crLHS159 + crLHS337);
rLHS(8,15)+=-gauss_weight*(-DN(3,0)*crLHS319 + DN(3,0)*crLHS321 + crLHS338 - crLHS339*crLHS90);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS163 + DN(2,2)*crLHS164 + crLHS118);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS165 + DN(2,2)*crLHS167 + crLHS196 + crLHS340);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS48 + DN(2,1)*crLHS169 + DN(2,2)*crLHS171 + crLHS236);
rLHS(9,3)+=gauss_weight*(DN(0,1)*crLHS319 - DN(0,1)*crLHS321 + crLHS202*crLHS90 - crLHS203);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS56 + DN(2,1)*crLHS174 + DN(2,2)*crLHS175 + crLHS272);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS177 + DN(2,2)*crLHS179 + crLHS290 + crLHS341);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS84 + DN(2,1)*crLHS183 + DN(2,2)*crLHS185 + crLHS305);
rLHS(9,7)+=gauss_weight*(DN(1,1)*crLHS319 - DN(1,1)*crLHS321 + crLHS293*crLHS90 - crLHS294);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS93 + DN(2,1)*crLHS189 + DN(2,2)*crLHS190 + crLHS329);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS192 + DN(2,2)*crLHS194 + crLHS12*crLHS342 + crLHS327);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS121 + DN(2,1)*crLHS198 + DN(2,2)*crLHS200 + crLHS344);
rLHS(9,11)+=DN(2,1)*crLHS331;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS129 + DN(2,1)*crLHS204 + DN(2,2)*crLHS205 + crLHS345);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS150 + DN(2,1)*crLHS207 + DN(2,2)*crLHS209 + crLHS347 + crLHS348);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS157 + DN(2,1)*crLHS213 + DN(2,2)*crLHS215 + crLHS349);
rLHS(9,15)+=-gauss_weight*(-DN(3,1)*crLHS319 + DN(3,1)*crLHS321 + crLHS350 - crLHS351*crLHS90);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS164 + DN(2,2)*crLHS219 + crLHS124);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS167 + DN(2,2)*crLHS220 + crLHS201);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS50 + DN(2,1)*crLHS171 + DN(2,2)*crLHS221 + crLHS239 + crLHS340);
rLHS(10,3)+=gauss_weight*(DN(0,2)*crLHS319 - DN(0,2)*crLHS321 + crLHS240*crLHS90 - crLHS241);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS58 + DN(2,1)*crLHS175 + DN(2,2)*crLHS223 + crLHS273);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS80 + DN(2,1)*crLHS179 + DN(2,2)*crLHS226 + crLHS292);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS86 + DN(2,1)*crLHS185 + DN(2,2)*crLHS228 + crLHS307 + crLHS341);
rLHS(10,7)+=gauss_weight*(DN(1,2)*crLHS319 - DN(1,2)*crLHS321 + crLHS308*crLHS90 - crLHS309);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS95 + DN(2,1)*crLHS190 + DN(2,2)*crLHS233 + crLHS330);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS117 + DN(2,1)*crLHS194 + DN(2,2)*crLHS235 + crLHS344);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS123 + DN(2,1)*crLHS200 + DN(2,2)*crLHS237 + crLHS12*crLHS352 + crLHS327);
rLHS(10,11)+=DN(2,2)*crLHS331;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS131 + DN(2,1)*crLHS205 + DN(2,2)*crLHS242 + crLHS354);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS153 + DN(2,1)*crLHS209 + DN(2,2)*crLHS244 + crLHS355);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS159 + DN(2,1)*crLHS215 + DN(2,2)*crLHS246 + crLHS348 + crLHS357);
rLHS(10,15)+=-gauss_weight*(-DN(3,2)*crLHS319 + DN(3,2)*crLHS321 + crLHS358 - crLHS359*crLHS90);
rLHS(11,0)+=gauss_weight*(DN(2,0)*crLHS34 + crLHS125);
rLHS(11,1)+=gauss_weight*(DN(2,1)*crLHS34 + crLHS202);
rLHS(11,2)+=gauss_weight*(DN(2,2)*crLHS34 + crLHS240);
rLHS(11,3)+=crLHS254;
rLHS(11,4)+=gauss_weight*(DN(2,0)*crLHS73 + crLHS274);
rLHS(11,5)+=gauss_weight*(DN(2,1)*crLHS73 + crLHS293);
rLHS(11,6)+=gauss_weight*(DN(2,2)*crLHS73 + crLHS308);
rLHS(11,7)+=crLHS317;
rLHS(11,8)+=DN(2,0)*crLHS360;
rLHS(11,9)+=DN(2,1)*crLHS360;
rLHS(11,10)+=DN(2,2)*crLHS360;
rLHS(11,11)+=crLHS252*(crLHS325 + crLHS342 + crLHS352);
rLHS(11,12)+=gauss_weight*(DN(2,0)*crLHS146 + crLHS339);
rLHS(11,13)+=gauss_weight*(DN(2,1)*crLHS146 + crLHS351);
rLHS(11,14)+=gauss_weight*(DN(2,2)*crLHS146 + crLHS359);
rLHS(11,15)+=crLHS361;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS134 + crLHS366);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS38 + DN(3,1)*crLHS40 + DN(3,2)*crLHS43 + crLHS206);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS46 + DN(3,1)*crLHS48 + DN(3,2)*crLHS50 + crLHS243);
rLHS(12,3)+=gauss_weight*(DN(0,0)*crLHS362 - DN(0,0)*crLHS364 + crLHS161*crLHS90 - crLHS162);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS54 + DN(3,1)*crLHS56 + DN(3,2)*crLHS58 + crLHS278 + crLHS367);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS75 + DN(3,1)*crLHS77 + DN(3,2)*crLHS80 + crLHS295);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS82 + DN(3,1)*crLHS84 + DN(3,2)*crLHS86 + crLHS310);
rLHS(12,7)+=gauss_weight*(DN(1,0)*crLHS362 - DN(1,0)*crLHS364 + crLHS282*crLHS90 - crLHS283);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS91 + DN(3,1)*crLHS93 + DN(3,2)*crLHS95 + crLHS334 + crLHS368);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS112 + DN(3,1)*crLHS114 + DN(3,2)*crLHS117 + crLHS345);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS119 + DN(3,1)*crLHS121 + DN(3,2)*crLHS123 + crLHS354);
rLHS(12,11)+=gauss_weight*(DN(2,0)*crLHS362 - DN(2,0)*crLHS364 + crLHS338*crLHS90 - crLHS339);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS129 + DN(3,2)*crLHS131 + crLHS12*crLHS369 + crLHS371);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS148 + DN(3,1)*crLHS150 + DN(3,2)*crLHS153 + crLHS373);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS155 + DN(3,1)*crLHS157 + DN(3,2)*crLHS159 + crLHS374);
rLHS(12,15)+=DN(3,0)*crLHS375;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS163 + DN(3,2)*crLHS164 + crLHS154);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS40 + DN(3,1)*crLHS165 + DN(3,2)*crLHS167 + crLHS211 + crLHS376);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS48 + DN(3,1)*crLHS169 + DN(3,2)*crLHS171 + crLHS245);
rLHS(13,3)+=gauss_weight*(DN(0,1)*crLHS362 - DN(0,1)*crLHS364 + crLHS217*crLHS90 - crLHS218);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS56 + DN(3,1)*crLHS174 + DN(3,2)*crLHS175 + crLHS280);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS177 + DN(3,2)*crLHS179 + crLHS297 + crLHS377);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS84 + DN(3,1)*crLHS183 + DN(3,2)*crLHS185 + crLHS311);
rLHS(13,7)+=gauss_weight*(DN(1,1)*crLHS362 - DN(1,1)*crLHS364 + crLHS300*crLHS90 - crLHS301);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS93 + DN(3,1)*crLHS189 + DN(3,2)*crLHS190 + crLHS336);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS114 + DN(3,1)*crLHS192 + DN(3,2)*crLHS194 + crLHS347 + crLHS378);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS121 + DN(3,1)*crLHS198 + DN(3,2)*crLHS200 + crLHS355);
rLHS(13,11)+=gauss_weight*(DN(2,1)*crLHS362 - DN(2,1)*crLHS364 + crLHS350*crLHS90 - crLHS351);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS129 + DN(3,1)*crLHS204 + DN(3,2)*crLHS205 + crLHS373);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS150 + DN(3,1)*crLHS207 + DN(3,2)*crLHS209 + crLHS12*crLHS379 + crLHS371);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS157 + DN(3,1)*crLHS213 + DN(3,2)*crLHS215 + crLHS380);
rLHS(13,15)+=DN(3,1)*crLHS375;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS164 + DN(3,2)*crLHS219 + crLHS160);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS43 + DN(3,1)*crLHS167 + DN(3,2)*crLHS220 + crLHS216);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS50 + DN(3,1)*crLHS171 + DN(3,2)*crLHS221 + crLHS248 + crLHS376);
rLHS(14,3)+=gauss_weight*(DN(0,2)*crLHS362 - DN(0,2)*crLHS364 + crLHS249*crLHS90 - crLHS250);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS58 + DN(3,1)*crLHS175 + DN(3,2)*crLHS223 + crLHS281);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS80 + DN(3,1)*crLHS179 + DN(3,2)*crLHS226 + crLHS299);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS86 + DN(3,1)*crLHS185 + DN(3,2)*crLHS228 + crLHS313 + crLHS377);
rLHS(14,7)+=gauss_weight*(DN(1,2)*crLHS362 - DN(1,2)*crLHS364 + crLHS314*crLHS90 - crLHS315);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS95 + DN(3,1)*crLHS190 + DN(3,2)*crLHS233 + crLHS337);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS117 + DN(3,1)*crLHS194 + DN(3,2)*crLHS235 + crLHS349);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS123 + DN(3,1)*crLHS200 + DN(3,2)*crLHS237 + crLHS357 + crLHS378);
rLHS(14,11)+=gauss_weight*(DN(2,2)*crLHS362 - DN(2,2)*crLHS364 + crLHS358*crLHS90 - crLHS359);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS131 + DN(3,1)*crLHS205 + DN(3,2)*crLHS242 + crLHS374);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS153 + DN(3,1)*crLHS209 + DN(3,2)*crLHS244 + crLHS380);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS159 + DN(3,1)*crLHS215 + DN(3,2)*crLHS246 + crLHS12*crLHS381 + crLHS371);
rLHS(14,15)+=DN(3,2)*crLHS375;
rLHS(15,0)+=gauss_weight*(DN(3,0)*crLHS34 + crLHS161);
rLHS(15,1)+=gauss_weight*(DN(3,1)*crLHS34 + crLHS217);
rLHS(15,2)+=gauss_weight*(DN(3,2)*crLHS34 + crLHS249);
rLHS(15,3)+=crLHS255;
rLHS(15,4)+=gauss_weight*(DN(3,0)*crLHS73 + crLHS282);
rLHS(15,5)+=gauss_weight*(DN(3,1)*crLHS73 + crLHS300);
rLHS(15,6)+=gauss_weight*(DN(3,2)*crLHS73 + crLHS314);
rLHS(15,7)+=crLHS318;
rLHS(15,8)+=gauss_weight*(DN(3,0)*crLHS110 + crLHS338);
rLHS(15,9)+=gauss_weight*(DN(3,1)*crLHS110 + crLHS350);
rLHS(15,10)+=gauss_weight*(DN(3,2)*crLHS110 + crLHS358);
rLHS(15,11)+=crLHS361;
rLHS(15,12)+=DN(3,0)*crLHS382;
rLHS(15,13)+=DN(3,1)*crLHS382;
rLHS(15,14)+=DN(3,2)*crLHS382;
rLHS(15,15)+=crLHS252*(crLHS369 + crLHS379 + crLHS381);

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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
const double crRHS4 = N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
const double crRHS5 = N[0]*rho;
const double crRHS6 = crRHS5*time_coeff;
const double crRHS7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS9 = crRHS7*crRHS8;
const double crRHS10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS11 = crRHS10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
const double crRHS12 = crRHS11 + crRHS9;
const double crRHS13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS14 = crRHS13 + crRHS7;
const double crRHS15 = crRHS2*stab_c3;
const double crRHS16 = rho*stab_c2*sqrt(pow(crRHS10, 2) + pow(crRHS8, 2));
const double crRHS17 = crRHS14*(h*(crRHS15*h + crRHS16)/stab_c1 + mu);
const double crRHS18 = rho*time_coeff;
const double crRHS19 = crRHS18*crRHS4;
const double crRHS20 = 1.0/(crRHS15 + crRHS16/h + crRHS18*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double crRHS21 = crRHS20*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS1 + crRHS11*rho + crRHS19 + crRHS3 + crRHS9*rho);
const double crRHS22 = N[0]*crRHS2;
const double crRHS23 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crRHS24 = crRHS23*crRHS5;
const double crRHS25 = rho*(DN(0,0)*crRHS8 + DN(0,1)*crRHS10);
const double crRHS26 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS27 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crRHS28 = N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
const double crRHS29 = crRHS8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
const double crRHS30 = crRHS10*crRHS13;
const double crRHS31 = crRHS29 + crRHS30;
const double crRHS32 = crRHS18*crRHS28;
const double crRHS33 = crRHS20*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS26 + crRHS27 + crRHS29*rho + crRHS30*rho + crRHS32);
const double crRHS34 = N[1]*rho;
const double crRHS35 = N[1]*crRHS2;
const double crRHS36 = crRHS23*crRHS34;
const double crRHS37 = rho*(DN(1,0)*crRHS8 + DN(1,1)*crRHS10);
const double crRHS38 = N[2]*rho;
const double crRHS39 = N[2]*crRHS2;
const double crRHS40 = crRHS23*crRHS38;
const double crRHS41 = rho*(DN(2,0)*crRHS8 + DN(2,1)*crRHS10);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS17 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS1 + N[0]*crRHS3 + crRHS12*crRHS5 - crRHS21*crRHS22 + crRHS21*crRHS24 + crRHS21*crRHS25 + crRHS4*crRHS6);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS17 + DN(0,1)*stress[1] - N[0]*crRHS26 + N[0]*crRHS27 - crRHS22*crRHS33 + crRHS24*crRHS33 + crRHS25*crRHS33 + crRHS28*crRHS6 + crRHS31*crRHS5);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS21 + DN(0,1)*crRHS33 + N[0]*crRHS14);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS17 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS1 + N[1]*crRHS19 + N[1]*crRHS3 + crRHS12*crRHS34 - crRHS21*crRHS35 + crRHS21*crRHS36 + crRHS21*crRHS37);
rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS17 + DN(1,1)*stress[1] - N[1]*crRHS26 + N[1]*crRHS27 + N[1]*crRHS32 + crRHS31*crRHS34 - crRHS33*crRHS35 + crRHS33*crRHS36 + crRHS33*crRHS37);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS21 + DN(1,1)*crRHS33 + N[1]*crRHS14);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS17 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS1 + N[2]*crRHS19 + N[2]*crRHS3 + crRHS12*crRHS38 - crRHS21*crRHS39 + crRHS21*crRHS40 + crRHS21*crRHS41);
rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS17 + DN(2,1)*stress[1] - N[2]*crRHS26 + N[2]*crRHS27 + N[2]*crRHS32 + crRHS31*crRHS38 - crRHS33*crRHS39 + crRHS33*crRHS40 + crRHS33*crRHS41);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS21 + DN(2,1)*crRHS33 + N[2]*crRHS14);

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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
const double crRHS4 = N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
const double crRHS5 = N[0]*rho;
const double crRHS6 = crRHS5*time_coeff;
const double crRHS7 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS9 = crRHS7*crRHS8;
const double crRHS10 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS11 = crRHS10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0));
const double crRHS12 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS13 = crRHS12*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
const double crRHS14 = crRHS11 + crRHS13 + crRHS9;
const double crRHS15 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS16 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS17 = crRHS15 + crRHS16 + crRHS7;
const double crRHS18 = crRHS2*stab_c3;
const double crRHS19 = rho*stab_c2*sqrt(pow(crRHS10, 2) + pow(crRHS12, 2) + pow(crRHS8, 2));
const double crRHS20 = crRHS17*(h*(crRHS18*h + crRHS19)/stab_c1 + mu);
const double crRHS21 = rho*time_coeff;
const double crRHS22 = crRHS21*crRHS4;
const double crRHS23 = 1.0/(crRHS18 + crRHS19/h + crRHS21*dyn_tau/dt + mu*stab_c1/pow(h, 2));
const double crRHS24 = crRHS23*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crRHS1 + crRHS11*rho + crRHS13*rho + crRHS22 + crRHS3 + crRHS9*rho);
const double crRHS25 = N[0]*crRHS2;
const double crRHS26 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crRHS27 = crRHS26*crRHS5;
const double crRHS28 = rho*(DN(0,0)*crRHS8 + DN(0,1)*crRHS10 + DN(0,2)*crRHS12);
const double crRHS29 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS30 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crRHS31 = N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
const double crRHS32 = crRHS8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1));
const double crRHS33 = crRHS10*crRHS15;
const double crRHS34 = crRHS12*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
const double crRHS35 = crRHS32 + crRHS33 + crRHS34;
const double crRHS36 = crRHS21*crRHS31;
const double crRHS37 = crRHS23*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS29 + crRHS30 + crRHS32*rho + crRHS33*rho + crRHS34*rho + crRHS36);
const double crRHS38 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS39 = crRHS2*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crRHS40 = N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
const double crRHS41 = crRHS8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2));
const double crRHS42 = crRHS10*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
const double crRHS43 = crRHS12*crRHS16;
const double crRHS44 = crRHS41 + crRHS42 + crRHS43;
const double crRHS45 = crRHS21*crRHS40;
const double crRHS46 = crRHS23*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS38 + crRHS39 + crRHS41*rho + crRHS42*rho + crRHS43*rho + crRHS45);
const double crRHS47 = N[1]*rho;
const double crRHS48 = N[1]*crRHS2;
const double crRHS49 = crRHS26*crRHS47;
const double crRHS50 = rho*(DN(1,0)*crRHS8 + DN(1,1)*crRHS10 + DN(1,2)*crRHS12);
const double crRHS51 = N[2]*rho;
const double crRHS52 = N[2]*crRHS2;
const double crRHS53 = crRHS26*crRHS51;
const double crRHS54 = rho*(DN(2,0)*crRHS8 + DN(2,1)*crRHS10 + DN(2,2)*crRHS12);
const double crRHS55 = N[3]*rho;
const double crRHS56 = N[3]*crRHS2;
const double crRHS57 = crRHS26*crRHS55;
const double crRHS58 = rho*(DN(3,0)*crRHS8 + DN(3,1)*crRHS10 + DN(3,2)*crRHS12);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS20 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] - N[0]*crRHS1 + N[0]*crRHS3 + crRHS14*crRHS5 - crRHS24*crRHS25 + crRHS24*crRHS27 + crRHS24*crRHS28 + crRHS4*crRHS6);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS20 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS29 + N[0]*crRHS30 - crRHS25*crRHS37 + crRHS27*crRHS37 + crRHS28*crRHS37 + crRHS31*crRHS6 + crRHS35*crRHS5);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS20 + DN(0,2)*stress[2] - N[0]*crRHS38 + N[0]*crRHS39 - crRHS25*crRHS46 + crRHS27*crRHS46 + crRHS28*crRHS46 + crRHS40*crRHS6 + crRHS44*crRHS5);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS24 + DN(0,1)*crRHS37 + DN(0,2)*crRHS46 + N[0]*crRHS17);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS20 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] - N[1]*crRHS1 + N[1]*crRHS22 + N[1]*crRHS3 + crRHS14*crRHS47 - crRHS24*crRHS48 + crRHS24*crRHS49 + crRHS24*crRHS50);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS20 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS29 + N[1]*crRHS30 + N[1]*crRHS36 + crRHS35*crRHS47 - crRHS37*crRHS48 + crRHS37*crRHS49 + crRHS37*crRHS50);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS20 + DN(1,2)*stress[2] - N[1]*crRHS38 + N[1]*crRHS39 + N[1]*crRHS45 + crRHS44*crRHS47 - crRHS46*crRHS48 + crRHS46*crRHS49 + crRHS46*crRHS50);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS24 + DN(1,1)*crRHS37 + DN(1,2)*crRHS46 + N[1]*crRHS17);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS20 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] - N[2]*crRHS1 + N[2]*crRHS22 + N[2]*crRHS3 + crRHS14*crRHS51 - crRHS24*crRHS52 + crRHS24*crRHS53 + crRHS24*crRHS54);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS20 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS29 + N[2]*crRHS30 + N[2]*crRHS36 + crRHS35*crRHS51 - crRHS37*crRHS52 + crRHS37*crRHS53 + crRHS37*crRHS54);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS20 + DN(2,2)*stress[2] - N[2]*crRHS38 + N[2]*crRHS39 + N[2]*crRHS45 + crRHS44*crRHS51 - crRHS46*crRHS52 + crRHS46*crRHS53 + crRHS46*crRHS54);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS24 + DN(2,1)*crRHS37 + DN(2,2)*crRHS46 + N[2]*crRHS17);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS20 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] - N[3]*crRHS1 + N[3]*crRHS22 + N[3]*crRHS3 + crRHS14*crRHS55 - crRHS24*crRHS56 + crRHS24*crRHS57 + crRHS24*crRHS58);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS20 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS29 + N[3]*crRHS30 + N[3]*crRHS36 + crRHS35*crRHS55 - crRHS37*crRHS56 + crRHS37*crRHS57 + crRHS37*crRHS58);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS20 + DN(3,2)*stress[2] - N[3]*crRHS38 + N[3]*crRHS39 + N[3]*crRHS45 + crRHS44*crRHS55 - crRHS46*crRHS56 + crRHS46*crRHS57 + crRHS46*crRHS58);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS24 + DN(3,1)*crRHS37 + DN(3,2)*crRHS46 + N[3]*crRHS17);

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
    const double time_coeff = rData.TopOptTimeCoefficient;

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

    const BoundedMatrix<double,2,3> v_conv_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = pow(DN(0,0), 2);
const double crLHS4 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crLHS5 = crLHS4*stab_c3;
const double crLHS6 = N[0]*v_conv_ns(0,0) + N[1]*v_conv_ns(1,0) + N[2]*v_conv_ns(2,0);
const double crLHS7 = N[0]*v_conv_ns(0,1) + N[1]*v_conv_ns(1,1) + N[2]*v_conv_ns(2,1);
const double crLHS8 = rho*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2));
const double crLHS9 = DN(0,0)*v_conv_ns(0,0) + DN(1,0)*v_conv_ns(1,0) + DN(2,0)*v_conv_ns(2,0);
const double crLHS10 = DN(0,0)*v_conv_ns(0,1);
const double crLHS11 = DN(1,0)*v_conv_ns(1,1);
const double crLHS12 = DN(2,0)*v_conv_ns(2,1);
const double crLHS13 = crLHS10 + crLHS11 + crLHS12;
const double crLHS14 = DN(0,1)*v_conv_ns(0,0);
const double crLHS15 = DN(1,1)*v_conv_ns(1,0);
const double crLHS16 = DN(2,1)*v_conv_ns(2,0);
const double crLHS17 = crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = DN(0,1)*v_conv_ns(0,1) + DN(1,1)*v_conv_ns(1,1) + DN(2,1)*v_conv_ns(2,1);
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
const double crLHS42 = bdf0*time_coeff;
const double crLHS43 = crLHS22*crLHS42 + crLHS24*crLHS35 + crLHS41;
const double crLHS44 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS45 = C(1,2)*DN(0,1);
const double crLHS46 = C(2,2)*DN(0,0) + crLHS45;
const double crLHS47 = DN(0,0)*crLHS20;
const double crLHS48 = DN(0,1)*crLHS47;
const double crLHS49 = crLHS13*crLHS33;
const double crLHS50 = crLHS41*rho;
const double crLHS51 = pow(rho, 2);
const double crLHS52 = crLHS35*crLHS51;
const double crLHS53 = crLHS49*crLHS52;
const double crLHS54 = crLHS18*crLHS24;
const double crLHS55 = crLHS23 + crLHS25 - crLHS27 - crLHS29 + crLHS54;
const double crLHS56 = 1.0*crLHS10 + 1.0*crLHS11 + 1.0*crLHS12;
const double crLHS57 = crLHS24*crLHS32;
const double crLHS58 = crLHS56*crLHS57;
const double crLHS59 = DN(0,0)*N[0];
const double crLHS60 = DN(0,0)*crLHS33;
const double crLHS61 = DN(0,0)*crLHS9 + DN(0,1)*crLHS13;
const double crLHS62 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS63 = C(0,2)*DN(1,0);
const double crLHS64 = C(2,2)*DN(1,1) + crLHS63;
const double crLHS65 = N[1]*rho;
const double crLHS66 = crLHS65*crLHS9;
const double crLHS67 = N[1]*crLHS4;
const double crLHS68 = DN(1,0)*crLHS6;
const double crLHS69 = crLHS68*rho;
const double crLHS70 = DN(1,1)*crLHS7;
const double crLHS71 = crLHS70*rho;
const double crLHS72 = -crLHS67 + crLHS69 + crLHS71;
const double crLHS73 = -crLHS66 + crLHS72;
const double crLHS74 = crLHS33*crLHS73;
const double crLHS75 = crLHS37*crLHS65;
const double crLHS76 = -crLHS73*crLHS9 + crLHS75;
const double crLHS77 = N[1]*crLHS23;
const double crLHS78 = crLHS24*crLHS42;
const double crLHS79 = N[1]*crLHS78 + crLHS77;
const double crLHS80 = crLHS35*crLHS65 + crLHS79;
const double crLHS81 = DN(0,0)*DN(1,0);
const double crLHS82 = N[1]*crLHS25 + crLHS20*crLHS81;
const double crLHS83 = C(0,1)*DN(1,1) + crLHS63;
const double crLHS84 = C(1,2)*DN(1,1);
const double crLHS85 = C(2,2)*DN(1,0) + crLHS84;
const double crLHS86 = DN(1,1)*crLHS47;
const double crLHS87 = crLHS18*crLHS65;
const double crLHS88 = crLHS66 + crLHS67 - crLHS69 - crLHS71 + crLHS87;
const double crLHS89 = crLHS13*crLHS24;
const double crLHS90 = crLHS49*rho;
const double crLHS91 = N[1]*crLHS89 - crLHS77*crLHS90;
const double crLHS92 = DN(0,0)*N[1];
const double crLHS93 = DN(1,0)*crLHS33;
const double crLHS94 = DN(1,0)*crLHS9 + DN(1,1)*crLHS13;
const double crLHS95 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS96 = C(0,2)*DN(2,0);
const double crLHS97 = C(2,2)*DN(2,1) + crLHS96;
const double crLHS98 = N[2]*rho;
const double crLHS99 = crLHS9*crLHS98;
const double crLHS100 = N[2]*crLHS4;
const double crLHS101 = DN(2,0)*crLHS6;
const double crLHS102 = crLHS101*rho;
const double crLHS103 = DN(2,1)*crLHS7;
const double crLHS104 = crLHS103*rho;
const double crLHS105 = -crLHS100 + crLHS102 + crLHS104;
const double crLHS106 = crLHS105 - crLHS99;
const double crLHS107 = crLHS106*crLHS33;
const double crLHS108 = crLHS37*crLHS98;
const double crLHS109 = -crLHS106*crLHS9 + crLHS108;
const double crLHS110 = N[2]*crLHS23;
const double crLHS111 = N[2]*crLHS78 + crLHS110;
const double crLHS112 = crLHS111 + crLHS35*crLHS98;
const double crLHS113 = DN(0,0)*DN(2,0);
const double crLHS114 = N[2]*crLHS25 + crLHS113*crLHS20;
const double crLHS115 = C(0,1)*DN(2,1) + crLHS96;
const double crLHS116 = C(1,2)*DN(2,1);
const double crLHS117 = C(2,2)*DN(2,0) + crLHS116;
const double crLHS118 = DN(2,1)*crLHS47;
const double crLHS119 = crLHS18*crLHS98;
const double crLHS120 = crLHS100 - crLHS102 - crLHS104 + crLHS119 + crLHS99;
const double crLHS121 = N[2]*crLHS89 - crLHS110*crLHS90;
const double crLHS122 = DN(0,0)*N[2];
const double crLHS123 = DN(2,0)*crLHS33;
const double crLHS124 = DN(2,0)*crLHS9 + DN(2,1)*crLHS13;
const double crLHS125 = C(0,1)*DN(0,0) + crLHS45;
const double crLHS126 = crLHS17*crLHS33;
const double crLHS127 = crLHS126*crLHS52;
const double crLHS128 = 1.0*crLHS14 + 1.0*crLHS15 + 1.0*crLHS16;
const double crLHS129 = crLHS128*crLHS57;
const double crLHS130 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS131 = pow(DN(0,1), 2);
const double crLHS132 = crLHS30 - crLHS54;
const double crLHS133 = crLHS132*crLHS33;
const double crLHS134 = -crLHS132*crLHS18 + crLHS38;
const double crLHS135 = DN(0,1)*N[0];
const double crLHS136 = DN(0,1)*crLHS33;
const double crLHS137 = DN(0,0)*crLHS17 + DN(0,1)*crLHS18;
const double crLHS138 = C(0,1)*DN(1,0) + crLHS84;
const double crLHS139 = DN(0,1)*crLHS20;
const double crLHS140 = DN(1,0)*crLHS139;
const double crLHS141 = crLHS17*crLHS24;
const double crLHS142 = crLHS126*rho;
const double crLHS143 = N[1]*crLHS141 - crLHS142*crLHS77;
const double crLHS144 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS145 = crLHS72 - crLHS87;
const double crLHS146 = crLHS145*crLHS33;
const double crLHS147 = -crLHS145*crLHS18 + crLHS75;
const double crLHS148 = DN(0,1)*DN(1,1);
const double crLHS149 = N[1]*crLHS54 + crLHS148*crLHS20;
const double crLHS150 = DN(0,1)*N[1];
const double crLHS151 = DN(1,1)*crLHS33;
const double crLHS152 = DN(1,0)*crLHS17 + DN(1,1)*crLHS18;
const double crLHS153 = C(0,1)*DN(2,0) + crLHS116;
const double crLHS154 = DN(2,0)*crLHS139;
const double crLHS155 = N[2]*crLHS141 - crLHS110*crLHS142;
const double crLHS156 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS157 = crLHS105 - crLHS119;
const double crLHS158 = crLHS157*crLHS33;
const double crLHS159 = crLHS108 - crLHS157*crLHS18;
const double crLHS160 = DN(0,1)*DN(2,1);
const double crLHS161 = N[2]*crLHS54 + crLHS160*crLHS20;
const double crLHS162 = DN(0,1)*N[2];
const double crLHS163 = DN(2,1)*crLHS33;
const double crLHS164 = DN(2,0)*crLHS17 + DN(2,1)*crLHS18;
const double crLHS165 = crLHS33*gauss_weight;
const double crLHS166 = DN(1,0)*N[0];
const double crLHS167 = DN(1,1)*N[0];
const double crLHS168 = crLHS165*(crLHS148 + crLHS81);
const double crLHS169 = DN(2,0)*N[0];
const double crLHS170 = DN(2,1)*N[0];
const double crLHS171 = crLHS165*(crLHS113 + crLHS160);
const double crLHS172 = crLHS68 + crLHS70;
const double crLHS173 = crLHS172*rho;
const double crLHS174 = crLHS33*crLHS65;
const double crLHS175 = crLHS172*crLHS24 + crLHS79;
const double crLHS176 = crLHS32*crLHS65;
const double crLHS177 = crLHS176*crLHS56;
const double crLHS178 = crLHS172*crLHS51;
const double crLHS179 = crLHS178*crLHS49;
const double crLHS180 = pow(DN(1,0), 2);
const double crLHS181 = pow(N[1], 2);
const double crLHS182 = crLHS181*rho;
const double crLHS183 = crLHS181*crLHS4;
const double crLHS184 = crLHS172*crLHS65 + crLHS182*crLHS42 + crLHS183;
const double crLHS185 = DN(1,0)*crLHS20;
const double crLHS186 = DN(1,1)*crLHS185;
const double crLHS187 = DN(1,0)*N[1];
const double crLHS188 = N[2]*crLHS65;
const double crLHS189 = N[2]*crLHS67 + crLHS188*crLHS42;
const double crLHS190 = crLHS172*crLHS98 + crLHS189;
const double crLHS191 = DN(1,0)*DN(2,0);
const double crLHS192 = N[2]*crLHS66 + crLHS191*crLHS20;
const double crLHS193 = DN(2,1)*crLHS185;
const double crLHS194 = crLHS33*crLHS98;
const double crLHS195 = crLHS194*crLHS67;
const double crLHS196 = crLHS13*crLHS188 - crLHS13*crLHS195;
const double crLHS197 = DN(1,0)*N[2];
const double crLHS198 = crLHS128*crLHS176;
const double crLHS199 = crLHS126*crLHS178;
const double crLHS200 = pow(DN(1,1), 2);
const double crLHS201 = DN(1,1)*N[1];
const double crLHS202 = DN(2,0)*crLHS20;
const double crLHS203 = DN(1,1)*crLHS202;
const double crLHS204 = crLHS17*crLHS188 - crLHS17*crLHS195;
const double crLHS205 = DN(1,1)*DN(2,1);
const double crLHS206 = N[2]*crLHS87 + crLHS20*crLHS205;
const double crLHS207 = DN(1,1)*N[2];
const double crLHS208 = DN(2,0)*N[1];
const double crLHS209 = DN(2,1)*N[1];
const double crLHS210 = crLHS165*(crLHS191 + crLHS205);
const double crLHS211 = crLHS101 + crLHS103;
const double crLHS212 = crLHS211*rho;
const double crLHS213 = crLHS111 + crLHS211*crLHS24;
const double crLHS214 = crLHS32*crLHS98;
const double crLHS215 = crLHS214*crLHS56;
const double crLHS216 = crLHS211*crLHS51;
const double crLHS217 = crLHS216*crLHS49;
const double crLHS218 = crLHS189 + crLHS211*crLHS65;
const double crLHS219 = pow(DN(2,0), 2);
const double crLHS220 = pow(N[2], 2);
const double crLHS221 = crLHS220*rho;
const double crLHS222 = crLHS220*crLHS4;
const double crLHS223 = crLHS211*crLHS98 + crLHS221*crLHS42 + crLHS222;
const double crLHS224 = DN(2,1)*crLHS202;
const double crLHS225 = DN(2,0)*N[2];
const double crLHS226 = crLHS128*crLHS214;
const double crLHS227 = crLHS126*crLHS216;
const double crLHS228 = pow(DN(2,1), 2);
const double crLHS229 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS20*crLHS3 + crLHS22*crLHS9 + crLHS23*crLHS34 + crLHS34*crLHS36 - crLHS39*crLHS40 + crLHS43);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS44 + DN(0,1)*crLHS46 - N[0]*crLHS53 + crLHS13*crLHS22 + crLHS48 - crLHS49*crLHS50 - crLHS55*crLHS58);
rLHS(0,2)+=-gauss_weight*(crLHS23*crLHS60 + crLHS36*crLHS60 + crLHS40*crLHS61 + crLHS59);
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS64 + crLHS23*crLHS74 + crLHS36*crLHS74 - crLHS40*crLHS76 + crLHS80 + crLHS82);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS83 + DN(0,1)*crLHS85 - N[1]*crLHS53 - crLHS58*crLHS88 + crLHS86 + crLHS91);
rLHS(0,5)+=-gauss_weight*(crLHS23*crLHS93 + crLHS36*crLHS93 + crLHS40*crLHS94 + crLHS92);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS95 + DN(0,1)*crLHS97 + crLHS107*crLHS23 + crLHS107*crLHS36 - crLHS109*crLHS40 + crLHS112 + crLHS114);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS115 + DN(0,1)*crLHS117 - N[2]*crLHS53 + crLHS118 - crLHS120*crLHS58 + crLHS121);
rLHS(0,8)+=-gauss_weight*(crLHS122 + crLHS123*crLHS23 + crLHS123*crLHS36 + crLHS124*crLHS40);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS125 - N[0]*crLHS127 - crLHS126*crLHS50 - crLHS129*crLHS55 + crLHS17*crLHS22 + crLHS48);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS46 + DN(0,1)*crLHS130 + crLHS131*crLHS20 + crLHS133*crLHS23 + crLHS133*crLHS36 - crLHS134*crLHS40 + crLHS18*crLHS22 + crLHS43);
rLHS(1,2)+=-gauss_weight*(crLHS135 + crLHS136*crLHS23 + crLHS136*crLHS36 + crLHS137*crLHS40);
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS64 + DN(0,1)*crLHS138 - N[1]*crLHS127 - crLHS129*crLHS88 + crLHS140 + crLHS143);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS85 + DN(0,1)*crLHS144 + crLHS146*crLHS23 + crLHS146*crLHS36 - crLHS147*crLHS40 + crLHS149 + crLHS80);
rLHS(1,5)+=-gauss_weight*(crLHS150 + crLHS151*crLHS23 + crLHS151*crLHS36 + crLHS152*crLHS40);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS97 + DN(0,1)*crLHS153 - N[2]*crLHS127 - crLHS120*crLHS129 + crLHS154 + crLHS155);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS117 + DN(0,1)*crLHS156 + crLHS112 + crLHS158*crLHS23 + crLHS158*crLHS36 - crLHS159*crLHS40 + crLHS161);
rLHS(1,8)+=-gauss_weight*(crLHS162 + crLHS163*crLHS23 + crLHS163*crLHS36 + crLHS164*crLHS40);
rLHS(2,0)+=gauss_weight*(crLHS135*crLHS142 - crLHS31*crLHS60 + crLHS59);
rLHS(2,1)+=gauss_weight*(-crLHS132*crLHS136 + crLHS135 + crLHS59*crLHS90);
rLHS(2,2)+=crLHS165*(crLHS131 + crLHS3);
rLHS(2,3)+=gauss_weight*(crLHS142*crLHS150 + crLHS166 - crLHS60*crLHS73);
rLHS(2,4)+=gauss_weight*(-crLHS136*crLHS145 + crLHS167 + crLHS90*crLHS92);
rLHS(2,5)+=crLHS168;
rLHS(2,6)+=gauss_weight*(-crLHS106*crLHS60 + crLHS142*crLHS162 + crLHS169);
rLHS(2,7)+=gauss_weight*(crLHS122*crLHS90 - crLHS136*crLHS157 + crLHS170);
rLHS(2,8)+=crLHS171;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS173*crLHS34 - crLHS174*crLHS39 + crLHS175 + crLHS34*crLHS67 + crLHS82);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS44 + DN(1,1)*crLHS46 - N[0]*crLHS179 + crLHS140 - crLHS177*crLHS55 + crLHS91);
rLHS(3,2)+=-gauss_weight*(crLHS166 + crLHS173*crLHS60 + crLHS174*crLHS61 + crLHS60*crLHS67);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS64 + crLHS173*crLHS74 - crLHS174*crLHS76 + crLHS180*crLHS20 + crLHS182*crLHS9 + crLHS184 + crLHS67*crLHS74);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS83 + DN(1,1)*crLHS85 - N[1]*crLHS179 + crLHS13*crLHS182 - crLHS177*crLHS88 - crLHS183*crLHS90 + crLHS186);
rLHS(3,5)+=-gauss_weight*(crLHS173*crLHS93 + crLHS174*crLHS94 + crLHS187 + crLHS67*crLHS93);
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS95 + DN(1,1)*crLHS97 + crLHS107*crLHS173 + crLHS107*crLHS67 - crLHS109*crLHS174 + crLHS190 + crLHS192);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS115 + DN(1,1)*crLHS117 - N[2]*crLHS179 - crLHS120*crLHS177 + crLHS193 + crLHS196);
rLHS(3,8)+=-gauss_weight*(crLHS123*crLHS173 + crLHS123*crLHS67 + crLHS124*crLHS174 + crLHS197);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS125 - N[0]*crLHS199 + crLHS143 - crLHS198*crLHS55 + crLHS86);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS46 + DN(1,1)*crLHS130 + crLHS133*crLHS173 + crLHS133*crLHS67 - crLHS134*crLHS174 + crLHS149 + crLHS175);
rLHS(4,2)+=-gauss_weight*(crLHS136*crLHS173 + crLHS136*crLHS67 + crLHS137*crLHS174 + crLHS167);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS64 + DN(1,1)*crLHS138 - N[1]*crLHS199 - crLHS142*crLHS183 + crLHS17*crLHS182 + crLHS186 - crLHS198*crLHS88);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS85 + DN(1,1)*crLHS144 + crLHS146*crLHS173 + crLHS146*crLHS67 - crLHS147*crLHS174 + crLHS18*crLHS182 + crLHS184 + crLHS20*crLHS200);
rLHS(4,5)+=-gauss_weight*(crLHS151*crLHS173 + crLHS151*crLHS67 + crLHS152*crLHS174 + crLHS201);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS97 + DN(1,1)*crLHS153 - N[2]*crLHS199 - crLHS120*crLHS198 + crLHS203 + crLHS204);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS117 + DN(1,1)*crLHS156 + crLHS158*crLHS173 + crLHS158*crLHS67 - crLHS159*crLHS174 + crLHS190 + crLHS206);
rLHS(4,8)+=-gauss_weight*(crLHS163*crLHS173 + crLHS163*crLHS67 + crLHS164*crLHS174 + crLHS207);
rLHS(5,0)+=gauss_weight*(crLHS142*crLHS167 - crLHS31*crLHS93 + crLHS92);
rLHS(5,1)+=gauss_weight*(-crLHS132*crLHS151 + crLHS150 + crLHS166*crLHS90);
rLHS(5,2)+=crLHS168;
rLHS(5,3)+=gauss_weight*(crLHS142*crLHS201 + crLHS187 - crLHS73*crLHS93);
rLHS(5,4)+=gauss_weight*(-crLHS145*crLHS151 + crLHS187*crLHS90 + crLHS201);
rLHS(5,5)+=crLHS165*(crLHS180 + crLHS200);
rLHS(5,6)+=gauss_weight*(-crLHS106*crLHS93 + crLHS142*crLHS207 + crLHS208);
rLHS(5,7)+=gauss_weight*(-crLHS151*crLHS157 + crLHS197*crLHS90 + crLHS209);
rLHS(5,8)+=crLHS210;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS100*crLHS34 + crLHS114 - crLHS194*crLHS39 + crLHS212*crLHS34 + crLHS213);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS44 + DN(2,1)*crLHS46 - N[0]*crLHS217 + crLHS121 + crLHS154 - crLHS215*crLHS55);
rLHS(6,2)+=-gauss_weight*(crLHS100*crLHS60 + crLHS169 + crLHS194*crLHS61 + crLHS212*crLHS60);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS64 + crLHS100*crLHS74 + crLHS192 - crLHS194*crLHS76 + crLHS212*crLHS74 + crLHS218);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS83 + DN(2,1)*crLHS85 - N[1]*crLHS217 + crLHS196 + crLHS203 - crLHS215*crLHS88);
rLHS(6,5)+=-gauss_weight*(crLHS100*crLHS93 + crLHS194*crLHS94 + crLHS208 + crLHS212*crLHS93);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS95 + DN(2,1)*crLHS97 + crLHS100*crLHS107 + crLHS107*crLHS212 - crLHS109*crLHS194 + crLHS20*crLHS219 + crLHS221*crLHS9 + crLHS223);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS115 + DN(2,1)*crLHS117 - N[2]*crLHS217 - crLHS120*crLHS215 + crLHS13*crLHS221 - crLHS222*crLHS90 + crLHS224);
rLHS(6,8)+=-gauss_weight*(crLHS100*crLHS123 + crLHS123*crLHS212 + crLHS124*crLHS194 + crLHS225);
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS125 - N[0]*crLHS227 + crLHS118 + crLHS155 - crLHS226*crLHS55);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS46 + DN(2,1)*crLHS130 + crLHS100*crLHS133 + crLHS133*crLHS212 - crLHS134*crLHS194 + crLHS161 + crLHS213);
rLHS(7,2)+=-gauss_weight*(crLHS100*crLHS136 + crLHS136*crLHS212 + crLHS137*crLHS194 + crLHS170);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS64 + DN(2,1)*crLHS138 - N[1]*crLHS227 + crLHS193 + crLHS204 - crLHS226*crLHS88);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS85 + DN(2,1)*crLHS144 + crLHS100*crLHS146 + crLHS146*crLHS212 - crLHS147*crLHS194 + crLHS206 + crLHS218);
rLHS(7,5)+=-gauss_weight*(crLHS100*crLHS151 + crLHS151*crLHS212 + crLHS152*crLHS194 + crLHS209);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS97 + DN(2,1)*crLHS153 - N[2]*crLHS227 - crLHS120*crLHS226 - crLHS142*crLHS222 + crLHS17*crLHS221 + crLHS224);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS117 + DN(2,1)*crLHS156 + crLHS100*crLHS158 + crLHS158*crLHS212 - crLHS159*crLHS194 + crLHS18*crLHS221 + crLHS20*crLHS228 + crLHS223);
rLHS(7,8)+=-gauss_weight*(crLHS100*crLHS163 + crLHS163*crLHS212 + crLHS164*crLHS194 + crLHS229);
rLHS(8,0)+=gauss_weight*(crLHS122 - crLHS123*crLHS31 + crLHS142*crLHS170);
rLHS(8,1)+=gauss_weight*(-crLHS132*crLHS163 + crLHS162 + crLHS169*crLHS90);
rLHS(8,2)+=crLHS171;
rLHS(8,3)+=gauss_weight*(-crLHS123*crLHS73 + crLHS142*crLHS209 + crLHS197);
rLHS(8,4)+=gauss_weight*(-crLHS145*crLHS163 + crLHS207 + crLHS208*crLHS90);
rLHS(8,5)+=crLHS210;
rLHS(8,6)+=gauss_weight*(-crLHS106*crLHS123 + crLHS142*crLHS229 + crLHS225);
rLHS(8,7)+=gauss_weight*(-crLHS157*crLHS163 + crLHS225*crLHS90 + crLHS229);
rLHS(8,8)+=crLHS165*(crLHS219 + crLHS228);
    
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
    const double time_coeff = rData.TopOptTimeCoefficient;

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

    const BoundedMatrix<double,3,4> v_conv_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

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
const double crLHS8 = N[0]*v_conv_ns(0,0) + N[1]*v_conv_ns(1,0) + N[2]*v_conv_ns(2,0) + N[3]*v_conv_ns(3,0);
const double crLHS9 = N[0]*v_conv_ns(0,1) + N[1]*v_conv_ns(1,1) + N[2]*v_conv_ns(2,1) + N[3]*v_conv_ns(3,1);
const double crLHS10 = N[0]*v_conv_ns(0,2) + N[1]*v_conv_ns(1,2) + N[2]*v_conv_ns(2,2) + N[3]*v_conv_ns(3,2);
const double crLHS11 = rho*stab_c2*sqrt(pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2));
const double crLHS12 = DN(0,0)*v_conv_ns(0,0) + DN(1,0)*v_conv_ns(1,0) + DN(2,0)*v_conv_ns(2,0) + DN(3,0)*v_conv_ns(3,0);
const double crLHS13 = DN(0,0)*v_conv_ns(0,1) + DN(1,0)*v_conv_ns(1,1) + DN(2,0)*v_conv_ns(2,1) + DN(3,0)*v_conv_ns(3,1);
const double crLHS14 = DN(0,0)*v_conv_ns(0,2) + DN(1,0)*v_conv_ns(1,2) + DN(2,0)*v_conv_ns(2,2) + DN(3,0)*v_conv_ns(3,2);
const double crLHS15 = DN(0,1)*v_conv_ns(0,0) + DN(1,1)*v_conv_ns(1,0) + DN(2,1)*v_conv_ns(2,0) + DN(3,1)*v_conv_ns(3,0);
const double crLHS16 = DN(0,1)*v_conv_ns(0,1) + DN(1,1)*v_conv_ns(1,1) + DN(2,1)*v_conv_ns(2,1) + DN(3,1)*v_conv_ns(3,1);
const double crLHS17 = DN(0,1)*v_conv_ns(0,2) + DN(1,1)*v_conv_ns(1,2) + DN(2,1)*v_conv_ns(2,2) + DN(3,1)*v_conv_ns(3,2);
const double crLHS18 = DN(0,2)*v_conv_ns(0,0) + DN(1,2)*v_conv_ns(1,0) + DN(2,2)*v_conv_ns(2,0) + DN(3,2)*v_conv_ns(3,0);
const double crLHS19 = DN(0,2)*v_conv_ns(0,1) + DN(1,2)*v_conv_ns(1,1) + DN(2,2)*v_conv_ns(2,1) + DN(3,2)*v_conv_ns(3,1);
const double crLHS20 = DN(0,2)*v_conv_ns(0,2) + DN(1,2)*v_conv_ns(1,2) + DN(2,2)*v_conv_ns(2,2) + DN(3,2)*v_conv_ns(3,2);
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
const double crLHS44 = bdf0*time_coeff;
const double crLHS45 = crLHS24*crLHS44 + crLHS26*crLHS35 + crLHS43;
const double crLHS46 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS47 = C(1,3)*DN(0,1);
const double crLHS48 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS47;
const double crLHS49 = C(3,5)*DN(0,0);
const double crLHS50 = C(4,5)*DN(0,2);
const double crLHS51 = C(1,5)*DN(0,1) + crLHS49 + crLHS50;
const double crLHS52 = DN(0,0)*crLHS22;
const double crLHS53 = DN(0,1)*crLHS52;
const double crLHS54 = crLHS13*crLHS33;
const double crLHS55 = crLHS43*rho;
const double crLHS56 = pow(rho, 2);
const double crLHS57 = crLHS35*crLHS56;
const double crLHS58 = N[0]*crLHS57;
const double crLHS59 = crLHS16*crLHS26;
const double crLHS60 = crLHS31 - crLHS59;
const double crLHS61 = crLHS13*crLHS27 - crLHS13*crLHS60 + crLHS19*crLHS39;
const double crLHS62 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS63 = C(3,4)*DN(0,1);
const double crLHS64 = C(2,3)*DN(0,2) + crLHS49 + crLHS63;
const double crLHS65 = C(2,5)*DN(0,2);
const double crLHS66 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS65;
const double crLHS67 = DN(0,2)*crLHS52;
const double crLHS68 = crLHS14*crLHS33;
const double crLHS69 = crLHS20*crLHS26;
const double crLHS70 = crLHS31 - crLHS69;
const double crLHS71 = crLHS14*crLHS27 - crLHS14*crLHS70 + crLHS17*crLHS37;
const double crLHS72 = DN(0,0)*N[0];
const double crLHS73 = DN(0,0)*crLHS33;
const double crLHS74 = DN(0,0)*crLHS12 + DN(0,1)*crLHS13 + DN(0,2)*crLHS14;
const double crLHS75 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS76 = C(0,3)*DN(1,0);
const double crLHS77 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS76;
const double crLHS78 = C(0,5)*DN(1,0);
const double crLHS79 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS78;
const double crLHS80 = N[1]*rho;
const double crLHS81 = crLHS12*crLHS80;
const double crLHS82 = N[1]*crLHS6;
const double crLHS83 = DN(1,0)*crLHS8;
const double crLHS84 = DN(1,1)*crLHS9;
const double crLHS85 = DN(1,2)*crLHS10;
const double crLHS86 = -crLHS82 + crLHS83*rho + crLHS84*rho + crLHS85*rho;
const double crLHS87 = -crLHS81 + crLHS86;
const double crLHS88 = crLHS33*crLHS87;
const double crLHS89 = crLHS13*crLHS80;
const double crLHS90 = crLHS15*crLHS89;
const double crLHS91 = crLHS14*crLHS80;
const double crLHS92 = crLHS18*crLHS91;
const double crLHS93 = -crLHS12*crLHS87 + crLHS90 + crLHS92;
const double crLHS94 = N[1]*crLHS25;
const double crLHS95 = N[1]*crLHS26;
const double crLHS96 = crLHS44*crLHS95 + crLHS94;
const double crLHS97 = crLHS35*crLHS80 + crLHS96;
const double crLHS98 = DN(0,0)*DN(1,0);
const double crLHS99 = N[1]*crLHS27 + crLHS22*crLHS98;
const double crLHS100 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS76;
const double crLHS101 = C(1,3)*DN(1,1);
const double crLHS102 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS101;
const double crLHS103 = C(3,5)*DN(1,0);
const double crLHS104 = C(4,5)*DN(1,2);
const double crLHS105 = C(1,5)*DN(1,1) + crLHS103 + crLHS104;
const double crLHS106 = DN(1,1)*crLHS52;
const double crLHS107 = crLHS16*crLHS80;
const double crLHS108 = -crLHS107 + crLHS86;
const double crLHS109 = -crLHS108*crLHS13 + crLHS13*crLHS81 + crLHS19*crLHS91;
const double crLHS110 = N[1]*crLHS57;
const double crLHS111 = crLHS94*rho;
const double crLHS112 = N[1]*crLHS37 - crLHS111*crLHS54;
const double crLHS113 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS78;
const double crLHS114 = C(3,4)*DN(1,1);
const double crLHS115 = C(2,3)*DN(1,2) + crLHS103 + crLHS114;
const double crLHS116 = C(2,5)*DN(1,2);
const double crLHS117 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS116;
const double crLHS118 = DN(1,2)*crLHS52;
const double crLHS119 = crLHS20*crLHS80;
const double crLHS120 = -crLHS119 + crLHS86;
const double crLHS121 = -crLHS120*crLHS14 + crLHS14*crLHS81 + crLHS17*crLHS89;
const double crLHS122 = N[1]*crLHS39 - crLHS111*crLHS68;
const double crLHS123 = DN(0,0)*N[1];
const double crLHS124 = DN(1,0)*crLHS33;
const double crLHS125 = DN(1,0)*crLHS12 + DN(1,1)*crLHS13 + DN(1,2)*crLHS14;
const double crLHS126 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS127 = C(0,3)*DN(2,0);
const double crLHS128 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS127;
const double crLHS129 = C(0,5)*DN(2,0);
const double crLHS130 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS129;
const double crLHS131 = N[2]*rho;
const double crLHS132 = crLHS12*crLHS131;
const double crLHS133 = N[2]*crLHS6;
const double crLHS134 = DN(2,0)*crLHS8;
const double crLHS135 = DN(2,1)*crLHS9;
const double crLHS136 = DN(2,2)*crLHS10;
const double crLHS137 = -crLHS133 + crLHS134*rho + crLHS135*rho + crLHS136*rho;
const double crLHS138 = -crLHS132 + crLHS137;
const double crLHS139 = crLHS138*crLHS33;
const double crLHS140 = crLHS13*crLHS131;
const double crLHS141 = crLHS140*crLHS15;
const double crLHS142 = crLHS131*crLHS14;
const double crLHS143 = crLHS142*crLHS18;
const double crLHS144 = -crLHS12*crLHS138 + crLHS141 + crLHS143;
const double crLHS145 = N[2]*crLHS25;
const double crLHS146 = N[2]*crLHS26;
const double crLHS147 = crLHS145 + crLHS146*crLHS44;
const double crLHS148 = crLHS131*crLHS35 + crLHS147;
const double crLHS149 = DN(0,0)*DN(2,0);
const double crLHS150 = N[2]*crLHS27 + crLHS149*crLHS22;
const double crLHS151 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS127;
const double crLHS152 = C(1,3)*DN(2,1);
const double crLHS153 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS152;
const double crLHS154 = C(3,5)*DN(2,0);
const double crLHS155 = C(4,5)*DN(2,2);
const double crLHS156 = C(1,5)*DN(2,1) + crLHS154 + crLHS155;
const double crLHS157 = DN(2,1)*crLHS52;
const double crLHS158 = crLHS131*crLHS16;
const double crLHS159 = crLHS137 - crLHS158;
const double crLHS160 = crLHS13*crLHS132 - crLHS13*crLHS159 + crLHS142*crLHS19;
const double crLHS161 = N[2]*crLHS57;
const double crLHS162 = crLHS145*rho;
const double crLHS163 = N[2]*crLHS37 - crLHS162*crLHS54;
const double crLHS164 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS129;
const double crLHS165 = C(3,4)*DN(2,1);
const double crLHS166 = C(2,3)*DN(2,2) + crLHS154 + crLHS165;
const double crLHS167 = C(2,5)*DN(2,2);
const double crLHS168 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS167;
const double crLHS169 = DN(2,2)*crLHS52;
const double crLHS170 = crLHS131*crLHS20;
const double crLHS171 = crLHS137 - crLHS170;
const double crLHS172 = crLHS132*crLHS14 - crLHS14*crLHS171 + crLHS140*crLHS17;
const double crLHS173 = N[2]*crLHS39 - crLHS162*crLHS68;
const double crLHS174 = DN(0,0)*N[2];
const double crLHS175 = DN(2,0)*crLHS33;
const double crLHS176 = DN(2,0)*crLHS12 + DN(2,1)*crLHS13 + DN(2,2)*crLHS14;
const double crLHS177 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS178 = C(0,3)*DN(3,0);
const double crLHS179 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS178;
const double crLHS180 = C(0,5)*DN(3,0);
const double crLHS181 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS180;
const double crLHS182 = N[3]*rho;
const double crLHS183 = crLHS12*crLHS182;
const double crLHS184 = N[3]*crLHS6;
const double crLHS185 = DN(3,0)*crLHS8;
const double crLHS186 = DN(3,1)*crLHS9;
const double crLHS187 = DN(3,2)*crLHS10;
const double crLHS188 = -crLHS184 + crLHS185*rho + crLHS186*rho + crLHS187*rho;
const double crLHS189 = -crLHS183 + crLHS188;
const double crLHS190 = crLHS189*crLHS33;
const double crLHS191 = crLHS13*crLHS182;
const double crLHS192 = crLHS15*crLHS191;
const double crLHS193 = crLHS14*crLHS182;
const double crLHS194 = crLHS18*crLHS193;
const double crLHS195 = -crLHS12*crLHS189 + crLHS192 + crLHS194;
const double crLHS196 = N[3]*crLHS25;
const double crLHS197 = N[3]*crLHS26;
const double crLHS198 = crLHS196 + crLHS197*crLHS44;
const double crLHS199 = crLHS182*crLHS35 + crLHS198;
const double crLHS200 = DN(0,0)*DN(3,0);
const double crLHS201 = N[3]*crLHS27 + crLHS200*crLHS22;
const double crLHS202 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS178;
const double crLHS203 = C(1,3)*DN(3,1);
const double crLHS204 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS203;
const double crLHS205 = C(3,5)*DN(3,0);
const double crLHS206 = C(4,5)*DN(3,2);
const double crLHS207 = C(1,5)*DN(3,1) + crLHS205 + crLHS206;
const double crLHS208 = DN(3,1)*crLHS52;
const double crLHS209 = crLHS16*crLHS182;
const double crLHS210 = crLHS188 - crLHS209;
const double crLHS211 = crLHS13*crLHS183 - crLHS13*crLHS210 + crLHS19*crLHS193;
const double crLHS212 = N[3]*crLHS57;
const double crLHS213 = crLHS196*rho;
const double crLHS214 = N[3]*crLHS37 - crLHS213*crLHS54;
const double crLHS215 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS180;
const double crLHS216 = C(3,4)*DN(3,1);
const double crLHS217 = C(2,3)*DN(3,2) + crLHS205 + crLHS216;
const double crLHS218 = C(2,5)*DN(3,2);
const double crLHS219 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS218;
const double crLHS220 = DN(3,2)*crLHS52;
const double crLHS221 = crLHS182*crLHS20;
const double crLHS222 = crLHS188 - crLHS221;
const double crLHS223 = crLHS14*crLHS183 - crLHS14*crLHS222 + crLHS17*crLHS191;
const double crLHS224 = N[3]*crLHS39 - crLHS213*crLHS68;
const double crLHS225 = DN(0,0)*N[3];
const double crLHS226 = DN(3,0)*crLHS33;
const double crLHS227 = DN(3,0)*crLHS12 + DN(3,1)*crLHS13 + DN(3,2)*crLHS14;
const double crLHS228 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS47;
const double crLHS229 = C(0,4)*DN(0,0) + crLHS50 + crLHS63;
const double crLHS230 = crLHS15*crLHS33;
const double crLHS231 = crLHS17*crLHS26;
const double crLHS232 = -crLHS15*crLHS32 + crLHS15*crLHS59 + crLHS18*crLHS231;
const double crLHS233 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS234 = C(1,4)*DN(0,1);
const double crLHS235 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS234;
const double crLHS236 = pow(DN(0,1), 2);
const double crLHS237 = crLHS33*crLHS60;
const double crLHS238 = crLHS19*crLHS231;
const double crLHS239 = -crLHS16*crLHS60 + crLHS238 + crLHS38;
const double crLHS240 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS234;
const double crLHS241 = C(2,4)*DN(0,2);
const double crLHS242 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS241;
const double crLHS243 = DN(0,1)*crLHS22;
const double crLHS244 = DN(0,2)*crLHS243;
const double crLHS245 = crLHS17*crLHS33;
const double crLHS246 = crLHS15*crLHS39 + crLHS17*crLHS59 - crLHS17*crLHS70;
const double crLHS247 = DN(0,1)*N[0];
const double crLHS248 = DN(0,1)*crLHS33;
const double crLHS249 = DN(0,0)*crLHS15 + DN(0,1)*crLHS16 + DN(0,2)*crLHS17;
const double crLHS250 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS101;
const double crLHS251 = C(0,4)*DN(1,0) + crLHS104 + crLHS114;
const double crLHS252 = DN(1,0)*crLHS243;
const double crLHS253 = crLHS17*crLHS80;
const double crLHS254 = crLHS107*crLHS15 - crLHS15*crLHS87 + crLHS18*crLHS253;
const double crLHS255 = crLHS15*crLHS26;
const double crLHS256 = N[1]*crLHS255 - crLHS111*crLHS230;
const double crLHS257 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS258 = C(1,4)*DN(1,1);
const double crLHS259 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS258;
const double crLHS260 = crLHS108*crLHS33;
const double crLHS261 = crLHS19*crLHS253;
const double crLHS262 = -crLHS108*crLHS16 + crLHS261 + crLHS90;
const double crLHS263 = DN(0,1)*DN(1,1);
const double crLHS264 = N[1]*crLHS59 + crLHS22*crLHS263;
const double crLHS265 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS258;
const double crLHS266 = C(2,4)*DN(1,2);
const double crLHS267 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS266;
const double crLHS268 = DN(1,2)*crLHS243;
const double crLHS269 = crLHS107*crLHS17 - crLHS120*crLHS17 + crLHS15*crLHS91;
const double crLHS270 = N[1]*crLHS231 - crLHS111*crLHS245;
const double crLHS271 = DN(0,1)*N[1];
const double crLHS272 = DN(1,1)*crLHS33;
const double crLHS273 = DN(1,0)*crLHS15 + DN(1,1)*crLHS16 + DN(1,2)*crLHS17;
const double crLHS274 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS152;
const double crLHS275 = C(0,4)*DN(2,0) + crLHS155 + crLHS165;
const double crLHS276 = DN(2,0)*crLHS243;
const double crLHS277 = crLHS131*crLHS17;
const double crLHS278 = -crLHS138*crLHS15 + crLHS15*crLHS158 + crLHS18*crLHS277;
const double crLHS279 = N[2]*crLHS255 - crLHS162*crLHS230;
const double crLHS280 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS281 = C(1,4)*DN(2,1);
const double crLHS282 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS281;
const double crLHS283 = crLHS159*crLHS33;
const double crLHS284 = crLHS19*crLHS277;
const double crLHS285 = crLHS141 - crLHS159*crLHS16 + crLHS284;
const double crLHS286 = DN(0,1)*DN(2,1);
const double crLHS287 = N[2]*crLHS59 + crLHS22*crLHS286;
const double crLHS288 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS281;
const double crLHS289 = C(2,4)*DN(2,2);
const double crLHS290 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS289;
const double crLHS291 = DN(2,2)*crLHS243;
const double crLHS292 = crLHS142*crLHS15 + crLHS158*crLHS17 - crLHS17*crLHS171;
const double crLHS293 = N[2]*crLHS231 - crLHS162*crLHS245;
const double crLHS294 = DN(0,1)*N[2];
const double crLHS295 = DN(2,1)*crLHS33;
const double crLHS296 = DN(2,0)*crLHS15 + DN(2,1)*crLHS16 + DN(2,2)*crLHS17;
const double crLHS297 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS203;
const double crLHS298 = C(0,4)*DN(3,0) + crLHS206 + crLHS216;
const double crLHS299 = DN(3,0)*crLHS243;
const double crLHS300 = crLHS17*crLHS182;
const double crLHS301 = -crLHS15*crLHS189 + crLHS15*crLHS209 + crLHS18*crLHS300;
const double crLHS302 = N[3]*crLHS255 - crLHS213*crLHS230;
const double crLHS303 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS304 = C(1,4)*DN(3,1);
const double crLHS305 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS304;
const double crLHS306 = crLHS210*crLHS33;
const double crLHS307 = crLHS19*crLHS300;
const double crLHS308 = -crLHS16*crLHS210 + crLHS192 + crLHS307;
const double crLHS309 = DN(0,1)*DN(3,1);
const double crLHS310 = N[3]*crLHS59 + crLHS22*crLHS309;
const double crLHS311 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS304;
const double crLHS312 = C(2,4)*DN(3,2);
const double crLHS313 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS312;
const double crLHS314 = DN(3,2)*crLHS243;
const double crLHS315 = crLHS15*crLHS193 + crLHS17*crLHS209 - crLHS17*crLHS222;
const double crLHS316 = N[3]*crLHS231 - crLHS213*crLHS245;
const double crLHS317 = DN(0,1)*N[3];
const double crLHS318 = DN(3,1)*crLHS33;
const double crLHS319 = DN(3,0)*crLHS15 + DN(3,1)*crLHS16 + DN(3,2)*crLHS17;
const double crLHS320 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS65;
const double crLHS321 = crLHS18*crLHS33;
const double crLHS322 = -crLHS18*crLHS32 + crLHS18*crLHS69 + crLHS19*crLHS255;
const double crLHS323 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS241;
const double crLHS324 = crLHS19*crLHS33;
const double crLHS325 = crLHS18*crLHS37 - crLHS19*crLHS60 + crLHS19*crLHS69;
const double crLHS326 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS327 = pow(DN(0,2), 2);
const double crLHS328 = crLHS33*crLHS70;
const double crLHS329 = -crLHS20*crLHS70 + crLHS238 + crLHS40;
const double crLHS330 = DN(0,2)*N[0];
const double crLHS331 = DN(0,2)*crLHS33;
const double crLHS332 = DN(0,0)*crLHS18 + DN(0,1)*crLHS19 + DN(0,2)*crLHS20;
const double crLHS333 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS116;
const double crLHS334 = DN(0,2)*crLHS22;
const double crLHS335 = DN(1,0)*crLHS334;
const double crLHS336 = crLHS15*crLHS19;
const double crLHS337 = crLHS119*crLHS18 - crLHS18*crLHS87 + crLHS336*crLHS80;
const double crLHS338 = -crLHS111*crLHS321 + crLHS18*crLHS95;
const double crLHS339 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS266;
const double crLHS340 = DN(1,1)*crLHS334;
const double crLHS341 = -crLHS108*crLHS19 + crLHS119*crLHS19 + crLHS18*crLHS89;
const double crLHS342 = -crLHS111*crLHS324 + crLHS19*crLHS95;
const double crLHS343 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS344 = crLHS120*crLHS33;
const double crLHS345 = -crLHS120*crLHS20 + crLHS261 + crLHS92;
const double crLHS346 = DN(0,2)*DN(1,2);
const double crLHS347 = N[1]*crLHS69 + crLHS22*crLHS346;
const double crLHS348 = DN(0,2)*N[1];
const double crLHS349 = DN(1,2)*crLHS33;
const double crLHS350 = DN(1,0)*crLHS18 + DN(1,1)*crLHS19 + DN(1,2)*crLHS20;
const double crLHS351 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS167;
const double crLHS352 = DN(2,0)*crLHS334;
const double crLHS353 = crLHS131*crLHS336 - crLHS138*crLHS18 + crLHS170*crLHS18;
const double crLHS354 = crLHS146*crLHS18 - crLHS162*crLHS321;
const double crLHS355 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS289;
const double crLHS356 = DN(2,1)*crLHS334;
const double crLHS357 = crLHS140*crLHS18 - crLHS159*crLHS19 + crLHS170*crLHS19;
const double crLHS358 = crLHS146*crLHS19 - crLHS162*crLHS324;
const double crLHS359 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS360 = crLHS171*crLHS33;
const double crLHS361 = crLHS143 - crLHS171*crLHS20 + crLHS284;
const double crLHS362 = DN(0,2)*DN(2,2);
const double crLHS363 = N[2]*crLHS69 + crLHS22*crLHS362;
const double crLHS364 = DN(0,2)*N[2];
const double crLHS365 = DN(2,2)*crLHS33;
const double crLHS366 = DN(2,0)*crLHS18 + DN(2,1)*crLHS19 + DN(2,2)*crLHS20;
const double crLHS367 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS218;
const double crLHS368 = DN(3,0)*crLHS334;
const double crLHS369 = -crLHS18*crLHS189 + crLHS18*crLHS221 + crLHS182*crLHS336;
const double crLHS370 = crLHS18*crLHS197 - crLHS213*crLHS321;
const double crLHS371 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS312;
const double crLHS372 = DN(3,1)*crLHS334;
const double crLHS373 = crLHS18*crLHS191 - crLHS19*crLHS210 + crLHS19*crLHS221;
const double crLHS374 = crLHS19*crLHS197 - crLHS213*crLHS324;
const double crLHS375 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS376 = crLHS222*crLHS33;
const double crLHS377 = crLHS194 - crLHS20*crLHS222 + crLHS307;
const double crLHS378 = DN(0,2)*DN(3,2);
const double crLHS379 = N[3]*crLHS69 + crLHS22*crLHS378;
const double crLHS380 = DN(0,2)*N[3];
const double crLHS381 = DN(3,2)*crLHS33;
const double crLHS382 = DN(3,0)*crLHS18 + DN(3,1)*crLHS19 + DN(3,2)*crLHS20;
const double crLHS383 = crLHS247*rho;
const double crLHS384 = crLHS330*rho;
const double crLHS385 = crLHS72*rho;
const double crLHS386 = crLHS33*gauss_weight;
const double crLHS387 = DN(1,0)*N[0];
const double crLHS388 = crLHS271*rho;
const double crLHS389 = crLHS348*rho;
const double crLHS390 = DN(1,1)*N[0];
const double crLHS391 = crLHS123*rho;
const double crLHS392 = DN(1,2)*N[0];
const double crLHS393 = crLHS386*(crLHS263 + crLHS346 + crLHS98);
const double crLHS394 = DN(2,0)*N[0];
const double crLHS395 = crLHS294*rho;
const double crLHS396 = crLHS364*rho;
const double crLHS397 = DN(2,1)*N[0];
const double crLHS398 = crLHS174*rho;
const double crLHS399 = DN(2,2)*N[0];
const double crLHS400 = crLHS386*(crLHS149 + crLHS286 + crLHS362);
const double crLHS401 = DN(3,0)*N[0];
const double crLHS402 = crLHS317*rho;
const double crLHS403 = crLHS380*rho;
const double crLHS404 = DN(3,1)*N[0];
const double crLHS405 = crLHS225*rho;
const double crLHS406 = DN(3,2)*N[0];
const double crLHS407 = crLHS386*(crLHS200 + crLHS309 + crLHS378);
const double crLHS408 = crLHS83 + crLHS84 + crLHS85;
const double crLHS409 = crLHS408*rho;
const double crLHS410 = crLHS33*crLHS80;
const double crLHS411 = crLHS26*crLHS408 + crLHS96;
const double crLHS412 = crLHS408*crLHS56;
const double crLHS413 = N[0]*crLHS412;
const double crLHS414 = pow(DN(1,0), 2);
const double crLHS415 = pow(N[1], 2);
const double crLHS416 = crLHS415*rho;
const double crLHS417 = crLHS415*crLHS6;
const double crLHS418 = crLHS408*crLHS80 + crLHS416*crLHS44 + crLHS417;
const double crLHS419 = DN(1,0)*crLHS22;
const double crLHS420 = DN(1,1)*crLHS419;
const double crLHS421 = crLHS417*rho;
const double crLHS422 = N[1]*crLHS412;
const double crLHS423 = DN(1,2)*crLHS419;
const double crLHS424 = DN(1,0)*N[1];
const double crLHS425 = N[2]*crLHS80;
const double crLHS426 = N[2]*crLHS82 + crLHS425*crLHS44;
const double crLHS427 = crLHS131*crLHS408 + crLHS426;
const double crLHS428 = DN(1,0)*DN(2,0);
const double crLHS429 = N[2]*crLHS81 + crLHS22*crLHS428;
const double crLHS430 = DN(2,1)*crLHS419;
const double crLHS431 = N[2]*crLHS412;
const double crLHS432 = crLHS33*crLHS82;
const double crLHS433 = N[2]*crLHS89 - crLHS140*crLHS432;
const double crLHS434 = DN(2,2)*crLHS419;
const double crLHS435 = N[2]*crLHS91 - crLHS142*crLHS432;
const double crLHS436 = DN(1,0)*N[2];
const double crLHS437 = N[3]*crLHS80;
const double crLHS438 = N[3]*crLHS82 + crLHS437*crLHS44;
const double crLHS439 = crLHS182*crLHS408 + crLHS438;
const double crLHS440 = DN(1,0)*DN(3,0);
const double crLHS441 = N[3]*crLHS81 + crLHS22*crLHS440;
const double crLHS442 = DN(3,1)*crLHS419;
const double crLHS443 = N[3]*crLHS412;
const double crLHS444 = N[3]*crLHS89 - crLHS191*crLHS432;
const double crLHS445 = DN(3,2)*crLHS419;
const double crLHS446 = N[3]*crLHS91 - crLHS193*crLHS432;
const double crLHS447 = DN(1,0)*N[3];
const double crLHS448 = pow(DN(1,1), 2);
const double crLHS449 = DN(1,1)*crLHS22;
const double crLHS450 = DN(1,2)*crLHS449;
const double crLHS451 = DN(1,1)*N[1];
const double crLHS452 = DN(2,0)*crLHS449;
const double crLHS453 = crLHS15*crLHS80;
const double crLHS454 = crLHS131*crLHS33;
const double crLHS455 = crLHS15*crLHS82;
const double crLHS456 = N[2]*crLHS453 - crLHS454*crLHS455;
const double crLHS457 = DN(1,1)*DN(2,1);
const double crLHS458 = N[2]*crLHS107 + crLHS22*crLHS457;
const double crLHS459 = DN(2,2)*crLHS449;
const double crLHS460 = N[2]*crLHS253 - crLHS277*crLHS432;
const double crLHS461 = DN(1,1)*N[2];
const double crLHS462 = DN(3,0)*crLHS449;
const double crLHS463 = crLHS182*crLHS33;
const double crLHS464 = N[3]*crLHS453 - crLHS455*crLHS463;
const double crLHS465 = DN(1,1)*DN(3,1);
const double crLHS466 = N[3]*crLHS107 + crLHS22*crLHS465;
const double crLHS467 = DN(3,2)*crLHS449;
const double crLHS468 = N[3]*crLHS253 - crLHS300*crLHS432;
const double crLHS469 = DN(1,1)*N[3];
const double crLHS470 = pow(DN(1,2), 2);
const double crLHS471 = DN(1,2)*N[1];
const double crLHS472 = DN(1,2)*crLHS22;
const double crLHS473 = DN(2,0)*crLHS472;
const double crLHS474 = crLHS454*crLHS82;
const double crLHS475 = crLHS18*crLHS425 - crLHS18*crLHS474;
const double crLHS476 = DN(2,1)*crLHS472;
const double crLHS477 = crLHS19*crLHS425 - crLHS19*crLHS474;
const double crLHS478 = DN(1,2)*DN(2,2);
const double crLHS479 = N[2]*crLHS119 + crLHS22*crLHS478;
const double crLHS480 = DN(1,2)*N[2];
const double crLHS481 = DN(3,0)*crLHS472;
const double crLHS482 = crLHS463*crLHS82;
const double crLHS483 = crLHS18*crLHS437 - crLHS18*crLHS482;
const double crLHS484 = DN(3,1)*crLHS472;
const double crLHS485 = crLHS19*crLHS437 - crLHS19*crLHS482;
const double crLHS486 = DN(1,2)*DN(3,2);
const double crLHS487 = N[3]*crLHS119 + crLHS22*crLHS486;
const double crLHS488 = DN(1,2)*N[3];
const double crLHS489 = crLHS390*rho;
const double crLHS490 = crLHS392*rho;
const double crLHS491 = crLHS387*rho;
const double crLHS492 = crLHS451*rho;
const double crLHS493 = crLHS471*rho;
const double crLHS494 = crLHS424*rho;
const double crLHS495 = DN(2,0)*N[1];
const double crLHS496 = crLHS461*rho;
const double crLHS497 = crLHS480*rho;
const double crLHS498 = DN(2,1)*N[1];
const double crLHS499 = crLHS436*rho;
const double crLHS500 = DN(2,2)*N[1];
const double crLHS501 = crLHS386*(crLHS428 + crLHS457 + crLHS478);
const double crLHS502 = DN(3,0)*N[1];
const double crLHS503 = crLHS469*rho;
const double crLHS504 = crLHS488*rho;
const double crLHS505 = DN(3,1)*N[1];
const double crLHS506 = crLHS447*rho;
const double crLHS507 = DN(3,2)*N[1];
const double crLHS508 = crLHS386*(crLHS440 + crLHS465 + crLHS486);
const double crLHS509 = crLHS134 + crLHS135 + crLHS136;
const double crLHS510 = crLHS509*rho;
const double crLHS511 = crLHS147 + crLHS26*crLHS509;
const double crLHS512 = crLHS509*crLHS56;
const double crLHS513 = N[0]*crLHS512;
const double crLHS514 = crLHS426 + crLHS509*crLHS80;
const double crLHS515 = N[1]*crLHS512;
const double crLHS516 = pow(DN(2,0), 2);
const double crLHS517 = pow(N[2], 2);
const double crLHS518 = crLHS517*rho;
const double crLHS519 = crLHS517*crLHS6;
const double crLHS520 = crLHS131*crLHS509 + crLHS44*crLHS518 + crLHS519;
const double crLHS521 = DN(2,0)*crLHS22;
const double crLHS522 = DN(2,1)*crLHS521;
const double crLHS523 = crLHS519*rho;
const double crLHS524 = N[2]*crLHS512;
const double crLHS525 = DN(2,2)*crLHS521;
const double crLHS526 = DN(2,0)*N[2];
const double crLHS527 = N[3]*crLHS131;
const double crLHS528 = N[3]*crLHS133 + crLHS44*crLHS527;
const double crLHS529 = crLHS182*crLHS509 + crLHS528;
const double crLHS530 = DN(2,0)*DN(3,0);
const double crLHS531 = N[3]*crLHS132 + crLHS22*crLHS530;
const double crLHS532 = DN(3,1)*crLHS521;
const double crLHS533 = N[3]*crLHS512;
const double crLHS534 = crLHS133*crLHS33;
const double crLHS535 = N[3]*crLHS140 - crLHS191*crLHS534;
const double crLHS536 = DN(3,2)*crLHS521;
const double crLHS537 = N[3]*crLHS142 - crLHS193*crLHS534;
const double crLHS538 = DN(2,0)*N[3];
const double crLHS539 = pow(DN(2,1), 2);
const double crLHS540 = DN(2,1)*crLHS22;
const double crLHS541 = DN(2,2)*crLHS540;
const double crLHS542 = DN(2,1)*N[2];
const double crLHS543 = DN(3,0)*crLHS540;
const double crLHS544 = crLHS133*crLHS463;
const double crLHS545 = crLHS15*crLHS527 - crLHS15*crLHS544;
const double crLHS546 = DN(2,1)*DN(3,1);
const double crLHS547 = N[3]*crLHS158 + crLHS22*crLHS546;
const double crLHS548 = DN(3,2)*crLHS540;
const double crLHS549 = N[3]*crLHS277 - crLHS300*crLHS534;
const double crLHS550 = DN(2,1)*N[3];
const double crLHS551 = pow(DN(2,2), 2);
const double crLHS552 = DN(2,2)*N[2];
const double crLHS553 = DN(2,2)*crLHS22;
const double crLHS554 = DN(3,0)*crLHS553;
const double crLHS555 = crLHS18*crLHS527 - crLHS18*crLHS544;
const double crLHS556 = DN(3,1)*crLHS553;
const double crLHS557 = crLHS19*crLHS527 - crLHS19*crLHS544;
const double crLHS558 = DN(2,2)*DN(3,2);
const double crLHS559 = N[3]*crLHS170 + crLHS22*crLHS558;
const double crLHS560 = DN(2,2)*N[3];
const double crLHS561 = crLHS397*rho;
const double crLHS562 = crLHS399*rho;
const double crLHS563 = crLHS394*rho;
const double crLHS564 = crLHS498*rho;
const double crLHS565 = crLHS500*rho;
const double crLHS566 = crLHS495*rho;
const double crLHS567 = crLHS542*rho;
const double crLHS568 = crLHS552*rho;
const double crLHS569 = crLHS526*rho;
const double crLHS570 = DN(3,0)*N[2];
const double crLHS571 = crLHS550*rho;
const double crLHS572 = crLHS560*rho;
const double crLHS573 = DN(3,1)*N[2];
const double crLHS574 = crLHS538*rho;
const double crLHS575 = DN(3,2)*N[2];
const double crLHS576 = crLHS386*(crLHS530 + crLHS546 + crLHS558);
const double crLHS577 = crLHS185 + crLHS186 + crLHS187;
const double crLHS578 = crLHS577*rho;
const double crLHS579 = crLHS198 + crLHS26*crLHS577;
const double crLHS580 = crLHS56*crLHS577;
const double crLHS581 = N[0]*crLHS580;
const double crLHS582 = crLHS438 + crLHS577*crLHS80;
const double crLHS583 = N[1]*crLHS580;
const double crLHS584 = crLHS131*crLHS577 + crLHS528;
const double crLHS585 = N[2]*crLHS580;
const double crLHS586 = pow(DN(3,0), 2);
const double crLHS587 = pow(N[3], 2);
const double crLHS588 = crLHS587*rho;
const double crLHS589 = crLHS587*crLHS6;
const double crLHS590 = crLHS182*crLHS577 + crLHS44*crLHS588 + crLHS589;
const double crLHS591 = DN(3,0)*crLHS22;
const double crLHS592 = DN(3,1)*crLHS591;
const double crLHS593 = crLHS589*rho;
const double crLHS594 = N[3]*crLHS580;
const double crLHS595 = DN(3,2)*crLHS591;
const double crLHS596 = DN(3,0)*N[3];
const double crLHS597 = pow(DN(3,1), 2);
const double crLHS598 = DN(3,1)*DN(3,2)*crLHS22;
const double crLHS599 = DN(3,1)*N[3];
const double crLHS600 = pow(DN(3,2), 2);
const double crLHS601 = DN(3,2)*N[3];
const double crLHS602 = crLHS404*rho;
const double crLHS603 = crLHS406*rho;
const double crLHS604 = crLHS401*rho;
const double crLHS605 = crLHS505*rho;
const double crLHS606 = crLHS507*rho;
const double crLHS607 = crLHS502*rho;
const double crLHS608 = crLHS573*rho;
const double crLHS609 = crLHS575*rho;
const double crLHS610 = crLHS570*rho;
const double crLHS611 = crLHS599*rho;
const double crLHS612 = crLHS601*rho;
const double crLHS613 = crLHS596*rho;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS24 + crLHS22*crLHS5 + crLHS25*crLHS34 + crLHS34*crLHS36 - crLHS41*crLHS42 + crLHS45);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS46 + DN(0,1)*crLHS48 + DN(0,2)*crLHS51 + crLHS13*crLHS24 - crLHS42*crLHS61 + crLHS53 - crLHS54*crLHS55 - crLHS54*crLHS58);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS64 + DN(0,2)*crLHS66 + crLHS14*crLHS24 - crLHS42*crLHS71 - crLHS55*crLHS68 - crLHS58*crLHS68 + crLHS67);
rLHS(0,3)+=-gauss_weight*(crLHS25*crLHS73 + crLHS36*crLHS73 + crLHS42*crLHS74 + crLHS72);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS77 + DN(0,2)*crLHS79 + crLHS25*crLHS88 + crLHS36*crLHS88 - crLHS42*crLHS93 + crLHS97 + crLHS99);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS100 + DN(0,1)*crLHS102 + DN(0,2)*crLHS105 + crLHS106 - crLHS109*crLHS42 - crLHS110*crLHS54 + crLHS112);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS113 + DN(0,1)*crLHS115 + DN(0,2)*crLHS117 - crLHS110*crLHS68 + crLHS118 - crLHS121*crLHS42 + crLHS122);
rLHS(0,7)+=-gauss_weight*(crLHS123 + crLHS124*crLHS25 + crLHS124*crLHS36 + crLHS125*crLHS42);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS126 + DN(0,1)*crLHS128 + DN(0,2)*crLHS130 + crLHS139*crLHS25 + crLHS139*crLHS36 - crLHS144*crLHS42 + crLHS148 + crLHS150);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS151 + DN(0,1)*crLHS153 + DN(0,2)*crLHS156 + crLHS157 - crLHS160*crLHS42 - crLHS161*crLHS54 + crLHS163);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS164 + DN(0,1)*crLHS166 + DN(0,2)*crLHS168 - crLHS161*crLHS68 + crLHS169 - crLHS172*crLHS42 + crLHS173);
rLHS(0,11)+=-gauss_weight*(crLHS174 + crLHS175*crLHS25 + crLHS175*crLHS36 + crLHS176*crLHS42);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS177 + DN(0,1)*crLHS179 + DN(0,2)*crLHS181 + crLHS190*crLHS25 + crLHS190*crLHS36 - crLHS195*crLHS42 + crLHS199 + crLHS201);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS202 + DN(0,1)*crLHS204 + DN(0,2)*crLHS207 + crLHS208 - crLHS211*crLHS42 - crLHS212*crLHS54 + crLHS214);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS215 + DN(0,1)*crLHS217 + DN(0,2)*crLHS219 - crLHS212*crLHS68 + crLHS220 - crLHS223*crLHS42 + crLHS224);
rLHS(0,15)+=-gauss_weight*(crLHS225 + crLHS226*crLHS25 + crLHS226*crLHS36 + crLHS227*crLHS42);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS228 + DN(0,2)*crLHS229 + crLHS15*crLHS24 - crLHS230*crLHS55 - crLHS230*crLHS58 - crLHS232*crLHS42 + crLHS53);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS48 + DN(0,1)*crLHS233 + DN(0,2)*crLHS235 + crLHS16*crLHS24 + crLHS22*crLHS236 + crLHS237*crLHS25 + crLHS237*crLHS36 - crLHS239*crLHS42 + crLHS45);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS64 + DN(0,1)*crLHS240 + DN(0,2)*crLHS242 + crLHS17*crLHS24 + crLHS244 - crLHS245*crLHS55 - crLHS245*crLHS58 - crLHS246*crLHS42);
rLHS(1,3)+=-gauss_weight*(crLHS247 + crLHS248*crLHS25 + crLHS248*crLHS36 + crLHS249*crLHS42);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS250 + DN(0,2)*crLHS251 - crLHS110*crLHS230 + crLHS252 - crLHS254*crLHS42 + crLHS256);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS102 + DN(0,1)*crLHS257 + DN(0,2)*crLHS259 + crLHS25*crLHS260 + crLHS260*crLHS36 - crLHS262*crLHS42 + crLHS264 + crLHS97);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS115 + DN(0,1)*crLHS265 + DN(0,2)*crLHS267 - crLHS110*crLHS245 + crLHS268 - crLHS269*crLHS42 + crLHS270);
rLHS(1,7)+=-gauss_weight*(crLHS25*crLHS272 + crLHS271 + crLHS272*crLHS36 + crLHS273*crLHS42);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS128 + DN(0,1)*crLHS274 + DN(0,2)*crLHS275 - crLHS161*crLHS230 + crLHS276 - crLHS278*crLHS42 + crLHS279);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS153 + DN(0,1)*crLHS280 + DN(0,2)*crLHS282 + crLHS148 + crLHS25*crLHS283 + crLHS283*crLHS36 - crLHS285*crLHS42 + crLHS287);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS166 + DN(0,1)*crLHS288 + DN(0,2)*crLHS290 - crLHS161*crLHS245 + crLHS291 - crLHS292*crLHS42 + crLHS293);
rLHS(1,11)+=-gauss_weight*(crLHS25*crLHS295 + crLHS294 + crLHS295*crLHS36 + crLHS296*crLHS42);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS179 + DN(0,1)*crLHS297 + DN(0,2)*crLHS298 - crLHS212*crLHS230 + crLHS299 - crLHS301*crLHS42 + crLHS302);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS204 + DN(0,1)*crLHS303 + DN(0,2)*crLHS305 + crLHS199 + crLHS25*crLHS306 + crLHS306*crLHS36 - crLHS308*crLHS42 + crLHS310);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS217 + DN(0,1)*crLHS311 + DN(0,2)*crLHS313 - crLHS212*crLHS245 + crLHS314 - crLHS315*crLHS42 + crLHS316);
rLHS(1,15)+=-gauss_weight*(crLHS25*crLHS318 + crLHS317 + crLHS318*crLHS36 + crLHS319*crLHS42);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS229 + DN(0,2)*crLHS320 + crLHS18*crLHS24 - crLHS321*crLHS55 - crLHS321*crLHS58 - crLHS322*crLHS42 + crLHS67);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS51 + DN(0,1)*crLHS235 + DN(0,2)*crLHS323 + crLHS19*crLHS24 + crLHS244 - crLHS324*crLHS55 - crLHS324*crLHS58 - crLHS325*crLHS42);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS66 + DN(0,1)*crLHS242 + DN(0,2)*crLHS326 + crLHS20*crLHS24 + crLHS22*crLHS327 + crLHS25*crLHS328 + crLHS328*crLHS36 - crLHS329*crLHS42 + crLHS45);
rLHS(2,3)+=-gauss_weight*(crLHS25*crLHS331 + crLHS330 + crLHS331*crLHS36 + crLHS332*crLHS42);
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS251 + DN(0,2)*crLHS333 - crLHS110*crLHS321 + crLHS335 - crLHS337*crLHS42 + crLHS338);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS105 + DN(0,1)*crLHS259 + DN(0,2)*crLHS339 - crLHS110*crLHS324 + crLHS340 - crLHS341*crLHS42 + crLHS342);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS117 + DN(0,1)*crLHS267 + DN(0,2)*crLHS343 + crLHS25*crLHS344 + crLHS344*crLHS36 - crLHS345*crLHS42 + crLHS347 + crLHS97);
rLHS(2,7)+=-gauss_weight*(crLHS25*crLHS349 + crLHS348 + crLHS349*crLHS36 + crLHS350*crLHS42);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS130 + DN(0,1)*crLHS275 + DN(0,2)*crLHS351 - crLHS161*crLHS321 + crLHS352 - crLHS353*crLHS42 + crLHS354);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS156 + DN(0,1)*crLHS282 + DN(0,2)*crLHS355 - crLHS161*crLHS324 + crLHS356 - crLHS357*crLHS42 + crLHS358);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS168 + DN(0,1)*crLHS290 + DN(0,2)*crLHS359 + crLHS148 + crLHS25*crLHS360 + crLHS36*crLHS360 - crLHS361*crLHS42 + crLHS363);
rLHS(2,11)+=-gauss_weight*(crLHS25*crLHS365 + crLHS36*crLHS365 + crLHS364 + crLHS366*crLHS42);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS181 + DN(0,1)*crLHS298 + DN(0,2)*crLHS367 - crLHS212*crLHS321 + crLHS368 - crLHS369*crLHS42 + crLHS370);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS207 + DN(0,1)*crLHS305 + DN(0,2)*crLHS371 - crLHS212*crLHS324 + crLHS372 - crLHS373*crLHS42 + crLHS374);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS219 + DN(0,1)*crLHS313 + DN(0,2)*crLHS375 + crLHS199 + crLHS25*crLHS376 + crLHS36*crLHS376 - crLHS377*crLHS42 + crLHS379);
rLHS(2,15)+=-gauss_weight*(crLHS25*crLHS381 + crLHS36*crLHS381 + crLHS380 + crLHS382*crLHS42);
rLHS(3,0)+=gauss_weight*(crLHS230*crLHS383 - crLHS32*crLHS73 + crLHS321*crLHS384 + crLHS72);
rLHS(3,1)+=gauss_weight*(crLHS247 - crLHS248*crLHS60 + crLHS324*crLHS384 + crLHS385*crLHS54);
rLHS(3,2)+=gauss_weight*(crLHS245*crLHS383 + crLHS330 - crLHS331*crLHS70 + crLHS385*crLHS68);
rLHS(3,3)+=crLHS386*(crLHS236 + crLHS327 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS230*crLHS388 + crLHS321*crLHS389 + crLHS387 - crLHS73*crLHS87);
rLHS(3,5)+=gauss_weight*(-crLHS108*crLHS248 + crLHS324*crLHS389 + crLHS390 + crLHS391*crLHS54);
rLHS(3,6)+=gauss_weight*(-crLHS120*crLHS331 + crLHS245*crLHS388 + crLHS391*crLHS68 + crLHS392);
rLHS(3,7)+=crLHS393;
rLHS(3,8)+=gauss_weight*(-crLHS138*crLHS73 + crLHS230*crLHS395 + crLHS321*crLHS396 + crLHS394);
rLHS(3,9)+=gauss_weight*(-crLHS159*crLHS248 + crLHS324*crLHS396 + crLHS397 + crLHS398*crLHS54);
rLHS(3,10)+=gauss_weight*(-crLHS171*crLHS331 + crLHS245*crLHS395 + crLHS398*crLHS68 + crLHS399);
rLHS(3,11)+=crLHS400;
rLHS(3,12)+=gauss_weight*(-crLHS189*crLHS73 + crLHS230*crLHS402 + crLHS321*crLHS403 + crLHS401);
rLHS(3,13)+=gauss_weight*(-crLHS210*crLHS248 + crLHS324*crLHS403 + crLHS404 + crLHS405*crLHS54);
rLHS(3,14)+=gauss_weight*(-crLHS222*crLHS331 + crLHS245*crLHS402 + crLHS405*crLHS68 + crLHS406);
rLHS(3,15)+=crLHS407;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS34*crLHS409 + crLHS34*crLHS82 - crLHS41*crLHS410 + crLHS411 + crLHS99);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS46 + DN(1,1)*crLHS48 + DN(1,2)*crLHS51 + crLHS112 + crLHS252 - crLHS410*crLHS61 - crLHS413*crLHS54);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS64 + DN(1,2)*crLHS66 + crLHS122 + crLHS335 - crLHS410*crLHS71 - crLHS413*crLHS68);
rLHS(4,3)+=-gauss_weight*(crLHS387 + crLHS409*crLHS73 + crLHS410*crLHS74 + crLHS73*crLHS82);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS77 + DN(1,2)*crLHS79 + crLHS12*crLHS416 + crLHS22*crLHS414 + crLHS409*crLHS88 - crLHS410*crLHS93 + crLHS418 + crLHS82*crLHS88);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS100 + DN(1,1)*crLHS102 + DN(1,2)*crLHS105 - crLHS109*crLHS410 + crLHS13*crLHS416 + crLHS420 - crLHS421*crLHS54 - crLHS422*crLHS54);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS113 + DN(1,1)*crLHS115 + DN(1,2)*crLHS117 - crLHS121*crLHS410 + crLHS14*crLHS416 - crLHS421*crLHS68 - crLHS422*crLHS68 + crLHS423);
rLHS(4,7)+=-gauss_weight*(crLHS124*crLHS409 + crLHS124*crLHS82 + crLHS125*crLHS410 + crLHS424);
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS126 + DN(1,1)*crLHS128 + DN(1,2)*crLHS130 + crLHS139*crLHS409 + crLHS139*crLHS82 - crLHS144*crLHS410 + crLHS427 + crLHS429);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS151 + DN(1,1)*crLHS153 + DN(1,2)*crLHS156 - crLHS160*crLHS410 + crLHS430 - crLHS431*crLHS54 + crLHS433);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS164 + DN(1,1)*crLHS166 + DN(1,2)*crLHS168 - crLHS172*crLHS410 - crLHS431*crLHS68 + crLHS434 + crLHS435);
rLHS(4,11)+=-gauss_weight*(crLHS175*crLHS409 + crLHS175*crLHS82 + crLHS176*crLHS410 + crLHS436);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS177 + DN(1,1)*crLHS179 + DN(1,2)*crLHS181 + crLHS190*crLHS409 + crLHS190*crLHS82 - crLHS195*crLHS410 + crLHS439 + crLHS441);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS202 + DN(1,1)*crLHS204 + DN(1,2)*crLHS207 - crLHS211*crLHS410 + crLHS442 - crLHS443*crLHS54 + crLHS444);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS215 + DN(1,1)*crLHS217 + DN(1,2)*crLHS219 - crLHS223*crLHS410 - crLHS443*crLHS68 + crLHS445 + crLHS446);
rLHS(4,15)+=-gauss_weight*(crLHS226*crLHS409 + crLHS226*crLHS82 + crLHS227*crLHS410 + crLHS447);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS228 + DN(1,2)*crLHS229 + crLHS106 - crLHS230*crLHS413 - crLHS232*crLHS410 + crLHS256);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS48 + DN(1,1)*crLHS233 + DN(1,2)*crLHS235 + crLHS237*crLHS409 + crLHS237*crLHS82 - crLHS239*crLHS410 + crLHS264 + crLHS411);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS64 + DN(1,1)*crLHS240 + DN(1,2)*crLHS242 - crLHS245*crLHS413 - crLHS246*crLHS410 + crLHS270 + crLHS340);
rLHS(5,3)+=-gauss_weight*(crLHS248*crLHS409 + crLHS248*crLHS82 + crLHS249*crLHS410 + crLHS390);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS250 + DN(1,2)*crLHS251 + crLHS15*crLHS416 - crLHS230*crLHS421 - crLHS230*crLHS422 - crLHS254*crLHS410 + crLHS420);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS102 + DN(1,1)*crLHS257 + DN(1,2)*crLHS259 + crLHS16*crLHS416 + crLHS22*crLHS448 + crLHS260*crLHS409 + crLHS260*crLHS82 - crLHS262*crLHS410 + crLHS418);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS115 + DN(1,1)*crLHS265 + DN(1,2)*crLHS267 + crLHS17*crLHS416 - crLHS245*crLHS421 - crLHS245*crLHS422 - crLHS269*crLHS410 + crLHS450);
rLHS(5,7)+=-gauss_weight*(crLHS272*crLHS409 + crLHS272*crLHS82 + crLHS273*crLHS410 + crLHS451);
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS128 + DN(1,1)*crLHS274 + DN(1,2)*crLHS275 - crLHS230*crLHS431 - crLHS278*crLHS410 + crLHS452 + crLHS456);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS153 + DN(1,1)*crLHS280 + DN(1,2)*crLHS282 + crLHS283*crLHS409 + crLHS283*crLHS82 - crLHS285*crLHS410 + crLHS427 + crLHS458);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS166 + DN(1,1)*crLHS288 + DN(1,2)*crLHS290 - crLHS245*crLHS431 - crLHS292*crLHS410 + crLHS459 + crLHS460);
rLHS(5,11)+=-gauss_weight*(crLHS295*crLHS409 + crLHS295*crLHS82 + crLHS296*crLHS410 + crLHS461);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS179 + DN(1,1)*crLHS297 + DN(1,2)*crLHS298 - crLHS230*crLHS443 - crLHS301*crLHS410 + crLHS462 + crLHS464);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS204 + DN(1,1)*crLHS303 + DN(1,2)*crLHS305 + crLHS306*crLHS409 + crLHS306*crLHS82 - crLHS308*crLHS410 + crLHS439 + crLHS466);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS217 + DN(1,1)*crLHS311 + DN(1,2)*crLHS313 - crLHS245*crLHS443 - crLHS315*crLHS410 + crLHS467 + crLHS468);
rLHS(5,15)+=-gauss_weight*(crLHS318*crLHS409 + crLHS318*crLHS82 + crLHS319*crLHS410 + crLHS469);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS229 + DN(1,2)*crLHS320 + crLHS118 - crLHS321*crLHS413 - crLHS322*crLHS410 + crLHS338);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS51 + DN(1,1)*crLHS235 + DN(1,2)*crLHS323 + crLHS268 - crLHS324*crLHS413 - crLHS325*crLHS410 + crLHS342);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS66 + DN(1,1)*crLHS242 + DN(1,2)*crLHS326 + crLHS328*crLHS409 + crLHS328*crLHS82 - crLHS329*crLHS410 + crLHS347 + crLHS411);
rLHS(6,3)+=-gauss_weight*(crLHS331*crLHS409 + crLHS331*crLHS82 + crLHS332*crLHS410 + crLHS392);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS251 + DN(1,2)*crLHS333 + crLHS18*crLHS416 - crLHS321*crLHS421 - crLHS321*crLHS422 - crLHS337*crLHS410 + crLHS423);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS105 + DN(1,1)*crLHS259 + DN(1,2)*crLHS339 + crLHS19*crLHS416 - crLHS324*crLHS421 - crLHS324*crLHS422 - crLHS341*crLHS410 + crLHS450);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS117 + DN(1,1)*crLHS267 + DN(1,2)*crLHS343 + crLHS20*crLHS416 + crLHS22*crLHS470 + crLHS344*crLHS409 + crLHS344*crLHS82 - crLHS345*crLHS410 + crLHS418);
rLHS(6,7)+=-gauss_weight*(crLHS349*crLHS409 + crLHS349*crLHS82 + crLHS350*crLHS410 + crLHS471);
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS130 + DN(1,1)*crLHS275 + DN(1,2)*crLHS351 - crLHS321*crLHS431 - crLHS353*crLHS410 + crLHS473 + crLHS475);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS156 + DN(1,1)*crLHS282 + DN(1,2)*crLHS355 - crLHS324*crLHS431 - crLHS357*crLHS410 + crLHS476 + crLHS477);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS168 + DN(1,1)*crLHS290 + DN(1,2)*crLHS359 + crLHS360*crLHS409 + crLHS360*crLHS82 - crLHS361*crLHS410 + crLHS427 + crLHS479);
rLHS(6,11)+=-gauss_weight*(crLHS365*crLHS409 + crLHS365*crLHS82 + crLHS366*crLHS410 + crLHS480);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS181 + DN(1,1)*crLHS298 + DN(1,2)*crLHS367 - crLHS321*crLHS443 - crLHS369*crLHS410 + crLHS481 + crLHS483);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS207 + DN(1,1)*crLHS305 + DN(1,2)*crLHS371 - crLHS324*crLHS443 - crLHS373*crLHS410 + crLHS484 + crLHS485);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS219 + DN(1,1)*crLHS313 + DN(1,2)*crLHS375 + crLHS376*crLHS409 + crLHS376*crLHS82 - crLHS377*crLHS410 + crLHS439 + crLHS487);
rLHS(6,15)+=-gauss_weight*(crLHS381*crLHS409 + crLHS381*crLHS82 + crLHS382*crLHS410 + crLHS488);
rLHS(7,0)+=gauss_weight*(crLHS123 - crLHS124*crLHS32 + crLHS230*crLHS489 + crLHS321*crLHS490);
rLHS(7,1)+=gauss_weight*(crLHS271 - crLHS272*crLHS60 + crLHS324*crLHS490 + crLHS491*crLHS54);
rLHS(7,2)+=gauss_weight*(crLHS245*crLHS489 + crLHS348 - crLHS349*crLHS70 + crLHS491*crLHS68);
rLHS(7,3)+=crLHS393;
rLHS(7,4)+=gauss_weight*(-crLHS124*crLHS87 + crLHS230*crLHS492 + crLHS321*crLHS493 + crLHS424);
rLHS(7,5)+=gauss_weight*(-crLHS108*crLHS272 + crLHS324*crLHS493 + crLHS451 + crLHS494*crLHS54);
rLHS(7,6)+=gauss_weight*(-crLHS120*crLHS349 + crLHS245*crLHS492 + crLHS471 + crLHS494*crLHS68);
rLHS(7,7)+=crLHS386*(crLHS414 + crLHS448 + crLHS470);
rLHS(7,8)+=gauss_weight*(-crLHS124*crLHS138 + crLHS230*crLHS496 + crLHS321*crLHS497 + crLHS495);
rLHS(7,9)+=gauss_weight*(-crLHS159*crLHS272 + crLHS324*crLHS497 + crLHS498 + crLHS499*crLHS54);
rLHS(7,10)+=gauss_weight*(-crLHS171*crLHS349 + crLHS245*crLHS496 + crLHS499*crLHS68 + crLHS500);
rLHS(7,11)+=crLHS501;
rLHS(7,12)+=gauss_weight*(-crLHS124*crLHS189 + crLHS230*crLHS503 + crLHS321*crLHS504 + crLHS502);
rLHS(7,13)+=gauss_weight*(-crLHS210*crLHS272 + crLHS324*crLHS504 + crLHS505 + crLHS506*crLHS54);
rLHS(7,14)+=gauss_weight*(-crLHS222*crLHS349 + crLHS245*crLHS503 + crLHS506*crLHS68 + crLHS507);
rLHS(7,15)+=crLHS508;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS133*crLHS34 + crLHS150 + crLHS34*crLHS510 - crLHS41*crLHS454 + crLHS511);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS46 + DN(2,1)*crLHS48 + DN(2,2)*crLHS51 + crLHS163 + crLHS276 - crLHS454*crLHS61 - crLHS513*crLHS54);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS64 + DN(2,2)*crLHS66 + crLHS173 + crLHS352 - crLHS454*crLHS71 - crLHS513*crLHS68);
rLHS(8,3)+=-gauss_weight*(crLHS133*crLHS73 + crLHS394 + crLHS454*crLHS74 + crLHS510*crLHS73);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS77 + DN(2,2)*crLHS79 + crLHS133*crLHS88 + crLHS429 - crLHS454*crLHS93 + crLHS510*crLHS88 + crLHS514);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS100 + DN(2,1)*crLHS102 + DN(2,2)*crLHS105 - crLHS109*crLHS454 + crLHS433 + crLHS452 - crLHS515*crLHS54);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS113 + DN(2,1)*crLHS115 + DN(2,2)*crLHS117 - crLHS121*crLHS454 + crLHS435 + crLHS473 - crLHS515*crLHS68);
rLHS(8,7)+=-gauss_weight*(crLHS124*crLHS133 + crLHS124*crLHS510 + crLHS125*crLHS454 + crLHS495);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS126 + DN(2,1)*crLHS128 + DN(2,2)*crLHS130 + crLHS12*crLHS518 + crLHS133*crLHS139 + crLHS139*crLHS510 - crLHS144*crLHS454 + crLHS22*crLHS516 + crLHS520);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS151 + DN(2,1)*crLHS153 + DN(2,2)*crLHS156 + crLHS13*crLHS518 - crLHS160*crLHS454 + crLHS522 - crLHS523*crLHS54 - crLHS524*crLHS54);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS164 + DN(2,1)*crLHS166 + DN(2,2)*crLHS168 + crLHS14*crLHS518 - crLHS172*crLHS454 - crLHS523*crLHS68 - crLHS524*crLHS68 + crLHS525);
rLHS(8,11)+=-gauss_weight*(crLHS133*crLHS175 + crLHS175*crLHS510 + crLHS176*crLHS454 + crLHS526);
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS177 + DN(2,1)*crLHS179 + DN(2,2)*crLHS181 + crLHS133*crLHS190 + crLHS190*crLHS510 - crLHS195*crLHS454 + crLHS529 + crLHS531);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS202 + DN(2,1)*crLHS204 + DN(2,2)*crLHS207 - crLHS211*crLHS454 + crLHS532 - crLHS533*crLHS54 + crLHS535);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS215 + DN(2,1)*crLHS217 + DN(2,2)*crLHS219 - crLHS223*crLHS454 - crLHS533*crLHS68 + crLHS536 + crLHS537);
rLHS(8,15)+=-gauss_weight*(crLHS133*crLHS226 + crLHS226*crLHS510 + crLHS227*crLHS454 + crLHS538);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS228 + DN(2,2)*crLHS229 + crLHS157 - crLHS230*crLHS513 - crLHS232*crLHS454 + crLHS279);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS48 + DN(2,1)*crLHS233 + DN(2,2)*crLHS235 + crLHS133*crLHS237 + crLHS237*crLHS510 - crLHS239*crLHS454 + crLHS287 + crLHS511);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS64 + DN(2,1)*crLHS240 + DN(2,2)*crLHS242 - crLHS245*crLHS513 - crLHS246*crLHS454 + crLHS293 + crLHS356);
rLHS(9,3)+=-gauss_weight*(crLHS133*crLHS248 + crLHS248*crLHS510 + crLHS249*crLHS454 + crLHS397);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS250 + DN(2,2)*crLHS251 - crLHS230*crLHS515 - crLHS254*crLHS454 + crLHS430 + crLHS456);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS102 + DN(2,1)*crLHS257 + DN(2,2)*crLHS259 + crLHS133*crLHS260 + crLHS260*crLHS510 - crLHS262*crLHS454 + crLHS458 + crLHS514);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS115 + DN(2,1)*crLHS265 + DN(2,2)*crLHS267 - crLHS245*crLHS515 - crLHS269*crLHS454 + crLHS460 + crLHS476);
rLHS(9,7)+=-gauss_weight*(crLHS133*crLHS272 + crLHS272*crLHS510 + crLHS273*crLHS454 + crLHS498);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS128 + DN(2,1)*crLHS274 + DN(2,2)*crLHS275 + crLHS15*crLHS518 - crLHS230*crLHS523 - crLHS230*crLHS524 - crLHS278*crLHS454 + crLHS522);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS153 + DN(2,1)*crLHS280 + DN(2,2)*crLHS282 + crLHS133*crLHS283 + crLHS16*crLHS518 + crLHS22*crLHS539 + crLHS283*crLHS510 - crLHS285*crLHS454 + crLHS520);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS166 + DN(2,1)*crLHS288 + DN(2,2)*crLHS290 + crLHS17*crLHS518 - crLHS245*crLHS523 - crLHS245*crLHS524 - crLHS292*crLHS454 + crLHS541);
rLHS(9,11)+=-gauss_weight*(crLHS133*crLHS295 + crLHS295*crLHS510 + crLHS296*crLHS454 + crLHS542);
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS179 + DN(2,1)*crLHS297 + DN(2,2)*crLHS298 - crLHS230*crLHS533 - crLHS301*crLHS454 + crLHS543 + crLHS545);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS204 + DN(2,1)*crLHS303 + DN(2,2)*crLHS305 + crLHS133*crLHS306 + crLHS306*crLHS510 - crLHS308*crLHS454 + crLHS529 + crLHS547);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS217 + DN(2,1)*crLHS311 + DN(2,2)*crLHS313 - crLHS245*crLHS533 - crLHS315*crLHS454 + crLHS548 + crLHS549);
rLHS(9,15)+=-gauss_weight*(crLHS133*crLHS318 + crLHS318*crLHS510 + crLHS319*crLHS454 + crLHS550);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS229 + DN(2,2)*crLHS320 + crLHS169 - crLHS321*crLHS513 - crLHS322*crLHS454 + crLHS354);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS51 + DN(2,1)*crLHS235 + DN(2,2)*crLHS323 + crLHS291 - crLHS324*crLHS513 - crLHS325*crLHS454 + crLHS358);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS66 + DN(2,1)*crLHS242 + DN(2,2)*crLHS326 + crLHS133*crLHS328 + crLHS328*crLHS510 - crLHS329*crLHS454 + crLHS363 + crLHS511);
rLHS(10,3)+=-gauss_weight*(crLHS133*crLHS331 + crLHS331*crLHS510 + crLHS332*crLHS454 + crLHS399);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS251 + DN(2,2)*crLHS333 - crLHS321*crLHS515 - crLHS337*crLHS454 + crLHS434 + crLHS475);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS105 + DN(2,1)*crLHS259 + DN(2,2)*crLHS339 - crLHS324*crLHS515 - crLHS341*crLHS454 + crLHS459 + crLHS477);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS117 + DN(2,1)*crLHS267 + DN(2,2)*crLHS343 + crLHS133*crLHS344 + crLHS344*crLHS510 - crLHS345*crLHS454 + crLHS479 + crLHS514);
rLHS(10,7)+=-gauss_weight*(crLHS133*crLHS349 + crLHS349*crLHS510 + crLHS350*crLHS454 + crLHS500);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS130 + DN(2,1)*crLHS275 + DN(2,2)*crLHS351 + crLHS18*crLHS518 - crLHS321*crLHS523 - crLHS321*crLHS524 - crLHS353*crLHS454 + crLHS525);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS156 + DN(2,1)*crLHS282 + DN(2,2)*crLHS355 + crLHS19*crLHS518 - crLHS324*crLHS523 - crLHS324*crLHS524 - crLHS357*crLHS454 + crLHS541);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS168 + DN(2,1)*crLHS290 + DN(2,2)*crLHS359 + crLHS133*crLHS360 + crLHS20*crLHS518 + crLHS22*crLHS551 + crLHS360*crLHS510 - crLHS361*crLHS454 + crLHS520);
rLHS(10,11)+=-gauss_weight*(crLHS133*crLHS365 + crLHS365*crLHS510 + crLHS366*crLHS454 + crLHS552);
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS181 + DN(2,1)*crLHS298 + DN(2,2)*crLHS367 - crLHS321*crLHS533 - crLHS369*crLHS454 + crLHS554 + crLHS555);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS207 + DN(2,1)*crLHS305 + DN(2,2)*crLHS371 - crLHS324*crLHS533 - crLHS373*crLHS454 + crLHS556 + crLHS557);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS219 + DN(2,1)*crLHS313 + DN(2,2)*crLHS375 + crLHS133*crLHS376 + crLHS376*crLHS510 - crLHS377*crLHS454 + crLHS529 + crLHS559);
rLHS(10,15)+=-gauss_weight*(crLHS133*crLHS381 + crLHS381*crLHS510 + crLHS382*crLHS454 + crLHS560);
rLHS(11,0)+=gauss_weight*(crLHS174 - crLHS175*crLHS32 + crLHS230*crLHS561 + crLHS321*crLHS562);
rLHS(11,1)+=gauss_weight*(crLHS294 - crLHS295*crLHS60 + crLHS324*crLHS562 + crLHS54*crLHS563);
rLHS(11,2)+=gauss_weight*(crLHS245*crLHS561 + crLHS364 - crLHS365*crLHS70 + crLHS563*crLHS68);
rLHS(11,3)+=crLHS400;
rLHS(11,4)+=gauss_weight*(-crLHS175*crLHS87 + crLHS230*crLHS564 + crLHS321*crLHS565 + crLHS436);
rLHS(11,5)+=gauss_weight*(-crLHS108*crLHS295 + crLHS324*crLHS565 + crLHS461 + crLHS54*crLHS566);
rLHS(11,6)+=gauss_weight*(-crLHS120*crLHS365 + crLHS245*crLHS564 + crLHS480 + crLHS566*crLHS68);
rLHS(11,7)+=crLHS501;
rLHS(11,8)+=gauss_weight*(-crLHS138*crLHS175 + crLHS230*crLHS567 + crLHS321*crLHS568 + crLHS526);
rLHS(11,9)+=gauss_weight*(-crLHS159*crLHS295 + crLHS324*crLHS568 + crLHS54*crLHS569 + crLHS542);
rLHS(11,10)+=gauss_weight*(-crLHS171*crLHS365 + crLHS245*crLHS567 + crLHS552 + crLHS569*crLHS68);
rLHS(11,11)+=crLHS386*(crLHS516 + crLHS539 + crLHS551);
rLHS(11,12)+=gauss_weight*(-crLHS175*crLHS189 + crLHS230*crLHS571 + crLHS321*crLHS572 + crLHS570);
rLHS(11,13)+=gauss_weight*(-crLHS210*crLHS295 + crLHS324*crLHS572 + crLHS54*crLHS574 + crLHS573);
rLHS(11,14)+=gauss_weight*(-crLHS222*crLHS365 + crLHS245*crLHS571 + crLHS574*crLHS68 + crLHS575);
rLHS(11,15)+=crLHS576;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS184*crLHS34 + crLHS201 + crLHS34*crLHS578 - crLHS41*crLHS463 + crLHS579);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS46 + DN(3,1)*crLHS48 + DN(3,2)*crLHS51 + crLHS214 + crLHS299 - crLHS463*crLHS61 - crLHS54*crLHS581);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS62 + DN(3,1)*crLHS64 + DN(3,2)*crLHS66 + crLHS224 + crLHS368 - crLHS463*crLHS71 - crLHS581*crLHS68);
rLHS(12,3)+=-gauss_weight*(crLHS184*crLHS73 + crLHS401 + crLHS463*crLHS74 + crLHS578*crLHS73);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS75 + DN(3,1)*crLHS77 + DN(3,2)*crLHS79 + crLHS184*crLHS88 + crLHS441 - crLHS463*crLHS93 + crLHS578*crLHS88 + crLHS582);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS100 + DN(3,1)*crLHS102 + DN(3,2)*crLHS105 - crLHS109*crLHS463 + crLHS444 + crLHS462 - crLHS54*crLHS583);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS113 + DN(3,1)*crLHS115 + DN(3,2)*crLHS117 - crLHS121*crLHS463 + crLHS446 + crLHS481 - crLHS583*crLHS68);
rLHS(12,7)+=-gauss_weight*(crLHS124*crLHS184 + crLHS124*crLHS578 + crLHS125*crLHS463 + crLHS502);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS126 + DN(3,1)*crLHS128 + DN(3,2)*crLHS130 + crLHS139*crLHS184 + crLHS139*crLHS578 - crLHS144*crLHS463 + crLHS531 + crLHS584);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS151 + DN(3,1)*crLHS153 + DN(3,2)*crLHS156 - crLHS160*crLHS463 + crLHS535 - crLHS54*crLHS585 + crLHS543);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS164 + DN(3,1)*crLHS166 + DN(3,2)*crLHS168 - crLHS172*crLHS463 + crLHS537 + crLHS554 - crLHS585*crLHS68);
rLHS(12,11)+=-gauss_weight*(crLHS175*crLHS184 + crLHS175*crLHS578 + crLHS176*crLHS463 + crLHS570);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS177 + DN(3,1)*crLHS179 + DN(3,2)*crLHS181 + crLHS12*crLHS588 + crLHS184*crLHS190 + crLHS190*crLHS578 - crLHS195*crLHS463 + crLHS22*crLHS586 + crLHS590);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS202 + DN(3,1)*crLHS204 + DN(3,2)*crLHS207 + crLHS13*crLHS588 - crLHS211*crLHS463 - crLHS54*crLHS593 - crLHS54*crLHS594 + crLHS592);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS215 + DN(3,1)*crLHS217 + DN(3,2)*crLHS219 + crLHS14*crLHS588 - crLHS223*crLHS463 - crLHS593*crLHS68 - crLHS594*crLHS68 + crLHS595);
rLHS(12,15)+=-gauss_weight*(crLHS184*crLHS226 + crLHS226*crLHS578 + crLHS227*crLHS463 + crLHS596);
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS228 + DN(3,2)*crLHS229 + crLHS208 - crLHS230*crLHS581 - crLHS232*crLHS463 + crLHS302);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS48 + DN(3,1)*crLHS233 + DN(3,2)*crLHS235 + crLHS184*crLHS237 + crLHS237*crLHS578 - crLHS239*crLHS463 + crLHS310 + crLHS579);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS64 + DN(3,1)*crLHS240 + DN(3,2)*crLHS242 - crLHS245*crLHS581 - crLHS246*crLHS463 + crLHS316 + crLHS372);
rLHS(13,3)+=-gauss_weight*(crLHS184*crLHS248 + crLHS248*crLHS578 + crLHS249*crLHS463 + crLHS404);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS250 + DN(3,2)*crLHS251 - crLHS230*crLHS583 - crLHS254*crLHS463 + crLHS442 + crLHS464);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS102 + DN(3,1)*crLHS257 + DN(3,2)*crLHS259 + crLHS184*crLHS260 + crLHS260*crLHS578 - crLHS262*crLHS463 + crLHS466 + crLHS582);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS115 + DN(3,1)*crLHS265 + DN(3,2)*crLHS267 - crLHS245*crLHS583 - crLHS269*crLHS463 + crLHS468 + crLHS484);
rLHS(13,7)+=-gauss_weight*(crLHS184*crLHS272 + crLHS272*crLHS578 + crLHS273*crLHS463 + crLHS505);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS128 + DN(3,1)*crLHS274 + DN(3,2)*crLHS275 - crLHS230*crLHS585 - crLHS278*crLHS463 + crLHS532 + crLHS545);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS153 + DN(3,1)*crLHS280 + DN(3,2)*crLHS282 + crLHS184*crLHS283 + crLHS283*crLHS578 - crLHS285*crLHS463 + crLHS547 + crLHS584);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS166 + DN(3,1)*crLHS288 + DN(3,2)*crLHS290 - crLHS245*crLHS585 - crLHS292*crLHS463 + crLHS549 + crLHS556);
rLHS(13,11)+=-gauss_weight*(crLHS184*crLHS295 + crLHS295*crLHS578 + crLHS296*crLHS463 + crLHS573);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS179 + DN(3,1)*crLHS297 + DN(3,2)*crLHS298 + crLHS15*crLHS588 - crLHS230*crLHS593 - crLHS230*crLHS594 - crLHS301*crLHS463 + crLHS592);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS204 + DN(3,1)*crLHS303 + DN(3,2)*crLHS305 + crLHS16*crLHS588 + crLHS184*crLHS306 + crLHS22*crLHS597 + crLHS306*crLHS578 - crLHS308*crLHS463 + crLHS590);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS217 + DN(3,1)*crLHS311 + DN(3,2)*crLHS313 + crLHS17*crLHS588 - crLHS245*crLHS593 - crLHS245*crLHS594 - crLHS315*crLHS463 + crLHS598);
rLHS(13,15)+=-gauss_weight*(crLHS184*crLHS318 + crLHS318*crLHS578 + crLHS319*crLHS463 + crLHS599);
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS229 + DN(3,2)*crLHS320 + crLHS220 - crLHS321*crLHS581 - crLHS322*crLHS463 + crLHS370);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS51 + DN(3,1)*crLHS235 + DN(3,2)*crLHS323 + crLHS314 - crLHS324*crLHS581 - crLHS325*crLHS463 + crLHS374);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS66 + DN(3,1)*crLHS242 + DN(3,2)*crLHS326 + crLHS184*crLHS328 + crLHS328*crLHS578 - crLHS329*crLHS463 + crLHS379 + crLHS579);
rLHS(14,3)+=-gauss_weight*(crLHS184*crLHS331 + crLHS331*crLHS578 + crLHS332*crLHS463 + crLHS406);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS79 + DN(3,1)*crLHS251 + DN(3,2)*crLHS333 - crLHS321*crLHS583 - crLHS337*crLHS463 + crLHS445 + crLHS483);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS105 + DN(3,1)*crLHS259 + DN(3,2)*crLHS339 - crLHS324*crLHS583 - crLHS341*crLHS463 + crLHS467 + crLHS485);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS117 + DN(3,1)*crLHS267 + DN(3,2)*crLHS343 + crLHS184*crLHS344 + crLHS344*crLHS578 - crLHS345*crLHS463 + crLHS487 + crLHS582);
rLHS(14,7)+=-gauss_weight*(crLHS184*crLHS349 + crLHS349*crLHS578 + crLHS350*crLHS463 + crLHS507);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS130 + DN(3,1)*crLHS275 + DN(3,2)*crLHS351 - crLHS321*crLHS585 - crLHS353*crLHS463 + crLHS536 + crLHS555);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS156 + DN(3,1)*crLHS282 + DN(3,2)*crLHS355 - crLHS324*crLHS585 - crLHS357*crLHS463 + crLHS548 + crLHS557);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS168 + DN(3,1)*crLHS290 + DN(3,2)*crLHS359 + crLHS184*crLHS360 + crLHS360*crLHS578 - crLHS361*crLHS463 + crLHS559 + crLHS584);
rLHS(14,11)+=-gauss_weight*(crLHS184*crLHS365 + crLHS365*crLHS578 + crLHS366*crLHS463 + crLHS575);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS181 + DN(3,1)*crLHS298 + DN(3,2)*crLHS367 + crLHS18*crLHS588 - crLHS321*crLHS593 - crLHS321*crLHS594 - crLHS369*crLHS463 + crLHS595);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS207 + DN(3,1)*crLHS305 + DN(3,2)*crLHS371 + crLHS19*crLHS588 - crLHS324*crLHS593 - crLHS324*crLHS594 - crLHS373*crLHS463 + crLHS598);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS219 + DN(3,1)*crLHS313 + DN(3,2)*crLHS375 + crLHS184*crLHS376 + crLHS20*crLHS588 + crLHS22*crLHS600 + crLHS376*crLHS578 - crLHS377*crLHS463 + crLHS590);
rLHS(14,15)+=-gauss_weight*(crLHS184*crLHS381 + crLHS381*crLHS578 + crLHS382*crLHS463 + crLHS601);
rLHS(15,0)+=gauss_weight*(crLHS225 - crLHS226*crLHS32 + crLHS230*crLHS602 + crLHS321*crLHS603);
rLHS(15,1)+=gauss_weight*(crLHS317 - crLHS318*crLHS60 + crLHS324*crLHS603 + crLHS54*crLHS604);
rLHS(15,2)+=gauss_weight*(crLHS245*crLHS602 + crLHS380 - crLHS381*crLHS70 + crLHS604*crLHS68);
rLHS(15,3)+=crLHS407;
rLHS(15,4)+=gauss_weight*(-crLHS226*crLHS87 + crLHS230*crLHS605 + crLHS321*crLHS606 + crLHS447);
rLHS(15,5)+=gauss_weight*(-crLHS108*crLHS318 + crLHS324*crLHS606 + crLHS469 + crLHS54*crLHS607);
rLHS(15,6)+=gauss_weight*(-crLHS120*crLHS381 + crLHS245*crLHS605 + crLHS488 + crLHS607*crLHS68);
rLHS(15,7)+=crLHS508;
rLHS(15,8)+=gauss_weight*(-crLHS138*crLHS226 + crLHS230*crLHS608 + crLHS321*crLHS609 + crLHS538);
rLHS(15,9)+=gauss_weight*(-crLHS159*crLHS318 + crLHS324*crLHS609 + crLHS54*crLHS610 + crLHS550);
rLHS(15,10)+=gauss_weight*(-crLHS171*crLHS381 + crLHS245*crLHS608 + crLHS560 + crLHS610*crLHS68);
rLHS(15,11)+=crLHS576;
rLHS(15,12)+=gauss_weight*(-crLHS189*crLHS226 + crLHS230*crLHS611 + crLHS321*crLHS612 + crLHS596);
rLHS(15,13)+=gauss_weight*(-crLHS210*crLHS318 + crLHS324*crLHS612 + crLHS54*crLHS613 + crLHS599);
rLHS(15,14)+=gauss_weight*(-crLHS222*crLHS381 + crLHS245*crLHS611 + crLHS601 + crLHS613*crLHS68);
rLHS(15,15)+=crLHS386*(crLHS586 + crLHS597 + crLHS600);

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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
    const BoundedMatrix<double,2,3> v_conv_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
    const BoundedMatrix<double,2,3>& f_adj = rData.BodyForce_adj;
    const array_1d<double,3>& p_adj = rData.Pressure_adj;
    const BoundedMatrix<double,2,3> functional_v = rData.Functional_derivative_velocity;
    const array_1d<double,3>& functional_t = rData.Functional_derivative_transport_scalar;
    const array_1d<double,3>& functional_t_adj = rData.Functional_derivative_transport_scalar_adj;
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
const double crRHS5 = N[0]*functional_v(0,0) + N[1]*functional_v(1,0) + N[2]*functional_v(2,0);
const double crRHS6 = 2.0*crRHS2*functional_weights[0];
const double crRHS7 = N[0]*crRHS6;
const double crRHS8 = DN(0,1)*functional_v(0,0);
const double crRHS9 = DN(1,1)*functional_v(1,0);
const double crRHS10 = DN(2,1)*functional_v(2,0);
const double crRHS11 = DN(0,0)*functional_v(0,1) + DN(1,0)*functional_v(1,1) + DN(2,0)*functional_v(2,1);
const double crRHS12 = 2.0*functional_weights[2]*mu*(-crRHS10 + crRHS11 - crRHS8 - crRHS9);
const double crRHS13 = DN(0,0)*functional_t[0] + DN(1,0)*functional_t[1] + DN(2,0)*functional_t[2];
const double crRHS14 = N[0]*t_ConvCoeff[0] + N[1]*t_ConvCoeff[1] + N[2]*t_ConvCoeff[2];
const double crRHS15 = crRHS14*(N[0]*functional_t_adj[0] + N[1]*functional_t_adj[1] + N[2]*functional_t_adj[2]);
const double crRHS16 = crRHS13*crRHS15;
const double crRHS17 = crRHS14*functional_weights[6]*(N[0]*functional_t[0] + N[1]*functional_t[1] + N[2]*functional_t[2]);
const double crRHS18 = crRHS13*crRHS17;
const double crRHS19 = N[0]*v_conv_ns(0,0) + N[1]*v_conv_ns(1,0) + N[2]*v_conv_ns(2,0);
const double crRHS20 = N[0]*v_conv_ns(0,1) + N[1]*v_conv_ns(1,1) + N[2]*v_conv_ns(2,1);
const double crRHS21 = rho*(DN(0,0)*crRHS19 + DN(0,1)*crRHS20);
const double crRHS22 = 1.0*DN(0,0)*functional_v(0,0) + 1.0*DN(1,0)*functional_v(1,0) + 1.0*DN(2,0)*functional_v(2,0);
const double crRHS23 = 0.5*crRHS10 + 0.5*crRHS11 + 0.5*crRHS8 + 0.5*crRHS9;
const double crRHS24 = 4.0*functional_weights[1]*mu;
const double crRHS25 = N[0]*(bdf0*v_adj(0,0) + bdf1*vn_adj(0,0) + bdf2*vnn_adj(0,0)) + N[1]*(bdf0*v_adj(1,0) + bdf1*vn_adj(1,0) + bdf2*vnn_adj(1,0)) + N[2]*(bdf0*v_adj(2,0) + bdf1*vn_adj(2,0) + bdf2*vnn_adj(2,0));
const double crRHS26 = N[0]*rho;
const double crRHS27 = crRHS26*time_coeff;
const double crRHS28 = DN(0,0)*v_conv_ns(0,0) + DN(1,0)*v_conv_ns(1,0) + DN(2,0)*v_conv_ns(2,0);
const double crRHS29 = crRHS28*crRHS3;
const double crRHS30 = DN(0,0)*v_conv_ns(0,1) + DN(1,0)*v_conv_ns(1,1) + DN(2,0)*v_conv_ns(2,1);
const double crRHS31 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1);
const double crRHS32 = crRHS30*crRHS31;
const double crRHS33 = crRHS29 + crRHS32;
const double crRHS34 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0);
const double crRHS35 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1);
const double crRHS36 = crRHS34 + crRHS35;
const double crRHS37 = crRHS2*stab_c3;
const double crRHS38 = rho*stab_c2*sqrt(pow(crRHS19, 2) + pow(crRHS20, 2));
const double crRHS39 = DN(0,1)*v_conv_ns(0,0) + DN(1,1)*v_conv_ns(1,0) + DN(2,1)*v_conv_ns(2,0);
const double crRHS40 = DN(0,1)*v_conv_ns(0,1) + DN(1,1)*v_conv_ns(1,1) + DN(2,1)*v_conv_ns(2,1);
const double crRHS41 = rho*stab_c3*sqrt(pow(crRHS28, 2) + pow(crRHS30, 2) + pow(crRHS39, 2) + pow(crRHS40, 2));
const double crRHS42 = crRHS36*(h*(crRHS37*h + crRHS38 + crRHS41*h)/stab_c1 + mu);
const double crRHS43 = crRHS19*rho;
const double crRHS44 = crRHS20*rho;
const double crRHS45 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS16 + crRHS18 + crRHS19*crRHS6 + crRHS29*rho + crRHS32*rho - crRHS34*crRHS43 + crRHS4 - crRHS44*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0));
const double crRHS46 = 1.0/(crRHS37 + crRHS38/h + crRHS41 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS47 = crRHS45*crRHS46;
const double crRHS48 = N[0]*crRHS2;
const double crRHS49 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS50 = crRHS2*crRHS31;
const double crRHS51 = crRHS3*crRHS39;
const double crRHS52 = crRHS31*crRHS40;
const double crRHS53 = DN(0,1)*functional_t[0] + DN(1,1)*functional_t[1] + DN(2,1)*functional_t[2];
const double crRHS54 = crRHS15*crRHS53;
const double crRHS55 = crRHS17*crRHS53;
const double crRHS56 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + crRHS20*crRHS6 - crRHS35*crRHS44 - crRHS43*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS49 + crRHS50 + crRHS51*rho + crRHS52*rho + crRHS54 + crRHS55;
const double crRHS57 = crRHS28*crRHS45 + crRHS30*crRHS56;
const double crRHS58 = crRHS26*crRHS46;
const double crRHS59 = N[0]*functional_v(0,1) + N[1]*functional_v(1,1) + N[2]*functional_v(2,1);
const double crRHS60 = 1.0*DN(0,1)*functional_v(0,1) + 1.0*DN(1,1)*functional_v(1,1) + 1.0*DN(2,1)*functional_v(2,1);
const double crRHS61 = N[0]*(bdf0*v_adj(0,1) + bdf1*vn_adj(0,1) + bdf2*vnn_adj(0,1)) + N[1]*(bdf0*v_adj(1,1) + bdf1*vn_adj(1,1) + bdf2*vnn_adj(1,1)) + N[2]*(bdf0*v_adj(2,1) + bdf1*vn_adj(2,1) + bdf2*vnn_adj(2,1));
const double crRHS62 = crRHS51 + crRHS52;
const double crRHS63 = crRHS46*crRHS56;
const double crRHS64 = crRHS39*crRHS45 + crRHS40*crRHS56;
const double crRHS65 = N[1]*crRHS6;
const double crRHS66 = rho*(DN(1,0)*crRHS19 + DN(1,1)*crRHS20);
const double crRHS67 = N[1]*rho;
const double crRHS68 = crRHS67*time_coeff;
const double crRHS69 = N[1]*crRHS2;
const double crRHS70 = crRHS46*crRHS67;
const double crRHS71 = N[2]*crRHS6;
const double crRHS72 = rho*(DN(2,0)*crRHS19 + DN(2,1)*crRHS20);
const double crRHS73 = N[2]*rho;
const double crRHS74 = crRHS73*time_coeff;
const double crRHS75 = N[2]*crRHS2;
const double crRHS76 = crRHS46*crRHS73;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS42 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS12 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS16 + N[0]*crRHS18 + N[0]*crRHS4 + crRHS21*crRHS3 - crRHS21*crRHS47 + crRHS24*(DN(0,0)*crRHS22 + DN(0,1)*crRHS23) + crRHS25*crRHS27 + crRHS26*crRHS33 - crRHS47*crRHS48 + crRHS5*crRHS7 - crRHS57*crRHS58);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS12 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS42 + DN(0,1)*stress_adj[1] - N[0]*crRHS49 + N[0]*crRHS50 + N[0]*crRHS54 + N[0]*crRHS55 + crRHS21*crRHS31 - crRHS21*crRHS63 + crRHS24*(DN(0,0)*crRHS23 + DN(0,1)*crRHS60) + crRHS26*crRHS62 + crRHS27*crRHS61 - crRHS48*crRHS63 - crRHS58*crRHS64 + crRHS59*crRHS7);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS47 + DN(0,1)*crRHS63 + N[0]*crRHS36);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS42 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS12 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS16 + N[1]*crRHS18 + N[1]*crRHS4 + crRHS24*(DN(1,0)*crRHS22 + DN(1,1)*crRHS23) + crRHS25*crRHS68 + crRHS3*crRHS66 + crRHS33*crRHS67 - crRHS47*crRHS66 - crRHS47*crRHS69 + crRHS5*crRHS65 - crRHS57*crRHS70);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS12 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS42 + DN(1,1)*stress_adj[1] - N[1]*crRHS49 + N[1]*crRHS50 + N[1]*crRHS54 + N[1]*crRHS55 + crRHS24*(DN(1,0)*crRHS23 + DN(1,1)*crRHS60) + crRHS31*crRHS66 + crRHS59*crRHS65 + crRHS61*crRHS68 + crRHS62*crRHS67 - crRHS63*crRHS66 - crRHS63*crRHS69 - crRHS64*crRHS70);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS47 + DN(1,1)*crRHS63 + N[1]*crRHS36);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS42 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS12 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS16 + N[2]*crRHS18 + N[2]*crRHS4 + crRHS24*(DN(2,0)*crRHS22 + DN(2,1)*crRHS23) + crRHS25*crRHS74 + crRHS3*crRHS72 + crRHS33*crRHS73 - crRHS47*crRHS72 - crRHS47*crRHS75 + crRHS5*crRHS71 - crRHS57*crRHS76);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS12 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS42 + DN(2,1)*stress_adj[1] - N[2]*crRHS49 + N[2]*crRHS50 + N[2]*crRHS54 + N[2]*crRHS55 + crRHS24*(DN(2,0)*crRHS23 + DN(2,1)*crRHS60) + crRHS31*crRHS72 + crRHS59*crRHS71 + crRHS61*crRHS74 + crRHS62*crRHS73 - crRHS63*crRHS72 - crRHS63*crRHS75 - crRHS64*crRHS76);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS47 + DN(2,1)*crRHS63 + N[2]*crRHS36);

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
    const double time_coeff = rData.TopOptTimeCoefficient;

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
    const BoundedMatrix<double,3,4> v_conv_ns = rData.Convection_velocity_adj; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
    const BoundedMatrix<double,3,4>& f_adj = rData.BodyForce_adj;
    const array_1d<double,4>& p_adj = rData.Pressure_adj;
    // const array_1d<double,4>& pn_adj = rData.Pressure_adj_OldStep1;
    // const array_1d<double,4>& pnn_adj = rData.Pressure_adj_OldStep2;
    const BoundedMatrix<double,3,4> functional_v = rData.Functional_derivative_velocity;
    const array_1d<double,4>& functional_t = rData.Functional_derivative_transport_scalar;
    const array_1d<double,4>& functional_t_adj = rData.Functional_derivative_transport_scalar_adj;
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
const double crRHS5 = N[0]*functional_v(0,0) + N[1]*functional_v(1,0) + N[2]*functional_v(2,0) + N[3]*functional_v(3,0);
const double crRHS6 = 2.0*crRHS2*functional_weights[0];
const double crRHS7 = N[0]*crRHS6;
const double crRHS8 = DN(0,0)*functional_t[0] + DN(1,0)*functional_t[1] + DN(2,0)*functional_t[2] + DN(3,0)*functional_t[3];
const double crRHS9 = N[0]*t_ConvCoeff[0] + N[1]*t_ConvCoeff[1] + N[2]*t_ConvCoeff[2] + N[3]*t_ConvCoeff[3];
const double crRHS10 = crRHS9*(N[0]*functional_t_adj[0] + N[1]*functional_t_adj[1] + N[2]*functional_t_adj[2] + N[3]*functional_t_adj[3]);
const double crRHS11 = crRHS10*crRHS8;
const double crRHS12 = crRHS9*functional_weights[6]*(N[0]*functional_t[0] + N[1]*functional_t[1] + N[2]*functional_t[2] + N[3]*functional_t[3]);
const double crRHS13 = crRHS12*crRHS8;
const double crRHS14 = N[0]*(bdf0*v_adj(0,0) + bdf1*vn_adj(0,0) + bdf2*vnn_adj(0,0)) + N[1]*(bdf0*v_adj(1,0) + bdf1*vn_adj(1,0) + bdf2*vnn_adj(1,0)) + N[2]*(bdf0*v_adj(2,0) + bdf1*vn_adj(2,0) + bdf2*vnn_adj(2,0)) + N[3]*(bdf0*v_adj(3,0) + bdf1*vn_adj(3,0) + bdf2*vnn_adj(3,0));
const double crRHS15 = N[0]*rho;
const double crRHS16 = crRHS15*time_coeff;
const double crRHS17 = N[0]*v_conv_ns(0,0) + N[1]*v_conv_ns(1,0) + N[2]*v_conv_ns(2,0) + N[3]*v_conv_ns(3,0);
const double crRHS18 = N[0]*v_conv_ns(0,1) + N[1]*v_conv_ns(1,1) + N[2]*v_conv_ns(2,1) + N[3]*v_conv_ns(3,1);
const double crRHS19 = N[0]*v_conv_ns(0,2) + N[1]*v_conv_ns(1,2) + N[2]*v_conv_ns(2,2) + N[3]*v_conv_ns(3,2);
const double crRHS20 = rho*(DN(0,0)*crRHS17 + DN(0,1)*crRHS18 + DN(0,2)*crRHS19);
const double crRHS21 = DN(0,1)*functional_v(0,0);
const double crRHS22 = DN(1,1)*functional_v(1,0);
const double crRHS23 = DN(2,1)*functional_v(2,0);
const double crRHS24 = DN(3,1)*functional_v(3,0);
const double crRHS25 = DN(0,0)*functional_v(0,1) + DN(1,0)*functional_v(1,1) + DN(2,0)*functional_v(2,1) + DN(3,0)*functional_v(3,1);
const double crRHS26 = -crRHS21 - crRHS22 - crRHS23 - crRHS24 + crRHS25;
const double crRHS27 = DN(0,2)*functional_v(0,0);
const double crRHS28 = DN(1,2)*functional_v(1,0);
const double crRHS29 = DN(2,2)*functional_v(2,0);
const double crRHS30 = DN(3,2)*functional_v(3,0);
const double crRHS31 = DN(0,0)*functional_v(0,2) + DN(1,0)*functional_v(1,2) + DN(2,0)*functional_v(2,2) + DN(3,0)*functional_v(3,2);
const double crRHS32 = -crRHS27 - crRHS28 - crRHS29 - crRHS30 + crRHS31;
const double crRHS33 = 2.0*functional_weights[2]*mu;
const double crRHS34 = 1.0*DN(0,0)*functional_v(0,0) + 1.0*DN(1,0)*functional_v(1,0) + 1.0*DN(2,0)*functional_v(2,0) + 1.0*DN(3,0)*functional_v(3,0);
const double crRHS35 = 0.5*crRHS21 + 0.5*crRHS22 + 0.5*crRHS23 + 0.5*crRHS24 + 0.5*crRHS25;
const double crRHS36 = crRHS27 + crRHS28 + crRHS29 + crRHS30 + crRHS31;
const double crRHS37 = 0.5*DN(0,2);
const double crRHS38 = 4.0*functional_weights[1]*mu;
const double crRHS39 = DN(0,0)*v_conv_ns(0,0) + DN(1,0)*v_conv_ns(1,0) + DN(2,0)*v_conv_ns(2,0) + DN(3,0)*v_conv_ns(3,0);
const double crRHS40 = crRHS3*crRHS39;
const double crRHS41 = DN(0,0)*v_conv_ns(0,1) + DN(1,0)*v_conv_ns(1,1) + DN(2,0)*v_conv_ns(2,1) + DN(3,0)*v_conv_ns(3,1);
const double crRHS42 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1) + N[3]*v_adj(3,1);
const double crRHS43 = crRHS41*crRHS42;
const double crRHS44 = DN(0,0)*v_conv_ns(0,2) + DN(1,0)*v_conv_ns(1,2) + DN(2,0)*v_conv_ns(2,2) + DN(3,0)*v_conv_ns(3,2);
const double crRHS45 = N[0]*v_adj(0,2) + N[1]*v_adj(1,2) + N[2]*v_adj(2,2) + N[3]*v_adj(3,2);
const double crRHS46 = crRHS44*crRHS45;
const double crRHS47 = crRHS40 + crRHS43 + crRHS46;
const double crRHS48 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0) + DN(3,0)*v_adj(3,0);
const double crRHS49 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1) + DN(3,1)*v_adj(3,1);
const double crRHS50 = DN(0,2)*v_adj(0,2) + DN(1,2)*v_adj(1,2) + DN(2,2)*v_adj(2,2) + DN(3,2)*v_adj(3,2);
const double crRHS51 = crRHS48 + crRHS49 + crRHS50;
const double crRHS52 = crRHS2*stab_c3;
const double crRHS53 = rho*stab_c2*sqrt(pow(crRHS17, 2) + pow(crRHS18, 2) + pow(crRHS19, 2));
const double crRHS54 = DN(0,1)*v_conv_ns(0,0) + DN(1,1)*v_conv_ns(1,0) + DN(2,1)*v_conv_ns(2,0) + DN(3,1)*v_conv_ns(3,0);
const double crRHS55 = DN(0,1)*v_conv_ns(0,1) + DN(1,1)*v_conv_ns(1,1) + DN(2,1)*v_conv_ns(2,1) + DN(3,1)*v_conv_ns(3,1);
const double crRHS56 = DN(0,1)*v_conv_ns(0,2) + DN(1,1)*v_conv_ns(1,2) + DN(2,1)*v_conv_ns(2,2) + DN(3,1)*v_conv_ns(3,2);
const double crRHS57 = DN(0,2)*v_conv_ns(0,0) + DN(1,2)*v_conv_ns(1,0) + DN(2,2)*v_conv_ns(2,0) + DN(3,2)*v_conv_ns(3,0);
const double crRHS58 = DN(0,2)*v_conv_ns(0,1) + DN(1,2)*v_conv_ns(1,1) + DN(2,2)*v_conv_ns(2,1) + DN(3,2)*v_conv_ns(3,1);
const double crRHS59 = DN(0,2)*v_conv_ns(0,2) + DN(1,2)*v_conv_ns(1,2) + DN(2,2)*v_conv_ns(2,2) + DN(3,2)*v_conv_ns(3,2);
const double crRHS60 = rho*stab_c3*sqrt(pow(crRHS39, 2) + pow(crRHS41, 2) + pow(crRHS44, 2) + pow(crRHS54, 2) + pow(crRHS55, 2) + pow(crRHS56, 2) + pow(crRHS57, 2) + pow(crRHS58, 2) + pow(crRHS59, 2));
const double crRHS61 = crRHS51*(h*(crRHS52*h + crRHS53 + crRHS60*h)/stab_c1 + mu);
const double crRHS62 = crRHS17*rho;
const double crRHS63 = crRHS18*rho;
const double crRHS64 = crRHS19*rho;
const double crRHS65 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 + crRHS11 + crRHS13 + crRHS17*crRHS6 + crRHS4 + crRHS40*rho + crRHS43*rho + crRHS46*rho - crRHS48*crRHS62 - crRHS63*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS64*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0));
const double crRHS66 = 1.0/(crRHS52 + crRHS53/h + crRHS60 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS67 = crRHS65*crRHS66;
const double crRHS68 = N[0]*crRHS2;
const double crRHS69 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS70 = crRHS2*crRHS42;
const double crRHS71 = crRHS3*crRHS54;
const double crRHS72 = crRHS42*crRHS55;
const double crRHS73 = crRHS45*crRHS56;
const double crRHS74 = DN(0,1)*functional_t[0] + DN(1,1)*functional_t[1] + DN(2,1)*functional_t[2] + DN(3,1)*functional_t[3];
const double crRHS75 = crRHS10*crRHS74;
const double crRHS76 = crRHS12*crRHS74;
const double crRHS77 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] + crRHS18*crRHS6 - crRHS49*crRHS63 - crRHS62*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS64*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS69 + crRHS70 + crRHS71*rho + crRHS72*rho + crRHS73*rho + crRHS75 + crRHS76;
const double crRHS78 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS79 = crRHS2*crRHS45;
const double crRHS80 = crRHS3*crRHS57;
const double crRHS81 = crRHS42*crRHS58;
const double crRHS82 = crRHS45*crRHS59;
const double crRHS83 = DN(0,2)*functional_t[0] + DN(1,2)*functional_t[1] + DN(2,2)*functional_t[2] + DN(3,2)*functional_t[3];
const double crRHS84 = crRHS10*crRHS83;
const double crRHS85 = crRHS12*crRHS83;
const double crRHS86 = DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] + crRHS19*crRHS6 - crRHS50*crRHS64 - crRHS62*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS63*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS78 + crRHS79 + crRHS80*rho + crRHS81*rho + crRHS82*rho + crRHS84 + crRHS85;
const double crRHS87 = crRHS39*crRHS65 + crRHS41*crRHS77 + crRHS44*crRHS86;
const double crRHS88 = crRHS15*crRHS66;
const double crRHS89 = N[0]*functional_v(0,1) + N[1]*functional_v(1,1) + N[2]*functional_v(2,1) + N[3]*functional_v(3,1);
const double crRHS90 = N[0]*(bdf0*v_adj(0,1) + bdf1*vn_adj(0,1) + bdf2*vnn_adj(0,1)) + N[1]*(bdf0*v_adj(1,1) + bdf1*vn_adj(1,1) + bdf2*vnn_adj(1,1)) + N[2]*(bdf0*v_adj(2,1) + bdf1*vn_adj(2,1) + bdf2*vnn_adj(2,1)) + N[3]*(bdf0*v_adj(3,1) + bdf1*vn_adj(3,1) + bdf2*vnn_adj(3,1));
const double crRHS91 = DN(0,2)*functional_v(0,1);
const double crRHS92 = DN(1,2)*functional_v(1,1);
const double crRHS93 = DN(2,2)*functional_v(2,1);
const double crRHS94 = DN(3,2)*functional_v(3,1);
const double crRHS95 = DN(0,1)*functional_v(0,2) + DN(1,1)*functional_v(1,2) + DN(2,1)*functional_v(2,2) + DN(3,1)*functional_v(3,2);
const double crRHS96 = -crRHS91 - crRHS92 - crRHS93 - crRHS94 + crRHS95;
const double crRHS97 = 1.0*DN(0,1)*functional_v(0,1) + 1.0*DN(1,1)*functional_v(1,1) + 1.0*DN(2,1)*functional_v(2,1) + 1.0*DN(3,1)*functional_v(3,1);
const double crRHS98 = crRHS91 + crRHS92 + crRHS93 + crRHS94 + crRHS95;
const double crRHS99 = crRHS71 + crRHS72 + crRHS73;
const double crRHS100 = crRHS66*crRHS77;
const double crRHS101 = crRHS54*crRHS65 + crRHS55*crRHS77 + crRHS56*crRHS86;
const double crRHS102 = N[0]*functional_v(0,2) + N[1]*functional_v(1,2) + N[2]*functional_v(2,2) + N[3]*functional_v(3,2);
const double crRHS103 = N[0]*(bdf0*v_adj(0,2) + bdf1*vn_adj(0,2) + bdf2*vnn_adj(0,2)) + N[1]*(bdf0*v_adj(1,2) + bdf1*vn_adj(1,2) + bdf2*vnn_adj(1,2)) + N[2]*(bdf0*v_adj(2,2) + bdf1*vn_adj(2,2) + bdf2*vnn_adj(2,2)) + N[3]*(bdf0*v_adj(3,2) + bdf1*vn_adj(3,2) + bdf2*vnn_adj(3,2));
const double crRHS104 = 1.0*DN(0,2)*functional_v(0,2) + 1.0*DN(1,2)*functional_v(1,2) + 1.0*DN(2,2)*functional_v(2,2) + 1.0*DN(3,2)*functional_v(3,2);
const double crRHS105 = 0.5*crRHS36;
const double crRHS106 = 0.5*crRHS98;
const double crRHS107 = crRHS80 + crRHS81 + crRHS82;
const double crRHS108 = crRHS66*crRHS86;
const double crRHS109 = crRHS57*crRHS65 + crRHS58*crRHS77 + crRHS59*crRHS86;
const double crRHS110 = N[1]*crRHS6;
const double crRHS111 = N[1]*rho;
const double crRHS112 = crRHS111*time_coeff;
const double crRHS113 = rho*(DN(1,0)*crRHS17 + DN(1,1)*crRHS18 + DN(1,2)*crRHS19);
const double crRHS114 = N[1]*crRHS2;
const double crRHS115 = crRHS111*crRHS66;
const double crRHS116 = N[2]*crRHS6;
const double crRHS117 = N[2]*rho;
const double crRHS118 = crRHS117*time_coeff;
const double crRHS119 = rho*(DN(2,0)*crRHS17 + DN(2,1)*crRHS18 + DN(2,2)*crRHS19);
const double crRHS120 = N[2]*crRHS2;
const double crRHS121 = crRHS117*crRHS66;
const double crRHS122 = N[3]*crRHS6;
const double crRHS123 = N[3]*rho;
const double crRHS124 = crRHS123*time_coeff;
const double crRHS125 = rho*(DN(3,0)*crRHS17 + DN(3,1)*crRHS18 + DN(3,2)*crRHS19);
const double crRHS126 = N[3]*crRHS2;
const double crRHS127 = crRHS123*crRHS66;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS61 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS11 + N[0]*crRHS13 + N[0]*crRHS4 + crRHS14*crRHS16 + crRHS15*crRHS47 + crRHS20*crRHS3 - crRHS20*crRHS67 - crRHS33*(DN(0,1)*crRHS26 + DN(0,2)*crRHS32) + crRHS38*(DN(0,0)*crRHS34 + DN(0,1)*crRHS35 + crRHS36*crRHS37) + crRHS5*crRHS7 - crRHS67*crRHS68 - crRHS87*crRHS88);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS61 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS69 + N[0]*crRHS70 + N[0]*crRHS75 + N[0]*crRHS76 - crRHS100*crRHS20 - crRHS100*crRHS68 - crRHS101*crRHS88 + crRHS15*crRHS99 + crRHS16*crRHS90 + crRHS20*crRHS42 + crRHS33*(DN(0,0)*crRHS26 - DN(0,2)*crRHS96) + crRHS38*(DN(0,0)*crRHS35 + DN(0,1)*crRHS97 + crRHS37*crRHS98) + crRHS7*crRHS89);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS61 + DN(0,2)*stress_adj[2] - N[0]*crRHS78 + N[0]*crRHS79 + N[0]*crRHS84 + N[0]*crRHS85 + crRHS102*crRHS7 + crRHS103*crRHS16 + crRHS107*crRHS15 - crRHS108*crRHS20 - crRHS108*crRHS68 - crRHS109*crRHS88 + crRHS20*crRHS45 + crRHS33*(DN(0,0)*crRHS32 + DN(0,1)*crRHS96) + crRHS38*(DN(0,0)*crRHS105 + DN(0,1)*crRHS106 + DN(0,2)*crRHS104));
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS67 + DN(0,1)*crRHS100 + DN(0,2)*crRHS108 + N[0]*crRHS51);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS61 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS11 + N[1]*crRHS13 + N[1]*crRHS4 + crRHS110*crRHS5 + crRHS111*crRHS47 + crRHS112*crRHS14 + crRHS113*crRHS3 - crRHS113*crRHS67 - crRHS114*crRHS67 - crRHS115*crRHS87 - crRHS33*(DN(1,1)*crRHS26 + DN(1,2)*crRHS32) + crRHS38*(DN(1,0)*crRHS34 + DN(1,1)*crRHS35 + DN(1,2)*crRHS105));
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS61 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS69 + N[1]*crRHS70 + N[1]*crRHS75 + N[1]*crRHS76 - crRHS100*crRHS113 - crRHS100*crRHS114 - crRHS101*crRHS115 + crRHS110*crRHS89 + crRHS111*crRHS99 + crRHS112*crRHS90 + crRHS113*crRHS42 + crRHS33*(DN(1,0)*crRHS26 - DN(1,2)*crRHS96) + crRHS38*(DN(1,0)*crRHS35 + DN(1,1)*crRHS97 + DN(1,2)*crRHS106));
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS61 + DN(1,2)*stress_adj[2] - N[1]*crRHS78 + N[1]*crRHS79 + N[1]*crRHS84 + N[1]*crRHS85 + crRHS102*crRHS110 + crRHS103*crRHS112 + crRHS107*crRHS111 - crRHS108*crRHS113 - crRHS108*crRHS114 - crRHS109*crRHS115 + crRHS113*crRHS45 + crRHS33*(DN(1,0)*crRHS32 + DN(1,1)*crRHS96) + crRHS38*(DN(1,0)*crRHS105 + DN(1,1)*crRHS106 + DN(1,2)*crRHS104));
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS67 + DN(1,1)*crRHS100 + DN(1,2)*crRHS108 + N[1]*crRHS51);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS61 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS11 + N[2]*crRHS13 + N[2]*crRHS4 + crRHS116*crRHS5 + crRHS117*crRHS47 + crRHS118*crRHS14 + crRHS119*crRHS3 - crRHS119*crRHS67 - crRHS120*crRHS67 - crRHS121*crRHS87 - crRHS33*(DN(2,1)*crRHS26 + DN(2,2)*crRHS32) + crRHS38*(DN(2,0)*crRHS34 + DN(2,1)*crRHS35 + DN(2,2)*crRHS105));
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS61 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS69 + N[2]*crRHS70 + N[2]*crRHS75 + N[2]*crRHS76 - crRHS100*crRHS119 - crRHS100*crRHS120 - crRHS101*crRHS121 + crRHS116*crRHS89 + crRHS117*crRHS99 + crRHS118*crRHS90 + crRHS119*crRHS42 + crRHS33*(DN(2,0)*crRHS26 - DN(2,2)*crRHS96) + crRHS38*(DN(2,0)*crRHS35 + DN(2,1)*crRHS97 + DN(2,2)*crRHS106));
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS61 + DN(2,2)*stress_adj[2] - N[2]*crRHS78 + N[2]*crRHS79 + N[2]*crRHS84 + N[2]*crRHS85 + crRHS102*crRHS116 + crRHS103*crRHS118 + crRHS107*crRHS117 - crRHS108*crRHS119 - crRHS108*crRHS120 - crRHS109*crRHS121 + crRHS119*crRHS45 + crRHS33*(DN(2,0)*crRHS32 + DN(2,1)*crRHS96) + crRHS38*(DN(2,0)*crRHS105 + DN(2,1)*crRHS106 + DN(2,2)*crRHS104));
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS67 + DN(2,1)*crRHS100 + DN(2,2)*crRHS108 + N[2]*crRHS51);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS61 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS11 + N[3]*crRHS13 + N[3]*crRHS4 + crRHS122*crRHS5 + crRHS123*crRHS47 + crRHS124*crRHS14 + crRHS125*crRHS3 - crRHS125*crRHS67 - crRHS126*crRHS67 - crRHS127*crRHS87 - crRHS33*(DN(3,1)*crRHS26 + DN(3,2)*crRHS32) + crRHS38*(DN(3,0)*crRHS34 + DN(3,1)*crRHS35 + DN(3,2)*crRHS105));
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS61 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS69 + N[3]*crRHS70 + N[3]*crRHS75 + N[3]*crRHS76 - crRHS100*crRHS125 - crRHS100*crRHS126 - crRHS101*crRHS127 + crRHS122*crRHS89 + crRHS123*crRHS99 + crRHS124*crRHS90 + crRHS125*crRHS42 + crRHS33*(DN(3,0)*crRHS26 - DN(3,2)*crRHS96) + crRHS38*(DN(3,0)*crRHS35 + DN(3,1)*crRHS97 + DN(3,2)*crRHS106));
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS61 + DN(3,2)*stress_adj[2] - N[3]*crRHS78 + N[3]*crRHS79 + N[3]*crRHS84 + N[3]*crRHS85 + crRHS102*crRHS122 + crRHS103*crRHS124 + crRHS107*crRHS123 - crRHS108*crRHS125 - crRHS108*crRHS126 - crRHS109*crRHS127 + crRHS125*crRHS45 + crRHS33*(DN(3,0)*crRHS32 + DN(3,1)*crRHS96) + crRHS38*(DN(3,0)*crRHS105 + DN(3,1)*crRHS106 + DN(3,2)*crRHS104));
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS67 + DN(3,1)*crRHS100 + DN(3,2)*crRHS108 + N[3]*crRHS51);

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