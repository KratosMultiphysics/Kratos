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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    
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
const double crLHS10 = DN(0,0)*crLHS6;
const double crLHS11 = DN(0,1)*crLHS7;
const double crLHS12 = crLHS10 + crLHS11;
const double crLHS13 = N[0]*rho;
const double crLHS14 = N[0]*crLHS4;
const double crLHS15 = crLHS10*rho;
const double crLHS16 = crLHS11*rho;
const double crLHS17 = crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = 1.0/(crLHS5 + crLHS8/h + mu*stab_c1/pow(h, 2));
const double crLHS19 = 1.0*crLHS18;
const double crLHS20 = crLHS19*rho;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = 1.0*crLHS14;
const double crLHS23 = crLHS18*crLHS22;
const double crLHS24 = crLHS17*crLHS19;
const double crLHS25 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS26 = crLHS13*crLHS25;
const double crLHS27 = pow(N[0], 2)*crLHS4 + crLHS12*crLHS13 + crLHS17*crLHS21 - crLHS17*crLHS23 + crLHS24*crLHS26;
const double crLHS28 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS29 = C(1,2)*DN(0,1);
const double crLHS30 = C(2,2)*DN(0,0) + crLHS29;
const double crLHS31 = DN(0,0)*crLHS9;
const double crLHS32 = DN(0,1)*crLHS31;
const double crLHS33 = crLHS19*crLHS25;
const double crLHS34 = gauss_weight*(-N[0] + crLHS13*crLHS33 + crLHS21 - crLHS23);
const double crLHS35 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS36 = C(0,2)*DN(1,0);
const double crLHS37 = C(2,2)*DN(1,1) + crLHS36;
const double crLHS38 = DN(0,0)*DN(1,0);
const double crLHS39 = N[1]*crLHS14;
const double crLHS40 = crLHS38*crLHS9 + crLHS39;
const double crLHS41 = DN(1,0)*crLHS6;
const double crLHS42 = DN(1,1)*crLHS7;
const double crLHS43 = crLHS41 + crLHS42;
const double crLHS44 = N[1]*crLHS4;
const double crLHS45 = crLHS41*rho;
const double crLHS46 = crLHS42*rho;
const double crLHS47 = crLHS44 + crLHS45 + crLHS46;
const double crLHS48 = crLHS19*crLHS47;
const double crLHS49 = crLHS13*crLHS43 + crLHS21*crLHS47 - crLHS23*crLHS47 + crLHS26*crLHS48;
const double crLHS50 = C(0,1)*DN(1,1) + crLHS36;
const double crLHS51 = C(1,2)*DN(1,1);
const double crLHS52 = C(2,2)*DN(1,0) + crLHS51;
const double crLHS53 = DN(1,1)*crLHS31;
const double crLHS54 = DN(0,0)*N[1];
const double crLHS55 = DN(1,0)*N[0];
const double crLHS56 = crLHS20*crLHS25;
const double crLHS57 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS58 = C(0,2)*DN(2,0);
const double crLHS59 = C(2,2)*DN(2,1) + crLHS58;
const double crLHS60 = DN(0,0)*DN(2,0);
const double crLHS61 = N[2]*crLHS14;
const double crLHS62 = crLHS60*crLHS9 + crLHS61;
const double crLHS63 = DN(2,0)*crLHS6;
const double crLHS64 = DN(2,1)*crLHS7;
const double crLHS65 = crLHS63 + crLHS64;
const double crLHS66 = N[2]*crLHS4;
const double crLHS67 = crLHS63*rho;
const double crLHS68 = crLHS64*rho;
const double crLHS69 = crLHS66 + crLHS67 + crLHS68;
const double crLHS70 = crLHS19*crLHS69;
const double crLHS71 = crLHS13*crLHS65 + crLHS21*crLHS69 - crLHS23*crLHS69 + crLHS26*crLHS70;
const double crLHS72 = C(0,1)*DN(2,1) + crLHS58;
const double crLHS73 = C(1,2)*DN(2,1);
const double crLHS74 = C(2,2)*DN(2,0) + crLHS73;
const double crLHS75 = DN(2,1)*crLHS31;
const double crLHS76 = DN(0,0)*N[2];
const double crLHS77 = DN(2,0)*N[0];
const double crLHS78 = C(0,1)*DN(0,0) + crLHS29;
const double crLHS79 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS80 = pow(DN(0,1), 2);
const double crLHS81 = C(0,1)*DN(1,0) + crLHS51;
const double crLHS82 = DN(0,1)*crLHS9;
const double crLHS83 = DN(1,0)*crLHS82;
const double crLHS84 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS85 = DN(0,1)*DN(1,1);
const double crLHS86 = crLHS39 + crLHS85*crLHS9;
const double crLHS87 = DN(0,1)*N[1];
const double crLHS88 = DN(1,1)*N[0];
const double crLHS89 = C(0,1)*DN(2,0) + crLHS73;
const double crLHS90 = DN(2,0)*crLHS82;
const double crLHS91 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS92 = DN(0,1)*DN(2,1);
const double crLHS93 = crLHS61 + crLHS9*crLHS92;
const double crLHS94 = DN(0,1)*N[2];
const double crLHS95 = DN(2,1)*N[0];
const double crLHS96 = gauss_weight*(N[0] + crLHS18*(1.0*crLHS15 + 1.0*crLHS16 + crLHS22));
const double crLHS97 = crLHS19*gauss_weight;
const double crLHS98 = crLHS97*(crLHS38 + crLHS85);
const double crLHS99 = crLHS97*(crLHS60 + crLHS92);
const double crLHS100 = N[1]*rho;
const double crLHS101 = crLHS20*crLHS43;
const double crLHS102 = 1.0*crLHS44;
const double crLHS103 = crLHS102*crLHS18;
const double crLHS104 = crLHS100*crLHS25;
const double crLHS105 = crLHS100*crLHS12 + crLHS101*crLHS17 - crLHS103*crLHS17 + crLHS104*crLHS24;
const double crLHS106 = pow(DN(1,0), 2);
const double crLHS107 = pow(N[1], 2)*crLHS4 + crLHS100*crLHS43 + crLHS101*crLHS47 - crLHS103*crLHS47 + crLHS104*crLHS48;
const double crLHS108 = DN(1,0)*crLHS9;
const double crLHS109 = DN(1,1)*crLHS108;
const double crLHS110 = gauss_weight*(-N[1] + crLHS100*crLHS33 + crLHS101 - crLHS103);
const double crLHS111 = DN(1,0)*DN(2,0);
const double crLHS112 = N[2]*crLHS44;
const double crLHS113 = crLHS111*crLHS9 + crLHS112;
const double crLHS114 = crLHS100*crLHS65 + crLHS101*crLHS69 - crLHS103*crLHS69 + crLHS104*crLHS70;
const double crLHS115 = DN(2,1)*crLHS108;
const double crLHS116 = DN(1,0)*N[2];
const double crLHS117 = DN(2,0)*N[1];
const double crLHS118 = pow(DN(1,1), 2);
const double crLHS119 = DN(2,0)*crLHS9;
const double crLHS120 = DN(1,1)*crLHS119;
const double crLHS121 = DN(1,1)*DN(2,1);
const double crLHS122 = crLHS112 + crLHS121*crLHS9;
const double crLHS123 = DN(1,1)*N[2];
const double crLHS124 = DN(2,1)*N[1];
const double crLHS125 = gauss_weight*(N[1] + crLHS18*(crLHS102 + 1.0*crLHS45 + 1.0*crLHS46));
const double crLHS126 = crLHS97*(crLHS111 + crLHS121);
const double crLHS127 = N[2]*rho;
const double crLHS128 = crLHS20*crLHS65;
const double crLHS129 = 1.0*crLHS66;
const double crLHS130 = crLHS129*crLHS18;
const double crLHS131 = crLHS127*crLHS25;
const double crLHS132 = crLHS12*crLHS127 + crLHS128*crLHS17 - crLHS130*crLHS17 + crLHS131*crLHS24;
const double crLHS133 = crLHS127*crLHS43 + crLHS128*crLHS47 - crLHS130*crLHS47 + crLHS131*crLHS48;
const double crLHS134 = pow(DN(2,0), 2);
const double crLHS135 = pow(N[2], 2)*crLHS4 + crLHS127*crLHS65 + crLHS128*crLHS69 - crLHS130*crLHS69 + crLHS131*crLHS70;
const double crLHS136 = DN(2,1)*crLHS119;
const double crLHS137 = gauss_weight*(-N[2] + crLHS127*crLHS33 + crLHS128 - crLHS130);
const double crLHS138 = pow(DN(2,1), 2);
const double crLHS139 = gauss_weight*(N[2] + crLHS18*(crLHS129 + 1.0*crLHS67 + 1.0*crLHS68));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS27 + crLHS3*crLHS9);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS28 + DN(0,1)*crLHS30 + crLHS32);
rLHS(0,2)+=DN(0,0)*crLHS34;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS35 + DN(0,1)*crLHS37 + crLHS40 + crLHS49);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS50 + DN(0,1)*crLHS52 + crLHS53);
rLHS(0,5)+=-gauss_weight*(-DN(1,0)*crLHS21 + DN(1,0)*crLHS23 + crLHS54 - crLHS55*crLHS56);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS57 + DN(0,1)*crLHS59 + crLHS62 + crLHS71);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS72 + DN(0,1)*crLHS74 + crLHS75);
rLHS(0,8)+=-gauss_weight*(-DN(2,0)*crLHS21 + DN(2,0)*crLHS23 - crLHS56*crLHS77 + crLHS76);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS78 + crLHS32);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS30 + DN(0,1)*crLHS79 + crLHS27 + crLHS80*crLHS9);
rLHS(1,2)+=DN(0,1)*crLHS34;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS37 + DN(0,1)*crLHS81 + crLHS83);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS52 + DN(0,1)*crLHS84 + crLHS49 + crLHS86);
rLHS(1,5)+=-gauss_weight*(-DN(1,1)*crLHS21 + DN(1,1)*crLHS23 - crLHS56*crLHS88 + crLHS87);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS59 + DN(0,1)*crLHS89 + crLHS90);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS91 + crLHS71 + crLHS93);
rLHS(1,8)+=-gauss_weight*(-DN(2,1)*crLHS21 + DN(2,1)*crLHS23 - crLHS56*crLHS95 + crLHS94);
rLHS(2,0)+=DN(0,0)*crLHS96;
rLHS(2,1)+=DN(0,1)*crLHS96;
rLHS(2,2)+=crLHS97*(crLHS3 + crLHS80);
rLHS(2,3)+=gauss_weight*(DN(0,0)*crLHS48 + crLHS55);
rLHS(2,4)+=gauss_weight*(DN(0,1)*crLHS48 + crLHS88);
rLHS(2,5)+=crLHS98;
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS70 + crLHS77);
rLHS(2,7)+=gauss_weight*(DN(0,1)*crLHS70 + crLHS95);
rLHS(2,8)+=crLHS99;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS105 + crLHS40);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS28 + DN(1,1)*crLHS30 + crLHS83);
rLHS(3,2)+=gauss_weight*(DN(0,0)*crLHS101 - DN(0,0)*crLHS103 + crLHS54*crLHS56 - crLHS55);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS35 + DN(1,1)*crLHS37 + crLHS106*crLHS9 + crLHS107);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS50 + DN(1,1)*crLHS52 + crLHS109);
rLHS(3,5)+=DN(1,0)*crLHS110;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS57 + DN(1,1)*crLHS59 + crLHS113 + crLHS114);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS72 + DN(1,1)*crLHS74 + crLHS115);
rLHS(3,8)+=-gauss_weight*(-DN(2,0)*crLHS101 + DN(2,0)*crLHS103 + crLHS116 - crLHS117*crLHS56);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS78 + crLHS53);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS30 + DN(1,1)*crLHS79 + crLHS105 + crLHS86);
rLHS(4,2)+=gauss_weight*(DN(0,1)*crLHS101 - DN(0,1)*crLHS103 + crLHS56*crLHS87 - crLHS88);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS37 + DN(1,1)*crLHS81 + crLHS109);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS52 + DN(1,1)*crLHS84 + crLHS107 + crLHS118*crLHS9);
rLHS(4,5)+=DN(1,1)*crLHS110;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS59 + DN(1,1)*crLHS89 + crLHS120);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS91 + crLHS114 + crLHS122);
rLHS(4,8)+=-gauss_weight*(-DN(2,1)*crLHS101 + DN(2,1)*crLHS103 + crLHS123 - crLHS124*crLHS56);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS24 + crLHS54);
rLHS(5,1)+=gauss_weight*(DN(1,1)*crLHS24 + crLHS87);
rLHS(5,2)+=crLHS98;
rLHS(5,3)+=DN(1,0)*crLHS125;
rLHS(5,4)+=DN(1,1)*crLHS125;
rLHS(5,5)+=crLHS97*(crLHS106 + crLHS118);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS70 + crLHS117);
rLHS(5,7)+=gauss_weight*(DN(1,1)*crLHS70 + crLHS124);
rLHS(5,8)+=crLHS126;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS132 + crLHS62);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS28 + DN(2,1)*crLHS30 + crLHS90);
rLHS(6,2)+=gauss_weight*(DN(0,0)*crLHS128 - DN(0,0)*crLHS130 + crLHS56*crLHS76 - crLHS77);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS35 + DN(2,1)*crLHS37 + crLHS113 + crLHS133);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS50 + DN(2,1)*crLHS52 + crLHS120);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS128 - DN(1,0)*crLHS130 + crLHS116*crLHS56 - crLHS117);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS57 + DN(2,1)*crLHS59 + crLHS134*crLHS9 + crLHS135);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS72 + DN(2,1)*crLHS74 + crLHS136);
rLHS(6,8)+=DN(2,0)*crLHS137;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS78 + crLHS75);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS30 + DN(2,1)*crLHS79 + crLHS132 + crLHS93);
rLHS(7,2)+=gauss_weight*(DN(0,1)*crLHS128 - DN(0,1)*crLHS130 + crLHS56*crLHS94 - crLHS95);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS37 + DN(2,1)*crLHS81 + crLHS115);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS52 + DN(2,1)*crLHS84 + crLHS122 + crLHS133);
rLHS(7,5)+=gauss_weight*(DN(1,1)*crLHS128 - DN(1,1)*crLHS130 + crLHS123*crLHS56 - crLHS124);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS59 + DN(2,1)*crLHS89 + crLHS136);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS91 + crLHS135 + crLHS138*crLHS9);
rLHS(7,8)+=DN(2,1)*crLHS137;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS24 + crLHS76);
rLHS(8,1)+=gauss_weight*(DN(2,1)*crLHS24 + crLHS94);
rLHS(8,2)+=crLHS99;
rLHS(8,3)+=gauss_weight*(DN(2,0)*crLHS48 + crLHS116);
rLHS(8,4)+=gauss_weight*(DN(2,1)*crLHS48 + crLHS123);
rLHS(8,5)+=crLHS126;
rLHS(8,6)+=DN(2,0)*crLHS139;
rLHS(8,7)+=DN(2,1)*crLHS139;
rLHS(8,8)+=crLHS97*(crLHS134 + crLHS138);
    
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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    
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
const double crLHS13 = DN(0,0)*crLHS8;
const double crLHS14 = DN(0,1)*crLHS9;
const double crLHS15 = DN(0,2)*crLHS10;
const double crLHS16 = crLHS13 + crLHS14 + crLHS15;
const double crLHS17 = N[0]*rho;
const double crLHS18 = N[0]*crLHS6;
const double crLHS19 = crLHS13*rho;
const double crLHS20 = crLHS14*rho;
const double crLHS21 = crLHS15*rho;
const double crLHS22 = crLHS18 + crLHS19 + crLHS20 + crLHS21;
const double crLHS23 = 1.0/(crLHS11/h + crLHS7 + mu*stab_c1/pow(h, 2));
const double crLHS24 = 1.0*crLHS23;
const double crLHS25 = crLHS24*rho;
const double crLHS26 = crLHS16*crLHS25;
const double crLHS27 = 1.0*crLHS18;
const double crLHS28 = crLHS23*crLHS27;
const double crLHS29 = crLHS22*crLHS24;
const double crLHS30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS31 = crLHS17*crLHS30;
const double crLHS32 = pow(N[0], 2)*crLHS6 + crLHS16*crLHS17 + crLHS22*crLHS26 - crLHS22*crLHS28 + crLHS29*crLHS31;
const double crLHS33 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS34 = C(1,3)*DN(0,1);
const double crLHS35 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS34;
const double crLHS36 = C(3,5)*DN(0,0);
const double crLHS37 = C(4,5)*DN(0,2);
const double crLHS38 = C(1,5)*DN(0,1) + crLHS36 + crLHS37;
const double crLHS39 = DN(0,0)*crLHS12;
const double crLHS40 = DN(0,1)*crLHS39;
const double crLHS41 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS42 = C(3,4)*DN(0,1);
const double crLHS43 = C(2,3)*DN(0,2) + crLHS36 + crLHS42;
const double crLHS44 = C(2,5)*DN(0,2);
const double crLHS45 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS44;
const double crLHS46 = DN(0,2)*crLHS39;
const double crLHS47 = crLHS24*crLHS30;
const double crLHS48 = gauss_weight*(-N[0] + crLHS17*crLHS47 + crLHS26 - crLHS28);
const double crLHS49 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS50 = C(0,3)*DN(1,0);
const double crLHS51 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS50;
const double crLHS52 = C(0,5)*DN(1,0);
const double crLHS53 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS52;
const double crLHS54 = DN(0,0)*DN(1,0);
const double crLHS55 = N[1]*crLHS18;
const double crLHS56 = crLHS12*crLHS54 + crLHS55;
const double crLHS57 = DN(1,0)*crLHS8;
const double crLHS58 = DN(1,1)*crLHS9;
const double crLHS59 = DN(1,2)*crLHS10;
const double crLHS60 = crLHS57 + crLHS58 + crLHS59;
const double crLHS61 = N[1]*crLHS6;
const double crLHS62 = crLHS57*rho;
const double crLHS63 = crLHS58*rho;
const double crLHS64 = crLHS59*rho;
const double crLHS65 = crLHS61 + crLHS62 + crLHS63 + crLHS64;
const double crLHS66 = crLHS24*crLHS65;
const double crLHS67 = crLHS17*crLHS60 + crLHS26*crLHS65 - crLHS28*crLHS65 + crLHS31*crLHS66;
const double crLHS68 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS50;
const double crLHS69 = C(1,3)*DN(1,1);
const double crLHS70 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS69;
const double crLHS71 = C(3,5)*DN(1,0);
const double crLHS72 = C(4,5)*DN(1,2);
const double crLHS73 = C(1,5)*DN(1,1) + crLHS71 + crLHS72;
const double crLHS74 = DN(1,1)*crLHS39;
const double crLHS75 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS52;
const double crLHS76 = C(3,4)*DN(1,1);
const double crLHS77 = C(2,3)*DN(1,2) + crLHS71 + crLHS76;
const double crLHS78 = C(2,5)*DN(1,2);
const double crLHS79 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS78;
const double crLHS80 = DN(1,2)*crLHS39;
const double crLHS81 = DN(0,0)*N[1];
const double crLHS82 = DN(1,0)*N[0];
const double crLHS83 = crLHS25*crLHS30;
const double crLHS84 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS85 = C(0,3)*DN(2,0);
const double crLHS86 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS85;
const double crLHS87 = C(0,5)*DN(2,0);
const double crLHS88 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS87;
const double crLHS89 = DN(0,0)*DN(2,0);
const double crLHS90 = N[2]*crLHS18;
const double crLHS91 = crLHS12*crLHS89 + crLHS90;
const double crLHS92 = DN(2,0)*crLHS8;
const double crLHS93 = DN(2,1)*crLHS9;
const double crLHS94 = DN(2,2)*crLHS10;
const double crLHS95 = crLHS92 + crLHS93 + crLHS94;
const double crLHS96 = N[2]*crLHS6;
const double crLHS97 = crLHS92*rho;
const double crLHS98 = crLHS93*rho;
const double crLHS99 = crLHS94*rho;
const double crLHS100 = crLHS96 + crLHS97 + crLHS98 + crLHS99;
const double crLHS101 = crLHS100*crLHS24;
const double crLHS102 = crLHS100*crLHS26 - crLHS100*crLHS28 + crLHS101*crLHS31 + crLHS17*crLHS95;
const double crLHS103 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS85;
const double crLHS104 = C(1,3)*DN(2,1);
const double crLHS105 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS104;
const double crLHS106 = C(3,5)*DN(2,0);
const double crLHS107 = C(4,5)*DN(2,2);
const double crLHS108 = C(1,5)*DN(2,1) + crLHS106 + crLHS107;
const double crLHS109 = DN(2,1)*crLHS39;
const double crLHS110 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS87;
const double crLHS111 = C(3,4)*DN(2,1);
const double crLHS112 = C(2,3)*DN(2,2) + crLHS106 + crLHS111;
const double crLHS113 = C(2,5)*DN(2,2);
const double crLHS114 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS113;
const double crLHS115 = DN(2,2)*crLHS39;
const double crLHS116 = DN(0,0)*N[2];
const double crLHS117 = DN(2,0)*N[0];
const double crLHS118 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS119 = C(0,3)*DN(3,0);
const double crLHS120 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS119;
const double crLHS121 = C(0,5)*DN(3,0);
const double crLHS122 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS121;
const double crLHS123 = DN(0,0)*DN(3,0);
const double crLHS124 = N[3]*crLHS18;
const double crLHS125 = crLHS12*crLHS123 + crLHS124;
const double crLHS126 = DN(3,0)*crLHS8;
const double crLHS127 = DN(3,1)*crLHS9;
const double crLHS128 = DN(3,2)*crLHS10;
const double crLHS129 = crLHS126 + crLHS127 + crLHS128;
const double crLHS130 = N[3]*crLHS6;
const double crLHS131 = crLHS126*rho;
const double crLHS132 = crLHS127*rho;
const double crLHS133 = crLHS128*rho;
const double crLHS134 = crLHS130 + crLHS131 + crLHS132 + crLHS133;
const double crLHS135 = crLHS134*crLHS24;
const double crLHS136 = crLHS129*crLHS17 + crLHS134*crLHS26 - crLHS134*crLHS28 + crLHS135*crLHS31;
const double crLHS137 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS119;
const double crLHS138 = C(1,3)*DN(3,1);
const double crLHS139 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS138;
const double crLHS140 = C(3,5)*DN(3,0);
const double crLHS141 = C(4,5)*DN(3,2);
const double crLHS142 = C(1,5)*DN(3,1) + crLHS140 + crLHS141;
const double crLHS143 = DN(3,1)*crLHS39;
const double crLHS144 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS121;
const double crLHS145 = C(3,4)*DN(3,1);
const double crLHS146 = C(2,3)*DN(3,2) + crLHS140 + crLHS145;
const double crLHS147 = C(2,5)*DN(3,2);
const double crLHS148 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS147;
const double crLHS149 = DN(3,2)*crLHS39;
const double crLHS150 = DN(0,0)*N[3];
const double crLHS151 = DN(3,0)*N[0];
const double crLHS152 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS34;
const double crLHS153 = C(0,4)*DN(0,0) + crLHS37 + crLHS42;
const double crLHS154 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS155 = C(1,4)*DN(0,1);
const double crLHS156 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS155;
const double crLHS157 = pow(DN(0,1), 2);
const double crLHS158 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS155;
const double crLHS159 = C(2,4)*DN(0,2);
const double crLHS160 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS159;
const double crLHS161 = DN(0,1)*crLHS12;
const double crLHS162 = DN(0,2)*crLHS161;
const double crLHS163 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS69;
const double crLHS164 = C(0,4)*DN(1,0) + crLHS72 + crLHS76;
const double crLHS165 = DN(1,0)*crLHS161;
const double crLHS166 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS167 = C(1,4)*DN(1,1);
const double crLHS168 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS167;
const double crLHS169 = DN(0,1)*DN(1,1);
const double crLHS170 = crLHS12*crLHS169;
const double crLHS171 = crLHS55 + crLHS67;
const double crLHS172 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS167;
const double crLHS173 = C(2,4)*DN(1,2);
const double crLHS174 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS173;
const double crLHS175 = DN(1,2)*crLHS161;
const double crLHS176 = DN(0,1)*N[1];
const double crLHS177 = DN(1,1)*N[0];
const double crLHS178 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS104;
const double crLHS179 = C(0,4)*DN(2,0) + crLHS107 + crLHS111;
const double crLHS180 = DN(2,0)*crLHS161;
const double crLHS181 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS182 = C(1,4)*DN(2,1);
const double crLHS183 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS182;
const double crLHS184 = DN(0,1)*DN(2,1);
const double crLHS185 = crLHS12*crLHS184;
const double crLHS186 = crLHS102 + crLHS90;
const double crLHS187 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS182;
const double crLHS188 = C(2,4)*DN(2,2);
const double crLHS189 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS188;
const double crLHS190 = DN(2,2)*crLHS161;
const double crLHS191 = DN(0,1)*N[2];
const double crLHS192 = DN(2,1)*N[0];
const double crLHS193 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS138;
const double crLHS194 = C(0,4)*DN(3,0) + crLHS141 + crLHS145;
const double crLHS195 = DN(3,0)*crLHS161;
const double crLHS196 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS197 = C(1,4)*DN(3,1);
const double crLHS198 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS197;
const double crLHS199 = DN(0,1)*DN(3,1);
const double crLHS200 = crLHS12*crLHS199;
const double crLHS201 = crLHS124 + crLHS136;
const double crLHS202 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS197;
const double crLHS203 = C(2,4)*DN(3,2);
const double crLHS204 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS203;
const double crLHS205 = DN(3,2)*crLHS161;
const double crLHS206 = DN(0,1)*N[3];
const double crLHS207 = DN(3,1)*N[0];
const double crLHS208 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS44;
const double crLHS209 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS159;
const double crLHS210 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS211 = pow(DN(0,2), 2);
const double crLHS212 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS78;
const double crLHS213 = DN(0,2)*crLHS12;
const double crLHS214 = DN(1,0)*crLHS213;
const double crLHS215 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS173;
const double crLHS216 = DN(1,1)*crLHS213;
const double crLHS217 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS218 = DN(0,2)*DN(1,2);
const double crLHS219 = crLHS12*crLHS218;
const double crLHS220 = DN(0,2)*N[1];
const double crLHS221 = DN(1,2)*N[0];
const double crLHS222 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS113;
const double crLHS223 = DN(2,0)*crLHS213;
const double crLHS224 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS188;
const double crLHS225 = DN(2,1)*crLHS213;
const double crLHS226 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS227 = DN(0,2)*DN(2,2);
const double crLHS228 = crLHS12*crLHS227;
const double crLHS229 = DN(0,2)*N[2];
const double crLHS230 = DN(2,2)*N[0];
const double crLHS231 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS147;
const double crLHS232 = DN(3,0)*crLHS213;
const double crLHS233 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS203;
const double crLHS234 = DN(3,1)*crLHS213;
const double crLHS235 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS236 = DN(0,2)*DN(3,2);
const double crLHS237 = crLHS12*crLHS236;
const double crLHS238 = DN(0,2)*N[3];
const double crLHS239 = DN(3,2)*N[0];
const double crLHS240 = gauss_weight*(N[0] + crLHS23*(1.0*crLHS19 + 1.0*crLHS20 + 1.0*crLHS21 + crLHS27));
const double crLHS241 = crLHS24*gauss_weight;
const double crLHS242 = crLHS241*(crLHS169 + crLHS218 + crLHS54);
const double crLHS243 = crLHS241*(crLHS184 + crLHS227 + crLHS89);
const double crLHS244 = crLHS241*(crLHS123 + crLHS199 + crLHS236);
const double crLHS245 = N[1]*rho;
const double crLHS246 = crLHS25*crLHS60;
const double crLHS247 = 1.0*crLHS61;
const double crLHS248 = crLHS23*crLHS247;
const double crLHS249 = crLHS245*crLHS30;
const double crLHS250 = crLHS16*crLHS245 + crLHS22*crLHS246 - crLHS22*crLHS248 + crLHS249*crLHS29;
const double crLHS251 = pow(DN(1,0), 2);
const double crLHS252 = pow(N[1], 2)*crLHS6 + crLHS245*crLHS60 + crLHS246*crLHS65 - crLHS248*crLHS65 + crLHS249*crLHS66;
const double crLHS253 = DN(1,0)*crLHS12;
const double crLHS254 = DN(1,1)*crLHS253;
const double crLHS255 = DN(1,2)*crLHS253;
const double crLHS256 = gauss_weight*(-N[1] + crLHS245*crLHS47 + crLHS246 - crLHS248);
const double crLHS257 = DN(1,0)*DN(2,0);
const double crLHS258 = N[2]*crLHS61;
const double crLHS259 = crLHS12*crLHS257 + crLHS258;
const double crLHS260 = crLHS100*crLHS246 - crLHS100*crLHS248 + crLHS101*crLHS249 + crLHS245*crLHS95;
const double crLHS261 = DN(2,1)*crLHS253;
const double crLHS262 = DN(2,2)*crLHS253;
const double crLHS263 = DN(1,0)*N[2];
const double crLHS264 = DN(2,0)*N[1];
const double crLHS265 = DN(1,0)*DN(3,0);
const double crLHS266 = N[3]*crLHS61;
const double crLHS267 = crLHS12*crLHS265 + crLHS266;
const double crLHS268 = crLHS129*crLHS245 + crLHS134*crLHS246 - crLHS134*crLHS248 + crLHS135*crLHS249;
const double crLHS269 = DN(3,1)*crLHS253;
const double crLHS270 = DN(3,2)*crLHS253;
const double crLHS271 = DN(1,0)*N[3];
const double crLHS272 = DN(3,0)*N[1];
const double crLHS273 = crLHS250 + crLHS55;
const double crLHS274 = pow(DN(1,1), 2);
const double crLHS275 = DN(1,1)*crLHS12;
const double crLHS276 = DN(1,2)*crLHS275;
const double crLHS277 = DN(2,0)*crLHS275;
const double crLHS278 = DN(1,1)*DN(2,1);
const double crLHS279 = crLHS12*crLHS278;
const double crLHS280 = crLHS258 + crLHS260;
const double crLHS281 = DN(2,2)*crLHS275;
const double crLHS282 = DN(1,1)*N[2];
const double crLHS283 = DN(2,1)*N[1];
const double crLHS284 = DN(3,0)*crLHS275;
const double crLHS285 = DN(1,1)*DN(3,1);
const double crLHS286 = crLHS12*crLHS285;
const double crLHS287 = crLHS266 + crLHS268;
const double crLHS288 = DN(3,2)*crLHS275;
const double crLHS289 = DN(1,1)*N[3];
const double crLHS290 = DN(3,1)*N[1];
const double crLHS291 = pow(DN(1,2), 2);
const double crLHS292 = DN(1,2)*crLHS12;
const double crLHS293 = DN(2,0)*crLHS292;
const double crLHS294 = DN(2,1)*crLHS292;
const double crLHS295 = DN(1,2)*DN(2,2);
const double crLHS296 = crLHS12*crLHS295;
const double crLHS297 = DN(1,2)*N[2];
const double crLHS298 = DN(2,2)*N[1];
const double crLHS299 = DN(3,0)*crLHS292;
const double crLHS300 = DN(3,1)*crLHS292;
const double crLHS301 = DN(1,2)*DN(3,2);
const double crLHS302 = crLHS12*crLHS301;
const double crLHS303 = DN(1,2)*N[3];
const double crLHS304 = DN(3,2)*N[1];
const double crLHS305 = gauss_weight*(N[1] + crLHS23*(crLHS247 + 1.0*crLHS62 + 1.0*crLHS63 + 1.0*crLHS64));
const double crLHS306 = crLHS241*(crLHS257 + crLHS278 + crLHS295);
const double crLHS307 = crLHS241*(crLHS265 + crLHS285 + crLHS301);
const double crLHS308 = N[2]*rho;
const double crLHS309 = crLHS25*crLHS95;
const double crLHS310 = 1.0*crLHS96;
const double crLHS311 = crLHS23*crLHS310;
const double crLHS312 = crLHS30*crLHS308;
const double crLHS313 = crLHS16*crLHS308 + crLHS22*crLHS309 - crLHS22*crLHS311 + crLHS29*crLHS312;
const double crLHS314 = crLHS308*crLHS60 + crLHS309*crLHS65 - crLHS311*crLHS65 + crLHS312*crLHS66;
const double crLHS315 = pow(DN(2,0), 2);
const double crLHS316 = pow(N[2], 2)*crLHS6 + crLHS100*crLHS309 - crLHS100*crLHS311 + crLHS101*crLHS312 + crLHS308*crLHS95;
const double crLHS317 = DN(2,0)*crLHS12;
const double crLHS318 = DN(2,1)*crLHS317;
const double crLHS319 = DN(2,2)*crLHS317;
const double crLHS320 = gauss_weight*(-N[2] + crLHS308*crLHS47 + crLHS309 - crLHS311);
const double crLHS321 = DN(2,0)*DN(3,0);
const double crLHS322 = N[3]*crLHS96;
const double crLHS323 = crLHS12*crLHS321 + crLHS322;
const double crLHS324 = crLHS129*crLHS308 + crLHS134*crLHS309 - crLHS134*crLHS311 + crLHS135*crLHS312;
const double crLHS325 = DN(3,1)*crLHS317;
const double crLHS326 = DN(3,2)*crLHS317;
const double crLHS327 = DN(2,0)*N[3];
const double crLHS328 = DN(3,0)*N[2];
const double crLHS329 = crLHS313 + crLHS90;
const double crLHS330 = crLHS258 + crLHS314;
const double crLHS331 = pow(DN(2,1), 2);
const double crLHS332 = DN(2,1)*crLHS12;
const double crLHS333 = DN(2,2)*crLHS332;
const double crLHS334 = DN(3,0)*crLHS332;
const double crLHS335 = DN(2,1)*DN(3,1);
const double crLHS336 = crLHS12*crLHS335;
const double crLHS337 = crLHS322 + crLHS324;
const double crLHS338 = DN(3,2)*crLHS332;
const double crLHS339 = DN(2,1)*N[3];
const double crLHS340 = DN(3,1)*N[2];
const double crLHS341 = pow(DN(2,2), 2);
const double crLHS342 = DN(2,2)*crLHS12;
const double crLHS343 = DN(3,0)*crLHS342;
const double crLHS344 = DN(3,1)*crLHS342;
const double crLHS345 = DN(2,2)*DN(3,2);
const double crLHS346 = crLHS12*crLHS345;
const double crLHS347 = DN(2,2)*N[3];
const double crLHS348 = DN(3,2)*N[2];
const double crLHS349 = gauss_weight*(N[2] + crLHS23*(crLHS310 + 1.0*crLHS97 + 1.0*crLHS98 + 1.0*crLHS99));
const double crLHS350 = crLHS241*(crLHS321 + crLHS335 + crLHS345);
const double crLHS351 = N[3]*rho;
const double crLHS352 = crLHS129*crLHS25;
const double crLHS353 = 1.0*crLHS130;
const double crLHS354 = crLHS23*crLHS353;
const double crLHS355 = crLHS30*crLHS351;
const double crLHS356 = crLHS16*crLHS351 + crLHS22*crLHS352 - crLHS22*crLHS354 + crLHS29*crLHS355;
const double crLHS357 = crLHS351*crLHS60 + crLHS352*crLHS65 - crLHS354*crLHS65 + crLHS355*crLHS66;
const double crLHS358 = crLHS100*crLHS352 - crLHS100*crLHS354 + crLHS101*crLHS355 + crLHS351*crLHS95;
const double crLHS359 = pow(DN(3,0), 2);
const double crLHS360 = pow(N[3], 2)*crLHS6 + crLHS129*crLHS351 + crLHS134*crLHS352 - crLHS134*crLHS354 + crLHS135*crLHS355;
const double crLHS361 = DN(3,0)*crLHS12;
const double crLHS362 = DN(3,1)*crLHS361;
const double crLHS363 = DN(3,2)*crLHS361;
const double crLHS364 = gauss_weight*(-N[3] + crLHS351*crLHS47 + crLHS352 - crLHS354);
const double crLHS365 = crLHS124 + crLHS356;
const double crLHS366 = crLHS266 + crLHS357;
const double crLHS367 = crLHS322 + crLHS358;
const double crLHS368 = pow(DN(3,1), 2);
const double crLHS369 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS370 = pow(DN(3,2), 2);
const double crLHS371 = gauss_weight*(N[3] + crLHS23*(1.0*crLHS131 + 1.0*crLHS132 + 1.0*crLHS133 + crLHS353));
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS32);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS33 + DN(0,1)*crLHS35 + DN(0,2)*crLHS38 + crLHS40);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS41 + DN(0,1)*crLHS43 + DN(0,2)*crLHS45 + crLHS46);
rLHS(0,3)+=DN(0,0)*crLHS48;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS49 + DN(0,1)*crLHS51 + DN(0,2)*crLHS53 + crLHS56 + crLHS67);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS68 + DN(0,1)*crLHS70 + DN(0,2)*crLHS73 + crLHS74);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS77 + DN(0,2)*crLHS79 + crLHS80);
rLHS(0,7)+=-gauss_weight*(-DN(1,0)*crLHS26 + DN(1,0)*crLHS28 + crLHS81 - crLHS82*crLHS83);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS84 + DN(0,1)*crLHS86 + DN(0,2)*crLHS88 + crLHS102 + crLHS91);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS103 + DN(0,1)*crLHS105 + DN(0,2)*crLHS108 + crLHS109);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS112 + DN(0,2)*crLHS114 + crLHS115);
rLHS(0,11)+=-gauss_weight*(-DN(2,0)*crLHS26 + DN(2,0)*crLHS28 + crLHS116 - crLHS117*crLHS83);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS118 + DN(0,1)*crLHS120 + DN(0,2)*crLHS122 + crLHS125 + crLHS136);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS137 + DN(0,1)*crLHS139 + DN(0,2)*crLHS142 + crLHS143);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS144 + DN(0,1)*crLHS146 + DN(0,2)*crLHS148 + crLHS149);
rLHS(0,15)+=-gauss_weight*(-DN(3,0)*crLHS26 + DN(3,0)*crLHS28 + crLHS150 - crLHS151*crLHS83);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS152 + DN(0,2)*crLHS153 + crLHS40);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS35 + DN(0,1)*crLHS154 + DN(0,2)*crLHS156 + crLHS12*crLHS157 + crLHS32);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS158 + DN(0,2)*crLHS160 + crLHS162);
rLHS(1,3)+=DN(0,1)*crLHS48;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS51 + DN(0,1)*crLHS163 + DN(0,2)*crLHS164 + crLHS165);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS70 + DN(0,1)*crLHS166 + DN(0,2)*crLHS168 + crLHS170 + crLHS171);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS172 + DN(0,2)*crLHS174 + crLHS175);
rLHS(1,7)+=-gauss_weight*(-DN(1,1)*crLHS26 + DN(1,1)*crLHS28 + crLHS176 - crLHS177*crLHS83);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS86 + DN(0,1)*crLHS178 + DN(0,2)*crLHS179 + crLHS180);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS105 + DN(0,1)*crLHS181 + DN(0,2)*crLHS183 + crLHS185 + crLHS186);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS187 + DN(0,2)*crLHS189 + crLHS190);
rLHS(1,11)+=-gauss_weight*(-DN(2,1)*crLHS26 + DN(2,1)*crLHS28 + crLHS191 - crLHS192*crLHS83);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS120 + DN(0,1)*crLHS193 + DN(0,2)*crLHS194 + crLHS195);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS139 + DN(0,1)*crLHS196 + DN(0,2)*crLHS198 + crLHS200 + crLHS201);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS146 + DN(0,1)*crLHS202 + DN(0,2)*crLHS204 + crLHS205);
rLHS(1,15)+=-gauss_weight*(-DN(3,1)*crLHS26 + DN(3,1)*crLHS28 + crLHS206 - crLHS207*crLHS83);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS153 + DN(0,2)*crLHS208 + crLHS46);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS156 + DN(0,2)*crLHS209 + crLHS162);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS160 + DN(0,2)*crLHS210 + crLHS12*crLHS211 + crLHS32);
rLHS(2,3)+=DN(0,2)*crLHS48;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS53 + DN(0,1)*crLHS164 + DN(0,2)*crLHS212 + crLHS214);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS73 + DN(0,1)*crLHS168 + DN(0,2)*crLHS215 + crLHS216);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS174 + DN(0,2)*crLHS217 + crLHS171 + crLHS219);
rLHS(2,7)+=-gauss_weight*(-DN(1,2)*crLHS26 + DN(1,2)*crLHS28 + crLHS220 - crLHS221*crLHS83);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS88 + DN(0,1)*crLHS179 + DN(0,2)*crLHS222 + crLHS223);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS108 + DN(0,1)*crLHS183 + DN(0,2)*crLHS224 + crLHS225);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS189 + DN(0,2)*crLHS226 + crLHS186 + crLHS228);
rLHS(2,11)+=-gauss_weight*(-DN(2,2)*crLHS26 + DN(2,2)*crLHS28 + crLHS229 - crLHS230*crLHS83);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS122 + DN(0,1)*crLHS194 + DN(0,2)*crLHS231 + crLHS232);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS142 + DN(0,1)*crLHS198 + DN(0,2)*crLHS233 + crLHS234);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS148 + DN(0,1)*crLHS204 + DN(0,2)*crLHS235 + crLHS201 + crLHS237);
rLHS(2,15)+=-gauss_weight*(-DN(3,2)*crLHS26 + DN(3,2)*crLHS28 + crLHS238 - crLHS239*crLHS83);
rLHS(3,0)+=DN(0,0)*crLHS240;
rLHS(3,1)+=DN(0,1)*crLHS240;
rLHS(3,2)+=DN(0,2)*crLHS240;
rLHS(3,3)+=crLHS241*(crLHS157 + crLHS211 + crLHS5);
rLHS(3,4)+=gauss_weight*(DN(0,0)*crLHS66 + crLHS82);
rLHS(3,5)+=gauss_weight*(DN(0,1)*crLHS66 + crLHS177);
rLHS(3,6)+=gauss_weight*(DN(0,2)*crLHS66 + crLHS221);
rLHS(3,7)+=crLHS242;
rLHS(3,8)+=gauss_weight*(DN(0,0)*crLHS101 + crLHS117);
rLHS(3,9)+=gauss_weight*(DN(0,1)*crLHS101 + crLHS192);
rLHS(3,10)+=gauss_weight*(DN(0,2)*crLHS101 + crLHS230);
rLHS(3,11)+=crLHS243;
rLHS(3,12)+=gauss_weight*(DN(0,0)*crLHS135 + crLHS151);
rLHS(3,13)+=gauss_weight*(DN(0,1)*crLHS135 + crLHS207);
rLHS(3,14)+=gauss_weight*(DN(0,2)*crLHS135 + crLHS239);
rLHS(3,15)+=crLHS244;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS250 + crLHS56);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS33 + DN(1,1)*crLHS35 + DN(1,2)*crLHS38 + crLHS165);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS41 + DN(1,1)*crLHS43 + DN(1,2)*crLHS45 + crLHS214);
rLHS(4,3)+=gauss_weight*(DN(0,0)*crLHS246 - DN(0,0)*crLHS248 + crLHS81*crLHS83 - crLHS82);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS49 + DN(1,1)*crLHS51 + DN(1,2)*crLHS53 + crLHS12*crLHS251 + crLHS252);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS68 + DN(1,1)*crLHS70 + DN(1,2)*crLHS73 + crLHS254);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS77 + DN(1,2)*crLHS79 + crLHS255);
rLHS(4,7)+=DN(1,0)*crLHS256;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS84 + DN(1,1)*crLHS86 + DN(1,2)*crLHS88 + crLHS259 + crLHS260);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS103 + DN(1,1)*crLHS105 + DN(1,2)*crLHS108 + crLHS261);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS112 + DN(1,2)*crLHS114 + crLHS262);
rLHS(4,11)+=-gauss_weight*(-DN(2,0)*crLHS246 + DN(2,0)*crLHS248 + crLHS263 - crLHS264*crLHS83);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS118 + DN(1,1)*crLHS120 + DN(1,2)*crLHS122 + crLHS267 + crLHS268);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS137 + DN(1,1)*crLHS139 + DN(1,2)*crLHS142 + crLHS269);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS144 + DN(1,1)*crLHS146 + DN(1,2)*crLHS148 + crLHS270);
rLHS(4,15)+=-gauss_weight*(-DN(3,0)*crLHS246 + DN(3,0)*crLHS248 + crLHS271 - crLHS272*crLHS83);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS152 + DN(1,2)*crLHS153 + crLHS74);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS35 + DN(1,1)*crLHS154 + DN(1,2)*crLHS156 + crLHS170 + crLHS273);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS158 + DN(1,2)*crLHS160 + crLHS216);
rLHS(5,3)+=gauss_weight*(DN(0,1)*crLHS246 - DN(0,1)*crLHS248 + crLHS176*crLHS83 - crLHS177);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS51 + DN(1,1)*crLHS163 + DN(1,2)*crLHS164 + crLHS254);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS70 + DN(1,1)*crLHS166 + DN(1,2)*crLHS168 + crLHS12*crLHS274 + crLHS252);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS172 + DN(1,2)*crLHS174 + crLHS276);
rLHS(5,7)+=DN(1,1)*crLHS256;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS86 + DN(1,1)*crLHS178 + DN(1,2)*crLHS179 + crLHS277);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS105 + DN(1,1)*crLHS181 + DN(1,2)*crLHS183 + crLHS279 + crLHS280);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS187 + DN(1,2)*crLHS189 + crLHS281);
rLHS(5,11)+=-gauss_weight*(-DN(2,1)*crLHS246 + DN(2,1)*crLHS248 + crLHS282 - crLHS283*crLHS83);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS120 + DN(1,1)*crLHS193 + DN(1,2)*crLHS194 + crLHS284);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS139 + DN(1,1)*crLHS196 + DN(1,2)*crLHS198 + crLHS286 + crLHS287);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS146 + DN(1,1)*crLHS202 + DN(1,2)*crLHS204 + crLHS288);
rLHS(5,15)+=-gauss_weight*(-DN(3,1)*crLHS246 + DN(3,1)*crLHS248 + crLHS289 - crLHS290*crLHS83);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS153 + DN(1,2)*crLHS208 + crLHS80);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS156 + DN(1,2)*crLHS209 + crLHS175);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS160 + DN(1,2)*crLHS210 + crLHS219 + crLHS273);
rLHS(6,3)+=gauss_weight*(DN(0,2)*crLHS246 - DN(0,2)*crLHS248 + crLHS220*crLHS83 - crLHS221);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS53 + DN(1,1)*crLHS164 + DN(1,2)*crLHS212 + crLHS255);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS73 + DN(1,1)*crLHS168 + DN(1,2)*crLHS215 + crLHS276);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS174 + DN(1,2)*crLHS217 + crLHS12*crLHS291 + crLHS252);
rLHS(6,7)+=DN(1,2)*crLHS256;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS88 + DN(1,1)*crLHS179 + DN(1,2)*crLHS222 + crLHS293);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS108 + DN(1,1)*crLHS183 + DN(1,2)*crLHS224 + crLHS294);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS189 + DN(1,2)*crLHS226 + crLHS280 + crLHS296);
rLHS(6,11)+=-gauss_weight*(-DN(2,2)*crLHS246 + DN(2,2)*crLHS248 + crLHS297 - crLHS298*crLHS83);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS122 + DN(1,1)*crLHS194 + DN(1,2)*crLHS231 + crLHS299);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS142 + DN(1,1)*crLHS198 + DN(1,2)*crLHS233 + crLHS300);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS148 + DN(1,1)*crLHS204 + DN(1,2)*crLHS235 + crLHS287 + crLHS302);
rLHS(6,15)+=-gauss_weight*(-DN(3,2)*crLHS246 + DN(3,2)*crLHS248 + crLHS303 - crLHS304*crLHS83);
rLHS(7,0)+=gauss_weight*(DN(1,0)*crLHS29 + crLHS81);
rLHS(7,1)+=gauss_weight*(DN(1,1)*crLHS29 + crLHS176);
rLHS(7,2)+=gauss_weight*(DN(1,2)*crLHS29 + crLHS220);
rLHS(7,3)+=crLHS242;
rLHS(7,4)+=DN(1,0)*crLHS305;
rLHS(7,5)+=DN(1,1)*crLHS305;
rLHS(7,6)+=DN(1,2)*crLHS305;
rLHS(7,7)+=crLHS241*(crLHS251 + crLHS274 + crLHS291);
rLHS(7,8)+=gauss_weight*(DN(1,0)*crLHS101 + crLHS264);
rLHS(7,9)+=gauss_weight*(DN(1,1)*crLHS101 + crLHS283);
rLHS(7,10)+=gauss_weight*(DN(1,2)*crLHS101 + crLHS298);
rLHS(7,11)+=crLHS306;
rLHS(7,12)+=gauss_weight*(DN(1,0)*crLHS135 + crLHS272);
rLHS(7,13)+=gauss_weight*(DN(1,1)*crLHS135 + crLHS290);
rLHS(7,14)+=gauss_weight*(DN(1,2)*crLHS135 + crLHS304);
rLHS(7,15)+=crLHS307;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS313 + crLHS91);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS33 + DN(2,1)*crLHS35 + DN(2,2)*crLHS38 + crLHS180);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS41 + DN(2,1)*crLHS43 + DN(2,2)*crLHS45 + crLHS223);
rLHS(8,3)+=gauss_weight*(DN(0,0)*crLHS309 - DN(0,0)*crLHS311 + crLHS116*crLHS83 - crLHS117);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS49 + DN(2,1)*crLHS51 + DN(2,2)*crLHS53 + crLHS259 + crLHS314);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS68 + DN(2,1)*crLHS70 + DN(2,2)*crLHS73 + crLHS277);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS77 + DN(2,2)*crLHS79 + crLHS293);
rLHS(8,7)+=gauss_weight*(DN(1,0)*crLHS309 - DN(1,0)*crLHS311 + crLHS263*crLHS83 - crLHS264);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS84 + DN(2,1)*crLHS86 + DN(2,2)*crLHS88 + crLHS12*crLHS315 + crLHS316);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS103 + DN(2,1)*crLHS105 + DN(2,2)*crLHS108 + crLHS318);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS112 + DN(2,2)*crLHS114 + crLHS319);
rLHS(8,11)+=DN(2,0)*crLHS320;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS118 + DN(2,1)*crLHS120 + DN(2,2)*crLHS122 + crLHS323 + crLHS324);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS137 + DN(2,1)*crLHS139 + DN(2,2)*crLHS142 + crLHS325);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS144 + DN(2,1)*crLHS146 + DN(2,2)*crLHS148 + crLHS326);
rLHS(8,15)+=-gauss_weight*(-DN(3,0)*crLHS309 + DN(3,0)*crLHS311 + crLHS327 - crLHS328*crLHS83);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS152 + DN(2,2)*crLHS153 + crLHS109);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS35 + DN(2,1)*crLHS154 + DN(2,2)*crLHS156 + crLHS185 + crLHS329);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS158 + DN(2,2)*crLHS160 + crLHS225);
rLHS(9,3)+=gauss_weight*(DN(0,1)*crLHS309 - DN(0,1)*crLHS311 + crLHS191*crLHS83 - crLHS192);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS51 + DN(2,1)*crLHS163 + DN(2,2)*crLHS164 + crLHS261);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS70 + DN(2,1)*crLHS166 + DN(2,2)*crLHS168 + crLHS279 + crLHS330);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS172 + DN(2,2)*crLHS174 + crLHS294);
rLHS(9,7)+=gauss_weight*(DN(1,1)*crLHS309 - DN(1,1)*crLHS311 + crLHS282*crLHS83 - crLHS283);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS86 + DN(2,1)*crLHS178 + DN(2,2)*crLHS179 + crLHS318);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS105 + DN(2,1)*crLHS181 + DN(2,2)*crLHS183 + crLHS12*crLHS331 + crLHS316);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS187 + DN(2,2)*crLHS189 + crLHS333);
rLHS(9,11)+=DN(2,1)*crLHS320;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS120 + DN(2,1)*crLHS193 + DN(2,2)*crLHS194 + crLHS334);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS139 + DN(2,1)*crLHS196 + DN(2,2)*crLHS198 + crLHS336 + crLHS337);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS146 + DN(2,1)*crLHS202 + DN(2,2)*crLHS204 + crLHS338);
rLHS(9,15)+=-gauss_weight*(-DN(3,1)*crLHS309 + DN(3,1)*crLHS311 + crLHS339 - crLHS340*crLHS83);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS153 + DN(2,2)*crLHS208 + crLHS115);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS156 + DN(2,2)*crLHS209 + crLHS190);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS160 + DN(2,2)*crLHS210 + crLHS228 + crLHS329);
rLHS(10,3)+=gauss_weight*(DN(0,2)*crLHS309 - DN(0,2)*crLHS311 + crLHS229*crLHS83 - crLHS230);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS53 + DN(2,1)*crLHS164 + DN(2,2)*crLHS212 + crLHS262);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS73 + DN(2,1)*crLHS168 + DN(2,2)*crLHS215 + crLHS281);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS174 + DN(2,2)*crLHS217 + crLHS296 + crLHS330);
rLHS(10,7)+=gauss_weight*(DN(1,2)*crLHS309 - DN(1,2)*crLHS311 + crLHS297*crLHS83 - crLHS298);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS88 + DN(2,1)*crLHS179 + DN(2,2)*crLHS222 + crLHS319);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS108 + DN(2,1)*crLHS183 + DN(2,2)*crLHS224 + crLHS333);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS189 + DN(2,2)*crLHS226 + crLHS12*crLHS341 + crLHS316);
rLHS(10,11)+=DN(2,2)*crLHS320;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS122 + DN(2,1)*crLHS194 + DN(2,2)*crLHS231 + crLHS343);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS142 + DN(2,1)*crLHS198 + DN(2,2)*crLHS233 + crLHS344);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS148 + DN(2,1)*crLHS204 + DN(2,2)*crLHS235 + crLHS337 + crLHS346);
rLHS(10,15)+=-gauss_weight*(-DN(3,2)*crLHS309 + DN(3,2)*crLHS311 + crLHS347 - crLHS348*crLHS83);
rLHS(11,0)+=gauss_weight*(DN(2,0)*crLHS29 + crLHS116);
rLHS(11,1)+=gauss_weight*(DN(2,1)*crLHS29 + crLHS191);
rLHS(11,2)+=gauss_weight*(DN(2,2)*crLHS29 + crLHS229);
rLHS(11,3)+=crLHS243;
rLHS(11,4)+=gauss_weight*(DN(2,0)*crLHS66 + crLHS263);
rLHS(11,5)+=gauss_weight*(DN(2,1)*crLHS66 + crLHS282);
rLHS(11,6)+=gauss_weight*(DN(2,2)*crLHS66 + crLHS297);
rLHS(11,7)+=crLHS306;
rLHS(11,8)+=DN(2,0)*crLHS349;
rLHS(11,9)+=DN(2,1)*crLHS349;
rLHS(11,10)+=DN(2,2)*crLHS349;
rLHS(11,11)+=crLHS241*(crLHS315 + crLHS331 + crLHS341);
rLHS(11,12)+=gauss_weight*(DN(2,0)*crLHS135 + crLHS328);
rLHS(11,13)+=gauss_weight*(DN(2,1)*crLHS135 + crLHS340);
rLHS(11,14)+=gauss_weight*(DN(2,2)*crLHS135 + crLHS348);
rLHS(11,15)+=crLHS350;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS125 + crLHS356);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS33 + DN(3,1)*crLHS35 + DN(3,2)*crLHS38 + crLHS195);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS41 + DN(3,1)*crLHS43 + DN(3,2)*crLHS45 + crLHS232);
rLHS(12,3)+=gauss_weight*(DN(0,0)*crLHS352 - DN(0,0)*crLHS354 + crLHS150*crLHS83 - crLHS151);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS49 + DN(3,1)*crLHS51 + DN(3,2)*crLHS53 + crLHS267 + crLHS357);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS68 + DN(3,1)*crLHS70 + DN(3,2)*crLHS73 + crLHS284);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS75 + DN(3,1)*crLHS77 + DN(3,2)*crLHS79 + crLHS299);
rLHS(12,7)+=gauss_weight*(DN(1,0)*crLHS352 - DN(1,0)*crLHS354 + crLHS271*crLHS83 - crLHS272);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS84 + DN(3,1)*crLHS86 + DN(3,2)*crLHS88 + crLHS323 + crLHS358);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS103 + DN(3,1)*crLHS105 + DN(3,2)*crLHS108 + crLHS334);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS110 + DN(3,1)*crLHS112 + DN(3,2)*crLHS114 + crLHS343);
rLHS(12,11)+=gauss_weight*(DN(2,0)*crLHS352 - DN(2,0)*crLHS354 + crLHS327*crLHS83 - crLHS328);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS118 + DN(3,1)*crLHS120 + DN(3,2)*crLHS122 + crLHS12*crLHS359 + crLHS360);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS137 + DN(3,1)*crLHS139 + DN(3,2)*crLHS142 + crLHS362);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS144 + DN(3,1)*crLHS146 + DN(3,2)*crLHS148 + crLHS363);
rLHS(12,15)+=DN(3,0)*crLHS364;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS152 + DN(3,2)*crLHS153 + crLHS143);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS35 + DN(3,1)*crLHS154 + DN(3,2)*crLHS156 + crLHS200 + crLHS365);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS43 + DN(3,1)*crLHS158 + DN(3,2)*crLHS160 + crLHS234);
rLHS(13,3)+=gauss_weight*(DN(0,1)*crLHS352 - DN(0,1)*crLHS354 + crLHS206*crLHS83 - crLHS207);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS51 + DN(3,1)*crLHS163 + DN(3,2)*crLHS164 + crLHS269);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS70 + DN(3,1)*crLHS166 + DN(3,2)*crLHS168 + crLHS286 + crLHS366);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS172 + DN(3,2)*crLHS174 + crLHS300);
rLHS(13,7)+=gauss_weight*(DN(1,1)*crLHS352 - DN(1,1)*crLHS354 + crLHS289*crLHS83 - crLHS290);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS86 + DN(3,1)*crLHS178 + DN(3,2)*crLHS179 + crLHS325);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS105 + DN(3,1)*crLHS181 + DN(3,2)*crLHS183 + crLHS336 + crLHS367);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS112 + DN(3,1)*crLHS187 + DN(3,2)*crLHS189 + crLHS344);
rLHS(13,11)+=gauss_weight*(DN(2,1)*crLHS352 - DN(2,1)*crLHS354 + crLHS339*crLHS83 - crLHS340);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS120 + DN(3,1)*crLHS193 + DN(3,2)*crLHS194 + crLHS362);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS139 + DN(3,1)*crLHS196 + DN(3,2)*crLHS198 + crLHS12*crLHS368 + crLHS360);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS146 + DN(3,1)*crLHS202 + DN(3,2)*crLHS204 + crLHS369);
rLHS(13,15)+=DN(3,1)*crLHS364;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS153 + DN(3,2)*crLHS208 + crLHS149);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS38 + DN(3,1)*crLHS156 + DN(3,2)*crLHS209 + crLHS205);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS45 + DN(3,1)*crLHS160 + DN(3,2)*crLHS210 + crLHS237 + crLHS365);
rLHS(14,3)+=gauss_weight*(DN(0,2)*crLHS352 - DN(0,2)*crLHS354 + crLHS238*crLHS83 - crLHS239);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS53 + DN(3,1)*crLHS164 + DN(3,2)*crLHS212 + crLHS270);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS73 + DN(3,1)*crLHS168 + DN(3,2)*crLHS215 + crLHS288);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS79 + DN(3,1)*crLHS174 + DN(3,2)*crLHS217 + crLHS302 + crLHS366);
rLHS(14,7)+=gauss_weight*(DN(1,2)*crLHS352 - DN(1,2)*crLHS354 + crLHS303*crLHS83 - crLHS304);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS88 + DN(3,1)*crLHS179 + DN(3,2)*crLHS222 + crLHS326);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS108 + DN(3,1)*crLHS183 + DN(3,2)*crLHS224 + crLHS338);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS114 + DN(3,1)*crLHS189 + DN(3,2)*crLHS226 + crLHS346 + crLHS367);
rLHS(14,11)+=gauss_weight*(DN(2,2)*crLHS352 - DN(2,2)*crLHS354 + crLHS347*crLHS83 - crLHS348);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS122 + DN(3,1)*crLHS194 + DN(3,2)*crLHS231 + crLHS363);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS142 + DN(3,1)*crLHS198 + DN(3,2)*crLHS233 + crLHS369);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS148 + DN(3,1)*crLHS204 + DN(3,2)*crLHS235 + crLHS12*crLHS370 + crLHS360);
rLHS(14,15)+=DN(3,2)*crLHS364;
rLHS(15,0)+=gauss_weight*(DN(3,0)*crLHS29 + crLHS150);
rLHS(15,1)+=gauss_weight*(DN(3,1)*crLHS29 + crLHS206);
rLHS(15,2)+=gauss_weight*(DN(3,2)*crLHS29 + crLHS238);
rLHS(15,3)+=crLHS244;
rLHS(15,4)+=gauss_weight*(DN(3,0)*crLHS66 + crLHS271);
rLHS(15,5)+=gauss_weight*(DN(3,1)*crLHS66 + crLHS289);
rLHS(15,6)+=gauss_weight*(DN(3,2)*crLHS66 + crLHS303);
rLHS(15,7)+=crLHS307;
rLHS(15,8)+=gauss_weight*(DN(3,0)*crLHS101 + crLHS327);
rLHS(15,9)+=gauss_weight*(DN(3,1)*crLHS101 + crLHS339);
rLHS(15,10)+=gauss_weight*(DN(3,2)*crLHS101 + crLHS347);
rLHS(15,11)+=crLHS350;
rLHS(15,12)+=DN(3,0)*crLHS371;
rLHS(15,13)+=DN(3,1)*crLHS371;
rLHS(15,14)+=DN(3,2)*crLHS371;
rLHS(15,15)+=crLHS241*(crLHS359 + crLHS368 + crLHS370);

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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    // const double bdf1 = rData.bdf1;
    // const double bdf2 = rData.bdf2;

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
const double crRHS4 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crRHS5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crRHS6 = crRHS4*crRHS5;
const double crRHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crRHS8 = crRHS7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
const double crRHS9 = crRHS6 + crRHS8;
const double crRHS10 = N[0]*rho;
const double crRHS11 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crRHS12 = crRHS11 + crRHS4;
const double crRHS13 = crRHS2*stab_c3;
const double crRHS14 = rho*stab_c2*sqrt(pow(crRHS5, 2) + pow(crRHS7, 2));
const double crRHS15 = crRHS12*(h*(crRHS13*h + crRHS14)/stab_c1 + mu);
const double crRHS16 = 1.0/(crRHS13 + crRHS14/h + mu*stab_c1/pow(h, 2));
const double crRHS17 = crRHS16*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crRHS1 + crRHS3 + crRHS6*rho + crRHS8*rho);
const double crRHS18 = N[0]*crRHS2;
const double crRHS19 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crRHS20 = crRHS10*crRHS19;
const double crRHS21 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS7);
const double crRHS22 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crRHS23 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
const double crRHS24 = crRHS5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
const double crRHS25 = crRHS11*crRHS7;
const double crRHS26 = crRHS24 + crRHS25;
const double crRHS27 = crRHS16*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crRHS22 + crRHS23 + crRHS24*rho + crRHS25*rho);
const double crRHS28 = N[1]*rho;
const double crRHS29 = N[1]*crRHS2;
const double crRHS30 = crRHS19*crRHS28;
const double crRHS31 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS7);
const double crRHS32 = N[2]*rho;
const double crRHS33 = N[2]*crRHS2;
const double crRHS34 = crRHS19*crRHS32;
const double crRHS35 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS7);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS15 + DN(0,0)*stress[0] + DN(0,1)*stress[2] - N[0]*crRHS1 + N[0]*crRHS3 + crRHS10*crRHS9 - crRHS17*crRHS18 + crRHS17*crRHS20 + crRHS17*crRHS21);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS15 + DN(0,1)*stress[1] - N[0]*crRHS22 + N[0]*crRHS23 + crRHS10*crRHS26 - crRHS18*crRHS27 + crRHS20*crRHS27 + crRHS21*crRHS27);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS17 + DN(0,1)*crRHS27 + N[0]*crRHS12);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS15 + DN(1,0)*stress[0] + DN(1,1)*stress[2] - N[1]*crRHS1 + N[1]*crRHS3 - crRHS17*crRHS29 + crRHS17*crRHS30 + crRHS17*crRHS31 + crRHS28*crRHS9);
rRHS[4]+=-gauss_weight*(DN(1,0)*stress[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS15 + DN(1,1)*stress[1] - N[1]*crRHS22 + N[1]*crRHS23 + crRHS26*crRHS28 - crRHS27*crRHS29 + crRHS27*crRHS30 + crRHS27*crRHS31);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS17 + DN(1,1)*crRHS27 + N[1]*crRHS12);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS15 + DN(2,0)*stress[0] + DN(2,1)*stress[2] - N[2]*crRHS1 + N[2]*crRHS3 - crRHS17*crRHS33 + crRHS17*crRHS34 + crRHS17*crRHS35 + crRHS32*crRHS9);
rRHS[7]+=-gauss_weight*(DN(2,0)*stress[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS15 + DN(2,1)*stress[1] - N[2]*crRHS22 + N[2]*crRHS23 + crRHS26*crRHS32 - crRHS27*crRHS33 + crRHS27*crRHS34 + crRHS27*crRHS35);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS17 + DN(2,1)*crRHS27 + N[2]*crRHS12);

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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    // const double bdf1 = rData.bdf1;
    // const double bdf2 = rData.bdf2;

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
    const array_1d<double,4>& stress = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // NAVIER-STOKES ELEMENTAL RHS VECTOR
    const double crRHS0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crRHS1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crRHS2 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crRHS3 = crRHS2*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
const double crRHS4 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crRHS5 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crRHS6 = crRHS4*crRHS5;
const double crRHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crRHS8 = crRHS7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0));
const double crRHS9 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crRHS10 = crRHS9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
const double crRHS11 = crRHS10 + crRHS6 + crRHS8;
const double crRHS12 = N[0]*rho;
const double crRHS13 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crRHS14 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crRHS15 = crRHS13 + crRHS14 + crRHS4;
const double crRHS16 = crRHS2*stab_c3;
const double crRHS17 = rho*stab_c2*sqrt(pow(crRHS5, 2) + pow(crRHS7, 2) + pow(crRHS9, 2));
const double crRHS18 = crRHS15*(h*(crRHS16*h + crRHS17)/stab_c1 + mu);
const double crRHS19 = 1.0/(crRHS16 + crRHS17/h + mu*stab_c1/pow(h, 2));
const double crRHS20 = crRHS19*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crRHS1 + crRHS10*rho + crRHS3 + crRHS6*rho + crRHS8*rho);
const double crRHS21 = N[0]*crRHS2;
const double crRHS22 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crRHS23 = crRHS12*crRHS22;
const double crRHS24 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS7 + DN(0,2)*crRHS9);
const double crRHS25 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crRHS26 = crRHS2*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
const double crRHS27 = crRHS5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1));
const double crRHS28 = crRHS13*crRHS7;
const double crRHS29 = crRHS9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
const double crRHS30 = crRHS27 + crRHS28 + crRHS29;
const double crRHS31 = crRHS19*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crRHS25 + crRHS26 + crRHS27*rho + crRHS28*rho + crRHS29*rho);
const double crRHS32 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crRHS33 = crRHS2*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
const double crRHS34 = crRHS5*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2));
const double crRHS35 = crRHS7*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
const double crRHS36 = crRHS14*crRHS9;
const double crRHS37 = crRHS34 + crRHS35 + crRHS36;
const double crRHS38 = crRHS19*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crRHS32 + crRHS33 + crRHS34*rho + crRHS35*rho + crRHS36*rho);
const double crRHS39 = N[1]*rho;
const double crRHS40 = N[1]*crRHS2;
const double crRHS41 = crRHS22*crRHS39;
const double crRHS42 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS7 + DN(1,2)*crRHS9);
const double crRHS43 = N[2]*rho;
const double crRHS44 = N[2]*crRHS2;
const double crRHS45 = crRHS22*crRHS43;
const double crRHS46 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS7 + DN(2,2)*crRHS9);
const double crRHS47 = N[3]*rho;
const double crRHS48 = N[3]*crRHS2;
const double crRHS49 = crRHS22*crRHS47;
const double crRHS50 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS7 + DN(3,2)*crRHS9);
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS18 + DN(0,0)*stress[0] + DN(0,1)*stress[3] + DN(0,2)*stress[5] - N[0]*crRHS1 + N[0]*crRHS3 + crRHS11*crRHS12 - crRHS20*crRHS21 + crRHS20*crRHS23 + crRHS20*crRHS24);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS18 + DN(0,1)*stress[1] + DN(0,2)*stress[4] - N[0]*crRHS25 + N[0]*crRHS26 + crRHS12*crRHS30 - crRHS21*crRHS31 + crRHS23*crRHS31 + crRHS24*crRHS31);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress[5] + DN(0,1)*stress[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS18 + DN(0,2)*stress[2] - N[0]*crRHS32 + N[0]*crRHS33 + crRHS12*crRHS37 - crRHS21*crRHS38 + crRHS23*crRHS38 + crRHS24*crRHS38);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS20 + DN(0,1)*crRHS31 + DN(0,2)*crRHS38 + N[0]*crRHS15);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS18 + DN(1,0)*stress[0] + DN(1,1)*stress[3] + DN(1,2)*stress[5] - N[1]*crRHS1 + N[1]*crRHS3 + crRHS11*crRHS39 - crRHS20*crRHS40 + crRHS20*crRHS41 + crRHS20*crRHS42);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS18 + DN(1,1)*stress[1] + DN(1,2)*stress[4] - N[1]*crRHS25 + N[1]*crRHS26 + crRHS30*crRHS39 - crRHS31*crRHS40 + crRHS31*crRHS41 + crRHS31*crRHS42);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress[5] + DN(1,1)*stress[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS18 + DN(1,2)*stress[2] - N[1]*crRHS32 + N[1]*crRHS33 + crRHS37*crRHS39 - crRHS38*crRHS40 + crRHS38*crRHS41 + crRHS38*crRHS42);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS20 + DN(1,1)*crRHS31 + DN(1,2)*crRHS38 + N[1]*crRHS15);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS18 + DN(2,0)*stress[0] + DN(2,1)*stress[3] + DN(2,2)*stress[5] - N[2]*crRHS1 + N[2]*crRHS3 + crRHS11*crRHS43 - crRHS20*crRHS44 + crRHS20*crRHS45 + crRHS20*crRHS46);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS18 + DN(2,1)*stress[1] + DN(2,2)*stress[4] - N[2]*crRHS25 + N[2]*crRHS26 + crRHS30*crRHS43 - crRHS31*crRHS44 + crRHS31*crRHS45 + crRHS31*crRHS46);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress[5] + DN(2,1)*stress[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS18 + DN(2,2)*stress[2] - N[2]*crRHS32 + N[2]*crRHS33 + crRHS37*crRHS43 - crRHS38*crRHS44 + crRHS38*crRHS45 + crRHS38*crRHS46);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS20 + DN(2,1)*crRHS31 + DN(2,2)*crRHS38 + N[2]*crRHS15);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS18 + DN(3,0)*stress[0] + DN(3,1)*stress[3] + DN(3,2)*stress[5] - N[3]*crRHS1 + N[3]*crRHS3 + crRHS11*crRHS47 - crRHS20*crRHS48 + crRHS20*crRHS49 + crRHS20*crRHS50);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS18 + DN(3,1)*stress[1] + DN(3,2)*stress[4] - N[3]*crRHS25 + N[3]*crRHS26 + crRHS30*crRHS47 - crRHS31*crRHS48 + crRHS31*crRHS49 + crRHS31*crRHS50);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress[5] + DN(3,1)*stress[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS18 + DN(3,2)*stress[2] - N[3]*crRHS32 + N[3]*crRHS33 + crRHS37*crRHS47 - crRHS38*crRHS48 + crRHS38*crRHS49 + crRHS38*crRHS50);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS20 + DN(3,1)*crRHS31 + DN(3,2)*crRHS38 + N[3]*crRHS15);

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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    
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
const double crLHS9 = h*(crLHS5*h + crLHS8)/stab_c1 + mu;
const double crLHS10 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crLHS11 = pow(N[0], 2);
const double crLHS12 = crLHS11*rho;
const double crLHS13 = N[0]*rho;
const double crLHS14 = crLHS10*crLHS13;
const double crLHS15 = N[0]*crLHS4;
const double crLHS16 = DN(0,0)*crLHS6;
const double crLHS17 = DN(0,1)*crLHS7;
const double crLHS18 = -crLHS15 + crLHS16*rho + crLHS17*rho;
const double crLHS19 = -crLHS14 + crLHS18;
const double crLHS20 = 1.0/(crLHS5 + crLHS8/h + mu*stab_c1/pow(h, 2));
const double crLHS21 = crLHS15*crLHS20;
const double crLHS22 = crLHS16 + crLHS17;
const double crLHS23 = crLHS20*rho;
const double crLHS24 = crLHS22*crLHS23;
const double crLHS25 = crLHS11*crLHS4;
const double crLHS26 = crLHS13*crLHS22 + crLHS25;
const double crLHS27 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS28 = C(1,2)*DN(0,1);
const double crLHS29 = C(2,2)*DN(0,0) + crLHS28;
const double crLHS30 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1);
const double crLHS31 = DN(0,0)*crLHS9;
const double crLHS32 = DN(0,1)*crLHS31;
const double crLHS33 = crLHS23*crLHS25;
const double crLHS34 = crLHS20*pow(rho, 2);
const double crLHS35 = crLHS22*crLHS34;
const double crLHS36 = crLHS30*crLHS35;
const double crLHS37 = gauss_weight*(N[0] + crLHS21 - crLHS24);
const double crLHS38 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS39 = C(0,2)*DN(1,0);
const double crLHS40 = C(2,2)*DN(1,1) + crLHS39;
const double crLHS41 = N[1]*rho;
const double crLHS42 = crLHS10*crLHS41;
const double crLHS43 = N[1]*crLHS4;
const double crLHS44 = DN(1,0)*crLHS6;
const double crLHS45 = DN(1,1)*crLHS7;
const double crLHS46 = -crLHS43 + crLHS44*rho + crLHS45*rho;
const double crLHS47 = -crLHS42 + crLHS46;
const double crLHS48 = N[1]*crLHS15;
const double crLHS49 = crLHS22*crLHS41 + crLHS48;
const double crLHS50 = DN(0,0)*DN(1,0);
const double crLHS51 = N[1]*crLHS14 + crLHS50*crLHS9;
const double crLHS52 = C(0,1)*DN(1,1) + crLHS39;
const double crLHS53 = C(1,2)*DN(1,1);
const double crLHS54 = C(2,2)*DN(1,0) + crLHS53;
const double crLHS55 = DN(1,1)*crLHS31;
const double crLHS56 = crLHS13*crLHS30;
const double crLHS57 = crLHS23*crLHS30;
const double crLHS58 = N[1]*crLHS56 - crLHS48*crLHS57;
const double crLHS59 = DN(0,0)*N[1];
const double crLHS60 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS61 = C(0,2)*DN(2,0);
const double crLHS62 = C(2,2)*DN(2,1) + crLHS61;
const double crLHS63 = N[2]*rho;
const double crLHS64 = N[2]*crLHS4;
const double crLHS65 = DN(2,0)*crLHS6;
const double crLHS66 = DN(2,1)*crLHS7;
const double crLHS67 = -crLHS64 + crLHS65*rho + crLHS66*rho;
const double crLHS68 = -crLHS10*crLHS63 + crLHS67;
const double crLHS69 = N[2]*crLHS15;
const double crLHS70 = crLHS22*crLHS63 + crLHS69;
const double crLHS71 = DN(0,0)*DN(2,0);
const double crLHS72 = N[2]*crLHS14 + crLHS71*crLHS9;
const double crLHS73 = C(0,1)*DN(2,1) + crLHS61;
const double crLHS74 = C(1,2)*DN(2,1);
const double crLHS75 = C(2,2)*DN(2,0) + crLHS74;
const double crLHS76 = DN(2,1)*crLHS31;
const double crLHS77 = N[2]*crLHS56 - crLHS57*crLHS69;
const double crLHS78 = DN(0,0)*N[2];
const double crLHS79 = C(0,1)*DN(0,0) + crLHS28;
const double crLHS80 = DN(0,1)*v_ns(0,0) + DN(1,1)*v_ns(1,0) + DN(2,1)*v_ns(2,0);
const double crLHS81 = crLHS35*crLHS80;
const double crLHS82 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS83 = pow(DN(0,1), 2);
const double crLHS84 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crLHS85 = crLHS13*crLHS84;
const double crLHS86 = crLHS18 - crLHS85;
const double crLHS87 = C(0,1)*DN(1,0) + crLHS53;
const double crLHS88 = DN(0,1)*crLHS9;
const double crLHS89 = DN(1,0)*crLHS88;
const double crLHS90 = crLHS13*crLHS80;
const double crLHS91 = crLHS23*crLHS80;
const double crLHS92 = N[1]*crLHS90 - crLHS48*crLHS91;
const double crLHS93 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS94 = crLHS41*crLHS84;
const double crLHS95 = crLHS46 - crLHS94;
const double crLHS96 = DN(0,1)*DN(1,1);
const double crLHS97 = N[1]*crLHS85 + crLHS9*crLHS96;
const double crLHS98 = DN(0,1)*N[1];
const double crLHS99 = C(0,1)*DN(2,0) + crLHS74;
const double crLHS100 = DN(2,0)*crLHS88;
const double crLHS101 = N[2]*crLHS90 - crLHS69*crLHS91;
const double crLHS102 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS103 = -crLHS63*crLHS84 + crLHS67;
const double crLHS104 = DN(0,1)*DN(2,1);
const double crLHS105 = N[2]*crLHS85 + crLHS104*crLHS9;
const double crLHS106 = DN(0,1)*N[2];
const double crLHS107 = DN(0,0)*N[0];
const double crLHS108 = DN(0,1)*N[0];
const double crLHS109 = DN(0,0)*crLHS20;
const double crLHS110 = DN(0,1)*crLHS20;
const double crLHS111 = crLHS20*gauss_weight;
const double crLHS112 = DN(1,0)*N[0];
const double crLHS113 = DN(1,1)*N[0];
const double crLHS114 = crLHS111*(crLHS50 + crLHS96);
const double crLHS115 = DN(2,0)*N[0];
const double crLHS116 = DN(2,1)*N[0];
const double crLHS117 = crLHS111*(crLHS104 + crLHS71);
const double crLHS118 = crLHS20*crLHS43;
const double crLHS119 = crLHS44 + crLHS45;
const double crLHS120 = crLHS119*crLHS23;
const double crLHS121 = crLHS119*crLHS13 + crLHS48;
const double crLHS122 = crLHS119*crLHS34;
const double crLHS123 = crLHS122*crLHS30;
const double crLHS124 = pow(DN(1,0), 2);
const double crLHS125 = pow(N[1], 2);
const double crLHS126 = crLHS125*rho;
const double crLHS127 = crLHS125*crLHS4;
const double crLHS128 = crLHS119*crLHS41 + crLHS127;
const double crLHS129 = DN(1,0)*crLHS9;
const double crLHS130 = DN(1,1)*crLHS129;
const double crLHS131 = gauss_weight*(N[1] + crLHS118 - crLHS120);
const double crLHS132 = N[2]*crLHS43;
const double crLHS133 = crLHS119*crLHS63 + crLHS132;
const double crLHS134 = DN(1,0)*DN(2,0);
const double crLHS135 = N[2]*crLHS42 + crLHS134*crLHS9;
const double crLHS136 = DN(2,1)*crLHS129;
const double crLHS137 = N[2]*crLHS41;
const double crLHS138 = crLHS118*crLHS63;
const double crLHS139 = crLHS137*crLHS30 - crLHS138*crLHS30;
const double crLHS140 = DN(1,0)*N[2];
const double crLHS141 = crLHS122*crLHS80;
const double crLHS142 = pow(DN(1,1), 2);
const double crLHS143 = DN(2,0)*crLHS9;
const double crLHS144 = DN(1,1)*crLHS143;
const double crLHS145 = crLHS137*crLHS80 - crLHS138*crLHS80;
const double crLHS146 = DN(1,1)*DN(2,1);
const double crLHS147 = N[2]*crLHS94 + crLHS146*crLHS9;
const double crLHS148 = DN(1,1)*N[2];
const double crLHS149 = DN(1,0)*crLHS20;
const double crLHS150 = DN(1,1)*crLHS20;
const double crLHS151 = DN(1,0)*N[1];
const double crLHS152 = DN(1,1)*N[1];
const double crLHS153 = DN(2,0)*N[1];
const double crLHS154 = DN(2,1)*N[1];
const double crLHS155 = crLHS111*(crLHS134 + crLHS146);
const double crLHS156 = crLHS20*crLHS64;
const double crLHS157 = crLHS65 + crLHS66;
const double crLHS158 = crLHS157*crLHS23;
const double crLHS159 = crLHS13*crLHS157 + crLHS69;
const double crLHS160 = crLHS157*crLHS34;
const double crLHS161 = crLHS160*crLHS30;
const double crLHS162 = crLHS132 + crLHS157*crLHS41;
const double crLHS163 = pow(DN(2,0), 2);
const double crLHS164 = pow(N[2], 2);
const double crLHS165 = crLHS164*rho;
const double crLHS166 = crLHS164*crLHS4;
const double crLHS167 = crLHS157*crLHS63 + crLHS166;
const double crLHS168 = DN(2,1)*crLHS143;
const double crLHS169 = gauss_weight*(N[2] + crLHS156 - crLHS158);
const double crLHS170 = crLHS160*crLHS80;
const double crLHS171 = pow(DN(2,1), 2);
const double crLHS172 = DN(2,0)*crLHS20;
const double crLHS173 = DN(2,1)*crLHS20;
const double crLHS174 = DN(2,0)*N[2];
const double crLHS175 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS10*crLHS12 + crLHS19*crLHS21 - crLHS19*crLHS24 + crLHS26 + crLHS3*crLHS9);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS27 + DN(0,1)*crLHS29 + N[0]*crLHS36 + crLHS12*crLHS30 - crLHS30*crLHS33 + crLHS32);
rLHS(0,2)+=-DN(0,0)*crLHS37;
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS38 + DN(0,1)*crLHS40 + crLHS21*crLHS47 - crLHS24*crLHS47 + crLHS49 + crLHS51);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS52 + DN(0,1)*crLHS54 + N[1]*crLHS36 + crLHS55 + crLHS58);
rLHS(0,5)+=-gauss_weight*(DN(1,0)*crLHS21 - DN(1,0)*crLHS24 + crLHS59);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS60 + DN(0,1)*crLHS62 + crLHS21*crLHS68 - crLHS24*crLHS68 + crLHS70 + crLHS72);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS73 + DN(0,1)*crLHS75 + N[2]*crLHS36 + crLHS76 + crLHS77);
rLHS(0,8)+=-gauss_weight*(DN(2,0)*crLHS21 - DN(2,0)*crLHS24 + crLHS78);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS79 + N[0]*crLHS81 + crLHS12*crLHS80 + crLHS32 - crLHS33*crLHS80);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS29 + DN(0,1)*crLHS82 + crLHS12*crLHS84 + crLHS21*crLHS86 - crLHS24*crLHS86 + crLHS26 + crLHS83*crLHS9);
rLHS(1,2)+=-DN(0,1)*crLHS37;
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS87 + N[1]*crLHS81 + crLHS89 + crLHS92);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS93 + crLHS21*crLHS95 - crLHS24*crLHS95 + crLHS49 + crLHS97);
rLHS(1,5)+=-gauss_weight*(DN(1,1)*crLHS21 - DN(1,1)*crLHS24 + crLHS98);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS99 + N[2]*crLHS81 + crLHS100 + crLHS101);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS102 + crLHS103*crLHS21 - crLHS103*crLHS24 + crLHS105 + crLHS70);
rLHS(1,8)+=-gauss_weight*(DN(2,1)*crLHS21 - DN(2,1)*crLHS24 + crLHS106);
rLHS(2,0)+=gauss_weight*(crLHS107 + crLHS108*crLHS91 - crLHS109*crLHS19);
rLHS(2,1)+=gauss_weight*(crLHS107*crLHS57 + crLHS108 - crLHS110*crLHS86);
rLHS(2,2)+=crLHS111*(crLHS3 + crLHS83);
rLHS(2,3)+=gauss_weight*(-crLHS109*crLHS47 + crLHS112 + crLHS91*crLHS98);
rLHS(2,4)+=gauss_weight*(-crLHS110*crLHS95 + crLHS113 + crLHS57*crLHS59);
rLHS(2,5)+=crLHS114;
rLHS(2,6)+=gauss_weight*(crLHS106*crLHS91 - crLHS109*crLHS68 + crLHS115);
rLHS(2,7)+=gauss_weight*(-crLHS103*crLHS110 + crLHS116 + crLHS57*crLHS78);
rLHS(2,8)+=crLHS117;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS118*crLHS19 - crLHS120*crLHS19 + crLHS121 + crLHS51);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS27 + DN(1,1)*crLHS29 + N[0]*crLHS123 + crLHS58 + crLHS89);
rLHS(3,2)+=-gauss_weight*(-DN(0,0)*crLHS120 + crLHS109*crLHS43 + crLHS112);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS38 + DN(1,1)*crLHS40 + crLHS10*crLHS126 + crLHS118*crLHS47 - crLHS120*crLHS47 + crLHS124*crLHS9 + crLHS128);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS52 + DN(1,1)*crLHS54 + N[1]*crLHS123 + crLHS126*crLHS30 - crLHS127*crLHS57 + crLHS130);
rLHS(3,5)+=-DN(1,0)*crLHS131;
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS60 + DN(1,1)*crLHS62 + crLHS118*crLHS68 - crLHS120*crLHS68 + crLHS133 + crLHS135);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS73 + DN(1,1)*crLHS75 + N[2]*crLHS123 + crLHS136 + crLHS139);
rLHS(3,8)+=-gauss_weight*(DN(2,0)*crLHS118 - DN(2,0)*crLHS120 + crLHS140);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS79 + N[0]*crLHS141 + crLHS55 + crLHS92);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS29 + DN(1,1)*crLHS82 + crLHS118*crLHS86 - crLHS120*crLHS86 + crLHS121 + crLHS97);
rLHS(4,2)+=-gauss_weight*(-DN(0,1)*crLHS120 + crLHS110*crLHS43 + crLHS113);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS87 + N[1]*crLHS141 + crLHS126*crLHS80 - crLHS127*crLHS91 + crLHS130);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS93 + crLHS118*crLHS95 - crLHS120*crLHS95 + crLHS126*crLHS84 + crLHS128 + crLHS142*crLHS9);
rLHS(4,5)+=-DN(1,1)*crLHS131;
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS99 + N[2]*crLHS141 + crLHS144 + crLHS145);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS102 + crLHS103*crLHS118 - crLHS103*crLHS120 + crLHS133 + crLHS147);
rLHS(4,8)+=-gauss_weight*(DN(2,1)*crLHS118 - DN(2,1)*crLHS120 + crLHS148);
rLHS(5,0)+=gauss_weight*(crLHS113*crLHS91 - crLHS149*crLHS19 + crLHS59);
rLHS(5,1)+=gauss_weight*(crLHS112*crLHS57 - crLHS150*crLHS86 + crLHS98);
rLHS(5,2)+=crLHS114;
rLHS(5,3)+=gauss_weight*(-crLHS149*crLHS47 + crLHS151 + crLHS152*crLHS91);
rLHS(5,4)+=gauss_weight*(-crLHS150*crLHS95 + crLHS151*crLHS57 + crLHS152);
rLHS(5,5)+=crLHS111*(crLHS124 + crLHS142);
rLHS(5,6)+=gauss_weight*(crLHS148*crLHS91 - crLHS149*crLHS68 + crLHS153);
rLHS(5,7)+=gauss_weight*(-crLHS103*crLHS150 + crLHS140*crLHS57 + crLHS154);
rLHS(5,8)+=crLHS155;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS156*crLHS19 - crLHS158*crLHS19 + crLHS159 + crLHS72);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS27 + DN(2,1)*crLHS29 + N[0]*crLHS161 + crLHS100 + crLHS77);
rLHS(6,2)+=-gauss_weight*(-DN(0,0)*crLHS158 + crLHS109*crLHS64 + crLHS115);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS38 + DN(2,1)*crLHS40 + crLHS135 + crLHS156*crLHS47 - crLHS158*crLHS47 + crLHS162);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS52 + DN(2,1)*crLHS54 + N[1]*crLHS161 + crLHS139 + crLHS144);
rLHS(6,5)+=-gauss_weight*(-DN(1,0)*crLHS158 + crLHS149*crLHS64 + crLHS153);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS60 + DN(2,1)*crLHS62 + crLHS10*crLHS165 + crLHS156*crLHS68 - crLHS158*crLHS68 + crLHS163*crLHS9 + crLHS167);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS73 + DN(2,1)*crLHS75 + N[2]*crLHS161 + crLHS165*crLHS30 - crLHS166*crLHS57 + crLHS168);
rLHS(6,8)+=-DN(2,0)*crLHS169;
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS79 + N[0]*crLHS170 + crLHS101 + crLHS76);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS29 + DN(2,1)*crLHS82 + crLHS105 + crLHS156*crLHS86 - crLHS158*crLHS86 + crLHS159);
rLHS(7,2)+=-gauss_weight*(-DN(0,1)*crLHS158 + crLHS110*crLHS64 + crLHS116);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS87 + N[1]*crLHS170 + crLHS136 + crLHS145);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS93 + crLHS147 + crLHS156*crLHS95 - crLHS158*crLHS95 + crLHS162);
rLHS(7,5)+=-gauss_weight*(-DN(1,1)*crLHS158 + crLHS150*crLHS64 + crLHS154);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS99 + N[2]*crLHS170 + crLHS165*crLHS80 - crLHS166*crLHS91 + crLHS168);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS102 + crLHS103*crLHS156 - crLHS103*crLHS158 + crLHS165*crLHS84 + crLHS167 + crLHS171*crLHS9);
rLHS(7,8)+=-DN(2,1)*crLHS169;
rLHS(8,0)+=gauss_weight*(crLHS116*crLHS91 - crLHS172*crLHS19 + crLHS78);
rLHS(8,1)+=gauss_weight*(crLHS106 + crLHS115*crLHS57 - crLHS173*crLHS86);
rLHS(8,2)+=crLHS117;
rLHS(8,3)+=gauss_weight*(crLHS140 + crLHS154*crLHS91 - crLHS172*crLHS47);
rLHS(8,4)+=gauss_weight*(crLHS148 + crLHS153*crLHS57 - crLHS173*crLHS95);
rLHS(8,5)+=crLHS155;
rLHS(8,6)+=gauss_weight*(-crLHS172*crLHS68 + crLHS174 + crLHS175*crLHS91);
rLHS(8,7)+=gauss_weight*(-crLHS103*crLHS173 + crLHS174*crLHS57 + crLHS175);
rLHS(8,8)+=crLHS111*(crLHS163 + crLHS171);
    
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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    
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
const double crLHS12 = h*(crLHS11 + crLHS7*h)/stab_c1 + mu;
const double crLHS13 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crLHS14 = pow(N[0], 2);
const double crLHS15 = crLHS14*rho;
const double crLHS16 = N[0]*rho;
const double crLHS17 = crLHS13*crLHS16;
const double crLHS18 = N[0]*crLHS6;
const double crLHS19 = DN(0,0)*crLHS8;
const double crLHS20 = DN(0,1)*crLHS9;
const double crLHS21 = DN(0,2)*crLHS10;
const double crLHS22 = -crLHS18 + crLHS19*rho + crLHS20*rho + crLHS21*rho;
const double crLHS23 = -crLHS17 + crLHS22;
const double crLHS24 = 1.0/(crLHS11/h + crLHS7 + mu*stab_c1/pow(h, 2));
const double crLHS25 = crLHS18*crLHS24;
const double crLHS26 = crLHS19 + crLHS20 + crLHS21;
const double crLHS27 = crLHS24*rho;
const double crLHS28 = crLHS26*crLHS27;
const double crLHS29 = crLHS14*crLHS6;
const double crLHS30 = crLHS16*crLHS26 + crLHS29;
const double crLHS31 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS32 = C(1,3)*DN(0,1);
const double crLHS33 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS32;
const double crLHS34 = C(3,5)*DN(0,0);
const double crLHS35 = C(4,5)*DN(0,2);
const double crLHS36 = C(1,5)*DN(0,1) + crLHS34 + crLHS35;
const double crLHS37 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crLHS38 = DN(0,0)*crLHS12;
const double crLHS39 = DN(0,1)*crLHS38;
const double crLHS40 = crLHS27*crLHS29;
const double crLHS41 = crLHS24*pow(rho, 2);
const double crLHS42 = crLHS26*crLHS41;
const double crLHS43 = N[0]*crLHS42;
const double crLHS44 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS45 = C(3,4)*DN(0,1);
const double crLHS46 = C(2,3)*DN(0,2) + crLHS34 + crLHS45;
const double crLHS47 = C(2,5)*DN(0,2);
const double crLHS48 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS47;
const double crLHS49 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crLHS50 = DN(0,2)*crLHS38;
const double crLHS51 = gauss_weight*(N[0] + crLHS25 - crLHS28);
const double crLHS52 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS53 = C(0,3)*DN(1,0);
const double crLHS54 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS53;
const double crLHS55 = C(0,5)*DN(1,0);
const double crLHS56 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS55;
const double crLHS57 = N[1]*rho;
const double crLHS58 = crLHS13*crLHS57;
const double crLHS59 = N[1]*crLHS6;
const double crLHS60 = DN(1,0)*crLHS8;
const double crLHS61 = DN(1,1)*crLHS9;
const double crLHS62 = DN(1,2)*crLHS10;
const double crLHS63 = -crLHS59 + crLHS60*rho + crLHS61*rho + crLHS62*rho;
const double crLHS64 = -crLHS58 + crLHS63;
const double crLHS65 = N[1]*crLHS18;
const double crLHS66 = crLHS26*crLHS57 + crLHS65;
const double crLHS67 = DN(0,0)*DN(1,0);
const double crLHS68 = N[1]*crLHS17 + crLHS12*crLHS67;
const double crLHS69 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS53;
const double crLHS70 = C(1,3)*DN(1,1);
const double crLHS71 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS70;
const double crLHS72 = C(3,5)*DN(1,0);
const double crLHS73 = C(4,5)*DN(1,2);
const double crLHS74 = C(1,5)*DN(1,1) + crLHS72 + crLHS73;
const double crLHS75 = DN(1,1)*crLHS38;
const double crLHS76 = N[1]*crLHS42;
const double crLHS77 = N[1]*crLHS16;
const double crLHS78 = crLHS27*crLHS65;
const double crLHS79 = crLHS37*crLHS77 - crLHS37*crLHS78;
const double crLHS80 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS55;
const double crLHS81 = C(3,4)*DN(1,1);
const double crLHS82 = C(2,3)*DN(1,2) + crLHS72 + crLHS81;
const double crLHS83 = C(2,5)*DN(1,2);
const double crLHS84 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS83;
const double crLHS85 = DN(1,2)*crLHS38;
const double crLHS86 = crLHS49*crLHS77 - crLHS49*crLHS78;
const double crLHS87 = DN(0,0)*N[1];
const double crLHS88 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS89 = C(0,3)*DN(2,0);
const double crLHS90 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS89;
const double crLHS91 = C(0,5)*DN(2,0);
const double crLHS92 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS91;
const double crLHS93 = N[2]*rho;
const double crLHS94 = crLHS13*crLHS93;
const double crLHS95 = N[2]*crLHS6;
const double crLHS96 = DN(2,0)*crLHS8;
const double crLHS97 = DN(2,1)*crLHS9;
const double crLHS98 = DN(2,2)*crLHS10;
const double crLHS99 = -crLHS95 + crLHS96*rho + crLHS97*rho + crLHS98*rho;
const double crLHS100 = -crLHS94 + crLHS99;
const double crLHS101 = N[2]*crLHS18;
const double crLHS102 = crLHS101 + crLHS26*crLHS93;
const double crLHS103 = DN(0,0)*DN(2,0);
const double crLHS104 = N[2]*crLHS17 + crLHS103*crLHS12;
const double crLHS105 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS89;
const double crLHS106 = C(1,3)*DN(2,1);
const double crLHS107 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS106;
const double crLHS108 = C(3,5)*DN(2,0);
const double crLHS109 = C(4,5)*DN(2,2);
const double crLHS110 = C(1,5)*DN(2,1) + crLHS108 + crLHS109;
const double crLHS111 = DN(2,1)*crLHS38;
const double crLHS112 = N[2]*crLHS42;
const double crLHS113 = N[2]*crLHS16;
const double crLHS114 = crLHS101*crLHS27;
const double crLHS115 = crLHS113*crLHS37 - crLHS114*crLHS37;
const double crLHS116 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS91;
const double crLHS117 = C(3,4)*DN(2,1);
const double crLHS118 = C(2,3)*DN(2,2) + crLHS108 + crLHS117;
const double crLHS119 = C(2,5)*DN(2,2);
const double crLHS120 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS119;
const double crLHS121 = DN(2,2)*crLHS38;
const double crLHS122 = crLHS113*crLHS49 - crLHS114*crLHS49;
const double crLHS123 = DN(0,0)*N[2];
const double crLHS124 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS125 = C(0,3)*DN(3,0);
const double crLHS126 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS125;
const double crLHS127 = C(0,5)*DN(3,0);
const double crLHS128 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS127;
const double crLHS129 = N[3]*rho;
const double crLHS130 = N[3]*crLHS6;
const double crLHS131 = DN(3,0)*crLHS8;
const double crLHS132 = DN(3,1)*crLHS9;
const double crLHS133 = DN(3,2)*crLHS10;
const double crLHS134 = -crLHS130 + crLHS131*rho + crLHS132*rho + crLHS133*rho;
const double crLHS135 = -crLHS129*crLHS13 + crLHS134;
const double crLHS136 = N[3]*crLHS18;
const double crLHS137 = crLHS129*crLHS26 + crLHS136;
const double crLHS138 = DN(0,0)*DN(3,0);
const double crLHS139 = N[3]*crLHS17 + crLHS12*crLHS138;
const double crLHS140 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS125;
const double crLHS141 = C(1,3)*DN(3,1);
const double crLHS142 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS141;
const double crLHS143 = C(3,5)*DN(3,0);
const double crLHS144 = C(4,5)*DN(3,2);
const double crLHS145 = C(1,5)*DN(3,1) + crLHS143 + crLHS144;
const double crLHS146 = DN(3,1)*crLHS38;
const double crLHS147 = N[3]*crLHS42;
const double crLHS148 = N[3]*crLHS16;
const double crLHS149 = crLHS136*crLHS27;
const double crLHS150 = crLHS148*crLHS37 - crLHS149*crLHS37;
const double crLHS151 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS127;
const double crLHS152 = C(3,4)*DN(3,1);
const double crLHS153 = C(2,3)*DN(3,2) + crLHS143 + crLHS152;
const double crLHS154 = C(2,5)*DN(3,2);
const double crLHS155 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS154;
const double crLHS156 = DN(3,2)*crLHS38;
const double crLHS157 = crLHS148*crLHS49 - crLHS149*crLHS49;
const double crLHS158 = DN(0,0)*N[3];
const double crLHS159 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS32;
const double crLHS160 = C(0,4)*DN(0,0) + crLHS35 + crLHS45;
const double crLHS161 = DN(0,1)*v_ns(0,0) + DN(1,1)*v_ns(1,0) + DN(2,1)*v_ns(2,0) + DN(3,1)*v_ns(3,0);
const double crLHS162 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS163 = C(1,4)*DN(0,1);
const double crLHS164 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS163;
const double crLHS165 = pow(DN(0,1), 2);
const double crLHS166 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crLHS167 = crLHS16*crLHS166;
const double crLHS168 = -crLHS167 + crLHS22;
const double crLHS169 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS163;
const double crLHS170 = C(2,4)*DN(0,2);
const double crLHS171 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS170;
const double crLHS172 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crLHS173 = DN(0,1)*crLHS12;
const double crLHS174 = DN(0,2)*crLHS173;
const double crLHS175 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS70;
const double crLHS176 = C(0,4)*DN(1,0) + crLHS73 + crLHS81;
const double crLHS177 = DN(1,0)*crLHS173;
const double crLHS178 = crLHS161*crLHS77 - crLHS161*crLHS78;
const double crLHS179 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS180 = C(1,4)*DN(1,1);
const double crLHS181 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS180;
const double crLHS182 = crLHS166*crLHS57;
const double crLHS183 = -crLHS182 + crLHS63;
const double crLHS184 = DN(0,1)*DN(1,1);
const double crLHS185 = N[1]*crLHS167 + crLHS12*crLHS184;
const double crLHS186 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS180;
const double crLHS187 = C(2,4)*DN(1,2);
const double crLHS188 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS187;
const double crLHS189 = DN(1,2)*crLHS173;
const double crLHS190 = crLHS172*crLHS77 - crLHS172*crLHS78;
const double crLHS191 = DN(0,1)*N[1];
const double crLHS192 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS106;
const double crLHS193 = C(0,4)*DN(2,0) + crLHS109 + crLHS117;
const double crLHS194 = DN(2,0)*crLHS173;
const double crLHS195 = crLHS113*crLHS161 - crLHS114*crLHS161;
const double crLHS196 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS197 = C(1,4)*DN(2,1);
const double crLHS198 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS197;
const double crLHS199 = crLHS166*crLHS93;
const double crLHS200 = -crLHS199 + crLHS99;
const double crLHS201 = DN(0,1)*DN(2,1);
const double crLHS202 = N[2]*crLHS167 + crLHS12*crLHS201;
const double crLHS203 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS197;
const double crLHS204 = C(2,4)*DN(2,2);
const double crLHS205 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS204;
const double crLHS206 = DN(2,2)*crLHS173;
const double crLHS207 = crLHS113*crLHS172 - crLHS114*crLHS172;
const double crLHS208 = DN(0,1)*N[2];
const double crLHS209 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS141;
const double crLHS210 = C(0,4)*DN(3,0) + crLHS144 + crLHS152;
const double crLHS211 = DN(3,0)*crLHS173;
const double crLHS212 = crLHS148*crLHS161 - crLHS149*crLHS161;
const double crLHS213 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS214 = C(1,4)*DN(3,1);
const double crLHS215 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS214;
const double crLHS216 = -crLHS129*crLHS166 + crLHS134;
const double crLHS217 = DN(0,1)*DN(3,1);
const double crLHS218 = N[3]*crLHS167 + crLHS12*crLHS217;
const double crLHS219 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS214;
const double crLHS220 = C(2,4)*DN(3,2);
const double crLHS221 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS220;
const double crLHS222 = DN(3,2)*crLHS173;
const double crLHS223 = crLHS148*crLHS172 - crLHS149*crLHS172;
const double crLHS224 = DN(0,1)*N[3];
const double crLHS225 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS47;
const double crLHS226 = DN(0,2)*v_ns(0,0) + DN(1,2)*v_ns(1,0) + DN(2,2)*v_ns(2,0) + DN(3,2)*v_ns(3,0);
const double crLHS227 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS170;
const double crLHS228 = DN(0,2)*v_ns(0,1) + DN(1,2)*v_ns(1,1) + DN(2,2)*v_ns(2,1) + DN(3,2)*v_ns(3,1);
const double crLHS229 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS230 = pow(DN(0,2), 2);
const double crLHS231 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crLHS232 = crLHS16*crLHS231;
const double crLHS233 = crLHS22 - crLHS232;
const double crLHS234 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS83;
const double crLHS235 = DN(0,2)*crLHS12;
const double crLHS236 = DN(1,0)*crLHS235;
const double crLHS237 = crLHS226*crLHS77 - crLHS226*crLHS78;
const double crLHS238 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS187;
const double crLHS239 = DN(1,1)*crLHS235;
const double crLHS240 = crLHS228*crLHS77 - crLHS228*crLHS78;
const double crLHS241 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS242 = crLHS231*crLHS57;
const double crLHS243 = -crLHS242 + crLHS63;
const double crLHS244 = DN(0,2)*DN(1,2);
const double crLHS245 = N[1]*crLHS232 + crLHS12*crLHS244;
const double crLHS246 = DN(0,2)*N[1];
const double crLHS247 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS119;
const double crLHS248 = DN(2,0)*crLHS235;
const double crLHS249 = crLHS113*crLHS226 - crLHS114*crLHS226;
const double crLHS250 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS204;
const double crLHS251 = DN(2,1)*crLHS235;
const double crLHS252 = crLHS113*crLHS228 - crLHS114*crLHS228;
const double crLHS253 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS254 = crLHS231*crLHS93;
const double crLHS255 = -crLHS254 + crLHS99;
const double crLHS256 = DN(0,2)*DN(2,2);
const double crLHS257 = N[2]*crLHS232 + crLHS12*crLHS256;
const double crLHS258 = DN(0,2)*N[2];
const double crLHS259 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS154;
const double crLHS260 = DN(3,0)*crLHS235;
const double crLHS261 = crLHS148*crLHS226 - crLHS149*crLHS226;
const double crLHS262 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS220;
const double crLHS263 = DN(3,1)*crLHS235;
const double crLHS264 = crLHS148*crLHS228 - crLHS149*crLHS228;
const double crLHS265 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS266 = -crLHS129*crLHS231 + crLHS134;
const double crLHS267 = DN(0,2)*DN(3,2);
const double crLHS268 = N[3]*crLHS232 + crLHS12*crLHS267;
const double crLHS269 = DN(0,2)*N[3];
const double crLHS270 = DN(0,0)*N[0];
const double crLHS271 = DN(0,1)*N[0];
const double crLHS272 = crLHS27*crLHS271;
const double crLHS273 = DN(0,2)*N[0];
const double crLHS274 = crLHS27*crLHS273;
const double crLHS275 = DN(0,0)*crLHS24;
const double crLHS276 = crLHS27*crLHS270;
const double crLHS277 = DN(0,1)*crLHS24;
const double crLHS278 = DN(0,2)*crLHS24;
const double crLHS279 = crLHS24*gauss_weight;
const double crLHS280 = DN(1,0)*N[0];
const double crLHS281 = crLHS191*crLHS27;
const double crLHS282 = crLHS246*crLHS27;
const double crLHS283 = DN(1,1)*N[0];
const double crLHS284 = crLHS27*crLHS87;
const double crLHS285 = DN(1,2)*N[0];
const double crLHS286 = crLHS279*(crLHS184 + crLHS244 + crLHS67);
const double crLHS287 = DN(2,0)*N[0];
const double crLHS288 = crLHS208*crLHS27;
const double crLHS289 = crLHS258*crLHS27;
const double crLHS290 = DN(2,1)*N[0];
const double crLHS291 = crLHS123*crLHS27;
const double crLHS292 = DN(2,2)*N[0];
const double crLHS293 = crLHS279*(crLHS103 + crLHS201 + crLHS256);
const double crLHS294 = DN(3,0)*N[0];
const double crLHS295 = crLHS224*crLHS27;
const double crLHS296 = crLHS269*crLHS27;
const double crLHS297 = DN(3,1)*N[0];
const double crLHS298 = crLHS158*crLHS27;
const double crLHS299 = DN(3,2)*N[0];
const double crLHS300 = crLHS279*(crLHS138 + crLHS217 + crLHS267);
const double crLHS301 = crLHS24*crLHS59;
const double crLHS302 = crLHS60 + crLHS61 + crLHS62;
const double crLHS303 = crLHS27*crLHS302;
const double crLHS304 = crLHS16*crLHS302 + crLHS65;
const double crLHS305 = crLHS302*crLHS41;
const double crLHS306 = N[0]*crLHS305;
const double crLHS307 = pow(DN(1,0), 2);
const double crLHS308 = pow(N[1], 2);
const double crLHS309 = crLHS308*rho;
const double crLHS310 = crLHS308*crLHS6;
const double crLHS311 = crLHS302*crLHS57 + crLHS310;
const double crLHS312 = DN(1,0)*crLHS12;
const double crLHS313 = DN(1,1)*crLHS312;
const double crLHS314 = crLHS27*crLHS310;
const double crLHS315 = N[1]*crLHS305;
const double crLHS316 = DN(1,2)*crLHS312;
const double crLHS317 = gauss_weight*(N[1] + crLHS301 - crLHS303);
const double crLHS318 = N[2]*crLHS59;
const double crLHS319 = crLHS302*crLHS93 + crLHS318;
const double crLHS320 = DN(1,0)*DN(2,0);
const double crLHS321 = N[2]*crLHS58 + crLHS12*crLHS320;
const double crLHS322 = DN(2,1)*crLHS312;
const double crLHS323 = N[2]*crLHS305;
const double crLHS324 = N[2]*crLHS57;
const double crLHS325 = crLHS301*crLHS93;
const double crLHS326 = crLHS324*crLHS37 - crLHS325*crLHS37;
const double crLHS327 = DN(2,2)*crLHS312;
const double crLHS328 = crLHS324*crLHS49 - crLHS325*crLHS49;
const double crLHS329 = DN(1,0)*N[2];
const double crLHS330 = N[3]*crLHS59;
const double crLHS331 = crLHS129*crLHS302 + crLHS330;
const double crLHS332 = DN(1,0)*DN(3,0);
const double crLHS333 = N[3]*crLHS58 + crLHS12*crLHS332;
const double crLHS334 = DN(3,1)*crLHS312;
const double crLHS335 = N[3]*crLHS305;
const double crLHS336 = N[3]*crLHS57;
const double crLHS337 = crLHS129*crLHS301;
const double crLHS338 = crLHS336*crLHS37 - crLHS337*crLHS37;
const double crLHS339 = DN(3,2)*crLHS312;
const double crLHS340 = crLHS336*crLHS49 - crLHS337*crLHS49;
const double crLHS341 = DN(1,0)*N[3];
const double crLHS342 = pow(DN(1,1), 2);
const double crLHS343 = DN(1,1)*crLHS12;
const double crLHS344 = DN(1,2)*crLHS343;
const double crLHS345 = DN(2,0)*crLHS343;
const double crLHS346 = crLHS161*crLHS324 - crLHS161*crLHS325;
const double crLHS347 = DN(1,1)*DN(2,1);
const double crLHS348 = N[2]*crLHS182 + crLHS12*crLHS347;
const double crLHS349 = DN(2,2)*crLHS343;
const double crLHS350 = crLHS172*crLHS324 - crLHS172*crLHS325;
const double crLHS351 = DN(1,1)*N[2];
const double crLHS352 = DN(3,0)*crLHS343;
const double crLHS353 = crLHS161*crLHS336 - crLHS161*crLHS337;
const double crLHS354 = DN(1,1)*DN(3,1);
const double crLHS355 = N[3]*crLHS182 + crLHS12*crLHS354;
const double crLHS356 = DN(3,2)*crLHS343;
const double crLHS357 = crLHS172*crLHS336 - crLHS172*crLHS337;
const double crLHS358 = DN(1,1)*N[3];
const double crLHS359 = pow(DN(1,2), 2);
const double crLHS360 = DN(1,2)*crLHS12;
const double crLHS361 = DN(2,0)*crLHS360;
const double crLHS362 = crLHS226*crLHS324 - crLHS226*crLHS325;
const double crLHS363 = DN(2,1)*crLHS360;
const double crLHS364 = crLHS228*crLHS324 - crLHS228*crLHS325;
const double crLHS365 = DN(1,2)*DN(2,2);
const double crLHS366 = N[2]*crLHS242 + crLHS12*crLHS365;
const double crLHS367 = DN(1,2)*N[2];
const double crLHS368 = DN(3,0)*crLHS360;
const double crLHS369 = crLHS226*crLHS336 - crLHS226*crLHS337;
const double crLHS370 = DN(3,1)*crLHS360;
const double crLHS371 = crLHS228*crLHS336 - crLHS228*crLHS337;
const double crLHS372 = DN(1,2)*DN(3,2);
const double crLHS373 = N[3]*crLHS242 + crLHS12*crLHS372;
const double crLHS374 = DN(1,2)*N[3];
const double crLHS375 = crLHS27*crLHS283;
const double crLHS376 = crLHS27*crLHS285;
const double crLHS377 = DN(1,0)*crLHS24;
const double crLHS378 = crLHS27*crLHS280;
const double crLHS379 = DN(1,1)*crLHS24;
const double crLHS380 = DN(1,2)*crLHS24;
const double crLHS381 = DN(1,0)*N[1];
const double crLHS382 = DN(1,1)*N[1];
const double crLHS383 = crLHS27*crLHS382;
const double crLHS384 = DN(1,2)*N[1];
const double crLHS385 = crLHS27*crLHS384;
const double crLHS386 = crLHS27*crLHS381;
const double crLHS387 = DN(2,0)*N[1];
const double crLHS388 = crLHS27*crLHS351;
const double crLHS389 = crLHS27*crLHS367;
const double crLHS390 = DN(2,1)*N[1];
const double crLHS391 = crLHS27*crLHS329;
const double crLHS392 = DN(2,2)*N[1];
const double crLHS393 = crLHS279*(crLHS320 + crLHS347 + crLHS365);
const double crLHS394 = DN(3,0)*N[1];
const double crLHS395 = crLHS27*crLHS358;
const double crLHS396 = crLHS27*crLHS374;
const double crLHS397 = DN(3,1)*N[1];
const double crLHS398 = crLHS27*crLHS341;
const double crLHS399 = DN(3,2)*N[1];
const double crLHS400 = crLHS279*(crLHS332 + crLHS354 + crLHS372);
const double crLHS401 = crLHS24*crLHS95;
const double crLHS402 = crLHS96 + crLHS97 + crLHS98;
const double crLHS403 = crLHS27*crLHS402;
const double crLHS404 = crLHS101 + crLHS16*crLHS402;
const double crLHS405 = crLHS402*crLHS41;
const double crLHS406 = N[0]*crLHS405;
const double crLHS407 = crLHS318 + crLHS402*crLHS57;
const double crLHS408 = N[1]*crLHS405;
const double crLHS409 = pow(DN(2,0), 2);
const double crLHS410 = pow(N[2], 2);
const double crLHS411 = crLHS410*rho;
const double crLHS412 = crLHS410*crLHS6;
const double crLHS413 = crLHS402*crLHS93 + crLHS412;
const double crLHS414 = DN(2,0)*crLHS12;
const double crLHS415 = DN(2,1)*crLHS414;
const double crLHS416 = crLHS27*crLHS412;
const double crLHS417 = N[2]*crLHS405;
const double crLHS418 = DN(2,2)*crLHS414;
const double crLHS419 = gauss_weight*(N[2] + crLHS401 - crLHS403);
const double crLHS420 = N[3]*crLHS95;
const double crLHS421 = crLHS129*crLHS402 + crLHS420;
const double crLHS422 = DN(2,0)*DN(3,0);
const double crLHS423 = N[3]*crLHS94 + crLHS12*crLHS422;
const double crLHS424 = DN(3,1)*crLHS414;
const double crLHS425 = N[3]*crLHS405;
const double crLHS426 = N[3]*crLHS93;
const double crLHS427 = crLHS129*crLHS401;
const double crLHS428 = crLHS37*crLHS426 - crLHS37*crLHS427;
const double crLHS429 = DN(3,2)*crLHS414;
const double crLHS430 = crLHS426*crLHS49 - crLHS427*crLHS49;
const double crLHS431 = DN(2,0)*N[3];
const double crLHS432 = pow(DN(2,1), 2);
const double crLHS433 = DN(2,1)*crLHS12;
const double crLHS434 = DN(2,2)*crLHS433;
const double crLHS435 = DN(3,0)*crLHS433;
const double crLHS436 = crLHS161*crLHS426 - crLHS161*crLHS427;
const double crLHS437 = DN(2,1)*DN(3,1);
const double crLHS438 = N[3]*crLHS199 + crLHS12*crLHS437;
const double crLHS439 = DN(3,2)*crLHS433;
const double crLHS440 = crLHS172*crLHS426 - crLHS172*crLHS427;
const double crLHS441 = DN(2,1)*N[3];
const double crLHS442 = pow(DN(2,2), 2);
const double crLHS443 = DN(2,2)*crLHS12;
const double crLHS444 = DN(3,0)*crLHS443;
const double crLHS445 = crLHS226*crLHS426 - crLHS226*crLHS427;
const double crLHS446 = DN(3,1)*crLHS443;
const double crLHS447 = crLHS228*crLHS426 - crLHS228*crLHS427;
const double crLHS448 = DN(2,2)*DN(3,2);
const double crLHS449 = N[3]*crLHS254 + crLHS12*crLHS448;
const double crLHS450 = DN(2,2)*N[3];
const double crLHS451 = crLHS27*crLHS290;
const double crLHS452 = crLHS27*crLHS292;
const double crLHS453 = DN(2,0)*crLHS24;
const double crLHS454 = crLHS27*crLHS287;
const double crLHS455 = DN(2,1)*crLHS24;
const double crLHS456 = DN(2,2)*crLHS24;
const double crLHS457 = crLHS27*crLHS390;
const double crLHS458 = crLHS27*crLHS392;
const double crLHS459 = crLHS27*crLHS387;
const double crLHS460 = DN(2,0)*N[2];
const double crLHS461 = DN(2,1)*N[2];
const double crLHS462 = crLHS27*crLHS461;
const double crLHS463 = DN(2,2)*N[2];
const double crLHS464 = crLHS27*crLHS463;
const double crLHS465 = crLHS27*crLHS460;
const double crLHS466 = DN(3,0)*N[2];
const double crLHS467 = crLHS27*crLHS441;
const double crLHS468 = crLHS27*crLHS450;
const double crLHS469 = DN(3,1)*N[2];
const double crLHS470 = crLHS27*crLHS431;
const double crLHS471 = DN(3,2)*N[2];
const double crLHS472 = crLHS279*(crLHS422 + crLHS437 + crLHS448);
const double crLHS473 = crLHS130*crLHS24;
const double crLHS474 = crLHS131 + crLHS132 + crLHS133;
const double crLHS475 = crLHS27*crLHS474;
const double crLHS476 = crLHS136 + crLHS16*crLHS474;
const double crLHS477 = crLHS41*crLHS474;
const double crLHS478 = N[0]*crLHS477;
const double crLHS479 = crLHS330 + crLHS474*crLHS57;
const double crLHS480 = N[1]*crLHS477;
const double crLHS481 = crLHS420 + crLHS474*crLHS93;
const double crLHS482 = N[2]*crLHS477;
const double crLHS483 = pow(DN(3,0), 2);
const double crLHS484 = pow(N[3], 2);
const double crLHS485 = crLHS484*rho;
const double crLHS486 = crLHS484*crLHS6;
const double crLHS487 = crLHS129*crLHS474 + crLHS486;
const double crLHS488 = DN(3,0)*crLHS12;
const double crLHS489 = DN(3,1)*crLHS488;
const double crLHS490 = crLHS27*crLHS486;
const double crLHS491 = N[3]*crLHS477;
const double crLHS492 = DN(3,2)*crLHS488;
const double crLHS493 = gauss_weight*(N[3] + crLHS473 - crLHS475);
const double crLHS494 = pow(DN(3,1), 2);
const double crLHS495 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS496 = pow(DN(3,2), 2);
const double crLHS497 = crLHS27*crLHS297;
const double crLHS498 = crLHS27*crLHS299;
const double crLHS499 = DN(3,0)*crLHS24;
const double crLHS500 = crLHS27*crLHS294;
const double crLHS501 = DN(3,1)*crLHS24;
const double crLHS502 = DN(3,2)*crLHS24;
const double crLHS503 = crLHS27*crLHS397;
const double crLHS504 = crLHS27*crLHS399;
const double crLHS505 = crLHS27*crLHS394;
const double crLHS506 = crLHS27*crLHS469;
const double crLHS507 = crLHS27*crLHS471;
const double crLHS508 = crLHS27*crLHS466;
const double crLHS509 = DN(3,0)*N[3];
const double crLHS510 = DN(3,1)*N[3];
const double crLHS511 = crLHS27*crLHS510;
const double crLHS512 = DN(3,2)*N[3];
const double crLHS513 = crLHS27*crLHS512;
const double crLHS514 = crLHS27*crLHS509;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS5 + crLHS13*crLHS15 + crLHS23*crLHS25 - crLHS23*crLHS28 + crLHS30);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS31 + DN(0,1)*crLHS33 + DN(0,2)*crLHS36 + crLHS15*crLHS37 - crLHS37*crLHS40 + crLHS37*crLHS43 + crLHS39);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS44 + DN(0,1)*crLHS46 + DN(0,2)*crLHS48 + crLHS15*crLHS49 - crLHS40*crLHS49 + crLHS43*crLHS49 + crLHS50);
rLHS(0,3)+=-DN(0,0)*crLHS51;
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS52 + DN(0,1)*crLHS54 + DN(0,2)*crLHS56 + crLHS25*crLHS64 - crLHS28*crLHS64 + crLHS66 + crLHS68);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS69 + DN(0,1)*crLHS71 + DN(0,2)*crLHS74 + crLHS37*crLHS76 + crLHS75 + crLHS79);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS80 + DN(0,1)*crLHS82 + DN(0,2)*crLHS84 + crLHS49*crLHS76 + crLHS85 + crLHS86);
rLHS(0,7)+=-gauss_weight*(DN(1,0)*crLHS25 - DN(1,0)*crLHS28 + crLHS87);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS88 + DN(0,1)*crLHS90 + DN(0,2)*crLHS92 + crLHS100*crLHS25 - crLHS100*crLHS28 + crLHS102 + crLHS104);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS105 + DN(0,1)*crLHS107 + DN(0,2)*crLHS110 + crLHS111 + crLHS112*crLHS37 + crLHS115);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS116 + DN(0,1)*crLHS118 + DN(0,2)*crLHS120 + crLHS112*crLHS49 + crLHS121 + crLHS122);
rLHS(0,11)+=-gauss_weight*(DN(2,0)*crLHS25 - DN(2,0)*crLHS28 + crLHS123);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS124 + DN(0,1)*crLHS126 + DN(0,2)*crLHS128 + crLHS135*crLHS25 - crLHS135*crLHS28 + crLHS137 + crLHS139);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS140 + DN(0,1)*crLHS142 + DN(0,2)*crLHS145 + crLHS146 + crLHS147*crLHS37 + crLHS150);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS151 + DN(0,1)*crLHS153 + DN(0,2)*crLHS155 + crLHS147*crLHS49 + crLHS156 + crLHS157);
rLHS(0,15)+=-gauss_weight*(DN(3,0)*crLHS25 - DN(3,0)*crLHS28 + crLHS158);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS159 + DN(0,2)*crLHS160 + crLHS15*crLHS161 - crLHS161*crLHS40 + crLHS161*crLHS43 + crLHS39);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS33 + DN(0,1)*crLHS162 + DN(0,2)*crLHS164 + crLHS12*crLHS165 + crLHS15*crLHS166 + crLHS168*crLHS25 - crLHS168*crLHS28 + crLHS30);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS46 + DN(0,1)*crLHS169 + DN(0,2)*crLHS171 + crLHS15*crLHS172 - crLHS172*crLHS40 + crLHS172*crLHS43 + crLHS174);
rLHS(1,3)+=-DN(0,1)*crLHS51;
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS54 + DN(0,1)*crLHS175 + DN(0,2)*crLHS176 + crLHS161*crLHS76 + crLHS177 + crLHS178);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS71 + DN(0,1)*crLHS179 + DN(0,2)*crLHS181 + crLHS183*crLHS25 - crLHS183*crLHS28 + crLHS185 + crLHS66);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS82 + DN(0,1)*crLHS186 + DN(0,2)*crLHS188 + crLHS172*crLHS76 + crLHS189 + crLHS190);
rLHS(1,7)+=-gauss_weight*(DN(1,1)*crLHS25 - DN(1,1)*crLHS28 + crLHS191);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS90 + DN(0,1)*crLHS192 + DN(0,2)*crLHS193 + crLHS112*crLHS161 + crLHS194 + crLHS195);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS107 + DN(0,1)*crLHS196 + DN(0,2)*crLHS198 + crLHS102 + crLHS200*crLHS25 - crLHS200*crLHS28 + crLHS202);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS118 + DN(0,1)*crLHS203 + DN(0,2)*crLHS205 + crLHS112*crLHS172 + crLHS206 + crLHS207);
rLHS(1,11)+=-gauss_weight*(DN(2,1)*crLHS25 - DN(2,1)*crLHS28 + crLHS208);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS126 + DN(0,1)*crLHS209 + DN(0,2)*crLHS210 + crLHS147*crLHS161 + crLHS211 + crLHS212);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS142 + DN(0,1)*crLHS213 + DN(0,2)*crLHS215 + crLHS137 + crLHS216*crLHS25 - crLHS216*crLHS28 + crLHS218);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS153 + DN(0,1)*crLHS219 + DN(0,2)*crLHS221 + crLHS147*crLHS172 + crLHS222 + crLHS223);
rLHS(1,15)+=-gauss_weight*(DN(3,1)*crLHS25 - DN(3,1)*crLHS28 + crLHS224);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS160 + DN(0,2)*crLHS225 + crLHS15*crLHS226 - crLHS226*crLHS40 + crLHS226*crLHS43 + crLHS50);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS36 + DN(0,1)*crLHS164 + DN(0,2)*crLHS227 + crLHS15*crLHS228 + crLHS174 - crLHS228*crLHS40 + crLHS228*crLHS43);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS48 + DN(0,1)*crLHS171 + DN(0,2)*crLHS229 + crLHS12*crLHS230 + crLHS15*crLHS231 + crLHS233*crLHS25 - crLHS233*crLHS28 + crLHS30);
rLHS(2,3)+=-DN(0,2)*crLHS51;
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS56 + DN(0,1)*crLHS176 + DN(0,2)*crLHS234 + crLHS226*crLHS76 + crLHS236 + crLHS237);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS181 + DN(0,2)*crLHS238 + crLHS228*crLHS76 + crLHS239 + crLHS240);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS84 + DN(0,1)*crLHS188 + DN(0,2)*crLHS241 + crLHS243*crLHS25 - crLHS243*crLHS28 + crLHS245 + crLHS66);
rLHS(2,7)+=-gauss_weight*(DN(1,2)*crLHS25 - DN(1,2)*crLHS28 + crLHS246);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS92 + DN(0,1)*crLHS193 + DN(0,2)*crLHS247 + crLHS112*crLHS226 + crLHS248 + crLHS249);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS198 + DN(0,2)*crLHS250 + crLHS112*crLHS228 + crLHS251 + crLHS252);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS120 + DN(0,1)*crLHS205 + DN(0,2)*crLHS253 + crLHS102 + crLHS25*crLHS255 - crLHS255*crLHS28 + crLHS257);
rLHS(2,11)+=-gauss_weight*(DN(2,2)*crLHS25 - DN(2,2)*crLHS28 + crLHS258);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS128 + DN(0,1)*crLHS210 + DN(0,2)*crLHS259 + crLHS147*crLHS226 + crLHS260 + crLHS261);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS145 + DN(0,1)*crLHS215 + DN(0,2)*crLHS262 + crLHS147*crLHS228 + crLHS263 + crLHS264);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS155 + DN(0,1)*crLHS221 + DN(0,2)*crLHS265 + crLHS137 + crLHS25*crLHS266 - crLHS266*crLHS28 + crLHS268);
rLHS(2,15)+=-gauss_weight*(DN(3,2)*crLHS25 - DN(3,2)*crLHS28 + crLHS269);
rLHS(3,0)+=gauss_weight*(crLHS161*crLHS272 + crLHS226*crLHS274 - crLHS23*crLHS275 + crLHS270);
rLHS(3,1)+=gauss_weight*(-crLHS168*crLHS277 + crLHS228*crLHS274 + crLHS271 + crLHS276*crLHS37);
rLHS(3,2)+=gauss_weight*(crLHS172*crLHS272 - crLHS233*crLHS278 + crLHS273 + crLHS276*crLHS49);
rLHS(3,3)+=crLHS279*(crLHS165 + crLHS230 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS161*crLHS281 + crLHS226*crLHS282 - crLHS275*crLHS64 + crLHS280);
rLHS(3,5)+=gauss_weight*(-crLHS183*crLHS277 + crLHS228*crLHS282 + crLHS283 + crLHS284*crLHS37);
rLHS(3,6)+=gauss_weight*(crLHS172*crLHS281 - crLHS243*crLHS278 + crLHS284*crLHS49 + crLHS285);
rLHS(3,7)+=crLHS286;
rLHS(3,8)+=gauss_weight*(-crLHS100*crLHS275 + crLHS161*crLHS288 + crLHS226*crLHS289 + crLHS287);
rLHS(3,9)+=gauss_weight*(-crLHS200*crLHS277 + crLHS228*crLHS289 + crLHS290 + crLHS291*crLHS37);
rLHS(3,10)+=gauss_weight*(crLHS172*crLHS288 - crLHS255*crLHS278 + crLHS291*crLHS49 + crLHS292);
rLHS(3,11)+=crLHS293;
rLHS(3,12)+=gauss_weight*(-crLHS135*crLHS275 + crLHS161*crLHS295 + crLHS226*crLHS296 + crLHS294);
rLHS(3,13)+=gauss_weight*(-crLHS216*crLHS277 + crLHS228*crLHS296 + crLHS297 + crLHS298*crLHS37);
rLHS(3,14)+=gauss_weight*(crLHS172*crLHS295 - crLHS266*crLHS278 + crLHS298*crLHS49 + crLHS299);
rLHS(3,15)+=crLHS300;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS23*crLHS301 - crLHS23*crLHS303 + crLHS304 + crLHS68);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS31 + DN(1,1)*crLHS33 + DN(1,2)*crLHS36 + crLHS177 + crLHS306*crLHS37 + crLHS79);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS44 + DN(1,1)*crLHS46 + DN(1,2)*crLHS48 + crLHS236 + crLHS306*crLHS49 + crLHS86);
rLHS(4,3)+=-gauss_weight*(-DN(0,0)*crLHS303 + crLHS275*crLHS59 + crLHS280);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS52 + DN(1,1)*crLHS54 + DN(1,2)*crLHS56 + crLHS12*crLHS307 + crLHS13*crLHS309 + crLHS301*crLHS64 - crLHS303*crLHS64 + crLHS311);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS69 + DN(1,1)*crLHS71 + DN(1,2)*crLHS74 + crLHS309*crLHS37 + crLHS313 - crLHS314*crLHS37 + crLHS315*crLHS37);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS80 + DN(1,1)*crLHS82 + DN(1,2)*crLHS84 + crLHS309*crLHS49 - crLHS314*crLHS49 + crLHS315*crLHS49 + crLHS316);
rLHS(4,7)+=-DN(1,0)*crLHS317;
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS88 + DN(1,1)*crLHS90 + DN(1,2)*crLHS92 + crLHS100*crLHS301 - crLHS100*crLHS303 + crLHS319 + crLHS321);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS105 + DN(1,1)*crLHS107 + DN(1,2)*crLHS110 + crLHS322 + crLHS323*crLHS37 + crLHS326);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS116 + DN(1,1)*crLHS118 + DN(1,2)*crLHS120 + crLHS323*crLHS49 + crLHS327 + crLHS328);
rLHS(4,11)+=-gauss_weight*(DN(2,0)*crLHS301 - DN(2,0)*crLHS303 + crLHS329);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS124 + DN(1,1)*crLHS126 + DN(1,2)*crLHS128 + crLHS135*crLHS301 - crLHS135*crLHS303 + crLHS331 + crLHS333);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS140 + DN(1,1)*crLHS142 + DN(1,2)*crLHS145 + crLHS334 + crLHS335*crLHS37 + crLHS338);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS151 + DN(1,1)*crLHS153 + DN(1,2)*crLHS155 + crLHS335*crLHS49 + crLHS339 + crLHS340);
rLHS(4,15)+=-gauss_weight*(DN(3,0)*crLHS301 - DN(3,0)*crLHS303 + crLHS341);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS159 + DN(1,2)*crLHS160 + crLHS161*crLHS306 + crLHS178 + crLHS75);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS33 + DN(1,1)*crLHS162 + DN(1,2)*crLHS164 + crLHS168*crLHS301 - crLHS168*crLHS303 + crLHS185 + crLHS304);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS46 + DN(1,1)*crLHS169 + DN(1,2)*crLHS171 + crLHS172*crLHS306 + crLHS190 + crLHS239);
rLHS(5,3)+=-gauss_weight*(-DN(0,1)*crLHS303 + crLHS277*crLHS59 + crLHS283);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS54 + DN(1,1)*crLHS175 + DN(1,2)*crLHS176 + crLHS161*crLHS309 - crLHS161*crLHS314 + crLHS161*crLHS315 + crLHS313);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS71 + DN(1,1)*crLHS179 + DN(1,2)*crLHS181 + crLHS12*crLHS342 + crLHS166*crLHS309 + crLHS183*crLHS301 - crLHS183*crLHS303 + crLHS311);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS82 + DN(1,1)*crLHS186 + DN(1,2)*crLHS188 + crLHS172*crLHS309 - crLHS172*crLHS314 + crLHS172*crLHS315 + crLHS344);
rLHS(5,7)+=-DN(1,1)*crLHS317;
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS90 + DN(1,1)*crLHS192 + DN(1,2)*crLHS193 + crLHS161*crLHS323 + crLHS345 + crLHS346);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS107 + DN(1,1)*crLHS196 + DN(1,2)*crLHS198 + crLHS200*crLHS301 - crLHS200*crLHS303 + crLHS319 + crLHS348);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS118 + DN(1,1)*crLHS203 + DN(1,2)*crLHS205 + crLHS172*crLHS323 + crLHS349 + crLHS350);
rLHS(5,11)+=-gauss_weight*(DN(2,1)*crLHS301 - DN(2,1)*crLHS303 + crLHS351);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS126 + DN(1,1)*crLHS209 + DN(1,2)*crLHS210 + crLHS161*crLHS335 + crLHS352 + crLHS353);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS142 + DN(1,1)*crLHS213 + DN(1,2)*crLHS215 + crLHS216*crLHS301 - crLHS216*crLHS303 + crLHS331 + crLHS355);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS153 + DN(1,1)*crLHS219 + DN(1,2)*crLHS221 + crLHS172*crLHS335 + crLHS356 + crLHS357);
rLHS(5,15)+=-gauss_weight*(DN(3,1)*crLHS301 - DN(3,1)*crLHS303 + crLHS358);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS160 + DN(1,2)*crLHS225 + crLHS226*crLHS306 + crLHS237 + crLHS85);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS36 + DN(1,1)*crLHS164 + DN(1,2)*crLHS227 + crLHS189 + crLHS228*crLHS306 + crLHS240);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS48 + DN(1,1)*crLHS171 + DN(1,2)*crLHS229 + crLHS233*crLHS301 - crLHS233*crLHS303 + crLHS245 + crLHS304);
rLHS(6,3)+=-gauss_weight*(-DN(0,2)*crLHS303 + crLHS278*crLHS59 + crLHS285);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS56 + DN(1,1)*crLHS176 + DN(1,2)*crLHS234 + crLHS226*crLHS309 - crLHS226*crLHS314 + crLHS226*crLHS315 + crLHS316);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS181 + DN(1,2)*crLHS238 + crLHS228*crLHS309 - crLHS228*crLHS314 + crLHS228*crLHS315 + crLHS344);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS84 + DN(1,1)*crLHS188 + DN(1,2)*crLHS241 + crLHS12*crLHS359 + crLHS231*crLHS309 + crLHS243*crLHS301 - crLHS243*crLHS303 + crLHS311);
rLHS(6,7)+=-DN(1,2)*crLHS317;
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS92 + DN(1,1)*crLHS193 + DN(1,2)*crLHS247 + crLHS226*crLHS323 + crLHS361 + crLHS362);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS198 + DN(1,2)*crLHS250 + crLHS228*crLHS323 + crLHS363 + crLHS364);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS120 + DN(1,1)*crLHS205 + DN(1,2)*crLHS253 + crLHS255*crLHS301 - crLHS255*crLHS303 + crLHS319 + crLHS366);
rLHS(6,11)+=-gauss_weight*(DN(2,2)*crLHS301 - DN(2,2)*crLHS303 + crLHS367);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS128 + DN(1,1)*crLHS210 + DN(1,2)*crLHS259 + crLHS226*crLHS335 + crLHS368 + crLHS369);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS145 + DN(1,1)*crLHS215 + DN(1,2)*crLHS262 + crLHS228*crLHS335 + crLHS370 + crLHS371);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS155 + DN(1,1)*crLHS221 + DN(1,2)*crLHS265 + crLHS266*crLHS301 - crLHS266*crLHS303 + crLHS331 + crLHS373);
rLHS(6,15)+=-gauss_weight*(DN(3,2)*crLHS301 - DN(3,2)*crLHS303 + crLHS374);
rLHS(7,0)+=gauss_weight*(crLHS161*crLHS375 + crLHS226*crLHS376 - crLHS23*crLHS377 + crLHS87);
rLHS(7,1)+=gauss_weight*(-crLHS168*crLHS379 + crLHS191 + crLHS228*crLHS376 + crLHS37*crLHS378);
rLHS(7,2)+=gauss_weight*(crLHS172*crLHS375 - crLHS233*crLHS380 + crLHS246 + crLHS378*crLHS49);
rLHS(7,3)+=crLHS286;
rLHS(7,4)+=gauss_weight*(crLHS161*crLHS383 + crLHS226*crLHS385 - crLHS377*crLHS64 + crLHS381);
rLHS(7,5)+=gauss_weight*(-crLHS183*crLHS379 + crLHS228*crLHS385 + crLHS37*crLHS386 + crLHS382);
rLHS(7,6)+=gauss_weight*(crLHS172*crLHS383 - crLHS243*crLHS380 + crLHS384 + crLHS386*crLHS49);
rLHS(7,7)+=crLHS279*(crLHS307 + crLHS342 + crLHS359);
rLHS(7,8)+=gauss_weight*(-crLHS100*crLHS377 + crLHS161*crLHS388 + crLHS226*crLHS389 + crLHS387);
rLHS(7,9)+=gauss_weight*(-crLHS200*crLHS379 + crLHS228*crLHS389 + crLHS37*crLHS391 + crLHS390);
rLHS(7,10)+=gauss_weight*(crLHS172*crLHS388 - crLHS255*crLHS380 + crLHS391*crLHS49 + crLHS392);
rLHS(7,11)+=crLHS393;
rLHS(7,12)+=gauss_weight*(-crLHS135*crLHS377 + crLHS161*crLHS395 + crLHS226*crLHS396 + crLHS394);
rLHS(7,13)+=gauss_weight*(-crLHS216*crLHS379 + crLHS228*crLHS396 + crLHS37*crLHS398 + crLHS397);
rLHS(7,14)+=gauss_weight*(crLHS172*crLHS395 - crLHS266*crLHS380 + crLHS398*crLHS49 + crLHS399);
rLHS(7,15)+=crLHS400;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS104 + crLHS23*crLHS401 - crLHS23*crLHS403 + crLHS404);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS31 + DN(2,1)*crLHS33 + DN(2,2)*crLHS36 + crLHS115 + crLHS194 + crLHS37*crLHS406);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS44 + DN(2,1)*crLHS46 + DN(2,2)*crLHS48 + crLHS122 + crLHS248 + crLHS406*crLHS49);
rLHS(8,3)+=-gauss_weight*(-DN(0,0)*crLHS403 + crLHS275*crLHS95 + crLHS287);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS52 + DN(2,1)*crLHS54 + DN(2,2)*crLHS56 + crLHS321 + crLHS401*crLHS64 - crLHS403*crLHS64 + crLHS407);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS69 + DN(2,1)*crLHS71 + DN(2,2)*crLHS74 + crLHS326 + crLHS345 + crLHS37*crLHS408);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS80 + DN(2,1)*crLHS82 + DN(2,2)*crLHS84 + crLHS328 + crLHS361 + crLHS408*crLHS49);
rLHS(8,7)+=-gauss_weight*(-DN(1,0)*crLHS403 + crLHS377*crLHS95 + crLHS387);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS88 + DN(2,1)*crLHS90 + DN(2,2)*crLHS92 + crLHS100*crLHS401 - crLHS100*crLHS403 + crLHS12*crLHS409 + crLHS13*crLHS411 + crLHS413);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS105 + DN(2,1)*crLHS107 + DN(2,2)*crLHS110 + crLHS37*crLHS411 - crLHS37*crLHS416 + crLHS37*crLHS417 + crLHS415);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS116 + DN(2,1)*crLHS118 + DN(2,2)*crLHS120 + crLHS411*crLHS49 - crLHS416*crLHS49 + crLHS417*crLHS49 + crLHS418);
rLHS(8,11)+=-DN(2,0)*crLHS419;
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS124 + DN(2,1)*crLHS126 + DN(2,2)*crLHS128 + crLHS135*crLHS401 - crLHS135*crLHS403 + crLHS421 + crLHS423);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS140 + DN(2,1)*crLHS142 + DN(2,2)*crLHS145 + crLHS37*crLHS425 + crLHS424 + crLHS428);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS151 + DN(2,1)*crLHS153 + DN(2,2)*crLHS155 + crLHS425*crLHS49 + crLHS429 + crLHS430);
rLHS(8,15)+=-gauss_weight*(DN(3,0)*crLHS401 - DN(3,0)*crLHS403 + crLHS431);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS159 + DN(2,2)*crLHS160 + crLHS111 + crLHS161*crLHS406 + crLHS195);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS33 + DN(2,1)*crLHS162 + DN(2,2)*crLHS164 + crLHS168*crLHS401 - crLHS168*crLHS403 + crLHS202 + crLHS404);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS46 + DN(2,1)*crLHS169 + DN(2,2)*crLHS171 + crLHS172*crLHS406 + crLHS207 + crLHS251);
rLHS(9,3)+=-gauss_weight*(-DN(0,1)*crLHS403 + crLHS277*crLHS95 + crLHS290);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS54 + DN(2,1)*crLHS175 + DN(2,2)*crLHS176 + crLHS161*crLHS408 + crLHS322 + crLHS346);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS71 + DN(2,1)*crLHS179 + DN(2,2)*crLHS181 + crLHS183*crLHS401 - crLHS183*crLHS403 + crLHS348 + crLHS407);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS82 + DN(2,1)*crLHS186 + DN(2,2)*crLHS188 + crLHS172*crLHS408 + crLHS350 + crLHS363);
rLHS(9,7)+=-gauss_weight*(-DN(1,1)*crLHS403 + crLHS379*crLHS95 + crLHS390);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS90 + DN(2,1)*crLHS192 + DN(2,2)*crLHS193 + crLHS161*crLHS411 - crLHS161*crLHS416 + crLHS161*crLHS417 + crLHS415);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS107 + DN(2,1)*crLHS196 + DN(2,2)*crLHS198 + crLHS12*crLHS432 + crLHS166*crLHS411 + crLHS200*crLHS401 - crLHS200*crLHS403 + crLHS413);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS118 + DN(2,1)*crLHS203 + DN(2,2)*crLHS205 + crLHS172*crLHS411 - crLHS172*crLHS416 + crLHS172*crLHS417 + crLHS434);
rLHS(9,11)+=-DN(2,1)*crLHS419;
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS126 + DN(2,1)*crLHS209 + DN(2,2)*crLHS210 + crLHS161*crLHS425 + crLHS435 + crLHS436);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS142 + DN(2,1)*crLHS213 + DN(2,2)*crLHS215 + crLHS216*crLHS401 - crLHS216*crLHS403 + crLHS421 + crLHS438);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS153 + DN(2,1)*crLHS219 + DN(2,2)*crLHS221 + crLHS172*crLHS425 + crLHS439 + crLHS440);
rLHS(9,15)+=-gauss_weight*(DN(3,1)*crLHS401 - DN(3,1)*crLHS403 + crLHS441);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS160 + DN(2,2)*crLHS225 + crLHS121 + crLHS226*crLHS406 + crLHS249);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS36 + DN(2,1)*crLHS164 + DN(2,2)*crLHS227 + crLHS206 + crLHS228*crLHS406 + crLHS252);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS48 + DN(2,1)*crLHS171 + DN(2,2)*crLHS229 + crLHS233*crLHS401 - crLHS233*crLHS403 + crLHS257 + crLHS404);
rLHS(10,3)+=-gauss_weight*(-DN(0,2)*crLHS403 + crLHS278*crLHS95 + crLHS292);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS56 + DN(2,1)*crLHS176 + DN(2,2)*crLHS234 + crLHS226*crLHS408 + crLHS327 + crLHS362);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS181 + DN(2,2)*crLHS238 + crLHS228*crLHS408 + crLHS349 + crLHS364);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS84 + DN(2,1)*crLHS188 + DN(2,2)*crLHS241 + crLHS243*crLHS401 - crLHS243*crLHS403 + crLHS366 + crLHS407);
rLHS(10,7)+=-gauss_weight*(-DN(1,2)*crLHS403 + crLHS380*crLHS95 + crLHS392);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS92 + DN(2,1)*crLHS193 + DN(2,2)*crLHS247 + crLHS226*crLHS411 - crLHS226*crLHS416 + crLHS226*crLHS417 + crLHS418);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS198 + DN(2,2)*crLHS250 + crLHS228*crLHS411 - crLHS228*crLHS416 + crLHS228*crLHS417 + crLHS434);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS120 + DN(2,1)*crLHS205 + DN(2,2)*crLHS253 + crLHS12*crLHS442 + crLHS231*crLHS411 + crLHS255*crLHS401 - crLHS255*crLHS403 + crLHS413);
rLHS(10,11)+=-DN(2,2)*crLHS419;
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS128 + DN(2,1)*crLHS210 + DN(2,2)*crLHS259 + crLHS226*crLHS425 + crLHS444 + crLHS445);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS145 + DN(2,1)*crLHS215 + DN(2,2)*crLHS262 + crLHS228*crLHS425 + crLHS446 + crLHS447);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS155 + DN(2,1)*crLHS221 + DN(2,2)*crLHS265 + crLHS266*crLHS401 - crLHS266*crLHS403 + crLHS421 + crLHS449);
rLHS(10,15)+=-gauss_weight*(DN(3,2)*crLHS401 - DN(3,2)*crLHS403 + crLHS450);
rLHS(11,0)+=gauss_weight*(crLHS123 + crLHS161*crLHS451 + crLHS226*crLHS452 - crLHS23*crLHS453);
rLHS(11,1)+=gauss_weight*(-crLHS168*crLHS455 + crLHS208 + crLHS228*crLHS452 + crLHS37*crLHS454);
rLHS(11,2)+=gauss_weight*(crLHS172*crLHS451 - crLHS233*crLHS456 + crLHS258 + crLHS454*crLHS49);
rLHS(11,3)+=crLHS293;
rLHS(11,4)+=gauss_weight*(crLHS161*crLHS457 + crLHS226*crLHS458 + crLHS329 - crLHS453*crLHS64);
rLHS(11,5)+=gauss_weight*(-crLHS183*crLHS455 + crLHS228*crLHS458 + crLHS351 + crLHS37*crLHS459);
rLHS(11,6)+=gauss_weight*(crLHS172*crLHS457 - crLHS243*crLHS456 + crLHS367 + crLHS459*crLHS49);
rLHS(11,7)+=crLHS393;
rLHS(11,8)+=gauss_weight*(-crLHS100*crLHS453 + crLHS161*crLHS462 + crLHS226*crLHS464 + crLHS460);
rLHS(11,9)+=gauss_weight*(-crLHS200*crLHS455 + crLHS228*crLHS464 + crLHS37*crLHS465 + crLHS461);
rLHS(11,10)+=gauss_weight*(crLHS172*crLHS462 - crLHS255*crLHS456 + crLHS463 + crLHS465*crLHS49);
rLHS(11,11)+=crLHS279*(crLHS409 + crLHS432 + crLHS442);
rLHS(11,12)+=gauss_weight*(-crLHS135*crLHS453 + crLHS161*crLHS467 + crLHS226*crLHS468 + crLHS466);
rLHS(11,13)+=gauss_weight*(-crLHS216*crLHS455 + crLHS228*crLHS468 + crLHS37*crLHS470 + crLHS469);
rLHS(11,14)+=gauss_weight*(crLHS172*crLHS467 - crLHS266*crLHS456 + crLHS470*crLHS49 + crLHS471);
rLHS(11,15)+=crLHS472;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS139 + crLHS23*crLHS473 - crLHS23*crLHS475 + crLHS476);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS31 + DN(3,1)*crLHS33 + DN(3,2)*crLHS36 + crLHS150 + crLHS211 + crLHS37*crLHS478);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS44 + DN(3,1)*crLHS46 + DN(3,2)*crLHS48 + crLHS157 + crLHS260 + crLHS478*crLHS49);
rLHS(12,3)+=-gauss_weight*(-DN(0,0)*crLHS475 + crLHS130*crLHS275 + crLHS294);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS52 + DN(3,1)*crLHS54 + DN(3,2)*crLHS56 + crLHS333 + crLHS473*crLHS64 - crLHS475*crLHS64 + crLHS479);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS69 + DN(3,1)*crLHS71 + DN(3,2)*crLHS74 + crLHS338 + crLHS352 + crLHS37*crLHS480);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS80 + DN(3,1)*crLHS82 + DN(3,2)*crLHS84 + crLHS340 + crLHS368 + crLHS480*crLHS49);
rLHS(12,7)+=-gauss_weight*(-DN(1,0)*crLHS475 + crLHS130*crLHS377 + crLHS394);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS88 + DN(3,1)*crLHS90 + DN(3,2)*crLHS92 + crLHS100*crLHS473 - crLHS100*crLHS475 + crLHS423 + crLHS481);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS105 + DN(3,1)*crLHS107 + DN(3,2)*crLHS110 + crLHS37*crLHS482 + crLHS428 + crLHS435);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS116 + DN(3,1)*crLHS118 + DN(3,2)*crLHS120 + crLHS430 + crLHS444 + crLHS482*crLHS49);
rLHS(12,11)+=-gauss_weight*(-DN(2,0)*crLHS475 + crLHS130*crLHS453 + crLHS466);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS124 + DN(3,1)*crLHS126 + DN(3,2)*crLHS128 + crLHS12*crLHS483 + crLHS13*crLHS485 + crLHS135*crLHS473 - crLHS135*crLHS475 + crLHS487);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS140 + DN(3,1)*crLHS142 + DN(3,2)*crLHS145 + crLHS37*crLHS485 - crLHS37*crLHS490 + crLHS37*crLHS491 + crLHS489);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS151 + DN(3,1)*crLHS153 + DN(3,2)*crLHS155 + crLHS485*crLHS49 - crLHS49*crLHS490 + crLHS49*crLHS491 + crLHS492);
rLHS(12,15)+=-DN(3,0)*crLHS493;
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS159 + DN(3,2)*crLHS160 + crLHS146 + crLHS161*crLHS478 + crLHS212);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS33 + DN(3,1)*crLHS162 + DN(3,2)*crLHS164 + crLHS168*crLHS473 - crLHS168*crLHS475 + crLHS218 + crLHS476);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS46 + DN(3,1)*crLHS169 + DN(3,2)*crLHS171 + crLHS172*crLHS478 + crLHS223 + crLHS263);
rLHS(13,3)+=-gauss_weight*(-DN(0,1)*crLHS475 + crLHS130*crLHS277 + crLHS297);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS54 + DN(3,1)*crLHS175 + DN(3,2)*crLHS176 + crLHS161*crLHS480 + crLHS334 + crLHS353);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS71 + DN(3,1)*crLHS179 + DN(3,2)*crLHS181 + crLHS183*crLHS473 - crLHS183*crLHS475 + crLHS355 + crLHS479);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS82 + DN(3,1)*crLHS186 + DN(3,2)*crLHS188 + crLHS172*crLHS480 + crLHS357 + crLHS370);
rLHS(13,7)+=-gauss_weight*(-DN(1,1)*crLHS475 + crLHS130*crLHS379 + crLHS397);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS90 + DN(3,1)*crLHS192 + DN(3,2)*crLHS193 + crLHS161*crLHS482 + crLHS424 + crLHS436);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS107 + DN(3,1)*crLHS196 + DN(3,2)*crLHS198 + crLHS200*crLHS473 - crLHS200*crLHS475 + crLHS438 + crLHS481);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS118 + DN(3,1)*crLHS203 + DN(3,2)*crLHS205 + crLHS172*crLHS482 + crLHS440 + crLHS446);
rLHS(13,11)+=-gauss_weight*(-DN(2,1)*crLHS475 + crLHS130*crLHS455 + crLHS469);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS126 + DN(3,1)*crLHS209 + DN(3,2)*crLHS210 + crLHS161*crLHS485 - crLHS161*crLHS490 + crLHS161*crLHS491 + crLHS489);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS142 + DN(3,1)*crLHS213 + DN(3,2)*crLHS215 + crLHS12*crLHS494 + crLHS166*crLHS485 + crLHS216*crLHS473 - crLHS216*crLHS475 + crLHS487);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS153 + DN(3,1)*crLHS219 + DN(3,2)*crLHS221 + crLHS172*crLHS485 - crLHS172*crLHS490 + crLHS172*crLHS491 + crLHS495);
rLHS(13,15)+=-DN(3,1)*crLHS493;
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS160 + DN(3,2)*crLHS225 + crLHS156 + crLHS226*crLHS478 + crLHS261);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS36 + DN(3,1)*crLHS164 + DN(3,2)*crLHS227 + crLHS222 + crLHS228*crLHS478 + crLHS264);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS48 + DN(3,1)*crLHS171 + DN(3,2)*crLHS229 + crLHS233*crLHS473 - crLHS233*crLHS475 + crLHS268 + crLHS476);
rLHS(14,3)+=-gauss_weight*(-DN(0,2)*crLHS475 + crLHS130*crLHS278 + crLHS299);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS56 + DN(3,1)*crLHS176 + DN(3,2)*crLHS234 + crLHS226*crLHS480 + crLHS339 + crLHS369);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS74 + DN(3,1)*crLHS181 + DN(3,2)*crLHS238 + crLHS228*crLHS480 + crLHS356 + crLHS371);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS84 + DN(3,1)*crLHS188 + DN(3,2)*crLHS241 + crLHS243*crLHS473 - crLHS243*crLHS475 + crLHS373 + crLHS479);
rLHS(14,7)+=-gauss_weight*(-DN(1,2)*crLHS475 + crLHS130*crLHS380 + crLHS399);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS92 + DN(3,1)*crLHS193 + DN(3,2)*crLHS247 + crLHS226*crLHS482 + crLHS429 + crLHS445);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS110 + DN(3,1)*crLHS198 + DN(3,2)*crLHS250 + crLHS228*crLHS482 + crLHS439 + crLHS447);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS120 + DN(3,1)*crLHS205 + DN(3,2)*crLHS253 + crLHS255*crLHS473 - crLHS255*crLHS475 + crLHS449 + crLHS481);
rLHS(14,11)+=-gauss_weight*(-DN(2,2)*crLHS475 + crLHS130*crLHS456 + crLHS471);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS128 + DN(3,1)*crLHS210 + DN(3,2)*crLHS259 + crLHS226*crLHS485 - crLHS226*crLHS490 + crLHS226*crLHS491 + crLHS492);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS145 + DN(3,1)*crLHS215 + DN(3,2)*crLHS262 + crLHS228*crLHS485 - crLHS228*crLHS490 + crLHS228*crLHS491 + crLHS495);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS155 + DN(3,1)*crLHS221 + DN(3,2)*crLHS265 + crLHS12*crLHS496 + crLHS231*crLHS485 + crLHS266*crLHS473 - crLHS266*crLHS475 + crLHS487);
rLHS(14,15)+=-DN(3,2)*crLHS493;
rLHS(15,0)+=gauss_weight*(crLHS158 + crLHS161*crLHS497 + crLHS226*crLHS498 - crLHS23*crLHS499);
rLHS(15,1)+=gauss_weight*(-crLHS168*crLHS501 + crLHS224 + crLHS228*crLHS498 + crLHS37*crLHS500);
rLHS(15,2)+=gauss_weight*(crLHS172*crLHS497 - crLHS233*crLHS502 + crLHS269 + crLHS49*crLHS500);
rLHS(15,3)+=crLHS300;
rLHS(15,4)+=gauss_weight*(crLHS161*crLHS503 + crLHS226*crLHS504 + crLHS341 - crLHS499*crLHS64);
rLHS(15,5)+=gauss_weight*(-crLHS183*crLHS501 + crLHS228*crLHS504 + crLHS358 + crLHS37*crLHS505);
rLHS(15,6)+=gauss_weight*(crLHS172*crLHS503 - crLHS243*crLHS502 + crLHS374 + crLHS49*crLHS505);
rLHS(15,7)+=crLHS400;
rLHS(15,8)+=gauss_weight*(-crLHS100*crLHS499 + crLHS161*crLHS506 + crLHS226*crLHS507 + crLHS431);
rLHS(15,9)+=gauss_weight*(-crLHS200*crLHS501 + crLHS228*crLHS507 + crLHS37*crLHS508 + crLHS441);
rLHS(15,10)+=gauss_weight*(crLHS172*crLHS506 - crLHS255*crLHS502 + crLHS450 + crLHS49*crLHS508);
rLHS(15,11)+=crLHS472;
rLHS(15,12)+=gauss_weight*(-crLHS135*crLHS499 + crLHS161*crLHS511 + crLHS226*crLHS513 + crLHS509);
rLHS(15,13)+=gauss_weight*(-crLHS216*crLHS501 + crLHS228*crLHS513 + crLHS37*crLHS514 + crLHS510);
rLHS(15,14)+=gauss_weight*(crLHS172*crLHS511 - crLHS266*crLHS502 + crLHS49*crLHS514 + crLHS512);
rLHS(15,15)+=crLHS279*(crLHS483 + crLHS494 + crLHS496);

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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    // const double bdf1 = rData.bdf1;
    // const double bdf2 = rData.bdf2;

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
const double crRHS23 = crRHS18*crRHS3;
const double crRHS24 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1);
const double crRHS25 = crRHS11*crRHS24;
const double crRHS26 = crRHS23 + crRHS25;
const double crRHS27 = N[0]*rho;
const double crRHS28 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0);
const double crRHS29 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1);
const double crRHS30 = crRHS28 + crRHS29;
const double crRHS31 = crRHS2*stab_c3;
const double crRHS32 = rho*stab_c2*sqrt(pow(crRHS16, 2) + pow(crRHS5, 2));
const double crRHS33 = crRHS30*(h*(crRHS31*h + crRHS32)/stab_c1 + mu);
const double crRHS34 = crRHS5*rho;
const double crRHS35 = crRHS16*rho;
const double crRHS36 = crRHS14*functional_weights[6];
const double crRHS37 = 1.0/(crRHS31 + crRHS32/h + mu*stab_c1/pow(h, 2));
const double crRHS38 = crRHS37*(DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS14 + crRHS23*rho + crRHS25*rho - crRHS28*crRHS34 - crRHS35*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0)) + crRHS36 + crRHS4 + crRHS7);
const double crRHS39 = N[0]*crRHS2;
const double crRHS40 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS41 = crRHS2*crRHS24;
const double crRHS42 = crRHS16*crRHS6;
const double crRHS43 = crRHS13*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2]);
const double crRHS44 = N[0]*crRHS43;
const double crRHS45 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crRHS46 = 1.0*crRHS45;
const double crRHS47 = crRHS20*crRHS3;
const double crRHS48 = crRHS24*crRHS45;
const double crRHS49 = crRHS47 + crRHS48;
const double crRHS50 = crRHS43*functional_weights[6];
const double crRHS51 = crRHS37*(DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] - crRHS29*crRHS35 - crRHS34*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS40 + crRHS41 + crRHS42 + crRHS43 + crRHS47*rho + crRHS48*rho + crRHS50);
const double crRHS52 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS16);
const double crRHS53 = N[1]*rho;
const double crRHS54 = N[1]*crRHS2;
const double crRHS55 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS16);
const double crRHS56 = N[2]*rho;
const double crRHS57 = N[2]*crRHS2;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS33 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS12 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS15*functional_weights[6] + crRHS15 + crRHS17*crRHS3 + crRHS17*crRHS38 + crRHS22*(DN(0,0)*crRHS19 + DN(0,1)*crRHS21) + crRHS26*crRHS27 - crRHS38*crRHS39);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS12 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS33 + DN(0,1)*stress_adj[1] - N[0]*crRHS40 + N[0]*crRHS41 + N[0]*crRHS42 + crRHS17*crRHS24 + crRHS17*crRHS51 + crRHS22*(DN(0,0)*crRHS21 + DN(0,1)*crRHS46) + crRHS27*crRHS49 - crRHS39*crRHS51 + crRHS44*functional_weights[6] + crRHS44);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS38 + DN(0,1)*crRHS51 + N[0]*crRHS30);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS33 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS12 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS14 + N[1]*crRHS36 + N[1]*crRHS4 + N[1]*crRHS7 + crRHS22*(DN(1,0)*crRHS19 + DN(1,1)*crRHS21) + crRHS26*crRHS53 + crRHS3*crRHS52 + crRHS38*crRHS52 - crRHS38*crRHS54);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS12 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS33 + DN(1,1)*stress_adj[1] - N[1]*crRHS40 + N[1]*crRHS41 + N[1]*crRHS42 + N[1]*crRHS43 + N[1]*crRHS50 + crRHS22*(DN(1,0)*crRHS21 + DN(1,1)*crRHS46) + crRHS24*crRHS52 + crRHS49*crRHS53 + crRHS51*crRHS52 - crRHS51*crRHS54);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS38 + DN(1,1)*crRHS51 + N[1]*crRHS30);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS33 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS12 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS14 + N[2]*crRHS36 + N[2]*crRHS4 + N[2]*crRHS7 + crRHS22*(DN(2,0)*crRHS19 + DN(2,1)*crRHS21) + crRHS26*crRHS56 + crRHS3*crRHS55 + crRHS38*crRHS55 - crRHS38*crRHS57);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS12 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS33 + DN(2,1)*stress_adj[1] - N[2]*crRHS40 + N[2]*crRHS41 + N[2]*crRHS42 + N[2]*crRHS43 + N[2]*crRHS50 + crRHS22*(DN(2,0)*crRHS21 + DN(2,1)*crRHS46) + crRHS24*crRHS55 + crRHS49*crRHS56 + crRHS51*crRHS55 - crRHS51*crRHS57);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS38 + DN(2,1)*crRHS51 + N[2]*crRHS30);

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
    // const double dt = rData.DeltaTime;
    // const double bdf0 = rData.bdf0;
    // const double bdf1 = rData.bdf1;
    // const double bdf2 = rData.bdf2;

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
const double crRHS11 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crRHS12 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crRHS13 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS11 + DN(0,2)*crRHS12);
const double crRHS14 = DN(0,1)*v_ns(0,0);
const double crRHS15 = DN(1,1)*v_ns(1,0);
const double crRHS16 = DN(2,1)*v_ns(2,0);
const double crRHS17 = DN(3,1)*v_ns(3,0);
const double crRHS18 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crRHS19 = -crRHS14 - crRHS15 - crRHS16 - crRHS17 + crRHS18;
const double crRHS20 = DN(0,2)*v_ns(0,0);
const double crRHS21 = DN(1,2)*v_ns(1,0);
const double crRHS22 = DN(2,2)*v_ns(2,0);
const double crRHS23 = DN(3,2)*v_ns(3,0);
const double crRHS24 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crRHS25 = -crRHS20 - crRHS21 - crRHS22 - crRHS23 + crRHS24;
const double crRHS26 = 2.0*functional_weights[2]*mu;
const double crRHS27 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crRHS28 = 1.0*crRHS27;
const double crRHS29 = crRHS14 + crRHS15 + crRHS16 + crRHS17;
const double crRHS30 = 0.5*crRHS18 + 0.5*crRHS29;
const double crRHS31 = crRHS20 + crRHS21 + crRHS22 + crRHS23;
const double crRHS32 = crRHS24 + crRHS31;
const double crRHS33 = 0.5*DN(0,2);
const double crRHS34 = 4.0*functional_weights[1]*mu;
const double crRHS35 = crRHS27*crRHS3;
const double crRHS36 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1) + N[3]*v_adj(3,1);
const double crRHS37 = crRHS18*crRHS36;
const double crRHS38 = N[0]*v_adj(0,2) + N[1]*v_adj(1,2) + N[2]*v_adj(2,2) + N[3]*v_adj(3,2);
const double crRHS39 = crRHS24*crRHS38;
const double crRHS40 = crRHS35 + crRHS37 + crRHS39;
const double crRHS41 = N[0]*rho;
const double crRHS42 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0) + DN(3,0)*v_adj(3,0);
const double crRHS43 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1) + DN(3,1)*v_adj(3,1);
const double crRHS44 = DN(0,2)*v_adj(0,2) + DN(1,2)*v_adj(1,2) + DN(2,2)*v_adj(2,2) + DN(3,2)*v_adj(3,2);
const double crRHS45 = crRHS42 + crRHS43 + crRHS44;
const double crRHS46 = crRHS2*stab_c3;
const double crRHS47 = rho*stab_c2*sqrt(pow(crRHS11, 2) + pow(crRHS12, 2) + pow(crRHS5, 2));
const double crRHS48 = crRHS45*(h*(crRHS46*h + crRHS47)/stab_c1 + mu);
const double crRHS49 = crRHS5*rho;
const double crRHS50 = crRHS11*rho;
const double crRHS51 = crRHS12*rho;
const double crRHS52 = crRHS9*functional_weights[6];
const double crRHS53 = 1.0/(crRHS46 + crRHS47/h + mu*stab_c1/pow(h, 2));
const double crRHS54 = crRHS53*(DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 + crRHS35*rho + crRHS37*rho + crRHS39*rho + crRHS4 - crRHS42*crRHS49 - crRHS50*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS51*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0)) + crRHS52 + crRHS7 + crRHS9);
const double crRHS55 = N[0]*crRHS2;
const double crRHS56 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS57 = crRHS2*crRHS36;
const double crRHS58 = crRHS11*crRHS6;
const double crRHS59 = crRHS8*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3]);
const double crRHS60 = N[0]*crRHS59;
const double crRHS61 = DN(0,2)*v_ns(0,1);
const double crRHS62 = DN(1,2)*v_ns(1,1);
const double crRHS63 = DN(2,2)*v_ns(2,1);
const double crRHS64 = DN(3,2)*v_ns(3,1);
const double crRHS65 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crRHS66 = -crRHS61 - crRHS62 - crRHS63 - crRHS64 + crRHS65;
const double crRHS67 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crRHS68 = 1.0*crRHS67;
const double crRHS69 = crRHS61 + crRHS62 + crRHS63 + crRHS64;
const double crRHS70 = crRHS65 + crRHS69;
const double crRHS71 = crRHS29*crRHS3;
const double crRHS72 = crRHS36*crRHS67;
const double crRHS73 = crRHS38*crRHS65;
const double crRHS74 = crRHS71 + crRHS72 + crRHS73;
const double crRHS75 = crRHS59*functional_weights[6];
const double crRHS76 = crRHS53*(DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] - crRHS43*crRHS50 - crRHS49*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS51*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS56 + crRHS57 + crRHS58 + crRHS59 + crRHS71*rho + crRHS72*rho + crRHS73*rho + crRHS75);
const double crRHS77 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS78 = crRHS2*crRHS38;
const double crRHS79 = crRHS12*crRHS6;
const double crRHS80 = crRHS8*(DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3]);
const double crRHS81 = N[0]*crRHS80;
const double crRHS82 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crRHS83 = 1.0*crRHS82;
const double crRHS84 = 0.5*crRHS32;
const double crRHS85 = 0.5*crRHS70;
const double crRHS86 = crRHS3*crRHS31;
const double crRHS87 = crRHS36*crRHS69;
const double crRHS88 = crRHS38*crRHS82;
const double crRHS89 = crRHS86 + crRHS87 + crRHS88;
const double crRHS90 = crRHS80*functional_weights[6];
const double crRHS91 = crRHS53*(DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] - crRHS44*crRHS51 - crRHS49*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS50*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS77 + crRHS78 + crRHS79 + crRHS80 + crRHS86*rho + crRHS87*rho + crRHS88*rho + crRHS90);
const double crRHS92 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS11 + DN(1,2)*crRHS12);
const double crRHS93 = N[1]*rho;
const double crRHS94 = N[1]*crRHS2;
const double crRHS95 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS11 + DN(2,2)*crRHS12);
const double crRHS96 = N[2]*rho;
const double crRHS97 = N[2]*crRHS2;
const double crRHS98 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS11 + DN(3,2)*crRHS12);
const double crRHS99 = N[3]*rho;
const double crRHS100 = N[3]*crRHS2;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS48 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS10*functional_weights[6] + crRHS10 + crRHS13*crRHS3 + crRHS13*crRHS54 - crRHS26*(DN(0,1)*crRHS19 + DN(0,2)*crRHS25) + crRHS34*(DN(0,0)*crRHS28 + DN(0,1)*crRHS30 + crRHS32*crRHS33) + crRHS40*crRHS41 - crRHS54*crRHS55);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS48 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS56 + N[0]*crRHS57 + N[0]*crRHS58 + crRHS13*crRHS36 + crRHS13*crRHS76 + crRHS26*(DN(0,0)*crRHS19 - DN(0,2)*crRHS66) + crRHS34*(DN(0,0)*crRHS30 + DN(0,1)*crRHS68 + crRHS33*crRHS70) + crRHS41*crRHS74 - crRHS55*crRHS76 + crRHS60*functional_weights[6] + crRHS60);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS48 + DN(0,2)*stress_adj[2] - N[0]*crRHS77 + N[0]*crRHS78 + N[0]*crRHS79 + crRHS13*crRHS38 + crRHS13*crRHS91 + crRHS26*(DN(0,0)*crRHS25 + DN(0,1)*crRHS66) + crRHS34*(DN(0,0)*crRHS84 + DN(0,1)*crRHS85 + DN(0,2)*crRHS83) + crRHS41*crRHS89 - crRHS55*crRHS91 + crRHS81*functional_weights[6] + crRHS81);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS54 + DN(0,1)*crRHS76 + DN(0,2)*crRHS91 + N[0]*crRHS45);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS48 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS4 + N[1]*crRHS52 + N[1]*crRHS7 + N[1]*crRHS9 - crRHS26*(DN(1,1)*crRHS19 + DN(1,2)*crRHS25) + crRHS3*crRHS92 + crRHS34*(DN(1,0)*crRHS28 + DN(1,1)*crRHS30 + DN(1,2)*crRHS84) + crRHS40*crRHS93 + crRHS54*crRHS92 - crRHS54*crRHS94);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS48 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS56 + N[1]*crRHS57 + N[1]*crRHS58 + N[1]*crRHS59 + N[1]*crRHS75 + crRHS26*(DN(1,0)*crRHS19 - DN(1,2)*crRHS66) + crRHS34*(DN(1,0)*crRHS30 + DN(1,1)*crRHS68 + DN(1,2)*crRHS85) + crRHS36*crRHS92 + crRHS74*crRHS93 + crRHS76*crRHS92 - crRHS76*crRHS94);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS48 + DN(1,2)*stress_adj[2] - N[1]*crRHS77 + N[1]*crRHS78 + N[1]*crRHS79 + N[1]*crRHS80 + N[1]*crRHS90 + crRHS26*(DN(1,0)*crRHS25 + DN(1,1)*crRHS66) + crRHS34*(DN(1,0)*crRHS84 + DN(1,1)*crRHS85 + DN(1,2)*crRHS83) + crRHS38*crRHS92 + crRHS89*crRHS93 + crRHS91*crRHS92 - crRHS91*crRHS94);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS54 + DN(1,1)*crRHS76 + DN(1,2)*crRHS91 + N[1]*crRHS45);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS48 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS4 + N[2]*crRHS52 + N[2]*crRHS7 + N[2]*crRHS9 - crRHS26*(DN(2,1)*crRHS19 + DN(2,2)*crRHS25) + crRHS3*crRHS95 + crRHS34*(DN(2,0)*crRHS28 + DN(2,1)*crRHS30 + DN(2,2)*crRHS84) + crRHS40*crRHS96 + crRHS54*crRHS95 - crRHS54*crRHS97);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS48 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS56 + N[2]*crRHS57 + N[2]*crRHS58 + N[2]*crRHS59 + N[2]*crRHS75 + crRHS26*(DN(2,0)*crRHS19 - DN(2,2)*crRHS66) + crRHS34*(DN(2,0)*crRHS30 + DN(2,1)*crRHS68 + DN(2,2)*crRHS85) + crRHS36*crRHS95 + crRHS74*crRHS96 + crRHS76*crRHS95 - crRHS76*crRHS97);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS48 + DN(2,2)*stress_adj[2] - N[2]*crRHS77 + N[2]*crRHS78 + N[2]*crRHS79 + N[2]*crRHS80 + N[2]*crRHS90 + crRHS26*(DN(2,0)*crRHS25 + DN(2,1)*crRHS66) + crRHS34*(DN(2,0)*crRHS84 + DN(2,1)*crRHS85 + DN(2,2)*crRHS83) + crRHS38*crRHS95 + crRHS89*crRHS96 + crRHS91*crRHS95 - crRHS91*crRHS97);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS54 + DN(2,1)*crRHS76 + DN(2,2)*crRHS91 + N[2]*crRHS45);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS48 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS4 + N[3]*crRHS52 + N[3]*crRHS7 + N[3]*crRHS9 - crRHS100*crRHS54 - crRHS26*(DN(3,1)*crRHS19 + DN(3,2)*crRHS25) + crRHS3*crRHS98 + crRHS34*(DN(3,0)*crRHS28 + DN(3,1)*crRHS30 + DN(3,2)*crRHS84) + crRHS40*crRHS99 + crRHS54*crRHS98);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS48 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS56 + N[3]*crRHS57 + N[3]*crRHS58 + N[3]*crRHS59 + N[3]*crRHS75 - crRHS100*crRHS76 + crRHS26*(DN(3,0)*crRHS19 - DN(3,2)*crRHS66) + crRHS34*(DN(3,0)*crRHS30 + DN(3,1)*crRHS68 + DN(3,2)*crRHS85) + crRHS36*crRHS98 + crRHS74*crRHS99 + crRHS76*crRHS98);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS48 + DN(3,2)*stress_adj[2] - N[3]*crRHS77 + N[3]*crRHS78 + N[3]*crRHS79 + N[3]*crRHS80 + N[3]*crRHS90 - crRHS100*crRHS91 + crRHS26*(DN(3,0)*crRHS25 + DN(3,1)*crRHS66) + crRHS34*(DN(3,0)*crRHS84 + DN(3,1)*crRHS85 + DN(3,2)*crRHS83) + crRHS38*crRHS98 + crRHS89*crRHS99 + crRHS91*crRHS98);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS54 + DN(3,1)*crRHS76 + DN(3,2)*crRHS91 + N[3]*crRHS45);

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