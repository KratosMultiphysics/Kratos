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
const double crLHS3 = DN(0,0)*DN(0,0);
const double crLHS4 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crLHS5 = crLHS4*stab_c3;
const double crLHS6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crLHS7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crLHS8 = rho*stab_c2*sqrt(crLHS6*crLHS6 + crLHS7*crLHS7);
const double crLHS9 = h*(crLHS5*h + crLHS8)*1.0/stab_c1 + mu;
const double crLHS10 = DN(0,0)*crLHS6;
const double crLHS11 = DN(0,1)*crLHS7;
const double crLHS12 = crLHS10 + crLHS11;
const double crLHS13 = N[0]*rho;
const double crLHS14 = N[0]*crLHS4;
const double crLHS15 = crLHS10*rho;
const double crLHS16 = crLHS11*rho;
const double crLHS17 = crLHS14 + crLHS15 + crLHS16;
const double crLHS18 = 1.0/(crLHS5 + crLHS8*1.0/h + mu*stab_c1*1.0/(h*h));
const double crLHS19 = 1.0*crLHS18;
const double crLHS20 = crLHS19*rho;
const double crLHS21 = crLHS12*crLHS20;
const double crLHS22 = 1.0*crLHS14;
const double crLHS23 = crLHS18*crLHS22;
const double crLHS24 = crLHS17*crLHS19;
const double crLHS25 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crLHS26 = crLHS13*crLHS25;
const double crLHS27 = crLHS12*crLHS13 + crLHS17*crLHS21 - crLHS17*crLHS23 + crLHS24*crLHS26 + crLHS4*(N[0]*N[0]);
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
const double crLHS80 = DN(0,1)*DN(0,1);
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
const double crLHS106 = DN(1,0)*DN(1,0);
const double crLHS107 = crLHS100*crLHS43 + crLHS101*crLHS47 - crLHS103*crLHS47 + crLHS104*crLHS48 + crLHS4*(N[1]*N[1]);
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
const double crLHS118 = DN(1,1)*DN(1,1);
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
const double crLHS134 = DN(2,0)*DN(2,0);
const double crLHS135 = crLHS127*crLHS65 + crLHS128*crLHS69 - crLHS130*crLHS69 + crLHS131*crLHS70 + crLHS4*(N[2]*N[2]);
const double crLHS136 = DN(2,1)*crLHS119;
const double crLHS137 = gauss_weight*(-N[2] + crLHS127*crLHS33 + crLHS128 - crLHS130);
const double crLHS138 = DN(2,1)*DN(2,1);
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
    const BoundedMatrix<double,4,4>& C = rData.C;

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
const double crLHS5 = DN(0,0)*DN(0,0);
const double crLHS6 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crLHS7 = crLHS6*stab_c3;
const double crLHS8 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crLHS9 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crLHS10 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crLHS11 = rho*stab_c2*sqrt(crLHS10*crLHS10 + crLHS8*crLHS8 + crLHS9*crLHS9);
const double crLHS12 = h*(crLHS11 + crLHS7*h)*1.0/stab_c1 + mu;
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
const double crLHS23 = 1.0/(crLHS11*1.0/h + crLHS7 + mu*stab_c1*1.0/(h*h));
const double crLHS24 = 1.0*crLHS23;
const double crLHS25 = crLHS24*rho;
const double crLHS26 = crLHS16*crLHS25;
const double crLHS27 = 1.0*crLHS18;
const double crLHS28 = crLHS23*crLHS27;
const double crLHS29 = crLHS22*crLHS24;
const double crLHS30 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crLHS31 = crLHS17*crLHS30;
const double crLHS32 = crLHS16*crLHS17 + crLHS22*crLHS26 - crLHS22*crLHS28 + crLHS29*crLHS31 + crLHS6*(N[0]*N[0]);
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
const double crLHS157 = DN(0,1)*DN(0,1);
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
const double crLHS211 = DN(0,2)*DN(0,2);
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
const double crLHS251 = DN(1,0)*DN(1,0);
const double crLHS252 = crLHS245*crLHS60 + crLHS246*crLHS65 - crLHS248*crLHS65 + crLHS249*crLHS66 + crLHS6*(N[1]*N[1]);
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
const double crLHS274 = DN(1,1)*DN(1,1);
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
const double crLHS291 = DN(1,2)*DN(1,2);
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
const double crLHS315 = DN(2,0)*DN(2,0);
const double crLHS316 = crLHS100*crLHS309 - crLHS100*crLHS311 + crLHS101*crLHS312 + crLHS308*crLHS95 + crLHS6*(N[2]*N[2]);
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
const double crLHS331 = DN(2,1)*DN(2,1);
const double crLHS332 = DN(2,1)*crLHS12;
const double crLHS333 = DN(2,2)*crLHS332;
const double crLHS334 = DN(3,0)*crLHS332;
const double crLHS335 = DN(2,1)*DN(3,1);
const double crLHS336 = crLHS12*crLHS335;
const double crLHS337 = crLHS322 + crLHS324;
const double crLHS338 = DN(3,2)*crLHS332;
const double crLHS339 = DN(2,1)*N[3];
const double crLHS340 = DN(3,1)*N[2];
const double crLHS341 = DN(2,2)*DN(2,2);
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
const double crLHS359 = DN(3,0)*DN(3,0);
const double crLHS360 = crLHS129*crLHS351 + crLHS134*crLHS352 - crLHS134*crLHS354 + crLHS135*crLHS355 + crLHS6*(N[3]*N[3]);
const double crLHS361 = DN(3,0)*crLHS12;
const double crLHS362 = DN(3,1)*crLHS361;
const double crLHS363 = DN(3,2)*crLHS361;
const double crLHS364 = gauss_weight*(-N[3] + crLHS351*crLHS47 + crLHS352 - crLHS354);
const double crLHS365 = crLHS124 + crLHS356;
const double crLHS366 = crLHS266 + crLHS357;
const double crLHS367 = crLHS322 + crLHS358;
const double crLHS368 = DN(3,1)*DN(3,1);
const double crLHS369 = DN(3,1)*DN(3,2)*crLHS12;
const double crLHS370 = DN(3,2)*DN(3,2);
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
const double crRHS14 = rho*stab_c2*sqrt(crRHS5*crRHS5 + crRHS7*crRHS7);
const double crRHS15 = crRHS12*(h*(crRHS13*h + crRHS14)*1.0/stab_c1 + mu);
const double crRHS16 = 1.0*1.0/(crRHS13 + crRHS14*1.0/h + mu*stab_c1*1.0/(h*h));
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
const double crRHS17 = rho*stab_c2*sqrt(crRHS5*crRHS5 + crRHS7*crRHS7 + crRHS9*crRHS9);
const double crRHS18 = crRHS15*(h*(crRHS16*h + crRHS17)*1.0/stab_c1 + mu);
const double crRHS19 = 1.0*1.0/(crRHS16 + crRHS17*1.0/h + mu*stab_c1*1.0/(h*h));
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

    const BoundedMatrix<double,2,3> v_ns = rData.Velocity; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double crLHS1 = C(0,2)*DN(0,0);
const double crLHS2 = C(2,2)*DN(0,1) + crLHS1;
const double crLHS3 = DN(0,0)*DN(0,0);
const double crLHS4 = h*h;
const double crLHS5 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2];
const double crLHS6 = crLHS5*stab_c3;
const double crLHS7 = crLHS4*crLHS6*1.0/stab_c1 + mu;
const double crLHS8 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crLHS9 = N[0]*N[0];
const double crLHS10 = crLHS9*rho;
const double crLHS11 = N[0]*crLHS5;
const double crLHS12 = N[0]*rho;
const double crLHS13 = crLHS12*crLHS8;
const double crLHS14 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0);
const double crLHS15 = DN(0,0)*crLHS14;
const double crLHS16 = crLHS15*rho;
const double crLHS17 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1);
const double crLHS18 = DN(0,1)*crLHS17;
const double crLHS19 = crLHS18*rho;
const double crLHS20 = -crLHS11 + crLHS16 + crLHS19;
const double crLHS21 = -crLHS13 + crLHS20;
const double crLHS22 = 1.0/(crLHS6 + mu*stab_c1*1.0/crLHS4);
const double crLHS23 = 1.0*crLHS22;
const double crLHS24 = crLHS21*crLHS23;
const double crLHS25 = crLHS15 + crLHS18;
const double crLHS26 = crLHS25*rho;
const double crLHS27 = DN(0,0)*v_ns(0,1);
const double crLHS28 = DN(1,0)*v_ns(1,1);
const double crLHS29 = DN(2,0)*v_ns(2,1);
const double crLHS30 = crLHS27 + crLHS28 + crLHS29;
const double crLHS31 = DN(0,1)*v_ns(0,0);
const double crLHS32 = DN(1,1)*v_ns(1,0);
const double crLHS33 = DN(2,1)*v_ns(2,0);
const double crLHS34 = crLHS31 + crLHS32 + crLHS33;
const double crLHS35 = crLHS30*crLHS34;
const double crLHS36 = crLHS12*crLHS35;
const double crLHS37 = -crLHS21*crLHS8 + crLHS36;
const double crLHS38 = crLHS12*crLHS23;
const double crLHS39 = crLHS5*crLHS9;
const double crLHS40 = crLHS12*crLHS25 + crLHS39;
const double crLHS41 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS42 = C(1,2)*DN(0,1);
const double crLHS43 = C(2,2)*DN(0,0) + crLHS42;
const double crLHS44 = DN(0,0)*crLHS7;
const double crLHS45 = DN(0,1)*crLHS44;
const double crLHS46 = crLHS23*crLHS30;
const double crLHS47 = crLHS39*rho;
const double crLHS48 = rho*rho;
const double crLHS49 = crLHS25*crLHS48;
const double crLHS50 = crLHS46*crLHS49;
const double crLHS51 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crLHS52 = crLHS12*crLHS51;
const double crLHS53 = crLHS11 + crLHS13 - crLHS16 - crLHS19 + crLHS52;
const double crLHS54 = 1.0*crLHS27 + 1.0*crLHS28 + 1.0*crLHS29;
const double crLHS55 = crLHS12*crLHS22;
const double crLHS56 = crLHS54*crLHS55;
const double crLHS57 = DN(0,0)*N[0];
const double crLHS58 = DN(0,0)*crLHS23;
const double crLHS59 = DN(0,0)*crLHS8 + DN(0,1)*crLHS30;
const double crLHS60 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS61 = C(0,2)*DN(1,0);
const double crLHS62 = C(2,2)*DN(1,1) + crLHS61;
const double crLHS63 = N[1]*rho;
const double crLHS64 = crLHS63*crLHS8;
const double crLHS65 = N[1]*crLHS5;
const double crLHS66 = DN(1,0)*crLHS14;
const double crLHS67 = crLHS66*rho;
const double crLHS68 = DN(1,1)*crLHS17;
const double crLHS69 = crLHS68*rho;
const double crLHS70 = -crLHS65 + crLHS67 + crLHS69;
const double crLHS71 = -crLHS64 + crLHS70;
const double crLHS72 = crLHS23*crLHS71;
const double crLHS73 = crLHS35*crLHS63;
const double crLHS74 = -crLHS71*crLHS8 + crLHS73;
const double crLHS75 = N[1]*crLHS11;
const double crLHS76 = crLHS25*crLHS63 + crLHS75;
const double crLHS77 = DN(0,0)*DN(1,0);
const double crLHS78 = N[1]*crLHS13 + crLHS7*crLHS77;
const double crLHS79 = C(0,1)*DN(1,1) + crLHS61;
const double crLHS80 = C(1,2)*DN(1,1);
const double crLHS81 = C(2,2)*DN(1,0) + crLHS80;
const double crLHS82 = DN(1,1)*crLHS44;
const double crLHS83 = crLHS51*crLHS63;
const double crLHS84 = crLHS64 + crLHS65 - crLHS67 - crLHS69 + crLHS83;
const double crLHS85 = crLHS12*crLHS30;
const double crLHS86 = crLHS46*rho;
const double crLHS87 = N[1]*crLHS85 - crLHS75*crLHS86;
const double crLHS88 = DN(0,0)*N[1];
const double crLHS89 = DN(1,0)*crLHS23;
const double crLHS90 = DN(1,0)*crLHS8 + DN(1,1)*crLHS30;
const double crLHS91 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS92 = C(0,2)*DN(2,0);
const double crLHS93 = C(2,2)*DN(2,1) + crLHS92;
const double crLHS94 = N[2]*rho;
const double crLHS95 = crLHS8*crLHS94;
const double crLHS96 = N[2]*crLHS5;
const double crLHS97 = DN(2,0)*crLHS14;
const double crLHS98 = crLHS97*rho;
const double crLHS99 = DN(2,1)*crLHS17;
const double crLHS100 = crLHS99*rho;
const double crLHS101 = crLHS100 - crLHS96 + crLHS98;
const double crLHS102 = crLHS101 - crLHS95;
const double crLHS103 = crLHS102*crLHS23;
const double crLHS104 = crLHS35*crLHS94;
const double crLHS105 = -crLHS102*crLHS8 + crLHS104;
const double crLHS106 = N[2]*crLHS11;
const double crLHS107 = crLHS106 + crLHS25*crLHS94;
const double crLHS108 = DN(0,0)*DN(2,0);
const double crLHS109 = N[2]*crLHS13 + crLHS108*crLHS7;
const double crLHS110 = C(0,1)*DN(2,1) + crLHS92;
const double crLHS111 = C(1,2)*DN(2,1);
const double crLHS112 = C(2,2)*DN(2,0) + crLHS111;
const double crLHS113 = DN(2,1)*crLHS44;
const double crLHS114 = crLHS51*crLHS94;
const double crLHS115 = -crLHS100 + crLHS114 + crLHS95 + crLHS96 - crLHS98;
const double crLHS116 = N[2]*crLHS85 - crLHS106*crLHS86;
const double crLHS117 = DN(0,0)*N[2];
const double crLHS118 = DN(2,0)*crLHS23;
const double crLHS119 = DN(2,0)*crLHS8 + DN(2,1)*crLHS30;
const double crLHS120 = C(0,1)*DN(0,0) + crLHS42;
const double crLHS121 = crLHS23*crLHS34;
const double crLHS122 = crLHS121*crLHS49;
const double crLHS123 = 1.0*crLHS31 + 1.0*crLHS32 + 1.0*crLHS33;
const double crLHS124 = crLHS123*crLHS55;
const double crLHS125 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS126 = DN(0,1)*DN(0,1);
const double crLHS127 = crLHS20 - crLHS52;
const double crLHS128 = crLHS127*crLHS23;
const double crLHS129 = -crLHS127*crLHS51 + crLHS36;
const double crLHS130 = DN(0,1)*N[0];
const double crLHS131 = DN(0,1)*crLHS23;
const double crLHS132 = DN(0,0)*crLHS34 + DN(0,1)*crLHS51;
const double crLHS133 = C(0,1)*DN(1,0) + crLHS80;
const double crLHS134 = DN(0,1)*crLHS7;
const double crLHS135 = DN(1,0)*crLHS134;
const double crLHS136 = crLHS12*crLHS34;
const double crLHS137 = crLHS121*rho;
const double crLHS138 = N[1]*crLHS136 - crLHS137*crLHS75;
const double crLHS139 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS140 = crLHS70 - crLHS83;
const double crLHS141 = crLHS140*crLHS23;
const double crLHS142 = -crLHS140*crLHS51 + crLHS73;
const double crLHS143 = DN(0,1)*DN(1,1);
const double crLHS144 = N[1]*crLHS52 + crLHS143*crLHS7;
const double crLHS145 = DN(0,1)*N[1];
const double crLHS146 = DN(1,1)*crLHS23;
const double crLHS147 = DN(1,0)*crLHS34 + DN(1,1)*crLHS51;
const double crLHS148 = C(0,1)*DN(2,0) + crLHS111;
const double crLHS149 = DN(2,0)*crLHS134;
const double crLHS150 = N[2]*crLHS136 - crLHS106*crLHS137;
const double crLHS151 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS152 = crLHS101 - crLHS114;
const double crLHS153 = crLHS152*crLHS23;
const double crLHS154 = crLHS104 - crLHS152*crLHS51;
const double crLHS155 = DN(0,1)*DN(2,1);
const double crLHS156 = N[2]*crLHS52 + crLHS155*crLHS7;
const double crLHS157 = DN(0,1)*N[2];
const double crLHS158 = DN(2,1)*crLHS23;
const double crLHS159 = DN(2,0)*crLHS34 + DN(2,1)*crLHS51;
const double crLHS160 = crLHS23*gauss_weight;
const double crLHS161 = DN(1,0)*N[0];
const double crLHS162 = DN(1,1)*N[0];
const double crLHS163 = crLHS160*(crLHS143 + crLHS77);
const double crLHS164 = DN(2,0)*N[0];
const double crLHS165 = DN(2,1)*N[0];
const double crLHS166 = crLHS160*(crLHS108 + crLHS155);
const double crLHS167 = crLHS66 + crLHS68;
const double crLHS168 = crLHS167*rho;
const double crLHS169 = crLHS23*crLHS63;
const double crLHS170 = crLHS12*crLHS167 + crLHS75;
const double crLHS171 = crLHS22*crLHS63;
const double crLHS172 = crLHS171*crLHS54;
const double crLHS173 = crLHS167*crLHS48;
const double crLHS174 = crLHS173*crLHS46;
const double crLHS175 = DN(1,0)*DN(1,0);
const double crLHS176 = N[1]*N[1];
const double crLHS177 = crLHS176*rho;
const double crLHS178 = crLHS176*crLHS5;
const double crLHS179 = crLHS167*crLHS63 + crLHS178;
const double crLHS180 = DN(1,0)*crLHS7;
const double crLHS181 = DN(1,1)*crLHS180;
const double crLHS182 = DN(1,0)*N[1];
const double crLHS183 = N[2]*crLHS65;
const double crLHS184 = crLHS167*crLHS94 + crLHS183;
const double crLHS185 = DN(1,0)*DN(2,0);
const double crLHS186 = N[2]*crLHS64 + crLHS185*crLHS7;
const double crLHS187 = DN(2,1)*crLHS180;
const double crLHS188 = N[2]*crLHS63;
const double crLHS189 = crLHS23*crLHS94;
const double crLHS190 = crLHS189*crLHS65;
const double crLHS191 = crLHS188*crLHS30 - crLHS190*crLHS30;
const double crLHS192 = DN(1,0)*N[2];
const double crLHS193 = crLHS123*crLHS171;
const double crLHS194 = crLHS121*crLHS173;
const double crLHS195 = DN(1,1)*DN(1,1);
const double crLHS196 = DN(1,1)*N[1];
const double crLHS197 = DN(2,0)*crLHS7;
const double crLHS198 = DN(1,1)*crLHS197;
const double crLHS199 = crLHS188*crLHS34 - crLHS190*crLHS34;
const double crLHS200 = DN(1,1)*DN(2,1);
const double crLHS201 = N[2]*crLHS83 + crLHS200*crLHS7;
const double crLHS202 = DN(1,1)*N[2];
const double crLHS203 = DN(2,0)*N[1];
const double crLHS204 = DN(2,1)*N[1];
const double crLHS205 = crLHS160*(crLHS185 + crLHS200);
const double crLHS206 = crLHS97 + crLHS99;
const double crLHS207 = crLHS206*rho;
const double crLHS208 = crLHS106 + crLHS12*crLHS206;
const double crLHS209 = crLHS22*crLHS94;
const double crLHS210 = crLHS209*crLHS54;
const double crLHS211 = crLHS206*crLHS48;
const double crLHS212 = crLHS211*crLHS46;
const double crLHS213 = crLHS183 + crLHS206*crLHS63;
const double crLHS214 = DN(2,0)*DN(2,0);
const double crLHS215 = N[2]*N[2];
const double crLHS216 = crLHS215*rho;
const double crLHS217 = crLHS215*crLHS5;
const double crLHS218 = crLHS206*crLHS94 + crLHS217;
const double crLHS219 = DN(2,1)*crLHS197;
const double crLHS220 = DN(2,0)*N[2];
const double crLHS221 = crLHS123*crLHS209;
const double crLHS222 = crLHS121*crLHS211;
const double crLHS223 = DN(2,1)*DN(2,1);
const double crLHS224 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS10*crLHS8 + crLHS11*crLHS24 + crLHS24*crLHS26 + crLHS3*crLHS7 - crLHS37*crLHS38 + crLHS40);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS41 + DN(0,1)*crLHS43 - N[0]*crLHS50 + crLHS10*crLHS30 + crLHS45 - crLHS46*crLHS47 - crLHS53*crLHS56);
rLHS(0,2)+=-gauss_weight*(crLHS11*crLHS58 + crLHS26*crLHS58 + crLHS38*crLHS59 + crLHS57);
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS60 + DN(0,1)*crLHS62 + crLHS11*crLHS72 + crLHS26*crLHS72 - crLHS38*crLHS74 + crLHS76 + crLHS78);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS81 - N[1]*crLHS50 - crLHS56*crLHS84 + crLHS82 + crLHS87);
rLHS(0,5)+=-gauss_weight*(crLHS11*crLHS89 + crLHS26*crLHS89 + crLHS38*crLHS90 + crLHS88);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS91 + DN(0,1)*crLHS93 + crLHS103*crLHS11 + crLHS103*crLHS26 - crLHS105*crLHS38 + crLHS107 + crLHS109);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS112 - N[2]*crLHS50 + crLHS113 - crLHS115*crLHS56 + crLHS116);
rLHS(0,8)+=-gauss_weight*(crLHS11*crLHS118 + crLHS117 + crLHS118*crLHS26 + crLHS119*crLHS38);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS120 - N[0]*crLHS122 + crLHS10*crLHS34 - crLHS121*crLHS47 - crLHS124*crLHS53 + crLHS45);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS125 + crLHS10*crLHS51 + crLHS11*crLHS128 + crLHS126*crLHS7 + crLHS128*crLHS26 - crLHS129*crLHS38 + crLHS40);
rLHS(1,2)+=-gauss_weight*(crLHS11*crLHS131 + crLHS130 + crLHS131*crLHS26 + crLHS132*crLHS38);
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS133 - N[1]*crLHS122 - crLHS124*crLHS84 + crLHS135 + crLHS138);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS81 + DN(0,1)*crLHS139 + crLHS11*crLHS141 + crLHS141*crLHS26 - crLHS142*crLHS38 + crLHS144 + crLHS76);
rLHS(1,5)+=-gauss_weight*(crLHS11*crLHS146 + crLHS145 + crLHS146*crLHS26 + crLHS147*crLHS38);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS93 + DN(0,1)*crLHS148 - N[2]*crLHS122 - crLHS115*crLHS124 + crLHS149 + crLHS150);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS151 + crLHS107 + crLHS11*crLHS153 + crLHS153*crLHS26 - crLHS154*crLHS38 + crLHS156);
rLHS(1,8)+=-gauss_weight*(crLHS11*crLHS158 + crLHS157 + crLHS158*crLHS26 + crLHS159*crLHS38);
rLHS(2,0)+=gauss_weight*(crLHS130*crLHS137 - crLHS21*crLHS58 + crLHS57);
rLHS(2,1)+=gauss_weight*(-crLHS127*crLHS131 + crLHS130 + crLHS57*crLHS86);
rLHS(2,2)+=crLHS160*(crLHS126 + crLHS3);
rLHS(2,3)+=gauss_weight*(crLHS137*crLHS145 + crLHS161 - crLHS58*crLHS71);
rLHS(2,4)+=gauss_weight*(-crLHS131*crLHS140 + crLHS162 + crLHS86*crLHS88);
rLHS(2,5)+=crLHS163;
rLHS(2,6)+=gauss_weight*(-crLHS102*crLHS58 + crLHS137*crLHS157 + crLHS164);
rLHS(2,7)+=gauss_weight*(crLHS117*crLHS86 - crLHS131*crLHS152 + crLHS165);
rLHS(2,8)+=crLHS166;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS168*crLHS24 - crLHS169*crLHS37 + crLHS170 + crLHS24*crLHS65 + crLHS78);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS41 + DN(1,1)*crLHS43 - N[0]*crLHS174 + crLHS135 - crLHS172*crLHS53 + crLHS87);
rLHS(3,2)+=-gauss_weight*(crLHS161 + crLHS168*crLHS58 + crLHS169*crLHS59 + crLHS58*crLHS65);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS60 + DN(1,1)*crLHS62 + crLHS168*crLHS72 - crLHS169*crLHS74 + crLHS175*crLHS7 + crLHS177*crLHS8 + crLHS179 + crLHS65*crLHS72);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS81 - N[1]*crLHS174 - crLHS172*crLHS84 + crLHS177*crLHS30 - crLHS178*crLHS86 + crLHS181);
rLHS(3,5)+=-gauss_weight*(crLHS168*crLHS89 + crLHS169*crLHS90 + crLHS182 + crLHS65*crLHS89);
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS91 + DN(1,1)*crLHS93 + crLHS103*crLHS168 + crLHS103*crLHS65 - crLHS105*crLHS169 + crLHS184 + crLHS186);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS112 - N[2]*crLHS174 - crLHS115*crLHS172 + crLHS187 + crLHS191);
rLHS(3,8)+=-gauss_weight*(crLHS118*crLHS168 + crLHS118*crLHS65 + crLHS119*crLHS169 + crLHS192);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS120 - N[0]*crLHS194 + crLHS138 - crLHS193*crLHS53 + crLHS82);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS125 + crLHS128*crLHS168 + crLHS128*crLHS65 - crLHS129*crLHS169 + crLHS144 + crLHS170);
rLHS(4,2)+=-gauss_weight*(crLHS131*crLHS168 + crLHS131*crLHS65 + crLHS132*crLHS169 + crLHS162);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS133 - N[1]*crLHS194 - crLHS137*crLHS178 + crLHS177*crLHS34 + crLHS181 - crLHS193*crLHS84);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS81 + DN(1,1)*crLHS139 + crLHS141*crLHS168 + crLHS141*crLHS65 - crLHS142*crLHS169 + crLHS177*crLHS51 + crLHS179 + crLHS195*crLHS7);
rLHS(4,5)+=-gauss_weight*(crLHS146*crLHS168 + crLHS146*crLHS65 + crLHS147*crLHS169 + crLHS196);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS93 + DN(1,1)*crLHS148 - N[2]*crLHS194 - crLHS115*crLHS193 + crLHS198 + crLHS199);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS151 + crLHS153*crLHS168 + crLHS153*crLHS65 - crLHS154*crLHS169 + crLHS184 + crLHS201);
rLHS(4,8)+=-gauss_weight*(crLHS158*crLHS168 + crLHS158*crLHS65 + crLHS159*crLHS169 + crLHS202);
rLHS(5,0)+=gauss_weight*(crLHS137*crLHS162 - crLHS21*crLHS89 + crLHS88);
rLHS(5,1)+=gauss_weight*(-crLHS127*crLHS146 + crLHS145 + crLHS161*crLHS86);
rLHS(5,2)+=crLHS163;
rLHS(5,3)+=gauss_weight*(crLHS137*crLHS196 + crLHS182 - crLHS71*crLHS89);
rLHS(5,4)+=gauss_weight*(-crLHS140*crLHS146 + crLHS182*crLHS86 + crLHS196);
rLHS(5,5)+=crLHS160*(crLHS175 + crLHS195);
rLHS(5,6)+=gauss_weight*(-crLHS102*crLHS89 + crLHS137*crLHS202 + crLHS203);
rLHS(5,7)+=gauss_weight*(-crLHS146*crLHS152 + crLHS192*crLHS86 + crLHS204);
rLHS(5,8)+=crLHS205;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS109 - crLHS189*crLHS37 + crLHS207*crLHS24 + crLHS208 + crLHS24*crLHS96);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS41 + DN(2,1)*crLHS43 - N[0]*crLHS212 + crLHS116 + crLHS149 - crLHS210*crLHS53);
rLHS(6,2)+=-gauss_weight*(crLHS164 + crLHS189*crLHS59 + crLHS207*crLHS58 + crLHS58*crLHS96);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS60 + DN(2,1)*crLHS62 + crLHS186 - crLHS189*crLHS74 + crLHS207*crLHS72 + crLHS213 + crLHS72*crLHS96);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS81 - N[1]*crLHS212 + crLHS191 + crLHS198 - crLHS210*crLHS84);
rLHS(6,5)+=-gauss_weight*(crLHS189*crLHS90 + crLHS203 + crLHS207*crLHS89 + crLHS89*crLHS96);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS91 + DN(2,1)*crLHS93 + crLHS103*crLHS207 + crLHS103*crLHS96 - crLHS105*crLHS189 + crLHS214*crLHS7 + crLHS216*crLHS8 + crLHS218);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS112 - N[2]*crLHS212 - crLHS115*crLHS210 + crLHS216*crLHS30 - crLHS217*crLHS86 + crLHS219);
rLHS(6,8)+=-gauss_weight*(crLHS118*crLHS207 + crLHS118*crLHS96 + crLHS119*crLHS189 + crLHS220);
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS120 - N[0]*crLHS222 + crLHS113 + crLHS150 - crLHS221*crLHS53);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS125 + crLHS128*crLHS207 + crLHS128*crLHS96 - crLHS129*crLHS189 + crLHS156 + crLHS208);
rLHS(7,2)+=-gauss_weight*(crLHS131*crLHS207 + crLHS131*crLHS96 + crLHS132*crLHS189 + crLHS165);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS133 - N[1]*crLHS222 + crLHS187 + crLHS199 - crLHS221*crLHS84);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS81 + DN(2,1)*crLHS139 + crLHS141*crLHS207 + crLHS141*crLHS96 - crLHS142*crLHS189 + crLHS201 + crLHS213);
rLHS(7,5)+=-gauss_weight*(crLHS146*crLHS207 + crLHS146*crLHS96 + crLHS147*crLHS189 + crLHS204);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS93 + DN(2,1)*crLHS148 - N[2]*crLHS222 - crLHS115*crLHS221 - crLHS137*crLHS217 + crLHS216*crLHS34 + crLHS219);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS151 + crLHS153*crLHS207 + crLHS153*crLHS96 - crLHS154*crLHS189 + crLHS216*crLHS51 + crLHS218 + crLHS223*crLHS7);
rLHS(7,8)+=-gauss_weight*(crLHS158*crLHS207 + crLHS158*crLHS96 + crLHS159*crLHS189 + crLHS224);
rLHS(8,0)+=gauss_weight*(crLHS117 - crLHS118*crLHS21 + crLHS137*crLHS165);
rLHS(8,1)+=gauss_weight*(-crLHS127*crLHS158 + crLHS157 + crLHS164*crLHS86);
rLHS(8,2)+=crLHS166;
rLHS(8,3)+=gauss_weight*(-crLHS118*crLHS71 + crLHS137*crLHS204 + crLHS192);
rLHS(8,4)+=gauss_weight*(-crLHS140*crLHS158 + crLHS202 + crLHS203*crLHS86);
rLHS(8,5)+=crLHS205;
rLHS(8,6)+=gauss_weight*(-crLHS102*crLHS118 + crLHS137*crLHS224 + crLHS220);
rLHS(8,7)+=gauss_weight*(-crLHS152*crLHS158 + crLHS220*crLHS86 + crLHS224);
rLHS(8,8)+=crLHS160*(crLHS214 + crLHS223);
    
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
    const BoundedMatrix<double,4,4>& C = rData.C;

    // Get shape function values
    const array_1d<double,4>& N = rData.N;
    const BoundedMatrix<double,4,3>& DN = rData.DN_DX;

    // const double dyn_tau = rData.DynamicTau;
    // Stabilization parameters 
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;
    constexpr double stab_c3 = 2.0;

    const BoundedMatrix<double,3,4> v_ns = rData.Velocity; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL LHS MATRIX
    const double crLHS0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double crLHS1 = C(0,3)*DN(0,0);
const double crLHS2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crLHS1;
const double crLHS3 = C(0,5)*DN(0,0);
const double crLHS4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crLHS3;
const double crLHS5 = DN(0,0)*DN(0,0);
const double crLHS6 = h*h;
const double crLHS7 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crLHS8 = crLHS7*stab_c3;
const double crLHS9 = crLHS6*crLHS8*1.0/stab_c1 + mu;
const double crLHS10 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crLHS11 = N[0]*N[0];
const double crLHS12 = crLHS11*rho;
const double crLHS13 = N[0]*crLHS7;
const double crLHS14 = N[0]*rho;
const double crLHS15 = crLHS10*crLHS14;
const double crLHS16 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0) + N[3]*v_ns(3,0);
const double crLHS17 = DN(0,0)*crLHS16;
const double crLHS18 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crLHS19 = DN(0,1)*crLHS18;
const double crLHS20 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crLHS21 = DN(0,2)*crLHS20;
const double crLHS22 = -crLHS13 + crLHS17*rho + crLHS19*rho + crLHS21*rho;
const double crLHS23 = -crLHS15 + crLHS22;
const double crLHS24 = 1.0*1.0/(crLHS8 + mu*stab_c1*1.0/crLHS6);
const double crLHS25 = crLHS23*crLHS24;
const double crLHS26 = crLHS17 + crLHS19 + crLHS21;
const double crLHS27 = crLHS26*rho;
const double crLHS28 = DN(0,1)*v_ns(0,0) + DN(1,1)*v_ns(1,0) + DN(2,1)*v_ns(2,0) + DN(3,1)*v_ns(3,0);
const double crLHS29 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crLHS30 = crLHS14*crLHS29;
const double crLHS31 = crLHS28*crLHS30;
const double crLHS32 = DN(0,2)*v_ns(0,0) + DN(1,2)*v_ns(1,0) + DN(2,2)*v_ns(2,0) + DN(3,2)*v_ns(3,0);
const double crLHS33 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crLHS34 = crLHS14*crLHS33;
const double crLHS35 = crLHS32*crLHS34;
const double crLHS36 = -crLHS10*crLHS23 + crLHS31 + crLHS35;
const double crLHS37 = crLHS14*crLHS24;
const double crLHS38 = crLHS11*crLHS7;
const double crLHS39 = crLHS14*crLHS26 + crLHS38;
const double crLHS40 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS41 = C(1,3)*DN(0,1);
const double crLHS42 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS41;
const double crLHS43 = C(3,5)*DN(0,0);
const double crLHS44 = C(4,5)*DN(0,2);
const double crLHS45 = C(1,5)*DN(0,1) + crLHS43 + crLHS44;
const double crLHS46 = DN(0,0)*crLHS9;
const double crLHS47 = DN(0,1)*crLHS46;
const double crLHS48 = crLHS24*crLHS29;
const double crLHS49 = crLHS38*rho;
const double crLHS50 = rho*rho;
const double crLHS51 = crLHS26*crLHS50;
const double crLHS52 = N[0]*crLHS51;
const double crLHS53 = DN(0,2)*v_ns(0,1) + DN(1,2)*v_ns(1,1) + DN(2,2)*v_ns(2,1) + DN(3,2)*v_ns(3,1);
const double crLHS54 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crLHS55 = crLHS14*crLHS54;
const double crLHS56 = crLHS22 - crLHS55;
const double crLHS57 = crLHS15*crLHS29 - crLHS29*crLHS56 + crLHS34*crLHS53;
const double crLHS58 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS59 = C(3,4)*DN(0,1);
const double crLHS60 = C(2,3)*DN(0,2) + crLHS43 + crLHS59;
const double crLHS61 = C(2,5)*DN(0,2);
const double crLHS62 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS61;
const double crLHS63 = DN(0,2)*crLHS46;
const double crLHS64 = crLHS24*crLHS33;
const double crLHS65 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crLHS66 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crLHS67 = crLHS14*crLHS66;
const double crLHS68 = crLHS22 - crLHS67;
const double crLHS69 = crLHS15*crLHS33 + crLHS30*crLHS65 - crLHS33*crLHS68;
const double crLHS70 = DN(0,0)*N[0];
const double crLHS71 = DN(0,0)*crLHS24;
const double crLHS72 = DN(0,0)*crLHS10 + DN(0,1)*crLHS29 + DN(0,2)*crLHS33;
const double crLHS73 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS74 = C(0,3)*DN(1,0);
const double crLHS75 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS74;
const double crLHS76 = C(0,5)*DN(1,0);
const double crLHS77 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS76;
const double crLHS78 = N[1]*rho;
const double crLHS79 = crLHS10*crLHS78;
const double crLHS80 = N[1]*crLHS7;
const double crLHS81 = DN(1,0)*crLHS16;
const double crLHS82 = DN(1,1)*crLHS18;
const double crLHS83 = DN(1,2)*crLHS20;
const double crLHS84 = -crLHS80 + crLHS81*rho + crLHS82*rho + crLHS83*rho;
const double crLHS85 = -crLHS79 + crLHS84;
const double crLHS86 = crLHS24*crLHS85;
const double crLHS87 = crLHS29*crLHS78;
const double crLHS88 = crLHS28*crLHS87;
const double crLHS89 = crLHS33*crLHS78;
const double crLHS90 = crLHS32*crLHS89;
const double crLHS91 = -crLHS10*crLHS85 + crLHS88 + crLHS90;
const double crLHS92 = N[1]*crLHS13;
const double crLHS93 = crLHS26*crLHS78 + crLHS92;
const double crLHS94 = DN(0,0)*DN(1,0);
const double crLHS95 = N[1]*crLHS15 + crLHS9*crLHS94;
const double crLHS96 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS74;
const double crLHS97 = C(1,3)*DN(1,1);
const double crLHS98 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS97;
const double crLHS99 = C(3,5)*DN(1,0);
const double crLHS100 = C(4,5)*DN(1,2);
const double crLHS101 = C(1,5)*DN(1,1) + crLHS100 + crLHS99;
const double crLHS102 = DN(1,1)*crLHS46;
const double crLHS103 = crLHS54*crLHS78;
const double crLHS104 = -crLHS103 + crLHS84;
const double crLHS105 = -crLHS104*crLHS29 + crLHS29*crLHS79 + crLHS53*crLHS89;
const double crLHS106 = N[1]*crLHS51;
const double crLHS107 = crLHS92*rho;
const double crLHS108 = N[1]*crLHS30 - crLHS107*crLHS48;
const double crLHS109 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS76;
const double crLHS110 = C(3,4)*DN(1,1);
const double crLHS111 = C(2,3)*DN(1,2) + crLHS110 + crLHS99;
const double crLHS112 = C(2,5)*DN(1,2);
const double crLHS113 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS112;
const double crLHS114 = DN(1,2)*crLHS46;
const double crLHS115 = crLHS66*crLHS78;
const double crLHS116 = -crLHS115 + crLHS84;
const double crLHS117 = -crLHS116*crLHS33 + crLHS33*crLHS79 + crLHS65*crLHS87;
const double crLHS118 = N[1]*crLHS34 - crLHS107*crLHS64;
const double crLHS119 = DN(0,0)*N[1];
const double crLHS120 = DN(1,0)*crLHS24;
const double crLHS121 = DN(1,0)*crLHS10 + DN(1,1)*crLHS29 + DN(1,2)*crLHS33;
const double crLHS122 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS123 = C(0,3)*DN(2,0);
const double crLHS124 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS123;
const double crLHS125 = C(0,5)*DN(2,0);
const double crLHS126 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS125;
const double crLHS127 = N[2]*rho;
const double crLHS128 = crLHS10*crLHS127;
const double crLHS129 = N[2]*crLHS7;
const double crLHS130 = DN(2,0)*crLHS16;
const double crLHS131 = DN(2,1)*crLHS18;
const double crLHS132 = DN(2,2)*crLHS20;
const double crLHS133 = -crLHS129 + crLHS130*rho + crLHS131*rho + crLHS132*rho;
const double crLHS134 = -crLHS128 + crLHS133;
const double crLHS135 = crLHS134*crLHS24;
const double crLHS136 = crLHS127*crLHS29;
const double crLHS137 = crLHS136*crLHS28;
const double crLHS138 = crLHS127*crLHS33;
const double crLHS139 = crLHS138*crLHS32;
const double crLHS140 = -crLHS10*crLHS134 + crLHS137 + crLHS139;
const double crLHS141 = N[2]*crLHS13;
const double crLHS142 = crLHS127*crLHS26 + crLHS141;
const double crLHS143 = DN(0,0)*DN(2,0);
const double crLHS144 = N[2]*crLHS15 + crLHS143*crLHS9;
const double crLHS145 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS123;
const double crLHS146 = C(1,3)*DN(2,1);
const double crLHS147 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS146;
const double crLHS148 = C(3,5)*DN(2,0);
const double crLHS149 = C(4,5)*DN(2,2);
const double crLHS150 = C(1,5)*DN(2,1) + crLHS148 + crLHS149;
const double crLHS151 = DN(2,1)*crLHS46;
const double crLHS152 = crLHS127*crLHS54;
const double crLHS153 = crLHS133 - crLHS152;
const double crLHS154 = crLHS128*crLHS29 + crLHS138*crLHS53 - crLHS153*crLHS29;
const double crLHS155 = N[2]*crLHS51;
const double crLHS156 = crLHS141*rho;
const double crLHS157 = N[2]*crLHS30 - crLHS156*crLHS48;
const double crLHS158 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS125;
const double crLHS159 = C(3,4)*DN(2,1);
const double crLHS160 = C(2,3)*DN(2,2) + crLHS148 + crLHS159;
const double crLHS161 = C(2,5)*DN(2,2);
const double crLHS162 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS161;
const double crLHS163 = DN(2,2)*crLHS46;
const double crLHS164 = crLHS127*crLHS66;
const double crLHS165 = crLHS133 - crLHS164;
const double crLHS166 = crLHS128*crLHS33 + crLHS136*crLHS65 - crLHS165*crLHS33;
const double crLHS167 = N[2]*crLHS34 - crLHS156*crLHS64;
const double crLHS168 = DN(0,0)*N[2];
const double crLHS169 = DN(2,0)*crLHS24;
const double crLHS170 = DN(2,0)*crLHS10 + DN(2,1)*crLHS29 + DN(2,2)*crLHS33;
const double crLHS171 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS172 = C(0,3)*DN(3,0);
const double crLHS173 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS172;
const double crLHS174 = C(0,5)*DN(3,0);
const double crLHS175 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS174;
const double crLHS176 = N[3]*rho;
const double crLHS177 = crLHS10*crLHS176;
const double crLHS178 = N[3]*crLHS7;
const double crLHS179 = DN(3,0)*crLHS16;
const double crLHS180 = DN(3,1)*crLHS18;
const double crLHS181 = DN(3,2)*crLHS20;
const double crLHS182 = -crLHS178 + crLHS179*rho + crLHS180*rho + crLHS181*rho;
const double crLHS183 = -crLHS177 + crLHS182;
const double crLHS184 = crLHS183*crLHS24;
const double crLHS185 = crLHS176*crLHS29;
const double crLHS186 = crLHS185*crLHS28;
const double crLHS187 = crLHS176*crLHS33;
const double crLHS188 = crLHS187*crLHS32;
const double crLHS189 = -crLHS10*crLHS183 + crLHS186 + crLHS188;
const double crLHS190 = N[3]*crLHS13;
const double crLHS191 = crLHS176*crLHS26 + crLHS190;
const double crLHS192 = DN(0,0)*DN(3,0);
const double crLHS193 = N[3]*crLHS15 + crLHS192*crLHS9;
const double crLHS194 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS172;
const double crLHS195 = C(1,3)*DN(3,1);
const double crLHS196 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS195;
const double crLHS197 = C(3,5)*DN(3,0);
const double crLHS198 = C(4,5)*DN(3,2);
const double crLHS199 = C(1,5)*DN(3,1) + crLHS197 + crLHS198;
const double crLHS200 = DN(3,1)*crLHS46;
const double crLHS201 = crLHS176*crLHS54;
const double crLHS202 = crLHS182 - crLHS201;
const double crLHS203 = crLHS177*crLHS29 + crLHS187*crLHS53 - crLHS202*crLHS29;
const double crLHS204 = N[3]*crLHS51;
const double crLHS205 = crLHS190*rho;
const double crLHS206 = N[3]*crLHS30 - crLHS205*crLHS48;
const double crLHS207 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS174;
const double crLHS208 = C(3,4)*DN(3,1);
const double crLHS209 = C(2,3)*DN(3,2) + crLHS197 + crLHS208;
const double crLHS210 = C(2,5)*DN(3,2);
const double crLHS211 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS210;
const double crLHS212 = DN(3,2)*crLHS46;
const double crLHS213 = crLHS176*crLHS66;
const double crLHS214 = crLHS182 - crLHS213;
const double crLHS215 = crLHS177*crLHS33 + crLHS185*crLHS65 - crLHS214*crLHS33;
const double crLHS216 = N[3]*crLHS34 - crLHS205*crLHS64;
const double crLHS217 = DN(0,0)*N[3];
const double crLHS218 = DN(3,0)*crLHS24;
const double crLHS219 = DN(3,0)*crLHS10 + DN(3,1)*crLHS29 + DN(3,2)*crLHS33;
const double crLHS220 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS41;
const double crLHS221 = C(0,4)*DN(0,0) + crLHS44 + crLHS59;
const double crLHS222 = crLHS24*crLHS28;
const double crLHS223 = crLHS14*crLHS65;
const double crLHS224 = crLHS223*crLHS32 - crLHS23*crLHS28 + crLHS28*crLHS55;
const double crLHS225 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS226 = C(1,4)*DN(0,1);
const double crLHS227 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS226;
const double crLHS228 = DN(0,1)*DN(0,1);
const double crLHS229 = crLHS24*crLHS56;
const double crLHS230 = crLHS223*crLHS53;
const double crLHS231 = crLHS230 + crLHS31 - crLHS54*crLHS56;
const double crLHS232 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS226;
const double crLHS233 = C(2,4)*DN(0,2);
const double crLHS234 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS233;
const double crLHS235 = DN(0,1)*crLHS9;
const double crLHS236 = DN(0,2)*crLHS235;
const double crLHS237 = crLHS24*crLHS65;
const double crLHS238 = crLHS28*crLHS34 + crLHS55*crLHS65 - crLHS65*crLHS68;
const double crLHS239 = DN(0,1)*N[0];
const double crLHS240 = DN(0,1)*crLHS24;
const double crLHS241 = DN(0,0)*crLHS28 + DN(0,1)*crLHS54 + DN(0,2)*crLHS65;
const double crLHS242 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS97;
const double crLHS243 = C(0,4)*DN(1,0) + crLHS100 + crLHS110;
const double crLHS244 = DN(1,0)*crLHS235;
const double crLHS245 = crLHS65*crLHS78;
const double crLHS246 = crLHS103*crLHS28 + crLHS245*crLHS32 - crLHS28*crLHS85;
const double crLHS247 = crLHS14*crLHS28;
const double crLHS248 = N[1]*crLHS247 - crLHS107*crLHS222;
const double crLHS249 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS250 = C(1,4)*DN(1,1);
const double crLHS251 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS250;
const double crLHS252 = crLHS104*crLHS24;
const double crLHS253 = crLHS245*crLHS53;
const double crLHS254 = -crLHS104*crLHS54 + crLHS253 + crLHS88;
const double crLHS255 = DN(0,1)*DN(1,1);
const double crLHS256 = N[1]*crLHS55 + crLHS255*crLHS9;
const double crLHS257 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS250;
const double crLHS258 = C(2,4)*DN(1,2);
const double crLHS259 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS258;
const double crLHS260 = DN(1,2)*crLHS235;
const double crLHS261 = crLHS103*crLHS65 - crLHS116*crLHS65 + crLHS28*crLHS89;
const double crLHS262 = N[1]*crLHS223 - crLHS107*crLHS237;
const double crLHS263 = DN(0,1)*N[1];
const double crLHS264 = DN(1,1)*crLHS24;
const double crLHS265 = DN(1,0)*crLHS28 + DN(1,1)*crLHS54 + DN(1,2)*crLHS65;
const double crLHS266 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS146;
const double crLHS267 = C(0,4)*DN(2,0) + crLHS149 + crLHS159;
const double crLHS268 = DN(2,0)*crLHS235;
const double crLHS269 = crLHS127*crLHS65;
const double crLHS270 = -crLHS134*crLHS28 + crLHS152*crLHS28 + crLHS269*crLHS32;
const double crLHS271 = N[2]*crLHS247 - crLHS156*crLHS222;
const double crLHS272 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS273 = C(1,4)*DN(2,1);
const double crLHS274 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS273;
const double crLHS275 = crLHS153*crLHS24;
const double crLHS276 = crLHS269*crLHS53;
const double crLHS277 = crLHS137 - crLHS153*crLHS54 + crLHS276;
const double crLHS278 = DN(0,1)*DN(2,1);
const double crLHS279 = N[2]*crLHS55 + crLHS278*crLHS9;
const double crLHS280 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS273;
const double crLHS281 = C(2,4)*DN(2,2);
const double crLHS282 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS281;
const double crLHS283 = DN(2,2)*crLHS235;
const double crLHS284 = crLHS138*crLHS28 + crLHS152*crLHS65 - crLHS165*crLHS65;
const double crLHS285 = N[2]*crLHS223 - crLHS156*crLHS237;
const double crLHS286 = DN(0,1)*N[2];
const double crLHS287 = DN(2,1)*crLHS24;
const double crLHS288 = DN(2,0)*crLHS28 + DN(2,1)*crLHS54 + DN(2,2)*crLHS65;
const double crLHS289 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS195;
const double crLHS290 = C(0,4)*DN(3,0) + crLHS198 + crLHS208;
const double crLHS291 = DN(3,0)*crLHS235;
const double crLHS292 = crLHS176*crLHS65;
const double crLHS293 = -crLHS183*crLHS28 + crLHS201*crLHS28 + crLHS292*crLHS32;
const double crLHS294 = N[3]*crLHS247 - crLHS205*crLHS222;
const double crLHS295 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS296 = C(1,4)*DN(3,1);
const double crLHS297 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS296;
const double crLHS298 = crLHS202*crLHS24;
const double crLHS299 = crLHS292*crLHS53;
const double crLHS300 = crLHS186 - crLHS202*crLHS54 + crLHS299;
const double crLHS301 = DN(0,1)*DN(3,1);
const double crLHS302 = N[3]*crLHS55 + crLHS301*crLHS9;
const double crLHS303 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS296;
const double crLHS304 = C(2,4)*DN(3,2);
const double crLHS305 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS304;
const double crLHS306 = DN(3,2)*crLHS235;
const double crLHS307 = crLHS187*crLHS28 + crLHS201*crLHS65 - crLHS214*crLHS65;
const double crLHS308 = N[3]*crLHS223 - crLHS205*crLHS237;
const double crLHS309 = DN(0,1)*N[3];
const double crLHS310 = DN(3,1)*crLHS24;
const double crLHS311 = DN(3,0)*crLHS28 + DN(3,1)*crLHS54 + DN(3,2)*crLHS65;
const double crLHS312 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS61;
const double crLHS313 = crLHS24*crLHS32;
const double crLHS314 = -crLHS23*crLHS32 + crLHS247*crLHS53 + crLHS32*crLHS67;
const double crLHS315 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS233;
const double crLHS316 = crLHS24*crLHS53;
const double crLHS317 = crLHS30*crLHS32 - crLHS53*crLHS56 + crLHS53*crLHS67;
const double crLHS318 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS319 = DN(0,2)*DN(0,2);
const double crLHS320 = crLHS24*crLHS68;
const double crLHS321 = crLHS230 + crLHS35 - crLHS66*crLHS68;
const double crLHS322 = DN(0,2)*N[0];
const double crLHS323 = DN(0,2)*crLHS24;
const double crLHS324 = DN(0,0)*crLHS32 + DN(0,1)*crLHS53 + DN(0,2)*crLHS66;
const double crLHS325 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS112;
const double crLHS326 = DN(0,2)*crLHS9;
const double crLHS327 = DN(1,0)*crLHS326;
const double crLHS328 = crLHS28*crLHS53;
const double crLHS329 = crLHS115*crLHS32 - crLHS32*crLHS85 + crLHS328*crLHS78;
const double crLHS330 = N[1]*crLHS14;
const double crLHS331 = -crLHS107*crLHS313 + crLHS32*crLHS330;
const double crLHS332 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS258;
const double crLHS333 = DN(1,1)*crLHS326;
const double crLHS334 = -crLHS104*crLHS53 + crLHS115*crLHS53 + crLHS32*crLHS87;
const double crLHS335 = -crLHS107*crLHS316 + crLHS330*crLHS53;
const double crLHS336 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS337 = crLHS116*crLHS24;
const double crLHS338 = -crLHS116*crLHS66 + crLHS253 + crLHS90;
const double crLHS339 = DN(0,2)*DN(1,2);
const double crLHS340 = N[1]*crLHS67 + crLHS339*crLHS9;
const double crLHS341 = DN(0,2)*N[1];
const double crLHS342 = DN(1,2)*crLHS24;
const double crLHS343 = DN(1,0)*crLHS32 + DN(1,1)*crLHS53 + DN(1,2)*crLHS66;
const double crLHS344 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS161;
const double crLHS345 = DN(2,0)*crLHS326;
const double crLHS346 = crLHS127*crLHS328 - crLHS134*crLHS32 + crLHS164*crLHS32;
const double crLHS347 = N[2]*crLHS14;
const double crLHS348 = -crLHS156*crLHS313 + crLHS32*crLHS347;
const double crLHS349 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS281;
const double crLHS350 = DN(2,1)*crLHS326;
const double crLHS351 = crLHS136*crLHS32 - crLHS153*crLHS53 + crLHS164*crLHS53;
const double crLHS352 = -crLHS156*crLHS316 + crLHS347*crLHS53;
const double crLHS353 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS354 = crLHS165*crLHS24;
const double crLHS355 = crLHS139 - crLHS165*crLHS66 + crLHS276;
const double crLHS356 = DN(0,2)*DN(2,2);
const double crLHS357 = N[2]*crLHS67 + crLHS356*crLHS9;
const double crLHS358 = DN(0,2)*N[2];
const double crLHS359 = DN(2,2)*crLHS24;
const double crLHS360 = DN(2,0)*crLHS32 + DN(2,1)*crLHS53 + DN(2,2)*crLHS66;
const double crLHS361 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS210;
const double crLHS362 = DN(3,0)*crLHS326;
const double crLHS363 = crLHS176*crLHS328 - crLHS183*crLHS32 + crLHS213*crLHS32;
const double crLHS364 = N[3]*crLHS14;
const double crLHS365 = -crLHS205*crLHS313 + crLHS32*crLHS364;
const double crLHS366 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS304;
const double crLHS367 = DN(3,1)*crLHS326;
const double crLHS368 = crLHS185*crLHS32 - crLHS202*crLHS53 + crLHS213*crLHS53;
const double crLHS369 = -crLHS205*crLHS316 + crLHS364*crLHS53;
const double crLHS370 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS371 = crLHS214*crLHS24;
const double crLHS372 = crLHS188 - crLHS214*crLHS66 + crLHS299;
const double crLHS373 = DN(0,2)*DN(3,2);
const double crLHS374 = N[3]*crLHS67 + crLHS373*crLHS9;
const double crLHS375 = DN(0,2)*N[3];
const double crLHS376 = DN(3,2)*crLHS24;
const double crLHS377 = DN(3,0)*crLHS32 + DN(3,1)*crLHS53 + DN(3,2)*crLHS66;
const double crLHS378 = crLHS239*rho;
const double crLHS379 = crLHS322*rho;
const double crLHS380 = crLHS70*rho;
const double crLHS381 = crLHS24*gauss_weight;
const double crLHS382 = DN(1,0)*N[0];
const double crLHS383 = crLHS263*rho;
const double crLHS384 = crLHS341*rho;
const double crLHS385 = DN(1,1)*N[0];
const double crLHS386 = crLHS119*rho;
const double crLHS387 = DN(1,2)*N[0];
const double crLHS388 = crLHS381*(crLHS255 + crLHS339 + crLHS94);
const double crLHS389 = DN(2,0)*N[0];
const double crLHS390 = crLHS286*rho;
const double crLHS391 = crLHS358*rho;
const double crLHS392 = DN(2,1)*N[0];
const double crLHS393 = crLHS168*rho;
const double crLHS394 = DN(2,2)*N[0];
const double crLHS395 = crLHS381*(crLHS143 + crLHS278 + crLHS356);
const double crLHS396 = DN(3,0)*N[0];
const double crLHS397 = crLHS309*rho;
const double crLHS398 = crLHS375*rho;
const double crLHS399 = DN(3,1)*N[0];
const double crLHS400 = crLHS217*rho;
const double crLHS401 = DN(3,2)*N[0];
const double crLHS402 = crLHS381*(crLHS192 + crLHS301 + crLHS373);
const double crLHS403 = crLHS81 + crLHS82 + crLHS83;
const double crLHS404 = crLHS403*rho;
const double crLHS405 = crLHS24*crLHS78;
const double crLHS406 = crLHS14*crLHS403 + crLHS92;
const double crLHS407 = crLHS403*crLHS50;
const double crLHS408 = N[0]*crLHS407;
const double crLHS409 = DN(1,0)*DN(1,0);
const double crLHS410 = N[1]*N[1];
const double crLHS411 = crLHS410*rho;
const double crLHS412 = crLHS410*crLHS7;
const double crLHS413 = crLHS403*crLHS78 + crLHS412;
const double crLHS414 = DN(1,0)*crLHS9;
const double crLHS415 = DN(1,1)*crLHS414;
const double crLHS416 = crLHS412*rho;
const double crLHS417 = N[1]*crLHS407;
const double crLHS418 = DN(1,2)*crLHS414;
const double crLHS419 = DN(1,0)*N[1];
const double crLHS420 = N[2]*crLHS80;
const double crLHS421 = crLHS127*crLHS403 + crLHS420;
const double crLHS422 = DN(1,0)*DN(2,0);
const double crLHS423 = N[2]*crLHS79 + crLHS422*crLHS9;
const double crLHS424 = DN(2,1)*crLHS414;
const double crLHS425 = N[2]*crLHS407;
const double crLHS426 = crLHS24*crLHS80;
const double crLHS427 = N[2]*crLHS87 - crLHS136*crLHS426;
const double crLHS428 = DN(2,2)*crLHS414;
const double crLHS429 = N[2]*crLHS89 - crLHS138*crLHS426;
const double crLHS430 = DN(1,0)*N[2];
const double crLHS431 = N[3]*crLHS80;
const double crLHS432 = crLHS176*crLHS403 + crLHS431;
const double crLHS433 = DN(1,0)*DN(3,0);
const double crLHS434 = N[3]*crLHS79 + crLHS433*crLHS9;
const double crLHS435 = DN(3,1)*crLHS414;
const double crLHS436 = N[3]*crLHS407;
const double crLHS437 = N[3]*crLHS87 - crLHS185*crLHS426;
const double crLHS438 = DN(3,2)*crLHS414;
const double crLHS439 = N[3]*crLHS89 - crLHS187*crLHS426;
const double crLHS440 = DN(1,0)*N[3];
const double crLHS441 = DN(1,1)*DN(1,1);
const double crLHS442 = DN(1,1)*crLHS9;
const double crLHS443 = DN(1,2)*crLHS442;
const double crLHS444 = DN(1,1)*N[1];
const double crLHS445 = DN(2,0)*crLHS442;
const double crLHS446 = crLHS28*crLHS78;
const double crLHS447 = crLHS127*crLHS24;
const double crLHS448 = crLHS28*crLHS80;
const double crLHS449 = N[2]*crLHS446 - crLHS447*crLHS448;
const double crLHS450 = DN(1,1)*DN(2,1);
const double crLHS451 = N[2]*crLHS103 + crLHS450*crLHS9;
const double crLHS452 = DN(2,2)*crLHS442;
const double crLHS453 = N[2]*crLHS245 - crLHS269*crLHS426;
const double crLHS454 = DN(1,1)*N[2];
const double crLHS455 = DN(3,0)*crLHS442;
const double crLHS456 = crLHS176*crLHS24;
const double crLHS457 = N[3]*crLHS446 - crLHS448*crLHS456;
const double crLHS458 = DN(1,1)*DN(3,1);
const double crLHS459 = N[3]*crLHS103 + crLHS458*crLHS9;
const double crLHS460 = DN(3,2)*crLHS442;
const double crLHS461 = N[3]*crLHS245 - crLHS292*crLHS426;
const double crLHS462 = DN(1,1)*N[3];
const double crLHS463 = DN(1,2)*DN(1,2);
const double crLHS464 = DN(1,2)*N[1];
const double crLHS465 = DN(1,2)*crLHS9;
const double crLHS466 = DN(2,0)*crLHS465;
const double crLHS467 = N[2]*crLHS78;
const double crLHS468 = crLHS447*crLHS80;
const double crLHS469 = crLHS32*crLHS467 - crLHS32*crLHS468;
const double crLHS470 = DN(2,1)*crLHS465;
const double crLHS471 = crLHS467*crLHS53 - crLHS468*crLHS53;
const double crLHS472 = DN(1,2)*DN(2,2);
const double crLHS473 = N[2]*crLHS115 + crLHS472*crLHS9;
const double crLHS474 = DN(1,2)*N[2];
const double crLHS475 = DN(3,0)*crLHS465;
const double crLHS476 = N[3]*crLHS78;
const double crLHS477 = crLHS456*crLHS80;
const double crLHS478 = crLHS32*crLHS476 - crLHS32*crLHS477;
const double crLHS479 = DN(3,1)*crLHS465;
const double crLHS480 = crLHS476*crLHS53 - crLHS477*crLHS53;
const double crLHS481 = DN(1,2)*DN(3,2);
const double crLHS482 = N[3]*crLHS115 + crLHS481*crLHS9;
const double crLHS483 = DN(1,2)*N[3];
const double crLHS484 = crLHS385*rho;
const double crLHS485 = crLHS387*rho;
const double crLHS486 = crLHS382*rho;
const double crLHS487 = crLHS444*rho;
const double crLHS488 = crLHS464*rho;
const double crLHS489 = crLHS419*rho;
const double crLHS490 = DN(2,0)*N[1];
const double crLHS491 = crLHS454*rho;
const double crLHS492 = crLHS474*rho;
const double crLHS493 = DN(2,1)*N[1];
const double crLHS494 = crLHS430*rho;
const double crLHS495 = DN(2,2)*N[1];
const double crLHS496 = crLHS381*(crLHS422 + crLHS450 + crLHS472);
const double crLHS497 = DN(3,0)*N[1];
const double crLHS498 = crLHS462*rho;
const double crLHS499 = crLHS483*rho;
const double crLHS500 = DN(3,1)*N[1];
const double crLHS501 = crLHS440*rho;
const double crLHS502 = DN(3,2)*N[1];
const double crLHS503 = crLHS381*(crLHS433 + crLHS458 + crLHS481);
const double crLHS504 = crLHS130 + crLHS131 + crLHS132;
const double crLHS505 = crLHS504*rho;
const double crLHS506 = crLHS14*crLHS504 + crLHS141;
const double crLHS507 = crLHS50*crLHS504;
const double crLHS508 = N[0]*crLHS507;
const double crLHS509 = crLHS420 + crLHS504*crLHS78;
const double crLHS510 = N[1]*crLHS507;
const double crLHS511 = DN(2,0)*DN(2,0);
const double crLHS512 = N[2]*N[2];
const double crLHS513 = crLHS512*rho;
const double crLHS514 = crLHS512*crLHS7;
const double crLHS515 = crLHS127*crLHS504 + crLHS514;
const double crLHS516 = DN(2,0)*crLHS9;
const double crLHS517 = DN(2,1)*crLHS516;
const double crLHS518 = crLHS514*rho;
const double crLHS519 = N[2]*crLHS507;
const double crLHS520 = DN(2,2)*crLHS516;
const double crLHS521 = DN(2,0)*N[2];
const double crLHS522 = N[3]*crLHS129;
const double crLHS523 = crLHS176*crLHS504 + crLHS522;
const double crLHS524 = DN(2,0)*DN(3,0);
const double crLHS525 = N[3]*crLHS128 + crLHS524*crLHS9;
const double crLHS526 = DN(3,1)*crLHS516;
const double crLHS527 = N[3]*crLHS507;
const double crLHS528 = crLHS129*crLHS24;
const double crLHS529 = N[3]*crLHS136 - crLHS185*crLHS528;
const double crLHS530 = DN(3,2)*crLHS516;
const double crLHS531 = N[3]*crLHS138 - crLHS187*crLHS528;
const double crLHS532 = DN(2,0)*N[3];
const double crLHS533 = DN(2,1)*DN(2,1);
const double crLHS534 = DN(2,1)*crLHS9;
const double crLHS535 = DN(2,2)*crLHS534;
const double crLHS536 = DN(2,1)*N[2];
const double crLHS537 = DN(3,0)*crLHS534;
const double crLHS538 = N[3]*crLHS127;
const double crLHS539 = crLHS129*crLHS456;
const double crLHS540 = crLHS28*crLHS538 - crLHS28*crLHS539;
const double crLHS541 = DN(2,1)*DN(3,1);
const double crLHS542 = N[3]*crLHS152 + crLHS541*crLHS9;
const double crLHS543 = DN(3,2)*crLHS534;
const double crLHS544 = N[3]*crLHS269 - crLHS292*crLHS528;
const double crLHS545 = DN(2,1)*N[3];
const double crLHS546 = DN(2,2)*DN(2,2);
const double crLHS547 = DN(2,2)*N[2];
const double crLHS548 = DN(2,2)*crLHS9;
const double crLHS549 = DN(3,0)*crLHS548;
const double crLHS550 = crLHS32*crLHS538 - crLHS32*crLHS539;
const double crLHS551 = DN(3,1)*crLHS548;
const double crLHS552 = crLHS53*crLHS538 - crLHS53*crLHS539;
const double crLHS553 = DN(2,2)*DN(3,2);
const double crLHS554 = N[3]*crLHS164 + crLHS553*crLHS9;
const double crLHS555 = DN(2,2)*N[3];
const double crLHS556 = crLHS392*rho;
const double crLHS557 = crLHS394*rho;
const double crLHS558 = crLHS389*rho;
const double crLHS559 = crLHS493*rho;
const double crLHS560 = crLHS495*rho;
const double crLHS561 = crLHS490*rho;
const double crLHS562 = crLHS536*rho;
const double crLHS563 = crLHS547*rho;
const double crLHS564 = crLHS521*rho;
const double crLHS565 = DN(3,0)*N[2];
const double crLHS566 = crLHS545*rho;
const double crLHS567 = crLHS555*rho;
const double crLHS568 = DN(3,1)*N[2];
const double crLHS569 = crLHS532*rho;
const double crLHS570 = DN(3,2)*N[2];
const double crLHS571 = crLHS381*(crLHS524 + crLHS541 + crLHS553);
const double crLHS572 = crLHS179 + crLHS180 + crLHS181;
const double crLHS573 = crLHS572*rho;
const double crLHS574 = crLHS14*crLHS572 + crLHS190;
const double crLHS575 = crLHS50*crLHS572;
const double crLHS576 = N[0]*crLHS575;
const double crLHS577 = crLHS431 + crLHS572*crLHS78;
const double crLHS578 = N[1]*crLHS575;
const double crLHS579 = crLHS127*crLHS572 + crLHS522;
const double crLHS580 = N[2]*crLHS575;
const double crLHS581 = DN(3,0)*DN(3,0);
const double crLHS582 = N[3]*N[3];
const double crLHS583 = crLHS582*rho;
const double crLHS584 = crLHS582*crLHS7;
const double crLHS585 = crLHS176*crLHS572 + crLHS584;
const double crLHS586 = DN(3,0)*crLHS9;
const double crLHS587 = DN(3,1)*crLHS586;
const double crLHS588 = crLHS584*rho;
const double crLHS589 = N[3]*crLHS575;
const double crLHS590 = DN(3,2)*crLHS586;
const double crLHS591 = DN(3,0)*N[3];
const double crLHS592 = DN(3,1)*DN(3,1);
const double crLHS593 = DN(3,1)*DN(3,2)*crLHS9;
const double crLHS594 = DN(3,1)*N[3];
const double crLHS595 = DN(3,2)*DN(3,2);
const double crLHS596 = DN(3,2)*N[3];
const double crLHS597 = crLHS399*rho;
const double crLHS598 = crLHS401*rho;
const double crLHS599 = crLHS396*rho;
const double crLHS600 = crLHS500*rho;
const double crLHS601 = crLHS502*rho;
const double crLHS602 = crLHS497*rho;
const double crLHS603 = crLHS568*rho;
const double crLHS604 = crLHS570*rho;
const double crLHS605 = crLHS565*rho;
const double crLHS606 = crLHS594*rho;
const double crLHS607 = crLHS596*rho;
const double crLHS608 = crLHS591*rho;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS10*crLHS12 + crLHS13*crLHS25 + crLHS25*crLHS27 - crLHS36*crLHS37 + crLHS39 + crLHS5*crLHS9);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS40 + DN(0,1)*crLHS42 + DN(0,2)*crLHS45 + crLHS12*crLHS29 - crLHS37*crLHS57 + crLHS47 - crLHS48*crLHS49 - crLHS48*crLHS52);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS58 + DN(0,1)*crLHS60 + DN(0,2)*crLHS62 + crLHS12*crLHS33 - crLHS37*crLHS69 - crLHS49*crLHS64 - crLHS52*crLHS64 + crLHS63);
rLHS(0,3)+=-gauss_weight*(crLHS13*crLHS71 + crLHS27*crLHS71 + crLHS37*crLHS72 + crLHS70);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS73 + DN(0,1)*crLHS75 + DN(0,2)*crLHS77 + crLHS13*crLHS86 + crLHS27*crLHS86 - crLHS37*crLHS91 + crLHS93 + crLHS95);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS96 + DN(0,1)*crLHS98 + DN(0,2)*crLHS101 + crLHS102 - crLHS105*crLHS37 - crLHS106*crLHS48 + crLHS108);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS109 + DN(0,1)*crLHS111 + DN(0,2)*crLHS113 - crLHS106*crLHS64 + crLHS114 - crLHS117*crLHS37 + crLHS118);
rLHS(0,7)+=-gauss_weight*(crLHS119 + crLHS120*crLHS13 + crLHS120*crLHS27 + crLHS121*crLHS37);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS122 + DN(0,1)*crLHS124 + DN(0,2)*crLHS126 + crLHS13*crLHS135 + crLHS135*crLHS27 - crLHS140*crLHS37 + crLHS142 + crLHS144);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS145 + DN(0,1)*crLHS147 + DN(0,2)*crLHS150 + crLHS151 - crLHS154*crLHS37 - crLHS155*crLHS48 + crLHS157);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS158 + DN(0,1)*crLHS160 + DN(0,2)*crLHS162 - crLHS155*crLHS64 + crLHS163 - crLHS166*crLHS37 + crLHS167);
rLHS(0,11)+=-gauss_weight*(crLHS13*crLHS169 + crLHS168 + crLHS169*crLHS27 + crLHS170*crLHS37);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS171 + DN(0,1)*crLHS173 + DN(0,2)*crLHS175 + crLHS13*crLHS184 + crLHS184*crLHS27 - crLHS189*crLHS37 + crLHS191 + crLHS193);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS194 + DN(0,1)*crLHS196 + DN(0,2)*crLHS199 + crLHS200 - crLHS203*crLHS37 - crLHS204*crLHS48 + crLHS206);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS207 + DN(0,1)*crLHS209 + DN(0,2)*crLHS211 - crLHS204*crLHS64 + crLHS212 - crLHS215*crLHS37 + crLHS216);
rLHS(0,15)+=-gauss_weight*(crLHS13*crLHS218 + crLHS217 + crLHS218*crLHS27 + crLHS219*crLHS37);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS220 + DN(0,2)*crLHS221 + crLHS12*crLHS28 - crLHS222*crLHS49 - crLHS222*crLHS52 - crLHS224*crLHS37 + crLHS47);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS42 + DN(0,1)*crLHS225 + DN(0,2)*crLHS227 + crLHS12*crLHS54 + crLHS13*crLHS229 + crLHS228*crLHS9 + crLHS229*crLHS27 - crLHS231*crLHS37 + crLHS39);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS60 + DN(0,1)*crLHS232 + DN(0,2)*crLHS234 + crLHS12*crLHS65 + crLHS236 - crLHS237*crLHS49 - crLHS237*crLHS52 - crLHS238*crLHS37);
rLHS(1,3)+=-gauss_weight*(crLHS13*crLHS240 + crLHS239 + crLHS240*crLHS27 + crLHS241*crLHS37);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS242 + DN(0,2)*crLHS243 - crLHS106*crLHS222 + crLHS244 - crLHS246*crLHS37 + crLHS248);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS98 + DN(0,1)*crLHS249 + DN(0,2)*crLHS251 + crLHS13*crLHS252 + crLHS252*crLHS27 - crLHS254*crLHS37 + crLHS256 + crLHS93);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS111 + DN(0,1)*crLHS257 + DN(0,2)*crLHS259 - crLHS106*crLHS237 + crLHS260 - crLHS261*crLHS37 + crLHS262);
rLHS(1,7)+=-gauss_weight*(crLHS13*crLHS264 + crLHS263 + crLHS264*crLHS27 + crLHS265*crLHS37);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS124 + DN(0,1)*crLHS266 + DN(0,2)*crLHS267 - crLHS155*crLHS222 + crLHS268 - crLHS270*crLHS37 + crLHS271);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS147 + DN(0,1)*crLHS272 + DN(0,2)*crLHS274 + crLHS13*crLHS275 + crLHS142 + crLHS27*crLHS275 - crLHS277*crLHS37 + crLHS279);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS160 + DN(0,1)*crLHS280 + DN(0,2)*crLHS282 - crLHS155*crLHS237 + crLHS283 - crLHS284*crLHS37 + crLHS285);
rLHS(1,11)+=-gauss_weight*(crLHS13*crLHS287 + crLHS27*crLHS287 + crLHS286 + crLHS288*crLHS37);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS173 + DN(0,1)*crLHS289 + DN(0,2)*crLHS290 - crLHS204*crLHS222 + crLHS291 - crLHS293*crLHS37 + crLHS294);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS196 + DN(0,1)*crLHS295 + DN(0,2)*crLHS297 + crLHS13*crLHS298 + crLHS191 + crLHS27*crLHS298 - crLHS300*crLHS37 + crLHS302);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS209 + DN(0,1)*crLHS303 + DN(0,2)*crLHS305 - crLHS204*crLHS237 + crLHS306 - crLHS307*crLHS37 + crLHS308);
rLHS(1,15)+=-gauss_weight*(crLHS13*crLHS310 + crLHS27*crLHS310 + crLHS309 + crLHS311*crLHS37);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS221 + DN(0,2)*crLHS312 + crLHS12*crLHS32 - crLHS313*crLHS49 - crLHS313*crLHS52 - crLHS314*crLHS37 + crLHS63);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS227 + DN(0,2)*crLHS315 + crLHS12*crLHS53 + crLHS236 - crLHS316*crLHS49 - crLHS316*crLHS52 - crLHS317*crLHS37);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS234 + DN(0,2)*crLHS318 + crLHS12*crLHS66 + crLHS13*crLHS320 + crLHS27*crLHS320 + crLHS319*crLHS9 - crLHS321*crLHS37 + crLHS39);
rLHS(2,3)+=-gauss_weight*(crLHS13*crLHS323 + crLHS27*crLHS323 + crLHS322 + crLHS324*crLHS37);
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS243 + DN(0,2)*crLHS325 - crLHS106*crLHS313 + crLHS327 - crLHS329*crLHS37 + crLHS331);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS101 + DN(0,1)*crLHS251 + DN(0,2)*crLHS332 - crLHS106*crLHS316 + crLHS333 - crLHS334*crLHS37 + crLHS335);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS113 + DN(0,1)*crLHS259 + DN(0,2)*crLHS336 + crLHS13*crLHS337 + crLHS27*crLHS337 - crLHS338*crLHS37 + crLHS340 + crLHS93);
rLHS(2,7)+=-gauss_weight*(crLHS13*crLHS342 + crLHS27*crLHS342 + crLHS341 + crLHS343*crLHS37);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS126 + DN(0,1)*crLHS267 + DN(0,2)*crLHS344 - crLHS155*crLHS313 + crLHS345 - crLHS346*crLHS37 + crLHS348);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS150 + DN(0,1)*crLHS274 + DN(0,2)*crLHS349 - crLHS155*crLHS316 + crLHS350 - crLHS351*crLHS37 + crLHS352);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS162 + DN(0,1)*crLHS282 + DN(0,2)*crLHS353 + crLHS13*crLHS354 + crLHS142 + crLHS27*crLHS354 - crLHS355*crLHS37 + crLHS357);
rLHS(2,11)+=-gauss_weight*(crLHS13*crLHS359 + crLHS27*crLHS359 + crLHS358 + crLHS360*crLHS37);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS175 + DN(0,1)*crLHS290 + DN(0,2)*crLHS361 - crLHS204*crLHS313 + crLHS362 - crLHS363*crLHS37 + crLHS365);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS199 + DN(0,1)*crLHS297 + DN(0,2)*crLHS366 - crLHS204*crLHS316 + crLHS367 - crLHS368*crLHS37 + crLHS369);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS211 + DN(0,1)*crLHS305 + DN(0,2)*crLHS370 + crLHS13*crLHS371 + crLHS191 + crLHS27*crLHS371 - crLHS37*crLHS372 + crLHS374);
rLHS(2,15)+=-gauss_weight*(crLHS13*crLHS376 + crLHS27*crLHS376 + crLHS37*crLHS377 + crLHS375);
rLHS(3,0)+=gauss_weight*(crLHS222*crLHS378 - crLHS23*crLHS71 + crLHS313*crLHS379 + crLHS70);
rLHS(3,1)+=gauss_weight*(crLHS239 - crLHS240*crLHS56 + crLHS316*crLHS379 + crLHS380*crLHS48);
rLHS(3,2)+=gauss_weight*(crLHS237*crLHS378 + crLHS322 - crLHS323*crLHS68 + crLHS380*crLHS64);
rLHS(3,3)+=crLHS381*(crLHS228 + crLHS319 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS222*crLHS383 + crLHS313*crLHS384 + crLHS382 - crLHS71*crLHS85);
rLHS(3,5)+=gauss_weight*(-crLHS104*crLHS240 + crLHS316*crLHS384 + crLHS385 + crLHS386*crLHS48);
rLHS(3,6)+=gauss_weight*(-crLHS116*crLHS323 + crLHS237*crLHS383 + crLHS386*crLHS64 + crLHS387);
rLHS(3,7)+=crLHS388;
rLHS(3,8)+=gauss_weight*(-crLHS134*crLHS71 + crLHS222*crLHS390 + crLHS313*crLHS391 + crLHS389);
rLHS(3,9)+=gauss_weight*(-crLHS153*crLHS240 + crLHS316*crLHS391 + crLHS392 + crLHS393*crLHS48);
rLHS(3,10)+=gauss_weight*(-crLHS165*crLHS323 + crLHS237*crLHS390 + crLHS393*crLHS64 + crLHS394);
rLHS(3,11)+=crLHS395;
rLHS(3,12)+=gauss_weight*(-crLHS183*crLHS71 + crLHS222*crLHS397 + crLHS313*crLHS398 + crLHS396);
rLHS(3,13)+=gauss_weight*(-crLHS202*crLHS240 + crLHS316*crLHS398 + crLHS399 + crLHS400*crLHS48);
rLHS(3,14)+=gauss_weight*(-crLHS214*crLHS323 + crLHS237*crLHS397 + crLHS400*crLHS64 + crLHS401);
rLHS(3,15)+=crLHS402;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS25*crLHS404 + crLHS25*crLHS80 - crLHS36*crLHS405 + crLHS406 + crLHS95);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS40 + DN(1,1)*crLHS42 + DN(1,2)*crLHS45 + crLHS108 + crLHS244 - crLHS405*crLHS57 - crLHS408*crLHS48);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS58 + DN(1,1)*crLHS60 + DN(1,2)*crLHS62 + crLHS118 + crLHS327 - crLHS405*crLHS69 - crLHS408*crLHS64);
rLHS(4,3)+=-gauss_weight*(crLHS382 + crLHS404*crLHS71 + crLHS405*crLHS72 + crLHS71*crLHS80);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS73 + DN(1,1)*crLHS75 + DN(1,2)*crLHS77 + crLHS10*crLHS411 + crLHS404*crLHS86 - crLHS405*crLHS91 + crLHS409*crLHS9 + crLHS413 + crLHS80*crLHS86);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS96 + DN(1,1)*crLHS98 + DN(1,2)*crLHS101 - crLHS105*crLHS405 + crLHS29*crLHS411 + crLHS415 - crLHS416*crLHS48 - crLHS417*crLHS48);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS109 + DN(1,1)*crLHS111 + DN(1,2)*crLHS113 - crLHS117*crLHS405 + crLHS33*crLHS411 - crLHS416*crLHS64 - crLHS417*crLHS64 + crLHS418);
rLHS(4,7)+=-gauss_weight*(crLHS120*crLHS404 + crLHS120*crLHS80 + crLHS121*crLHS405 + crLHS419);
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS122 + DN(1,1)*crLHS124 + DN(1,2)*crLHS126 + crLHS135*crLHS404 + crLHS135*crLHS80 - crLHS140*crLHS405 + crLHS421 + crLHS423);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS145 + DN(1,1)*crLHS147 + DN(1,2)*crLHS150 - crLHS154*crLHS405 + crLHS424 - crLHS425*crLHS48 + crLHS427);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS158 + DN(1,1)*crLHS160 + DN(1,2)*crLHS162 - crLHS166*crLHS405 - crLHS425*crLHS64 + crLHS428 + crLHS429);
rLHS(4,11)+=-gauss_weight*(crLHS169*crLHS404 + crLHS169*crLHS80 + crLHS170*crLHS405 + crLHS430);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS171 + DN(1,1)*crLHS173 + DN(1,2)*crLHS175 + crLHS184*crLHS404 + crLHS184*crLHS80 - crLHS189*crLHS405 + crLHS432 + crLHS434);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS194 + DN(1,1)*crLHS196 + DN(1,2)*crLHS199 - crLHS203*crLHS405 + crLHS435 - crLHS436*crLHS48 + crLHS437);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS207 + DN(1,1)*crLHS209 + DN(1,2)*crLHS211 - crLHS215*crLHS405 - crLHS436*crLHS64 + crLHS438 + crLHS439);
rLHS(4,15)+=-gauss_weight*(crLHS218*crLHS404 + crLHS218*crLHS80 + crLHS219*crLHS405 + crLHS440);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS220 + DN(1,2)*crLHS221 + crLHS102 - crLHS222*crLHS408 - crLHS224*crLHS405 + crLHS248);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS42 + DN(1,1)*crLHS225 + DN(1,2)*crLHS227 + crLHS229*crLHS404 + crLHS229*crLHS80 - crLHS231*crLHS405 + crLHS256 + crLHS406);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS60 + DN(1,1)*crLHS232 + DN(1,2)*crLHS234 - crLHS237*crLHS408 - crLHS238*crLHS405 + crLHS262 + crLHS333);
rLHS(5,3)+=-gauss_weight*(crLHS240*crLHS404 + crLHS240*crLHS80 + crLHS241*crLHS405 + crLHS385);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS242 + DN(1,2)*crLHS243 - crLHS222*crLHS416 - crLHS222*crLHS417 - crLHS246*crLHS405 + crLHS28*crLHS411 + crLHS415);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS98 + DN(1,1)*crLHS249 + DN(1,2)*crLHS251 + crLHS252*crLHS404 + crLHS252*crLHS80 - crLHS254*crLHS405 + crLHS411*crLHS54 + crLHS413 + crLHS441*crLHS9);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS111 + DN(1,1)*crLHS257 + DN(1,2)*crLHS259 - crLHS237*crLHS416 - crLHS237*crLHS417 - crLHS261*crLHS405 + crLHS411*crLHS65 + crLHS443);
rLHS(5,7)+=-gauss_weight*(crLHS264*crLHS404 + crLHS264*crLHS80 + crLHS265*crLHS405 + crLHS444);
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS124 + DN(1,1)*crLHS266 + DN(1,2)*crLHS267 - crLHS222*crLHS425 - crLHS270*crLHS405 + crLHS445 + crLHS449);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS147 + DN(1,1)*crLHS272 + DN(1,2)*crLHS274 + crLHS275*crLHS404 + crLHS275*crLHS80 - crLHS277*crLHS405 + crLHS421 + crLHS451);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS160 + DN(1,1)*crLHS280 + DN(1,2)*crLHS282 - crLHS237*crLHS425 - crLHS284*crLHS405 + crLHS452 + crLHS453);
rLHS(5,11)+=-gauss_weight*(crLHS287*crLHS404 + crLHS287*crLHS80 + crLHS288*crLHS405 + crLHS454);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS173 + DN(1,1)*crLHS289 + DN(1,2)*crLHS290 - crLHS222*crLHS436 - crLHS293*crLHS405 + crLHS455 + crLHS457);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS196 + DN(1,1)*crLHS295 + DN(1,2)*crLHS297 + crLHS298*crLHS404 + crLHS298*crLHS80 - crLHS300*crLHS405 + crLHS432 + crLHS459);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS209 + DN(1,1)*crLHS303 + DN(1,2)*crLHS305 - crLHS237*crLHS436 - crLHS307*crLHS405 + crLHS460 + crLHS461);
rLHS(5,15)+=-gauss_weight*(crLHS310*crLHS404 + crLHS310*crLHS80 + crLHS311*crLHS405 + crLHS462);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS221 + DN(1,2)*crLHS312 + crLHS114 - crLHS313*crLHS408 - crLHS314*crLHS405 + crLHS331);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS227 + DN(1,2)*crLHS315 + crLHS260 - crLHS316*crLHS408 - crLHS317*crLHS405 + crLHS335);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS234 + DN(1,2)*crLHS318 + crLHS320*crLHS404 + crLHS320*crLHS80 - crLHS321*crLHS405 + crLHS340 + crLHS406);
rLHS(6,3)+=-gauss_weight*(crLHS323*crLHS404 + crLHS323*crLHS80 + crLHS324*crLHS405 + crLHS387);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS243 + DN(1,2)*crLHS325 - crLHS313*crLHS416 - crLHS313*crLHS417 + crLHS32*crLHS411 - crLHS329*crLHS405 + crLHS418);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS101 + DN(1,1)*crLHS251 + DN(1,2)*crLHS332 - crLHS316*crLHS416 - crLHS316*crLHS417 - crLHS334*crLHS405 + crLHS411*crLHS53 + crLHS443);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS113 + DN(1,1)*crLHS259 + DN(1,2)*crLHS336 + crLHS337*crLHS404 + crLHS337*crLHS80 - crLHS338*crLHS405 + crLHS411*crLHS66 + crLHS413 + crLHS463*crLHS9);
rLHS(6,7)+=-gauss_weight*(crLHS342*crLHS404 + crLHS342*crLHS80 + crLHS343*crLHS405 + crLHS464);
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS126 + DN(1,1)*crLHS267 + DN(1,2)*crLHS344 - crLHS313*crLHS425 - crLHS346*crLHS405 + crLHS466 + crLHS469);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS150 + DN(1,1)*crLHS274 + DN(1,2)*crLHS349 - crLHS316*crLHS425 - crLHS351*crLHS405 + crLHS470 + crLHS471);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS162 + DN(1,1)*crLHS282 + DN(1,2)*crLHS353 + crLHS354*crLHS404 + crLHS354*crLHS80 - crLHS355*crLHS405 + crLHS421 + crLHS473);
rLHS(6,11)+=-gauss_weight*(crLHS359*crLHS404 + crLHS359*crLHS80 + crLHS360*crLHS405 + crLHS474);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS175 + DN(1,1)*crLHS290 + DN(1,2)*crLHS361 - crLHS313*crLHS436 - crLHS363*crLHS405 + crLHS475 + crLHS478);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS199 + DN(1,1)*crLHS297 + DN(1,2)*crLHS366 - crLHS316*crLHS436 - crLHS368*crLHS405 + crLHS479 + crLHS480);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS211 + DN(1,1)*crLHS305 + DN(1,2)*crLHS370 + crLHS371*crLHS404 + crLHS371*crLHS80 - crLHS372*crLHS405 + crLHS432 + crLHS482);
rLHS(6,15)+=-gauss_weight*(crLHS376*crLHS404 + crLHS376*crLHS80 + crLHS377*crLHS405 + crLHS483);
rLHS(7,0)+=gauss_weight*(crLHS119 - crLHS120*crLHS23 + crLHS222*crLHS484 + crLHS313*crLHS485);
rLHS(7,1)+=gauss_weight*(crLHS263 - crLHS264*crLHS56 + crLHS316*crLHS485 + crLHS48*crLHS486);
rLHS(7,2)+=gauss_weight*(crLHS237*crLHS484 + crLHS341 - crLHS342*crLHS68 + crLHS486*crLHS64);
rLHS(7,3)+=crLHS388;
rLHS(7,4)+=gauss_weight*(-crLHS120*crLHS85 + crLHS222*crLHS487 + crLHS313*crLHS488 + crLHS419);
rLHS(7,5)+=gauss_weight*(-crLHS104*crLHS264 + crLHS316*crLHS488 + crLHS444 + crLHS48*crLHS489);
rLHS(7,6)+=gauss_weight*(-crLHS116*crLHS342 + crLHS237*crLHS487 + crLHS464 + crLHS489*crLHS64);
rLHS(7,7)+=crLHS381*(crLHS409 + crLHS441 + crLHS463);
rLHS(7,8)+=gauss_weight*(-crLHS120*crLHS134 + crLHS222*crLHS491 + crLHS313*crLHS492 + crLHS490);
rLHS(7,9)+=gauss_weight*(-crLHS153*crLHS264 + crLHS316*crLHS492 + crLHS48*crLHS494 + crLHS493);
rLHS(7,10)+=gauss_weight*(-crLHS165*crLHS342 + crLHS237*crLHS491 + crLHS494*crLHS64 + crLHS495);
rLHS(7,11)+=crLHS496;
rLHS(7,12)+=gauss_weight*(-crLHS120*crLHS183 + crLHS222*crLHS498 + crLHS313*crLHS499 + crLHS497);
rLHS(7,13)+=gauss_weight*(-crLHS202*crLHS264 + crLHS316*crLHS499 + crLHS48*crLHS501 + crLHS500);
rLHS(7,14)+=gauss_weight*(-crLHS214*crLHS342 + crLHS237*crLHS498 + crLHS501*crLHS64 + crLHS502);
rLHS(7,15)+=crLHS503;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS129*crLHS25 + crLHS144 + crLHS25*crLHS505 - crLHS36*crLHS447 + crLHS506);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS40 + DN(2,1)*crLHS42 + DN(2,2)*crLHS45 + crLHS157 + crLHS268 - crLHS447*crLHS57 - crLHS48*crLHS508);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS58 + DN(2,1)*crLHS60 + DN(2,2)*crLHS62 + crLHS167 + crLHS345 - crLHS447*crLHS69 - crLHS508*crLHS64);
rLHS(8,3)+=-gauss_weight*(crLHS129*crLHS71 + crLHS389 + crLHS447*crLHS72 + crLHS505*crLHS71);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS73 + DN(2,1)*crLHS75 + DN(2,2)*crLHS77 + crLHS129*crLHS86 + crLHS423 - crLHS447*crLHS91 + crLHS505*crLHS86 + crLHS509);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS96 + DN(2,1)*crLHS98 + DN(2,2)*crLHS101 - crLHS105*crLHS447 + crLHS427 + crLHS445 - crLHS48*crLHS510);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS109 + DN(2,1)*crLHS111 + DN(2,2)*crLHS113 - crLHS117*crLHS447 + crLHS429 + crLHS466 - crLHS510*crLHS64);
rLHS(8,7)+=-gauss_weight*(crLHS120*crLHS129 + crLHS120*crLHS505 + crLHS121*crLHS447 + crLHS490);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS122 + DN(2,1)*crLHS124 + DN(2,2)*crLHS126 + crLHS10*crLHS513 + crLHS129*crLHS135 + crLHS135*crLHS505 - crLHS140*crLHS447 + crLHS511*crLHS9 + crLHS515);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS145 + DN(2,1)*crLHS147 + DN(2,2)*crLHS150 - crLHS154*crLHS447 + crLHS29*crLHS513 - crLHS48*crLHS518 - crLHS48*crLHS519 + crLHS517);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS158 + DN(2,1)*crLHS160 + DN(2,2)*crLHS162 - crLHS166*crLHS447 + crLHS33*crLHS513 - crLHS518*crLHS64 - crLHS519*crLHS64 + crLHS520);
rLHS(8,11)+=-gauss_weight*(crLHS129*crLHS169 + crLHS169*crLHS505 + crLHS170*crLHS447 + crLHS521);
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS171 + DN(2,1)*crLHS173 + DN(2,2)*crLHS175 + crLHS129*crLHS184 + crLHS184*crLHS505 - crLHS189*crLHS447 + crLHS523 + crLHS525);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS194 + DN(2,1)*crLHS196 + DN(2,2)*crLHS199 - crLHS203*crLHS447 - crLHS48*crLHS527 + crLHS526 + crLHS529);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS207 + DN(2,1)*crLHS209 + DN(2,2)*crLHS211 - crLHS215*crLHS447 - crLHS527*crLHS64 + crLHS530 + crLHS531);
rLHS(8,15)+=-gauss_weight*(crLHS129*crLHS218 + crLHS218*crLHS505 + crLHS219*crLHS447 + crLHS532);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS220 + DN(2,2)*crLHS221 + crLHS151 - crLHS222*crLHS508 - crLHS224*crLHS447 + crLHS271);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS42 + DN(2,1)*crLHS225 + DN(2,2)*crLHS227 + crLHS129*crLHS229 + crLHS229*crLHS505 - crLHS231*crLHS447 + crLHS279 + crLHS506);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS60 + DN(2,1)*crLHS232 + DN(2,2)*crLHS234 - crLHS237*crLHS508 - crLHS238*crLHS447 + crLHS285 + crLHS350);
rLHS(9,3)+=-gauss_weight*(crLHS129*crLHS240 + crLHS240*crLHS505 + crLHS241*crLHS447 + crLHS392);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS242 + DN(2,2)*crLHS243 - crLHS222*crLHS510 - crLHS246*crLHS447 + crLHS424 + crLHS449);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS98 + DN(2,1)*crLHS249 + DN(2,2)*crLHS251 + crLHS129*crLHS252 + crLHS252*crLHS505 - crLHS254*crLHS447 + crLHS451 + crLHS509);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS111 + DN(2,1)*crLHS257 + DN(2,2)*crLHS259 - crLHS237*crLHS510 - crLHS261*crLHS447 + crLHS453 + crLHS470);
rLHS(9,7)+=-gauss_weight*(crLHS129*crLHS264 + crLHS264*crLHS505 + crLHS265*crLHS447 + crLHS493);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS124 + DN(2,1)*crLHS266 + DN(2,2)*crLHS267 - crLHS222*crLHS518 - crLHS222*crLHS519 - crLHS270*crLHS447 + crLHS28*crLHS513 + crLHS517);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS147 + DN(2,1)*crLHS272 + DN(2,2)*crLHS274 + crLHS129*crLHS275 + crLHS275*crLHS505 - crLHS277*crLHS447 + crLHS513*crLHS54 + crLHS515 + crLHS533*crLHS9);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS160 + DN(2,1)*crLHS280 + DN(2,2)*crLHS282 - crLHS237*crLHS518 - crLHS237*crLHS519 - crLHS284*crLHS447 + crLHS513*crLHS65 + crLHS535);
rLHS(9,11)+=-gauss_weight*(crLHS129*crLHS287 + crLHS287*crLHS505 + crLHS288*crLHS447 + crLHS536);
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS173 + DN(2,1)*crLHS289 + DN(2,2)*crLHS290 - crLHS222*crLHS527 - crLHS293*crLHS447 + crLHS537 + crLHS540);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS196 + DN(2,1)*crLHS295 + DN(2,2)*crLHS297 + crLHS129*crLHS298 + crLHS298*crLHS505 - crLHS300*crLHS447 + crLHS523 + crLHS542);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS209 + DN(2,1)*crLHS303 + DN(2,2)*crLHS305 - crLHS237*crLHS527 - crLHS307*crLHS447 + crLHS543 + crLHS544);
rLHS(9,15)+=-gauss_weight*(crLHS129*crLHS310 + crLHS310*crLHS505 + crLHS311*crLHS447 + crLHS545);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS221 + DN(2,2)*crLHS312 + crLHS163 - crLHS313*crLHS508 - crLHS314*crLHS447 + crLHS348);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS227 + DN(2,2)*crLHS315 + crLHS283 - crLHS316*crLHS508 - crLHS317*crLHS447 + crLHS352);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS234 + DN(2,2)*crLHS318 + crLHS129*crLHS320 + crLHS320*crLHS505 - crLHS321*crLHS447 + crLHS357 + crLHS506);
rLHS(10,3)+=-gauss_weight*(crLHS129*crLHS323 + crLHS323*crLHS505 + crLHS324*crLHS447 + crLHS394);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS243 + DN(2,2)*crLHS325 - crLHS313*crLHS510 - crLHS329*crLHS447 + crLHS428 + crLHS469);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS101 + DN(2,1)*crLHS251 + DN(2,2)*crLHS332 - crLHS316*crLHS510 - crLHS334*crLHS447 + crLHS452 + crLHS471);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS113 + DN(2,1)*crLHS259 + DN(2,2)*crLHS336 + crLHS129*crLHS337 + crLHS337*crLHS505 - crLHS338*crLHS447 + crLHS473 + crLHS509);
rLHS(10,7)+=-gauss_weight*(crLHS129*crLHS342 + crLHS342*crLHS505 + crLHS343*crLHS447 + crLHS495);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS126 + DN(2,1)*crLHS267 + DN(2,2)*crLHS344 - crLHS313*crLHS518 - crLHS313*crLHS519 + crLHS32*crLHS513 - crLHS346*crLHS447 + crLHS520);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS150 + DN(2,1)*crLHS274 + DN(2,2)*crLHS349 - crLHS316*crLHS518 - crLHS316*crLHS519 - crLHS351*crLHS447 + crLHS513*crLHS53 + crLHS535);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS162 + DN(2,1)*crLHS282 + DN(2,2)*crLHS353 + crLHS129*crLHS354 + crLHS354*crLHS505 - crLHS355*crLHS447 + crLHS513*crLHS66 + crLHS515 + crLHS546*crLHS9);
rLHS(10,11)+=-gauss_weight*(crLHS129*crLHS359 + crLHS359*crLHS505 + crLHS360*crLHS447 + crLHS547);
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS175 + DN(2,1)*crLHS290 + DN(2,2)*crLHS361 - crLHS313*crLHS527 - crLHS363*crLHS447 + crLHS549 + crLHS550);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS199 + DN(2,1)*crLHS297 + DN(2,2)*crLHS366 - crLHS316*crLHS527 - crLHS368*crLHS447 + crLHS551 + crLHS552);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS211 + DN(2,1)*crLHS305 + DN(2,2)*crLHS370 + crLHS129*crLHS371 + crLHS371*crLHS505 - crLHS372*crLHS447 + crLHS523 + crLHS554);
rLHS(10,15)+=-gauss_weight*(crLHS129*crLHS376 + crLHS376*crLHS505 + crLHS377*crLHS447 + crLHS555);
rLHS(11,0)+=gauss_weight*(crLHS168 - crLHS169*crLHS23 + crLHS222*crLHS556 + crLHS313*crLHS557);
rLHS(11,1)+=gauss_weight*(crLHS286 - crLHS287*crLHS56 + crLHS316*crLHS557 + crLHS48*crLHS558);
rLHS(11,2)+=gauss_weight*(crLHS237*crLHS556 + crLHS358 - crLHS359*crLHS68 + crLHS558*crLHS64);
rLHS(11,3)+=crLHS395;
rLHS(11,4)+=gauss_weight*(-crLHS169*crLHS85 + crLHS222*crLHS559 + crLHS313*crLHS560 + crLHS430);
rLHS(11,5)+=gauss_weight*(-crLHS104*crLHS287 + crLHS316*crLHS560 + crLHS454 + crLHS48*crLHS561);
rLHS(11,6)+=gauss_weight*(-crLHS116*crLHS359 + crLHS237*crLHS559 + crLHS474 + crLHS561*crLHS64);
rLHS(11,7)+=crLHS496;
rLHS(11,8)+=gauss_weight*(-crLHS134*crLHS169 + crLHS222*crLHS562 + crLHS313*crLHS563 + crLHS521);
rLHS(11,9)+=gauss_weight*(-crLHS153*crLHS287 + crLHS316*crLHS563 + crLHS48*crLHS564 + crLHS536);
rLHS(11,10)+=gauss_weight*(-crLHS165*crLHS359 + crLHS237*crLHS562 + crLHS547 + crLHS564*crLHS64);
rLHS(11,11)+=crLHS381*(crLHS511 + crLHS533 + crLHS546);
rLHS(11,12)+=gauss_weight*(-crLHS169*crLHS183 + crLHS222*crLHS566 + crLHS313*crLHS567 + crLHS565);
rLHS(11,13)+=gauss_weight*(-crLHS202*crLHS287 + crLHS316*crLHS567 + crLHS48*crLHS569 + crLHS568);
rLHS(11,14)+=gauss_weight*(-crLHS214*crLHS359 + crLHS237*crLHS566 + crLHS569*crLHS64 + crLHS570);
rLHS(11,15)+=crLHS571;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS178*crLHS25 + crLHS193 + crLHS25*crLHS573 - crLHS36*crLHS456 + crLHS574);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS40 + DN(3,1)*crLHS42 + DN(3,2)*crLHS45 + crLHS206 + crLHS291 - crLHS456*crLHS57 - crLHS48*crLHS576);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS58 + DN(3,1)*crLHS60 + DN(3,2)*crLHS62 + crLHS216 + crLHS362 - crLHS456*crLHS69 - crLHS576*crLHS64);
rLHS(12,3)+=-gauss_weight*(crLHS178*crLHS71 + crLHS396 + crLHS456*crLHS72 + crLHS573*crLHS71);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS73 + DN(3,1)*crLHS75 + DN(3,2)*crLHS77 + crLHS178*crLHS86 + crLHS434 - crLHS456*crLHS91 + crLHS573*crLHS86 + crLHS577);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS96 + DN(3,1)*crLHS98 + DN(3,2)*crLHS101 - crLHS105*crLHS456 + crLHS437 + crLHS455 - crLHS48*crLHS578);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS109 + DN(3,1)*crLHS111 + DN(3,2)*crLHS113 - crLHS117*crLHS456 + crLHS439 + crLHS475 - crLHS578*crLHS64);
rLHS(12,7)+=-gauss_weight*(crLHS120*crLHS178 + crLHS120*crLHS573 + crLHS121*crLHS456 + crLHS497);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS122 + DN(3,1)*crLHS124 + DN(3,2)*crLHS126 + crLHS135*crLHS178 + crLHS135*crLHS573 - crLHS140*crLHS456 + crLHS525 + crLHS579);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS145 + DN(3,1)*crLHS147 + DN(3,2)*crLHS150 - crLHS154*crLHS456 - crLHS48*crLHS580 + crLHS529 + crLHS537);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS158 + DN(3,1)*crLHS160 + DN(3,2)*crLHS162 - crLHS166*crLHS456 + crLHS531 + crLHS549 - crLHS580*crLHS64);
rLHS(12,11)+=-gauss_weight*(crLHS169*crLHS178 + crLHS169*crLHS573 + crLHS170*crLHS456 + crLHS565);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS171 + DN(3,1)*crLHS173 + DN(3,2)*crLHS175 + crLHS10*crLHS583 + crLHS178*crLHS184 + crLHS184*crLHS573 - crLHS189*crLHS456 + crLHS581*crLHS9 + crLHS585);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS194 + DN(3,1)*crLHS196 + DN(3,2)*crLHS199 - crLHS203*crLHS456 + crLHS29*crLHS583 - crLHS48*crLHS588 - crLHS48*crLHS589 + crLHS587);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS207 + DN(3,1)*crLHS209 + DN(3,2)*crLHS211 - crLHS215*crLHS456 + crLHS33*crLHS583 - crLHS588*crLHS64 - crLHS589*crLHS64 + crLHS590);
rLHS(12,15)+=-gauss_weight*(crLHS178*crLHS218 + crLHS218*crLHS573 + crLHS219*crLHS456 + crLHS591);
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS220 + DN(3,2)*crLHS221 + crLHS200 - crLHS222*crLHS576 - crLHS224*crLHS456 + crLHS294);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS42 + DN(3,1)*crLHS225 + DN(3,2)*crLHS227 + crLHS178*crLHS229 + crLHS229*crLHS573 - crLHS231*crLHS456 + crLHS302 + crLHS574);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS60 + DN(3,1)*crLHS232 + DN(3,2)*crLHS234 - crLHS237*crLHS576 - crLHS238*crLHS456 + crLHS308 + crLHS367);
rLHS(13,3)+=-gauss_weight*(crLHS178*crLHS240 + crLHS240*crLHS573 + crLHS241*crLHS456 + crLHS399);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS75 + DN(3,1)*crLHS242 + DN(3,2)*crLHS243 - crLHS222*crLHS578 - crLHS246*crLHS456 + crLHS435 + crLHS457);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS98 + DN(3,1)*crLHS249 + DN(3,2)*crLHS251 + crLHS178*crLHS252 + crLHS252*crLHS573 - crLHS254*crLHS456 + crLHS459 + crLHS577);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS111 + DN(3,1)*crLHS257 + DN(3,2)*crLHS259 - crLHS237*crLHS578 - crLHS261*crLHS456 + crLHS461 + crLHS479);
rLHS(13,7)+=-gauss_weight*(crLHS178*crLHS264 + crLHS264*crLHS573 + crLHS265*crLHS456 + crLHS500);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS124 + DN(3,1)*crLHS266 + DN(3,2)*crLHS267 - crLHS222*crLHS580 - crLHS270*crLHS456 + crLHS526 + crLHS540);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS147 + DN(3,1)*crLHS272 + DN(3,2)*crLHS274 + crLHS178*crLHS275 + crLHS275*crLHS573 - crLHS277*crLHS456 + crLHS542 + crLHS579);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS160 + DN(3,1)*crLHS280 + DN(3,2)*crLHS282 - crLHS237*crLHS580 - crLHS284*crLHS456 + crLHS544 + crLHS551);
rLHS(13,11)+=-gauss_weight*(crLHS178*crLHS287 + crLHS287*crLHS573 + crLHS288*crLHS456 + crLHS568);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS173 + DN(3,1)*crLHS289 + DN(3,2)*crLHS290 - crLHS222*crLHS588 - crLHS222*crLHS589 + crLHS28*crLHS583 - crLHS293*crLHS456 + crLHS587);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS196 + DN(3,1)*crLHS295 + DN(3,2)*crLHS297 + crLHS178*crLHS298 + crLHS298*crLHS573 - crLHS300*crLHS456 + crLHS54*crLHS583 + crLHS585 + crLHS592*crLHS9);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS209 + DN(3,1)*crLHS303 + DN(3,2)*crLHS305 - crLHS237*crLHS588 - crLHS237*crLHS589 - crLHS307*crLHS456 + crLHS583*crLHS65 + crLHS593);
rLHS(13,15)+=-gauss_weight*(crLHS178*crLHS310 + crLHS310*crLHS573 + crLHS311*crLHS456 + crLHS594);
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS221 + DN(3,2)*crLHS312 + crLHS212 - crLHS313*crLHS576 - crLHS314*crLHS456 + crLHS365);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS45 + DN(3,1)*crLHS227 + DN(3,2)*crLHS315 + crLHS306 - crLHS316*crLHS576 - crLHS317*crLHS456 + crLHS369);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS62 + DN(3,1)*crLHS234 + DN(3,2)*crLHS318 + crLHS178*crLHS320 + crLHS320*crLHS573 - crLHS321*crLHS456 + crLHS374 + crLHS574);
rLHS(14,3)+=-gauss_weight*(crLHS178*crLHS323 + crLHS323*crLHS573 + crLHS324*crLHS456 + crLHS401);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS243 + DN(3,2)*crLHS325 - crLHS313*crLHS578 - crLHS329*crLHS456 + crLHS438 + crLHS478);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS101 + DN(3,1)*crLHS251 + DN(3,2)*crLHS332 - crLHS316*crLHS578 - crLHS334*crLHS456 + crLHS460 + crLHS480);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS113 + DN(3,1)*crLHS259 + DN(3,2)*crLHS336 + crLHS178*crLHS337 + crLHS337*crLHS573 - crLHS338*crLHS456 + crLHS482 + crLHS577);
rLHS(14,7)+=-gauss_weight*(crLHS178*crLHS342 + crLHS342*crLHS573 + crLHS343*crLHS456 + crLHS502);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS126 + DN(3,1)*crLHS267 + DN(3,2)*crLHS344 - crLHS313*crLHS580 - crLHS346*crLHS456 + crLHS530 + crLHS550);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS150 + DN(3,1)*crLHS274 + DN(3,2)*crLHS349 - crLHS316*crLHS580 - crLHS351*crLHS456 + crLHS543 + crLHS552);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS162 + DN(3,1)*crLHS282 + DN(3,2)*crLHS353 + crLHS178*crLHS354 + crLHS354*crLHS573 - crLHS355*crLHS456 + crLHS554 + crLHS579);
rLHS(14,11)+=-gauss_weight*(crLHS178*crLHS359 + crLHS359*crLHS573 + crLHS360*crLHS456 + crLHS570);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS175 + DN(3,1)*crLHS290 + DN(3,2)*crLHS361 - crLHS313*crLHS588 - crLHS313*crLHS589 + crLHS32*crLHS583 - crLHS363*crLHS456 + crLHS590);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS199 + DN(3,1)*crLHS297 + DN(3,2)*crLHS366 - crLHS316*crLHS588 - crLHS316*crLHS589 - crLHS368*crLHS456 + crLHS53*crLHS583 + crLHS593);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS211 + DN(3,1)*crLHS305 + DN(3,2)*crLHS370 + crLHS178*crLHS371 + crLHS371*crLHS573 - crLHS372*crLHS456 + crLHS583*crLHS66 + crLHS585 + crLHS595*crLHS9);
rLHS(14,15)+=-gauss_weight*(crLHS178*crLHS376 + crLHS376*crLHS573 + crLHS377*crLHS456 + crLHS596);
rLHS(15,0)+=gauss_weight*(crLHS217 - crLHS218*crLHS23 + crLHS222*crLHS597 + crLHS313*crLHS598);
rLHS(15,1)+=gauss_weight*(crLHS309 - crLHS310*crLHS56 + crLHS316*crLHS598 + crLHS48*crLHS599);
rLHS(15,2)+=gauss_weight*(crLHS237*crLHS597 + crLHS375 - crLHS376*crLHS68 + crLHS599*crLHS64);
rLHS(15,3)+=crLHS402;
rLHS(15,4)+=gauss_weight*(-crLHS218*crLHS85 + crLHS222*crLHS600 + crLHS313*crLHS601 + crLHS440);
rLHS(15,5)+=gauss_weight*(-crLHS104*crLHS310 + crLHS316*crLHS601 + crLHS462 + crLHS48*crLHS602);
rLHS(15,6)+=gauss_weight*(-crLHS116*crLHS376 + crLHS237*crLHS600 + crLHS483 + crLHS602*crLHS64);
rLHS(15,7)+=crLHS503;
rLHS(15,8)+=gauss_weight*(-crLHS134*crLHS218 + crLHS222*crLHS603 + crLHS313*crLHS604 + crLHS532);
rLHS(15,9)+=gauss_weight*(-crLHS153*crLHS310 + crLHS316*crLHS604 + crLHS48*crLHS605 + crLHS545);
rLHS(15,10)+=gauss_weight*(-crLHS165*crLHS376 + crLHS237*crLHS603 + crLHS555 + crLHS605*crLHS64);
rLHS(15,11)+=crLHS571;
rLHS(15,12)+=gauss_weight*(-crLHS183*crLHS218 + crLHS222*crLHS606 + crLHS313*crLHS607 + crLHS591);
rLHS(15,13)+=gauss_weight*(-crLHS202*crLHS310 + crLHS316*crLHS607 + crLHS48*crLHS608 + crLHS594);
rLHS(15,14)+=gauss_weight*(-crLHS214*crLHS376 + crLHS237*crLHS606 + crLHS596 + crLHS608*crLHS64);
rLHS(15,15)+=crLHS381*(crLHS581 + crLHS592 + crLHS595);

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
    const BoundedMatrix<double,2,3> v_ns = rData.Velocity; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
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
const double crRHS23 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0);
const double crRHS24 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1);
const double crRHS25 = crRHS23 + crRHS24;
const double crRHS26 = h*h;
const double crRHS27 = crRHS2*stab_c3;
const double crRHS28 = crRHS25*(crRHS26*crRHS27*1.0/stab_c1 + mu);
const double crRHS29 = crRHS18*crRHS3;
const double crRHS30 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1);
const double crRHS31 = crRHS11*crRHS30;
const double crRHS32 = crRHS29 + crRHS31;
const double crRHS33 = N[0]*rho;
const double crRHS34 = crRHS5*rho;
const double crRHS35 = crRHS16*rho;
const double crRHS36 = crRHS14*functional_weights[6];
const double crRHS37 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS14 - crRHS23*crRHS34 + crRHS29*rho + crRHS31*rho - crRHS35*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0)) + crRHS36 + crRHS4 + crRHS7;
const double crRHS38 = 1.0*1.0/(crRHS27 + mu*stab_c1*1.0/crRHS26);
const double crRHS39 = crRHS37*crRHS38;
const double crRHS40 = N[0]*crRHS2;
const double crRHS41 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS42 = crRHS2*crRHS30;
const double crRHS43 = crRHS20*crRHS3;
const double crRHS44 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crRHS45 = crRHS30*crRHS44;
const double crRHS46 = crRHS16*crRHS6;
const double crRHS47 = crRHS13*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2]);
const double crRHS48 = crRHS47*functional_weights[6];
const double crRHS49 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] - crRHS24*crRHS35 - crRHS34*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS41 + crRHS42 + crRHS43*rho + crRHS45*rho + crRHS46 + crRHS47 + crRHS48;
const double crRHS50 = crRHS11*crRHS49 + crRHS18*crRHS37;
const double crRHS51 = crRHS33*crRHS38;
const double crRHS52 = 1.0*crRHS44;
const double crRHS53 = crRHS43 + crRHS45;
const double crRHS54 = crRHS38*crRHS49;
const double crRHS55 = crRHS20*crRHS37 + crRHS44*crRHS49;
const double crRHS56 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS16);
const double crRHS57 = N[1]*rho;
const double crRHS58 = N[1]*crRHS2;
const double crRHS59 = crRHS38*crRHS57;
const double crRHS60 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS16);
const double crRHS61 = N[2]*rho;
const double crRHS62 = N[2]*crRHS2;
const double crRHS63 = crRHS38*crRHS61;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS28 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS12 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS15*functional_weights[6] + crRHS15 + crRHS17*crRHS3 - crRHS17*crRHS39 + crRHS22*(DN(0,0)*crRHS19 + DN(0,1)*crRHS21) + crRHS32*crRHS33 - crRHS39*crRHS40 - crRHS50*crRHS51);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS12 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS28 + DN(0,1)*stress_adj[1] - N[0]*crRHS41 + N[0]*crRHS42 + N[0]*crRHS46 + N[0]*crRHS47 + N[0]*crRHS48 + crRHS17*crRHS30 - crRHS17*crRHS54 + crRHS22*(DN(0,0)*crRHS21 + DN(0,1)*crRHS52) + crRHS33*crRHS53 - crRHS40*crRHS54 - crRHS51*crRHS55);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS39 + DN(0,1)*crRHS54 + N[0]*crRHS25);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS28 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS12 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS14 + N[1]*crRHS36 + N[1]*crRHS4 + N[1]*crRHS7 + crRHS22*(DN(1,0)*crRHS19 + DN(1,1)*crRHS21) + crRHS3*crRHS56 + crRHS32*crRHS57 - crRHS39*crRHS56 - crRHS39*crRHS58 - crRHS50*crRHS59);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS12 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS28 + DN(1,1)*stress_adj[1] - N[1]*crRHS41 + N[1]*crRHS42 + N[1]*crRHS46 + N[1]*crRHS47 + N[1]*crRHS48 + crRHS22*(DN(1,0)*crRHS21 + DN(1,1)*crRHS52) + crRHS30*crRHS56 + crRHS53*crRHS57 - crRHS54*crRHS56 - crRHS54*crRHS58 - crRHS55*crRHS59);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS39 + DN(1,1)*crRHS54 + N[1]*crRHS25);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS28 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS12 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS14 + N[2]*crRHS36 + N[2]*crRHS4 + N[2]*crRHS7 + crRHS22*(DN(2,0)*crRHS19 + DN(2,1)*crRHS21) + crRHS3*crRHS60 + crRHS32*crRHS61 - crRHS39*crRHS60 - crRHS39*crRHS62 - crRHS50*crRHS63);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS12 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS28 + DN(2,1)*stress_adj[1] - N[2]*crRHS41 + N[2]*crRHS42 + N[2]*crRHS46 + N[2]*crRHS47 + N[2]*crRHS48 + crRHS22*(DN(2,0)*crRHS21 + DN(2,1)*crRHS52) + crRHS30*crRHS60 + crRHS53*crRHS61 - crRHS54*crRHS60 - crRHS54*crRHS62 - crRHS55*crRHS63);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS39 + DN(2,1)*crRHS54 + N[2]*crRHS25);

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
    const BoundedMatrix<double,3,4> v_ns = rData.Velocity; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION
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
const double crRHS14 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0) + DN(3,0)*v_adj(3,0);
const double crRHS15 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1) + DN(3,1)*v_adj(3,1);
const double crRHS16 = DN(0,2)*v_adj(0,2) + DN(1,2)*v_adj(1,2) + DN(2,2)*v_adj(2,2) + DN(3,2)*v_adj(3,2);
const double crRHS17 = crRHS14 + crRHS15 + crRHS16;
const double crRHS18 = h*h;
const double crRHS19 = crRHS2*stab_c3;
const double crRHS20 = crRHS17*(crRHS18*crRHS19*1.0/stab_c1 + mu);
const double crRHS21 = DN(0,1)*v_ns(0,0);
const double crRHS22 = DN(1,1)*v_ns(1,0);
const double crRHS23 = DN(2,1)*v_ns(2,0);
const double crRHS24 = DN(3,1)*v_ns(3,0);
const double crRHS25 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crRHS26 = -crRHS21 - crRHS22 - crRHS23 - crRHS24 + crRHS25;
const double crRHS27 = DN(0,2)*v_ns(0,0);
const double crRHS28 = DN(1,2)*v_ns(1,0);
const double crRHS29 = DN(2,2)*v_ns(2,0);
const double crRHS30 = DN(3,2)*v_ns(3,0);
const double crRHS31 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crRHS32 = -crRHS27 - crRHS28 - crRHS29 - crRHS30 + crRHS31;
const double crRHS33 = 2.0*functional_weights[2]*mu;
const double crRHS34 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crRHS35 = 1.0*crRHS34;
const double crRHS36 = crRHS21 + crRHS22 + crRHS23 + crRHS24;
const double crRHS37 = 0.5*crRHS25 + 0.5*crRHS36;
const double crRHS38 = crRHS27 + crRHS28 + crRHS29 + crRHS30;
const double crRHS39 = crRHS31 + crRHS38;
const double crRHS40 = 0.5*DN(0,2);
const double crRHS41 = 4.0*functional_weights[1]*mu;
const double crRHS42 = crRHS3*crRHS34;
const double crRHS43 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1) + N[3]*v_adj(3,1);
const double crRHS44 = crRHS25*crRHS43;
const double crRHS45 = N[0]*v_adj(0,2) + N[1]*v_adj(1,2) + N[2]*v_adj(2,2) + N[3]*v_adj(3,2);
const double crRHS46 = crRHS31*crRHS45;
const double crRHS47 = crRHS42 + crRHS44 + crRHS46;
const double crRHS48 = N[0]*rho;
const double crRHS49 = crRHS5*rho;
const double crRHS50 = crRHS11*rho;
const double crRHS51 = crRHS12*rho;
const double crRHS52 = crRHS9*functional_weights[6];
const double crRHS53 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 - crRHS14*crRHS49 + crRHS4 + crRHS42*rho + crRHS44*rho + crRHS46*rho - crRHS50*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS51*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0)) + crRHS52 + crRHS7 + crRHS9;
const double crRHS54 = 1.0*1.0/(crRHS19 + mu*stab_c1*1.0/crRHS18);
const double crRHS55 = crRHS53*crRHS54;
const double crRHS56 = N[0]*crRHS2;
const double crRHS57 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS58 = crRHS2*crRHS43;
const double crRHS59 = crRHS3*crRHS36;
const double crRHS60 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crRHS61 = crRHS43*crRHS60;
const double crRHS62 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crRHS63 = crRHS45*crRHS62;
const double crRHS64 = crRHS11*crRHS6;
const double crRHS65 = crRHS8*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3]);
const double crRHS66 = crRHS65*functional_weights[6];
const double crRHS67 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] - crRHS15*crRHS50 - crRHS49*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS51*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS57 + crRHS58 + crRHS59*rho + crRHS61*rho + crRHS63*rho + crRHS64 + crRHS65 + crRHS66;
const double crRHS68 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS69 = crRHS2*crRHS45;
const double crRHS70 = crRHS3*crRHS38;
const double crRHS71 = DN(0,2)*v_ns(0,1);
const double crRHS72 = DN(1,2)*v_ns(1,1);
const double crRHS73 = DN(2,2)*v_ns(2,1);
const double crRHS74 = DN(3,2)*v_ns(3,1);
const double crRHS75 = crRHS71 + crRHS72 + crRHS73 + crRHS74;
const double crRHS76 = crRHS43*crRHS75;
const double crRHS77 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crRHS78 = crRHS45*crRHS77;
const double crRHS79 = crRHS12*crRHS6;
const double crRHS80 = crRHS8*(DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3]);
const double crRHS81 = crRHS80*functional_weights[6];
const double crRHS82 = DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] - crRHS16*crRHS51 - crRHS49*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS50*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS68 + crRHS69 + crRHS70*rho + crRHS76*rho + crRHS78*rho + crRHS79 + crRHS80 + crRHS81;
const double crRHS83 = crRHS25*crRHS67 + crRHS31*crRHS82 + crRHS34*crRHS53;
const double crRHS84 = crRHS48*crRHS54;
const double crRHS85 = crRHS62 - crRHS71 - crRHS72 - crRHS73 - crRHS74;
const double crRHS86 = 1.0*crRHS60;
const double crRHS87 = crRHS62 + crRHS75;
const double crRHS88 = crRHS59 + crRHS61 + crRHS63;
const double crRHS89 = crRHS54*crRHS67;
const double crRHS90 = crRHS36*crRHS53 + crRHS60*crRHS67 + crRHS62*crRHS82;
const double crRHS91 = 1.0*crRHS77;
const double crRHS92 = 0.5*crRHS39;
const double crRHS93 = 0.5*crRHS87;
const double crRHS94 = crRHS70 + crRHS76 + crRHS78;
const double crRHS95 = crRHS54*crRHS82;
const double crRHS96 = crRHS38*crRHS53 + crRHS67*crRHS75 + crRHS77*crRHS82;
const double crRHS97 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS11 + DN(1,2)*crRHS12);
const double crRHS98 = N[1]*rho;
const double crRHS99 = N[1]*crRHS2;
const double crRHS100 = crRHS54*crRHS98;
const double crRHS101 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS11 + DN(2,2)*crRHS12);
const double crRHS102 = N[2]*rho;
const double crRHS103 = N[2]*crRHS2;
const double crRHS104 = crRHS102*crRHS54;
const double crRHS105 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS11 + DN(3,2)*crRHS12);
const double crRHS106 = N[3]*rho;
const double crRHS107 = N[3]*crRHS2;
const double crRHS108 = crRHS106*crRHS54;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS20 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS10*functional_weights[6] + crRHS10 + crRHS13*crRHS3 - crRHS13*crRHS55 - crRHS33*(DN(0,1)*crRHS26 + DN(0,2)*crRHS32) + crRHS41*(DN(0,0)*crRHS35 + DN(0,1)*crRHS37 + crRHS39*crRHS40) + crRHS47*crRHS48 - crRHS55*crRHS56 - crRHS83*crRHS84);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS20 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS57 + N[0]*crRHS58 + N[0]*crRHS64 + N[0]*crRHS65 + N[0]*crRHS66 + crRHS13*crRHS43 - crRHS13*crRHS89 + crRHS33*(DN(0,0)*crRHS26 - DN(0,2)*crRHS85) + crRHS41*(DN(0,0)*crRHS37 + DN(0,1)*crRHS86 + crRHS40*crRHS87) + crRHS48*crRHS88 - crRHS56*crRHS89 - crRHS84*crRHS90);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS20 + DN(0,2)*stress_adj[2] - N[0]*crRHS68 + N[0]*crRHS69 + N[0]*crRHS79 + N[0]*crRHS80 + N[0]*crRHS81 + crRHS13*crRHS45 - crRHS13*crRHS95 + crRHS33*(DN(0,0)*crRHS32 + DN(0,1)*crRHS85) + crRHS41*(DN(0,0)*crRHS92 + DN(0,1)*crRHS93 + DN(0,2)*crRHS91) + crRHS48*crRHS94 - crRHS56*crRHS95 - crRHS84*crRHS96);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS55 + DN(0,1)*crRHS89 + DN(0,2)*crRHS95 + N[0]*crRHS17);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS20 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS4 + N[1]*crRHS52 + N[1]*crRHS7 + N[1]*crRHS9 - crRHS100*crRHS83 + crRHS3*crRHS97 - crRHS33*(DN(1,1)*crRHS26 + DN(1,2)*crRHS32) + crRHS41*(DN(1,0)*crRHS35 + DN(1,1)*crRHS37 + DN(1,2)*crRHS92) + crRHS47*crRHS98 - crRHS55*crRHS97 - crRHS55*crRHS99);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS20 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS57 + N[1]*crRHS58 + N[1]*crRHS64 + N[1]*crRHS65 + N[1]*crRHS66 - crRHS100*crRHS90 + crRHS33*(DN(1,0)*crRHS26 - DN(1,2)*crRHS85) + crRHS41*(DN(1,0)*crRHS37 + DN(1,1)*crRHS86 + DN(1,2)*crRHS93) + crRHS43*crRHS97 + crRHS88*crRHS98 - crRHS89*crRHS97 - crRHS89*crRHS99);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS20 + DN(1,2)*stress_adj[2] - N[1]*crRHS68 + N[1]*crRHS69 + N[1]*crRHS79 + N[1]*crRHS80 + N[1]*crRHS81 - crRHS100*crRHS96 + crRHS33*(DN(1,0)*crRHS32 + DN(1,1)*crRHS85) + crRHS41*(DN(1,0)*crRHS92 + DN(1,1)*crRHS93 + DN(1,2)*crRHS91) + crRHS45*crRHS97 + crRHS94*crRHS98 - crRHS95*crRHS97 - crRHS95*crRHS99);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS55 + DN(1,1)*crRHS89 + DN(1,2)*crRHS95 + N[1]*crRHS17);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS20 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS4 + N[2]*crRHS52 + N[2]*crRHS7 + N[2]*crRHS9 + crRHS101*crRHS3 - crRHS101*crRHS55 + crRHS102*crRHS47 - crRHS103*crRHS55 - crRHS104*crRHS83 - crRHS33*(DN(2,1)*crRHS26 + DN(2,2)*crRHS32) + crRHS41*(DN(2,0)*crRHS35 + DN(2,1)*crRHS37 + DN(2,2)*crRHS92));
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS20 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS57 + N[2]*crRHS58 + N[2]*crRHS64 + N[2]*crRHS65 + N[2]*crRHS66 + crRHS101*crRHS43 - crRHS101*crRHS89 + crRHS102*crRHS88 - crRHS103*crRHS89 - crRHS104*crRHS90 + crRHS33*(DN(2,0)*crRHS26 - DN(2,2)*crRHS85) + crRHS41*(DN(2,0)*crRHS37 + DN(2,1)*crRHS86 + DN(2,2)*crRHS93));
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS20 + DN(2,2)*stress_adj[2] - N[2]*crRHS68 + N[2]*crRHS69 + N[2]*crRHS79 + N[2]*crRHS80 + N[2]*crRHS81 + crRHS101*crRHS45 - crRHS101*crRHS95 + crRHS102*crRHS94 - crRHS103*crRHS95 - crRHS104*crRHS96 + crRHS33*(DN(2,0)*crRHS32 + DN(2,1)*crRHS85) + crRHS41*(DN(2,0)*crRHS92 + DN(2,1)*crRHS93 + DN(2,2)*crRHS91));
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS55 + DN(2,1)*crRHS89 + DN(2,2)*crRHS95 + N[2]*crRHS17);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS20 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS4 + N[3]*crRHS52 + N[3]*crRHS7 + N[3]*crRHS9 + crRHS105*crRHS3 - crRHS105*crRHS55 + crRHS106*crRHS47 - crRHS107*crRHS55 - crRHS108*crRHS83 - crRHS33*(DN(3,1)*crRHS26 + DN(3,2)*crRHS32) + crRHS41*(DN(3,0)*crRHS35 + DN(3,1)*crRHS37 + DN(3,2)*crRHS92));
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS20 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS57 + N[3]*crRHS58 + N[3]*crRHS64 + N[3]*crRHS65 + N[3]*crRHS66 + crRHS105*crRHS43 - crRHS105*crRHS89 + crRHS106*crRHS88 - crRHS107*crRHS89 - crRHS108*crRHS90 + crRHS33*(DN(3,0)*crRHS26 - DN(3,2)*crRHS85) + crRHS41*(DN(3,0)*crRHS37 + DN(3,1)*crRHS86 + DN(3,2)*crRHS93));
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS20 + DN(3,2)*stress_adj[2] - N[3]*crRHS68 + N[3]*crRHS69 + N[3]*crRHS79 + N[3]*crRHS80 + N[3]*crRHS81 + crRHS105*crRHS45 - crRHS105*crRHS95 + crRHS106*crRHS94 - crRHS107*crRHS95 - crRHS108*crRHS96 + crRHS33*(DN(3,0)*crRHS32 + DN(3,1)*crRHS85) + crRHS41*(DN(3,0)*crRHS92 + DN(3,1)*crRHS93 + DN(3,2)*crRHS91));
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS55 + DN(3,1)*crRHS89 + DN(3,2)*crRHS95 + N[3]*crRHS17);

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