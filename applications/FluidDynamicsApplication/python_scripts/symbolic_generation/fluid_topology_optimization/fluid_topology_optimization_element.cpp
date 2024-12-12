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

    const BoundedMatrix<double,2,3> v_ns = rData.Velocity; // CONVECTIVE VELOCITY FOR THE ADJOINT IS THE PRIMAL NS SOLUTION

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
const double crLHS8 = pow(crLHS6, 2) + pow(crLHS7, 2);
const double crLHS9 = pow(crLHS8, 1.0/4.0)*stab_c3;
const double crLHS10 = sqrt(crLHS8)*rho*stab_c2;
const double crLHS11 = h*(crLHS10 + crLHS5*h + crLHS9*h)/stab_c1 + mu;
const double crLHS12 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crLHS13 = pow(N[0], 2);
const double crLHS14 = crLHS13*rho;
const double crLHS15 = N[0]*crLHS4;
const double crLHS16 = N[0]*rho;
const double crLHS17 = crLHS12*crLHS16;
const double crLHS18 = DN(0,0)*crLHS6;
const double crLHS19 = crLHS18*rho;
const double crLHS20 = DN(0,1)*crLHS7;
const double crLHS21 = crLHS20*rho;
const double crLHS22 = -crLHS15 + crLHS19 + crLHS21;
const double crLHS23 = -crLHS17 + crLHS22;
const double crLHS24 = 1.0/(crLHS10/h + crLHS5 + crLHS9 + mu*stab_c1/pow(h, 2));
const double crLHS25 = 1.0*crLHS24;
const double crLHS26 = crLHS23*crLHS25;
const double crLHS27 = crLHS18 + crLHS20;
const double crLHS28 = crLHS27*rho;
const double crLHS29 = DN(0,1)*v_ns(0,0);
const double crLHS30 = DN(1,1)*v_ns(1,0);
const double crLHS31 = DN(2,1)*v_ns(2,0);
const double crLHS32 = crLHS29 + crLHS30 + crLHS31;
const double crLHS33 = DN(0,0)*v_ns(0,1);
const double crLHS34 = DN(1,0)*v_ns(1,1);
const double crLHS35 = DN(2,0)*v_ns(2,1);
const double crLHS36 = crLHS33 + crLHS34 + crLHS35;
const double crLHS37 = crLHS32*crLHS36;
const double crLHS38 = crLHS16*crLHS37;
const double crLHS39 = -crLHS12*crLHS23 + crLHS38;
const double crLHS40 = crLHS16*crLHS25;
const double crLHS41 = crLHS13*crLHS4;
const double crLHS42 = crLHS16*crLHS27 + crLHS41;
const double crLHS43 = C(0,1)*DN(0,1) + crLHS1;
const double crLHS44 = C(1,2)*DN(0,1);
const double crLHS45 = C(2,2)*DN(0,0) + crLHS44;
const double crLHS46 = DN(0,0)*crLHS11;
const double crLHS47 = DN(0,1)*crLHS46;
const double crLHS48 = crLHS25*crLHS32;
const double crLHS49 = crLHS41*rho;
const double crLHS50 = pow(rho, 2);
const double crLHS51 = crLHS27*crLHS50;
const double crLHS52 = crLHS48*crLHS51;
const double crLHS53 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crLHS54 = crLHS16*crLHS53;
const double crLHS55 = crLHS15 + crLHS17 - crLHS19 - crLHS21 + crLHS54;
const double crLHS56 = 1.0*crLHS29 + 1.0*crLHS30 + 1.0*crLHS31;
const double crLHS57 = crLHS16*crLHS24;
const double crLHS58 = crLHS56*crLHS57;
const double crLHS59 = DN(0,0)*N[0];
const double crLHS60 = DN(0,0)*crLHS25;
const double crLHS61 = DN(0,0)*crLHS12 + DN(0,1)*crLHS32;
const double crLHS62 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double crLHS63 = C(0,2)*DN(1,0);
const double crLHS64 = C(2,2)*DN(1,1) + crLHS63;
const double crLHS65 = N[1]*rho;
const double crLHS66 = crLHS12*crLHS65;
const double crLHS67 = N[1]*crLHS4;
const double crLHS68 = DN(1,0)*crLHS6;
const double crLHS69 = crLHS68*rho;
const double crLHS70 = DN(1,1)*crLHS7;
const double crLHS71 = crLHS70*rho;
const double crLHS72 = -crLHS67 + crLHS69 + crLHS71;
const double crLHS73 = -crLHS66 + crLHS72;
const double crLHS74 = crLHS25*crLHS73;
const double crLHS75 = crLHS37*crLHS65;
const double crLHS76 = -crLHS12*crLHS73 + crLHS75;
const double crLHS77 = N[1]*crLHS15;
const double crLHS78 = crLHS27*crLHS65 + crLHS77;
const double crLHS79 = DN(0,0)*DN(1,0);
const double crLHS80 = N[1]*crLHS17 + crLHS11*crLHS79;
const double crLHS81 = C(0,1)*DN(1,1) + crLHS63;
const double crLHS82 = C(1,2)*DN(1,1);
const double crLHS83 = C(2,2)*DN(1,0) + crLHS82;
const double crLHS84 = DN(1,1)*crLHS46;
const double crLHS85 = crLHS53*crLHS65;
const double crLHS86 = crLHS66 + crLHS67 - crLHS69 - crLHS71 + crLHS85;
const double crLHS87 = crLHS16*crLHS32;
const double crLHS88 = crLHS48*rho;
const double crLHS89 = N[1]*crLHS87 - crLHS77*crLHS88;
const double crLHS90 = DN(0,0)*N[1];
const double crLHS91 = DN(1,0)*crLHS25;
const double crLHS92 = DN(1,0)*crLHS12 + DN(1,1)*crLHS32;
const double crLHS93 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS94 = C(0,2)*DN(2,0);
const double crLHS95 = C(2,2)*DN(2,1) + crLHS94;
const double crLHS96 = N[2]*rho;
const double crLHS97 = crLHS12*crLHS96;
const double crLHS98 = N[2]*crLHS4;
const double crLHS99 = DN(2,0)*crLHS6;
const double crLHS100 = crLHS99*rho;
const double crLHS101 = DN(2,1)*crLHS7;
const double crLHS102 = crLHS101*rho;
const double crLHS103 = crLHS100 + crLHS102 - crLHS98;
const double crLHS104 = crLHS103 - crLHS97;
const double crLHS105 = crLHS104*crLHS25;
const double crLHS106 = crLHS37*crLHS96;
const double crLHS107 = -crLHS104*crLHS12 + crLHS106;
const double crLHS108 = N[2]*crLHS15;
const double crLHS109 = crLHS108 + crLHS27*crLHS96;
const double crLHS110 = DN(0,0)*DN(2,0);
const double crLHS111 = N[2]*crLHS17 + crLHS11*crLHS110;
const double crLHS112 = C(0,1)*DN(2,1) + crLHS94;
const double crLHS113 = C(1,2)*DN(2,1);
const double crLHS114 = C(2,2)*DN(2,0) + crLHS113;
const double crLHS115 = DN(2,1)*crLHS46;
const double crLHS116 = crLHS53*crLHS96;
const double crLHS117 = -crLHS100 - crLHS102 + crLHS116 + crLHS97 + crLHS98;
const double crLHS118 = N[2]*crLHS87 - crLHS108*crLHS88;
const double crLHS119 = DN(0,0)*N[2];
const double crLHS120 = DN(2,0)*crLHS25;
const double crLHS121 = DN(2,0)*crLHS12 + DN(2,1)*crLHS32;
const double crLHS122 = C(0,1)*DN(0,0) + crLHS44;
const double crLHS123 = crLHS25*crLHS36;
const double crLHS124 = crLHS123*crLHS51;
const double crLHS125 = 1.0*crLHS33 + 1.0*crLHS34 + 1.0*crLHS35;
const double crLHS126 = crLHS125*crLHS57;
const double crLHS127 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS128 = pow(DN(0,1), 2);
const double crLHS129 = crLHS22 - crLHS54;
const double crLHS130 = crLHS129*crLHS25;
const double crLHS131 = -crLHS129*crLHS53 + crLHS38;
const double crLHS132 = DN(0,1)*N[0];
const double crLHS133 = DN(0,1)*crLHS25;
const double crLHS134 = DN(0,0)*crLHS36 + DN(0,1)*crLHS53;
const double crLHS135 = C(0,1)*DN(1,0) + crLHS82;
const double crLHS136 = DN(0,1)*crLHS11;
const double crLHS137 = DN(1,0)*crLHS136;
const double crLHS138 = crLHS16*crLHS36;
const double crLHS139 = crLHS123*rho;
const double crLHS140 = N[1]*crLHS138 - crLHS139*crLHS77;
const double crLHS141 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS142 = crLHS72 - crLHS85;
const double crLHS143 = crLHS142*crLHS25;
const double crLHS144 = -crLHS142*crLHS53 + crLHS75;
const double crLHS145 = DN(0,1)*DN(1,1);
const double crLHS146 = N[1]*crLHS54 + crLHS11*crLHS145;
const double crLHS147 = DN(0,1)*N[1];
const double crLHS148 = DN(1,1)*crLHS25;
const double crLHS149 = DN(1,0)*crLHS36 + DN(1,1)*crLHS53;
const double crLHS150 = C(0,1)*DN(2,0) + crLHS113;
const double crLHS151 = DN(2,0)*crLHS136;
const double crLHS152 = N[2]*crLHS138 - crLHS108*crLHS139;
const double crLHS153 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS154 = crLHS103 - crLHS116;
const double crLHS155 = crLHS154*crLHS25;
const double crLHS156 = crLHS106 - crLHS154*crLHS53;
const double crLHS157 = DN(0,1)*DN(2,1);
const double crLHS158 = N[2]*crLHS54 + crLHS11*crLHS157;
const double crLHS159 = DN(0,1)*N[2];
const double crLHS160 = DN(2,1)*crLHS25;
const double crLHS161 = DN(2,0)*crLHS36 + DN(2,1)*crLHS53;
const double crLHS162 = crLHS25*gauss_weight;
const double crLHS163 = DN(1,0)*N[0];
const double crLHS164 = DN(1,1)*N[0];
const double crLHS165 = crLHS162*(crLHS145 + crLHS79);
const double crLHS166 = DN(2,0)*N[0];
const double crLHS167 = DN(2,1)*N[0];
const double crLHS168 = crLHS162*(crLHS110 + crLHS157);
const double crLHS169 = crLHS68 + crLHS70;
const double crLHS170 = crLHS169*rho;
const double crLHS171 = crLHS25*crLHS65;
const double crLHS172 = crLHS16*crLHS169 + crLHS77;
const double crLHS173 = crLHS24*crLHS65;
const double crLHS174 = crLHS173*crLHS56;
const double crLHS175 = crLHS169*crLHS50;
const double crLHS176 = crLHS175*crLHS48;
const double crLHS177 = pow(DN(1,0), 2);
const double crLHS178 = pow(N[1], 2);
const double crLHS179 = crLHS178*rho;
const double crLHS180 = crLHS178*crLHS4;
const double crLHS181 = crLHS169*crLHS65 + crLHS180;
const double crLHS182 = DN(1,0)*crLHS11;
const double crLHS183 = DN(1,1)*crLHS182;
const double crLHS184 = DN(1,0)*N[1];
const double crLHS185 = N[2]*crLHS67;
const double crLHS186 = crLHS169*crLHS96 + crLHS185;
const double crLHS187 = DN(1,0)*DN(2,0);
const double crLHS188 = N[2]*crLHS66 + crLHS11*crLHS187;
const double crLHS189 = DN(2,1)*crLHS182;
const double crLHS190 = N[2]*crLHS65;
const double crLHS191 = crLHS25*crLHS96;
const double crLHS192 = crLHS191*crLHS67;
const double crLHS193 = crLHS190*crLHS32 - crLHS192*crLHS32;
const double crLHS194 = DN(1,0)*N[2];
const double crLHS195 = crLHS125*crLHS173;
const double crLHS196 = crLHS123*crLHS175;
const double crLHS197 = pow(DN(1,1), 2);
const double crLHS198 = DN(1,1)*N[1];
const double crLHS199 = DN(2,0)*crLHS11;
const double crLHS200 = DN(1,1)*crLHS199;
const double crLHS201 = crLHS190*crLHS36 - crLHS192*crLHS36;
const double crLHS202 = DN(1,1)*DN(2,1);
const double crLHS203 = N[2]*crLHS85 + crLHS11*crLHS202;
const double crLHS204 = DN(1,1)*N[2];
const double crLHS205 = DN(2,0)*N[1];
const double crLHS206 = DN(2,1)*N[1];
const double crLHS207 = crLHS162*(crLHS187 + crLHS202);
const double crLHS208 = crLHS101 + crLHS99;
const double crLHS209 = crLHS208*rho;
const double crLHS210 = crLHS108 + crLHS16*crLHS208;
const double crLHS211 = crLHS24*crLHS96;
const double crLHS212 = crLHS211*crLHS56;
const double crLHS213 = crLHS208*crLHS50;
const double crLHS214 = crLHS213*crLHS48;
const double crLHS215 = crLHS185 + crLHS208*crLHS65;
const double crLHS216 = pow(DN(2,0), 2);
const double crLHS217 = pow(N[2], 2);
const double crLHS218 = crLHS217*rho;
const double crLHS219 = crLHS217*crLHS4;
const double crLHS220 = crLHS208*crLHS96 + crLHS219;
const double crLHS221 = DN(2,1)*crLHS199;
const double crLHS222 = DN(2,0)*N[2];
const double crLHS223 = crLHS125*crLHS211;
const double crLHS224 = crLHS123*crLHS213;
const double crLHS225 = pow(DN(2,1), 2);
const double crLHS226 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS11*crLHS3 + crLHS12*crLHS14 + crLHS15*crLHS26 + crLHS26*crLHS28 - crLHS39*crLHS40 + crLHS42);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS45 - N[0]*crLHS52 + crLHS14*crLHS32 + crLHS47 - crLHS48*crLHS49 - crLHS55*crLHS58);
rLHS(0,2)+=-gauss_weight*(crLHS15*crLHS60 + crLHS28*crLHS60 + crLHS40*crLHS61 + crLHS59);
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS64 + crLHS15*crLHS74 + crLHS28*crLHS74 - crLHS40*crLHS76 + crLHS78 + crLHS80);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS81 + DN(0,1)*crLHS83 - N[1]*crLHS52 - crLHS58*crLHS86 + crLHS84 + crLHS89);
rLHS(0,5)+=-gauss_weight*(crLHS15*crLHS91 + crLHS28*crLHS91 + crLHS40*crLHS92 + crLHS90);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS93 + DN(0,1)*crLHS95 + crLHS105*crLHS15 + crLHS105*crLHS28 - crLHS107*crLHS40 + crLHS109 + crLHS111);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS114 - N[2]*crLHS52 + crLHS115 - crLHS117*crLHS58 + crLHS118);
rLHS(0,8)+=-gauss_weight*(crLHS119 + crLHS120*crLHS15 + crLHS120*crLHS28 + crLHS121*crLHS40);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS122 - N[0]*crLHS124 - crLHS123*crLHS49 - crLHS126*crLHS55 + crLHS14*crLHS36 + crLHS47);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS127 + crLHS11*crLHS128 + crLHS130*crLHS15 + crLHS130*crLHS28 - crLHS131*crLHS40 + crLHS14*crLHS53 + crLHS42);
rLHS(1,2)+=-gauss_weight*(crLHS132 + crLHS133*crLHS15 + crLHS133*crLHS28 + crLHS134*crLHS40);
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS64 + DN(0,1)*crLHS135 - N[1]*crLHS124 - crLHS126*crLHS86 + crLHS137 + crLHS140);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS83 + DN(0,1)*crLHS141 + crLHS143*crLHS15 + crLHS143*crLHS28 - crLHS144*crLHS40 + crLHS146 + crLHS78);
rLHS(1,5)+=-gauss_weight*(crLHS147 + crLHS148*crLHS15 + crLHS148*crLHS28 + crLHS149*crLHS40);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS95 + DN(0,1)*crLHS150 - N[2]*crLHS124 - crLHS117*crLHS126 + crLHS151 + crLHS152);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS153 + crLHS109 + crLHS15*crLHS155 + crLHS155*crLHS28 - crLHS156*crLHS40 + crLHS158);
rLHS(1,8)+=-gauss_weight*(crLHS15*crLHS160 + crLHS159 + crLHS160*crLHS28 + crLHS161*crLHS40);
rLHS(2,0)+=gauss_weight*(crLHS132*crLHS139 - crLHS23*crLHS60 + crLHS59);
rLHS(2,1)+=gauss_weight*(-crLHS129*crLHS133 + crLHS132 + crLHS59*crLHS88);
rLHS(2,2)+=crLHS162*(crLHS128 + crLHS3);
rLHS(2,3)+=gauss_weight*(crLHS139*crLHS147 + crLHS163 - crLHS60*crLHS73);
rLHS(2,4)+=gauss_weight*(-crLHS133*crLHS142 + crLHS164 + crLHS88*crLHS90);
rLHS(2,5)+=crLHS165;
rLHS(2,6)+=gauss_weight*(-crLHS104*crLHS60 + crLHS139*crLHS159 + crLHS166);
rLHS(2,7)+=gauss_weight*(crLHS119*crLHS88 - crLHS133*crLHS154 + crLHS167);
rLHS(2,8)+=crLHS168;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS170*crLHS26 - crLHS171*crLHS39 + crLHS172 + crLHS26*crLHS67 + crLHS80);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS45 - N[0]*crLHS176 + crLHS137 - crLHS174*crLHS55 + crLHS89);
rLHS(3,2)+=-gauss_weight*(crLHS163 + crLHS170*crLHS60 + crLHS171*crLHS61 + crLHS60*crLHS67);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS64 + crLHS11*crLHS177 + crLHS12*crLHS179 + crLHS170*crLHS74 - crLHS171*crLHS76 + crLHS181 + crLHS67*crLHS74);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS81 + DN(1,1)*crLHS83 - N[1]*crLHS176 - crLHS174*crLHS86 + crLHS179*crLHS32 - crLHS180*crLHS88 + crLHS183);
rLHS(3,5)+=-gauss_weight*(crLHS170*crLHS91 + crLHS171*crLHS92 + crLHS184 + crLHS67*crLHS91);
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS93 + DN(1,1)*crLHS95 + crLHS105*crLHS170 + crLHS105*crLHS67 - crLHS107*crLHS171 + crLHS186 + crLHS188);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS114 - N[2]*crLHS176 - crLHS117*crLHS174 + crLHS189 + crLHS193);
rLHS(3,8)+=-gauss_weight*(crLHS120*crLHS170 + crLHS120*crLHS67 + crLHS121*crLHS171 + crLHS194);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS122 - N[0]*crLHS196 + crLHS140 - crLHS195*crLHS55 + crLHS84);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS127 + crLHS130*crLHS170 + crLHS130*crLHS67 - crLHS131*crLHS171 + crLHS146 + crLHS172);
rLHS(4,2)+=-gauss_weight*(crLHS133*crLHS170 + crLHS133*crLHS67 + crLHS134*crLHS171 + crLHS164);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS64 + DN(1,1)*crLHS135 - N[1]*crLHS196 - crLHS139*crLHS180 + crLHS179*crLHS36 + crLHS183 - crLHS195*crLHS86);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS83 + DN(1,1)*crLHS141 + crLHS11*crLHS197 + crLHS143*crLHS170 + crLHS143*crLHS67 - crLHS144*crLHS171 + crLHS179*crLHS53 + crLHS181);
rLHS(4,5)+=-gauss_weight*(crLHS148*crLHS170 + crLHS148*crLHS67 + crLHS149*crLHS171 + crLHS198);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS95 + DN(1,1)*crLHS150 - N[2]*crLHS196 - crLHS117*crLHS195 + crLHS200 + crLHS201);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS153 + crLHS155*crLHS170 + crLHS155*crLHS67 - crLHS156*crLHS171 + crLHS186 + crLHS203);
rLHS(4,8)+=-gauss_weight*(crLHS160*crLHS170 + crLHS160*crLHS67 + crLHS161*crLHS171 + crLHS204);
rLHS(5,0)+=gauss_weight*(crLHS139*crLHS164 - crLHS23*crLHS91 + crLHS90);
rLHS(5,1)+=gauss_weight*(-crLHS129*crLHS148 + crLHS147 + crLHS163*crLHS88);
rLHS(5,2)+=crLHS165;
rLHS(5,3)+=gauss_weight*(crLHS139*crLHS198 + crLHS184 - crLHS73*crLHS91);
rLHS(5,4)+=gauss_weight*(-crLHS142*crLHS148 + crLHS184*crLHS88 + crLHS198);
rLHS(5,5)+=crLHS162*(crLHS177 + crLHS197);
rLHS(5,6)+=gauss_weight*(-crLHS104*crLHS91 + crLHS139*crLHS204 + crLHS205);
rLHS(5,7)+=gauss_weight*(-crLHS148*crLHS154 + crLHS194*crLHS88 + crLHS206);
rLHS(5,8)+=crLHS207;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS111 - crLHS191*crLHS39 + crLHS209*crLHS26 + crLHS210 + crLHS26*crLHS98);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS45 - N[0]*crLHS214 + crLHS118 + crLHS151 - crLHS212*crLHS55);
rLHS(6,2)+=-gauss_weight*(crLHS166 + crLHS191*crLHS61 + crLHS209*crLHS60 + crLHS60*crLHS98);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS64 + crLHS188 - crLHS191*crLHS76 + crLHS209*crLHS74 + crLHS215 + crLHS74*crLHS98);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS81 + DN(2,1)*crLHS83 - N[1]*crLHS214 + crLHS193 + crLHS200 - crLHS212*crLHS86);
rLHS(6,5)+=-gauss_weight*(crLHS191*crLHS92 + crLHS205 + crLHS209*crLHS91 + crLHS91*crLHS98);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS93 + DN(2,1)*crLHS95 + crLHS105*crLHS209 + crLHS105*crLHS98 - crLHS107*crLHS191 + crLHS11*crLHS216 + crLHS12*crLHS218 + crLHS220);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS114 - N[2]*crLHS214 - crLHS117*crLHS212 + crLHS218*crLHS32 - crLHS219*crLHS88 + crLHS221);
rLHS(6,8)+=-gauss_weight*(crLHS120*crLHS209 + crLHS120*crLHS98 + crLHS121*crLHS191 + crLHS222);
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS122 - N[0]*crLHS224 + crLHS115 + crLHS152 - crLHS223*crLHS55);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS127 + crLHS130*crLHS209 + crLHS130*crLHS98 - crLHS131*crLHS191 + crLHS158 + crLHS210);
rLHS(7,2)+=-gauss_weight*(crLHS133*crLHS209 + crLHS133*crLHS98 + crLHS134*crLHS191 + crLHS167);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS64 + DN(2,1)*crLHS135 - N[1]*crLHS224 + crLHS189 + crLHS201 - crLHS223*crLHS86);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS83 + DN(2,1)*crLHS141 + crLHS143*crLHS209 + crLHS143*crLHS98 - crLHS144*crLHS191 + crLHS203 + crLHS215);
rLHS(7,5)+=-gauss_weight*(crLHS148*crLHS209 + crLHS148*crLHS98 + crLHS149*crLHS191 + crLHS206);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS95 + DN(2,1)*crLHS150 - N[2]*crLHS224 - crLHS117*crLHS223 - crLHS139*crLHS219 + crLHS218*crLHS36 + crLHS221);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS153 + crLHS11*crLHS225 + crLHS155*crLHS209 + crLHS155*crLHS98 - crLHS156*crLHS191 + crLHS218*crLHS53 + crLHS220);
rLHS(7,8)+=-gauss_weight*(crLHS160*crLHS209 + crLHS160*crLHS98 + crLHS161*crLHS191 + crLHS226);
rLHS(8,0)+=gauss_weight*(crLHS119 - crLHS120*crLHS23 + crLHS139*crLHS167);
rLHS(8,1)+=gauss_weight*(-crLHS129*crLHS160 + crLHS159 + crLHS166*crLHS88);
rLHS(8,2)+=crLHS168;
rLHS(8,3)+=gauss_weight*(-crLHS120*crLHS73 + crLHS139*crLHS206 + crLHS194);
rLHS(8,4)+=gauss_weight*(-crLHS142*crLHS160 + crLHS204 + crLHS205*crLHS88);
rLHS(8,5)+=crLHS207;
rLHS(8,6)+=gauss_weight*(-crLHS104*crLHS120 + crLHS139*crLHS226 + crLHS222);
rLHS(8,7)+=gauss_weight*(-crLHS154*crLHS160 + crLHS222*crLHS88 + crLHS226);
rLHS(8,8)+=crLHS162*(crLHS216 + crLHS225);
    
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
const double crLHS5 = pow(DN(0,0), 2);
const double crLHS6 = N[0]*alpha[0] + N[1]*alpha[1] + N[2]*alpha[2] + N[3]*alpha[3];
const double crLHS7 = crLHS6*stab_c3;
const double crLHS8 = N[0]*v_ns(0,0) + N[1]*v_ns(1,0) + N[2]*v_ns(2,0) + N[3]*v_ns(3,0);
const double crLHS9 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crLHS10 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crLHS11 = pow(crLHS10, 2) + pow(crLHS8, 2) + pow(crLHS9, 2);
const double crLHS12 = pow(crLHS11, 1.0/4.0)*stab_c3;
const double crLHS13 = sqrt(crLHS11)*rho*stab_c2;
const double crLHS14 = h*(crLHS12*h + crLHS13 + crLHS7*h)/stab_c1 + mu;
const double crLHS15 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crLHS16 = pow(N[0], 2);
const double crLHS17 = crLHS16*rho;
const double crLHS18 = N[0]*crLHS6;
const double crLHS19 = N[0]*rho;
const double crLHS20 = crLHS15*crLHS19;
const double crLHS21 = DN(0,0)*crLHS8;
const double crLHS22 = DN(0,1)*crLHS9;
const double crLHS23 = DN(0,2)*crLHS10;
const double crLHS24 = -crLHS18 + crLHS21*rho + crLHS22*rho + crLHS23*rho;
const double crLHS25 = -crLHS20 + crLHS24;
const double crLHS26 = 1.0/(crLHS12 + crLHS13/h + crLHS7 + mu*stab_c1/pow(h, 2));
const double crLHS27 = crLHS25*crLHS26;
const double crLHS28 = crLHS21 + crLHS22 + crLHS23;
const double crLHS29 = crLHS28*rho;
const double crLHS30 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crLHS31 = DN(0,1)*v_ns(0,0) + DN(1,1)*v_ns(1,0) + DN(2,1)*v_ns(2,0) + DN(3,1)*v_ns(3,0);
const double crLHS32 = crLHS19*crLHS31;
const double crLHS33 = crLHS30*crLHS32;
const double crLHS34 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crLHS35 = DN(0,2)*v_ns(0,0) + DN(1,2)*v_ns(1,0) + DN(2,2)*v_ns(2,0) + DN(3,2)*v_ns(3,0);
const double crLHS36 = crLHS19*crLHS35;
const double crLHS37 = crLHS34*crLHS36;
const double crLHS38 = -crLHS15*crLHS25 + crLHS33 + crLHS37;
const double crLHS39 = crLHS19*crLHS26;
const double crLHS40 = crLHS16*crLHS6;
const double crLHS41 = crLHS19*crLHS28 + crLHS40;
const double crLHS42 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crLHS1;
const double crLHS43 = C(1,3)*DN(0,1);
const double crLHS44 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crLHS43;
const double crLHS45 = C(3,5)*DN(0,0);
const double crLHS46 = C(4,5)*DN(0,2);
const double crLHS47 = C(1,5)*DN(0,1) + crLHS45 + crLHS46;
const double crLHS48 = DN(0,0)*crLHS14;
const double crLHS49 = DN(0,1)*crLHS48;
const double crLHS50 = crLHS26*crLHS31;
const double crLHS51 = crLHS40*rho;
const double crLHS52 = pow(rho, 2);
const double crLHS53 = crLHS28*crLHS52;
const double crLHS54 = N[0]*crLHS53;
const double crLHS55 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crLHS56 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crLHS57 = crLHS19*crLHS56;
const double crLHS58 = crLHS24 - crLHS57;
const double crLHS59 = crLHS20*crLHS31 - crLHS31*crLHS58 + crLHS36*crLHS55;
const double crLHS60 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crLHS3;
const double crLHS61 = C(3,4)*DN(0,1);
const double crLHS62 = C(2,3)*DN(0,2) + crLHS45 + crLHS61;
const double crLHS63 = C(2,5)*DN(0,2);
const double crLHS64 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crLHS63;
const double crLHS65 = DN(0,2)*crLHS48;
const double crLHS66 = crLHS26*crLHS35;
const double crLHS67 = DN(0,2)*v_ns(0,1) + DN(1,2)*v_ns(1,1) + DN(2,2)*v_ns(2,1) + DN(3,2)*v_ns(3,1);
const double crLHS68 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crLHS69 = crLHS19*crLHS68;
const double crLHS70 = crLHS24 - crLHS69;
const double crLHS71 = crLHS20*crLHS35 + crLHS32*crLHS67 - crLHS35*crLHS70;
const double crLHS72 = DN(0,0)*N[0];
const double crLHS73 = DN(0,0)*crLHS26;
const double crLHS74 = DN(0,0)*crLHS15 + DN(0,1)*crLHS31 + DN(0,2)*crLHS35;
const double crLHS75 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double crLHS76 = C(0,3)*DN(1,0);
const double crLHS77 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crLHS76;
const double crLHS78 = C(0,5)*DN(1,0);
const double crLHS79 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crLHS78;
const double crLHS80 = N[1]*rho;
const double crLHS81 = crLHS15*crLHS80;
const double crLHS82 = N[1]*crLHS6;
const double crLHS83 = DN(1,0)*crLHS8;
const double crLHS84 = DN(1,1)*crLHS9;
const double crLHS85 = DN(1,2)*crLHS10;
const double crLHS86 = -crLHS82 + crLHS83*rho + crLHS84*rho + crLHS85*rho;
const double crLHS87 = -crLHS81 + crLHS86;
const double crLHS88 = crLHS26*crLHS87;
const double crLHS89 = crLHS31*crLHS80;
const double crLHS90 = crLHS30*crLHS89;
const double crLHS91 = crLHS35*crLHS80;
const double crLHS92 = crLHS34*crLHS91;
const double crLHS93 = -crLHS15*crLHS87 + crLHS90 + crLHS92;
const double crLHS94 = N[1]*crLHS18;
const double crLHS95 = crLHS28*crLHS80 + crLHS94;
const double crLHS96 = DN(0,0)*DN(1,0);
const double crLHS97 = N[1]*crLHS20 + crLHS14*crLHS96;
const double crLHS98 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS76;
const double crLHS99 = C(1,3)*DN(1,1);
const double crLHS100 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS99;
const double crLHS101 = C(3,5)*DN(1,0);
const double crLHS102 = C(4,5)*DN(1,2);
const double crLHS103 = C(1,5)*DN(1,1) + crLHS101 + crLHS102;
const double crLHS104 = DN(1,1)*crLHS48;
const double crLHS105 = crLHS56*crLHS80;
const double crLHS106 = -crLHS105 + crLHS86;
const double crLHS107 = -crLHS106*crLHS31 + crLHS31*crLHS81 + crLHS55*crLHS91;
const double crLHS108 = N[1]*crLHS53;
const double crLHS109 = crLHS94*rho;
const double crLHS110 = N[1]*crLHS32 - crLHS109*crLHS50;
const double crLHS111 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS78;
const double crLHS112 = C(3,4)*DN(1,1);
const double crLHS113 = C(2,3)*DN(1,2) + crLHS101 + crLHS112;
const double crLHS114 = C(2,5)*DN(1,2);
const double crLHS115 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS114;
const double crLHS116 = DN(1,2)*crLHS48;
const double crLHS117 = crLHS68*crLHS80;
const double crLHS118 = -crLHS117 + crLHS86;
const double crLHS119 = -crLHS118*crLHS35 + crLHS35*crLHS81 + crLHS67*crLHS89;
const double crLHS120 = N[1]*crLHS36 - crLHS109*crLHS66;
const double crLHS121 = DN(0,0)*N[1];
const double crLHS122 = DN(1,0)*crLHS26;
const double crLHS123 = DN(1,0)*crLHS15 + DN(1,1)*crLHS31 + DN(1,2)*crLHS35;
const double crLHS124 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS125 = C(0,3)*DN(2,0);
const double crLHS126 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS125;
const double crLHS127 = C(0,5)*DN(2,0);
const double crLHS128 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS127;
const double crLHS129 = N[2]*rho;
const double crLHS130 = crLHS129*crLHS15;
const double crLHS131 = N[2]*crLHS6;
const double crLHS132 = DN(2,0)*crLHS8;
const double crLHS133 = DN(2,1)*crLHS9;
const double crLHS134 = DN(2,2)*crLHS10;
const double crLHS135 = -crLHS131 + crLHS132*rho + crLHS133*rho + crLHS134*rho;
const double crLHS136 = -crLHS130 + crLHS135;
const double crLHS137 = crLHS136*crLHS26;
const double crLHS138 = crLHS129*crLHS31;
const double crLHS139 = crLHS138*crLHS30;
const double crLHS140 = crLHS129*crLHS35;
const double crLHS141 = crLHS140*crLHS34;
const double crLHS142 = -crLHS136*crLHS15 + crLHS139 + crLHS141;
const double crLHS143 = N[2]*crLHS18;
const double crLHS144 = crLHS129*crLHS28 + crLHS143;
const double crLHS145 = DN(0,0)*DN(2,0);
const double crLHS146 = N[2]*crLHS20 + crLHS14*crLHS145;
const double crLHS147 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS125;
const double crLHS148 = C(1,3)*DN(2,1);
const double crLHS149 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS148;
const double crLHS150 = C(3,5)*DN(2,0);
const double crLHS151 = C(4,5)*DN(2,2);
const double crLHS152 = C(1,5)*DN(2,1) + crLHS150 + crLHS151;
const double crLHS153 = DN(2,1)*crLHS48;
const double crLHS154 = crLHS129*crLHS56;
const double crLHS155 = crLHS135 - crLHS154;
const double crLHS156 = crLHS130*crLHS31 + crLHS140*crLHS55 - crLHS155*crLHS31;
const double crLHS157 = N[2]*crLHS53;
const double crLHS158 = crLHS143*rho;
const double crLHS159 = N[2]*crLHS32 - crLHS158*crLHS50;
const double crLHS160 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS127;
const double crLHS161 = C(3,4)*DN(2,1);
const double crLHS162 = C(2,3)*DN(2,2) + crLHS150 + crLHS161;
const double crLHS163 = C(2,5)*DN(2,2);
const double crLHS164 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS163;
const double crLHS165 = DN(2,2)*crLHS48;
const double crLHS166 = crLHS129*crLHS68;
const double crLHS167 = crLHS135 - crLHS166;
const double crLHS168 = crLHS130*crLHS35 + crLHS138*crLHS67 - crLHS167*crLHS35;
const double crLHS169 = N[2]*crLHS36 - crLHS158*crLHS66;
const double crLHS170 = DN(0,0)*N[2];
const double crLHS171 = DN(2,0)*crLHS26;
const double crLHS172 = DN(2,0)*crLHS15 + DN(2,1)*crLHS31 + DN(2,2)*crLHS35;
const double crLHS173 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS174 = C(0,3)*DN(3,0);
const double crLHS175 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS174;
const double crLHS176 = C(0,5)*DN(3,0);
const double crLHS177 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS176;
const double crLHS178 = N[3]*rho;
const double crLHS179 = crLHS15*crLHS178;
const double crLHS180 = N[3]*crLHS6;
const double crLHS181 = DN(3,0)*crLHS8;
const double crLHS182 = DN(3,1)*crLHS9;
const double crLHS183 = DN(3,2)*crLHS10;
const double crLHS184 = -crLHS180 + crLHS181*rho + crLHS182*rho + crLHS183*rho;
const double crLHS185 = -crLHS179 + crLHS184;
const double crLHS186 = crLHS185*crLHS26;
const double crLHS187 = crLHS178*crLHS31;
const double crLHS188 = crLHS187*crLHS30;
const double crLHS189 = crLHS178*crLHS35;
const double crLHS190 = crLHS189*crLHS34;
const double crLHS191 = -crLHS15*crLHS185 + crLHS188 + crLHS190;
const double crLHS192 = N[3]*crLHS18;
const double crLHS193 = crLHS178*crLHS28 + crLHS192;
const double crLHS194 = DN(0,0)*DN(3,0);
const double crLHS195 = N[3]*crLHS20 + crLHS14*crLHS194;
const double crLHS196 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS174;
const double crLHS197 = C(1,3)*DN(3,1);
const double crLHS198 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS197;
const double crLHS199 = C(3,5)*DN(3,0);
const double crLHS200 = C(4,5)*DN(3,2);
const double crLHS201 = C(1,5)*DN(3,1) + crLHS199 + crLHS200;
const double crLHS202 = DN(3,1)*crLHS48;
const double crLHS203 = crLHS178*crLHS56;
const double crLHS204 = crLHS184 - crLHS203;
const double crLHS205 = crLHS179*crLHS31 + crLHS189*crLHS55 - crLHS204*crLHS31;
const double crLHS206 = N[3]*crLHS53;
const double crLHS207 = crLHS192*rho;
const double crLHS208 = N[3]*crLHS32 - crLHS207*crLHS50;
const double crLHS209 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS176;
const double crLHS210 = C(3,4)*DN(3,1);
const double crLHS211 = C(2,3)*DN(3,2) + crLHS199 + crLHS210;
const double crLHS212 = C(2,5)*DN(3,2);
const double crLHS213 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS212;
const double crLHS214 = DN(3,2)*crLHS48;
const double crLHS215 = crLHS178*crLHS68;
const double crLHS216 = crLHS184 - crLHS215;
const double crLHS217 = crLHS179*crLHS35 + crLHS187*crLHS67 - crLHS216*crLHS35;
const double crLHS218 = N[3]*crLHS36 - crLHS207*crLHS66;
const double crLHS219 = DN(0,0)*N[3];
const double crLHS220 = DN(3,0)*crLHS26;
const double crLHS221 = DN(3,0)*crLHS15 + DN(3,1)*crLHS31 + DN(3,2)*crLHS35;
const double crLHS222 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS43;
const double crLHS223 = C(0,4)*DN(0,0) + crLHS46 + crLHS61;
const double crLHS224 = crLHS26*crLHS30;
const double crLHS225 = crLHS19*crLHS67;
const double crLHS226 = crLHS225*crLHS34 - crLHS25*crLHS30 + crLHS30*crLHS57;
const double crLHS227 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS228 = C(1,4)*DN(0,1);
const double crLHS229 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS228;
const double crLHS230 = pow(DN(0,1), 2);
const double crLHS231 = crLHS26*crLHS58;
const double crLHS232 = crLHS225*crLHS55;
const double crLHS233 = crLHS232 + crLHS33 - crLHS56*crLHS58;
const double crLHS234 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS228;
const double crLHS235 = C(2,4)*DN(0,2);
const double crLHS236 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS235;
const double crLHS237 = DN(0,1)*crLHS14;
const double crLHS238 = DN(0,2)*crLHS237;
const double crLHS239 = crLHS26*crLHS67;
const double crLHS240 = crLHS30*crLHS36 + crLHS57*crLHS67 - crLHS67*crLHS70;
const double crLHS241 = DN(0,1)*N[0];
const double crLHS242 = DN(0,1)*crLHS26;
const double crLHS243 = DN(0,0)*crLHS30 + DN(0,1)*crLHS56 + DN(0,2)*crLHS67;
const double crLHS244 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS99;
const double crLHS245 = C(0,4)*DN(1,0) + crLHS102 + crLHS112;
const double crLHS246 = DN(1,0)*crLHS237;
const double crLHS247 = crLHS67*crLHS80;
const double crLHS248 = crLHS105*crLHS30 + crLHS247*crLHS34 - crLHS30*crLHS87;
const double crLHS249 = crLHS19*crLHS30;
const double crLHS250 = N[1]*crLHS249 - crLHS109*crLHS224;
const double crLHS251 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS252 = C(1,4)*DN(1,1);
const double crLHS253 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS252;
const double crLHS254 = crLHS106*crLHS26;
const double crLHS255 = crLHS247*crLHS55;
const double crLHS256 = -crLHS106*crLHS56 + crLHS255 + crLHS90;
const double crLHS257 = DN(0,1)*DN(1,1);
const double crLHS258 = N[1]*crLHS57 + crLHS14*crLHS257;
const double crLHS259 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS252;
const double crLHS260 = C(2,4)*DN(1,2);
const double crLHS261 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS260;
const double crLHS262 = DN(1,2)*crLHS237;
const double crLHS263 = crLHS105*crLHS67 - crLHS118*crLHS67 + crLHS30*crLHS91;
const double crLHS264 = N[1]*crLHS225 - crLHS109*crLHS239;
const double crLHS265 = DN(0,1)*N[1];
const double crLHS266 = DN(1,1)*crLHS26;
const double crLHS267 = DN(1,0)*crLHS30 + DN(1,1)*crLHS56 + DN(1,2)*crLHS67;
const double crLHS268 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS148;
const double crLHS269 = C(0,4)*DN(2,0) + crLHS151 + crLHS161;
const double crLHS270 = DN(2,0)*crLHS237;
const double crLHS271 = crLHS129*crLHS67;
const double crLHS272 = -crLHS136*crLHS30 + crLHS154*crLHS30 + crLHS271*crLHS34;
const double crLHS273 = N[2]*crLHS249 - crLHS158*crLHS224;
const double crLHS274 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS275 = C(1,4)*DN(2,1);
const double crLHS276 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS275;
const double crLHS277 = crLHS155*crLHS26;
const double crLHS278 = crLHS271*crLHS55;
const double crLHS279 = crLHS139 - crLHS155*crLHS56 + crLHS278;
const double crLHS280 = DN(0,1)*DN(2,1);
const double crLHS281 = N[2]*crLHS57 + crLHS14*crLHS280;
const double crLHS282 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS275;
const double crLHS283 = C(2,4)*DN(2,2);
const double crLHS284 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS283;
const double crLHS285 = DN(2,2)*crLHS237;
const double crLHS286 = crLHS140*crLHS30 + crLHS154*crLHS67 - crLHS167*crLHS67;
const double crLHS287 = N[2]*crLHS225 - crLHS158*crLHS239;
const double crLHS288 = DN(0,1)*N[2];
const double crLHS289 = DN(2,1)*crLHS26;
const double crLHS290 = DN(2,0)*crLHS30 + DN(2,1)*crLHS56 + DN(2,2)*crLHS67;
const double crLHS291 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS197;
const double crLHS292 = C(0,4)*DN(3,0) + crLHS200 + crLHS210;
const double crLHS293 = DN(3,0)*crLHS237;
const double crLHS294 = crLHS178*crLHS67;
const double crLHS295 = -crLHS185*crLHS30 + crLHS203*crLHS30 + crLHS294*crLHS34;
const double crLHS296 = N[3]*crLHS249 - crLHS207*crLHS224;
const double crLHS297 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS298 = C(1,4)*DN(3,1);
const double crLHS299 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS298;
const double crLHS300 = crLHS204*crLHS26;
const double crLHS301 = crLHS294*crLHS55;
const double crLHS302 = crLHS188 - crLHS204*crLHS56 + crLHS301;
const double crLHS303 = DN(0,1)*DN(3,1);
const double crLHS304 = N[3]*crLHS57 + crLHS14*crLHS303;
const double crLHS305 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS298;
const double crLHS306 = C(2,4)*DN(3,2);
const double crLHS307 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS306;
const double crLHS308 = DN(3,2)*crLHS237;
const double crLHS309 = crLHS189*crLHS30 + crLHS203*crLHS67 - crLHS216*crLHS67;
const double crLHS310 = N[3]*crLHS225 - crLHS207*crLHS239;
const double crLHS311 = DN(0,1)*N[3];
const double crLHS312 = DN(3,1)*crLHS26;
const double crLHS313 = DN(3,0)*crLHS30 + DN(3,1)*crLHS56 + DN(3,2)*crLHS67;
const double crLHS314 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS63;
const double crLHS315 = crLHS26*crLHS34;
const double crLHS316 = crLHS249*crLHS55 - crLHS25*crLHS34 + crLHS34*crLHS69;
const double crLHS317 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS235;
const double crLHS318 = crLHS26*crLHS55;
const double crLHS319 = crLHS32*crLHS34 - crLHS55*crLHS58 + crLHS55*crLHS69;
const double crLHS320 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS321 = pow(DN(0,2), 2);
const double crLHS322 = crLHS26*crLHS70;
const double crLHS323 = crLHS232 + crLHS37 - crLHS68*crLHS70;
const double crLHS324 = DN(0,2)*N[0];
const double crLHS325 = DN(0,2)*crLHS26;
const double crLHS326 = DN(0,0)*crLHS34 + DN(0,1)*crLHS55 + DN(0,2)*crLHS68;
const double crLHS327 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS114;
const double crLHS328 = DN(0,2)*crLHS14;
const double crLHS329 = DN(1,0)*crLHS328;
const double crLHS330 = crLHS30*crLHS55;
const double crLHS331 = crLHS117*crLHS34 + crLHS330*crLHS80 - crLHS34*crLHS87;
const double crLHS332 = N[1]*crLHS19;
const double crLHS333 = -crLHS109*crLHS315 + crLHS332*crLHS34;
const double crLHS334 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS260;
const double crLHS335 = DN(1,1)*crLHS328;
const double crLHS336 = -crLHS106*crLHS55 + crLHS117*crLHS55 + crLHS34*crLHS89;
const double crLHS337 = -crLHS109*crLHS318 + crLHS332*crLHS55;
const double crLHS338 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS339 = crLHS118*crLHS26;
const double crLHS340 = -crLHS118*crLHS68 + crLHS255 + crLHS92;
const double crLHS341 = DN(0,2)*DN(1,2);
const double crLHS342 = N[1]*crLHS69 + crLHS14*crLHS341;
const double crLHS343 = DN(0,2)*N[1];
const double crLHS344 = DN(1,2)*crLHS26;
const double crLHS345 = DN(1,0)*crLHS34 + DN(1,1)*crLHS55 + DN(1,2)*crLHS68;
const double crLHS346 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS163;
const double crLHS347 = DN(2,0)*crLHS328;
const double crLHS348 = crLHS129*crLHS330 - crLHS136*crLHS34 + crLHS166*crLHS34;
const double crLHS349 = N[2]*crLHS19;
const double crLHS350 = -crLHS158*crLHS315 + crLHS34*crLHS349;
const double crLHS351 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS283;
const double crLHS352 = DN(2,1)*crLHS328;
const double crLHS353 = crLHS138*crLHS34 - crLHS155*crLHS55 + crLHS166*crLHS55;
const double crLHS354 = -crLHS158*crLHS318 + crLHS349*crLHS55;
const double crLHS355 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS356 = crLHS167*crLHS26;
const double crLHS357 = crLHS141 - crLHS167*crLHS68 + crLHS278;
const double crLHS358 = DN(0,2)*DN(2,2);
const double crLHS359 = N[2]*crLHS69 + crLHS14*crLHS358;
const double crLHS360 = DN(0,2)*N[2];
const double crLHS361 = DN(2,2)*crLHS26;
const double crLHS362 = DN(2,0)*crLHS34 + DN(2,1)*crLHS55 + DN(2,2)*crLHS68;
const double crLHS363 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS212;
const double crLHS364 = DN(3,0)*crLHS328;
const double crLHS365 = crLHS178*crLHS330 - crLHS185*crLHS34 + crLHS215*crLHS34;
const double crLHS366 = N[3]*crLHS19;
const double crLHS367 = -crLHS207*crLHS315 + crLHS34*crLHS366;
const double crLHS368 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS306;
const double crLHS369 = DN(3,1)*crLHS328;
const double crLHS370 = crLHS187*crLHS34 - crLHS204*crLHS55 + crLHS215*crLHS55;
const double crLHS371 = -crLHS207*crLHS318 + crLHS366*crLHS55;
const double crLHS372 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS373 = crLHS216*crLHS26;
const double crLHS374 = crLHS190 - crLHS216*crLHS68 + crLHS301;
const double crLHS375 = DN(0,2)*DN(3,2);
const double crLHS376 = N[3]*crLHS69 + crLHS14*crLHS375;
const double crLHS377 = DN(0,2)*N[3];
const double crLHS378 = DN(3,2)*crLHS26;
const double crLHS379 = DN(3,0)*crLHS34 + DN(3,1)*crLHS55 + DN(3,2)*crLHS68;
const double crLHS380 = crLHS241*rho;
const double crLHS381 = crLHS324*rho;
const double crLHS382 = crLHS72*rho;
const double crLHS383 = crLHS26*gauss_weight;
const double crLHS384 = DN(1,0)*N[0];
const double crLHS385 = crLHS265*rho;
const double crLHS386 = crLHS343*rho;
const double crLHS387 = DN(1,1)*N[0];
const double crLHS388 = crLHS121*rho;
const double crLHS389 = DN(1,2)*N[0];
const double crLHS390 = crLHS383*(crLHS257 + crLHS341 + crLHS96);
const double crLHS391 = DN(2,0)*N[0];
const double crLHS392 = crLHS288*rho;
const double crLHS393 = crLHS360*rho;
const double crLHS394 = DN(2,1)*N[0];
const double crLHS395 = crLHS170*rho;
const double crLHS396 = DN(2,2)*N[0];
const double crLHS397 = crLHS383*(crLHS145 + crLHS280 + crLHS358);
const double crLHS398 = DN(3,0)*N[0];
const double crLHS399 = crLHS311*rho;
const double crLHS400 = crLHS377*rho;
const double crLHS401 = DN(3,1)*N[0];
const double crLHS402 = crLHS219*rho;
const double crLHS403 = DN(3,2)*N[0];
const double crLHS404 = crLHS383*(crLHS194 + crLHS303 + crLHS375);
const double crLHS405 = crLHS83 + crLHS84 + crLHS85;
const double crLHS406 = crLHS405*rho;
const double crLHS407 = crLHS26*crLHS80;
const double crLHS408 = crLHS19*crLHS405 + crLHS94;
const double crLHS409 = crLHS405*crLHS52;
const double crLHS410 = N[0]*crLHS409;
const double crLHS411 = pow(DN(1,0), 2);
const double crLHS412 = pow(N[1], 2);
const double crLHS413 = crLHS412*rho;
const double crLHS414 = crLHS412*crLHS6;
const double crLHS415 = crLHS405*crLHS80 + crLHS414;
const double crLHS416 = DN(1,0)*crLHS14;
const double crLHS417 = DN(1,1)*crLHS416;
const double crLHS418 = crLHS414*rho;
const double crLHS419 = N[1]*crLHS409;
const double crLHS420 = DN(1,2)*crLHS416;
const double crLHS421 = DN(1,0)*N[1];
const double crLHS422 = N[2]*crLHS82;
const double crLHS423 = crLHS129*crLHS405 + crLHS422;
const double crLHS424 = DN(1,0)*DN(2,0);
const double crLHS425 = N[2]*crLHS81 + crLHS14*crLHS424;
const double crLHS426 = DN(2,1)*crLHS416;
const double crLHS427 = N[2]*crLHS409;
const double crLHS428 = crLHS26*crLHS82;
const double crLHS429 = N[2]*crLHS89 - crLHS138*crLHS428;
const double crLHS430 = DN(2,2)*crLHS416;
const double crLHS431 = N[2]*crLHS91 - crLHS140*crLHS428;
const double crLHS432 = DN(1,0)*N[2];
const double crLHS433 = N[3]*crLHS82;
const double crLHS434 = crLHS178*crLHS405 + crLHS433;
const double crLHS435 = DN(1,0)*DN(3,0);
const double crLHS436 = N[3]*crLHS81 + crLHS14*crLHS435;
const double crLHS437 = DN(3,1)*crLHS416;
const double crLHS438 = N[3]*crLHS409;
const double crLHS439 = N[3]*crLHS89 - crLHS187*crLHS428;
const double crLHS440 = DN(3,2)*crLHS416;
const double crLHS441 = N[3]*crLHS91 - crLHS189*crLHS428;
const double crLHS442 = DN(1,0)*N[3];
const double crLHS443 = pow(DN(1,1), 2);
const double crLHS444 = DN(1,1)*crLHS14;
const double crLHS445 = DN(1,2)*crLHS444;
const double crLHS446 = DN(1,1)*N[1];
const double crLHS447 = DN(2,0)*crLHS444;
const double crLHS448 = crLHS30*crLHS80;
const double crLHS449 = crLHS129*crLHS26;
const double crLHS450 = crLHS30*crLHS82;
const double crLHS451 = N[2]*crLHS448 - crLHS449*crLHS450;
const double crLHS452 = DN(1,1)*DN(2,1);
const double crLHS453 = N[2]*crLHS105 + crLHS14*crLHS452;
const double crLHS454 = DN(2,2)*crLHS444;
const double crLHS455 = N[2]*crLHS247 - crLHS271*crLHS428;
const double crLHS456 = DN(1,1)*N[2];
const double crLHS457 = DN(3,0)*crLHS444;
const double crLHS458 = crLHS178*crLHS26;
const double crLHS459 = N[3]*crLHS448 - crLHS450*crLHS458;
const double crLHS460 = DN(1,1)*DN(3,1);
const double crLHS461 = N[3]*crLHS105 + crLHS14*crLHS460;
const double crLHS462 = DN(3,2)*crLHS444;
const double crLHS463 = N[3]*crLHS247 - crLHS294*crLHS428;
const double crLHS464 = DN(1,1)*N[3];
const double crLHS465 = pow(DN(1,2), 2);
const double crLHS466 = DN(1,2)*N[1];
const double crLHS467 = DN(1,2)*crLHS14;
const double crLHS468 = DN(2,0)*crLHS467;
const double crLHS469 = N[2]*crLHS80;
const double crLHS470 = crLHS449*crLHS82;
const double crLHS471 = crLHS34*crLHS469 - crLHS34*crLHS470;
const double crLHS472 = DN(2,1)*crLHS467;
const double crLHS473 = crLHS469*crLHS55 - crLHS470*crLHS55;
const double crLHS474 = DN(1,2)*DN(2,2);
const double crLHS475 = N[2]*crLHS117 + crLHS14*crLHS474;
const double crLHS476 = DN(1,2)*N[2];
const double crLHS477 = DN(3,0)*crLHS467;
const double crLHS478 = N[3]*crLHS80;
const double crLHS479 = crLHS458*crLHS82;
const double crLHS480 = crLHS34*crLHS478 - crLHS34*crLHS479;
const double crLHS481 = DN(3,1)*crLHS467;
const double crLHS482 = crLHS478*crLHS55 - crLHS479*crLHS55;
const double crLHS483 = DN(1,2)*DN(3,2);
const double crLHS484 = N[3]*crLHS117 + crLHS14*crLHS483;
const double crLHS485 = DN(1,2)*N[3];
const double crLHS486 = crLHS387*rho;
const double crLHS487 = crLHS389*rho;
const double crLHS488 = crLHS384*rho;
const double crLHS489 = crLHS446*rho;
const double crLHS490 = crLHS466*rho;
const double crLHS491 = crLHS421*rho;
const double crLHS492 = DN(2,0)*N[1];
const double crLHS493 = crLHS456*rho;
const double crLHS494 = crLHS476*rho;
const double crLHS495 = DN(2,1)*N[1];
const double crLHS496 = crLHS432*rho;
const double crLHS497 = DN(2,2)*N[1];
const double crLHS498 = crLHS383*(crLHS424 + crLHS452 + crLHS474);
const double crLHS499 = DN(3,0)*N[1];
const double crLHS500 = crLHS464*rho;
const double crLHS501 = crLHS485*rho;
const double crLHS502 = DN(3,1)*N[1];
const double crLHS503 = crLHS442*rho;
const double crLHS504 = DN(3,2)*N[1];
const double crLHS505 = crLHS383*(crLHS435 + crLHS460 + crLHS483);
const double crLHS506 = crLHS132 + crLHS133 + crLHS134;
const double crLHS507 = crLHS506*rho;
const double crLHS508 = crLHS143 + crLHS19*crLHS506;
const double crLHS509 = crLHS506*crLHS52;
const double crLHS510 = N[0]*crLHS509;
const double crLHS511 = crLHS422 + crLHS506*crLHS80;
const double crLHS512 = N[1]*crLHS509;
const double crLHS513 = pow(DN(2,0), 2);
const double crLHS514 = pow(N[2], 2);
const double crLHS515 = crLHS514*rho;
const double crLHS516 = crLHS514*crLHS6;
const double crLHS517 = crLHS129*crLHS506 + crLHS516;
const double crLHS518 = DN(2,0)*crLHS14;
const double crLHS519 = DN(2,1)*crLHS518;
const double crLHS520 = crLHS516*rho;
const double crLHS521 = N[2]*crLHS509;
const double crLHS522 = DN(2,2)*crLHS518;
const double crLHS523 = DN(2,0)*N[2];
const double crLHS524 = N[3]*crLHS131;
const double crLHS525 = crLHS178*crLHS506 + crLHS524;
const double crLHS526 = DN(2,0)*DN(3,0);
const double crLHS527 = N[3]*crLHS130 + crLHS14*crLHS526;
const double crLHS528 = DN(3,1)*crLHS518;
const double crLHS529 = N[3]*crLHS509;
const double crLHS530 = crLHS131*crLHS26;
const double crLHS531 = N[3]*crLHS138 - crLHS187*crLHS530;
const double crLHS532 = DN(3,2)*crLHS518;
const double crLHS533 = N[3]*crLHS140 - crLHS189*crLHS530;
const double crLHS534 = DN(2,0)*N[3];
const double crLHS535 = pow(DN(2,1), 2);
const double crLHS536 = DN(2,1)*crLHS14;
const double crLHS537 = DN(2,2)*crLHS536;
const double crLHS538 = DN(2,1)*N[2];
const double crLHS539 = DN(3,0)*crLHS536;
const double crLHS540 = N[3]*crLHS129;
const double crLHS541 = crLHS131*crLHS458;
const double crLHS542 = crLHS30*crLHS540 - crLHS30*crLHS541;
const double crLHS543 = DN(2,1)*DN(3,1);
const double crLHS544 = N[3]*crLHS154 + crLHS14*crLHS543;
const double crLHS545 = DN(3,2)*crLHS536;
const double crLHS546 = N[3]*crLHS271 - crLHS294*crLHS530;
const double crLHS547 = DN(2,1)*N[3];
const double crLHS548 = pow(DN(2,2), 2);
const double crLHS549 = DN(2,2)*N[2];
const double crLHS550 = DN(2,2)*crLHS14;
const double crLHS551 = DN(3,0)*crLHS550;
const double crLHS552 = crLHS34*crLHS540 - crLHS34*crLHS541;
const double crLHS553 = DN(3,1)*crLHS550;
const double crLHS554 = crLHS540*crLHS55 - crLHS541*crLHS55;
const double crLHS555 = DN(2,2)*DN(3,2);
const double crLHS556 = N[3]*crLHS166 + crLHS14*crLHS555;
const double crLHS557 = DN(2,2)*N[3];
const double crLHS558 = crLHS394*rho;
const double crLHS559 = crLHS396*rho;
const double crLHS560 = crLHS391*rho;
const double crLHS561 = crLHS495*rho;
const double crLHS562 = crLHS497*rho;
const double crLHS563 = crLHS492*rho;
const double crLHS564 = crLHS538*rho;
const double crLHS565 = crLHS549*rho;
const double crLHS566 = crLHS523*rho;
const double crLHS567 = DN(3,0)*N[2];
const double crLHS568 = crLHS547*rho;
const double crLHS569 = crLHS557*rho;
const double crLHS570 = DN(3,1)*N[2];
const double crLHS571 = crLHS534*rho;
const double crLHS572 = DN(3,2)*N[2];
const double crLHS573 = crLHS383*(crLHS526 + crLHS543 + crLHS555);
const double crLHS574 = crLHS181 + crLHS182 + crLHS183;
const double crLHS575 = crLHS574*rho;
const double crLHS576 = crLHS19*crLHS574 + crLHS192;
const double crLHS577 = crLHS52*crLHS574;
const double crLHS578 = N[0]*crLHS577;
const double crLHS579 = crLHS433 + crLHS574*crLHS80;
const double crLHS580 = N[1]*crLHS577;
const double crLHS581 = crLHS129*crLHS574 + crLHS524;
const double crLHS582 = N[2]*crLHS577;
const double crLHS583 = pow(DN(3,0), 2);
const double crLHS584 = pow(N[3], 2);
const double crLHS585 = crLHS584*rho;
const double crLHS586 = crLHS584*crLHS6;
const double crLHS587 = crLHS178*crLHS574 + crLHS586;
const double crLHS588 = DN(3,0)*crLHS14;
const double crLHS589 = DN(3,1)*crLHS588;
const double crLHS590 = crLHS586*rho;
const double crLHS591 = N[3]*crLHS577;
const double crLHS592 = DN(3,2)*crLHS588;
const double crLHS593 = DN(3,0)*N[3];
const double crLHS594 = pow(DN(3,1), 2);
const double crLHS595 = DN(3,1)*DN(3,2)*crLHS14;
const double crLHS596 = DN(3,1)*N[3];
const double crLHS597 = pow(DN(3,2), 2);
const double crLHS598 = DN(3,2)*N[3];
const double crLHS599 = crLHS401*rho;
const double crLHS600 = crLHS403*rho;
const double crLHS601 = crLHS398*rho;
const double crLHS602 = crLHS502*rho;
const double crLHS603 = crLHS504*rho;
const double crLHS604 = crLHS499*rho;
const double crLHS605 = crLHS570*rho;
const double crLHS606 = crLHS572*rho;
const double crLHS607 = crLHS567*rho;
const double crLHS608 = crLHS596*rho;
const double crLHS609 = crLHS598*rho;
const double crLHS610 = crLHS593*rho;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS14*crLHS5 + crLHS15*crLHS17 + crLHS18*crLHS27 + crLHS27*crLHS29 - crLHS38*crLHS39 + crLHS41);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS42 + DN(0,1)*crLHS44 + DN(0,2)*crLHS47 + crLHS17*crLHS31 - crLHS39*crLHS59 + crLHS49 - crLHS50*crLHS51 - crLHS50*crLHS54);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS60 + DN(0,1)*crLHS62 + DN(0,2)*crLHS64 + crLHS17*crLHS35 - crLHS39*crLHS71 - crLHS51*crLHS66 - crLHS54*crLHS66 + crLHS65);
rLHS(0,3)+=-gauss_weight*(crLHS18*crLHS73 + crLHS29*crLHS73 + crLHS39*crLHS74 + crLHS72);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS75 + DN(0,1)*crLHS77 + DN(0,2)*crLHS79 + crLHS18*crLHS88 + crLHS29*crLHS88 - crLHS39*crLHS93 + crLHS95 + crLHS97);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS98 + DN(0,1)*crLHS100 + DN(0,2)*crLHS103 + crLHS104 - crLHS107*crLHS39 - crLHS108*crLHS50 + crLHS110);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS111 + DN(0,1)*crLHS113 + DN(0,2)*crLHS115 - crLHS108*crLHS66 + crLHS116 - crLHS119*crLHS39 + crLHS120);
rLHS(0,7)+=-gauss_weight*(crLHS121 + crLHS122*crLHS18 + crLHS122*crLHS29 + crLHS123*crLHS39);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS124 + DN(0,1)*crLHS126 + DN(0,2)*crLHS128 + crLHS137*crLHS18 + crLHS137*crLHS29 - crLHS142*crLHS39 + crLHS144 + crLHS146);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS147 + DN(0,1)*crLHS149 + DN(0,2)*crLHS152 + crLHS153 - crLHS156*crLHS39 - crLHS157*crLHS50 + crLHS159);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS160 + DN(0,1)*crLHS162 + DN(0,2)*crLHS164 - crLHS157*crLHS66 + crLHS165 - crLHS168*crLHS39 + crLHS169);
rLHS(0,11)+=-gauss_weight*(crLHS170 + crLHS171*crLHS18 + crLHS171*crLHS29 + crLHS172*crLHS39);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS173 + DN(0,1)*crLHS175 + DN(0,2)*crLHS177 + crLHS18*crLHS186 + crLHS186*crLHS29 - crLHS191*crLHS39 + crLHS193 + crLHS195);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS196 + DN(0,1)*crLHS198 + DN(0,2)*crLHS201 + crLHS202 - crLHS205*crLHS39 - crLHS206*crLHS50 + crLHS208);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS209 + DN(0,1)*crLHS211 + DN(0,2)*crLHS213 - crLHS206*crLHS66 + crLHS214 - crLHS217*crLHS39 + crLHS218);
rLHS(0,15)+=-gauss_weight*(crLHS18*crLHS220 + crLHS219 + crLHS220*crLHS29 + crLHS221*crLHS39);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS222 + DN(0,2)*crLHS223 + crLHS17*crLHS30 - crLHS224*crLHS51 - crLHS224*crLHS54 - crLHS226*crLHS39 + crLHS49);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS44 + DN(0,1)*crLHS227 + DN(0,2)*crLHS229 + crLHS14*crLHS230 + crLHS17*crLHS56 + crLHS18*crLHS231 + crLHS231*crLHS29 - crLHS233*crLHS39 + crLHS41);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS62 + DN(0,1)*crLHS234 + DN(0,2)*crLHS236 + crLHS17*crLHS67 + crLHS238 - crLHS239*crLHS51 - crLHS239*crLHS54 - crLHS240*crLHS39);
rLHS(1,3)+=-gauss_weight*(crLHS18*crLHS242 + crLHS241 + crLHS242*crLHS29 + crLHS243*crLHS39);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS77 + DN(0,1)*crLHS244 + DN(0,2)*crLHS245 - crLHS108*crLHS224 + crLHS246 - crLHS248*crLHS39 + crLHS250);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS100 + DN(0,1)*crLHS251 + DN(0,2)*crLHS253 + crLHS18*crLHS254 + crLHS254*crLHS29 - crLHS256*crLHS39 + crLHS258 + crLHS95);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS113 + DN(0,1)*crLHS259 + DN(0,2)*crLHS261 - crLHS108*crLHS239 + crLHS262 - crLHS263*crLHS39 + crLHS264);
rLHS(1,7)+=-gauss_weight*(crLHS18*crLHS266 + crLHS265 + crLHS266*crLHS29 + crLHS267*crLHS39);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS126 + DN(0,1)*crLHS268 + DN(0,2)*crLHS269 - crLHS157*crLHS224 + crLHS270 - crLHS272*crLHS39 + crLHS273);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS149 + DN(0,1)*crLHS274 + DN(0,2)*crLHS276 + crLHS144 + crLHS18*crLHS277 + crLHS277*crLHS29 - crLHS279*crLHS39 + crLHS281);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS162 + DN(0,1)*crLHS282 + DN(0,2)*crLHS284 - crLHS157*crLHS239 + crLHS285 - crLHS286*crLHS39 + crLHS287);
rLHS(1,11)+=-gauss_weight*(crLHS18*crLHS289 + crLHS288 + crLHS289*crLHS29 + crLHS290*crLHS39);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS175 + DN(0,1)*crLHS291 + DN(0,2)*crLHS292 - crLHS206*crLHS224 + crLHS293 - crLHS295*crLHS39 + crLHS296);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS198 + DN(0,1)*crLHS297 + DN(0,2)*crLHS299 + crLHS18*crLHS300 + crLHS193 + crLHS29*crLHS300 - crLHS302*crLHS39 + crLHS304);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS211 + DN(0,1)*crLHS305 + DN(0,2)*crLHS307 - crLHS206*crLHS239 + crLHS308 - crLHS309*crLHS39 + crLHS310);
rLHS(1,15)+=-gauss_weight*(crLHS18*crLHS312 + crLHS29*crLHS312 + crLHS311 + crLHS313*crLHS39);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS223 + DN(0,2)*crLHS314 + crLHS17*crLHS34 - crLHS315*crLHS51 - crLHS315*crLHS54 - crLHS316*crLHS39 + crLHS65);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS47 + DN(0,1)*crLHS229 + DN(0,2)*crLHS317 + crLHS17*crLHS55 + crLHS238 - crLHS318*crLHS51 - crLHS318*crLHS54 - crLHS319*crLHS39);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS64 + DN(0,1)*crLHS236 + DN(0,2)*crLHS320 + crLHS14*crLHS321 + crLHS17*crLHS68 + crLHS18*crLHS322 + crLHS29*crLHS322 - crLHS323*crLHS39 + crLHS41);
rLHS(2,3)+=-gauss_weight*(crLHS18*crLHS325 + crLHS29*crLHS325 + crLHS324 + crLHS326*crLHS39);
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS79 + DN(0,1)*crLHS245 + DN(0,2)*crLHS327 - crLHS108*crLHS315 + crLHS329 - crLHS331*crLHS39 + crLHS333);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS103 + DN(0,1)*crLHS253 + DN(0,2)*crLHS334 - crLHS108*crLHS318 + crLHS335 - crLHS336*crLHS39 + crLHS337);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS115 + DN(0,1)*crLHS261 + DN(0,2)*crLHS338 + crLHS18*crLHS339 + crLHS29*crLHS339 - crLHS340*crLHS39 + crLHS342 + crLHS95);
rLHS(2,7)+=-gauss_weight*(crLHS18*crLHS344 + crLHS29*crLHS344 + crLHS343 + crLHS345*crLHS39);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS128 + DN(0,1)*crLHS269 + DN(0,2)*crLHS346 - crLHS157*crLHS315 + crLHS347 - crLHS348*crLHS39 + crLHS350);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS152 + DN(0,1)*crLHS276 + DN(0,2)*crLHS351 - crLHS157*crLHS318 + crLHS352 - crLHS353*crLHS39 + crLHS354);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS164 + DN(0,1)*crLHS284 + DN(0,2)*crLHS355 + crLHS144 + crLHS18*crLHS356 + crLHS29*crLHS356 - crLHS357*crLHS39 + crLHS359);
rLHS(2,11)+=-gauss_weight*(crLHS18*crLHS361 + crLHS29*crLHS361 + crLHS360 + crLHS362*crLHS39);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS177 + DN(0,1)*crLHS292 + DN(0,2)*crLHS363 - crLHS206*crLHS315 + crLHS364 - crLHS365*crLHS39 + crLHS367);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS201 + DN(0,1)*crLHS299 + DN(0,2)*crLHS368 - crLHS206*crLHS318 + crLHS369 - crLHS370*crLHS39 + crLHS371);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS213 + DN(0,1)*crLHS307 + DN(0,2)*crLHS372 + crLHS18*crLHS373 + crLHS193 + crLHS29*crLHS373 - crLHS374*crLHS39 + crLHS376);
rLHS(2,15)+=-gauss_weight*(crLHS18*crLHS378 + crLHS29*crLHS378 + crLHS377 + crLHS379*crLHS39);
rLHS(3,0)+=gauss_weight*(crLHS224*crLHS380 - crLHS25*crLHS73 + crLHS315*crLHS381 + crLHS72);
rLHS(3,1)+=gauss_weight*(crLHS241 - crLHS242*crLHS58 + crLHS318*crLHS381 + crLHS382*crLHS50);
rLHS(3,2)+=gauss_weight*(crLHS239*crLHS380 + crLHS324 - crLHS325*crLHS70 + crLHS382*crLHS66);
rLHS(3,3)+=crLHS383*(crLHS230 + crLHS321 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS224*crLHS385 + crLHS315*crLHS386 + crLHS384 - crLHS73*crLHS87);
rLHS(3,5)+=gauss_weight*(-crLHS106*crLHS242 + crLHS318*crLHS386 + crLHS387 + crLHS388*crLHS50);
rLHS(3,6)+=gauss_weight*(-crLHS118*crLHS325 + crLHS239*crLHS385 + crLHS388*crLHS66 + crLHS389);
rLHS(3,7)+=crLHS390;
rLHS(3,8)+=gauss_weight*(-crLHS136*crLHS73 + crLHS224*crLHS392 + crLHS315*crLHS393 + crLHS391);
rLHS(3,9)+=gauss_weight*(-crLHS155*crLHS242 + crLHS318*crLHS393 + crLHS394 + crLHS395*crLHS50);
rLHS(3,10)+=gauss_weight*(-crLHS167*crLHS325 + crLHS239*crLHS392 + crLHS395*crLHS66 + crLHS396);
rLHS(3,11)+=crLHS397;
rLHS(3,12)+=gauss_weight*(-crLHS185*crLHS73 + crLHS224*crLHS399 + crLHS315*crLHS400 + crLHS398);
rLHS(3,13)+=gauss_weight*(-crLHS204*crLHS242 + crLHS318*crLHS400 + crLHS401 + crLHS402*crLHS50);
rLHS(3,14)+=gauss_weight*(-crLHS216*crLHS325 + crLHS239*crLHS399 + crLHS402*crLHS66 + crLHS403);
rLHS(3,15)+=crLHS404;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS27*crLHS406 + crLHS27*crLHS82 - crLHS38*crLHS407 + crLHS408 + crLHS97);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS42 + DN(1,1)*crLHS44 + DN(1,2)*crLHS47 + crLHS110 + crLHS246 - crLHS407*crLHS59 - crLHS410*crLHS50);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS60 + DN(1,1)*crLHS62 + DN(1,2)*crLHS64 + crLHS120 + crLHS329 - crLHS407*crLHS71 - crLHS410*crLHS66);
rLHS(4,3)+=-gauss_weight*(crLHS384 + crLHS406*crLHS73 + crLHS407*crLHS74 + crLHS73*crLHS82);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS75 + DN(1,1)*crLHS77 + DN(1,2)*crLHS79 + crLHS14*crLHS411 + crLHS15*crLHS413 + crLHS406*crLHS88 - crLHS407*crLHS93 + crLHS415 + crLHS82*crLHS88);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS98 + DN(1,1)*crLHS100 + DN(1,2)*crLHS103 - crLHS107*crLHS407 + crLHS31*crLHS413 + crLHS417 - crLHS418*crLHS50 - crLHS419*crLHS50);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS111 + DN(1,1)*crLHS113 + DN(1,2)*crLHS115 - crLHS119*crLHS407 + crLHS35*crLHS413 - crLHS418*crLHS66 - crLHS419*crLHS66 + crLHS420);
rLHS(4,7)+=-gauss_weight*(crLHS122*crLHS406 + crLHS122*crLHS82 + crLHS123*crLHS407 + crLHS421);
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS124 + DN(1,1)*crLHS126 + DN(1,2)*crLHS128 + crLHS137*crLHS406 + crLHS137*crLHS82 - crLHS142*crLHS407 + crLHS423 + crLHS425);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS147 + DN(1,1)*crLHS149 + DN(1,2)*crLHS152 - crLHS156*crLHS407 + crLHS426 - crLHS427*crLHS50 + crLHS429);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS160 + DN(1,1)*crLHS162 + DN(1,2)*crLHS164 - crLHS168*crLHS407 - crLHS427*crLHS66 + crLHS430 + crLHS431);
rLHS(4,11)+=-gauss_weight*(crLHS171*crLHS406 + crLHS171*crLHS82 + crLHS172*crLHS407 + crLHS432);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS173 + DN(1,1)*crLHS175 + DN(1,2)*crLHS177 + crLHS186*crLHS406 + crLHS186*crLHS82 - crLHS191*crLHS407 + crLHS434 + crLHS436);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS196 + DN(1,1)*crLHS198 + DN(1,2)*crLHS201 - crLHS205*crLHS407 + crLHS437 - crLHS438*crLHS50 + crLHS439);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS209 + DN(1,1)*crLHS211 + DN(1,2)*crLHS213 - crLHS217*crLHS407 - crLHS438*crLHS66 + crLHS440 + crLHS441);
rLHS(4,15)+=-gauss_weight*(crLHS220*crLHS406 + crLHS220*crLHS82 + crLHS221*crLHS407 + crLHS442);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS222 + DN(1,2)*crLHS223 + crLHS104 - crLHS224*crLHS410 - crLHS226*crLHS407 + crLHS250);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS44 + DN(1,1)*crLHS227 + DN(1,2)*crLHS229 + crLHS231*crLHS406 + crLHS231*crLHS82 - crLHS233*crLHS407 + crLHS258 + crLHS408);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS62 + DN(1,1)*crLHS234 + DN(1,2)*crLHS236 - crLHS239*crLHS410 - crLHS240*crLHS407 + crLHS264 + crLHS335);
rLHS(5,3)+=-gauss_weight*(crLHS242*crLHS406 + crLHS242*crLHS82 + crLHS243*crLHS407 + crLHS387);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS77 + DN(1,1)*crLHS244 + DN(1,2)*crLHS245 - crLHS224*crLHS418 - crLHS224*crLHS419 - crLHS248*crLHS407 + crLHS30*crLHS413 + crLHS417);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS100 + DN(1,1)*crLHS251 + DN(1,2)*crLHS253 + crLHS14*crLHS443 + crLHS254*crLHS406 + crLHS254*crLHS82 - crLHS256*crLHS407 + crLHS413*crLHS56 + crLHS415);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS113 + DN(1,1)*crLHS259 + DN(1,2)*crLHS261 - crLHS239*crLHS418 - crLHS239*crLHS419 - crLHS263*crLHS407 + crLHS413*crLHS67 + crLHS445);
rLHS(5,7)+=-gauss_weight*(crLHS266*crLHS406 + crLHS266*crLHS82 + crLHS267*crLHS407 + crLHS446);
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS126 + DN(1,1)*crLHS268 + DN(1,2)*crLHS269 - crLHS224*crLHS427 - crLHS272*crLHS407 + crLHS447 + crLHS451);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS149 + DN(1,1)*crLHS274 + DN(1,2)*crLHS276 + crLHS277*crLHS406 + crLHS277*crLHS82 - crLHS279*crLHS407 + crLHS423 + crLHS453);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS162 + DN(1,1)*crLHS282 + DN(1,2)*crLHS284 - crLHS239*crLHS427 - crLHS286*crLHS407 + crLHS454 + crLHS455);
rLHS(5,11)+=-gauss_weight*(crLHS289*crLHS406 + crLHS289*crLHS82 + crLHS290*crLHS407 + crLHS456);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS175 + DN(1,1)*crLHS291 + DN(1,2)*crLHS292 - crLHS224*crLHS438 - crLHS295*crLHS407 + crLHS457 + crLHS459);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS198 + DN(1,1)*crLHS297 + DN(1,2)*crLHS299 + crLHS300*crLHS406 + crLHS300*crLHS82 - crLHS302*crLHS407 + crLHS434 + crLHS461);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS211 + DN(1,1)*crLHS305 + DN(1,2)*crLHS307 - crLHS239*crLHS438 - crLHS309*crLHS407 + crLHS462 + crLHS463);
rLHS(5,15)+=-gauss_weight*(crLHS312*crLHS406 + crLHS312*crLHS82 + crLHS313*crLHS407 + crLHS464);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS223 + DN(1,2)*crLHS314 + crLHS116 - crLHS315*crLHS410 - crLHS316*crLHS407 + crLHS333);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS47 + DN(1,1)*crLHS229 + DN(1,2)*crLHS317 + crLHS262 - crLHS318*crLHS410 - crLHS319*crLHS407 + crLHS337);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS64 + DN(1,1)*crLHS236 + DN(1,2)*crLHS320 + crLHS322*crLHS406 + crLHS322*crLHS82 - crLHS323*crLHS407 + crLHS342 + crLHS408);
rLHS(6,3)+=-gauss_weight*(crLHS325*crLHS406 + crLHS325*crLHS82 + crLHS326*crLHS407 + crLHS389);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS79 + DN(1,1)*crLHS245 + DN(1,2)*crLHS327 - crLHS315*crLHS418 - crLHS315*crLHS419 - crLHS331*crLHS407 + crLHS34*crLHS413 + crLHS420);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS103 + DN(1,1)*crLHS253 + DN(1,2)*crLHS334 - crLHS318*crLHS418 - crLHS318*crLHS419 - crLHS336*crLHS407 + crLHS413*crLHS55 + crLHS445);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS115 + DN(1,1)*crLHS261 + DN(1,2)*crLHS338 + crLHS14*crLHS465 + crLHS339*crLHS406 + crLHS339*crLHS82 - crLHS340*crLHS407 + crLHS413*crLHS68 + crLHS415);
rLHS(6,7)+=-gauss_weight*(crLHS344*crLHS406 + crLHS344*crLHS82 + crLHS345*crLHS407 + crLHS466);
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS128 + DN(1,1)*crLHS269 + DN(1,2)*crLHS346 - crLHS315*crLHS427 - crLHS348*crLHS407 + crLHS468 + crLHS471);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS152 + DN(1,1)*crLHS276 + DN(1,2)*crLHS351 - crLHS318*crLHS427 - crLHS353*crLHS407 + crLHS472 + crLHS473);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS164 + DN(1,1)*crLHS284 + DN(1,2)*crLHS355 + crLHS356*crLHS406 + crLHS356*crLHS82 - crLHS357*crLHS407 + crLHS423 + crLHS475);
rLHS(6,11)+=-gauss_weight*(crLHS361*crLHS406 + crLHS361*crLHS82 + crLHS362*crLHS407 + crLHS476);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS177 + DN(1,1)*crLHS292 + DN(1,2)*crLHS363 - crLHS315*crLHS438 - crLHS365*crLHS407 + crLHS477 + crLHS480);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS201 + DN(1,1)*crLHS299 + DN(1,2)*crLHS368 - crLHS318*crLHS438 - crLHS370*crLHS407 + crLHS481 + crLHS482);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS213 + DN(1,1)*crLHS307 + DN(1,2)*crLHS372 + crLHS373*crLHS406 + crLHS373*crLHS82 - crLHS374*crLHS407 + crLHS434 + crLHS484);
rLHS(6,15)+=-gauss_weight*(crLHS378*crLHS406 + crLHS378*crLHS82 + crLHS379*crLHS407 + crLHS485);
rLHS(7,0)+=gauss_weight*(crLHS121 - crLHS122*crLHS25 + crLHS224*crLHS486 + crLHS315*crLHS487);
rLHS(7,1)+=gauss_weight*(crLHS265 - crLHS266*crLHS58 + crLHS318*crLHS487 + crLHS488*crLHS50);
rLHS(7,2)+=gauss_weight*(crLHS239*crLHS486 + crLHS343 - crLHS344*crLHS70 + crLHS488*crLHS66);
rLHS(7,3)+=crLHS390;
rLHS(7,4)+=gauss_weight*(-crLHS122*crLHS87 + crLHS224*crLHS489 + crLHS315*crLHS490 + crLHS421);
rLHS(7,5)+=gauss_weight*(-crLHS106*crLHS266 + crLHS318*crLHS490 + crLHS446 + crLHS491*crLHS50);
rLHS(7,6)+=gauss_weight*(-crLHS118*crLHS344 + crLHS239*crLHS489 + crLHS466 + crLHS491*crLHS66);
rLHS(7,7)+=crLHS383*(crLHS411 + crLHS443 + crLHS465);
rLHS(7,8)+=gauss_weight*(-crLHS122*crLHS136 + crLHS224*crLHS493 + crLHS315*crLHS494 + crLHS492);
rLHS(7,9)+=gauss_weight*(-crLHS155*crLHS266 + crLHS318*crLHS494 + crLHS495 + crLHS496*crLHS50);
rLHS(7,10)+=gauss_weight*(-crLHS167*crLHS344 + crLHS239*crLHS493 + crLHS496*crLHS66 + crLHS497);
rLHS(7,11)+=crLHS498;
rLHS(7,12)+=gauss_weight*(-crLHS122*crLHS185 + crLHS224*crLHS500 + crLHS315*crLHS501 + crLHS499);
rLHS(7,13)+=gauss_weight*(-crLHS204*crLHS266 + crLHS318*crLHS501 + crLHS50*crLHS503 + crLHS502);
rLHS(7,14)+=gauss_weight*(-crLHS216*crLHS344 + crLHS239*crLHS500 + crLHS503*crLHS66 + crLHS504);
rLHS(7,15)+=crLHS505;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS131*crLHS27 + crLHS146 + crLHS27*crLHS507 - crLHS38*crLHS449 + crLHS508);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS42 + DN(2,1)*crLHS44 + DN(2,2)*crLHS47 + crLHS159 + crLHS270 - crLHS449*crLHS59 - crLHS50*crLHS510);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS60 + DN(2,1)*crLHS62 + DN(2,2)*crLHS64 + crLHS169 + crLHS347 - crLHS449*crLHS71 - crLHS510*crLHS66);
rLHS(8,3)+=-gauss_weight*(crLHS131*crLHS73 + crLHS391 + crLHS449*crLHS74 + crLHS507*crLHS73);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS75 + DN(2,1)*crLHS77 + DN(2,2)*crLHS79 + crLHS131*crLHS88 + crLHS425 - crLHS449*crLHS93 + crLHS507*crLHS88 + crLHS511);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS98 + DN(2,1)*crLHS100 + DN(2,2)*crLHS103 - crLHS107*crLHS449 + crLHS429 + crLHS447 - crLHS50*crLHS512);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS111 + DN(2,1)*crLHS113 + DN(2,2)*crLHS115 - crLHS119*crLHS449 + crLHS431 + crLHS468 - crLHS512*crLHS66);
rLHS(8,7)+=-gauss_weight*(crLHS122*crLHS131 + crLHS122*crLHS507 + crLHS123*crLHS449 + crLHS492);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS124 + DN(2,1)*crLHS126 + DN(2,2)*crLHS128 + crLHS131*crLHS137 + crLHS137*crLHS507 + crLHS14*crLHS513 - crLHS142*crLHS449 + crLHS15*crLHS515 + crLHS517);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS147 + DN(2,1)*crLHS149 + DN(2,2)*crLHS152 - crLHS156*crLHS449 + crLHS31*crLHS515 - crLHS50*crLHS520 - crLHS50*crLHS521 + crLHS519);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS160 + DN(2,1)*crLHS162 + DN(2,2)*crLHS164 - crLHS168*crLHS449 + crLHS35*crLHS515 - crLHS520*crLHS66 - crLHS521*crLHS66 + crLHS522);
rLHS(8,11)+=-gauss_weight*(crLHS131*crLHS171 + crLHS171*crLHS507 + crLHS172*crLHS449 + crLHS523);
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS173 + DN(2,1)*crLHS175 + DN(2,2)*crLHS177 + crLHS131*crLHS186 + crLHS186*crLHS507 - crLHS191*crLHS449 + crLHS525 + crLHS527);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS196 + DN(2,1)*crLHS198 + DN(2,2)*crLHS201 - crLHS205*crLHS449 - crLHS50*crLHS529 + crLHS528 + crLHS531);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS209 + DN(2,1)*crLHS211 + DN(2,2)*crLHS213 - crLHS217*crLHS449 - crLHS529*crLHS66 + crLHS532 + crLHS533);
rLHS(8,15)+=-gauss_weight*(crLHS131*crLHS220 + crLHS220*crLHS507 + crLHS221*crLHS449 + crLHS534);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS222 + DN(2,2)*crLHS223 + crLHS153 - crLHS224*crLHS510 - crLHS226*crLHS449 + crLHS273);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS44 + DN(2,1)*crLHS227 + DN(2,2)*crLHS229 + crLHS131*crLHS231 + crLHS231*crLHS507 - crLHS233*crLHS449 + crLHS281 + crLHS508);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS62 + DN(2,1)*crLHS234 + DN(2,2)*crLHS236 - crLHS239*crLHS510 - crLHS240*crLHS449 + crLHS287 + crLHS352);
rLHS(9,3)+=-gauss_weight*(crLHS131*crLHS242 + crLHS242*crLHS507 + crLHS243*crLHS449 + crLHS394);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS77 + DN(2,1)*crLHS244 + DN(2,2)*crLHS245 - crLHS224*crLHS512 - crLHS248*crLHS449 + crLHS426 + crLHS451);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS100 + DN(2,1)*crLHS251 + DN(2,2)*crLHS253 + crLHS131*crLHS254 + crLHS254*crLHS507 - crLHS256*crLHS449 + crLHS453 + crLHS511);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS113 + DN(2,1)*crLHS259 + DN(2,2)*crLHS261 - crLHS239*crLHS512 - crLHS263*crLHS449 + crLHS455 + crLHS472);
rLHS(9,7)+=-gauss_weight*(crLHS131*crLHS266 + crLHS266*crLHS507 + crLHS267*crLHS449 + crLHS495);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS126 + DN(2,1)*crLHS268 + DN(2,2)*crLHS269 - crLHS224*crLHS520 - crLHS224*crLHS521 - crLHS272*crLHS449 + crLHS30*crLHS515 + crLHS519);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS149 + DN(2,1)*crLHS274 + DN(2,2)*crLHS276 + crLHS131*crLHS277 + crLHS14*crLHS535 + crLHS277*crLHS507 - crLHS279*crLHS449 + crLHS515*crLHS56 + crLHS517);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS162 + DN(2,1)*crLHS282 + DN(2,2)*crLHS284 - crLHS239*crLHS520 - crLHS239*crLHS521 - crLHS286*crLHS449 + crLHS515*crLHS67 + crLHS537);
rLHS(9,11)+=-gauss_weight*(crLHS131*crLHS289 + crLHS289*crLHS507 + crLHS290*crLHS449 + crLHS538);
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS175 + DN(2,1)*crLHS291 + DN(2,2)*crLHS292 - crLHS224*crLHS529 - crLHS295*crLHS449 + crLHS539 + crLHS542);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS198 + DN(2,1)*crLHS297 + DN(2,2)*crLHS299 + crLHS131*crLHS300 + crLHS300*crLHS507 - crLHS302*crLHS449 + crLHS525 + crLHS544);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS211 + DN(2,1)*crLHS305 + DN(2,2)*crLHS307 - crLHS239*crLHS529 - crLHS309*crLHS449 + crLHS545 + crLHS546);
rLHS(9,15)+=-gauss_weight*(crLHS131*crLHS312 + crLHS312*crLHS507 + crLHS313*crLHS449 + crLHS547);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS223 + DN(2,2)*crLHS314 + crLHS165 - crLHS315*crLHS510 - crLHS316*crLHS449 + crLHS350);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS47 + DN(2,1)*crLHS229 + DN(2,2)*crLHS317 + crLHS285 - crLHS318*crLHS510 - crLHS319*crLHS449 + crLHS354);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS64 + DN(2,1)*crLHS236 + DN(2,2)*crLHS320 + crLHS131*crLHS322 + crLHS322*crLHS507 - crLHS323*crLHS449 + crLHS359 + crLHS508);
rLHS(10,3)+=-gauss_weight*(crLHS131*crLHS325 + crLHS325*crLHS507 + crLHS326*crLHS449 + crLHS396);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS79 + DN(2,1)*crLHS245 + DN(2,2)*crLHS327 - crLHS315*crLHS512 - crLHS331*crLHS449 + crLHS430 + crLHS471);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS103 + DN(2,1)*crLHS253 + DN(2,2)*crLHS334 - crLHS318*crLHS512 - crLHS336*crLHS449 + crLHS454 + crLHS473);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS115 + DN(2,1)*crLHS261 + DN(2,2)*crLHS338 + crLHS131*crLHS339 + crLHS339*crLHS507 - crLHS340*crLHS449 + crLHS475 + crLHS511);
rLHS(10,7)+=-gauss_weight*(crLHS131*crLHS344 + crLHS344*crLHS507 + crLHS345*crLHS449 + crLHS497);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS128 + DN(2,1)*crLHS269 + DN(2,2)*crLHS346 - crLHS315*crLHS520 - crLHS315*crLHS521 + crLHS34*crLHS515 - crLHS348*crLHS449 + crLHS522);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS152 + DN(2,1)*crLHS276 + DN(2,2)*crLHS351 - crLHS318*crLHS520 - crLHS318*crLHS521 - crLHS353*crLHS449 + crLHS515*crLHS55 + crLHS537);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS164 + DN(2,1)*crLHS284 + DN(2,2)*crLHS355 + crLHS131*crLHS356 + crLHS14*crLHS548 + crLHS356*crLHS507 - crLHS357*crLHS449 + crLHS515*crLHS68 + crLHS517);
rLHS(10,11)+=-gauss_weight*(crLHS131*crLHS361 + crLHS361*crLHS507 + crLHS362*crLHS449 + crLHS549);
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS177 + DN(2,1)*crLHS292 + DN(2,2)*crLHS363 - crLHS315*crLHS529 - crLHS365*crLHS449 + crLHS551 + crLHS552);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS201 + DN(2,1)*crLHS299 + DN(2,2)*crLHS368 - crLHS318*crLHS529 - crLHS370*crLHS449 + crLHS553 + crLHS554);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS213 + DN(2,1)*crLHS307 + DN(2,2)*crLHS372 + crLHS131*crLHS373 + crLHS373*crLHS507 - crLHS374*crLHS449 + crLHS525 + crLHS556);
rLHS(10,15)+=-gauss_weight*(crLHS131*crLHS378 + crLHS378*crLHS507 + crLHS379*crLHS449 + crLHS557);
rLHS(11,0)+=gauss_weight*(crLHS170 - crLHS171*crLHS25 + crLHS224*crLHS558 + crLHS315*crLHS559);
rLHS(11,1)+=gauss_weight*(crLHS288 - crLHS289*crLHS58 + crLHS318*crLHS559 + crLHS50*crLHS560);
rLHS(11,2)+=gauss_weight*(crLHS239*crLHS558 + crLHS360 - crLHS361*crLHS70 + crLHS560*crLHS66);
rLHS(11,3)+=crLHS397;
rLHS(11,4)+=gauss_weight*(-crLHS171*crLHS87 + crLHS224*crLHS561 + crLHS315*crLHS562 + crLHS432);
rLHS(11,5)+=gauss_weight*(-crLHS106*crLHS289 + crLHS318*crLHS562 + crLHS456 + crLHS50*crLHS563);
rLHS(11,6)+=gauss_weight*(-crLHS118*crLHS361 + crLHS239*crLHS561 + crLHS476 + crLHS563*crLHS66);
rLHS(11,7)+=crLHS498;
rLHS(11,8)+=gauss_weight*(-crLHS136*crLHS171 + crLHS224*crLHS564 + crLHS315*crLHS565 + crLHS523);
rLHS(11,9)+=gauss_weight*(-crLHS155*crLHS289 + crLHS318*crLHS565 + crLHS50*crLHS566 + crLHS538);
rLHS(11,10)+=gauss_weight*(-crLHS167*crLHS361 + crLHS239*crLHS564 + crLHS549 + crLHS566*crLHS66);
rLHS(11,11)+=crLHS383*(crLHS513 + crLHS535 + crLHS548);
rLHS(11,12)+=gauss_weight*(-crLHS171*crLHS185 + crLHS224*crLHS568 + crLHS315*crLHS569 + crLHS567);
rLHS(11,13)+=gauss_weight*(-crLHS204*crLHS289 + crLHS318*crLHS569 + crLHS50*crLHS571 + crLHS570);
rLHS(11,14)+=gauss_weight*(-crLHS216*crLHS361 + crLHS239*crLHS568 + crLHS571*crLHS66 + crLHS572);
rLHS(11,15)+=crLHS573;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS180*crLHS27 + crLHS195 + crLHS27*crLHS575 - crLHS38*crLHS458 + crLHS576);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS42 + DN(3,1)*crLHS44 + DN(3,2)*crLHS47 + crLHS208 + crLHS293 - crLHS458*crLHS59 - crLHS50*crLHS578);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS60 + DN(3,1)*crLHS62 + DN(3,2)*crLHS64 + crLHS218 + crLHS364 - crLHS458*crLHS71 - crLHS578*crLHS66);
rLHS(12,3)+=-gauss_weight*(crLHS180*crLHS73 + crLHS398 + crLHS458*crLHS74 + crLHS575*crLHS73);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS75 + DN(3,1)*crLHS77 + DN(3,2)*crLHS79 + crLHS180*crLHS88 + crLHS436 - crLHS458*crLHS93 + crLHS575*crLHS88 + crLHS579);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS98 + DN(3,1)*crLHS100 + DN(3,2)*crLHS103 - crLHS107*crLHS458 + crLHS439 + crLHS457 - crLHS50*crLHS580);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS111 + DN(3,1)*crLHS113 + DN(3,2)*crLHS115 - crLHS119*crLHS458 + crLHS441 + crLHS477 - crLHS580*crLHS66);
rLHS(12,7)+=-gauss_weight*(crLHS122*crLHS180 + crLHS122*crLHS575 + crLHS123*crLHS458 + crLHS499);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS124 + DN(3,1)*crLHS126 + DN(3,2)*crLHS128 + crLHS137*crLHS180 + crLHS137*crLHS575 - crLHS142*crLHS458 + crLHS527 + crLHS581);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS147 + DN(3,1)*crLHS149 + DN(3,2)*crLHS152 - crLHS156*crLHS458 - crLHS50*crLHS582 + crLHS531 + crLHS539);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS160 + DN(3,1)*crLHS162 + DN(3,2)*crLHS164 - crLHS168*crLHS458 + crLHS533 + crLHS551 - crLHS582*crLHS66);
rLHS(12,11)+=-gauss_weight*(crLHS171*crLHS180 + crLHS171*crLHS575 + crLHS172*crLHS458 + crLHS567);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS173 + DN(3,1)*crLHS175 + DN(3,2)*crLHS177 + crLHS14*crLHS583 + crLHS15*crLHS585 + crLHS180*crLHS186 + crLHS186*crLHS575 - crLHS191*crLHS458 + crLHS587);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS196 + DN(3,1)*crLHS198 + DN(3,2)*crLHS201 - crLHS205*crLHS458 + crLHS31*crLHS585 - crLHS50*crLHS590 - crLHS50*crLHS591 + crLHS589);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS209 + DN(3,1)*crLHS211 + DN(3,2)*crLHS213 - crLHS217*crLHS458 + crLHS35*crLHS585 - crLHS590*crLHS66 - crLHS591*crLHS66 + crLHS592);
rLHS(12,15)+=-gauss_weight*(crLHS180*crLHS220 + crLHS220*crLHS575 + crLHS221*crLHS458 + crLHS593);
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS222 + DN(3,2)*crLHS223 + crLHS202 - crLHS224*crLHS578 - crLHS226*crLHS458 + crLHS296);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS44 + DN(3,1)*crLHS227 + DN(3,2)*crLHS229 + crLHS180*crLHS231 + crLHS231*crLHS575 - crLHS233*crLHS458 + crLHS304 + crLHS576);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS62 + DN(3,1)*crLHS234 + DN(3,2)*crLHS236 - crLHS239*crLHS578 - crLHS240*crLHS458 + crLHS310 + crLHS369);
rLHS(13,3)+=-gauss_weight*(crLHS180*crLHS242 + crLHS242*crLHS575 + crLHS243*crLHS458 + crLHS401);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS77 + DN(3,1)*crLHS244 + DN(3,2)*crLHS245 - crLHS224*crLHS580 - crLHS248*crLHS458 + crLHS437 + crLHS459);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS100 + DN(3,1)*crLHS251 + DN(3,2)*crLHS253 + crLHS180*crLHS254 + crLHS254*crLHS575 - crLHS256*crLHS458 + crLHS461 + crLHS579);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS113 + DN(3,1)*crLHS259 + DN(3,2)*crLHS261 - crLHS239*crLHS580 - crLHS263*crLHS458 + crLHS463 + crLHS481);
rLHS(13,7)+=-gauss_weight*(crLHS180*crLHS266 + crLHS266*crLHS575 + crLHS267*crLHS458 + crLHS502);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS126 + DN(3,1)*crLHS268 + DN(3,2)*crLHS269 - crLHS224*crLHS582 - crLHS272*crLHS458 + crLHS528 + crLHS542);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS149 + DN(3,1)*crLHS274 + DN(3,2)*crLHS276 + crLHS180*crLHS277 + crLHS277*crLHS575 - crLHS279*crLHS458 + crLHS544 + crLHS581);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS162 + DN(3,1)*crLHS282 + DN(3,2)*crLHS284 - crLHS239*crLHS582 - crLHS286*crLHS458 + crLHS546 + crLHS553);
rLHS(13,11)+=-gauss_weight*(crLHS180*crLHS289 + crLHS289*crLHS575 + crLHS290*crLHS458 + crLHS570);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS175 + DN(3,1)*crLHS291 + DN(3,2)*crLHS292 - crLHS224*crLHS590 - crLHS224*crLHS591 - crLHS295*crLHS458 + crLHS30*crLHS585 + crLHS589);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS198 + DN(3,1)*crLHS297 + DN(3,2)*crLHS299 + crLHS14*crLHS594 + crLHS180*crLHS300 + crLHS300*crLHS575 - crLHS302*crLHS458 + crLHS56*crLHS585 + crLHS587);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS211 + DN(3,1)*crLHS305 + DN(3,2)*crLHS307 - crLHS239*crLHS590 - crLHS239*crLHS591 - crLHS309*crLHS458 + crLHS585*crLHS67 + crLHS595);
rLHS(13,15)+=-gauss_weight*(crLHS180*crLHS312 + crLHS312*crLHS575 + crLHS313*crLHS458 + crLHS596);
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS223 + DN(3,2)*crLHS314 + crLHS214 - crLHS315*crLHS578 - crLHS316*crLHS458 + crLHS367);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS47 + DN(3,1)*crLHS229 + DN(3,2)*crLHS317 + crLHS308 - crLHS318*crLHS578 - crLHS319*crLHS458 + crLHS371);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS64 + DN(3,1)*crLHS236 + DN(3,2)*crLHS320 + crLHS180*crLHS322 + crLHS322*crLHS575 - crLHS323*crLHS458 + crLHS376 + crLHS576);
rLHS(14,3)+=-gauss_weight*(crLHS180*crLHS325 + crLHS325*crLHS575 + crLHS326*crLHS458 + crLHS403);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS79 + DN(3,1)*crLHS245 + DN(3,2)*crLHS327 - crLHS315*crLHS580 - crLHS331*crLHS458 + crLHS440 + crLHS480);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS103 + DN(3,1)*crLHS253 + DN(3,2)*crLHS334 - crLHS318*crLHS580 - crLHS336*crLHS458 + crLHS462 + crLHS482);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS115 + DN(3,1)*crLHS261 + DN(3,2)*crLHS338 + crLHS180*crLHS339 + crLHS339*crLHS575 - crLHS340*crLHS458 + crLHS484 + crLHS579);
rLHS(14,7)+=-gauss_weight*(crLHS180*crLHS344 + crLHS344*crLHS575 + crLHS345*crLHS458 + crLHS504);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS128 + DN(3,1)*crLHS269 + DN(3,2)*crLHS346 - crLHS315*crLHS582 - crLHS348*crLHS458 + crLHS532 + crLHS552);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS152 + DN(3,1)*crLHS276 + DN(3,2)*crLHS351 - crLHS318*crLHS582 - crLHS353*crLHS458 + crLHS545 + crLHS554);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS164 + DN(3,1)*crLHS284 + DN(3,2)*crLHS355 + crLHS180*crLHS356 + crLHS356*crLHS575 - crLHS357*crLHS458 + crLHS556 + crLHS581);
rLHS(14,11)+=-gauss_weight*(crLHS180*crLHS361 + crLHS361*crLHS575 + crLHS362*crLHS458 + crLHS572);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS177 + DN(3,1)*crLHS292 + DN(3,2)*crLHS363 - crLHS315*crLHS590 - crLHS315*crLHS591 + crLHS34*crLHS585 - crLHS365*crLHS458 + crLHS592);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS201 + DN(3,1)*crLHS299 + DN(3,2)*crLHS368 - crLHS318*crLHS590 - crLHS318*crLHS591 - crLHS370*crLHS458 + crLHS55*crLHS585 + crLHS595);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS213 + DN(3,1)*crLHS307 + DN(3,2)*crLHS372 + crLHS14*crLHS597 + crLHS180*crLHS373 + crLHS373*crLHS575 - crLHS374*crLHS458 + crLHS585*crLHS68 + crLHS587);
rLHS(14,15)+=-gauss_weight*(crLHS180*crLHS378 + crLHS378*crLHS575 + crLHS379*crLHS458 + crLHS598);
rLHS(15,0)+=gauss_weight*(crLHS219 - crLHS220*crLHS25 + crLHS224*crLHS599 + crLHS315*crLHS600);
rLHS(15,1)+=gauss_weight*(crLHS311 - crLHS312*crLHS58 + crLHS318*crLHS600 + crLHS50*crLHS601);
rLHS(15,2)+=gauss_weight*(crLHS239*crLHS599 + crLHS377 - crLHS378*crLHS70 + crLHS601*crLHS66);
rLHS(15,3)+=crLHS404;
rLHS(15,4)+=gauss_weight*(-crLHS220*crLHS87 + crLHS224*crLHS602 + crLHS315*crLHS603 + crLHS442);
rLHS(15,5)+=gauss_weight*(-crLHS106*crLHS312 + crLHS318*crLHS603 + crLHS464 + crLHS50*crLHS604);
rLHS(15,6)+=gauss_weight*(-crLHS118*crLHS378 + crLHS239*crLHS602 + crLHS485 + crLHS604*crLHS66);
rLHS(15,7)+=crLHS505;
rLHS(15,8)+=gauss_weight*(-crLHS136*crLHS220 + crLHS224*crLHS605 + crLHS315*crLHS606 + crLHS534);
rLHS(15,9)+=gauss_weight*(-crLHS155*crLHS312 + crLHS318*crLHS606 + crLHS50*crLHS607 + crLHS547);
rLHS(15,10)+=gauss_weight*(-crLHS167*crLHS378 + crLHS239*crLHS605 + crLHS557 + crLHS607*crLHS66);
rLHS(15,11)+=crLHS573;
rLHS(15,12)+=gauss_weight*(-crLHS185*crLHS220 + crLHS224*crLHS608 + crLHS315*crLHS609 + crLHS593);
rLHS(15,13)+=gauss_weight*(-crLHS204*crLHS312 + crLHS318*crLHS609 + crLHS50*crLHS610 + crLHS596);
rLHS(15,14)+=gauss_weight*(-crLHS216*crLHS378 + crLHS239*crLHS608 + crLHS598 + crLHS610*crLHS66);
rLHS(15,15)+=crLHS383*(crLHS583 + crLHS594 + crLHS597);

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
const double crRHS6 = N[0]*crRHS2;
const double crRHS7 = 2.0*functional_weights[0];
const double crRHS8 = crRHS6*crRHS7;
const double crRHS9 = DN(0,1)*v_ns(0,0);
const double crRHS10 = DN(1,1)*v_ns(1,0);
const double crRHS11 = DN(2,1)*v_ns(2,0);
const double crRHS12 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1);
const double crRHS13 = 2.0*functional_weights[2]*mu*(-crRHS10 - crRHS11 + crRHS12 - crRHS9);
const double crRHS14 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1);
const double crRHS15 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS14);
const double crRHS16 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0);
const double crRHS17 = 1.0*crRHS16;
const double crRHS18 = crRHS10 + crRHS11 + crRHS9;
const double crRHS19 = 0.5*crRHS12 + 0.5*crRHS18;
const double crRHS20 = 4.0*functional_weights[1]*mu;
const double crRHS21 = crRHS16*crRHS3;
const double crRHS22 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1);
const double crRHS23 = crRHS18*crRHS22;
const double crRHS24 = crRHS21 + crRHS23;
const double crRHS25 = N[0]*rho;
const double crRHS26 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0);
const double crRHS27 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1);
const double crRHS28 = crRHS26 + crRHS27;
const double crRHS29 = crRHS2*stab_c3;
const double crRHS30 = pow(crRHS14, 2) + pow(crRHS5, 2);
const double crRHS31 = pow(crRHS30, 1.0/4.0)*stab_c3;
const double crRHS32 = sqrt(crRHS30)*rho*stab_c2;
const double crRHS33 = crRHS28*(h*(crRHS29*h + crRHS31*h + crRHS32)/stab_c1 + mu);
const double crRHS34 = crRHS5*rho;
const double crRHS35 = crRHS14*rho;
const double crRHS36 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS21*rho + crRHS23*rho - crRHS26*crRHS34 - crRHS35*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0)) + crRHS4;
const double crRHS37 = 1.0/(crRHS29 + crRHS31 + crRHS32/h + mu*stab_c1/pow(h, 2));
const double crRHS38 = crRHS36*crRHS37;
const double crRHS39 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS40 = crRHS2*crRHS22;
const double crRHS41 = crRHS12*crRHS3;
const double crRHS42 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crRHS43 = crRHS22*crRHS42;
const double crRHS44 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] - crRHS27*crRHS35 - crRHS34*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS39 + crRHS40 + crRHS41*rho + crRHS43*rho;
const double crRHS45 = crRHS16*crRHS36 + crRHS18*crRHS44;
const double crRHS46 = crRHS25*crRHS37;
const double crRHS47 = 1.0*crRHS42;
const double crRHS48 = crRHS41 + crRHS43;
const double crRHS49 = crRHS37*crRHS44;
const double crRHS50 = crRHS12*crRHS36 + crRHS42*crRHS44;
const double crRHS51 = N[1]*crRHS2;
const double crRHS52 = crRHS51*crRHS7;
const double crRHS53 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS14);
const double crRHS54 = N[1]*rho;
const double crRHS55 = crRHS37*crRHS54;
const double crRHS56 = N[2]*crRHS2;
const double crRHS57 = crRHS56*crRHS7;
const double crRHS58 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS14);
const double crRHS59 = N[2]*rho;
const double crRHS60 = crRHS37*crRHS59;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS33 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS13 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS4 + crRHS15*crRHS3 - crRHS15*crRHS38 + crRHS20*(DN(0,0)*crRHS17 + DN(0,1)*crRHS19) + crRHS24*crRHS25 - crRHS38*crRHS6 - crRHS45*crRHS46 + crRHS5*crRHS8);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS13 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS33 + DN(0,1)*stress_adj[1] - N[0]*crRHS39 + N[0]*crRHS40 + crRHS14*crRHS8 + crRHS15*crRHS22 - crRHS15*crRHS49 + crRHS20*(DN(0,0)*crRHS19 + DN(0,1)*crRHS47) + crRHS25*crRHS48 - crRHS46*crRHS50 - crRHS49*crRHS6);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS38 + DN(0,1)*crRHS49 + N[0]*crRHS28);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS33 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS13 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS4 + crRHS20*(DN(1,0)*crRHS17 + DN(1,1)*crRHS19) + crRHS24*crRHS54 + crRHS3*crRHS53 - crRHS38*crRHS51 - crRHS38*crRHS53 - crRHS45*crRHS55 + crRHS5*crRHS52);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS13 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS33 + DN(1,1)*stress_adj[1] - N[1]*crRHS39 + N[1]*crRHS40 + crRHS14*crRHS52 + crRHS20*(DN(1,0)*crRHS19 + DN(1,1)*crRHS47) + crRHS22*crRHS53 + crRHS48*crRHS54 - crRHS49*crRHS51 - crRHS49*crRHS53 - crRHS50*crRHS55);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS38 + DN(1,1)*crRHS49 + N[1]*crRHS28);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS33 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS13 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS4 + crRHS20*(DN(2,0)*crRHS17 + DN(2,1)*crRHS19) + crRHS24*crRHS59 + crRHS3*crRHS58 - crRHS38*crRHS56 - crRHS38*crRHS58 - crRHS45*crRHS60 + crRHS5*crRHS57);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS13 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS33 + DN(2,1)*stress_adj[1] - N[2]*crRHS39 + N[2]*crRHS40 + crRHS14*crRHS57 + crRHS20*(DN(2,0)*crRHS19 + DN(2,1)*crRHS47) + crRHS22*crRHS58 + crRHS48*crRHS59 - crRHS49*crRHS56 - crRHS49*crRHS58 - crRHS50*crRHS60);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS38 + DN(2,1)*crRHS49 + N[2]*crRHS28);

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
const double crRHS6 = N[0]*crRHS2;
const double crRHS7 = 2.0*functional_weights[0];
const double crRHS8 = crRHS6*crRHS7;
const double crRHS9 = N[0]*v_ns(0,1) + N[1]*v_ns(1,1) + N[2]*v_ns(2,1) + N[3]*v_ns(3,1);
const double crRHS10 = N[0]*v_ns(0,2) + N[1]*v_ns(1,2) + N[2]*v_ns(2,2) + N[3]*v_ns(3,2);
const double crRHS11 = rho*(DN(0,0)*crRHS5 + DN(0,1)*crRHS9 + DN(0,2)*crRHS10);
const double crRHS12 = DN(0,1)*v_ns(0,0);
const double crRHS13 = DN(1,1)*v_ns(1,0);
const double crRHS14 = DN(2,1)*v_ns(2,0);
const double crRHS15 = DN(3,1)*v_ns(3,0);
const double crRHS16 = DN(0,0)*v_ns(0,1) + DN(1,0)*v_ns(1,1) + DN(2,0)*v_ns(2,1) + DN(3,0)*v_ns(3,1);
const double crRHS17 = -crRHS12 - crRHS13 - crRHS14 - crRHS15 + crRHS16;
const double crRHS18 = DN(0,2)*v_ns(0,0);
const double crRHS19 = DN(1,2)*v_ns(1,0);
const double crRHS20 = DN(2,2)*v_ns(2,0);
const double crRHS21 = DN(3,2)*v_ns(3,0);
const double crRHS22 = DN(0,0)*v_ns(0,2) + DN(1,0)*v_ns(1,2) + DN(2,0)*v_ns(2,2) + DN(3,0)*v_ns(3,2);
const double crRHS23 = -crRHS18 - crRHS19 - crRHS20 - crRHS21 + crRHS22;
const double crRHS24 = 2.0*functional_weights[2]*mu;
const double crRHS25 = DN(0,0)*v_ns(0,0) + DN(1,0)*v_ns(1,0) + DN(2,0)*v_ns(2,0) + DN(3,0)*v_ns(3,0);
const double crRHS26 = 1.0*crRHS25;
const double crRHS27 = crRHS12 + crRHS13 + crRHS14 + crRHS15;
const double crRHS28 = 0.5*crRHS16 + 0.5*crRHS27;
const double crRHS29 = crRHS18 + crRHS19 + crRHS20 + crRHS21;
const double crRHS30 = crRHS22 + crRHS29;
const double crRHS31 = 0.5*DN(0,2);
const double crRHS32 = 4.0*functional_weights[1]*mu;
const double crRHS33 = crRHS25*crRHS3;
const double crRHS34 = N[0]*v_adj(0,1) + N[1]*v_adj(1,1) + N[2]*v_adj(2,1) + N[3]*v_adj(3,1);
const double crRHS35 = crRHS27*crRHS34;
const double crRHS36 = N[0]*v_adj(0,2) + N[1]*v_adj(1,2) + N[2]*v_adj(2,2) + N[3]*v_adj(3,2);
const double crRHS37 = crRHS29*crRHS36;
const double crRHS38 = crRHS33 + crRHS35 + crRHS37;
const double crRHS39 = N[0]*rho;
const double crRHS40 = DN(0,0)*v_adj(0,0) + DN(1,0)*v_adj(1,0) + DN(2,0)*v_adj(2,0) + DN(3,0)*v_adj(3,0);
const double crRHS41 = DN(0,1)*v_adj(0,1) + DN(1,1)*v_adj(1,1) + DN(2,1)*v_adj(2,1) + DN(3,1)*v_adj(3,1);
const double crRHS42 = DN(0,2)*v_adj(0,2) + DN(1,2)*v_adj(1,2) + DN(2,2)*v_adj(2,2) + DN(3,2)*v_adj(3,2);
const double crRHS43 = crRHS40 + crRHS41 + crRHS42;
const double crRHS44 = crRHS2*stab_c3;
const double crRHS45 = pow(crRHS10, 2) + pow(crRHS5, 2) + pow(crRHS9, 2);
const double crRHS46 = pow(crRHS45, 1.0/4.0)*stab_c3;
const double crRHS47 = sqrt(crRHS45)*rho*stab_c2;
const double crRHS48 = crRHS43*(h*(crRHS44*h + crRHS46*h + crRHS47)/stab_c1 + mu);
const double crRHS49 = crRHS5*rho;
const double crRHS50 = crRHS9*rho;
const double crRHS51 = crRHS10*rho;
const double crRHS52 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 + crRHS33*rho + crRHS35*rho + crRHS37*rho + crRHS4 - crRHS40*crRHS49 - crRHS50*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS51*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0));
const double crRHS53 = 1.0/(crRHS44 + crRHS46 + crRHS47/h + mu*stab_c1/pow(h, 2));
const double crRHS54 = crRHS52*crRHS53;
const double crRHS55 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS56 = crRHS2*crRHS34;
const double crRHS57 = crRHS16*crRHS3;
const double crRHS58 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crRHS59 = crRHS34*crRHS58;
const double crRHS60 = DN(0,2)*v_ns(0,1);
const double crRHS61 = DN(1,2)*v_ns(1,1);
const double crRHS62 = DN(2,2)*v_ns(2,1);
const double crRHS63 = DN(3,2)*v_ns(3,1);
const double crRHS64 = crRHS60 + crRHS61 + crRHS62 + crRHS63;
const double crRHS65 = crRHS36*crRHS64;
const double crRHS66 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] - crRHS41*crRHS50 - crRHS49*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS51*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS55 + crRHS56 + crRHS57*rho + crRHS59*rho + crRHS65*rho;
const double crRHS67 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS68 = crRHS2*crRHS36;
const double crRHS69 = crRHS22*crRHS3;
const double crRHS70 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crRHS71 = crRHS34*crRHS70;
const double crRHS72 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crRHS73 = crRHS36*crRHS72;
const double crRHS74 = DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] - crRHS42*crRHS51 - crRHS49*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS50*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS67 + crRHS68 + crRHS69*rho + crRHS71*rho + crRHS73*rho;
const double crRHS75 = crRHS25*crRHS52 + crRHS27*crRHS66 + crRHS29*crRHS74;
const double crRHS76 = crRHS39*crRHS53;
const double crRHS77 = -crRHS60 - crRHS61 - crRHS62 - crRHS63 + crRHS70;
const double crRHS78 = 1.0*crRHS58;
const double crRHS79 = crRHS64 + crRHS70;
const double crRHS80 = crRHS57 + crRHS59 + crRHS65;
const double crRHS81 = crRHS53*crRHS66;
const double crRHS82 = crRHS16*crRHS52 + crRHS58*crRHS66 + crRHS64*crRHS74;
const double crRHS83 = 1.0*crRHS72;
const double crRHS84 = 0.5*crRHS30;
const double crRHS85 = 0.5*crRHS79;
const double crRHS86 = crRHS69 + crRHS71 + crRHS73;
const double crRHS87 = crRHS53*crRHS74;
const double crRHS88 = crRHS22*crRHS52 + crRHS66*crRHS70 + crRHS72*crRHS74;
const double crRHS89 = N[1]*crRHS2;
const double crRHS90 = crRHS7*crRHS89;
const double crRHS91 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS9 + DN(1,2)*crRHS10);
const double crRHS92 = N[1]*rho;
const double crRHS93 = crRHS53*crRHS92;
const double crRHS94 = N[2]*crRHS2;
const double crRHS95 = crRHS7*crRHS94;
const double crRHS96 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS9 + DN(2,2)*crRHS10);
const double crRHS97 = N[2]*rho;
const double crRHS98 = crRHS53*crRHS97;
const double crRHS99 = N[3]*crRHS2;
const double crRHS100 = crRHS7*crRHS99;
const double crRHS101 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS9 + DN(3,2)*crRHS10);
const double crRHS102 = N[3]*rho;
const double crRHS103 = crRHS102*crRHS53;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS48 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS4 + crRHS11*crRHS3 - crRHS11*crRHS54 - crRHS24*(DN(0,1)*crRHS17 + DN(0,2)*crRHS23) + crRHS32*(DN(0,0)*crRHS26 + DN(0,1)*crRHS28 + crRHS30*crRHS31) + crRHS38*crRHS39 + crRHS5*crRHS8 - crRHS54*crRHS6 - crRHS75*crRHS76);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS48 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS55 + N[0]*crRHS56 + crRHS11*crRHS34 - crRHS11*crRHS81 + crRHS24*(DN(0,0)*crRHS17 - DN(0,2)*crRHS77) + crRHS32*(DN(0,0)*crRHS28 + DN(0,1)*crRHS78 + crRHS31*crRHS79) + crRHS39*crRHS80 - crRHS6*crRHS81 - crRHS76*crRHS82 + crRHS8*crRHS9);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS48 + DN(0,2)*stress_adj[2] - N[0]*crRHS67 + N[0]*crRHS68 + crRHS10*crRHS8 + crRHS11*crRHS36 - crRHS11*crRHS87 + crRHS24*(DN(0,0)*crRHS23 + DN(0,1)*crRHS77) + crRHS32*(DN(0,0)*crRHS84 + DN(0,1)*crRHS85 + DN(0,2)*crRHS83) + crRHS39*crRHS86 - crRHS6*crRHS87 - crRHS76*crRHS88);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS54 + DN(0,1)*crRHS81 + DN(0,2)*crRHS87 + N[0]*crRHS43);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS48 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS4 - crRHS24*(DN(1,1)*crRHS17 + DN(1,2)*crRHS23) + crRHS3*crRHS91 + crRHS32*(DN(1,0)*crRHS26 + DN(1,1)*crRHS28 + DN(1,2)*crRHS84) + crRHS38*crRHS92 + crRHS5*crRHS90 - crRHS54*crRHS89 - crRHS54*crRHS91 - crRHS75*crRHS93);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS48 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS55 + N[1]*crRHS56 + crRHS24*(DN(1,0)*crRHS17 - DN(1,2)*crRHS77) + crRHS32*(DN(1,0)*crRHS28 + DN(1,1)*crRHS78 + DN(1,2)*crRHS85) + crRHS34*crRHS91 + crRHS80*crRHS92 - crRHS81*crRHS89 - crRHS81*crRHS91 - crRHS82*crRHS93 + crRHS9*crRHS90);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS48 + DN(1,2)*stress_adj[2] - N[1]*crRHS67 + N[1]*crRHS68 + crRHS10*crRHS90 + crRHS24*(DN(1,0)*crRHS23 + DN(1,1)*crRHS77) + crRHS32*(DN(1,0)*crRHS84 + DN(1,1)*crRHS85 + DN(1,2)*crRHS83) + crRHS36*crRHS91 + crRHS86*crRHS92 - crRHS87*crRHS89 - crRHS87*crRHS91 - crRHS88*crRHS93);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS54 + DN(1,1)*crRHS81 + DN(1,2)*crRHS87 + N[1]*crRHS43);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS48 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS4 - crRHS24*(DN(2,1)*crRHS17 + DN(2,2)*crRHS23) + crRHS3*crRHS96 + crRHS32*(DN(2,0)*crRHS26 + DN(2,1)*crRHS28 + DN(2,2)*crRHS84) + crRHS38*crRHS97 + crRHS5*crRHS95 - crRHS54*crRHS94 - crRHS54*crRHS96 - crRHS75*crRHS98);
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS48 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS55 + N[2]*crRHS56 + crRHS24*(DN(2,0)*crRHS17 - DN(2,2)*crRHS77) + crRHS32*(DN(2,0)*crRHS28 + DN(2,1)*crRHS78 + DN(2,2)*crRHS85) + crRHS34*crRHS96 + crRHS80*crRHS97 - crRHS81*crRHS94 - crRHS81*crRHS96 - crRHS82*crRHS98 + crRHS9*crRHS95);
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS48 + DN(2,2)*stress_adj[2] - N[2]*crRHS67 + N[2]*crRHS68 + crRHS10*crRHS95 + crRHS24*(DN(2,0)*crRHS23 + DN(2,1)*crRHS77) + crRHS32*(DN(2,0)*crRHS84 + DN(2,1)*crRHS85 + DN(2,2)*crRHS83) + crRHS36*crRHS96 + crRHS86*crRHS97 - crRHS87*crRHS94 - crRHS87*crRHS96 - crRHS88*crRHS98);
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS54 + DN(2,1)*crRHS81 + DN(2,2)*crRHS87 + N[2]*crRHS43);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS48 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS4 + crRHS100*crRHS5 + crRHS101*crRHS3 - crRHS101*crRHS54 + crRHS102*crRHS38 - crRHS103*crRHS75 - crRHS24*(DN(3,1)*crRHS17 + DN(3,2)*crRHS23) + crRHS32*(DN(3,0)*crRHS26 + DN(3,1)*crRHS28 + DN(3,2)*crRHS84) - crRHS54*crRHS99);
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS48 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS55 + N[3]*crRHS56 + crRHS100*crRHS9 + crRHS101*crRHS34 - crRHS101*crRHS81 + crRHS102*crRHS80 - crRHS103*crRHS82 + crRHS24*(DN(3,0)*crRHS17 - DN(3,2)*crRHS77) + crRHS32*(DN(3,0)*crRHS28 + DN(3,1)*crRHS78 + DN(3,2)*crRHS85) - crRHS81*crRHS99);
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS48 + DN(3,2)*stress_adj[2] - N[3]*crRHS67 + N[3]*crRHS68 + crRHS10*crRHS100 + crRHS101*crRHS36 - crRHS101*crRHS87 + crRHS102*crRHS86 - crRHS103*crRHS88 + crRHS24*(DN(3,0)*crRHS23 + DN(3,1)*crRHS77) + crRHS32*(DN(3,0)*crRHS84 + DN(3,1)*crRHS85 + DN(3,2)*crRHS83) - crRHS87*crRHS99);
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS54 + DN(3,1)*crRHS81 + DN(3,2)*crRHS87 + N[3]*crRHS43);

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