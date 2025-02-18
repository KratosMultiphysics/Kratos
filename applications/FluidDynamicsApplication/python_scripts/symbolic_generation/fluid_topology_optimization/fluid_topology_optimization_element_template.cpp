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
    //substitute_lhs_2D3N_NS    
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
    //substitute_lhs_3D4N_NS
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
    //substitute_rhs_2D3N_NS
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
    //substitute_rhs_3D4N_NS
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
    //substitute_lhs_2D3N_ADJ_NS    
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
    //substitute_lhs_3D4N_ADJ_NS
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
    // const array_1d<double,3>& pn_adj = rData.Pressure_adj_OldStep1;
    // const array_1d<double,3>& pnn_adj = rData.Pressure_adj_OldStep2;
    const array_1d<double,3>& stress_adj = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL RHS VECTOR
    //substitute_rhs_2D3N_ADJ_NS
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
    const array_1d<double,6>& stress_adj = rData.ShearStress;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;

    // ADJOINT NAVIER-STOKES ELEMENTAL RHS VECTOR
    //substitute_rhs_3D4N_ADJ_NS
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