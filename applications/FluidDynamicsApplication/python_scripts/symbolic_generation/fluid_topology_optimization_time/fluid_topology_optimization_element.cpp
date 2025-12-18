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
const double crLHS42 = crLHS24*crLHS35 + crLHS41;
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
const double crLHS77 = crLHS35*crLHS64 + crLHS76;
const double crLHS78 = DN(0,0)*DN(1,0);
const double crLHS79 = N[1]*crLHS25 + crLHS20*crLHS78;
const double crLHS80 = C(0,1)*DN(1,1) + crLHS62;
const double crLHS81 = C(1,2)*DN(1,1);
const double crLHS82 = C(2,2)*DN(1,0) + crLHS81;
const double crLHS83 = DN(1,1)*crLHS46;
const double crLHS84 = crLHS18*crLHS64;
const double crLHS85 = crLHS65 + crLHS66 - crLHS68 - crLHS70 + crLHS84;
const double crLHS86 = crLHS13*crLHS24;
const double crLHS87 = crLHS48*rho;
const double crLHS88 = N[1]*crLHS86 - crLHS76*crLHS87;
const double crLHS89 = DN(0,0)*N[1];
const double crLHS90 = DN(1,0)*crLHS33;
const double crLHS91 = DN(1,0)*crLHS9 + DN(1,1)*crLHS13;
const double crLHS92 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double crLHS93 = C(0,2)*DN(2,0);
const double crLHS94 = C(2,2)*DN(2,1) + crLHS93;
const double crLHS95 = N[2]*rho;
const double crLHS96 = crLHS9*crLHS95;
const double crLHS97 = N[2]*crLHS4;
const double crLHS98 = DN(2,0)*crLHS6;
const double crLHS99 = crLHS98*rho;
const double crLHS100 = DN(2,1)*crLHS7;
const double crLHS101 = crLHS100*rho;
const double crLHS102 = crLHS101 - crLHS97 + crLHS99;
const double crLHS103 = crLHS102 - crLHS96;
const double crLHS104 = crLHS103*crLHS33;
const double crLHS105 = crLHS37*crLHS95;
const double crLHS106 = -crLHS103*crLHS9 + crLHS105;
const double crLHS107 = N[2]*crLHS23;
const double crLHS108 = crLHS107 + crLHS35*crLHS95;
const double crLHS109 = DN(0,0)*DN(2,0);
const double crLHS110 = N[2]*crLHS25 + crLHS109*crLHS20;
const double crLHS111 = C(0,1)*DN(2,1) + crLHS93;
const double crLHS112 = C(1,2)*DN(2,1);
const double crLHS113 = C(2,2)*DN(2,0) + crLHS112;
const double crLHS114 = DN(2,1)*crLHS46;
const double crLHS115 = crLHS18*crLHS95;
const double crLHS116 = -crLHS101 + crLHS115 + crLHS96 + crLHS97 - crLHS99;
const double crLHS117 = N[2]*crLHS86 - crLHS107*crLHS87;
const double crLHS118 = DN(0,0)*N[2];
const double crLHS119 = DN(2,0)*crLHS33;
const double crLHS120 = DN(2,0)*crLHS9 + DN(2,1)*crLHS13;
const double crLHS121 = C(0,1)*DN(0,0) + crLHS44;
const double crLHS122 = crLHS17*crLHS33;
const double crLHS123 = crLHS122*crLHS51;
const double crLHS124 = 1.0*crLHS14 + 1.0*crLHS15 + 1.0*crLHS16;
const double crLHS125 = crLHS124*crLHS56;
const double crLHS126 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double crLHS127 = pow(DN(0,1), 2);
const double crLHS128 = crLHS30 - crLHS53;
const double crLHS129 = crLHS128*crLHS33;
const double crLHS130 = -crLHS128*crLHS18 + crLHS38;
const double crLHS131 = DN(0,1)*N[0];
const double crLHS132 = DN(0,1)*crLHS33;
const double crLHS133 = DN(0,0)*crLHS17 + DN(0,1)*crLHS18;
const double crLHS134 = C(0,1)*DN(1,0) + crLHS81;
const double crLHS135 = DN(0,1)*crLHS20;
const double crLHS136 = DN(1,0)*crLHS135;
const double crLHS137 = crLHS17*crLHS24;
const double crLHS138 = crLHS122*rho;
const double crLHS139 = N[1]*crLHS137 - crLHS138*crLHS76;
const double crLHS140 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double crLHS141 = crLHS71 - crLHS84;
const double crLHS142 = crLHS141*crLHS33;
const double crLHS143 = -crLHS141*crLHS18 + crLHS74;
const double crLHS144 = DN(0,1)*DN(1,1);
const double crLHS145 = N[1]*crLHS53 + crLHS144*crLHS20;
const double crLHS146 = DN(0,1)*N[1];
const double crLHS147 = DN(1,1)*crLHS33;
const double crLHS148 = DN(1,0)*crLHS17 + DN(1,1)*crLHS18;
const double crLHS149 = C(0,1)*DN(2,0) + crLHS112;
const double crLHS150 = DN(2,0)*crLHS135;
const double crLHS151 = N[2]*crLHS137 - crLHS107*crLHS138;
const double crLHS152 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double crLHS153 = crLHS102 - crLHS115;
const double crLHS154 = crLHS153*crLHS33;
const double crLHS155 = crLHS105 - crLHS153*crLHS18;
const double crLHS156 = DN(0,1)*DN(2,1);
const double crLHS157 = N[2]*crLHS53 + crLHS156*crLHS20;
const double crLHS158 = DN(0,1)*N[2];
const double crLHS159 = DN(2,1)*crLHS33;
const double crLHS160 = DN(2,0)*crLHS17 + DN(2,1)*crLHS18;
const double crLHS161 = crLHS33*gauss_weight;
const double crLHS162 = DN(1,0)*N[0];
const double crLHS163 = DN(1,1)*N[0];
const double crLHS164 = crLHS161*(crLHS144 + crLHS78);
const double crLHS165 = DN(2,0)*N[0];
const double crLHS166 = DN(2,1)*N[0];
const double crLHS167 = crLHS161*(crLHS109 + crLHS156);
const double crLHS168 = crLHS67 + crLHS69;
const double crLHS169 = crLHS168*rho;
const double crLHS170 = crLHS33*crLHS64;
const double crLHS171 = crLHS168*crLHS24 + crLHS76;
const double crLHS172 = crLHS32*crLHS64;
const double crLHS173 = crLHS172*crLHS55;
const double crLHS174 = crLHS168*crLHS50;
const double crLHS175 = crLHS174*crLHS48;
const double crLHS176 = pow(DN(1,0), 2);
const double crLHS177 = pow(N[1], 2);
const double crLHS178 = crLHS177*rho;
const double crLHS179 = crLHS177*crLHS4;
const double crLHS180 = crLHS168*crLHS64 + crLHS179;
const double crLHS181 = DN(1,0)*crLHS20;
const double crLHS182 = DN(1,1)*crLHS181;
const double crLHS183 = DN(1,0)*N[1];
const double crLHS184 = N[2]*crLHS66;
const double crLHS185 = crLHS168*crLHS95 + crLHS184;
const double crLHS186 = DN(1,0)*DN(2,0);
const double crLHS187 = N[2]*crLHS65 + crLHS186*crLHS20;
const double crLHS188 = DN(2,1)*crLHS181;
const double crLHS189 = N[2]*crLHS64;
const double crLHS190 = crLHS33*crLHS95;
const double crLHS191 = crLHS190*crLHS66;
const double crLHS192 = crLHS13*crLHS189 - crLHS13*crLHS191;
const double crLHS193 = DN(1,0)*N[2];
const double crLHS194 = crLHS124*crLHS172;
const double crLHS195 = crLHS122*crLHS174;
const double crLHS196 = pow(DN(1,1), 2);
const double crLHS197 = DN(1,1)*N[1];
const double crLHS198 = DN(2,0)*crLHS20;
const double crLHS199 = DN(1,1)*crLHS198;
const double crLHS200 = crLHS17*crLHS189 - crLHS17*crLHS191;
const double crLHS201 = DN(1,1)*DN(2,1);
const double crLHS202 = N[2]*crLHS84 + crLHS20*crLHS201;
const double crLHS203 = DN(1,1)*N[2];
const double crLHS204 = DN(2,0)*N[1];
const double crLHS205 = DN(2,1)*N[1];
const double crLHS206 = crLHS161*(crLHS186 + crLHS201);
const double crLHS207 = crLHS100 + crLHS98;
const double crLHS208 = crLHS207*rho;
const double crLHS209 = crLHS107 + crLHS207*crLHS24;
const double crLHS210 = crLHS32*crLHS95;
const double crLHS211 = crLHS210*crLHS55;
const double crLHS212 = crLHS207*crLHS50;
const double crLHS213 = crLHS212*crLHS48;
const double crLHS214 = crLHS184 + crLHS207*crLHS64;
const double crLHS215 = pow(DN(2,0), 2);
const double crLHS216 = pow(N[2], 2);
const double crLHS217 = crLHS216*rho;
const double crLHS218 = crLHS216*crLHS4;
const double crLHS219 = crLHS207*crLHS95 + crLHS218;
const double crLHS220 = DN(2,1)*crLHS198;
const double crLHS221 = DN(2,0)*N[2];
const double crLHS222 = crLHS124*crLHS210;
const double crLHS223 = crLHS122*crLHS212;
const double crLHS224 = pow(DN(2,1), 2);
const double crLHS225 = DN(2,1)*N[2];
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + crLHS20*crLHS3 + crLHS22*crLHS9 + crLHS23*crLHS34 + crLHS34*crLHS36 - crLHS39*crLHS40 + crLHS42);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS43 + DN(0,1)*crLHS45 - N[0]*crLHS52 + crLHS13*crLHS22 + crLHS47 - crLHS48*crLHS49 - crLHS54*crLHS57);
rLHS(0,2)+=-gauss_weight*(crLHS23*crLHS59 + crLHS36*crLHS59 + crLHS40*crLHS60 + crLHS58);
rLHS(0,3)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + crLHS23*crLHS73 + crLHS36*crLHS73 - crLHS40*crLHS75 + crLHS77 + crLHS79);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS80 + DN(0,1)*crLHS82 - N[1]*crLHS52 - crLHS57*crLHS85 + crLHS83 + crLHS88);
rLHS(0,5)+=-gauss_weight*(crLHS23*crLHS90 + crLHS36*crLHS90 + crLHS40*crLHS91 + crLHS89);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS92 + DN(0,1)*crLHS94 + crLHS104*crLHS23 + crLHS104*crLHS36 - crLHS106*crLHS40 + crLHS108 + crLHS110);
rLHS(0,7)+=gauss_weight*(DN(0,0)*crLHS111 + DN(0,1)*crLHS113 - N[2]*crLHS52 + crLHS114 - crLHS116*crLHS57 + crLHS117);
rLHS(0,8)+=-gauss_weight*(crLHS118 + crLHS119*crLHS23 + crLHS119*crLHS36 + crLHS120*crLHS40);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS121 - N[0]*crLHS123 - crLHS122*crLHS49 - crLHS125*crLHS54 + crLHS17*crLHS22 + crLHS47);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS126 + crLHS127*crLHS20 + crLHS129*crLHS23 + crLHS129*crLHS36 - crLHS130*crLHS40 + crLHS18*crLHS22 + crLHS42);
rLHS(1,2)+=-gauss_weight*(crLHS131 + crLHS132*crLHS23 + crLHS132*crLHS36 + crLHS133*crLHS40);
rLHS(1,3)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS134 - N[1]*crLHS123 - crLHS125*crLHS85 + crLHS136 + crLHS139);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS82 + DN(0,1)*crLHS140 + crLHS142*crLHS23 + crLHS142*crLHS36 - crLHS143*crLHS40 + crLHS145 + crLHS77);
rLHS(1,5)+=-gauss_weight*(crLHS146 + crLHS147*crLHS23 + crLHS147*crLHS36 + crLHS148*crLHS40);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS94 + DN(0,1)*crLHS149 - N[2]*crLHS123 - crLHS116*crLHS125 + crLHS150 + crLHS151);
rLHS(1,7)+=gauss_weight*(DN(0,0)*crLHS113 + DN(0,1)*crLHS152 + crLHS108 + crLHS154*crLHS23 + crLHS154*crLHS36 - crLHS155*crLHS40 + crLHS157);
rLHS(1,8)+=-gauss_weight*(crLHS158 + crLHS159*crLHS23 + crLHS159*crLHS36 + crLHS160*crLHS40);
rLHS(2,0)+=gauss_weight*(crLHS131*crLHS138 - crLHS31*crLHS59 + crLHS58);
rLHS(2,1)+=gauss_weight*(-crLHS128*crLHS132 + crLHS131 + crLHS58*crLHS87);
rLHS(2,2)+=crLHS161*(crLHS127 + crLHS3);
rLHS(2,3)+=gauss_weight*(crLHS138*crLHS146 + crLHS162 - crLHS59*crLHS72);
rLHS(2,4)+=gauss_weight*(-crLHS132*crLHS141 + crLHS163 + crLHS87*crLHS89);
rLHS(2,5)+=crLHS164;
rLHS(2,6)+=gauss_weight*(-crLHS103*crLHS59 + crLHS138*crLHS158 + crLHS165);
rLHS(2,7)+=gauss_weight*(crLHS118*crLHS87 - crLHS132*crLHS153 + crLHS166);
rLHS(2,8)+=crLHS167;
rLHS(3,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + crLHS169*crLHS34 - crLHS170*crLHS39 + crLHS171 + crLHS34*crLHS66 + crLHS79);
rLHS(3,1)+=gauss_weight*(DN(1,0)*crLHS43 + DN(1,1)*crLHS45 - N[0]*crLHS175 + crLHS136 - crLHS173*crLHS54 + crLHS88);
rLHS(3,2)+=-gauss_weight*(crLHS162 + crLHS169*crLHS59 + crLHS170*crLHS60 + crLHS59*crLHS66);
rLHS(3,3)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + crLHS169*crLHS73 - crLHS170*crLHS75 + crLHS176*crLHS20 + crLHS178*crLHS9 + crLHS180 + crLHS66*crLHS73);
rLHS(3,4)+=gauss_weight*(DN(1,0)*crLHS80 + DN(1,1)*crLHS82 - N[1]*crLHS175 + crLHS13*crLHS178 - crLHS173*crLHS85 - crLHS179*crLHS87 + crLHS182);
rLHS(3,5)+=-gauss_weight*(crLHS169*crLHS90 + crLHS170*crLHS91 + crLHS183 + crLHS66*crLHS90);
rLHS(3,6)+=gauss_weight*(DN(1,0)*crLHS92 + DN(1,1)*crLHS94 + crLHS104*crLHS169 + crLHS104*crLHS66 - crLHS106*crLHS170 + crLHS185 + crLHS187);
rLHS(3,7)+=gauss_weight*(DN(1,0)*crLHS111 + DN(1,1)*crLHS113 - N[2]*crLHS175 - crLHS116*crLHS173 + crLHS188 + crLHS192);
rLHS(3,8)+=-gauss_weight*(crLHS119*crLHS169 + crLHS119*crLHS66 + crLHS120*crLHS170 + crLHS193);
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS121 - N[0]*crLHS195 + crLHS139 - crLHS194*crLHS54 + crLHS83);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS126 + crLHS129*crLHS169 + crLHS129*crLHS66 - crLHS130*crLHS170 + crLHS145 + crLHS171);
rLHS(4,2)+=-gauss_weight*(crLHS132*crLHS169 + crLHS132*crLHS66 + crLHS133*crLHS170 + crLHS163);
rLHS(4,3)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS134 - N[1]*crLHS195 - crLHS138*crLHS179 + crLHS17*crLHS178 + crLHS182 - crLHS194*crLHS85);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS82 + DN(1,1)*crLHS140 + crLHS142*crLHS169 + crLHS142*crLHS66 - crLHS143*crLHS170 + crLHS178*crLHS18 + crLHS180 + crLHS196*crLHS20);
rLHS(4,5)+=-gauss_weight*(crLHS147*crLHS169 + crLHS147*crLHS66 + crLHS148*crLHS170 + crLHS197);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS94 + DN(1,1)*crLHS149 - N[2]*crLHS195 - crLHS116*crLHS194 + crLHS199 + crLHS200);
rLHS(4,7)+=gauss_weight*(DN(1,0)*crLHS113 + DN(1,1)*crLHS152 + crLHS154*crLHS169 + crLHS154*crLHS66 - crLHS155*crLHS170 + crLHS185 + crLHS202);
rLHS(4,8)+=-gauss_weight*(crLHS159*crLHS169 + crLHS159*crLHS66 + crLHS160*crLHS170 + crLHS203);
rLHS(5,0)+=gauss_weight*(crLHS138*crLHS163 - crLHS31*crLHS90 + crLHS89);
rLHS(5,1)+=gauss_weight*(-crLHS128*crLHS147 + crLHS146 + crLHS162*crLHS87);
rLHS(5,2)+=crLHS164;
rLHS(5,3)+=gauss_weight*(crLHS138*crLHS197 + crLHS183 - crLHS72*crLHS90);
rLHS(5,4)+=gauss_weight*(-crLHS141*crLHS147 + crLHS183*crLHS87 + crLHS197);
rLHS(5,5)+=crLHS161*(crLHS176 + crLHS196);
rLHS(5,6)+=gauss_weight*(-crLHS103*crLHS90 + crLHS138*crLHS203 + crLHS204);
rLHS(5,7)+=gauss_weight*(-crLHS147*crLHS153 + crLHS193*crLHS87 + crLHS205);
rLHS(5,8)+=crLHS206;
rLHS(6,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + crLHS110 - crLHS190*crLHS39 + crLHS208*crLHS34 + crLHS209 + crLHS34*crLHS97);
rLHS(6,1)+=gauss_weight*(DN(2,0)*crLHS43 + DN(2,1)*crLHS45 - N[0]*crLHS213 + crLHS117 + crLHS150 - crLHS211*crLHS54);
rLHS(6,2)+=-gauss_weight*(crLHS165 + crLHS190*crLHS60 + crLHS208*crLHS59 + crLHS59*crLHS97);
rLHS(6,3)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + crLHS187 - crLHS190*crLHS75 + crLHS208*crLHS73 + crLHS214 + crLHS73*crLHS97);
rLHS(6,4)+=gauss_weight*(DN(2,0)*crLHS80 + DN(2,1)*crLHS82 - N[1]*crLHS213 + crLHS192 + crLHS199 - crLHS211*crLHS85);
rLHS(6,5)+=-gauss_weight*(crLHS190*crLHS91 + crLHS204 + crLHS208*crLHS90 + crLHS90*crLHS97);
rLHS(6,6)+=gauss_weight*(DN(2,0)*crLHS92 + DN(2,1)*crLHS94 + crLHS104*crLHS208 + crLHS104*crLHS97 - crLHS106*crLHS190 + crLHS20*crLHS215 + crLHS217*crLHS9 + crLHS219);
rLHS(6,7)+=gauss_weight*(DN(2,0)*crLHS111 + DN(2,1)*crLHS113 - N[2]*crLHS213 - crLHS116*crLHS211 + crLHS13*crLHS217 - crLHS218*crLHS87 + crLHS220);
rLHS(6,8)+=-gauss_weight*(crLHS119*crLHS208 + crLHS119*crLHS97 + crLHS120*crLHS190 + crLHS221);
rLHS(7,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS121 - N[0]*crLHS223 + crLHS114 + crLHS151 - crLHS222*crLHS54);
rLHS(7,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS126 + crLHS129*crLHS208 + crLHS129*crLHS97 - crLHS130*crLHS190 + crLHS157 + crLHS209);
rLHS(7,2)+=-gauss_weight*(crLHS132*crLHS208 + crLHS132*crLHS97 + crLHS133*crLHS190 + crLHS166);
rLHS(7,3)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS134 - N[1]*crLHS223 + crLHS188 + crLHS200 - crLHS222*crLHS85);
rLHS(7,4)+=gauss_weight*(DN(2,0)*crLHS82 + DN(2,1)*crLHS140 + crLHS142*crLHS208 + crLHS142*crLHS97 - crLHS143*crLHS190 + crLHS202 + crLHS214);
rLHS(7,5)+=-gauss_weight*(crLHS147*crLHS208 + crLHS147*crLHS97 + crLHS148*crLHS190 + crLHS205);
rLHS(7,6)+=gauss_weight*(DN(2,0)*crLHS94 + DN(2,1)*crLHS149 - N[2]*crLHS223 - crLHS116*crLHS222 - crLHS138*crLHS218 + crLHS17*crLHS217 + crLHS220);
rLHS(7,7)+=gauss_weight*(DN(2,0)*crLHS113 + DN(2,1)*crLHS152 + crLHS154*crLHS208 + crLHS154*crLHS97 - crLHS155*crLHS190 + crLHS18*crLHS217 + crLHS20*crLHS224 + crLHS219);
rLHS(7,8)+=-gauss_weight*(crLHS159*crLHS208 + crLHS159*crLHS97 + crLHS160*crLHS190 + crLHS225);
rLHS(8,0)+=gauss_weight*(crLHS118 - crLHS119*crLHS31 + crLHS138*crLHS166);
rLHS(8,1)+=gauss_weight*(-crLHS128*crLHS159 + crLHS158 + crLHS165*crLHS87);
rLHS(8,2)+=crLHS167;
rLHS(8,3)+=gauss_weight*(-crLHS119*crLHS72 + crLHS138*crLHS205 + crLHS193);
rLHS(8,4)+=gauss_weight*(-crLHS141*crLHS159 + crLHS203 + crLHS204*crLHS87);
rLHS(8,5)+=crLHS206;
rLHS(8,6)+=gauss_weight*(-crLHS103*crLHS119 + crLHS138*crLHS225 + crLHS221);
rLHS(8,7)+=gauss_weight*(-crLHS153*crLHS159 + crLHS221*crLHS87 + crLHS225);
rLHS(8,8)+=crLHS161*(crLHS215 + crLHS224);
    
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
const double crLHS44 = crLHS26*crLHS35 + crLHS43;
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
const double crLHS94 = crLHS35*crLHS79 + crLHS93;
const double crLHS95 = DN(0,0)*DN(1,0);
const double crLHS96 = N[1]*crLHS27 + crLHS22*crLHS95;
const double crLHS97 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crLHS75;
const double crLHS98 = C(1,3)*DN(1,1);
const double crLHS99 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crLHS98;
const double crLHS100 = C(3,5)*DN(1,0);
const double crLHS101 = C(4,5)*DN(1,2);
const double crLHS102 = C(1,5)*DN(1,1) + crLHS100 + crLHS101;
const double crLHS103 = DN(1,1)*crLHS51;
const double crLHS104 = crLHS16*crLHS79;
const double crLHS105 = -crLHS104 + crLHS85;
const double crLHS106 = -crLHS105*crLHS13 + crLHS13*crLHS80 + crLHS19*crLHS90;
const double crLHS107 = N[1]*crLHS56;
const double crLHS108 = crLHS93*rho;
const double crLHS109 = N[1]*crLHS37 - crLHS108*crLHS53;
const double crLHS110 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crLHS77;
const double crLHS111 = C(3,4)*DN(1,1);
const double crLHS112 = C(2,3)*DN(1,2) + crLHS100 + crLHS111;
const double crLHS113 = C(2,5)*DN(1,2);
const double crLHS114 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crLHS113;
const double crLHS115 = DN(1,2)*crLHS51;
const double crLHS116 = crLHS20*crLHS79;
const double crLHS117 = -crLHS116 + crLHS85;
const double crLHS118 = -crLHS117*crLHS14 + crLHS14*crLHS80 + crLHS17*crLHS88;
const double crLHS119 = N[1]*crLHS39 - crLHS108*crLHS67;
const double crLHS120 = DN(0,0)*N[1];
const double crLHS121 = DN(1,0)*crLHS33;
const double crLHS122 = DN(1,0)*crLHS12 + DN(1,1)*crLHS13 + DN(1,2)*crLHS14;
const double crLHS123 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double crLHS124 = C(0,3)*DN(2,0);
const double crLHS125 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crLHS124;
const double crLHS126 = C(0,5)*DN(2,0);
const double crLHS127 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crLHS126;
const double crLHS128 = N[2]*rho;
const double crLHS129 = crLHS12*crLHS128;
const double crLHS130 = N[2]*crLHS6;
const double crLHS131 = DN(2,0)*crLHS8;
const double crLHS132 = DN(2,1)*crLHS9;
const double crLHS133 = DN(2,2)*crLHS10;
const double crLHS134 = -crLHS130 + crLHS131*rho + crLHS132*rho + crLHS133*rho;
const double crLHS135 = -crLHS129 + crLHS134;
const double crLHS136 = crLHS135*crLHS33;
const double crLHS137 = crLHS128*crLHS13;
const double crLHS138 = crLHS137*crLHS15;
const double crLHS139 = crLHS128*crLHS14;
const double crLHS140 = crLHS139*crLHS18;
const double crLHS141 = -crLHS12*crLHS135 + crLHS138 + crLHS140;
const double crLHS142 = N[2]*crLHS25;
const double crLHS143 = crLHS128*crLHS35 + crLHS142;
const double crLHS144 = DN(0,0)*DN(2,0);
const double crLHS145 = N[2]*crLHS27 + crLHS144*crLHS22;
const double crLHS146 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crLHS124;
const double crLHS147 = C(1,3)*DN(2,1);
const double crLHS148 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crLHS147;
const double crLHS149 = C(3,5)*DN(2,0);
const double crLHS150 = C(4,5)*DN(2,2);
const double crLHS151 = C(1,5)*DN(2,1) + crLHS149 + crLHS150;
const double crLHS152 = DN(2,1)*crLHS51;
const double crLHS153 = crLHS128*crLHS16;
const double crLHS154 = crLHS134 - crLHS153;
const double crLHS155 = crLHS129*crLHS13 - crLHS13*crLHS154 + crLHS139*crLHS19;
const double crLHS156 = N[2]*crLHS56;
const double crLHS157 = crLHS142*rho;
const double crLHS158 = N[2]*crLHS37 - crLHS157*crLHS53;
const double crLHS159 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crLHS126;
const double crLHS160 = C(3,4)*DN(2,1);
const double crLHS161 = C(2,3)*DN(2,2) + crLHS149 + crLHS160;
const double crLHS162 = C(2,5)*DN(2,2);
const double crLHS163 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crLHS162;
const double crLHS164 = DN(2,2)*crLHS51;
const double crLHS165 = crLHS128*crLHS20;
const double crLHS166 = crLHS134 - crLHS165;
const double crLHS167 = crLHS129*crLHS14 + crLHS137*crLHS17 - crLHS14*crLHS166;
const double crLHS168 = N[2]*crLHS39 - crLHS157*crLHS67;
const double crLHS169 = DN(0,0)*N[2];
const double crLHS170 = DN(2,0)*crLHS33;
const double crLHS171 = DN(2,0)*crLHS12 + DN(2,1)*crLHS13 + DN(2,2)*crLHS14;
const double crLHS172 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double crLHS173 = C(0,3)*DN(3,0);
const double crLHS174 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crLHS173;
const double crLHS175 = C(0,5)*DN(3,0);
const double crLHS176 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crLHS175;
const double crLHS177 = N[3]*rho;
const double crLHS178 = crLHS12*crLHS177;
const double crLHS179 = N[3]*crLHS6;
const double crLHS180 = DN(3,0)*crLHS8;
const double crLHS181 = DN(3,1)*crLHS9;
const double crLHS182 = DN(3,2)*crLHS10;
const double crLHS183 = -crLHS179 + crLHS180*rho + crLHS181*rho + crLHS182*rho;
const double crLHS184 = -crLHS178 + crLHS183;
const double crLHS185 = crLHS184*crLHS33;
const double crLHS186 = crLHS13*crLHS177;
const double crLHS187 = crLHS15*crLHS186;
const double crLHS188 = crLHS14*crLHS177;
const double crLHS189 = crLHS18*crLHS188;
const double crLHS190 = -crLHS12*crLHS184 + crLHS187 + crLHS189;
const double crLHS191 = N[3]*crLHS25;
const double crLHS192 = crLHS177*crLHS35 + crLHS191;
const double crLHS193 = DN(0,0)*DN(3,0);
const double crLHS194 = N[3]*crLHS27 + crLHS193*crLHS22;
const double crLHS195 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crLHS173;
const double crLHS196 = C(1,3)*DN(3,1);
const double crLHS197 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crLHS196;
const double crLHS198 = C(3,5)*DN(3,0);
const double crLHS199 = C(4,5)*DN(3,2);
const double crLHS200 = C(1,5)*DN(3,1) + crLHS198 + crLHS199;
const double crLHS201 = DN(3,1)*crLHS51;
const double crLHS202 = crLHS16*crLHS177;
const double crLHS203 = crLHS183 - crLHS202;
const double crLHS204 = crLHS13*crLHS178 - crLHS13*crLHS203 + crLHS188*crLHS19;
const double crLHS205 = N[3]*crLHS56;
const double crLHS206 = crLHS191*rho;
const double crLHS207 = N[3]*crLHS37 - crLHS206*crLHS53;
const double crLHS208 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crLHS175;
const double crLHS209 = C(3,4)*DN(3,1);
const double crLHS210 = C(2,3)*DN(3,2) + crLHS198 + crLHS209;
const double crLHS211 = C(2,5)*DN(3,2);
const double crLHS212 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crLHS211;
const double crLHS213 = DN(3,2)*crLHS51;
const double crLHS214 = crLHS177*crLHS20;
const double crLHS215 = crLHS183 - crLHS214;
const double crLHS216 = crLHS14*crLHS178 - crLHS14*crLHS215 + crLHS17*crLHS186;
const double crLHS217 = N[3]*crLHS39 - crLHS206*crLHS67;
const double crLHS218 = DN(0,0)*N[3];
const double crLHS219 = DN(3,0)*crLHS33;
const double crLHS220 = DN(3,0)*crLHS12 + DN(3,1)*crLHS13 + DN(3,2)*crLHS14;
const double crLHS221 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crLHS46;
const double crLHS222 = C(0,4)*DN(0,0) + crLHS49 + crLHS62;
const double crLHS223 = crLHS15*crLHS33;
const double crLHS224 = crLHS17*crLHS26;
const double crLHS225 = -crLHS15*crLHS32 + crLHS15*crLHS58 + crLHS18*crLHS224;
const double crLHS226 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double crLHS227 = C(1,4)*DN(0,1);
const double crLHS228 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crLHS227;
const double crLHS229 = pow(DN(0,1), 2);
const double crLHS230 = crLHS33*crLHS59;
const double crLHS231 = crLHS19*crLHS224;
const double crLHS232 = -crLHS16*crLHS59 + crLHS231 + crLHS38;
const double crLHS233 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crLHS227;
const double crLHS234 = C(2,4)*DN(0,2);
const double crLHS235 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crLHS234;
const double crLHS236 = DN(0,1)*crLHS22;
const double crLHS237 = DN(0,2)*crLHS236;
const double crLHS238 = crLHS17*crLHS33;
const double crLHS239 = crLHS15*crLHS39 + crLHS17*crLHS58 - crLHS17*crLHS69;
const double crLHS240 = DN(0,1)*N[0];
const double crLHS241 = DN(0,1)*crLHS33;
const double crLHS242 = DN(0,0)*crLHS15 + DN(0,1)*crLHS16 + DN(0,2)*crLHS17;
const double crLHS243 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crLHS98;
const double crLHS244 = C(0,4)*DN(1,0) + crLHS101 + crLHS111;
const double crLHS245 = DN(1,0)*crLHS236;
const double crLHS246 = crLHS17*crLHS79;
const double crLHS247 = crLHS104*crLHS15 - crLHS15*crLHS86 + crLHS18*crLHS246;
const double crLHS248 = crLHS15*crLHS26;
const double crLHS249 = N[1]*crLHS248 - crLHS108*crLHS223;
const double crLHS250 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double crLHS251 = C(1,4)*DN(1,1);
const double crLHS252 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crLHS251;
const double crLHS253 = crLHS105*crLHS33;
const double crLHS254 = crLHS19*crLHS246;
const double crLHS255 = -crLHS105*crLHS16 + crLHS254 + crLHS89;
const double crLHS256 = DN(0,1)*DN(1,1);
const double crLHS257 = N[1]*crLHS58 + crLHS22*crLHS256;
const double crLHS258 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crLHS251;
const double crLHS259 = C(2,4)*DN(1,2);
const double crLHS260 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crLHS259;
const double crLHS261 = DN(1,2)*crLHS236;
const double crLHS262 = crLHS104*crLHS17 - crLHS117*crLHS17 + crLHS15*crLHS90;
const double crLHS263 = N[1]*crLHS224 - crLHS108*crLHS238;
const double crLHS264 = DN(0,1)*N[1];
const double crLHS265 = DN(1,1)*crLHS33;
const double crLHS266 = DN(1,0)*crLHS15 + DN(1,1)*crLHS16 + DN(1,2)*crLHS17;
const double crLHS267 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crLHS147;
const double crLHS268 = C(0,4)*DN(2,0) + crLHS150 + crLHS160;
const double crLHS269 = DN(2,0)*crLHS236;
const double crLHS270 = crLHS128*crLHS17;
const double crLHS271 = -crLHS135*crLHS15 + crLHS15*crLHS153 + crLHS18*crLHS270;
const double crLHS272 = N[2]*crLHS248 - crLHS157*crLHS223;
const double crLHS273 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double crLHS274 = C(1,4)*DN(2,1);
const double crLHS275 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crLHS274;
const double crLHS276 = crLHS154*crLHS33;
const double crLHS277 = crLHS19*crLHS270;
const double crLHS278 = crLHS138 - crLHS154*crLHS16 + crLHS277;
const double crLHS279 = DN(0,1)*DN(2,1);
const double crLHS280 = N[2]*crLHS58 + crLHS22*crLHS279;
const double crLHS281 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crLHS274;
const double crLHS282 = C(2,4)*DN(2,2);
const double crLHS283 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crLHS282;
const double crLHS284 = DN(2,2)*crLHS236;
const double crLHS285 = crLHS139*crLHS15 + crLHS153*crLHS17 - crLHS166*crLHS17;
const double crLHS286 = N[2]*crLHS224 - crLHS157*crLHS238;
const double crLHS287 = DN(0,1)*N[2];
const double crLHS288 = DN(2,1)*crLHS33;
const double crLHS289 = DN(2,0)*crLHS15 + DN(2,1)*crLHS16 + DN(2,2)*crLHS17;
const double crLHS290 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crLHS196;
const double crLHS291 = C(0,4)*DN(3,0) + crLHS199 + crLHS209;
const double crLHS292 = DN(3,0)*crLHS236;
const double crLHS293 = crLHS17*crLHS177;
const double crLHS294 = -crLHS15*crLHS184 + crLHS15*crLHS202 + crLHS18*crLHS293;
const double crLHS295 = N[3]*crLHS248 - crLHS206*crLHS223;
const double crLHS296 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double crLHS297 = C(1,4)*DN(3,1);
const double crLHS298 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crLHS297;
const double crLHS299 = crLHS203*crLHS33;
const double crLHS300 = crLHS19*crLHS293;
const double crLHS301 = -crLHS16*crLHS203 + crLHS187 + crLHS300;
const double crLHS302 = DN(0,1)*DN(3,1);
const double crLHS303 = N[3]*crLHS58 + crLHS22*crLHS302;
const double crLHS304 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crLHS297;
const double crLHS305 = C(2,4)*DN(3,2);
const double crLHS306 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crLHS305;
const double crLHS307 = DN(3,2)*crLHS236;
const double crLHS308 = crLHS15*crLHS188 + crLHS17*crLHS202 - crLHS17*crLHS215;
const double crLHS309 = N[3]*crLHS224 - crLHS206*crLHS238;
const double crLHS310 = DN(0,1)*N[3];
const double crLHS311 = DN(3,1)*crLHS33;
const double crLHS312 = DN(3,0)*crLHS15 + DN(3,1)*crLHS16 + DN(3,2)*crLHS17;
const double crLHS313 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crLHS64;
const double crLHS314 = crLHS18*crLHS33;
const double crLHS315 = -crLHS18*crLHS32 + crLHS18*crLHS68 + crLHS19*crLHS248;
const double crLHS316 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crLHS234;
const double crLHS317 = crLHS19*crLHS33;
const double crLHS318 = crLHS18*crLHS37 - crLHS19*crLHS59 + crLHS19*crLHS68;
const double crLHS319 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double crLHS320 = pow(DN(0,2), 2);
const double crLHS321 = crLHS33*crLHS69;
const double crLHS322 = -crLHS20*crLHS69 + crLHS231 + crLHS40;
const double crLHS323 = DN(0,2)*N[0];
const double crLHS324 = DN(0,2)*crLHS33;
const double crLHS325 = DN(0,0)*crLHS18 + DN(0,1)*crLHS19 + DN(0,2)*crLHS20;
const double crLHS326 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crLHS113;
const double crLHS327 = DN(0,2)*crLHS22;
const double crLHS328 = DN(1,0)*crLHS327;
const double crLHS329 = crLHS15*crLHS19;
const double crLHS330 = crLHS116*crLHS18 - crLHS18*crLHS86 + crLHS329*crLHS79;
const double crLHS331 = N[1]*crLHS26;
const double crLHS332 = -crLHS108*crLHS314 + crLHS18*crLHS331;
const double crLHS333 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crLHS259;
const double crLHS334 = DN(1,1)*crLHS327;
const double crLHS335 = -crLHS105*crLHS19 + crLHS116*crLHS19 + crLHS18*crLHS88;
const double crLHS336 = -crLHS108*crLHS317 + crLHS19*crLHS331;
const double crLHS337 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double crLHS338 = crLHS117*crLHS33;
const double crLHS339 = -crLHS117*crLHS20 + crLHS254 + crLHS91;
const double crLHS340 = DN(0,2)*DN(1,2);
const double crLHS341 = N[1]*crLHS68 + crLHS22*crLHS340;
const double crLHS342 = DN(0,2)*N[1];
const double crLHS343 = DN(1,2)*crLHS33;
const double crLHS344 = DN(1,0)*crLHS18 + DN(1,1)*crLHS19 + DN(1,2)*crLHS20;
const double crLHS345 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crLHS162;
const double crLHS346 = DN(2,0)*crLHS327;
const double crLHS347 = crLHS128*crLHS329 - crLHS135*crLHS18 + crLHS165*crLHS18;
const double crLHS348 = N[2]*crLHS26;
const double crLHS349 = -crLHS157*crLHS314 + crLHS18*crLHS348;
const double crLHS350 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crLHS282;
const double crLHS351 = DN(2,1)*crLHS327;
const double crLHS352 = crLHS137*crLHS18 - crLHS154*crLHS19 + crLHS165*crLHS19;
const double crLHS353 = -crLHS157*crLHS317 + crLHS19*crLHS348;
const double crLHS354 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double crLHS355 = crLHS166*crLHS33;
const double crLHS356 = crLHS140 - crLHS166*crLHS20 + crLHS277;
const double crLHS357 = DN(0,2)*DN(2,2);
const double crLHS358 = N[2]*crLHS68 + crLHS22*crLHS357;
const double crLHS359 = DN(0,2)*N[2];
const double crLHS360 = DN(2,2)*crLHS33;
const double crLHS361 = DN(2,0)*crLHS18 + DN(2,1)*crLHS19 + DN(2,2)*crLHS20;
const double crLHS362 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crLHS211;
const double crLHS363 = DN(3,0)*crLHS327;
const double crLHS364 = crLHS177*crLHS329 - crLHS18*crLHS184 + crLHS18*crLHS214;
const double crLHS365 = N[3]*crLHS26;
const double crLHS366 = crLHS18*crLHS365 - crLHS206*crLHS314;
const double crLHS367 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crLHS305;
const double crLHS368 = DN(3,1)*crLHS327;
const double crLHS369 = crLHS18*crLHS186 - crLHS19*crLHS203 + crLHS19*crLHS214;
const double crLHS370 = crLHS19*crLHS365 - crLHS206*crLHS317;
const double crLHS371 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double crLHS372 = crLHS215*crLHS33;
const double crLHS373 = crLHS189 - crLHS20*crLHS215 + crLHS300;
const double crLHS374 = DN(0,2)*DN(3,2);
const double crLHS375 = N[3]*crLHS68 + crLHS22*crLHS374;
const double crLHS376 = DN(0,2)*N[3];
const double crLHS377 = DN(3,2)*crLHS33;
const double crLHS378 = DN(3,0)*crLHS18 + DN(3,1)*crLHS19 + DN(3,2)*crLHS20;
const double crLHS379 = crLHS240*rho;
const double crLHS380 = crLHS323*rho;
const double crLHS381 = crLHS71*rho;
const double crLHS382 = crLHS33*gauss_weight;
const double crLHS383 = DN(1,0)*N[0];
const double crLHS384 = crLHS264*rho;
const double crLHS385 = crLHS342*rho;
const double crLHS386 = DN(1,1)*N[0];
const double crLHS387 = crLHS120*rho;
const double crLHS388 = DN(1,2)*N[0];
const double crLHS389 = crLHS382*(crLHS256 + crLHS340 + crLHS95);
const double crLHS390 = DN(2,0)*N[0];
const double crLHS391 = crLHS287*rho;
const double crLHS392 = crLHS359*rho;
const double crLHS393 = DN(2,1)*N[0];
const double crLHS394 = crLHS169*rho;
const double crLHS395 = DN(2,2)*N[0];
const double crLHS396 = crLHS382*(crLHS144 + crLHS279 + crLHS357);
const double crLHS397 = DN(3,0)*N[0];
const double crLHS398 = crLHS310*rho;
const double crLHS399 = crLHS376*rho;
const double crLHS400 = DN(3,1)*N[0];
const double crLHS401 = crLHS218*rho;
const double crLHS402 = DN(3,2)*N[0];
const double crLHS403 = crLHS382*(crLHS193 + crLHS302 + crLHS374);
const double crLHS404 = crLHS82 + crLHS83 + crLHS84;
const double crLHS405 = crLHS404*rho;
const double crLHS406 = crLHS33*crLHS79;
const double crLHS407 = crLHS26*crLHS404 + crLHS93;
const double crLHS408 = crLHS404*crLHS55;
const double crLHS409 = N[0]*crLHS408;
const double crLHS410 = pow(DN(1,0), 2);
const double crLHS411 = pow(N[1], 2);
const double crLHS412 = crLHS411*rho;
const double crLHS413 = crLHS411*crLHS6;
const double crLHS414 = crLHS404*crLHS79 + crLHS413;
const double crLHS415 = DN(1,0)*crLHS22;
const double crLHS416 = DN(1,1)*crLHS415;
const double crLHS417 = crLHS413*rho;
const double crLHS418 = N[1]*crLHS408;
const double crLHS419 = DN(1,2)*crLHS415;
const double crLHS420 = DN(1,0)*N[1];
const double crLHS421 = N[2]*crLHS81;
const double crLHS422 = crLHS128*crLHS404 + crLHS421;
const double crLHS423 = DN(1,0)*DN(2,0);
const double crLHS424 = N[2]*crLHS80 + crLHS22*crLHS423;
const double crLHS425 = DN(2,1)*crLHS415;
const double crLHS426 = N[2]*crLHS408;
const double crLHS427 = crLHS33*crLHS81;
const double crLHS428 = N[2]*crLHS88 - crLHS137*crLHS427;
const double crLHS429 = DN(2,2)*crLHS415;
const double crLHS430 = N[2]*crLHS90 - crLHS139*crLHS427;
const double crLHS431 = DN(1,0)*N[2];
const double crLHS432 = N[3]*crLHS81;
const double crLHS433 = crLHS177*crLHS404 + crLHS432;
const double crLHS434 = DN(1,0)*DN(3,0);
const double crLHS435 = N[3]*crLHS80 + crLHS22*crLHS434;
const double crLHS436 = DN(3,1)*crLHS415;
const double crLHS437 = N[3]*crLHS408;
const double crLHS438 = N[3]*crLHS88 - crLHS186*crLHS427;
const double crLHS439 = DN(3,2)*crLHS415;
const double crLHS440 = N[3]*crLHS90 - crLHS188*crLHS427;
const double crLHS441 = DN(1,0)*N[3];
const double crLHS442 = pow(DN(1,1), 2);
const double crLHS443 = DN(1,1)*crLHS22;
const double crLHS444 = DN(1,2)*crLHS443;
const double crLHS445 = DN(1,1)*N[1];
const double crLHS446 = DN(2,0)*crLHS443;
const double crLHS447 = crLHS15*crLHS79;
const double crLHS448 = crLHS128*crLHS33;
const double crLHS449 = crLHS15*crLHS81;
const double crLHS450 = N[2]*crLHS447 - crLHS448*crLHS449;
const double crLHS451 = DN(1,1)*DN(2,1);
const double crLHS452 = N[2]*crLHS104 + crLHS22*crLHS451;
const double crLHS453 = DN(2,2)*crLHS443;
const double crLHS454 = N[2]*crLHS246 - crLHS270*crLHS427;
const double crLHS455 = DN(1,1)*N[2];
const double crLHS456 = DN(3,0)*crLHS443;
const double crLHS457 = crLHS177*crLHS33;
const double crLHS458 = N[3]*crLHS447 - crLHS449*crLHS457;
const double crLHS459 = DN(1,1)*DN(3,1);
const double crLHS460 = N[3]*crLHS104 + crLHS22*crLHS459;
const double crLHS461 = DN(3,2)*crLHS443;
const double crLHS462 = N[3]*crLHS246 - crLHS293*crLHS427;
const double crLHS463 = DN(1,1)*N[3];
const double crLHS464 = pow(DN(1,2), 2);
const double crLHS465 = DN(1,2)*N[1];
const double crLHS466 = DN(1,2)*crLHS22;
const double crLHS467 = DN(2,0)*crLHS466;
const double crLHS468 = N[2]*crLHS79;
const double crLHS469 = crLHS448*crLHS81;
const double crLHS470 = crLHS18*crLHS468 - crLHS18*crLHS469;
const double crLHS471 = DN(2,1)*crLHS466;
const double crLHS472 = crLHS19*crLHS468 - crLHS19*crLHS469;
const double crLHS473 = DN(1,2)*DN(2,2);
const double crLHS474 = N[2]*crLHS116 + crLHS22*crLHS473;
const double crLHS475 = DN(1,2)*N[2];
const double crLHS476 = DN(3,0)*crLHS466;
const double crLHS477 = N[3]*crLHS79;
const double crLHS478 = crLHS457*crLHS81;
const double crLHS479 = crLHS18*crLHS477 - crLHS18*crLHS478;
const double crLHS480 = DN(3,1)*crLHS466;
const double crLHS481 = crLHS19*crLHS477 - crLHS19*crLHS478;
const double crLHS482 = DN(1,2)*DN(3,2);
const double crLHS483 = N[3]*crLHS116 + crLHS22*crLHS482;
const double crLHS484 = DN(1,2)*N[3];
const double crLHS485 = crLHS386*rho;
const double crLHS486 = crLHS388*rho;
const double crLHS487 = crLHS383*rho;
const double crLHS488 = crLHS445*rho;
const double crLHS489 = crLHS465*rho;
const double crLHS490 = crLHS420*rho;
const double crLHS491 = DN(2,0)*N[1];
const double crLHS492 = crLHS455*rho;
const double crLHS493 = crLHS475*rho;
const double crLHS494 = DN(2,1)*N[1];
const double crLHS495 = crLHS431*rho;
const double crLHS496 = DN(2,2)*N[1];
const double crLHS497 = crLHS382*(crLHS423 + crLHS451 + crLHS473);
const double crLHS498 = DN(3,0)*N[1];
const double crLHS499 = crLHS463*rho;
const double crLHS500 = crLHS484*rho;
const double crLHS501 = DN(3,1)*N[1];
const double crLHS502 = crLHS441*rho;
const double crLHS503 = DN(3,2)*N[1];
const double crLHS504 = crLHS382*(crLHS434 + crLHS459 + crLHS482);
const double crLHS505 = crLHS131 + crLHS132 + crLHS133;
const double crLHS506 = crLHS505*rho;
const double crLHS507 = crLHS142 + crLHS26*crLHS505;
const double crLHS508 = crLHS505*crLHS55;
const double crLHS509 = N[0]*crLHS508;
const double crLHS510 = crLHS421 + crLHS505*crLHS79;
const double crLHS511 = N[1]*crLHS508;
const double crLHS512 = pow(DN(2,0), 2);
const double crLHS513 = pow(N[2], 2);
const double crLHS514 = crLHS513*rho;
const double crLHS515 = crLHS513*crLHS6;
const double crLHS516 = crLHS128*crLHS505 + crLHS515;
const double crLHS517 = DN(2,0)*crLHS22;
const double crLHS518 = DN(2,1)*crLHS517;
const double crLHS519 = crLHS515*rho;
const double crLHS520 = N[2]*crLHS508;
const double crLHS521 = DN(2,2)*crLHS517;
const double crLHS522 = DN(2,0)*N[2];
const double crLHS523 = N[3]*crLHS130;
const double crLHS524 = crLHS177*crLHS505 + crLHS523;
const double crLHS525 = DN(2,0)*DN(3,0);
const double crLHS526 = N[3]*crLHS129 + crLHS22*crLHS525;
const double crLHS527 = DN(3,1)*crLHS517;
const double crLHS528 = N[3]*crLHS508;
const double crLHS529 = crLHS130*crLHS33;
const double crLHS530 = N[3]*crLHS137 - crLHS186*crLHS529;
const double crLHS531 = DN(3,2)*crLHS517;
const double crLHS532 = N[3]*crLHS139 - crLHS188*crLHS529;
const double crLHS533 = DN(2,0)*N[3];
const double crLHS534 = pow(DN(2,1), 2);
const double crLHS535 = DN(2,1)*crLHS22;
const double crLHS536 = DN(2,2)*crLHS535;
const double crLHS537 = DN(2,1)*N[2];
const double crLHS538 = DN(3,0)*crLHS535;
const double crLHS539 = N[3]*crLHS128;
const double crLHS540 = crLHS130*crLHS457;
const double crLHS541 = crLHS15*crLHS539 - crLHS15*crLHS540;
const double crLHS542 = DN(2,1)*DN(3,1);
const double crLHS543 = N[3]*crLHS153 + crLHS22*crLHS542;
const double crLHS544 = DN(3,2)*crLHS535;
const double crLHS545 = N[3]*crLHS270 - crLHS293*crLHS529;
const double crLHS546 = DN(2,1)*N[3];
const double crLHS547 = pow(DN(2,2), 2);
const double crLHS548 = DN(2,2)*N[2];
const double crLHS549 = DN(2,2)*crLHS22;
const double crLHS550 = DN(3,0)*crLHS549;
const double crLHS551 = crLHS18*crLHS539 - crLHS18*crLHS540;
const double crLHS552 = DN(3,1)*crLHS549;
const double crLHS553 = crLHS19*crLHS539 - crLHS19*crLHS540;
const double crLHS554 = DN(2,2)*DN(3,2);
const double crLHS555 = N[3]*crLHS165 + crLHS22*crLHS554;
const double crLHS556 = DN(2,2)*N[3];
const double crLHS557 = crLHS393*rho;
const double crLHS558 = crLHS395*rho;
const double crLHS559 = crLHS390*rho;
const double crLHS560 = crLHS494*rho;
const double crLHS561 = crLHS496*rho;
const double crLHS562 = crLHS491*rho;
const double crLHS563 = crLHS537*rho;
const double crLHS564 = crLHS548*rho;
const double crLHS565 = crLHS522*rho;
const double crLHS566 = DN(3,0)*N[2];
const double crLHS567 = crLHS546*rho;
const double crLHS568 = crLHS556*rho;
const double crLHS569 = DN(3,1)*N[2];
const double crLHS570 = crLHS533*rho;
const double crLHS571 = DN(3,2)*N[2];
const double crLHS572 = crLHS382*(crLHS525 + crLHS542 + crLHS554);
const double crLHS573 = crLHS180 + crLHS181 + crLHS182;
const double crLHS574 = crLHS573*rho;
const double crLHS575 = crLHS191 + crLHS26*crLHS573;
const double crLHS576 = crLHS55*crLHS573;
const double crLHS577 = N[0]*crLHS576;
const double crLHS578 = crLHS432 + crLHS573*crLHS79;
const double crLHS579 = N[1]*crLHS576;
const double crLHS580 = crLHS128*crLHS573 + crLHS523;
const double crLHS581 = N[2]*crLHS576;
const double crLHS582 = pow(DN(3,0), 2);
const double crLHS583 = pow(N[3], 2);
const double crLHS584 = crLHS583*rho;
const double crLHS585 = crLHS583*crLHS6;
const double crLHS586 = crLHS177*crLHS573 + crLHS585;
const double crLHS587 = DN(3,0)*crLHS22;
const double crLHS588 = DN(3,1)*crLHS587;
const double crLHS589 = crLHS585*rho;
const double crLHS590 = N[3]*crLHS576;
const double crLHS591 = DN(3,2)*crLHS587;
const double crLHS592 = DN(3,0)*N[3];
const double crLHS593 = pow(DN(3,1), 2);
const double crLHS594 = DN(3,1)*DN(3,2)*crLHS22;
const double crLHS595 = DN(3,1)*N[3];
const double crLHS596 = pow(DN(3,2), 2);
const double crLHS597 = DN(3,2)*N[3];
const double crLHS598 = crLHS400*rho;
const double crLHS599 = crLHS402*rho;
const double crLHS600 = crLHS397*rho;
const double crLHS601 = crLHS501*rho;
const double crLHS602 = crLHS503*rho;
const double crLHS603 = crLHS498*rho;
const double crLHS604 = crLHS569*rho;
const double crLHS605 = crLHS571*rho;
const double crLHS606 = crLHS566*rho;
const double crLHS607 = crLHS595*rho;
const double crLHS608 = crLHS597*rho;
const double crLHS609 = crLHS592*rho;
rLHS(0,0)+=gauss_weight*(DN(0,0)*crLHS0 + DN(0,1)*crLHS2 + DN(0,2)*crLHS4 + crLHS12*crLHS24 + crLHS22*crLHS5 + crLHS25*crLHS34 + crLHS34*crLHS36 - crLHS41*crLHS42 + crLHS44);
rLHS(0,1)+=gauss_weight*(DN(0,0)*crLHS45 + DN(0,1)*crLHS47 + DN(0,2)*crLHS50 + crLHS13*crLHS24 - crLHS42*crLHS60 + crLHS52 - crLHS53*crLHS54 - crLHS53*crLHS57);
rLHS(0,2)+=gauss_weight*(DN(0,0)*crLHS61 + DN(0,1)*crLHS63 + DN(0,2)*crLHS65 + crLHS14*crLHS24 - crLHS42*crLHS70 - crLHS54*crLHS67 - crLHS57*crLHS67 + crLHS66);
rLHS(0,3)+=-gauss_weight*(crLHS25*crLHS72 + crLHS36*crLHS72 + crLHS42*crLHS73 + crLHS71);
rLHS(0,4)+=gauss_weight*(DN(0,0)*crLHS74 + DN(0,1)*crLHS76 + DN(0,2)*crLHS78 + crLHS25*crLHS87 + crLHS36*crLHS87 - crLHS42*crLHS92 + crLHS94 + crLHS96);
rLHS(0,5)+=gauss_weight*(DN(0,0)*crLHS97 + DN(0,1)*crLHS99 + DN(0,2)*crLHS102 + crLHS103 - crLHS106*crLHS42 - crLHS107*crLHS53 + crLHS109);
rLHS(0,6)+=gauss_weight*(DN(0,0)*crLHS110 + DN(0,1)*crLHS112 + DN(0,2)*crLHS114 - crLHS107*crLHS67 + crLHS115 - crLHS118*crLHS42 + crLHS119);
rLHS(0,7)+=-gauss_weight*(crLHS120 + crLHS121*crLHS25 + crLHS121*crLHS36 + crLHS122*crLHS42);
rLHS(0,8)+=gauss_weight*(DN(0,0)*crLHS123 + DN(0,1)*crLHS125 + DN(0,2)*crLHS127 + crLHS136*crLHS25 + crLHS136*crLHS36 - crLHS141*crLHS42 + crLHS143 + crLHS145);
rLHS(0,9)+=gauss_weight*(DN(0,0)*crLHS146 + DN(0,1)*crLHS148 + DN(0,2)*crLHS151 + crLHS152 - crLHS155*crLHS42 - crLHS156*crLHS53 + crLHS158);
rLHS(0,10)+=gauss_weight*(DN(0,0)*crLHS159 + DN(0,1)*crLHS161 + DN(0,2)*crLHS163 - crLHS156*crLHS67 + crLHS164 - crLHS167*crLHS42 + crLHS168);
rLHS(0,11)+=-gauss_weight*(crLHS169 + crLHS170*crLHS25 + crLHS170*crLHS36 + crLHS171*crLHS42);
rLHS(0,12)+=gauss_weight*(DN(0,0)*crLHS172 + DN(0,1)*crLHS174 + DN(0,2)*crLHS176 + crLHS185*crLHS25 + crLHS185*crLHS36 - crLHS190*crLHS42 + crLHS192 + crLHS194);
rLHS(0,13)+=gauss_weight*(DN(0,0)*crLHS195 + DN(0,1)*crLHS197 + DN(0,2)*crLHS200 + crLHS201 - crLHS204*crLHS42 - crLHS205*crLHS53 + crLHS207);
rLHS(0,14)+=gauss_weight*(DN(0,0)*crLHS208 + DN(0,1)*crLHS210 + DN(0,2)*crLHS212 - crLHS205*crLHS67 + crLHS213 - crLHS216*crLHS42 + crLHS217);
rLHS(0,15)+=-gauss_weight*(crLHS218 + crLHS219*crLHS25 + crLHS219*crLHS36 + crLHS220*crLHS42);
rLHS(1,0)+=gauss_weight*(DN(0,0)*crLHS2 + DN(0,1)*crLHS221 + DN(0,2)*crLHS222 + crLHS15*crLHS24 - crLHS223*crLHS54 - crLHS223*crLHS57 - crLHS225*crLHS42 + crLHS52);
rLHS(1,1)+=gauss_weight*(DN(0,0)*crLHS47 + DN(0,1)*crLHS226 + DN(0,2)*crLHS228 + crLHS16*crLHS24 + crLHS22*crLHS229 + crLHS230*crLHS25 + crLHS230*crLHS36 - crLHS232*crLHS42 + crLHS44);
rLHS(1,2)+=gauss_weight*(DN(0,0)*crLHS63 + DN(0,1)*crLHS233 + DN(0,2)*crLHS235 + crLHS17*crLHS24 + crLHS237 - crLHS238*crLHS54 - crLHS238*crLHS57 - crLHS239*crLHS42);
rLHS(1,3)+=-gauss_weight*(crLHS240 + crLHS241*crLHS25 + crLHS241*crLHS36 + crLHS242*crLHS42);
rLHS(1,4)+=gauss_weight*(DN(0,0)*crLHS76 + DN(0,1)*crLHS243 + DN(0,2)*crLHS244 - crLHS107*crLHS223 + crLHS245 - crLHS247*crLHS42 + crLHS249);
rLHS(1,5)+=gauss_weight*(DN(0,0)*crLHS99 + DN(0,1)*crLHS250 + DN(0,2)*crLHS252 + crLHS25*crLHS253 + crLHS253*crLHS36 - crLHS255*crLHS42 + crLHS257 + crLHS94);
rLHS(1,6)+=gauss_weight*(DN(0,0)*crLHS112 + DN(0,1)*crLHS258 + DN(0,2)*crLHS260 - crLHS107*crLHS238 + crLHS261 - crLHS262*crLHS42 + crLHS263);
rLHS(1,7)+=-gauss_weight*(crLHS25*crLHS265 + crLHS264 + crLHS265*crLHS36 + crLHS266*crLHS42);
rLHS(1,8)+=gauss_weight*(DN(0,0)*crLHS125 + DN(0,1)*crLHS267 + DN(0,2)*crLHS268 - crLHS156*crLHS223 + crLHS269 - crLHS271*crLHS42 + crLHS272);
rLHS(1,9)+=gauss_weight*(DN(0,0)*crLHS148 + DN(0,1)*crLHS273 + DN(0,2)*crLHS275 + crLHS143 + crLHS25*crLHS276 + crLHS276*crLHS36 - crLHS278*crLHS42 + crLHS280);
rLHS(1,10)+=gauss_weight*(DN(0,0)*crLHS161 + DN(0,1)*crLHS281 + DN(0,2)*crLHS283 - crLHS156*crLHS238 + crLHS284 - crLHS285*crLHS42 + crLHS286);
rLHS(1,11)+=-gauss_weight*(crLHS25*crLHS288 + crLHS287 + crLHS288*crLHS36 + crLHS289*crLHS42);
rLHS(1,12)+=gauss_weight*(DN(0,0)*crLHS174 + DN(0,1)*crLHS290 + DN(0,2)*crLHS291 - crLHS205*crLHS223 + crLHS292 - crLHS294*crLHS42 + crLHS295);
rLHS(1,13)+=gauss_weight*(DN(0,0)*crLHS197 + DN(0,1)*crLHS296 + DN(0,2)*crLHS298 + crLHS192 + crLHS25*crLHS299 + crLHS299*crLHS36 - crLHS301*crLHS42 + crLHS303);
rLHS(1,14)+=gauss_weight*(DN(0,0)*crLHS210 + DN(0,1)*crLHS304 + DN(0,2)*crLHS306 - crLHS205*crLHS238 + crLHS307 - crLHS308*crLHS42 + crLHS309);
rLHS(1,15)+=-gauss_weight*(crLHS25*crLHS311 + crLHS310 + crLHS311*crLHS36 + crLHS312*crLHS42);
rLHS(2,0)+=gauss_weight*(DN(0,0)*crLHS4 + DN(0,1)*crLHS222 + DN(0,2)*crLHS313 + crLHS18*crLHS24 - crLHS314*crLHS54 - crLHS314*crLHS57 - crLHS315*crLHS42 + crLHS66);
rLHS(2,1)+=gauss_weight*(DN(0,0)*crLHS50 + DN(0,1)*crLHS228 + DN(0,2)*crLHS316 + crLHS19*crLHS24 + crLHS237 - crLHS317*crLHS54 - crLHS317*crLHS57 - crLHS318*crLHS42);
rLHS(2,2)+=gauss_weight*(DN(0,0)*crLHS65 + DN(0,1)*crLHS235 + DN(0,2)*crLHS319 + crLHS20*crLHS24 + crLHS22*crLHS320 + crLHS25*crLHS321 + crLHS321*crLHS36 - crLHS322*crLHS42 + crLHS44);
rLHS(2,3)+=-gauss_weight*(crLHS25*crLHS324 + crLHS323 + crLHS324*crLHS36 + crLHS325*crLHS42);
rLHS(2,4)+=gauss_weight*(DN(0,0)*crLHS78 + DN(0,1)*crLHS244 + DN(0,2)*crLHS326 - crLHS107*crLHS314 + crLHS328 - crLHS330*crLHS42 + crLHS332);
rLHS(2,5)+=gauss_weight*(DN(0,0)*crLHS102 + DN(0,1)*crLHS252 + DN(0,2)*crLHS333 - crLHS107*crLHS317 + crLHS334 - crLHS335*crLHS42 + crLHS336);
rLHS(2,6)+=gauss_weight*(DN(0,0)*crLHS114 + DN(0,1)*crLHS260 + DN(0,2)*crLHS337 + crLHS25*crLHS338 + crLHS338*crLHS36 - crLHS339*crLHS42 + crLHS341 + crLHS94);
rLHS(2,7)+=-gauss_weight*(crLHS25*crLHS343 + crLHS342 + crLHS343*crLHS36 + crLHS344*crLHS42);
rLHS(2,8)+=gauss_weight*(DN(0,0)*crLHS127 + DN(0,1)*crLHS268 + DN(0,2)*crLHS345 - crLHS156*crLHS314 + crLHS346 - crLHS347*crLHS42 + crLHS349);
rLHS(2,9)+=gauss_weight*(DN(0,0)*crLHS151 + DN(0,1)*crLHS275 + DN(0,2)*crLHS350 - crLHS156*crLHS317 + crLHS351 - crLHS352*crLHS42 + crLHS353);
rLHS(2,10)+=gauss_weight*(DN(0,0)*crLHS163 + DN(0,1)*crLHS283 + DN(0,2)*crLHS354 + crLHS143 + crLHS25*crLHS355 + crLHS355*crLHS36 - crLHS356*crLHS42 + crLHS358);
rLHS(2,11)+=-gauss_weight*(crLHS25*crLHS360 + crLHS359 + crLHS36*crLHS360 + crLHS361*crLHS42);
rLHS(2,12)+=gauss_weight*(DN(0,0)*crLHS176 + DN(0,1)*crLHS291 + DN(0,2)*crLHS362 - crLHS205*crLHS314 + crLHS363 - crLHS364*crLHS42 + crLHS366);
rLHS(2,13)+=gauss_weight*(DN(0,0)*crLHS200 + DN(0,1)*crLHS298 + DN(0,2)*crLHS367 - crLHS205*crLHS317 + crLHS368 - crLHS369*crLHS42 + crLHS370);
rLHS(2,14)+=gauss_weight*(DN(0,0)*crLHS212 + DN(0,1)*crLHS306 + DN(0,2)*crLHS371 + crLHS192 + crLHS25*crLHS372 + crLHS36*crLHS372 - crLHS373*crLHS42 + crLHS375);
rLHS(2,15)+=-gauss_weight*(crLHS25*crLHS377 + crLHS36*crLHS377 + crLHS376 + crLHS378*crLHS42);
rLHS(3,0)+=gauss_weight*(crLHS223*crLHS379 + crLHS314*crLHS380 - crLHS32*crLHS72 + crLHS71);
rLHS(3,1)+=gauss_weight*(crLHS240 - crLHS241*crLHS59 + crLHS317*crLHS380 + crLHS381*crLHS53);
rLHS(3,2)+=gauss_weight*(crLHS238*crLHS379 + crLHS323 - crLHS324*crLHS69 + crLHS381*crLHS67);
rLHS(3,3)+=crLHS382*(crLHS229 + crLHS320 + crLHS5);
rLHS(3,4)+=gauss_weight*(crLHS223*crLHS384 + crLHS314*crLHS385 + crLHS383 - crLHS72*crLHS86);
rLHS(3,5)+=gauss_weight*(-crLHS105*crLHS241 + crLHS317*crLHS385 + crLHS386 + crLHS387*crLHS53);
rLHS(3,6)+=gauss_weight*(-crLHS117*crLHS324 + crLHS238*crLHS384 + crLHS387*crLHS67 + crLHS388);
rLHS(3,7)+=crLHS389;
rLHS(3,8)+=gauss_weight*(-crLHS135*crLHS72 + crLHS223*crLHS391 + crLHS314*crLHS392 + crLHS390);
rLHS(3,9)+=gauss_weight*(-crLHS154*crLHS241 + crLHS317*crLHS392 + crLHS393 + crLHS394*crLHS53);
rLHS(3,10)+=gauss_weight*(-crLHS166*crLHS324 + crLHS238*crLHS391 + crLHS394*crLHS67 + crLHS395);
rLHS(3,11)+=crLHS396;
rLHS(3,12)+=gauss_weight*(-crLHS184*crLHS72 + crLHS223*crLHS398 + crLHS314*crLHS399 + crLHS397);
rLHS(3,13)+=gauss_weight*(-crLHS203*crLHS241 + crLHS317*crLHS399 + crLHS400 + crLHS401*crLHS53);
rLHS(3,14)+=gauss_weight*(-crLHS215*crLHS324 + crLHS238*crLHS398 + crLHS401*crLHS67 + crLHS402);
rLHS(3,15)+=crLHS403;
rLHS(4,0)+=gauss_weight*(DN(1,0)*crLHS0 + DN(1,1)*crLHS2 + DN(1,2)*crLHS4 + crLHS34*crLHS405 + crLHS34*crLHS81 - crLHS406*crLHS41 + crLHS407 + crLHS96);
rLHS(4,1)+=gauss_weight*(DN(1,0)*crLHS45 + DN(1,1)*crLHS47 + DN(1,2)*crLHS50 + crLHS109 + crLHS245 - crLHS406*crLHS60 - crLHS409*crLHS53);
rLHS(4,2)+=gauss_weight*(DN(1,0)*crLHS61 + DN(1,1)*crLHS63 + DN(1,2)*crLHS65 + crLHS119 + crLHS328 - crLHS406*crLHS70 - crLHS409*crLHS67);
rLHS(4,3)+=-gauss_weight*(crLHS383 + crLHS405*crLHS72 + crLHS406*crLHS73 + crLHS72*crLHS81);
rLHS(4,4)+=gauss_weight*(DN(1,0)*crLHS74 + DN(1,1)*crLHS76 + DN(1,2)*crLHS78 + crLHS12*crLHS412 + crLHS22*crLHS410 + crLHS405*crLHS87 - crLHS406*crLHS92 + crLHS414 + crLHS81*crLHS87);
rLHS(4,5)+=gauss_weight*(DN(1,0)*crLHS97 + DN(1,1)*crLHS99 + DN(1,2)*crLHS102 - crLHS106*crLHS406 + crLHS13*crLHS412 + crLHS416 - crLHS417*crLHS53 - crLHS418*crLHS53);
rLHS(4,6)+=gauss_weight*(DN(1,0)*crLHS110 + DN(1,1)*crLHS112 + DN(1,2)*crLHS114 - crLHS118*crLHS406 + crLHS14*crLHS412 - crLHS417*crLHS67 - crLHS418*crLHS67 + crLHS419);
rLHS(4,7)+=-gauss_weight*(crLHS121*crLHS405 + crLHS121*crLHS81 + crLHS122*crLHS406 + crLHS420);
rLHS(4,8)+=gauss_weight*(DN(1,0)*crLHS123 + DN(1,1)*crLHS125 + DN(1,2)*crLHS127 + crLHS136*crLHS405 + crLHS136*crLHS81 - crLHS141*crLHS406 + crLHS422 + crLHS424);
rLHS(4,9)+=gauss_weight*(DN(1,0)*crLHS146 + DN(1,1)*crLHS148 + DN(1,2)*crLHS151 - crLHS155*crLHS406 + crLHS425 - crLHS426*crLHS53 + crLHS428);
rLHS(4,10)+=gauss_weight*(DN(1,0)*crLHS159 + DN(1,1)*crLHS161 + DN(1,2)*crLHS163 - crLHS167*crLHS406 - crLHS426*crLHS67 + crLHS429 + crLHS430);
rLHS(4,11)+=-gauss_weight*(crLHS170*crLHS405 + crLHS170*crLHS81 + crLHS171*crLHS406 + crLHS431);
rLHS(4,12)+=gauss_weight*(DN(1,0)*crLHS172 + DN(1,1)*crLHS174 + DN(1,2)*crLHS176 + crLHS185*crLHS405 + crLHS185*crLHS81 - crLHS190*crLHS406 + crLHS433 + crLHS435);
rLHS(4,13)+=gauss_weight*(DN(1,0)*crLHS195 + DN(1,1)*crLHS197 + DN(1,2)*crLHS200 - crLHS204*crLHS406 + crLHS436 - crLHS437*crLHS53 + crLHS438);
rLHS(4,14)+=gauss_weight*(DN(1,0)*crLHS208 + DN(1,1)*crLHS210 + DN(1,2)*crLHS212 - crLHS216*crLHS406 - crLHS437*crLHS67 + crLHS439 + crLHS440);
rLHS(4,15)+=-gauss_weight*(crLHS219*crLHS405 + crLHS219*crLHS81 + crLHS220*crLHS406 + crLHS441);
rLHS(5,0)+=gauss_weight*(DN(1,0)*crLHS2 + DN(1,1)*crLHS221 + DN(1,2)*crLHS222 + crLHS103 - crLHS223*crLHS409 - crLHS225*crLHS406 + crLHS249);
rLHS(5,1)+=gauss_weight*(DN(1,0)*crLHS47 + DN(1,1)*crLHS226 + DN(1,2)*crLHS228 + crLHS230*crLHS405 + crLHS230*crLHS81 - crLHS232*crLHS406 + crLHS257 + crLHS407);
rLHS(5,2)+=gauss_weight*(DN(1,0)*crLHS63 + DN(1,1)*crLHS233 + DN(1,2)*crLHS235 - crLHS238*crLHS409 - crLHS239*crLHS406 + crLHS263 + crLHS334);
rLHS(5,3)+=-gauss_weight*(crLHS241*crLHS405 + crLHS241*crLHS81 + crLHS242*crLHS406 + crLHS386);
rLHS(5,4)+=gauss_weight*(DN(1,0)*crLHS76 + DN(1,1)*crLHS243 + DN(1,2)*crLHS244 + crLHS15*crLHS412 - crLHS223*crLHS417 - crLHS223*crLHS418 - crLHS247*crLHS406 + crLHS416);
rLHS(5,5)+=gauss_weight*(DN(1,0)*crLHS99 + DN(1,1)*crLHS250 + DN(1,2)*crLHS252 + crLHS16*crLHS412 + crLHS22*crLHS442 + crLHS253*crLHS405 + crLHS253*crLHS81 - crLHS255*crLHS406 + crLHS414);
rLHS(5,6)+=gauss_weight*(DN(1,0)*crLHS112 + DN(1,1)*crLHS258 + DN(1,2)*crLHS260 + crLHS17*crLHS412 - crLHS238*crLHS417 - crLHS238*crLHS418 - crLHS262*crLHS406 + crLHS444);
rLHS(5,7)+=-gauss_weight*(crLHS265*crLHS405 + crLHS265*crLHS81 + crLHS266*crLHS406 + crLHS445);
rLHS(5,8)+=gauss_weight*(DN(1,0)*crLHS125 + DN(1,1)*crLHS267 + DN(1,2)*crLHS268 - crLHS223*crLHS426 - crLHS271*crLHS406 + crLHS446 + crLHS450);
rLHS(5,9)+=gauss_weight*(DN(1,0)*crLHS148 + DN(1,1)*crLHS273 + DN(1,2)*crLHS275 + crLHS276*crLHS405 + crLHS276*crLHS81 - crLHS278*crLHS406 + crLHS422 + crLHS452);
rLHS(5,10)+=gauss_weight*(DN(1,0)*crLHS161 + DN(1,1)*crLHS281 + DN(1,2)*crLHS283 - crLHS238*crLHS426 - crLHS285*crLHS406 + crLHS453 + crLHS454);
rLHS(5,11)+=-gauss_weight*(crLHS288*crLHS405 + crLHS288*crLHS81 + crLHS289*crLHS406 + crLHS455);
rLHS(5,12)+=gauss_weight*(DN(1,0)*crLHS174 + DN(1,1)*crLHS290 + DN(1,2)*crLHS291 - crLHS223*crLHS437 - crLHS294*crLHS406 + crLHS456 + crLHS458);
rLHS(5,13)+=gauss_weight*(DN(1,0)*crLHS197 + DN(1,1)*crLHS296 + DN(1,2)*crLHS298 + crLHS299*crLHS405 + crLHS299*crLHS81 - crLHS301*crLHS406 + crLHS433 + crLHS460);
rLHS(5,14)+=gauss_weight*(DN(1,0)*crLHS210 + DN(1,1)*crLHS304 + DN(1,2)*crLHS306 - crLHS238*crLHS437 - crLHS308*crLHS406 + crLHS461 + crLHS462);
rLHS(5,15)+=-gauss_weight*(crLHS311*crLHS405 + crLHS311*crLHS81 + crLHS312*crLHS406 + crLHS463);
rLHS(6,0)+=gauss_weight*(DN(1,0)*crLHS4 + DN(1,1)*crLHS222 + DN(1,2)*crLHS313 + crLHS115 - crLHS314*crLHS409 - crLHS315*crLHS406 + crLHS332);
rLHS(6,1)+=gauss_weight*(DN(1,0)*crLHS50 + DN(1,1)*crLHS228 + DN(1,2)*crLHS316 + crLHS261 - crLHS317*crLHS409 - crLHS318*crLHS406 + crLHS336);
rLHS(6,2)+=gauss_weight*(DN(1,0)*crLHS65 + DN(1,1)*crLHS235 + DN(1,2)*crLHS319 + crLHS321*crLHS405 + crLHS321*crLHS81 - crLHS322*crLHS406 + crLHS341 + crLHS407);
rLHS(6,3)+=-gauss_weight*(crLHS324*crLHS405 + crLHS324*crLHS81 + crLHS325*crLHS406 + crLHS388);
rLHS(6,4)+=gauss_weight*(DN(1,0)*crLHS78 + DN(1,1)*crLHS244 + DN(1,2)*crLHS326 + crLHS18*crLHS412 - crLHS314*crLHS417 - crLHS314*crLHS418 - crLHS330*crLHS406 + crLHS419);
rLHS(6,5)+=gauss_weight*(DN(1,0)*crLHS102 + DN(1,1)*crLHS252 + DN(1,2)*crLHS333 + crLHS19*crLHS412 - crLHS317*crLHS417 - crLHS317*crLHS418 - crLHS335*crLHS406 + crLHS444);
rLHS(6,6)+=gauss_weight*(DN(1,0)*crLHS114 + DN(1,1)*crLHS260 + DN(1,2)*crLHS337 + crLHS20*crLHS412 + crLHS22*crLHS464 + crLHS338*crLHS405 + crLHS338*crLHS81 - crLHS339*crLHS406 + crLHS414);
rLHS(6,7)+=-gauss_weight*(crLHS343*crLHS405 + crLHS343*crLHS81 + crLHS344*crLHS406 + crLHS465);
rLHS(6,8)+=gauss_weight*(DN(1,0)*crLHS127 + DN(1,1)*crLHS268 + DN(1,2)*crLHS345 - crLHS314*crLHS426 - crLHS347*crLHS406 + crLHS467 + crLHS470);
rLHS(6,9)+=gauss_weight*(DN(1,0)*crLHS151 + DN(1,1)*crLHS275 + DN(1,2)*crLHS350 - crLHS317*crLHS426 - crLHS352*crLHS406 + crLHS471 + crLHS472);
rLHS(6,10)+=gauss_weight*(DN(1,0)*crLHS163 + DN(1,1)*crLHS283 + DN(1,2)*crLHS354 + crLHS355*crLHS405 + crLHS355*crLHS81 - crLHS356*crLHS406 + crLHS422 + crLHS474);
rLHS(6,11)+=-gauss_weight*(crLHS360*crLHS405 + crLHS360*crLHS81 + crLHS361*crLHS406 + crLHS475);
rLHS(6,12)+=gauss_weight*(DN(1,0)*crLHS176 + DN(1,1)*crLHS291 + DN(1,2)*crLHS362 - crLHS314*crLHS437 - crLHS364*crLHS406 + crLHS476 + crLHS479);
rLHS(6,13)+=gauss_weight*(DN(1,0)*crLHS200 + DN(1,1)*crLHS298 + DN(1,2)*crLHS367 - crLHS317*crLHS437 - crLHS369*crLHS406 + crLHS480 + crLHS481);
rLHS(6,14)+=gauss_weight*(DN(1,0)*crLHS212 + DN(1,1)*crLHS306 + DN(1,2)*crLHS371 + crLHS372*crLHS405 + crLHS372*crLHS81 - crLHS373*crLHS406 + crLHS433 + crLHS483);
rLHS(6,15)+=-gauss_weight*(crLHS377*crLHS405 + crLHS377*crLHS81 + crLHS378*crLHS406 + crLHS484);
rLHS(7,0)+=gauss_weight*(crLHS120 - crLHS121*crLHS32 + crLHS223*crLHS485 + crLHS314*crLHS486);
rLHS(7,1)+=gauss_weight*(crLHS264 - crLHS265*crLHS59 + crLHS317*crLHS486 + crLHS487*crLHS53);
rLHS(7,2)+=gauss_weight*(crLHS238*crLHS485 + crLHS342 - crLHS343*crLHS69 + crLHS487*crLHS67);
rLHS(7,3)+=crLHS389;
rLHS(7,4)+=gauss_weight*(-crLHS121*crLHS86 + crLHS223*crLHS488 + crLHS314*crLHS489 + crLHS420);
rLHS(7,5)+=gauss_weight*(-crLHS105*crLHS265 + crLHS317*crLHS489 + crLHS445 + crLHS490*crLHS53);
rLHS(7,6)+=gauss_weight*(-crLHS117*crLHS343 + crLHS238*crLHS488 + crLHS465 + crLHS490*crLHS67);
rLHS(7,7)+=crLHS382*(crLHS410 + crLHS442 + crLHS464);
rLHS(7,8)+=gauss_weight*(-crLHS121*crLHS135 + crLHS223*crLHS492 + crLHS314*crLHS493 + crLHS491);
rLHS(7,9)+=gauss_weight*(-crLHS154*crLHS265 + crLHS317*crLHS493 + crLHS494 + crLHS495*crLHS53);
rLHS(7,10)+=gauss_weight*(-crLHS166*crLHS343 + crLHS238*crLHS492 + crLHS495*crLHS67 + crLHS496);
rLHS(7,11)+=crLHS497;
rLHS(7,12)+=gauss_weight*(-crLHS121*crLHS184 + crLHS223*crLHS499 + crLHS314*crLHS500 + crLHS498);
rLHS(7,13)+=gauss_weight*(-crLHS203*crLHS265 + crLHS317*crLHS500 + crLHS501 + crLHS502*crLHS53);
rLHS(7,14)+=gauss_weight*(-crLHS215*crLHS343 + crLHS238*crLHS499 + crLHS502*crLHS67 + crLHS503);
rLHS(7,15)+=crLHS504;
rLHS(8,0)+=gauss_weight*(DN(2,0)*crLHS0 + DN(2,1)*crLHS2 + DN(2,2)*crLHS4 + crLHS130*crLHS34 + crLHS145 + crLHS34*crLHS506 - crLHS41*crLHS448 + crLHS507);
rLHS(8,1)+=gauss_weight*(DN(2,0)*crLHS45 + DN(2,1)*crLHS47 + DN(2,2)*crLHS50 + crLHS158 + crLHS269 - crLHS448*crLHS60 - crLHS509*crLHS53);
rLHS(8,2)+=gauss_weight*(DN(2,0)*crLHS61 + DN(2,1)*crLHS63 + DN(2,2)*crLHS65 + crLHS168 + crLHS346 - crLHS448*crLHS70 - crLHS509*crLHS67);
rLHS(8,3)+=-gauss_weight*(crLHS130*crLHS72 + crLHS390 + crLHS448*crLHS73 + crLHS506*crLHS72);
rLHS(8,4)+=gauss_weight*(DN(2,0)*crLHS74 + DN(2,1)*crLHS76 + DN(2,2)*crLHS78 + crLHS130*crLHS87 + crLHS424 - crLHS448*crLHS92 + crLHS506*crLHS87 + crLHS510);
rLHS(8,5)+=gauss_weight*(DN(2,0)*crLHS97 + DN(2,1)*crLHS99 + DN(2,2)*crLHS102 - crLHS106*crLHS448 + crLHS428 + crLHS446 - crLHS511*crLHS53);
rLHS(8,6)+=gauss_weight*(DN(2,0)*crLHS110 + DN(2,1)*crLHS112 + DN(2,2)*crLHS114 - crLHS118*crLHS448 + crLHS430 + crLHS467 - crLHS511*crLHS67);
rLHS(8,7)+=-gauss_weight*(crLHS121*crLHS130 + crLHS121*crLHS506 + crLHS122*crLHS448 + crLHS491);
rLHS(8,8)+=gauss_weight*(DN(2,0)*crLHS123 + DN(2,1)*crLHS125 + DN(2,2)*crLHS127 + crLHS12*crLHS514 + crLHS130*crLHS136 + crLHS136*crLHS506 - crLHS141*crLHS448 + crLHS22*crLHS512 + crLHS516);
rLHS(8,9)+=gauss_weight*(DN(2,0)*crLHS146 + DN(2,1)*crLHS148 + DN(2,2)*crLHS151 + crLHS13*crLHS514 - crLHS155*crLHS448 + crLHS518 - crLHS519*crLHS53 - crLHS520*crLHS53);
rLHS(8,10)+=gauss_weight*(DN(2,0)*crLHS159 + DN(2,1)*crLHS161 + DN(2,2)*crLHS163 + crLHS14*crLHS514 - crLHS167*crLHS448 - crLHS519*crLHS67 - crLHS520*crLHS67 + crLHS521);
rLHS(8,11)+=-gauss_weight*(crLHS130*crLHS170 + crLHS170*crLHS506 + crLHS171*crLHS448 + crLHS522);
rLHS(8,12)+=gauss_weight*(DN(2,0)*crLHS172 + DN(2,1)*crLHS174 + DN(2,2)*crLHS176 + crLHS130*crLHS185 + crLHS185*crLHS506 - crLHS190*crLHS448 + crLHS524 + crLHS526);
rLHS(8,13)+=gauss_weight*(DN(2,0)*crLHS195 + DN(2,1)*crLHS197 + DN(2,2)*crLHS200 - crLHS204*crLHS448 + crLHS527 - crLHS528*crLHS53 + crLHS530);
rLHS(8,14)+=gauss_weight*(DN(2,0)*crLHS208 + DN(2,1)*crLHS210 + DN(2,2)*crLHS212 - crLHS216*crLHS448 - crLHS528*crLHS67 + crLHS531 + crLHS532);
rLHS(8,15)+=-gauss_weight*(crLHS130*crLHS219 + crLHS219*crLHS506 + crLHS220*crLHS448 + crLHS533);
rLHS(9,0)+=gauss_weight*(DN(2,0)*crLHS2 + DN(2,1)*crLHS221 + DN(2,2)*crLHS222 + crLHS152 - crLHS223*crLHS509 - crLHS225*crLHS448 + crLHS272);
rLHS(9,1)+=gauss_weight*(DN(2,0)*crLHS47 + DN(2,1)*crLHS226 + DN(2,2)*crLHS228 + crLHS130*crLHS230 + crLHS230*crLHS506 - crLHS232*crLHS448 + crLHS280 + crLHS507);
rLHS(9,2)+=gauss_weight*(DN(2,0)*crLHS63 + DN(2,1)*crLHS233 + DN(2,2)*crLHS235 - crLHS238*crLHS509 - crLHS239*crLHS448 + crLHS286 + crLHS351);
rLHS(9,3)+=-gauss_weight*(crLHS130*crLHS241 + crLHS241*crLHS506 + crLHS242*crLHS448 + crLHS393);
rLHS(9,4)+=gauss_weight*(DN(2,0)*crLHS76 + DN(2,1)*crLHS243 + DN(2,2)*crLHS244 - crLHS223*crLHS511 - crLHS247*crLHS448 + crLHS425 + crLHS450);
rLHS(9,5)+=gauss_weight*(DN(2,0)*crLHS99 + DN(2,1)*crLHS250 + DN(2,2)*crLHS252 + crLHS130*crLHS253 + crLHS253*crLHS506 - crLHS255*crLHS448 + crLHS452 + crLHS510);
rLHS(9,6)+=gauss_weight*(DN(2,0)*crLHS112 + DN(2,1)*crLHS258 + DN(2,2)*crLHS260 - crLHS238*crLHS511 - crLHS262*crLHS448 + crLHS454 + crLHS471);
rLHS(9,7)+=-gauss_weight*(crLHS130*crLHS265 + crLHS265*crLHS506 + crLHS266*crLHS448 + crLHS494);
rLHS(9,8)+=gauss_weight*(DN(2,0)*crLHS125 + DN(2,1)*crLHS267 + DN(2,2)*crLHS268 + crLHS15*crLHS514 - crLHS223*crLHS519 - crLHS223*crLHS520 - crLHS271*crLHS448 + crLHS518);
rLHS(9,9)+=gauss_weight*(DN(2,0)*crLHS148 + DN(2,1)*crLHS273 + DN(2,2)*crLHS275 + crLHS130*crLHS276 + crLHS16*crLHS514 + crLHS22*crLHS534 + crLHS276*crLHS506 - crLHS278*crLHS448 + crLHS516);
rLHS(9,10)+=gauss_weight*(DN(2,0)*crLHS161 + DN(2,1)*crLHS281 + DN(2,2)*crLHS283 + crLHS17*crLHS514 - crLHS238*crLHS519 - crLHS238*crLHS520 - crLHS285*crLHS448 + crLHS536);
rLHS(9,11)+=-gauss_weight*(crLHS130*crLHS288 + crLHS288*crLHS506 + crLHS289*crLHS448 + crLHS537);
rLHS(9,12)+=gauss_weight*(DN(2,0)*crLHS174 + DN(2,1)*crLHS290 + DN(2,2)*crLHS291 - crLHS223*crLHS528 - crLHS294*crLHS448 + crLHS538 + crLHS541);
rLHS(9,13)+=gauss_weight*(DN(2,0)*crLHS197 + DN(2,1)*crLHS296 + DN(2,2)*crLHS298 + crLHS130*crLHS299 + crLHS299*crLHS506 - crLHS301*crLHS448 + crLHS524 + crLHS543);
rLHS(9,14)+=gauss_weight*(DN(2,0)*crLHS210 + DN(2,1)*crLHS304 + DN(2,2)*crLHS306 - crLHS238*crLHS528 - crLHS308*crLHS448 + crLHS544 + crLHS545);
rLHS(9,15)+=-gauss_weight*(crLHS130*crLHS311 + crLHS311*crLHS506 + crLHS312*crLHS448 + crLHS546);
rLHS(10,0)+=gauss_weight*(DN(2,0)*crLHS4 + DN(2,1)*crLHS222 + DN(2,2)*crLHS313 + crLHS164 - crLHS314*crLHS509 - crLHS315*crLHS448 + crLHS349);
rLHS(10,1)+=gauss_weight*(DN(2,0)*crLHS50 + DN(2,1)*crLHS228 + DN(2,2)*crLHS316 + crLHS284 - crLHS317*crLHS509 - crLHS318*crLHS448 + crLHS353);
rLHS(10,2)+=gauss_weight*(DN(2,0)*crLHS65 + DN(2,1)*crLHS235 + DN(2,2)*crLHS319 + crLHS130*crLHS321 + crLHS321*crLHS506 - crLHS322*crLHS448 + crLHS358 + crLHS507);
rLHS(10,3)+=-gauss_weight*(crLHS130*crLHS324 + crLHS324*crLHS506 + crLHS325*crLHS448 + crLHS395);
rLHS(10,4)+=gauss_weight*(DN(2,0)*crLHS78 + DN(2,1)*crLHS244 + DN(2,2)*crLHS326 - crLHS314*crLHS511 - crLHS330*crLHS448 + crLHS429 + crLHS470);
rLHS(10,5)+=gauss_weight*(DN(2,0)*crLHS102 + DN(2,1)*crLHS252 + DN(2,2)*crLHS333 - crLHS317*crLHS511 - crLHS335*crLHS448 + crLHS453 + crLHS472);
rLHS(10,6)+=gauss_weight*(DN(2,0)*crLHS114 + DN(2,1)*crLHS260 + DN(2,2)*crLHS337 + crLHS130*crLHS338 + crLHS338*crLHS506 - crLHS339*crLHS448 + crLHS474 + crLHS510);
rLHS(10,7)+=-gauss_weight*(crLHS130*crLHS343 + crLHS343*crLHS506 + crLHS344*crLHS448 + crLHS496);
rLHS(10,8)+=gauss_weight*(DN(2,0)*crLHS127 + DN(2,1)*crLHS268 + DN(2,2)*crLHS345 + crLHS18*crLHS514 - crLHS314*crLHS519 - crLHS314*crLHS520 - crLHS347*crLHS448 + crLHS521);
rLHS(10,9)+=gauss_weight*(DN(2,0)*crLHS151 + DN(2,1)*crLHS275 + DN(2,2)*crLHS350 + crLHS19*crLHS514 - crLHS317*crLHS519 - crLHS317*crLHS520 - crLHS352*crLHS448 + crLHS536);
rLHS(10,10)+=gauss_weight*(DN(2,0)*crLHS163 + DN(2,1)*crLHS283 + DN(2,2)*crLHS354 + crLHS130*crLHS355 + crLHS20*crLHS514 + crLHS22*crLHS547 + crLHS355*crLHS506 - crLHS356*crLHS448 + crLHS516);
rLHS(10,11)+=-gauss_weight*(crLHS130*crLHS360 + crLHS360*crLHS506 + crLHS361*crLHS448 + crLHS548);
rLHS(10,12)+=gauss_weight*(DN(2,0)*crLHS176 + DN(2,1)*crLHS291 + DN(2,2)*crLHS362 - crLHS314*crLHS528 - crLHS364*crLHS448 + crLHS550 + crLHS551);
rLHS(10,13)+=gauss_weight*(DN(2,0)*crLHS200 + DN(2,1)*crLHS298 + DN(2,2)*crLHS367 - crLHS317*crLHS528 - crLHS369*crLHS448 + crLHS552 + crLHS553);
rLHS(10,14)+=gauss_weight*(DN(2,0)*crLHS212 + DN(2,1)*crLHS306 + DN(2,2)*crLHS371 + crLHS130*crLHS372 + crLHS372*crLHS506 - crLHS373*crLHS448 + crLHS524 + crLHS555);
rLHS(10,15)+=-gauss_weight*(crLHS130*crLHS377 + crLHS377*crLHS506 + crLHS378*crLHS448 + crLHS556);
rLHS(11,0)+=gauss_weight*(crLHS169 - crLHS170*crLHS32 + crLHS223*crLHS557 + crLHS314*crLHS558);
rLHS(11,1)+=gauss_weight*(crLHS287 - crLHS288*crLHS59 + crLHS317*crLHS558 + crLHS53*crLHS559);
rLHS(11,2)+=gauss_weight*(crLHS238*crLHS557 + crLHS359 - crLHS360*crLHS69 + crLHS559*crLHS67);
rLHS(11,3)+=crLHS396;
rLHS(11,4)+=gauss_weight*(-crLHS170*crLHS86 + crLHS223*crLHS560 + crLHS314*crLHS561 + crLHS431);
rLHS(11,5)+=gauss_weight*(-crLHS105*crLHS288 + crLHS317*crLHS561 + crLHS455 + crLHS53*crLHS562);
rLHS(11,6)+=gauss_weight*(-crLHS117*crLHS360 + crLHS238*crLHS560 + crLHS475 + crLHS562*crLHS67);
rLHS(11,7)+=crLHS497;
rLHS(11,8)+=gauss_weight*(-crLHS135*crLHS170 + crLHS223*crLHS563 + crLHS314*crLHS564 + crLHS522);
rLHS(11,9)+=gauss_weight*(-crLHS154*crLHS288 + crLHS317*crLHS564 + crLHS53*crLHS565 + crLHS537);
rLHS(11,10)+=gauss_weight*(-crLHS166*crLHS360 + crLHS238*crLHS563 + crLHS548 + crLHS565*crLHS67);
rLHS(11,11)+=crLHS382*(crLHS512 + crLHS534 + crLHS547);
rLHS(11,12)+=gauss_weight*(-crLHS170*crLHS184 + crLHS223*crLHS567 + crLHS314*crLHS568 + crLHS566);
rLHS(11,13)+=gauss_weight*(-crLHS203*crLHS288 + crLHS317*crLHS568 + crLHS53*crLHS570 + crLHS569);
rLHS(11,14)+=gauss_weight*(-crLHS215*crLHS360 + crLHS238*crLHS567 + crLHS570*crLHS67 + crLHS571);
rLHS(11,15)+=crLHS572;
rLHS(12,0)+=gauss_weight*(DN(3,0)*crLHS0 + DN(3,1)*crLHS2 + DN(3,2)*crLHS4 + crLHS179*crLHS34 + crLHS194 + crLHS34*crLHS574 - crLHS41*crLHS457 + crLHS575);
rLHS(12,1)+=gauss_weight*(DN(3,0)*crLHS45 + DN(3,1)*crLHS47 + DN(3,2)*crLHS50 + crLHS207 + crLHS292 - crLHS457*crLHS60 - crLHS53*crLHS577);
rLHS(12,2)+=gauss_weight*(DN(3,0)*crLHS61 + DN(3,1)*crLHS63 + DN(3,2)*crLHS65 + crLHS217 + crLHS363 - crLHS457*crLHS70 - crLHS577*crLHS67);
rLHS(12,3)+=-gauss_weight*(crLHS179*crLHS72 + crLHS397 + crLHS457*crLHS73 + crLHS574*crLHS72);
rLHS(12,4)+=gauss_weight*(DN(3,0)*crLHS74 + DN(3,1)*crLHS76 + DN(3,2)*crLHS78 + crLHS179*crLHS87 + crLHS435 - crLHS457*crLHS92 + crLHS574*crLHS87 + crLHS578);
rLHS(12,5)+=gauss_weight*(DN(3,0)*crLHS97 + DN(3,1)*crLHS99 + DN(3,2)*crLHS102 - crLHS106*crLHS457 + crLHS438 + crLHS456 - crLHS53*crLHS579);
rLHS(12,6)+=gauss_weight*(DN(3,0)*crLHS110 + DN(3,1)*crLHS112 + DN(3,2)*crLHS114 - crLHS118*crLHS457 + crLHS440 + crLHS476 - crLHS579*crLHS67);
rLHS(12,7)+=-gauss_weight*(crLHS121*crLHS179 + crLHS121*crLHS574 + crLHS122*crLHS457 + crLHS498);
rLHS(12,8)+=gauss_weight*(DN(3,0)*crLHS123 + DN(3,1)*crLHS125 + DN(3,2)*crLHS127 + crLHS136*crLHS179 + crLHS136*crLHS574 - crLHS141*crLHS457 + crLHS526 + crLHS580);
rLHS(12,9)+=gauss_weight*(DN(3,0)*crLHS146 + DN(3,1)*crLHS148 + DN(3,2)*crLHS151 - crLHS155*crLHS457 - crLHS53*crLHS581 + crLHS530 + crLHS538);
rLHS(12,10)+=gauss_weight*(DN(3,0)*crLHS159 + DN(3,1)*crLHS161 + DN(3,2)*crLHS163 - crLHS167*crLHS457 + crLHS532 + crLHS550 - crLHS581*crLHS67);
rLHS(12,11)+=-gauss_weight*(crLHS170*crLHS179 + crLHS170*crLHS574 + crLHS171*crLHS457 + crLHS566);
rLHS(12,12)+=gauss_weight*(DN(3,0)*crLHS172 + DN(3,1)*crLHS174 + DN(3,2)*crLHS176 + crLHS12*crLHS584 + crLHS179*crLHS185 + crLHS185*crLHS574 - crLHS190*crLHS457 + crLHS22*crLHS582 + crLHS586);
rLHS(12,13)+=gauss_weight*(DN(3,0)*crLHS195 + DN(3,1)*crLHS197 + DN(3,2)*crLHS200 + crLHS13*crLHS584 - crLHS204*crLHS457 - crLHS53*crLHS589 - crLHS53*crLHS590 + crLHS588);
rLHS(12,14)+=gauss_weight*(DN(3,0)*crLHS208 + DN(3,1)*crLHS210 + DN(3,2)*crLHS212 + crLHS14*crLHS584 - crLHS216*crLHS457 - crLHS589*crLHS67 - crLHS590*crLHS67 + crLHS591);
rLHS(12,15)+=-gauss_weight*(crLHS179*crLHS219 + crLHS219*crLHS574 + crLHS220*crLHS457 + crLHS592);
rLHS(13,0)+=gauss_weight*(DN(3,0)*crLHS2 + DN(3,1)*crLHS221 + DN(3,2)*crLHS222 + crLHS201 - crLHS223*crLHS577 - crLHS225*crLHS457 + crLHS295);
rLHS(13,1)+=gauss_weight*(DN(3,0)*crLHS47 + DN(3,1)*crLHS226 + DN(3,2)*crLHS228 + crLHS179*crLHS230 + crLHS230*crLHS574 - crLHS232*crLHS457 + crLHS303 + crLHS575);
rLHS(13,2)+=gauss_weight*(DN(3,0)*crLHS63 + DN(3,1)*crLHS233 + DN(3,2)*crLHS235 - crLHS238*crLHS577 - crLHS239*crLHS457 + crLHS309 + crLHS368);
rLHS(13,3)+=-gauss_weight*(crLHS179*crLHS241 + crLHS241*crLHS574 + crLHS242*crLHS457 + crLHS400);
rLHS(13,4)+=gauss_weight*(DN(3,0)*crLHS76 + DN(3,1)*crLHS243 + DN(3,2)*crLHS244 - crLHS223*crLHS579 - crLHS247*crLHS457 + crLHS436 + crLHS458);
rLHS(13,5)+=gauss_weight*(DN(3,0)*crLHS99 + DN(3,1)*crLHS250 + DN(3,2)*crLHS252 + crLHS179*crLHS253 + crLHS253*crLHS574 - crLHS255*crLHS457 + crLHS460 + crLHS578);
rLHS(13,6)+=gauss_weight*(DN(3,0)*crLHS112 + DN(3,1)*crLHS258 + DN(3,2)*crLHS260 - crLHS238*crLHS579 - crLHS262*crLHS457 + crLHS462 + crLHS480);
rLHS(13,7)+=-gauss_weight*(crLHS179*crLHS265 + crLHS265*crLHS574 + crLHS266*crLHS457 + crLHS501);
rLHS(13,8)+=gauss_weight*(DN(3,0)*crLHS125 + DN(3,1)*crLHS267 + DN(3,2)*crLHS268 - crLHS223*crLHS581 - crLHS271*crLHS457 + crLHS527 + crLHS541);
rLHS(13,9)+=gauss_weight*(DN(3,0)*crLHS148 + DN(3,1)*crLHS273 + DN(3,2)*crLHS275 + crLHS179*crLHS276 + crLHS276*crLHS574 - crLHS278*crLHS457 + crLHS543 + crLHS580);
rLHS(13,10)+=gauss_weight*(DN(3,0)*crLHS161 + DN(3,1)*crLHS281 + DN(3,2)*crLHS283 - crLHS238*crLHS581 - crLHS285*crLHS457 + crLHS545 + crLHS552);
rLHS(13,11)+=-gauss_weight*(crLHS179*crLHS288 + crLHS288*crLHS574 + crLHS289*crLHS457 + crLHS569);
rLHS(13,12)+=gauss_weight*(DN(3,0)*crLHS174 + DN(3,1)*crLHS290 + DN(3,2)*crLHS291 + crLHS15*crLHS584 - crLHS223*crLHS589 - crLHS223*crLHS590 - crLHS294*crLHS457 + crLHS588);
rLHS(13,13)+=gauss_weight*(DN(3,0)*crLHS197 + DN(3,1)*crLHS296 + DN(3,2)*crLHS298 + crLHS16*crLHS584 + crLHS179*crLHS299 + crLHS22*crLHS593 + crLHS299*crLHS574 - crLHS301*crLHS457 + crLHS586);
rLHS(13,14)+=gauss_weight*(DN(3,0)*crLHS210 + DN(3,1)*crLHS304 + DN(3,2)*crLHS306 + crLHS17*crLHS584 - crLHS238*crLHS589 - crLHS238*crLHS590 - crLHS308*crLHS457 + crLHS594);
rLHS(13,15)+=-gauss_weight*(crLHS179*crLHS311 + crLHS311*crLHS574 + crLHS312*crLHS457 + crLHS595);
rLHS(14,0)+=gauss_weight*(DN(3,0)*crLHS4 + DN(3,1)*crLHS222 + DN(3,2)*crLHS313 + crLHS213 - crLHS314*crLHS577 - crLHS315*crLHS457 + crLHS366);
rLHS(14,1)+=gauss_weight*(DN(3,0)*crLHS50 + DN(3,1)*crLHS228 + DN(3,2)*crLHS316 + crLHS307 - crLHS317*crLHS577 - crLHS318*crLHS457 + crLHS370);
rLHS(14,2)+=gauss_weight*(DN(3,0)*crLHS65 + DN(3,1)*crLHS235 + DN(3,2)*crLHS319 + crLHS179*crLHS321 + crLHS321*crLHS574 - crLHS322*crLHS457 + crLHS375 + crLHS575);
rLHS(14,3)+=-gauss_weight*(crLHS179*crLHS324 + crLHS324*crLHS574 + crLHS325*crLHS457 + crLHS402);
rLHS(14,4)+=gauss_weight*(DN(3,0)*crLHS78 + DN(3,1)*crLHS244 + DN(3,2)*crLHS326 - crLHS314*crLHS579 - crLHS330*crLHS457 + crLHS439 + crLHS479);
rLHS(14,5)+=gauss_weight*(DN(3,0)*crLHS102 + DN(3,1)*crLHS252 + DN(3,2)*crLHS333 - crLHS317*crLHS579 - crLHS335*crLHS457 + crLHS461 + crLHS481);
rLHS(14,6)+=gauss_weight*(DN(3,0)*crLHS114 + DN(3,1)*crLHS260 + DN(3,2)*crLHS337 + crLHS179*crLHS338 + crLHS338*crLHS574 - crLHS339*crLHS457 + crLHS483 + crLHS578);
rLHS(14,7)+=-gauss_weight*(crLHS179*crLHS343 + crLHS343*crLHS574 + crLHS344*crLHS457 + crLHS503);
rLHS(14,8)+=gauss_weight*(DN(3,0)*crLHS127 + DN(3,1)*crLHS268 + DN(3,2)*crLHS345 - crLHS314*crLHS581 - crLHS347*crLHS457 + crLHS531 + crLHS551);
rLHS(14,9)+=gauss_weight*(DN(3,0)*crLHS151 + DN(3,1)*crLHS275 + DN(3,2)*crLHS350 - crLHS317*crLHS581 - crLHS352*crLHS457 + crLHS544 + crLHS553);
rLHS(14,10)+=gauss_weight*(DN(3,0)*crLHS163 + DN(3,1)*crLHS283 + DN(3,2)*crLHS354 + crLHS179*crLHS355 + crLHS355*crLHS574 - crLHS356*crLHS457 + crLHS555 + crLHS580);
rLHS(14,11)+=-gauss_weight*(crLHS179*crLHS360 + crLHS360*crLHS574 + crLHS361*crLHS457 + crLHS571);
rLHS(14,12)+=gauss_weight*(DN(3,0)*crLHS176 + DN(3,1)*crLHS291 + DN(3,2)*crLHS362 + crLHS18*crLHS584 - crLHS314*crLHS589 - crLHS314*crLHS590 - crLHS364*crLHS457 + crLHS591);
rLHS(14,13)+=gauss_weight*(DN(3,0)*crLHS200 + DN(3,1)*crLHS298 + DN(3,2)*crLHS367 + crLHS19*crLHS584 - crLHS317*crLHS589 - crLHS317*crLHS590 - crLHS369*crLHS457 + crLHS594);
rLHS(14,14)+=gauss_weight*(DN(3,0)*crLHS212 + DN(3,1)*crLHS306 + DN(3,2)*crLHS371 + crLHS179*crLHS372 + crLHS20*crLHS584 + crLHS22*crLHS596 + crLHS372*crLHS574 - crLHS373*crLHS457 + crLHS586);
rLHS(14,15)+=-gauss_weight*(crLHS179*crLHS377 + crLHS377*crLHS574 + crLHS378*crLHS457 + crLHS597);
rLHS(15,0)+=gauss_weight*(crLHS218 - crLHS219*crLHS32 + crLHS223*crLHS598 + crLHS314*crLHS599);
rLHS(15,1)+=gauss_weight*(crLHS310 - crLHS311*crLHS59 + crLHS317*crLHS599 + crLHS53*crLHS600);
rLHS(15,2)+=gauss_weight*(crLHS238*crLHS598 + crLHS376 - crLHS377*crLHS69 + crLHS600*crLHS67);
rLHS(15,3)+=crLHS403;
rLHS(15,4)+=gauss_weight*(-crLHS219*crLHS86 + crLHS223*crLHS601 + crLHS314*crLHS602 + crLHS441);
rLHS(15,5)+=gauss_weight*(-crLHS105*crLHS311 + crLHS317*crLHS602 + crLHS463 + crLHS53*crLHS603);
rLHS(15,6)+=gauss_weight*(-crLHS117*crLHS377 + crLHS238*crLHS601 + crLHS484 + crLHS603*crLHS67);
rLHS(15,7)+=crLHS504;
rLHS(15,8)+=gauss_weight*(-crLHS135*crLHS219 + crLHS223*crLHS604 + crLHS314*crLHS605 + crLHS533);
rLHS(15,9)+=gauss_weight*(-crLHS154*crLHS311 + crLHS317*crLHS605 + crLHS53*crLHS606 + crLHS546);
rLHS(15,10)+=gauss_weight*(-crLHS166*crLHS377 + crLHS238*crLHS604 + crLHS556 + crLHS606*crLHS67);
rLHS(15,11)+=crLHS572;
rLHS(15,12)+=gauss_weight*(-crLHS184*crLHS219 + crLHS223*crLHS607 + crLHS314*crLHS608 + crLHS592);
rLHS(15,13)+=gauss_weight*(-crLHS203*crLHS311 + crLHS317*crLHS608 + crLHS53*crLHS609 + crLHS595);
rLHS(15,14)+=gauss_weight*(-crLHS215*crLHS377 + crLHS238*crLHS607 + crLHS597 + crLHS609*crLHS67);
rLHS(15,15)+=crLHS382*(crLHS582 + crLHS593 + crLHS596);

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
const double crRHS33 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1);
const double crRHS34 = rho*stab_c3*sqrt(pow(crRHS11, 2) + pow(crRHS18, 2) + pow(crRHS20, 2) + pow(crRHS33, 2));
const double crRHS35 = crRHS30*(h*(crRHS31*h + crRHS32 + crRHS34*h)/stab_c1 + mu);
const double crRHS36 = crRHS5*rho;
const double crRHS37 = crRHS16*rho;
const double crRHS38 = crRHS14*functional_weights[6];
const double crRHS39 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] - crRHS1 + crRHS14 + crRHS23*rho + crRHS25*rho - crRHS28*crRHS36 - crRHS37*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0)) + crRHS38 + crRHS4 + crRHS7;
const double crRHS40 = 1.0/(crRHS31 + crRHS32/h + crRHS34 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS41 = crRHS39*crRHS40;
const double crRHS42 = N[0]*crRHS2;
const double crRHS43 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1));
const double crRHS44 = crRHS2*crRHS24;
const double crRHS45 = crRHS20*crRHS3;
const double crRHS46 = crRHS24*crRHS33;
const double crRHS47 = crRHS16*crRHS6;
const double crRHS48 = crRHS13*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2]);
const double crRHS49 = crRHS48*functional_weights[6];
const double crRHS50 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] - crRHS29*crRHS37 - crRHS36*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1)) - crRHS43 + crRHS44 + crRHS45*rho + crRHS46*rho + crRHS47 + crRHS48 + crRHS49;
const double crRHS51 = crRHS11*crRHS50 + crRHS18*crRHS39;
const double crRHS52 = crRHS27*crRHS40;
const double crRHS53 = 1.0*crRHS33;
const double crRHS54 = crRHS45 + crRHS46;
const double crRHS55 = crRHS40*crRHS50;
const double crRHS56 = crRHS20*crRHS39 + crRHS33*crRHS50;
const double crRHS57 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS16);
const double crRHS58 = N[1]*rho;
const double crRHS59 = N[1]*crRHS2;
const double crRHS60 = crRHS40*crRHS58;
const double crRHS61 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS16);
const double crRHS62 = N[2]*rho;
const double crRHS63 = N[2]*crRHS2;
const double crRHS64 = crRHS40*crRHS62;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS35 + DN(0,0)*stress_adj[0] - DN(0,1)*crRHS12 + DN(0,1)*stress_adj[2] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS15*functional_weights[6] + crRHS15 + crRHS17*crRHS3 - crRHS17*crRHS41 + crRHS22*(DN(0,0)*crRHS19 + DN(0,1)*crRHS21) + crRHS26*crRHS27 - crRHS41*crRHS42 - crRHS51*crRHS52);
rRHS[1]+=-gauss_weight*(DN(0,0)*crRHS12 + DN(0,0)*stress_adj[2] - DN(0,1)*crRHS0 + DN(0,1)*crRHS35 + DN(0,1)*stress_adj[1] - N[0]*crRHS43 + N[0]*crRHS44 + N[0]*crRHS47 + N[0]*crRHS48 + N[0]*crRHS49 + crRHS17*crRHS24 - crRHS17*crRHS55 + crRHS22*(DN(0,0)*crRHS21 + DN(0,1)*crRHS53) + crRHS27*crRHS54 - crRHS42*crRHS55 - crRHS52*crRHS56);
rRHS[2]+=-gauss_weight*(DN(0,0)*crRHS41 + DN(0,1)*crRHS55 + N[0]*crRHS30);
rRHS[3]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS35 + DN(1,0)*stress_adj[0] - DN(1,1)*crRHS12 + DN(1,1)*stress_adj[2] - N[1]*crRHS1 + N[1]*crRHS14 + N[1]*crRHS38 + N[1]*crRHS4 + N[1]*crRHS7 + crRHS22*(DN(1,0)*crRHS19 + DN(1,1)*crRHS21) + crRHS26*crRHS58 + crRHS3*crRHS57 - crRHS41*crRHS57 - crRHS41*crRHS59 - crRHS51*crRHS60);
rRHS[4]+=-gauss_weight*(DN(1,0)*crRHS12 + DN(1,0)*stress_adj[2] - DN(1,1)*crRHS0 + DN(1,1)*crRHS35 + DN(1,1)*stress_adj[1] - N[1]*crRHS43 + N[1]*crRHS44 + N[1]*crRHS47 + N[1]*crRHS48 + N[1]*crRHS49 + crRHS22*(DN(1,0)*crRHS21 + DN(1,1)*crRHS53) + crRHS24*crRHS57 + crRHS54*crRHS58 - crRHS55*crRHS57 - crRHS55*crRHS59 - crRHS56*crRHS60);
rRHS[5]+=-gauss_weight*(DN(1,0)*crRHS41 + DN(1,1)*crRHS55 + N[1]*crRHS30);
rRHS[6]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS35 + DN(2,0)*stress_adj[0] - DN(2,1)*crRHS12 + DN(2,1)*stress_adj[2] - N[2]*crRHS1 + N[2]*crRHS14 + N[2]*crRHS38 + N[2]*crRHS4 + N[2]*crRHS7 + crRHS22*(DN(2,0)*crRHS19 + DN(2,1)*crRHS21) + crRHS26*crRHS62 + crRHS3*crRHS61 - crRHS41*crRHS61 - crRHS41*crRHS63 - crRHS51*crRHS64);
rRHS[7]+=-gauss_weight*(DN(2,0)*crRHS12 + DN(2,0)*stress_adj[2] - DN(2,1)*crRHS0 + DN(2,1)*crRHS35 + DN(2,1)*stress_adj[1] - N[2]*crRHS43 + N[2]*crRHS44 + N[2]*crRHS47 + N[2]*crRHS48 + N[2]*crRHS49 + crRHS22*(DN(2,0)*crRHS21 + DN(2,1)*crRHS53) + crRHS24*crRHS61 + crRHS54*crRHS62 - crRHS55*crRHS61 - crRHS55*crRHS63 - crRHS56*crRHS64);
rRHS[8]+=-gauss_weight*(DN(2,0)*crRHS41 + DN(2,1)*crRHS55 + N[2]*crRHS30);

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
const double crRHS48 = DN(0,1)*v_ns(0,1) + DN(1,1)*v_ns(1,1) + DN(2,1)*v_ns(2,1) + DN(3,1)*v_ns(3,1);
const double crRHS49 = DN(0,1)*v_ns(0,2) + DN(1,1)*v_ns(1,2) + DN(2,1)*v_ns(2,2) + DN(3,1)*v_ns(3,2);
const double crRHS50 = DN(0,2)*v_ns(0,1);
const double crRHS51 = DN(1,2)*v_ns(1,1);
const double crRHS52 = DN(2,2)*v_ns(2,1);
const double crRHS53 = DN(3,2)*v_ns(3,1);
const double crRHS54 = crRHS50 + crRHS51 + crRHS52 + crRHS53;
const double crRHS55 = DN(0,2)*v_ns(0,2) + DN(1,2)*v_ns(1,2) + DN(2,2)*v_ns(2,2) + DN(3,2)*v_ns(3,2);
const double crRHS56 = rho*stab_c3*sqrt(pow(crRHS18, 2) + pow(crRHS24, 2) + pow(crRHS27, 2) + pow(crRHS29, 2) + pow(crRHS31, 2) + pow(crRHS48, 2) + pow(crRHS49, 2) + pow(crRHS54, 2) + pow(crRHS55, 2));
const double crRHS57 = crRHS45*(h*(crRHS46*h + crRHS47 + crRHS56*h)/stab_c1 + mu);
const double crRHS58 = crRHS5*rho;
const double crRHS59 = crRHS11*rho;
const double crRHS60 = crRHS12*rho;
const double crRHS61 = crRHS9*functional_weights[6];
const double crRHS62 = DN(0,0)*p_adj[0] + DN(1,0)*p_adj[1] + DN(2,0)*p_adj[2] + DN(3,0)*p_adj[3] - crRHS1 + crRHS35*rho + crRHS37*rho + crRHS39*rho + crRHS4 - crRHS42*crRHS58 - crRHS59*(DN(0,1)*v_adj(0,0) + DN(1,1)*v_adj(1,0) + DN(2,1)*v_adj(2,0) + DN(3,1)*v_adj(3,0)) - crRHS60*(DN(0,2)*v_adj(0,0) + DN(1,2)*v_adj(1,0) + DN(2,2)*v_adj(2,0) + DN(3,2)*v_adj(3,0)) + crRHS61 + crRHS7 + crRHS9;
const double crRHS63 = 1.0/(crRHS46 + crRHS47/h + crRHS56 + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crRHS64 = crRHS62*crRHS63;
const double crRHS65 = N[0]*crRHS2;
const double crRHS66 = rho*(N[0]*f_adj(0,1) + N[1]*f_adj(1,1) + N[2]*f_adj(2,1) + N[3]*f_adj(3,1));
const double crRHS67 = crRHS2*crRHS36;
const double crRHS68 = crRHS29*crRHS3;
const double crRHS69 = crRHS36*crRHS48;
const double crRHS70 = crRHS38*crRHS49;
const double crRHS71 = crRHS11*crRHS6;
const double crRHS72 = crRHS8*(DN(0,1)*t[0] + DN(1,1)*t[1] + DN(2,1)*t[2] + DN(3,1)*t[3]);
const double crRHS73 = crRHS72*functional_weights[6];
const double crRHS74 = DN(0,1)*p_adj[0] + DN(1,1)*p_adj[1] + DN(2,1)*p_adj[2] + DN(3,1)*p_adj[3] - crRHS43*crRHS59 - crRHS58*(DN(0,0)*v_adj(0,1) + DN(1,0)*v_adj(1,1) + DN(2,0)*v_adj(2,1) + DN(3,0)*v_adj(3,1)) - crRHS60*(DN(0,2)*v_adj(0,1) + DN(1,2)*v_adj(1,1) + DN(2,2)*v_adj(2,1) + DN(3,2)*v_adj(3,1)) - crRHS66 + crRHS67 + crRHS68*rho + crRHS69*rho + crRHS70*rho + crRHS71 + crRHS72 + crRHS73;
const double crRHS75 = rho*(N[0]*f_adj(0,2) + N[1]*f_adj(1,2) + N[2]*f_adj(2,2) + N[3]*f_adj(3,2));
const double crRHS76 = crRHS2*crRHS38;
const double crRHS77 = crRHS3*crRHS31;
const double crRHS78 = crRHS36*crRHS54;
const double crRHS79 = crRHS38*crRHS55;
const double crRHS80 = crRHS12*crRHS6;
const double crRHS81 = crRHS8*(DN(0,2)*t[0] + DN(1,2)*t[1] + DN(2,2)*t[2] + DN(3,2)*t[3]);
const double crRHS82 = crRHS81*functional_weights[6];
const double crRHS83 = DN(0,2)*p_adj[0] + DN(1,2)*p_adj[1] + DN(2,2)*p_adj[2] + DN(3,2)*p_adj[3] - crRHS44*crRHS60 - crRHS58*(DN(0,0)*v_adj(0,2) + DN(1,0)*v_adj(1,2) + DN(2,0)*v_adj(2,2) + DN(3,0)*v_adj(3,2)) - crRHS59*(DN(0,1)*v_adj(0,2) + DN(1,1)*v_adj(1,2) + DN(2,1)*v_adj(2,2) + DN(3,1)*v_adj(3,2)) - crRHS75 + crRHS76 + crRHS77*rho + crRHS78*rho + crRHS79*rho + crRHS80 + crRHS81 + crRHS82;
const double crRHS84 = crRHS18*crRHS74 + crRHS24*crRHS83 + crRHS27*crRHS62;
const double crRHS85 = crRHS41*crRHS63;
const double crRHS86 = crRHS49 - crRHS50 - crRHS51 - crRHS52 - crRHS53;
const double crRHS87 = 1.0*crRHS48;
const double crRHS88 = crRHS49 + crRHS54;
const double crRHS89 = crRHS68 + crRHS69 + crRHS70;
const double crRHS90 = crRHS63*crRHS74;
const double crRHS91 = crRHS29*crRHS62 + crRHS48*crRHS74 + crRHS49*crRHS83;
const double crRHS92 = 1.0*crRHS55;
const double crRHS93 = 0.5*crRHS32;
const double crRHS94 = 0.5*crRHS88;
const double crRHS95 = crRHS77 + crRHS78 + crRHS79;
const double crRHS96 = crRHS63*crRHS83;
const double crRHS97 = crRHS31*crRHS62 + crRHS54*crRHS74 + crRHS55*crRHS83;
const double crRHS98 = rho*(DN(1,0)*crRHS5 + DN(1,1)*crRHS11 + DN(1,2)*crRHS12);
const double crRHS99 = N[1]*rho;
const double crRHS100 = N[1]*crRHS2;
const double crRHS101 = crRHS63*crRHS99;
const double crRHS102 = rho*(DN(2,0)*crRHS5 + DN(2,1)*crRHS11 + DN(2,2)*crRHS12);
const double crRHS103 = N[2]*rho;
const double crRHS104 = N[2]*crRHS2;
const double crRHS105 = crRHS103*crRHS63;
const double crRHS106 = rho*(DN(3,0)*crRHS5 + DN(3,1)*crRHS11 + DN(3,2)*crRHS12);
const double crRHS107 = N[3]*rho;
const double crRHS108 = N[3]*crRHS2;
const double crRHS109 = crRHS107*crRHS63;
rRHS[0]+=-gauss_weight*(-DN(0,0)*crRHS0 + DN(0,0)*crRHS57 + DN(0,0)*stress_adj[0] + DN(0,1)*stress_adj[3] + DN(0,2)*stress_adj[5] - N[0]*crRHS1 + N[0]*crRHS4 + N[0]*crRHS7 + crRHS10*functional_weights[6] + crRHS10 + crRHS13*crRHS3 - crRHS13*crRHS64 - crRHS26*(DN(0,1)*crRHS19 + DN(0,2)*crRHS25) + crRHS34*(DN(0,0)*crRHS28 + DN(0,1)*crRHS30 + crRHS32*crRHS33) + crRHS40*crRHS41 - crRHS64*crRHS65 - crRHS84*crRHS85);
rRHS[1]+=-gauss_weight*(DN(0,0)*stress_adj[3] - DN(0,1)*crRHS0 + DN(0,1)*crRHS57 + DN(0,1)*stress_adj[1] + DN(0,2)*stress_adj[4] - N[0]*crRHS66 + N[0]*crRHS67 + N[0]*crRHS71 + N[0]*crRHS72 + N[0]*crRHS73 + crRHS13*crRHS36 - crRHS13*crRHS90 + crRHS26*(DN(0,0)*crRHS19 - DN(0,2)*crRHS86) + crRHS34*(DN(0,0)*crRHS30 + DN(0,1)*crRHS87 + crRHS33*crRHS88) + crRHS41*crRHS89 - crRHS65*crRHS90 - crRHS85*crRHS91);
rRHS[2]+=-gauss_weight*(DN(0,0)*stress_adj[5] + DN(0,1)*stress_adj[4] - DN(0,2)*crRHS0 + DN(0,2)*crRHS57 + DN(0,2)*stress_adj[2] - N[0]*crRHS75 + N[0]*crRHS76 + N[0]*crRHS80 + N[0]*crRHS81 + N[0]*crRHS82 + crRHS13*crRHS38 - crRHS13*crRHS96 + crRHS26*(DN(0,0)*crRHS25 + DN(0,1)*crRHS86) + crRHS34*(DN(0,0)*crRHS93 + DN(0,1)*crRHS94 + DN(0,2)*crRHS92) + crRHS41*crRHS95 - crRHS65*crRHS96 - crRHS85*crRHS97);
rRHS[3]+=-gauss_weight*(DN(0,0)*crRHS64 + DN(0,1)*crRHS90 + DN(0,2)*crRHS96 + N[0]*crRHS45);
rRHS[4]+=-gauss_weight*(-DN(1,0)*crRHS0 + DN(1,0)*crRHS57 + DN(1,0)*stress_adj[0] + DN(1,1)*stress_adj[3] + DN(1,2)*stress_adj[5] - N[1]*crRHS1 + N[1]*crRHS4 + N[1]*crRHS61 + N[1]*crRHS7 + N[1]*crRHS9 - crRHS100*crRHS64 - crRHS101*crRHS84 - crRHS26*(DN(1,1)*crRHS19 + DN(1,2)*crRHS25) + crRHS3*crRHS98 + crRHS34*(DN(1,0)*crRHS28 + DN(1,1)*crRHS30 + DN(1,2)*crRHS93) + crRHS40*crRHS99 - crRHS64*crRHS98);
rRHS[5]+=-gauss_weight*(DN(1,0)*stress_adj[3] - DN(1,1)*crRHS0 + DN(1,1)*crRHS57 + DN(1,1)*stress_adj[1] + DN(1,2)*stress_adj[4] - N[1]*crRHS66 + N[1]*crRHS67 + N[1]*crRHS71 + N[1]*crRHS72 + N[1]*crRHS73 - crRHS100*crRHS90 - crRHS101*crRHS91 + crRHS26*(DN(1,0)*crRHS19 - DN(1,2)*crRHS86) + crRHS34*(DN(1,0)*crRHS30 + DN(1,1)*crRHS87 + DN(1,2)*crRHS94) + crRHS36*crRHS98 + crRHS89*crRHS99 - crRHS90*crRHS98);
rRHS[6]+=-gauss_weight*(DN(1,0)*stress_adj[5] + DN(1,1)*stress_adj[4] - DN(1,2)*crRHS0 + DN(1,2)*crRHS57 + DN(1,2)*stress_adj[2] - N[1]*crRHS75 + N[1]*crRHS76 + N[1]*crRHS80 + N[1]*crRHS81 + N[1]*crRHS82 - crRHS100*crRHS96 - crRHS101*crRHS97 + crRHS26*(DN(1,0)*crRHS25 + DN(1,1)*crRHS86) + crRHS34*(DN(1,0)*crRHS93 + DN(1,1)*crRHS94 + DN(1,2)*crRHS92) + crRHS38*crRHS98 + crRHS95*crRHS99 - crRHS96*crRHS98);
rRHS[7]+=-gauss_weight*(DN(1,0)*crRHS64 + DN(1,1)*crRHS90 + DN(1,2)*crRHS96 + N[1]*crRHS45);
rRHS[8]+=-gauss_weight*(-DN(2,0)*crRHS0 + DN(2,0)*crRHS57 + DN(2,0)*stress_adj[0] + DN(2,1)*stress_adj[3] + DN(2,2)*stress_adj[5] - N[2]*crRHS1 + N[2]*crRHS4 + N[2]*crRHS61 + N[2]*crRHS7 + N[2]*crRHS9 + crRHS102*crRHS3 - crRHS102*crRHS64 + crRHS103*crRHS40 - crRHS104*crRHS64 - crRHS105*crRHS84 - crRHS26*(DN(2,1)*crRHS19 + DN(2,2)*crRHS25) + crRHS34*(DN(2,0)*crRHS28 + DN(2,1)*crRHS30 + DN(2,2)*crRHS93));
rRHS[9]+=-gauss_weight*(DN(2,0)*stress_adj[3] - DN(2,1)*crRHS0 + DN(2,1)*crRHS57 + DN(2,1)*stress_adj[1] + DN(2,2)*stress_adj[4] - N[2]*crRHS66 + N[2]*crRHS67 + N[2]*crRHS71 + N[2]*crRHS72 + N[2]*crRHS73 + crRHS102*crRHS36 - crRHS102*crRHS90 + crRHS103*crRHS89 - crRHS104*crRHS90 - crRHS105*crRHS91 + crRHS26*(DN(2,0)*crRHS19 - DN(2,2)*crRHS86) + crRHS34*(DN(2,0)*crRHS30 + DN(2,1)*crRHS87 + DN(2,2)*crRHS94));
rRHS[10]+=-gauss_weight*(DN(2,0)*stress_adj[5] + DN(2,1)*stress_adj[4] - DN(2,2)*crRHS0 + DN(2,2)*crRHS57 + DN(2,2)*stress_adj[2] - N[2]*crRHS75 + N[2]*crRHS76 + N[2]*crRHS80 + N[2]*crRHS81 + N[2]*crRHS82 + crRHS102*crRHS38 - crRHS102*crRHS96 + crRHS103*crRHS95 - crRHS104*crRHS96 - crRHS105*crRHS97 + crRHS26*(DN(2,0)*crRHS25 + DN(2,1)*crRHS86) + crRHS34*(DN(2,0)*crRHS93 + DN(2,1)*crRHS94 + DN(2,2)*crRHS92));
rRHS[11]+=-gauss_weight*(DN(2,0)*crRHS64 + DN(2,1)*crRHS90 + DN(2,2)*crRHS96 + N[2]*crRHS45);
rRHS[12]+=-gauss_weight*(-DN(3,0)*crRHS0 + DN(3,0)*crRHS57 + DN(3,0)*stress_adj[0] + DN(3,1)*stress_adj[3] + DN(3,2)*stress_adj[5] - N[3]*crRHS1 + N[3]*crRHS4 + N[3]*crRHS61 + N[3]*crRHS7 + N[3]*crRHS9 + crRHS106*crRHS3 - crRHS106*crRHS64 + crRHS107*crRHS40 - crRHS108*crRHS64 - crRHS109*crRHS84 - crRHS26*(DN(3,1)*crRHS19 + DN(3,2)*crRHS25) + crRHS34*(DN(3,0)*crRHS28 + DN(3,1)*crRHS30 + DN(3,2)*crRHS93));
rRHS[13]+=-gauss_weight*(DN(3,0)*stress_adj[3] - DN(3,1)*crRHS0 + DN(3,1)*crRHS57 + DN(3,1)*stress_adj[1] + DN(3,2)*stress_adj[4] - N[3]*crRHS66 + N[3]*crRHS67 + N[3]*crRHS71 + N[3]*crRHS72 + N[3]*crRHS73 + crRHS106*crRHS36 - crRHS106*crRHS90 + crRHS107*crRHS89 - crRHS108*crRHS90 - crRHS109*crRHS91 + crRHS26*(DN(3,0)*crRHS19 - DN(3,2)*crRHS86) + crRHS34*(DN(3,0)*crRHS30 + DN(3,1)*crRHS87 + DN(3,2)*crRHS94));
rRHS[14]+=-gauss_weight*(DN(3,0)*stress_adj[5] + DN(3,1)*stress_adj[4] - DN(3,2)*crRHS0 + DN(3,2)*crRHS57 + DN(3,2)*stress_adj[2] - N[3]*crRHS75 + N[3]*crRHS76 + N[3]*crRHS80 + N[3]*crRHS81 + N[3]*crRHS82 + crRHS106*crRHS38 - crRHS106*crRHS96 + crRHS107*crRHS95 - crRHS108*crRHS96 - crRHS109*crRHS97 + crRHS26*(DN(3,0)*crRHS25 + DN(3,1)*crRHS86) + crRHS34*(DN(3,0)*crRHS93 + DN(3,1)*crRHS94 + DN(3,2)*crRHS92));
rRHS[15]+=-gauss_weight*(DN(3,0)*crRHS64 + DN(3,1)*crRHS90 + DN(3,2)*crRHS96 + N[3]*crRHS45);

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