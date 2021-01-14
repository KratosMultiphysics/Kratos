//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//

// Application includes
#include "qs_navier_stokes_semiexplicit.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesSemiExplicit<TDim,TNumNodes>::QSNavierStokesSemiExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}

/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
QSNavierStokesSemiExplicit<TDim,TNumNodes>::QSNavierStokesSemiExplicit(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            // this->FractionalVelocityEquationIdVector(rResult,rCurrentProcessInfo);
            this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        case 3:
        {
            this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
        {
        case 1:
        {
            // this->GetFractionalVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 3:
        {
            this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // TODO: we should check and correct the shapes.

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 3:
        {
            this->CalculateLocalPressureSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
            break;
        }
        case 4:
        {
            this->CalculateLocalEndOfStepSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateLocalFractionalVelocitySystem(
    BoundedVector<double, 6> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    //substitute_rhs_momentum_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLocalFractionalVelocitySystem(
    BoundedVector<double, 12> &rRightHandSideBoundedVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    //substitute_rhs_momentum_3D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideBoundedVector *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 3;

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != n_nodes )
        rLeftHandSideMatrix.resize(n_nodes,n_nodes);
    rLeftHandSideMatrix = ZeroMatrix(n_nodes,n_nodes);
    if( rRightHandSideVector.size() != n_nodes )
        rRightHandSideVector.resize(n_nodes);
    rRightHandSideVector = ZeroVector(n_nodes);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);

    //substitute_rhs_mass_2D
    //substitute_lhs_mass_2D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLocalPressureSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    constexpr unsigned int n_nodes = 4;

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != n_nodes )
        rLeftHandSideMatrix.resize(n_nodes,n_nodes);
    rLeftHandSideMatrix = ZeroMatrix(n_nodes,n_nodes);
    if( rRightHandSideVector.size() != n_nodes )
        rRightHandSideVector.resize(n_nodes);
    rRightHandSideVector = ZeroVector(n_nodes);

    // Struct to pass around the data
    ElementDataStruct data;
    this->FillElementData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const auto &dt = data.dt;
    const auto &f = data.forcing;
    const auto &fracv = data.fractional_velocity;
    const auto &fracvconv = data.fractional_convective_velocity;

    const auto &gamma = data.gamma;
    const auto &h = data.h;
    const auto &nu = data.nu;
    const auto &p = data.pressure;
    const auto &pn = data.pressure_old;
    const auto &rho = data.rho;
    const auto &v = data.velocity;
    const auto &vn = data.velocity_old;
    const auto &vconv = data.velocity_convective;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    // Hardcoded shape functions gradients for linear triangular element
    // This is explicitly done to minimize the matrix acceses
    // The notation DN_i_j means shape function for node i in dimension j
    const double &DN_DX_0_0 = data.DN_DX(0, 0);
    const double &DN_DX_0_1 = data.DN_DX(0, 1);
    const double &DN_DX_0_2 = data.DN_DX(0, 2);
    const double &DN_DX_1_0 = data.DN_DX(1, 0);
    const double &DN_DX_1_1 = data.DN_DX(1, 1);
    const double &DN_DX_1_2 = data.DN_DX(1, 2);
    const double &DN_DX_2_0 = data.DN_DX(2, 0);
    const double &DN_DX_2_1 = data.DN_DX(2, 1);
    const double &DN_DX_2_2 = data.DN_DX(2, 2);
    const double &DN_DX_3_0 = data.DN_DX(3, 0);
    const double &DN_DX_3_1 = data.DN_DX(3, 1);
    const double &DN_DX_3_2 = data.DN_DX(3, 2);

    //substitute_rhs_mass_3D
    //substitute_lhs_mass_3D

    // Here we assume that all the weights of the gauss points are the same so we multiply at the end by Volume/n_nodes
    rRightHandSideVector *= data.volume / static_cast<double>(n_nodes);
    rLeftHandSideMatrix *= data.volume / static_cast<double>(n_nodes);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateLocalEndOfStepSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    KRATOS_ERROR << "Calling the CalculateLocalEndOfStepSystem() method for the semiexplicit Navier-Stokes element. You should call the Calculate method instead.";

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            constexpr IndexType dim = 2;
            constexpr IndexType n_nodes = 3;
            constexpr IndexType block_size = 2;

            // Calculate the explicit residual vector
            BoundedVector<double, 6> rhs;
            this->CalculateLocalFractionalVelocitySystem(rhs, rCurrentProcessInfo);

            // Add the residual contribution
            // Note that the reaction is indeed the formulation residual
            auto& r_geometry = GetGeometry();
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                const IndexType aux = i_node * block_size;
                auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
                for (IndexType d = 0; d < dim; ++d) {
                    #pragma omp atomic
                    r_mom[d] += rhs[aux + d];
                }
            }
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::AddExplicitContribution(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
    {
        case 1:
        {
            constexpr IndexType dim = 3;
            constexpr IndexType n_nodes = 4;
            constexpr IndexType block_size = 3;

            // Calculate the explicit residual vector
            BoundedVector<double, 12> rhs;
            this->CalculateLocalFractionalVelocitySystem(rhs, rCurrentProcessInfo);

            // Add the residual contribution
            // Note that the reaction is indeed the formulation residual
            auto& r_geometry = GetGeometry();
            for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
                const IndexType aux = i_node * block_size;
                auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
                for (IndexType d = 0; d < dim; ++d) {
                    #pragma omp atomic
                    r_mom[d] += rhs[aux + d];
                }
            }
            break;
        }
        default:
        {
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ", rCurrentProcessInfo[FRACTIONAL_STEP]);
        }
    }


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 2;

    // Initialize and fill the mass matrix values
    const double one_six = 1.0 / 6.0;
    const double one_twelve = 1.0 / 12.0;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_six;    rMassMatrix(0, 2) = one_twelve; rMassMatrix(0, 4) = one_twelve;
    rMassMatrix(1, 1) = one_six;    rMassMatrix(1, 3) = one_twelve; rMassMatrix(1, 5) = one_twelve;
    rMassMatrix(2, 0) = one_twelve; rMassMatrix(2, 2) = one_six;    rMassMatrix(2, 4) = one_twelve;
    rMassMatrix(3, 1) = one_twelve; rMassMatrix(3, 3) = one_six;    rMassMatrix(3, 5) = one_twelve;
    rMassMatrix(4, 0) = one_twelve; rMassMatrix(4, 2) = one_twelve; rMassMatrix(4, 4) = one_six;
    rMassMatrix(5, 1) = one_twelve; rMassMatrix(5, 3) = one_twelve; rMassMatrix(5, 5) = one_six;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateMassMatrix(
    MatrixType &rMassMatrix,
    const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType n_nodes = 4;
    constexpr IndexType block_size = 3;

    // Initialize and fill the mass matrix values
    const double one_ten = 0.1;
    const double one_twenty = 0.05;
    const unsigned int size = n_nodes * block_size;
    rMassMatrix = ZeroMatrix(size, size);
    rMassMatrix(0, 0) = one_ten;     rMassMatrix(0, 3) = one_twenty;  rMassMatrix(0, 6) = one_twenty;  rMassMatrix(0,9) = one_twenty;
    rMassMatrix(1, 1) = one_ten;     rMassMatrix(1, 4) = one_twenty;  rMassMatrix(1, 7) = one_twenty;  rMassMatrix(1,10) = one_twenty;
    rMassMatrix(2, 2) = one_ten;     rMassMatrix(2, 5) = one_twenty;  rMassMatrix(2, 8) = one_twenty;  rMassMatrix(2,11) = one_twenty;
    rMassMatrix(3, 0) = one_twenty;  rMassMatrix(3, 3) = one_ten;     rMassMatrix(3, 6) = one_twenty;  rMassMatrix(3,9) = one_twenty;
    rMassMatrix(4, 1) = one_twenty;  rMassMatrix(4, 4) = one_ten;     rMassMatrix(4, 7) = one_twenty;  rMassMatrix(4,10) = one_twenty;
    rMassMatrix(5, 2) = one_twenty;  rMassMatrix(5, 5) = one_ten;     rMassMatrix(5, 8) = one_twenty;  rMassMatrix(5,11) = one_twenty;
    rMassMatrix(6, 0) = one_twenty;  rMassMatrix(6, 3) = one_twenty;  rMassMatrix(6, 6) = one_ten;     rMassMatrix(6,9) = one_twenty;
    rMassMatrix(7, 1) = one_twenty;  rMassMatrix(7, 4) = one_twenty;  rMassMatrix(7, 7) = one_ten;     rMassMatrix(7,10) = one_twenty;
    rMassMatrix(8, 2) = one_twenty;  rMassMatrix(8, 5) = one_twenty;  rMassMatrix(8, 8) = one_ten;     rMassMatrix(8,11) = one_twenty;
    rMassMatrix(9, 0) = one_twenty;  rMassMatrix(9, 3) = one_twenty;  rMassMatrix(9, 6) = one_twenty;  rMassMatrix(9,9) = one_ten;
    rMassMatrix(10, 1) = one_twenty; rMassMatrix(10, 4) = one_twenty; rMassMatrix(10, 7) = one_twenty; rMassMatrix(10,10) = one_ten;
    rMassMatrix(11, 2) = one_twenty; rMassMatrix(11, 5) = one_twenty; rMassMatrix(11, 8) = one_twenty; rMassMatrix(11,11) = one_ten;

    // Here we assume that all the Gauss points have the same weight so we multiply by the volume
    rMassMatrix *= GetGeometry().Area();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Define local variables
    const unsigned int local_size = 6; // two components velocity * three nodes
    const double one_third_area = GetGeometry().Area() / 3.0;
    // Initialize and calculate elemental lumped mass vector
    if (rLumpedMassVector.size() != local_size) {
        rLumpedMassVector.resize(local_size, false);
    }
    for (IndexType i = 0; i < local_size; ++i) {
        rLumpedMassVector(i) = one_third_area;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::CalculateLumpedMassVector(
    VectorType& rLumpedMassVector,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Define local variables
    const unsigned int local_size = 12; // three components velocity * four nodes
    const double one_fourth_volume = GetGeometry().Volume() / 4.0;
    // Initialize and calculate elemental lumped mass vector
    if (rLumpedMassVector.size() != local_size) {
        rLumpedMassVector.resize(local_size, false);
    }
    for (IndexType i = 0; i < local_size; ++i) {
        rLumpedMassVector(i) = one_fourth_volume;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::Calculate(
    const Variable<array_1d<double,3> > &rVariable,
    array_1d<double,3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == VELOCITY)
    {
        constexpr unsigned int nodes_per_dimension = 6;

        // Check sizes and initialize
        VectorType rRightHandSideVector;
        if( rRightHandSideVector.size() != nodes_per_dimension )
            rRightHandSideVector.resize(nodes_per_dimension);
        rRightHandSideVector = ZeroVector(nodes_per_dimension);

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Substitute the formulation symbols by the data structure values
        const auto &dt = data.dt;
        const auto &f = data.forcing;
        const auto &fracv = data.fractional_velocity;
        const auto &fracvconv = data.fractional_convective_velocity;

        const auto &gamma = data.gamma;
        const auto &h = data.h;
        const auto &nu = data.nu;
        const auto &p = data.pressure;
        const auto &pn = data.pressure_old;
        const auto &rho = data.rho;
        const auto &v = data.velocity;
        const auto &vn = data.velocity_old;
        const auto &vconv = data.velocity_convective;

        // Stabilization parameters
        const double stab_c1 = 4.0;
        const double stab_c2 = 2.0;

        // Hardcoded shape functions gradients for linear triangular element
        // This is explicitly done to minimize the matrix acceses
        // The notation DN_i_j means shape function for node i in dimension j
        const double &DN_DX_0_0 = data.DN_DX(0, 0);
        const double &DN_DX_0_1 = data.DN_DX(0, 1);
        const double &DN_DX_1_0 = data.DN_DX(1, 0);
        const double &DN_DX_1_1 = data.DN_DX(1, 1);
        const double &DN_DX_2_0 = data.DN_DX(2, 0);
        const double &DN_DX_2_1 = data.DN_DX(2, 1);

        //substitute_rhs_endofstep_2D

        SizeType Index = 0;
        GeometryType& rGeom = this->GetGeometry();

        for (SizeType i = 0; i < 3; ++i) // loop over nodes
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < 2; ++d) // loop over dimensions
            {
                rTemp[d] += rRightHandSideVector[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::Calculate(
    const Variable<array_1d<double,3> > &rVariable,
    array_1d<double,3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == VELOCITY)
    {

        constexpr unsigned int nodes_per_dimension = 12;

        // Check sizes and initialize
        VectorType rRightHandSideVector;
        if( rRightHandSideVector.size() != nodes_per_dimension )
            rRightHandSideVector.resize(nodes_per_dimension);
        rRightHandSideVector = ZeroVector(nodes_per_dimension);

        // Struct to pass around the data
        ElementDataStruct data;
        this->FillElementData(data, rCurrentProcessInfo);

        // Substitute the formulation symbols by the data structure values
        const auto &dt = data.dt;
        const auto &f = data.forcing;
        const auto &fracv = data.fractional_velocity;
        const auto &fracvconv = data.fractional_convective_velocity;

        const auto &gamma = data.gamma;
        const auto &h = data.h;
        const auto &nu = data.nu;
        const auto &p = data.pressure;
        const auto &pn = data.pressure_old;
        const auto &rho = data.rho;
        const auto &v = data.velocity;
        const auto &vn = data.velocity_old;
        const auto &vconv = data.velocity_convective;

        // Stabilization parameters
        const double stab_c1 = 4.0;
        const double stab_c2 = 2.0;

        // Hardcoded shape functions gradients for linear triangular element
        // This is explicitly done to minimize the matrix acceses
        // The notation DN_i_j means shape function for node i in dimension j
        const double &DN_DX_0_0 = data.DN_DX(0, 0);
        const double &DN_DX_0_1 = data.DN_DX(0, 1);
        const double &DN_DX_0_2 = data.DN_DX(0, 2);
        const double &DN_DX_1_0 = data.DN_DX(1, 0);
        const double &DN_DX_1_1 = data.DN_DX(1, 1);
        const double &DN_DX_1_2 = data.DN_DX(1, 2);
        const double &DN_DX_2_0 = data.DN_DX(2, 0);
        const double &DN_DX_2_1 = data.DN_DX(2, 1);
        const double &DN_DX_2_2 = data.DN_DX(2, 2);
        const double &DN_DX_3_0 = data.DN_DX(3, 0);
        const double &DN_DX_3_1 = data.DN_DX(3, 1);
        const double &DN_DX_3_2 = data.DN_DX(3, 2);

        //substitute_rhs_endofstep_3D

        SizeType Index = 0;
        GeometryType& rGeom = this->GetGeometry();

        for (SizeType i = 0; i < 4; ++i) // loopover nodes
        {
            rGeom[i].SetLock(); // So it is safe to write in the node in OpenMP
            array_1d<double,3>& rTemp = rGeom[i].FastGetSolutionStepValue(FRACT_VEL);
            for (SizeType d = 0; d < 3; ++d) // loop over dimensions
            {
                rTemp[d] += rRightHandSideVector[Index++];
            }
            rGeom[i].UnSetLock(); // Free the node for other threads
        }

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
int QSNavierStokesSemiExplicit<TDim,TNumNodes>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;
    return 0;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::FillElementData(
    ElementDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Getting data for the given geometry and integration method
    const auto& r_geometry = GetGeometry();
    const unsigned int local_size = r_geometry.size();
    array_1d<double,TNumNodes> N_aux;
    GeometryUtils::CalculateGeometryData(r_geometry,rData.DN_DX,N_aux,rData.volume);
    rData.N_gausspoint = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Initialize some scalar data
    rData.lumping_factor = 1.00 / double(TNumNodes);
    rData.mu = 0.0;
    rData.dt = rCurrentProcessInfo[DELTA_TIME];
    // rData.dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    rData.rho = 0.0;
    // Commented code should become useful when adding stabilization and to replace theta computation
    double theta = rCurrentProcessInfo[RUNGE_KUTTA_STEP];
    if (theta == 1) {
        theta = 0.0;
    }
    else if (theta == 2 or theta == 3) {
        theta = 0.5;
    }
    else {
        theta = 1;
    }
    // double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
    // if (theta == 0.0) {
    //     rData.explicit_step_coefficient = 0.0;
    // }
    // else {
    //     rData.explicit_step_coefficient = 1.0/((theta)*r_process_info[DELTA_TIME]);
    // }

    for(unsigned int node_element = 0; node_element<local_size; node_element++) {
        // Assign variables to local variables
        const auto& r_body_force = r_geometry[node_element].FastGetSolutionStepValue(BODY_FORCE);
        const auto& r_body_force_old = r_geometry[node_element].FastGetSolutionStepValue(BODY_FORCE,1);
        const auto& r_velocity = r_geometry[node_element].FastGetSolutionStepValue(VELOCITY);
        const auto& r_velocity_old = r_geometry[node_element].FastGetSolutionStepValue(VELOCITY,1);
        // Observations
        // * unknown acceleration approximated as (u-u_old)*explicit_step_coefficient = (u-u_old)/((theta)*dt)
        //   observe that for theta = 0.0, u = u_old and explicit_step_coefficient = 0

        for (unsigned int k = 0; k < TDim; ++k) {
            // forcing term: interpolation exploiting theta
            rData.forcing(node_element,k) = (1-theta) * r_body_force_old[k] + theta * r_body_force[k];
            // velocity previous time step
            rData.velocity_old(node_element,k) = r_velocity_old[k];
            // fractional velocity current time step
            // observe that VELOCITY and not FRACT_VEL is being used as DoF, to avoid saving additional DoFs
            rData.fractional_velocity(node_element,k) = r_velocity[k];
            // convective fractional velocity is fractional velocity of current explicit step
            // observe that VELOCITY and not FRACT_VEL is being used as dof, to avoid saving additional dofs
            rData.fractional_convective_velocity(node_element,k) = r_velocity[k];
        }
        // pressure current and previous time step
        rData.pressure[node_element] = r_geometry[node_element].FastGetSolutionStepValue(PRESSURE);
        rData.pressure_old[node_element] = r_geometry[node_element].FastGetSolutionStepValue(PRESSURE,1);
        // rData.oss_projection[node_element] = r_geometry[node_element].FastGetSolutionStepValue();
        rData.mu += r_geometry[node_element].FastGetSolutionStepValue(VISCOSITY);
        rData.rho += r_geometry[node_element].FastGetSolutionStepValue(DENSITY);
    }
    // divide by number of nodes scalar data
    rData.mu *= rData.lumping_factor;
    rData.rho *= rData.lumping_factor;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
double QSNavierStokesSemiExplicit<TDim,TNumNodes>::CalculateElementSize(
    const BoundedMatrix<double,TNumNodes, TDim>& rDN_DX)
{
    KRATOS_TRY;

    double h = 0.0;
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        double h_inv = 0.0;
        for (unsigned int k = 0; k < TDim; ++k) {
            h_inv += rDN_DX(i,k) * rDN_DX(i,k);
        }
        h += 1.0/h_inv;
    }
    h = sqrt(h) / static_cast<double>(TNumNodes);
    return h;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::FractionalVelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACT_VEL_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y,xpos+1).EquationId();
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::FractionalVelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(FRACT_VEL_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Z,xpos+2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::VelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::VelocityEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::PressureEquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rResult.size() != NumNodes)
        rResult.resize(NumNodes);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

    for (SizeType i = 0; i < NumNodes; ++i)
        rResult[i] = rGeom[i].GetDof(PRESSURE,pos).EquationId();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::GetFractionalVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::GetFractionalVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<2,3>::GetVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
    }
}

/***********************************************************************************/

template<>
void QSNavierStokesSemiExplicit<3,4>::GetVelocityDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TNumNodes>
void QSNavierStokesSemiExplicit<TDim,TNumNodes>::GetPressureDofList(DofsVectorType& rElementalDofList,
                                              const ProcessInfo& rCurrentProcessInfo) const
{
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rElementalDofList.size() != NumNodes)
        rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
    {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class QSNavierStokesSemiExplicit<2,3>;
template class QSNavierStokesSemiExplicit<3,4>;

/***********************************************************************************/
/***********************************************************************************/

}