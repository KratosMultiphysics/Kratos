// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "adjoint_finite_difference_membrane_element.h"
#include "structural_mechanics_application_variables.h"
#include "custom_response_functions/response_utilities/stress_response_definitions.h"
#include "includes/checks.h"
#include "custom_elements/membrane_elements/membrane_element.hpp"
#include "custom_response_functions/response_utilities/finite_difference_utility.h"


namespace Kratos
{

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // Get perturbation size
    const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
    const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = this->mHasRotationDofs ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;

    if( rDesignVariable == TEMPERATURE )
    {
        rOutput.resize( number_of_nodes, local_size, false);

        Vector RHS;
        Vector derived_RHS;
        
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

        for(IndexType i_node = 0; i_node < this->mpPrimalElement->GetGeometry().PointsNumber(); ++i_node)
        {
            // Get pseudo-load contribution from utility
            FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, rDesignVariable,
                                                                      this->mpPrimalElement->GetGeometry()[i_node], delta, 
                                                                      derived_RHS, rCurrentProcessInfo);

            KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit! [ derived_RHS.size() = " << derived_RHS.size() << ", local_size = " << local_size << " ]." << std::endl;

            for(IndexType i = 0; i < derived_RHS.size(); ++i)
                rOutput(i_node, i) = derived_RHS[i];
        }
    }
    // else if( rDesignVariable == TRUSS_PRESTRESS_PK2 )
    // {
    //     Vector RHS;
    //     this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);

    //     // Get pseudo-load from utility
    //     //FiniteDifferenceUtility::CalculateRightHandSideDerivative(*pGetPrimalElement(), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);
    //     rOutput.resize(1,RHS.size(), false);

    //     std::stringstream filename;
    //     filename << "sensitivity_element_" << mpPrimalElement->Id() << ".dat";

    //     std::ifstream inFile(filename.str());
    //     if (!inFile.is_open()) {
    //         KRATOS_ERROR << "Could not open file " << filename.str() << " for reading" << std::endl;
    //     }

    //     std::vector<double> values;
    //     double val;
    //     while (inFile >> val) {
    //         values.push_back(val);
    //     }
    //     inFile.close();

    //     // Create 1-row matrix
        
    //     for (std::size_t j = 0; j < values.size(); ++j) {
    //         rOutput(0, j) = values[j];
    //     }

    
    // }
    else
    {
        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        //KRATOS_WATCH(RHS)
        //KRATOS_WATCH(this->pGetPrimalElement())

        // Get pseudo-load from utility
        FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, rDesignVariable, delta, rOutput, rCurrentProcessInfo);
    }

    if (rOutput.size1() == 0 || rOutput.size2() == 0)
    {
        //std::cout << "Passing 0 of size 1 to sensitivity matrix for variable: " << rDesignVariable.Name() << std::endl;
        rOutput = ZeroMatrix(1, local_size);
    }

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    
    const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
    const SizeType num_dofs_per_node = (this->mHasRotationDofs) ?  2 * dimension : dimension;
    const SizeType local_size = number_of_nodes * num_dofs_per_node;
    
    if( rDesignVariable == SHAPE_SENSITIVITY )
    {
        const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);
        const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
        Vector derived_RHS;

        if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
            rOutput.resize(dimension * number_of_nodes, local_size, false);

        IndexType index = 0;

        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        for(auto& node_i : this->mpPrimalElement->GetGeometry())
        {
            for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
            {
                // Get pseudo-load contribution from utility
                FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, *coord_directions[coord_dir_i],
                                                                            node_i, delta, derived_RHS, rCurrentProcessInfo);

                KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

                for(IndexType i = 0; i < derived_RHS.size(); ++i)
                    rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
            }
            index++;
        }
    }
    else if( rDesignVariable == PRE_STRESS )
    {
        //std::cout << "CalculateSensitivityMatrix in MembraneElement" << std::endl;
        const Vector delta_vector = this->GetPerturbationSize(PRESTRESS_VECTOR, rCurrentProcessInfo);
        //const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> components = {&PRE_STRESS_SENSITIVITY_XX, &PRE_STRESS_SENSITIVITY_YY, &PRE_STRESS_SENSITIVITY_XY};
        Matrix derived_RHS;
        const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
        
        KRATOS_ERROR_IF_NOT(dimension > 1) << "CalculateSensitivityMatrix for Vector variables is only available for 2 and 3D!" << std::endl;
        
        SizeType vector_size;
        if (dimension == 2){
            vector_size = 1;
        }
        else if (dimension == 3){
            vector_size = 3;
        }

        if ( (rOutput.size1() != vector_size) || (rOutput.size2() != local_size ) )
            rOutput.resize(vector_size , local_size, false);

        IndexType index = 0;

        Vector RHS;
        this->pGetPrimalElement()->CalculateRightHandSide(RHS, rCurrentProcessInfo);
        
        // Get pseudo-load contribution from utility
        FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->pGetPrimalElement()), RHS, PRESTRESS_VECTOR,
                                                                    delta_vector, derived_RHS, rCurrentProcessInfo);

        KRATOS_ERROR_IF_NOT(derived_RHS.size2() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

        for(IndexType i = 0; i < derived_RHS.size1(); ++i)
            for(IndexType j = 0; j < derived_RHS.size2(); ++j)
                if (i==2){
                    rOutput(i,j) = 0.0;
                }
                else {
                    rOutput( i , j) = derived_RHS(i,j);
                }

    }
    else{
        //std::cout << "Passing 0 of size 3 to sensitivity matrix for variable: " << rDesignVariable.Name() << std::endl;
        rOutput = ZeroMatrix(3, local_size);
    }
    KRATOS_CATCH("")
}

template <class TPrimalElement>
int AdjointFiniteDifferenceMembraneElement<TPrimalElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int return_value = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF(this->mHasRotationDofs) << "Adjoint membrane element have rotation dofs!" << std::endl;
    KRATOS_ERROR_IF_NOT(this->mpPrimalElement) << "Primal element pointer is nullptr!" << std::endl;

    //TODO: Check() of primal element should be called, but is not possible because of DOF check!
    this->CheckDofs();
    this->CheckProperties(rCurrentProcessInfo);

    KRATOS_ERROR_IF(this->GetGeometry().Area() < std::numeric_limits<double>::epsilon()*1000)
        << "Element #" << this->Id() << " has an Area of zero!" << std::endl;

    return return_value;

    KRATOS_CATCH("")
}

// template <class TPrimalElement>
// void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
//                                        const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY;

//     //const double delta = this->GetPerturbationSize(rDesignVariable, rCurrentProcessInfo);

//     const SizeType number_of_nodes = this->mpPrimalElement->GetGeometry().PointsNumber();
//     const SizeType dimension = rCurrentProcessInfo.GetValue(DOMAIN_SIZE);
//     const SizeType num_dofs_per_node = (this->mHasRotationDofs) ?  2 * dimension : dimension;
//     const SizeType local_size = number_of_nodes * num_dofs_per_node;

//     //const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};
//     Vector derived_RHS;

//     // if ( (rOutput.size1() != dimension * number_of_nodes) || (rOutput.size2() != local_size ) )
//     //     rOutput.resize(dimension * number_of_nodes, local_size, false);

//     // IndexType index = 0;

//     Vector RHS;
//     this->mpPrimalElement->CalculateRightHandSide(RHS, rCurrentProcessInfo);
//     // for(auto& node_i : mpPrimalElement->GetGeometry())
//     // {
//     //     for(IndexType coord_dir_i = 0; coord_dir_i < dimension; ++coord_dir_i)
//     //     {
//     //         // Get pseudo-load contribution from utility
//     //         FiniteDifferenceUtility::CalculateRightHandSideDerivative(*(this->mpPrimalElement), RHS, *coord_directions[coord_dir_i],
//     //                                                                     node_i, delta, derived_RHS, rCurrentProcessInfo);

//     //         KRATOS_ERROR_IF_NOT(derived_RHS.size() == local_size) << "Size of the pseudo-load does not fit!" << std::endl;

//     //         for(IndexType i = 0; i < derived_RHS.size(); ++i)
//     //             rOutput( (coord_dir_i + index*dimension), i) = derived_RHS[i];
//     //     }
//     //     index++;
//     // }

//     rLeftHandSideMatrix = ZeroMatrix(0, local_size);

//     KRATOS_CATCH("")
// }

// private

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::CheckDofs() const
{
    const GeometryType& r_geom = this->GetGeometry();
    // verify that the dofs exist
    for (IndexType i = 0; i < r_geom.size(); ++i)
    {
        const auto& r_node = r_geom[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, r_node);

        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_X, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Y, r_node);
        KRATOS_CHECK_DOF_IN_NODE(ADJOINT_DISPLACEMENT_Z, r_node);

        KRATOS_ERROR_IF(r_node.GetBufferSize() < 2) << "This Element needs "
            << "at least a buffer size = 2" << std::endl;
    }
}

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::CheckProperties(const ProcessInfo& rCurrentProcessInfo) const
{
    // check properties
    if(this->pGetProperties() == nullptr)
        KRATOS_ERROR << "Properties not provided for element " << this->Id() << std::endl;

    //const GeometryType& geom = this->GetGeometry(); // TODO check if this can be const

    const PropertiesType & r_props = this->GetProperties();

    if (!r_props.Has(CONSTITUTIVE_LAW))
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;
    const ConstitutiveLaw::Pointer& claw = r_props[CONSTITUTIVE_LAW];
    if (claw == nullptr)
        KRATOS_ERROR << "CONSTITUTIVE_LAW not provided for element " << this->Id() << std::endl;

    if(!r_props.Has(THICKNESS))
        KRATOS_ERROR << "THICKNESS not provided for element " << this->Id() << std::endl;
    if(r_props[THICKNESS] <= 0.0)
        KRATOS_ERROR << "wrong THICKNESS value provided for element " << this->Id() << std::endl;

    if(!r_props.Has(DENSITY))
        KRATOS_ERROR << "DENSITY not provided for element " << this->Id() << std::endl;
    if(r_props[DENSITY] < 0.0)
        KRATOS_ERROR << "wrong DENSITY value provided for element " << this->Id() << std::endl;

}

template <class TPrimalElement>
double AdjointFiniteDifferenceMembraneElement<TPrimalElement>::GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const
{
    KRATOS_TRY;

    if(rDesignVariable == SHAPE_SENSITIVITY)
    {
        double dx, dy, dz, L = 0.0;

        const GeometryType& geometry = this->mpPrimalElement->GetGeometry();

        dx = geometry[1].X0() - geometry[0].X0();
        dy = geometry[1].Y0() - geometry[0].Y0();
        dz = geometry[1].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[1].X0();
        dy = geometry[2].Y0() - geometry[1].Y0();
        dz = geometry[2].Z0() - geometry[1].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        dx = geometry[2].X0() - geometry[0].X0();
        dy = geometry[2].Y0() - geometry[0].Y0();
        dz = geometry[2].Z0() - geometry[0].Z0();
        L += sqrt(dx*dx + dy*dy + dz*dz);
        L /= 3.0;

        return L;
    }
    else
        return 1.0;

    KRATOS_CATCH("")
}

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

template <class TPrimalElement>
void AdjointFiniteDifferenceMembraneElement<TPrimalElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType);

}

template class AdjointFiniteDifferenceMembraneElement<MembraneElement>;

} // namespace Kratos

