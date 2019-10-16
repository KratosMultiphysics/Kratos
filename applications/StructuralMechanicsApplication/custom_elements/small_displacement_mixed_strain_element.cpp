// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/small_displacement_mixed_strain_element.h"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedStrainElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedStrainElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedStrainElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size  = n_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    if (dim == 2) {
        for(unsigned int i = 0; i < n_nodes; ++i) {
            rResult[i * (dim + 1)] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[i * (dim + 1) + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 2] = this->GetGeometry()[i].GetDof(VOLUMETRIC_STRAIN).EquationId();
        }
    } else if (dim == 3) {
        for(unsigned int i = 0; i < n_nodes; ++i){
            rResult[i * (dim + 1)] = this->GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[i * (dim + 1) + 1] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 2] = this->GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[i * (dim + 1) + 3] = this->GetGeometry()[i].GetDof(VOLUMETRIC_STRAIN).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto &r_geometry = GetGeometry();
    const unsigned int n_nodes = r_geometry.PointsNumber();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();
    const unsigned int dof_size  = n_nodes*(dim+1);

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(unsigned int i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    } else if (dim == 3) {
        for(unsigned int i = 0; i < n_nodes; ++i){
            rElementalDofList[i * (dim + 1)] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * (dim + 1) + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * (dim + 1) + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * (dim + 1) + 3] = this->GetGeometry()[i].pGetDof(VOLUMETRIC_STRAIN);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::Initialize()
{
    KRATOS_TRY

    // Integration method initialization
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
    const auto &r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Constitutive Law Vector initialisation
    if (mConstitutiveLawVector.size() != r_integration_points.size()) {
        mConstitutiveLawVector.resize(r_integration_points.size());
    }

    // Initialize material
    InitializeMaterial();

    KRATOS_CATCH( "" )
}

void SmallDisplacementMixedStrainElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

void SmallDisplacementMixedStrainElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Set te constitutive law values
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    Vector strain(strain_size);
    Vector stress(strain_size);
    Matrix cons_matrix(strain_size, strain_size);
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto &r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    cons_law_values.SetStrainVector(strain);
    cons_law_values.SetStressVector(stress);
    cons_law_values.SetConstitutiveMatrix(cons_matrix);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // TODO: IMPLEMENT IN A MORE EFFICIENT MANNER
    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Compute the geometry data
    Vector w_gauss_container;
    Matrix N_gauss_container;
    GeometryType::ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(r_geometry, w_gauss_container, N_gauss_container, DN_DX_container);

    // Calculate the RHS contributions
    rLeftHandSideMatrix.clear();

    Matrix B_mat(strain_size, n_nodes*dim);
    Matrix dev_strain_op;
    CalculateDeviatoricStrainOperator(dev_strain_op);

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Get Gauss pt. values
        const auto &rN = row(N_gauss_container, i_gauss);
        const double w_gauss = w_gauss_container[i_gauss];

        // Add the deviatoric strain contribution
        // TODO: MOVE TO A FUNCTION
        Vector tot_strain = ZeroVector(strain_size);
        CalculateB(B_mat, DN_DX_container[i_gauss]);
        for (unsigned int i = 0; i < strain_size; ++i) {
            for (unsigned int j = 0; j < strain_size ; ++j) {
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int d = 0; d < dim; ++d) {
                        const unsigned int aux = i_node * dim + d;
                        tot_strain[i] += dev_strain_op(i,j) * B_mat(j,aux) * r_disp[d];
                    }
                }
            }
        }

        // Interpolate and add the nodal volumetric strain
        double vol_strain = 0.0;
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            vol_strain += rN[i_node] * r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }
        for (unsigned int d = 0; d < dim; ++d) {
            tot_strain[d] += vol_strain;
        }

        // Get the stress from the constitutive law
        // TODO: MOVE TO A FUNCTION
        Vector tot_stress(strain_size);
        Matrix cons_matrix(strain_size, strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto &r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        cons_law_values.SetShapeFunctionsValues(rN);
        cons_law_values.SetStrainVector(tot_strain);
        cons_law_values.SetStressVector(tot_stress);
        cons_law_values.SetConstitutiveMatrix(cons_matrix);
        mConstitutiveLawVector[i_gauss]->CalculateMaterialResponseCauchy(cons_law_values);

        // Add momentum deviatoric stress contribution
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                for (unsigned int d1 = 0; d1 < dim; ++d1) {
                    const unsigned int aux_i = i * block_size + d1;
                    for (unsigned int d2 = 0; d2 < dim; ++d2) {
                        const unsigned int aux_j = j * block_size + d2;
                        for (unsigned int l = 0; l < strain_size; ++l) {
                            for (unsigned int k = 0; k < strain_size; ++k) {
                                for (unsigned int m = 0; m < strain_size; ++m) {
                                    for (unsigned int n = 0; n < strain_size; ++n) {
                                        rLeftHandSideMatrix(aux_i, aux_j) += w_gauss * B_mat(l,i*dim+d1) * dev_strain_op(k,l) * cons_law_values.GetConstitutiveMatrix()(k,m) * dev_strain_op(m,n) * B_mat(n,j*dim+d2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Add momentum volume stress contribution
        double hyd_stress = 0.0;
        for (unsigned int d = 0; d < dim; ++d) {
            hyd_stress += cons_law_values.GetStressVector()[d];
        }
        hyd_stress /= 3.0;
        const double bulk_modulus = vol_strain > 1.0e-12 ? hyd_stress * vol_strain / std::pow(vol_strain, 2) : 1.0e+12;

        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const unsigned int aux_j = j * block_size + dim;
                for (unsigned int d = 0; d < dim; ++d) {
                    const unsigned int aux_i = i * block_size + d;
                    for (unsigned int l = 0; l < strain_size; ++l) {
                        rLeftHandSideMatrix(aux_i, aux_j) += w_gauss * B_mat(l,i*dim+d) * bulk_modulus * rN[j];
                    }
                }
            }
        }

        // Add mass conservation divergence contribution
        for (unsigned int i = 0; i < n_nodes; ++i) {
            const unsigned int aux_i = i * block_size + dim;
            for (unsigned int j = 0; j < n_nodes; ++j) {
                for (unsigned int d = 0; d < dim; ++d) {
                    const unsigned int aux_j = j * block_size + d;
                    rLeftHandSideMatrix(aux_i, aux_j) += w_gauss * rN[i] * DN_DX_container[i_gauss](j,d);
                }
            }
        }

        // Add mass conservation volumetric strain contribution
        for (unsigned int i = 0; i < n_nodes; ++i) {
            const unsigned int aux_i = i * block_size + dim;
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const unsigned int aux_j = j * block_size + dim;
                rLeftHandSideMatrix(aux_i, aux_j) += w_gauss * rN[i] * rN[j];
            }
        }

        // Calculate stabilization constants
        const double h = ComputeElementSize(DN_DX_container[i_gauss]);
        const double shear_modulus = GetProperties()[YOUNG_MODULUS] / (2.0 * (1.0 + GetProperties()[POISSON_RATIO]));
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = 0.15;

        // Add the volumetric strain momentum stabilization term - term 1
        const double aux_1 = w_gauss * bulk_modulus * tau_2;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d1 = 0; d1 < dim; ++d1) {
                for (unsigned int j = 0; j < n_nodes; ++j) {
                    for (unsigned int d2 = 0; d2 < dim; ++d2) {
                        rLeftHandSideMatrix(i * block_size + d1, j * block_size + d2) += aux_1 * DN_DX_container[i_gauss](i, d1) * DN_DX_container[i_gauss](j, d2);
                    }
                }
            }
        }

        // Add the volumetric strain momentum stabilization term - term 2
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d = 0; d < dim; ++d) {
                for (unsigned int j = 0; j < n_nodes; ++j) {
                    rLeftHandSideMatrix(i * block_size + d, j * block_size + dim) += aux_1 * DN_DX_container[i_gauss](i,d) * rN(j);
                }
            }
        }

        // Add the volumetric strain mass stabilization term - term 2
        const double aux_2 = w_gauss * tau_1 * bulk_modulus;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                for (unsigned int d = 0; d < dim; ++d) {
                    rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim) += aux_2 * DN_DX_container[i_gauss](i, d) * DN_DX_container[i_gauss](j, d);
                }
            }
        }

        // Add the divergence mass stabilization term
        const double aux_3 = w_gauss * tau_2;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                for (unsigned int d = 0; d < dim; ++d) {
                    rLeftHandSideMatrix(i * block_size + dim, j * block_size + d) -= aux_3 * rN(i) * DN_DX_container[i_gauss](j,d);
                }
            }
        }

        // Add the volumetric strain stabilization term
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                rLeftHandSideMatrix(i * block_size + dim, j * block_size + dim) -= aux_3 * rN(i) * rN(j);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    const auto &r_geometry = GetGeometry();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Compute the geometry data
    Vector w_gauss_container;
    Matrix N_gauss_container;
    GeometryType::ShapeFunctionsGradientsType DN_DX_container;
    CalculateGeometryData(r_geometry, w_gauss_container, N_gauss_container, DN_DX_container);

    // Calculate the RHS contributions
    rRightHandSideVector.clear();

    Matrix B_mat(strain_size, n_nodes*dim);
    Matrix dev_strain_op;
    CalculateDeviatoricStrainOperator(dev_strain_op);

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Get Gauss pt. values
        const auto &rN = row(N_gauss_container, i_gauss);
        const double w_gauss = w_gauss_container[i_gauss];

        // Add the deviatoric strain contribution
        // TODO: MOVE TO A FUNCTION
        Vector tot_strain = ZeroVector(strain_size);
        CalculateB(B_mat, DN_DX_container[i_gauss]);
        for (unsigned int i = 0; i < strain_size; ++i) {
            for (unsigned int j = 0; j < strain_size ; ++j) {
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    const auto &r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int d = 0; d < dim; ++d) {
                        const unsigned int aux = i_node * dim + d;
                        tot_strain[i] += dev_strain_op(i,j) * B_mat(j,aux) * r_disp[d];
                    }
                }
            }
        }

        // Interpolate and add the nodal volumetric strain
        double vol_strain = 0.0;
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            vol_strain += rN[i_node] * r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }
        for (unsigned int d = 0; d < dim; ++d) {
            tot_strain[d] += vol_strain;
        }

        // Get the stress from the constitutive law
        // TODO: MOVE TO A FUNCTION
        Vector tot_stress(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto &r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        cons_law_values.SetShapeFunctionsValues(rN);
        cons_law_values.SetStrainVector(tot_strain);
        cons_law_values.SetStressVector(tot_stress);
        mConstitutiveLawVector[i_gauss]->CalculateMaterialResponseCauchy(cons_law_values);

        // Add momentum body force contribution
        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const double aux = rN[i] * rN[j];
                for (unsigned int d = 0; d < dim; ++d) {
                    rRightHandSideVector[i * block_size + d] += w_gauss * aux * body_force[d];
                }
            }
        }

        // Add momentum stress contribution
        // Note that this includes both the deviatoric and volumetric stress contributions
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d = 0; d < dim; ++d) {
                for (unsigned int j = 0; j < strain_size; ++j) {
                    for (unsigned int k = 0; k < strain_size; ++k) {
                        rRightHandSideVector(i * block_size + d) -= w_gauss * B_mat(j, i * dim + d) * dev_strain_op(k,j) * cons_law_values.GetStressVector()[k];
                    }
                }
            }
        }

        // Add mass conservation divergence contribution
        // TODO: These can be included in the body force contribution loop
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const auto &r_disp = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT);
                for (unsigned int d = 0; d < dim; ++d) {
                    rRightHandSideVector(i * block_size + dim) -= w_gauss * rN[i] * DN_DX_container[i_gauss](j,d) * r_disp[d];
                }
            }
        }

        // Add mass conservation volumetric strain contribution
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const double &r_vol_strain = r_geometry[j].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
                rRightHandSideVector(i * block_size + dim) -= w_gauss * rN[i] * rN[j] * r_vol_strain;
            }
        }

        // Calculate stabilization constants
        double hyd_stress = 0.0;
        for (unsigned int d = 0; d < dim; ++d) {
            hyd_stress += cons_law_values.GetStressVector()[d];
        }
        hyd_stress /= 3.0;
        const double bulk_modulus = vol_strain > 1.0e-12 ? hyd_stress * vol_strain / std::pow(vol_strain, 2) : 1.0e+12;
        const double shear_modulus = GetProperties()[YOUNG_MODULUS] / (2.0 * (1.0 + GetProperties()[POISSON_RATIO]));

        const double h = ComputeElementSize(DN_DX_container[i_gauss]);
        const double tau_1 = 2.0 * std::pow(h, 2) / (2.0 * shear_modulus);
        const double tau_2 = 0.15;

        // Add the volumetric strain momentum stabilization term - term 1
        const double aux_1 = w_gauss * bulk_modulus * tau_2;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d1 = 0; d1 < dim; ++d1) {
                for (unsigned int j = 0; j < n_nodes; ++j) {
                    const auto &r_disp = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT);
                    for (unsigned int d2 = 0; d2 < dim; ++d2) {
                        rRightHandSideVector(i * block_size + d1) += aux_1 * DN_DX_container[i_gauss](i,d1) * DN_DX_container[i_gauss](j,d2) * r_disp(d2);
                    }
                }
            }
        }

        // Add the volumetric strain momentum stabilization term - term 2
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d = 0; d < dim; ++d) {
                for (unsigned int j = 0; j < n_nodes; ++j) {
                    const double &r_vol_strain = r_geometry[j].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
                    rRightHandSideVector(i * block_size + d) -= aux_1 * DN_DX_container[i_gauss](i,d) * rN(j) * r_vol_strain;
                }
            }
        }

        // Add the volumetric strain mass stabilization term - term 1
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int d = 0; d < dim; ++d) {
                rRightHandSideVector(i * block_size + dim) += w_gauss * tau_1 * DN_DX_container[i_gauss](i,d) * body_force[d];
            }
        }

        // Add the volumetric strain mass stabilization term - term 2
        const double aux_2 = w_gauss * tau_1 * bulk_modulus;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const double &r_vol_strain = r_geometry[j].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
                for (unsigned int d = 0; d < dim; ++d) {
                    rRightHandSideVector(i * block_size + dim) += aux_2 * DN_DX_container[i_gauss](i, d) * DN_DX_container[i_gauss](j, d) * r_vol_strain;
                }
            }
        }

        // Add the divergence mass stabilization term
        const double aux_3 = w_gauss * tau_2;
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const auto &r_disp = r_geometry[j].FastGetSolutionStepValue(DISPLACEMENT);
                for (unsigned int d = 0; d < dim; ++d) {
                    rRightHandSideVector(i * block_size + dim) += aux_3 * rN(i) * DN_DX_container[i_gauss](j,d) * r_disp(d);
                }
            }
        }

        // Add the volumetric strain stabilization term
        for (unsigned int i = 0; i < n_nodes; ++i) {
            for (unsigned int j = 0; j < n_nodes; ++j) {
                const double &r_vol_strain = r_geometry[j].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
                rRightHandSideVector(i * block_size + dim) += aux_3 * rN(i) * rN(j) * r_vol_strain;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto &r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto &r_geometry = GetGeometry();
        const auto &r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector) {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedStrainElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight = 0.0;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        // CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        // CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        // int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        // if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        //     // Contributions to stiffness matrix calculated on the reference config
        //     this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        // }

        // if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        //     this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
        // }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> SmallDisplacementMixedStrainElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i) {
        body_force[i] = 0.0;
    }

    const auto &r_properties = GetProperties();
    const double density = r_properties.Has(DENSITY) ? r_properties[DENSITY] : 0.0;

    if (r_properties.Has(VOLUME_ACCELERATION)) {
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];
    }

    const auto &r_geometry = GetGeometry();
    if(r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION)) {
        Vector N;
        N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
        }
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateGeometryData(
    const GeometryType &rGeometry,
    Vector &rWeightsContainer,
    Matrix &rShapeFunctionsContainer,
    GeometryType::ShapeFunctionsGradientsType &rDNDXContainer) const
{
    // Compute the geometry data
    const SizeType n_nodes = rGeometry.PointsNumber();
    const auto r_integration_method = GetIntegrationMethod();
    const SizeType n_gauss = rGeometry.IntegrationPointsNumber(r_integration_method);
    const auto& r_integration_points = rGeometry.IntegrationPoints(r_integration_method);

    // Calculate the shape function values
    if (rShapeFunctionsContainer.size1() != n_gauss || rShapeFunctionsContainer.size2() != n_nodes) {
        rShapeFunctionsContainer.resize(n_gauss, n_nodes,false);
    }
    rShapeFunctionsContainer = rGeometry.ShapeFunctionsValues(r_integration_method);

    // Calculate the shape function gradients values
    Vector detJ_container;
    rGeometry.ShapeFunctionsIntegrationPointsGradients(rDNDXContainer, detJ_container, r_integration_method);

    // Calculate the Gauss point weights (already multiplied by det(J))
    if (rWeightsContainer.size() != n_gauss) {
        rWeightsContainer.resize(n_gauss, false);
    }
    for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        rWeightsContainer[i_gauss] = detJ_container[i_gauss] * r_integration_points[i_gauss].Weight();
    }

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*2    ) = rDN_DX(i, 0);
            rB(1, i*2 + 1) = rDN_DX(i, 1);
            rB(2, i*2    ) = rDN_DX(i, 1);
            rB(2, i*2 + 1) = rDN_DX(i, 0);
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*3    ) = rDN_DX(i, 0);
            rB(1, i*3 + 1) = rDN_DX(i, 1);
            rB(2, i*3 + 2) = rDN_DX(i, 2);
            rB(3, i*3    ) = rDN_DX(i, 1);
            rB(3, i*3 + 1) = rDN_DX(i, 0);
            rB(4, i*3 + 1) = rDN_DX(i, 2);
            rB(4, i*3 + 2) = rDN_DX(i, 1);
            rB(5, i*3    ) = rDN_DX(i, 2);
            rB(5, i*3 + 2) = rDN_DX(i, 0);
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::CalculateDeviatoricStrainOperator(Matrix& rDevStrainOp) const
{
    KRATOS_TRY;

    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    rDevStrainOp = ZeroMatrix(strain_size, strain_size);

    const double two_thirds = 2.0/3.0;
    const double minus_one_third = -1.0/3.0;

    if (strain_size == 3) {
        rDevStrainOp(0,0) = two_thirds;
        rDevStrainOp(0,1) = minus_one_third;
        rDevStrainOp(1,0) = minus_one_third;
        rDevStrainOp(1,1) = two_thirds;
        rDevStrainOp(2,2) = 1.0;
    } else if (strain_size == 6) {
        rDevStrainOp(0,0) = two_thirds;
        rDevStrainOp(0,1) = minus_one_third;
        rDevStrainOp(0,2) = minus_one_third;
        rDevStrainOp(1,0) = minus_one_third;
        rDevStrainOp(1,1) = two_thirds;
        rDevStrainOp(1,2) = minus_one_third;
        rDevStrainOp(2,0) = minus_one_third;
        rDevStrainOp(2,1) = minus_one_third;
        rDevStrainOp(2,2) = two_thirds;
        rDevStrainOp(3,3) = 1.0;
        rDevStrainOp(4,4) = 1.0;
        rDevStrainOp(5,5) = 1.0;
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

double SmallDisplacementMixedStrainElement::ComputeElementSize(const Matrix &rDN_DX) const
{
    double h = 0.0;
    for (unsigned int i_node = 0; i_node < GetGeometry().PointsNumber(); ++i_node) {
        double h_inv = 0.0;
        for (unsigned int d = 0; d < GetGeometry().WorkingSpaceDimension(); ++d) {
            h_inv += rDN_DX(i_node, d) * rDN_DX(i_node, d);
        }
        h += 1.0 / h_inv;
    }
    h = sqrt(h) / static_cast<double>(GetGeometry().PointsNumber());
    return h;
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedStrainElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int base_element_check = SmallDisplacementMixedStrainElement::BaseType::Check(rCurrentProcessInfo);

    return base_element_check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainElement::BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainElement::BaseType);
}

} // Namespace Kratos
