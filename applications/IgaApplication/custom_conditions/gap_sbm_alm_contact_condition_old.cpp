// //  KRATOS  _____________
// //         /  _/ ____/   |
// //         / // / __/ /| |
// //       _/ // /_/ / ___ |
// //      /___/\____/_/  |_| Application
// //
// //  License:         BSD License
// //                   Kratos default license: kratos/license.txt
// //
// //  Main authors:    Andrea Gorgi
// //

// // Project includes
// #include "custom_conditions/gap_sbm_alm_contact_condition.h"

// namespace Kratos
// {

// void GapSbmALMContactCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
// {
//     BaseType::Initialize(rCurrentProcessInfo);
// }

// void GapSbmALMContactCondition::CalculateRightHandSide(
//     VectorType& rRightHandSideVector,
//     const ProcessInfo& rCurrentProcessInfo)
// {
//     MatrixType left_hand_side_matrix;
//     CalculateAll(
//         left_hand_side_matrix,
//         rRightHandSideVector,
//         rCurrentProcessInfo,
//         false,
//         true);
// }

// void GapSbmALMContactCondition::CalculateLeftHandSide(
//     MatrixType& rLeftHandSideMatrix,
//     const ProcessInfo& rCurrentProcessInfo)
// {
//     VectorType right_hand_side_vector;
//     CalculateAll(
//         rLeftHandSideMatrix,
//         right_hand_side_vector,
//         rCurrentProcessInfo,
//         true,
//         false);
// }

// void GapSbmALMContactCondition::CalculateLocalSystem(
//     MatrixType& rLeftHandSideMatrix,
//     VectorType& rRightHandSideVector,
//     const ProcessInfo& rCurrentProcessInfo)
// {
//     CalculateAll(
//         rLeftHandSideMatrix,
//         rRightHandSideVector,
//         rCurrentProcessInfo,
//         true,
//         true);
// }

// void GapSbmALMContactCondition::CalculateAll(
//     MatrixType& rLeftHandSideMatrix,
//     VectorType& rRightHandSideVector,
//     const ProcessInfo& rCurrentProcessInfo,
//     const bool CalculateStiffnessMatrixFlag,
//     const bool CalculateResidualVectorFlag)
// {
//     KRATOS_TRY

//     mSideIdentifier = GetValue(IDENTIFIER);

//     if (mSideIdentifier == "SLAVE" || mSideIdentifier == "INACTIVE") {
//         if (CalculateResidualVectorFlag) {
//             if (rRightHandSideVector.size() != 0) {
//                 rRightHandSideVector.resize(0, false);
//             }
//         }
//         if (CalculateStiffnessMatrixFlag) {
//             if (rLeftHandSideMatrix.size1() != 0 || rLeftHandSideMatrix.size2() != 0) {
//                 rLeftHandSideMatrix.resize(0, 0, false);
//             }
//         }
//         return;
//     }

//     const auto& r_master_surrogate_geometry = GetMasterSurrogateGeometry();
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];

//     const SizeType number_of_control_points_master = r_master_surrogate_geometry.size();
//     const SizeType number_of_control_points_slave = r_slave_surrogate_geometry.size();

//     const SizeType disp_size_master = mDim * number_of_control_points_master;
//     const SizeType disp_size_slave = mDim * number_of_control_points_slave;
//     const SizeType lm_size_slave = number_of_control_points_slave;
//     const SizeType mat_size = disp_size_master + disp_size_slave + lm_size_slave;

//     if (CalculateStiffnessMatrixFlag) {
//         if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size) {
//             rLeftHandSideMatrix.resize(mat_size, mat_size, false);
//         }
//         noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
//     }

//     if (CalculateResidualVectorFlag) {
//         if (rRightHandSideVector.size() != mat_size) {
//             rRightHandSideVector.resize(mat_size, false);
//         }
//         noalias(rRightHandSideVector) = ZeroVector(mat_size);
//     }

//     ContactKinematics kinematics;
//     EvaluateContactKinematics(kinematics);
//     SetValue(NORMAL_GAP, kinematics.NormalGap);

//     const double integration_weight = GetValue(INTEGRATION_WEIGHT);
//     const double penalty = CalculateALMPenalty(rCurrentProcessInfo);
//     const double scale_factor = GetScaleFactor(rCurrentProcessInfo);
//     const double contact_traction = kinematics.Lambda + penalty * kinematics.NormalGap;
//     const bool is_active = GetValue(ACTIVATION_LEVEL) == 3;

//     const IndexType slave_disp_offset = disp_size_master;
//     const IndexType lm_offset = disp_size_master + disp_size_slave;

//     if (is_active) {
//         for (IndexType i = 0; i < number_of_control_points_master; ++i) {
//             for (IndexType idim = 0; idim < mDim; ++idim) {
//                 const double g_i = kinematics.MasterShapeFunctions(i) * kinematics.Normal[idim];
//                 const IndexType row = mDim * i + idim;

//                 if (CalculateResidualVectorFlag) {
//                     rRightHandSideVector(row) -= g_i * contact_traction * integration_weight;
//                 }

//                 if (CalculateStiffnessMatrixFlag) {
//                     for (IndexType j = 0; j < number_of_control_points_master; ++j) {
//                         for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                             const double g_j = kinematics.MasterShapeFunctions(j) * kinematics.Normal[jdim];
//                             const IndexType column = mDim * j + jdim;
//                             rLeftHandSideMatrix(row, column) += penalty * g_i * g_j * integration_weight;
//                         }
//                     }

//                     for (IndexType j = 0; j < number_of_control_points_slave; ++j) {
//                         for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                             const double g_j = -kinematics.SlaveShapeFunctions(j) * kinematics.Normal[jdim];
//                             const IndexType column = slave_disp_offset + mDim * j + jdim;
//                             rLeftHandSideMatrix(row, column) += penalty * g_i * g_j * integration_weight;
//                         }

//                         const IndexType lambda_column = lm_offset + j;
//                         rLeftHandSideMatrix(row, lambda_column) +=
//                             g_i * kinematics.SlaveShapeFunctions(j) * integration_weight;
//                     }
//                 }
//             }
//         }

//         for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//             for (IndexType idim = 0; idim < mDim; ++idim) {
//                 const double g_i = -kinematics.SlaveShapeFunctions(i) * kinematics.Normal[idim];
//                 const IndexType row = slave_disp_offset + mDim * i + idim;

//                 if (CalculateResidualVectorFlag) {
//                     rRightHandSideVector(row) -= g_i * contact_traction * integration_weight;
//                 }

//                 if (CalculateStiffnessMatrixFlag) {
//                     for (IndexType j = 0; j < number_of_control_points_master; ++j) {
//                         for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                             const double g_j = kinematics.MasterShapeFunctions(j) * kinematics.Normal[jdim];
//                             const IndexType column = mDim * j + jdim;
//                             rLeftHandSideMatrix(row, column) += penalty * g_i * g_j * integration_weight;
//                         }
//                     }

//                     for (IndexType j = 0; j < number_of_control_points_slave; ++j) {
//                         for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                             const double g_j = -kinematics.SlaveShapeFunctions(j) * kinematics.Normal[jdim];
//                             const IndexType column = slave_disp_offset + mDim * j + jdim;
//                             rLeftHandSideMatrix(row, column) += penalty * g_i * g_j * integration_weight;
//                         }

//                         const IndexType lambda_column = lm_offset + j;
//                         rLeftHandSideMatrix(row, lambda_column) +=
//                             g_i * kinematics.SlaveShapeFunctions(j) * integration_weight;
//                     }
//                 }
//             }
//         }

//         for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//             const IndexType row = lm_offset + i;

//             if (CalculateResidualVectorFlag) {
//                 rRightHandSideVector(row) -=
//                     kinematics.SlaveShapeFunctions(i) * kinematics.NormalGap * integration_weight;
//             }

//             if (CalculateStiffnessMatrixFlag) {
//                 for (IndexType j = 0; j < number_of_control_points_master; ++j) {
//                     for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                         const double g_j = kinematics.MasterShapeFunctions(j) * kinematics.Normal[jdim];
//                         const IndexType column = mDim * j + jdim;
//                         rLeftHandSideMatrix(row, column) +=
//                             kinematics.SlaveShapeFunctions(i) * g_j * integration_weight;
//                     }
//                 }

//                 for (IndexType j = 0; j < number_of_control_points_slave; ++j) {
//                     for (IndexType jdim = 0; jdim < mDim; ++jdim) {
//                         const double g_j = -kinematics.SlaveShapeFunctions(j) * kinematics.Normal[jdim];
//                         const IndexType column = slave_disp_offset + mDim * j + jdim;
//                         rLeftHandSideMatrix(row, column) +=
//                             kinematics.SlaveShapeFunctions(i) * g_j * integration_weight;
//                     }
//                 }
//             }
//         }
//     } else {
//         for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//             const IndexType row = lm_offset + i;

//             if (CalculateResidualVectorFlag) {
//                 rRightHandSideVector(row) -=
//                     scale_factor * kinematics.SlaveShapeFunctions(i) * kinematics.Lambda * integration_weight;
//             }

//             if (CalculateStiffnessMatrixFlag) {
//                 for (IndexType j = 0; j < number_of_control_points_slave; ++j) {
//                     const IndexType column = lm_offset + j;
//                     rLeftHandSideMatrix(row, column) +=
//                         scale_factor * kinematics.SlaveShapeFunctions(i) *
//                         kinematics.SlaveShapeFunctions(j) * integration_weight;
//                 }
//             }
//         }
//     }

//     KRATOS_CATCH("")
// }

// void GapSbmALMContactCondition::EvaluateContactKinematics(ContactKinematics& rKinematics) const
// {
//     const auto& r_master_surrogate_geometry = GetMasterSurrogateGeometry();
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];

//     ComputeTaylorExpansionContribution(
//         r_master_surrogate_geometry,
//         mDistanceVectorSkinReferenceMaster,
//         mBasisFunctionsOrderMaster,
//         rKinematics.MasterShapeFunctions);

//     ComputeTaylorExpansionContribution(
//         r_slave_surrogate_geometry,
//         mDistanceVectorSkinReferenceSlave,
//         mBasisFunctionsOrderSlave,
//         rKinematics.SlaveShapeFunctions);

//     GetSolutionCoefficientVector(rKinematics.MasterDisplacementCoefficients, 0);
//     GetSolutionCoefficientVector(rKinematics.SlaveDisplacementCoefficients, 1);
//     GetLambdaCoefficientVector(rKinematics.LambdaCoefficients);

//     rKinematics.Normal = mNormalPhysicalSpaceMaster;

//     const array_1d<double, 3> master_center = GetGeometry().Center().Coordinates();
//     const array_1d<double, 3> slave_center = pGetProjectionNode()->Coordinates();
//     rKinematics.InitialGap =
//         (master_center[0] - slave_center[0]) * rKinematics.Normal[0] +
//         (master_center[1] - slave_center[1]) * rKinematics.Normal[1];

//     double master_displacement_normal = 0.0;
//     for (IndexType i = 0; i < r_master_surrogate_geometry.size(); ++i) {
//         master_displacement_normal +=
//             rKinematics.MasterShapeFunctions(i) *
//             (rKinematics.MasterDisplacementCoefficients[mDim * i] * rKinematics.Normal[0] +
//              rKinematics.MasterDisplacementCoefficients[mDim * i + 1] * rKinematics.Normal[1]);
//     }

//     double slave_displacement_normal = 0.0;
//     double lambda = 0.0;
//     for (IndexType i = 0; i < r_slave_surrogate_geometry.size(); ++i) {
//         slave_displacement_normal +=
//             rKinematics.SlaveShapeFunctions(i) *
//             (rKinematics.SlaveDisplacementCoefficients[mDim * i] * rKinematics.Normal[0] +
//              rKinematics.SlaveDisplacementCoefficients[mDim * i + 1] * rKinematics.Normal[1]);
//         lambda += rKinematics.SlaveShapeFunctions(i) * rKinematics.LambdaCoefficients[i];
//     }

//     rKinematics.NormalGap =
//         rKinematics.InitialGap + master_displacement_normal - slave_displacement_normal;
//     rKinematics.Lambda = lambda;
// }

// void GapSbmALMContactCondition::GetLambdaCoefficientVector(Vector& rValues) const
// {
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];
//     const SizeType number_of_control_points_slave = r_slave_surrogate_geometry.size();

//     if (rValues.size() != number_of_control_points_slave) {
//         rValues.resize(number_of_control_points_slave, false);
//     }

//     for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//         rValues[i] = r_slave_surrogate_geometry[i].FastGetSolutionStepValue(SCALAR_LAGRANGE_MULTIPLIER);
//     }
// }

// void GapSbmALMContactCondition::EquationIdVector(
//     EquationIdVectorType& rResult,
//     const ProcessInfo& rCurrentProcessInfo) const
// {
//     const std::string& r_side_identifier = GetValue(IDENTIFIER);
//     if (r_side_identifier == "SLAVE" || r_side_identifier == "INACTIVE") {
//         if (rResult.size() != 0) {
//             rResult.resize(0, false);
//         }
//         return;
//     }

//     const auto& r_master_surrogate_geometry = GetMasterSurrogateGeometry();
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];

//     const SizeType number_of_control_points_master = r_master_surrogate_geometry.size();
//     const SizeType number_of_control_points_slave = r_slave_surrogate_geometry.size();

//     const SizeType mat_size = mDim * (number_of_control_points_master + number_of_control_points_slave)
//         + number_of_control_points_slave;

//     if (rResult.size() != mat_size) {
//         rResult.resize(mat_size, false);
//     }

//     IndexType index = 0;
//     for (IndexType i = 0; i < number_of_control_points_master; ++i) {
//         const auto& r_node = r_master_surrogate_geometry[i];
//         rResult[index++] = r_node.GetDof(DISPLACEMENT_X).EquationId();
//         rResult[index++] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
//     }

//     for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//         const auto& r_node = r_slave_surrogate_geometry[i];
//         rResult[index++] = r_node.GetDof(DISPLACEMENT_X).EquationId();
//         rResult[index++] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
//     }

//     for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//         rResult[index++] = r_slave_surrogate_geometry[i].GetDof(SCALAR_LAGRANGE_MULTIPLIER).EquationId();
//     }
// }

// void GapSbmALMContactCondition::GetDofList(
//     DofsVectorType& rElementalDofList,
//     const ProcessInfo& rCurrentProcessInfo) const
// {
//     const std::string& r_side_identifier = GetValue(IDENTIFIER);
//     if (r_side_identifier == "SLAVE" || r_side_identifier == "INACTIVE") {
//         if (rElementalDofList.size() != 0) {
//             rElementalDofList.resize(0);
//         }
//         return;
//     }

//     const auto& r_master_surrogate_geometry = GetMasterSurrogateGeometry();
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];

//     const SizeType number_of_control_points_master = r_master_surrogate_geometry.size();
//     const SizeType number_of_control_points_slave = r_slave_surrogate_geometry.size();

//     rElementalDofList.resize(0);
//     rElementalDofList.reserve(mDim * (number_of_control_points_master + number_of_control_points_slave)
//         + number_of_control_points_slave);

//     for (IndexType i = 0; i < number_of_control_points_master; ++i) {
//         const auto& r_node = r_master_surrogate_geometry[i];
//         rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
//         rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
//     }

//     for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//         const auto& r_node = r_slave_surrogate_geometry[i];
//         rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
//         rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
//     }

//     for (IndexType i = 0; i < number_of_control_points_slave; ++i) {
//         rElementalDofList.push_back(r_slave_surrogate_geometry[i].pGetDof(SCALAR_LAGRANGE_MULTIPLIER));
//     }
// }

// int GapSbmALMContactCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
// {
//     KRATOS_TRY

//     const std::string& r_side_identifier = GetValue(IDENTIFIER);
//     if (r_side_identifier == "SLAVE" || r_side_identifier == "INACTIVE") {
//         return 0;
//     }

//     const auto& r_master_surrogate_geometry = GetMasterSurrogateGeometry();
//     const auto& r_slave_surrogate_geometry = *pGetProjectionNode()->GetValue(NEIGHBOUR_GEOMETRIES)[0];

//     for (const auto& r_node : r_master_surrogate_geometry) {
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
//         KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
//         KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
//     }

//     for (const auto& r_node : r_slave_surrogate_geometry) {
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(SCALAR_LAGRANGE_MULTIPLIER, r_node);
//         KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
//         KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
//         KRATOS_CHECK_DOF_IN_NODE(SCALAR_LAGRANGE_MULTIPLIER, r_node);
//     }

//     KRATOS_ERROR_IF_NOT(mDim == 2)
//         << "\"GapSbmALMContactCondition\" #" << Id()
//         << " supports only 2D conditions, but the current dimension is " << mDim << std::endl;

//     return 0;

//     KRATOS_CATCH("")

//     return 0;
// }

// double GapSbmALMContactCondition::CalculateCurrentNormalGap() const
// {
//     ContactKinematics kinematics;
//     EvaluateContactKinematics(kinematics);
//     return kinematics.NormalGap;
// }

// double GapSbmALMContactCondition::CalculateCurrentLambda() const
// {
//     ContactKinematics kinematics;
//     EvaluateContactKinematics(kinematics);
//     return kinematics.Lambda;
// }

// double GapSbmALMContactCondition::CalculateALMPenalty(const ProcessInfo& rCurrentProcessInfo) const
// {
//     double reference_penalty = 0.0;
//     if (rCurrentProcessInfo.Has(INITIAL_PENALTY)) {
//         reference_penalty = rCurrentProcessInfo[INITIAL_PENALTY];
//     } else if (GetProperties().Has(PENALTY_FACTOR)) {
//         reference_penalty = GetProperties()[PENALTY_FACTOR];
//     } else {
//         reference_penalty = mPenalty;
//     }

//     if (reference_penalty <= 0.0) {
//         reference_penalty = 1.0;
//     }

//     double characteristic_length = 0.0;
//     if (this->Has(CHARACTERISTIC_GEOMETRY_LENGTH)) {
//         const Vector& r_characteristic_length = this->GetValue(CHARACTERISTIC_GEOMETRY_LENGTH);
//         if (r_characteristic_length.size() > 0) {
//             characteristic_length = r_characteristic_length[0];
//         }
//     }

//     if (characteristic_length <= 0.0 && this->Has(KNOT_SPAN_SIZES)) {
//         const Vector& r_mesh_size_uv = this->GetValue(KNOT_SPAN_SIZES);
//         if (r_mesh_size_uv.size() > 0) {
//             characteristic_length = r_mesh_size_uv[0];
//             for (IndexType i = 1; i < r_mesh_size_uv.size(); ++i) {
//                 characteristic_length = std::min(characteristic_length, r_mesh_size_uv[i]);
//             }
//         }
//     }

//     if (characteristic_length <= 0.0) {
//         characteristic_length = 1.0;
//     }

//     return reference_penalty / characteristic_length *
//         static_cast<double>(mBasisFunctionsOrderMaster * mBasisFunctionsOrderMaster) / 4.0;
// }

// double GapSbmALMContactCondition::CalculateAugmentedTraction(const ProcessInfo& rCurrentProcessInfo) const
// {
//     ContactKinematics kinematics;
//     EvaluateContactKinematics(kinematics);
//     return GetScaleFactor(rCurrentProcessInfo) * kinematics.Lambda
//         + CalculateALMPenalty(rCurrentProcessInfo) * kinematics.NormalGap;
// }

// double GapSbmALMContactCondition::GetScaleFactor(const ProcessInfo& rCurrentProcessInfo) const
// {
//     return rCurrentProcessInfo.Has(SCALE_FACTOR) ? rCurrentProcessInfo[SCALE_FACTOR] : 1.0;
// }

// } // namespace Kratos
