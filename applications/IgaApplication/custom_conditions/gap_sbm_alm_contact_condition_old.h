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

// #pragma once

// // Project includes
// #include "custom_conditions/gap_sbm_contact_condition.h"

// namespace Kratos
// {

// class KRATOS_API(IGA_APPLICATION) GapSbmALMContactCondition
//     : public GapSbmContactCondition
// {
// public:
//     KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmALMContactCondition);

//     using BaseType = GapSbmContactCondition;
//     using IndexType = BaseType::IndexType;
//     using SizeType = BaseType::SizeType;

//     GapSbmALMContactCondition(
//         IndexType NewId,
//         GeometryType::Pointer pGeometry)
//         : BaseType(NewId, pGeometry)
//     {
//     }

//     GapSbmALMContactCondition(
//         IndexType NewId,
//         GeometryType::Pointer pGeometry,
//         PropertiesType::Pointer pProperties)
//         : BaseType(NewId, pGeometry, pProperties)
//     {
//     }

//     GapSbmALMContactCondition()
//         : BaseType()
//     {
//     }

//     ~GapSbmALMContactCondition() override = default;

//     void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

//     Condition::Pointer Create(
//         IndexType NewId,
//         GeometryType::Pointer pGeom,
//         PropertiesType::Pointer pProperties) const override
//     {
//         return Kratos::make_intrusive<GapSbmALMContactCondition>(NewId, pGeom, pProperties);
//     }

//     Condition::Pointer Create(
//         IndexType NewId,
//         NodesArrayType const& ThisNodes,
//         PropertiesType::Pointer pProperties) const override
//     {
//         return Kratos::make_intrusive<GapSbmALMContactCondition>(
//             NewId,
//             GetGeometry().Create(ThisNodes),
//             pProperties);
//     }

//     void CalculateRightHandSide(
//         VectorType& rRightHandSideVector,
//         const ProcessInfo& rCurrentProcessInfo) override;

//     void CalculateLeftHandSide(
//         MatrixType& rLeftHandSideMatrix,
//         const ProcessInfo& rCurrentProcessInfo) override;

//     void CalculateLocalSystem(
//         MatrixType& rLeftHandSideMatrix,
//         VectorType& rRightHandSideVector,
//         const ProcessInfo& rCurrentProcessInfo) override;

//     void EquationIdVector(
//         EquationIdVectorType& rResult,
//         const ProcessInfo& rCurrentProcessInfo) const override;

//     void GetDofList(
//         DofsVectorType& rElementalDofList,
//         const ProcessInfo& rCurrentProcessInfo) const override;

//     int Check(const ProcessInfo& rCurrentProcessInfo) const override;

//     double CalculateCurrentNormalGap() const;

//     double CalculateCurrentLambda() const;

//     double CalculateALMPenalty(const ProcessInfo& rCurrentProcessInfo) const;

//     double CalculateAugmentedTraction(const ProcessInfo& rCurrentProcessInfo) const;

//     std::string Info() const override
//     {
//         std::stringstream buffer;
//         buffer << "\"GapSbmALMContactCondition\" #" << Id();
//         return buffer.str();
//     }

//     void PrintInfo(std::ostream& rOStream) const override
//     {
//         rOStream << "\"GapSbmALMContactCondition\" #" << Id();
//     }

// protected:
//     struct ContactKinematics
//     {
//         Vector MasterShapeFunctions;
//         Vector SlaveShapeFunctions;
//         Vector MasterDisplacementCoefficients;
//         Vector SlaveDisplacementCoefficients;
//         Vector LambdaCoefficients;
//         array_1d<double, 3> Normal = {0.0, 0.0, 0.0};
//         double InitialGap = 0.0;
//         double NormalGap = 0.0;
//         double Lambda = 0.0;
//     };

//     void CalculateAll(
//         MatrixType& rLeftHandSideMatrix,
//         VectorType& rRightHandSideVector,
//         const ProcessInfo& rCurrentProcessInfo,
//         const bool CalculateStiffnessMatrixFlag,
//         const bool CalculateResidualVectorFlag);

//     void EvaluateContactKinematics(ContactKinematics& rKinematics) const;

//     void GetLambdaCoefficientVector(Vector& rValues) const;

//     double GetScaleFactor(const ProcessInfo& rCurrentProcessInfo) const;

// private:
//     const GeometryType& GetMasterSurrogateGeometry() const
//     {
//         return *this->GetValue(NEIGHBOUR_GEOMETRIES)[0];
//     }
// };

// } // namespace Kratos
