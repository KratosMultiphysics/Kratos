//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once


// System includes
#include "includes/define.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "includes/variables.h"


// External includes

// Project includes
#include "iga_application_variables.h"

// Project includes
#include "includes/constitutive_law.h"

#include "geometries/quadrature_point_coupling_geometry_2d.h"

namespace Kratos
{
    /// Condition for penalty support condition
    class SbmContact2DCondition
        : public Condition
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Counted pointer definition of SbmContact2DCondition
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SbmContact2DCondition);

        /// Size types
        typedef std::size_t SizeType;
        typedef std::size_t IndexType;
        using NodeType = typename Condition::NodeType;

        ///@}
        ///@name Life Cycle
        ///@{
        
        void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

        /// Constructor with Id and geometry
        SbmContact2DCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry)
            : Condition(NewId, pGeometry)
        {};

        /// Constructor with Id, geometry and property
        SbmContact2DCondition(
            IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pPropMaster, 
            PropertiesType::Pointer pPropSlave)
            : Condition(NewId, pGeometry),
              mpPropMaster(pPropMaster),
              mpPropSlave(pPropSlave)
        {
            const int integrationDomain = GetGeometry().GetValue(ACTIVATION_LEVEL);

            if (integrationDomain == 1) {
                mpPropMaster = pPropSlave;
                mpPropSlave = pPropMaster;
            }
        };

        /// Default constructor
        SbmContact2DCondition() : Condition()
        {};

        /// Destructor
        // ~SbmContact2DCondition() override
        // {
        // };

        ///@}
        ///@name Life Cycle
        ///@{

        /// Create with Id, pointer to geometry and pointer to property
        Condition::Pointer Create(
            IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pPropMaster,
            PropertiesType::Pointer pPropSlave
        ) const
        {
            return Kratos::make_intrusive<SbmContact2DCondition>(
                NewId, pGeom, pPropMaster, pPropSlave);
        };

        GeometryType & GetSlaveGeometry() 
        {
            return this->GetGeometry().GetGeometryPart(QuadraturePointCouplingGeometry2D<Point>::Slave);
        };

        GeometryType const& GetSlaveGeometry() const
        {
            return this->GetGeometry().GetGeometryPart(QuadraturePointCouplingGeometry2D<Point>::Slave);
        };


        /**
         * @brief This method returns the paired geometry (constant version)
         * @return The master geometry (master in the definition of Popp which is the opposite of the standard)
         */
        GeometryType & GetMasterGeometry() 
        {
            return this->GetGeometry().GetGeometryPart(QuadraturePointCouplingGeometry2D<Point>::Master);
        };


        GeometryType const& GetMasterGeometry() const
        {
            return this->GetGeometry().GetGeometryPart(QuadraturePointCouplingGeometry2D<Point>::Master);
        };

        /// Create with Id, pointer to geometry and pointer to property
        // Condition::Pointer Create(
        //     IndexType NewId,
        //     NodesArrayType const& ThisNodes,
        //     PropertiesType::Pointer pProperties
        // ) const override
        // {
        //     return Kratos::make_intrusive<SbmContact2DCondition>(
        //         NewId, GetGeometry().Create(ThisNodes), pProperties);
        // };

        ///@}
        ///@name Operations
        ///@{

        /**
        * @brief This is called during the assembling process in order
        *        to calculate the condition right hand side matrix
        * @param rLeftHandSideMatrix the condition right hand side matrix
        * @param rCurrentProcessInfo the current process info
        */
        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

        /**
        * @brief This is called during the assembling process in order
        *        to calculate the condition left hand side matrix
        * @param rLeftHandSideMatrix the condition left hand side matrix
        * @param rCurrentProcessInfo the current process info
        */
        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo) override;

        /**
         * @brief This function provides a more general interface to the element.
         * @details It is designed so that rLHSvariables and rRHSvariables are
         *          passed to the element thus telling what is the desired output
         * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
         * @param rRightHandSideVector container for the desired RHS output
         * @param rCurrentProcessInfo the current process info instance
         */
        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

        /**
        * @brief Sets on rResult the ID's of the element degrees of freedom
        * @param rResult The vector containing the equation id
        * @param rCurrentProcessInfo The current process info instance
        */
        void EquationIdVector(
            EquationIdVectorType& rResult,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        /**
        * @brief Sets on rConditionDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the dof of the element
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
        ) const override;

        void GetValuesVector(Vector& rValues, IndexType index) const;

        void GetStrainVector(Vector& strainVector, IndexType index) const;

        void SetConstitutiveVariables(Vector& StrainVector, IndexType index, const Kratos::ProcessInfo& rCurrentProcessInfo); 

        void SetGap(); 

        void CalculateB(
            Matrix& rB, 
            Matrix& r_DN_DX,
            const SizeType number_of_control_points) const;
        ///@}
        ///@name Check
        ///@{

        /// Performs check if Penalty factor is provided.
        int Check(const ProcessInfo& rCurrentProcessInfo) const override;

        void CalculateOnIntegrationPoints(
            const Variable<double> &rVariable,
            std::vector<double> &rOutput,
            const ProcessInfo &rCurrentProcessInfo) override;
        
        void CalculateOnIntegrationPoints(
            const Variable<Vector>& rVariable,
            std::vector<Vector>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "\"SbmContact2DCondition\" #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "\"SbmContact2DCondition\" #" << Id();
        }

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override
        {
            pGetGeometry()->PrintData(rOStream);
        }

        ///@}

    protected:

    
    /**
     * Internal variables used in the constitutive calculations
     */
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);

            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);

            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };

    ///@name Protected static Member Variables
    ///@{
    void InitializeMaterial();
    void InitializeMemberVariables();
    void InitializeSbmMemberVariables();

    void ComputeGradientTaylorExpansionContribution(
        const GeometryType& rGeometry,
        const Vector& rDistanceVector,
        const SizeType BasisFunctionsOrder,
        Matrix& rGradHSum) const;

        void ComputeTaylorSumContribution(
            const GeometryType& rGeometry,
            const Vector& rDistanceVector,
            const SizeType BasisFunctionsOrder,
            Matrix& rHSum) const;

        void ApplyConstitutiveLaw(
            ConstitutiveLaw::Pointer pConstitutiveLaw,
            Vector& rStrain,
            ConstitutiveLaw::Parameters& rValues,
            ConstitutiveVariables& rConstitutiveVariables) const;

    PropertiesType& GetPropertiesMaster()
    {
        KRATOS_DEBUG_ERROR_IF(mpPropMaster == nullptr)
            << "Tryining to get the master properties of " << Info()
            << ", which are uninitialized." << std::endl;
        return *mpPropMaster;
    }

    PropertiesType& GetPropertiesSlave()
    {
        KRATOS_DEBUG_ERROR_IF(mpPropSlave == nullptr)
            << "Tryining to get the slave properties of " << Info()
            << ", which are uninitialized." << std::endl;
        return *mpPropSlave;
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    //@}
    ///@name Protected member Variables
    ///@{
    ConstitutiveLaw::Pointer mpConstitutiveLawSlave; /// The pointer containing the constitutive laws
    ConstitutiveLaw::Pointer mpConstitutiveLawMaster; /// The pointer containing the constitutive laws

    PropertiesType::Pointer mpPropMaster;
    PropertiesType::Pointer mpPropSlave;

    bool m_contact_is_active = false;

    SizeType mMasterDim = 0;
    SizeType mSlaveDim = 0;

    SizeType mMasterBasisFunctionsOrder = 0;
    SizeType mSlaveBasisFunctionsOrder = 0;

    double mMasterCharacteristicLength = 0.0;
    double mSlaveCharacteristicLength = 0.0;

    NodeType* mpProjectionNodeMaster = nullptr;
    NodeType* mpProjectionNodeSlave = nullptr;

    Vector mDistanceMaster = ZeroVector(2);
    Vector mDistanceSlave = ZeroVector(2);

    Vector mNormalMaster = ZeroVector(3);
    Vector mNormalSlave = ZeroVector(3);

    Vector mTrueNormalMaster = ZeroVector(3);
    Vector mTrueNormalSlave = ZeroVector(3);

    double mMasterTrueDotSurrogateNormal = 0.0;
    double mSlaveTrueDotSurrogateNormal = 0.0;

    double mIntegrationWeightMaster = 0.0;
    double mIntegrationWeightSlave = 0.0;

    ///@}

    private:
        ///@name Serialization
        ///@{

        // void GetDBMatrix(IndexType index, GeometryType::Pointer r_geometry, 
        //                                         Vector& old_displacement, Matrix& DB, const Kratos::ProcessInfo& rCurrentProcessInfo);

        const Matrix GetConstitutiveMatrix(IndexType index, Matrix& r_B, GeometryType r_geometry,
                                            Vector& old_displacement, const Kratos::ProcessInfo& rCurrentProcessInfo,
                                            Vector& stress_vector); 
        
        PropertiesType& GetProperty(IndexType index) {
            if (index == 0) {
                return *mpPropMaster;
            }
            else if (index == 1) {
                return *mpPropSlave;
            }
        }

        ConstitutiveLaw::Pointer GetConstitutiveLaw(IndexType index) {
            if (index == 0) {
                return mpConstitutiveLawMaster;
            }
            else if (index == 1) {
                return mpConstitutiveLawSlave;
            }
        }

        void GetDeformed(const Matrix& N, Vector& reference_position, Vector& displacement, Vector& deformed_position);

        double ComputeTaylorTerm(const double derivative, const double dx, const int n_k, const double dy, const int k) const;
                           
        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
        }

        virtual void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
        }

        ///@}

    }; // Class SupportPenaltyLaplacianCondition

}  // namespace Kratos.
