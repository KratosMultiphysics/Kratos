// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//

#pragma once

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

    template< unsigned int TDim, unsigned int TNumNodes >
    class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement :
        public Element
    {
    public:

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);

        using BaseType  = Element ;

        using IndexType = std::size_t;
        using PropertiesType = Properties;
        using NodeType = Node;
        using GeometryType = Geometry<NodeType>;
        using NodesArrayType = Geometry<NodeType>::PointsArrayType;
        using VectorType = Vector;
        using MatrixType = Matrix;
        using DofsVectorType = Element::DofsVectorType;
        using EquationIdVectorType = Element::EquationIdVectorType;

        /// The definition of the sizetype
        using SizeType = std::size_t;

    	std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
        bool mIsInitialised = false;
        bool mIsPressureCoupled = false;

        struct ElementVariables
        {
            ///Properties variables
            double WaterDensity;
            double SolidDensity;
            double WaterHeatCapacity;
            double SolidHeatCapacity;
            double WaterThermalConductivity;
            double SolidThermalConductivityXX;
            double SolidThermalConductivityXY;
            double SolidThermalConductivityYX;
            double SolidThermalConductivityYY;
            double Porosity;
            double Saturation;
            double LongitudinalDispersivity;
            double TransverseDispersivity;
            double SolidCompressibility;
            double DynamicViscosityInverse;

            ///ProcessInfo variables
            double DtTemperatureCoefficient;

            ///Nodal variables
            array_1d<double, TNumNodes> TemperatureVector;
            array_1d<double, TNumNodes> DtTemperatureVector;
            array_1d<double, TDim> DischargeVector;

            ///Constitutive Law parameters
            Matrix ConstitutiveMatrix;
            Vector N;
            Matrix GradNT;
            Matrix GradNTInitialConfiguration;
            BoundedMatrix<double, TDim, TDim> PermiabilityMatrix;

            Vector detJContainer;
            Matrix NContainer;
            GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

            // needed for updated Lagrangian:
            double detJ;
            double IntegrationCoefficient;

            //Auxiliary Variables
            BoundedMatrix<double, TNumNodes, TNumNodes> TMatrix;
            array_1d<double, TNumNodes> TVector;
            BoundedMatrix<double, TNumNodes, TDim> TDimMatrix;
        };

    	/// Default Constructor
        explicit TransientThermalElement(IndexType NewId = 0);

        /// Constructor using an array of nodes
        TransientThermalElement(IndexType NewId,
            const NodesArrayType& ThisNodes);

        /// Constructor using Geometry
        TransientThermalElement(IndexType NewId,
            GeometryType::Pointer pGeometry);

        /// Constructor using Properties
        TransientThermalElement(IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties);

        // Assignment operator.
        TransientThermalElement& operator=(TransientThermalElement const& rOther) = delete;

        // Copy constructor.
        TransientThermalElement(TransientThermalElement const& rOther) = delete;

        /// Destructor
        ~TransientThermalElement() override;

        Element::Pointer Create(IndexType NewId,
            NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const override;

        Element::Pointer Create(IndexType NewId,
            GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const override;

        int Check(const ProcessInfo & rCurrentProcessInfo) const override;

        void Initialize(const ProcessInfo & rCurrentProcessInfo) override;

        void GetDofList(DofsVectorType & rElementalDofList,
            const ProcessInfo & rCurrentProcessInfo) const override;

        void EquationIdVector(EquationIdVectorType & rResult,
            const ProcessInfo & rCurrentProcessInfo) const override;

        void CalculateDischargeVector(ElementVariables& rVariables);


    protected:

        void CalculateAll(MatrixType & rLeftHandSideMatrix,
            VectorType & rRightHandSideVector,
            const ProcessInfo & CurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag);

        void InitializeElementVariables(ElementVariables & rVariables,
            const ProcessInfo & CurrentProcessInfo);

        void CalculateAndAddLHS(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);

        void CalculateAndAddRHS(VectorType & rRightHandSideVector, ElementVariables & rVariables, unsigned int GPoint);

        void CalculateKinematics(ElementVariables & rVariables, unsigned int PointNumber);

        void CalculateAndAddConductivityMatrix(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);
        void CalculateAndAddCapacityMatrix(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);
        void CalculateAndAddConvectionMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
        void CalculateAndAddConductivityVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);
        void CalculateAndAddCapacityVector(VectorType& rRightHandSideVector, ElementVariables & rVariables);

        void InitializeNodalTemperatureVariables(ElementVariables& rVariables);

        unsigned int GetNumberOfDOF() const;

        virtual double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
            unsigned int PointNumber, double detJ);

        void InitializeProperties(ElementVariables& rVariables);

        virtual void CalculateConductivityMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            ElementVariables& rVariables);

        virtual void CalculateCapacityMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            ElementVariables& rVariables) const;

        void CalculateConvectionMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            ElementVariables& rVariables);
        
        virtual void CalculateCapacityVector(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            array_1d<double, TNumNodes>& TVector, ElementVariables& rVariables) const;

        virtual void  CalculateConductivityVector(BoundedMatrix<double, TNumNodes, TDim>& TDimMatrix,
            BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix, array_1d<double, TNumNodes>& TVector,
            const ElementVariables& rVariables);

        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculatePermiabilityMatrix(BoundedMatrix<double, TDim, TDim>& C, ElementVariables& rVariables);

        GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    private:

        /// Serialization

        friend class Serializer;

        void save(Serializer & rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        }

        void load(Serializer & rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
        }

    }; // Class TransientThermalElement

} // namespace Kratos