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
//                   John van Esch
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

        bool mIsInitialised = false;

        struct ElementVariables
        {
            double waterDensity;
            double solidDensity;
            double waterHeatCapacity;
            double solidHeatCapacity;
            double waterThermalConductivity;
            double solidThermalConductivityXX;
            double solidThermalConductivityXY;
            double solidThermalConductivityYX;
            double solidThermalConductivityYY;
            double porosity;
            double saturation;

            double dtTemperatureCoefficient;
            array_1d<double, TNumNodes> temperatureVector;
            array_1d<double, TNumNodes> dtTemperatureVector;
            Matrix constitutiveMatrix;
            Vector N;
            Matrix GradNT;
            Matrix GradNTInitialConfiguration;
            Vector detJContainer;
            Matrix NContainer;
            GeometryType::ShapeFunctionsGradientsType DN_DXContainer;
            double detJ;
            double IntegrationCoefficient;
            BoundedMatrix<double, TNumNodes, TNumNodes> conductivityMatrix;
            BoundedMatrix<double, TNumNodes, TNumNodes> capacityMatrix;
            array_1d<double, TNumNodes> conductivityVector;
            array_1d<double, TNumNodes> capacityVector;
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


    protected:

        void CalculateAll(MatrixType & rLeftHandSideMatrix,
            VectorType & rRightHandSideVector,
            const ProcessInfo & CurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag);

        void InitializeElementVariables(ElementVariables & rVariables,
            const ProcessInfo & CurrentProcessInfo);

        void CalculateAndAddLHS(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);

        void CalculateAndAddRHS(VectorType & rRightHandSideVector, ElementVariables & rVariables);

        void CalculateKinematics(ElementVariables & rVariables, unsigned int PointNumber);

        void CalculateAndAddConductivityMatrix(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);
        void CalculateAndAddCapacityMatrix(MatrixType & rLeftHandSideMatrix, ElementVariables & rVariables);
        void CalculateAndAddConductivityVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);
        void CalculateAndAddCapacityVector(VectorType& rRightHandSideVector, ElementVariables & rVariables);

        void InitializeNodalTemperatureVariables(ElementVariables& rVariables);

        unsigned int GetNumberOfDOF() const;

        virtual double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
            unsigned int PointNumber, double detJ);

        void InitializeProperties(ElementVariables& rVariables);

        virtual void CalculateConductivityMatrix(ElementVariables& rVariables);

        virtual void CalculateCapacityMatrix(ElementVariables& rVariables) const;
        
        virtual void CalculateCapacityVector(ElementVariables& rVariables) const;

        virtual void  CalculateConductivityVector(ElementVariables& rVariables);

        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

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