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

#if !defined(KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED )
#define  KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_utilities/stress_strain_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"
#include "custom_constitutive/thermal_dispersion_2D_law.h"

namespace Kratos
{

    template< unsigned int TDim, unsigned int TNumNodes >
    class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientThermalElement :
        public Element
    {

    public:

        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientThermalElement);
        //
        typedef Element BaseType;


        struct ElementVariables
        {
            ///Properties variables
            double WaterDensity;
            double SolidDensity;
            double Density;
            double WaterHeatCapacity;
            double SolidHeatCapacity;
            double HeatCapacity;
            double WaterThermalConductivity;
            double SolidThermalConductivity;
            double ThermalConductivity;
            double Porosity;

            ///Nodal variables
            array_1d<double, TNumNodes> TemperatureVector;

            ///General elemental variables
            //Vector VoigtVector;

            ///Variables computed at each GP
            Matrix B;
            BoundedMatrix<double, TDim, TNumNodes* TDim> Nu;

            ///Constitutive Law parameters
            Matrix ConstitutiveMatrix;
            Vector Np;
            Vector NpT;
            Matrix GradNp;
            Matrix GradNpT;
            Matrix GradNpTInitialConfiguration;

            Matrix F;
            double detF;
            Vector detJContainer;
            Matrix NContainer;
            GeometryType::ShapeFunctionsGradientsType DN_DXContainer;

            // needed for updated Lagrangian:
            double detJ;
            double detJInitialConfiguration;
            double IntegrationCoefficient;
            double IntegrationCoefficientInitialConfiguration;

            //Auxiliary Variables
            BoundedMatrix<double, TNumNodes, TNumNodes> TMatrix;
            array_1d<double, TNumNodes> TVector;
        };


        /// Default Constructor
        TransientThermalElement(IndexType NewId = 0) : BaseType(NewId) {}

        /// Constructor using an array of nodes
        TransientThermalElement(IndexType NewId, const NodesArrayType& ThisNodes) : BaseType(NewId, ThisNodes) {}

        /// Constructor using Geometry
        TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry) {}

        /// Constructor using Properties
        TransientThermalElement(IndexType NewId, GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties) {}

        /// Destructor
        ~TransientThermalElement() override {}

        ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
            PropertiesType::Pointer pProperties) const override;

        Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom,
            PropertiesType::Pointer pProperties) const override;

        ///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        int Check(const ProcessInfo& rCurrentProcessInfo) const override;


    protected:

        void CalculateAll(MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& CurrentProcessInfo,
            const bool CalculateStiffnessMatrixFlag,
            const bool CalculateResidualVectorFlag);

        void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& CurrentProcessInfo);

        void CalculateKinematics(ElementVariables& rVariables, unsigned int PointNumber);

        void InitializeProperties(ElementVariables& rVariables);

        void InitializeNodalTemperatureVariables(ElementVariables& rVariables);

        virtual double CalculateIntegrationCoefficient(const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
            unsigned int PointNumber, double detJ);

        void CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);

        void CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, unsigned int GPoint);

        void CalculateAndAddConductivityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
        void CalculateAndAddCapacityMatrix(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables);
        void CalculateAndAddCapacityVector(VectorType& rRightHandSideVector, ElementVariables& rVariables);

        virtual void CalculateConductivityMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            const ElementVariables& rVariables) const;

        virtual void CalculateCapacityMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            const ElementVariables& rVariables) const;

        virtual void CalculateCapacityVector(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            array_1d<double, TNumNodes>& TVector, const ElementVariables& rVariables) const;

    }; // Class TransientThermalElement

} // namespace Kratos

#endif // KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED  defined