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

        typedef Element BaseType;

        typedef std::size_t IndexType;
        typedef Properties PropertiesType;
        typedef Node <3> NodeType;
        typedef Geometry<NodeType> GeometryType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef Vector VectorType;
        typedef Matrix MatrixType;
        typedef Element::DofsVectorType DofsVectorType;
        typedef Element::EquationIdVectorType EquationIdVectorType;

        /// The definition of the sizetype
        typedef std::size_t SizeType;

    	std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
        bool mIsInitialised = false;

//        typedef typename BaseType::ElementVariables ElementVariables;

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
            double SolidThermalConductivityXX;
            double SolidThermalConductivityXY;
            double SolidThermalConductivityYX;
            double SolidThermalConductivityYY;
            double Porosity;
            double Saturation;

            ///ProcessInfo variables
            double DtTemperatureCoefficient;

            ///Nodal variables
            array_1d<double, TNumNodes> TemperatureVector;
            array_1d<double, TNumNodes> DtTemperatureVector;

            ///Variables computed at each GP
            Matrix B;

            ///Constitutive Law parameters
            BoundedMatrix<double, TDim, TDim> ConstitutiveMatrix;
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
            BoundedMatrix<double, TNumNodes, TDim> TDimMatrix;
        };

            	/// Default Constructor
        TransientThermalElement(IndexType NewId = 0) : BaseType(NewId) {}

        /// Constructor using an array of nodes
        TransientThermalElement(IndexType NewId,
            const NodesArrayType & ThisNodes) : BaseType(NewId, ThisNodes) {}

        /// Constructor using Geometry
        TransientThermalElement(IndexType NewId,
            GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry) {}

        /// Constructor using Properties
        TransientThermalElement(IndexType NewId,
            GeometryType::Pointer pGeometry,
            PropertiesType::Pointer pProperties) : BaseType(NewId, pGeometry, pProperties) {}

        /// Destructor
        ~TransientThermalElement() override {}

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

        void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;
        void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;
        void GetValuesVector(Vector& rValues, int Step) const override;
        void GetFirstDerivativesVector(Vector& rValues, int Step) const override;
        void GetSecondDerivativesVector(Vector& rValues, int Step) const override;


    protected:

        GeometryData::IntegrationMethod mThisIntegrationMethod;
        /// Member Variables

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
        
        virtual void CalculateCapacityVector(BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix,
            array_1d<double, TNumNodes>& TVector, ElementVariables& rVariables) const;

        virtual void  CalculateConductivityVector(BoundedMatrix<double, TNumNodes, TDim>& TDimMatrix,
            BoundedMatrix<double, TNumNodes, TNumNodes>& TMatrix, array_1d<double, TNumNodes>& TVector,
            const ElementVariables& rVariables);

        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo) override;

        void CalculateThermalDispersionMatrix(BoundedMatrix<double, TDim, TDim>& C, ElementVariables& rVariables);

        void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
        void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

        GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    private:

        /// Member Variables

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

        // Assignment operator.
        TransientThermalElement& operator=(TransientThermalElement const& rOther);

        // Copy constructor.
        TransientThermalElement(TransientThermalElement const& rOther);

        template <class TValueType>
        inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
        {
            rNode.SetLock();
            rNode.FastGetSolutionStepValue(Var) = Value;
            rNode.UnSetLock();
        }

    }; // Class TransientThermalElement

} // namespace Kratos

#endif // KRATOS_GEO_THERMAL_ELEMENT_H_INCLUDED  defined