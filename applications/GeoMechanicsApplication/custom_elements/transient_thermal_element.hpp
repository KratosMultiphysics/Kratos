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


    protected:



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