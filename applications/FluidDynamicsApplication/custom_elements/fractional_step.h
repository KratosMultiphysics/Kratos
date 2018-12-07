//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FRACTIONAL_STEP_H_INCLUDED )
#define  KRATOS_FRACTIONAL_STEP_H_INCLUDED

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

    ///@addtogroup FluidDynamicsApplication
    ///@{

    ///@name Kratos Globals
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name  Functions
    ///@{

    ///@}
    ///@name Kratos Classes
    ///@{

    /// A stabilized element for the incompressible Navier-Stokes equations.
    /**
     */
    template< unsigned int TDim >
    class FractionalStep : public Element
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of FractionalStep
        KRATOS_CLASS_POINTER_DEFINITION(FractionalStep);

        /// Node type (default is: Node<3>)
        typedef Node <3> NodeType;

        /// Geometry type (using with given NodeType)
        typedef Geometry<NodeType> GeometryType;

        /// Definition of nodes container type, redefined from GeometryType
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        /// Vector type for local contributions to the linear system
        typedef Vector VectorType;

        /// Matrix type for local contributions to the linear system
        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        /// Type for shape function values container
        typedef Kratos::Vector ShapeFunctionsType;

        /// Type for a matrix containing the shape function gradients
        typedef Kratos::Matrix ShapeFunctionDerivativesType;

        /// Type for an array of shape function gradient matrices
        typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

        ///@}
        ///@name Life Cycle
        ///@{

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        FractionalStep(IndexType NewId = 0) :
            Element(NewId)
        {}

        /// Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        FractionalStep(IndexType NewId, const NodesArrayType& ThisNodes) :
            Element(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        FractionalStep(IndexType NewId, GeometryType::Pointer pGeometry) :
            Element(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        FractionalStep(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            Element(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        ~FractionalStep() override
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new FractionalStep element, created using given input
         * @param NewId the ID of the new element
         * @param ThisNodes the nodes of the new element
         * @param pProperties the properties assigned to the new element
         * @return a Pointer to the new element
         */
	    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                Element::PropertiesType::Pointer pProperties) const override
	    {
            return Kratos::make_shared< FractionalStep<TDim> >(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
	    }
	
        /**
         * Returns a pointer to a new FractionalStep element, created using given input
         * @param NewId the ID of the new element
         * @param pGeom a pointer to the geometry
         * @param pProperties the properties assigned to the new element
         * @return a Pointer to the new element
         */
		
        Element::Pointer Create(IndexType NewId, Element::GeometryType::Pointer pGeom, Element::PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_shared< FractionalStep<TDim> >(NewId, pGeom, pProperties);
        }
        
        /**
         * Clones the selected element variables, creating a new one
         * @param NewId the ID of the new element
         * @param ThisNodes the nodes of the new element
         * @param pProperties the properties assigned to the new element
         * @return a Pointer to the new element
         */
        
        Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
        {
            Element::Pointer pNewElement = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );
            
            pNewElement->SetData(this->GetData());
            pNewElement->SetFlags(this->GetFlags());
            
            return pNewElement;
        }
        
        void Initialize() override;

        /// Initializes the element and all geometric information required for the problem.
        void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override;


        void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo) override;

        /// Calculate the element's local contribution to the system for the current step.
        void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          ProcessInfo& rCurrentProcessInfo) override;


        void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                           ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;
            KRATOS_THROW_ERROR(std::logic_error,"FractionalStep::CalculateLeftHandSide not implemented","");
            KRATOS_CATCH("");
        }

        void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override
        {
            KRATOS_TRY;
            KRATOS_THROW_ERROR(std::logic_error,"FractionalStep::CalculateRightHandSide not implemented","");
            KRATOS_CATCH("");
        }

        /**
         * @param rVariable Use DIVPROJ
         * @param rOutput (unused)
         * @param rCurrentProcessInfo Process info instance (unused)
         */
        void Calculate(const Variable<double>& rVariable,
                               double& rOutput,
                               const ProcessInfo& rCurrentProcessInfo) override;


        /**
         * @param rVariable Use ADVPROJ or VELOCITY
         * @param Output (unused)
         * @param rCurrentProcessInfo Process info instance (unused)
         */
        void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                               array_1d<double, 3 > & rOutput,
                               const ProcessInfo& rCurrentProcessInfo) override;


        // The following methods have different implementations depending on TDim
        /// Provides the global indices for each one of this element's local rows
        /**
         * this determines the elemental equation ID vector for all elemental
         * DOFs
         * @param rResult A vector containing the global Id of each row
         * @param rCurrentProcessInfo the current process info object (unused)
         */
        void EquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo) override;

        /// Returns a list of the element's Dofs
        /**
         * @param ElementalDofList the list of DOFs
         * @param rCurrentProcessInfo the current process info instance
         */
        void GetDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo) override;

    
        GeometryData::IntegrationMethod GetIntegrationMethod() const override;

        /// Obtain an array_1d<double,3> elemental variable, evaluated on gauss points.
        /**
         * @param rVariable Kratos vector variable to get
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                std::vector<array_1d<double, 3 > >& rValues,
                const ProcessInfo& rCurrentProcessInfo) override;

        /// Obtain a double elemental variable, evaluated on gauss points.
        /**
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                std::vector<double>& rValues,
                const ProcessInfo& rCurrentProcessInfo) override;

        /// Obtain an array_1d<double,6> elemental variable, evaluated on gauss points.
        /**
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                std::vector<array_1d<double, 6 > >& rValues,
                const ProcessInfo& rCurrentProcessInfo) override
        {
            this->GetElementalValueForOutput< array_1d<double,6> >(rVariable,rValues);
        }

        /// Obtain a Vector elemental variable, evaluated on gauss points.
        /**
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                std::vector<Vector>& rValues,
                const ProcessInfo& rCurrentProcessInfo) override
        {
            this->GetElementalValueForOutput< Vector >(rVariable,rValues);
        }

        /// Obtain a Matrix elemental variable, evaluated on gauss points.
        /**
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                std::vector<Matrix>& rValues,
                const ProcessInfo& rCurrentProcessInfo) override
        {
            this->GetElementalValueForOutput< Matrix >(rVariable,rValues);
        }

        ///@}
        ///@name Access
        ///@{

        ///@}
        ///@name Elemental Data
        ///@{

        /// Checks the input and that all required Kratos variables have been registered.
        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
         * @return 0 if no errors were found.
         */
        int Check(const ProcessInfo& rCurrentProcessInfo) override;

        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "FractionalStep #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override
        {
            rOStream << "FractionalStep" << TDim << "D";
        }

//        /// Print object's data.
//        virtual void PrintData(std::ostream& rOStream) const;

        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{

        void CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix,
                                                    VectorType& rRightHandSideVector,
                                                    const ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix,
                                          VectorType& rRightHandSideVector,
                                          const ProcessInfo& rCurrentProcessInfo);

        void CalculateEndOfStepSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo);

        void VelocityEquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        void PressureEquationIdVector(EquationIdVectorType& rResult,
                                      ProcessInfo& rCurrentProcessInfo);

        void GetVelocityDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        void GetPressureDofList(DofsVectorType& rElementalDofList,
                                ProcessInfo& rCurrentProcessInfo);

        void GetPressureValues(Vector& rValues,
                               const int Step = 0);

        void GetVelocityValues(Vector& rValues,
                               const int Step = 0);

        /// Determine integration point weights and shape funcition derivatives from the element's geometry.
        virtual void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
                                           Matrix& rNContainer,
                                           Vector& rGaussWeights);

        double ElementSize(/*ShapeFunctionDerivativesType& rDN_DX*/);


        /**
         * @brief EffectiveViscosity Calculate the viscosity at given integration point, using Smagorinsky if enabled.
         *
         * The Smagorinsky model is used only if the C_SMAGORINSKY is defined on the elemental data container.
         *
         * @note: This function is redefined when using Non-Newtonian constitutive models. It is important to keep its
         * signature, otherwise non-Newtonian models will stop working.
         *
         * @param Density The fluid's density at the integration point.
         * @param rN Nodal shape functions evaluated at the integration points (area coordinates for the point).
         * @param rDN_DX Shape function derivatives at the integration point.
         * @param ElemSize Representative length of the element (used only for Smagorinsky).
         * @param rProcessInfo ProcessInfo instance passed from the ModelPart.
         * @return Effective viscosity, in dynamic units (Pa*s or equivalent).
         */
        virtual double EffectiveViscosity(double Density,
                                          const ShapeFunctionsType &rN,
                                          const ShapeFunctionDerivativesType &rDN_DX,
                                          double ElemSize,
                                          const ProcessInfo &rProcessInfo);

        /**
         * @brief EquivalentStrainRate Calculate the second invariant of the strain rate tensor GammaDot = (2SijSij)^0.5.
         *
         * @note Our implementation of non-Newtonian consitutive models such as Bingham relies on this funcition being
         * defined on all fluid elements.
         *
         * @param rDN_DX Shape function derivatives at the integration point.
         * @return GammaDot = (2SijSij)^0.5.
         */
        double EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const;


        /// Add integration point contribution to the mass matrix.
        /**
         * A constistent mass matrix is used.
         * @param rMassMatrix The local matrix where the result will be added.
         * @param rN Elemental shape functions.
         * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
         */
        void AddMomentumMassTerm(Matrix& rMassMatrix,
                                 const ShapeFunctionsType& rN,
                                 const double Weight);

        virtual void AddMomentumSystemTerms(Matrix& rLHSMatrix,
                                    Vector& rRHSVector,
                                    const double Density,
                                    const Vector& rConvOperator,
                                    const array_1d<double,3>& rBodyForce,
                                    const double OldPressure,
                                    const double TauOne,
                                    const double TauTwo,
                                    const array_1d<double,3>& rMomentumProjection,
                                    const double MassProjection,
                                    const ShapeFunctionsType& rN,
                                    const ShapeFunctionDerivativesType& rDN_DX,
                                    const double Weight);

        void AddViscousTerm(MatrixType& rDampingMatrix,
                            const ShapeFunctionDerivativesType& rShapeDeriv,
                            const double Weight);


//        virtual void CalculateTau(double& TauOne,
//                                  double& TauTwo,
//                                  double ElemSize,
//                                  const ProcessInfo& rCurrentProcessInfo);

        /// Calculate Stabilization parameters.
        /**
         * Calculates both tau parameters based on a given advective velocity.
         * Takes time step and dynamic coefficient from given ProcessInfo instance.
         * ProcessInfo variables DELTA_TIME and DYNAMIC_TAU will be used.
         * @param TauOne First stabilization parameter (momentum equation)
         * @param TauTwo Second stabilization parameter (mass equation)
         * @param ElemSize Characteristic element size (h)
         * @param rAdvVel advection velocity
         * @param Area Elemental area
         * @param Density Density on integrartion point
         * @param Viscosity Dynamic viscosity (mu) on integrartion point
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void CalculateTau(double& TauOne,
                                  double& TauTwo,
                                  double ElemSize,
                                  const array_1d< double, 3 > & rAdvVel,
                                  const double Density,
                                  const double Viscosity,
                                  const ProcessInfo& rCurrentProcessInfo);

        virtual void CalculateProjectionRHS(VectorType& rMomentumRHS,
                                            VectorType& rMassRHS,
                                            const ShapeFunctionsType& rN,
                                            const ShapeFunctionDerivativesType& rDN_DX,
                                            const double Weight);

        virtual void CalculateProjectionRHS(VectorType& rConvTerm,
                                            VectorType& rPresTerm,
                                            VectorType& rDivTerm,
                                            const ShapeFunctionsType& rN,
                                            const ShapeFunctionDerivativesType& rDN_DX,
                                            const double Weight);


        void ModulatedGradientDiffusion(MatrixType& rDampingMatrix,
                                        const ShapeFunctionDerivativesType& rDN_DX,
                                        const double Weight);


        /// Write the convective operator evaluated at this point (for each nodal funciton) to an array
        /**
         * Evaluate the convective operator for each node's shape function at an arbitrary point
         * @param rResult Output vector
         * @param rVelocity Velocity evaluated at the integration point
         * @param rShapeDeriv Derivatives of shape functions evaluated at the integration point
         * @see GetAdvectiveVel provides rVelocity
         */
        void ConvectionOperator(Vector& rResult,
                                const array_1d<double,3>& rConvVel,
                                const ShapeFunctionDerivativesType& DN_DX);

        /// Evaluate ALE convective velocity (velocity-mesh velocity) at a given point.
        /**
         * @param rConvVel container for result.
         * @param N Shape functions at point of evaluation.
         */
        virtual void EvaluateConvVelocity(array_1d<double,3>& rConvVel,
                                          const ShapeFunctionsType& N);

        /// Write the value of a variable at a point inside the element to a double
        /**
         * Evaluate a nodal variable in the point where the form functions take the
         * values given by rShapeFunc and write the result to rResult.
         * This is an auxiliary function used to compute values in integration points.
         * @param rResult The variable where the value will be added to
         * @param rVariable The nodal variable to be read
         * @param rShapeFunc The values of the form functions in the point
         */
        template< class TVariableType >
        void EvaluateInPoint(TVariableType& rResult,
                             const Kratos::Variable<TVariableType>& Var,
                             const ShapeFunctionsType& rShapeFunc)
        {
            GeometryType& rGeom = this->GetGeometry();
            const SizeType NumNodes = rGeom.PointsNumber();

            rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

            for(SizeType i = 1; i < NumNodes; i++)
            {
                rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
            }
        }

        /// Write the value of a variable at a point inside the element to a double
        /**
         * Evaluate a nodal variable in the point where the form functions take the
         * values given by rShapeFunc and write the result to rResult.
         * This is an auxiliary function used to compute values in integration points.
         * @param rResult The variable where the value will be added to
         * @param rVariable The nodal variable to be read
         * @param rShapeFunc The values of the form functions in the point
         * @param Step Number of time steps back
         */
        template< class TVariableType >
        void EvaluateInPoint(TVariableType& rResult,
                             const Kratos::Variable<TVariableType>& Var,
                             const ShapeFunctionsType& rShapeFunc,
                             const IndexType Step)
        {
            GeometryType& rGeom = this->GetGeometry();
            const SizeType NumNodes = rGeom.PointsNumber();

            rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var,Step);

            for(SizeType i = 1; i < NumNodes; i++)
            {
                rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var,Step);
            }
        }

        void EvaluateGradientInPoint(array_1d<double,TDim>& rResult,
                                     const Kratos::Variable<double>& Var,
                                     const ShapeFunctionDerivativesType& rDN_DX)
        {
            GeometryType& rGeom = this->GetGeometry();
            const SizeType NumNodes = rGeom.PointsNumber();
	    
			const double& var = rGeom[0].FastGetSolutionStepValue(Var);
			for (SizeType d = 0; d < TDim; ++d)
				rResult[d] = rDN_DX(0,d) * var;
			
			for(SizeType i = 1; i < NumNodes; i++)
			{
			  const double& var = rGeom[i].FastGetSolutionStepValue(Var);
			  for (SizeType d = 0; d < TDim; ++d)
						rResult[d] += rDN_DX(i,d) * var;
				
			}
		}

        void EvaluateDivergenceInPoint(double& rResult,
                                       const Kratos::Variable< array_1d<double,3> >& Var,
                                       const ShapeFunctionDerivativesType& rDN_DX)
        {
            GeometryType& rGeom = this->GetGeometry();
            const SizeType NumNodes = rGeom.PointsNumber();

            rResult = 0.0;
            for(SizeType i = 0; i < NumNodes; i++)
            {
                const array_1d<double,3>& var = rGeom[i].FastGetSolutionStepValue(Var);
                for (SizeType d = 0; d < TDim; ++d)
                {
                    rResult += rDN_DX(i,d) * var[d];
                }
            }
        }

        /// Helper function to print results on gauss points
        /** Reads a variable from the element's database and returns it in a format
          * that can be used by GetValueOnIntegrationPoints functions.
          * @see GetValueOnIntegrationPoints
          */
        template<class TValueType>
        void GetElementalValueForOutput(const Kratos::Variable<TValueType>& rVariable,
                                        std::vector<TValueType>& rOutput)
        {
            unsigned int NumValues = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_2);
            rOutput.resize(NumValues);
            /*
             The cast is done to avoid modification of the element's data. Data modification
             would happen if rVariable is not stored now (would initialize a pointer to &rVariable
             with associated value of 0.0). This is catastrophic if the variable referenced
             goes out of scope.
             */
            const FractionalStep<TDim>* const_this = static_cast<const FractionalStep<TDim>*> (this);
            const TValueType& Val = const_this->GetValue(rVariable);

            for (unsigned int i = 0; i < NumValues; i++)
                rOutput[i] = Val;
        }

        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{

        ///@}
        ///@name Member Variables
        ///@{

        ///@}
        ///@name Serialization
        ///@{

        friend class Serializer;

        void save(Serializer& rSerializer) const override
        {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
        }

        void load(Serializer& rSerializer) override
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        }

        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        FractionalStep & operator=(FractionalStep const& rOther);

        /// Copy constructor.
        FractionalStep(FractionalStep const& rOther);

        ///@}

    }; // Class FractionalStep

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     FractionalStep<TDim>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const FractionalStep<TDim>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif // KRATOS_FRACTIONAL_STEP_H_INCLUDED  defined
