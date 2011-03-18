/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

#ifndef KRATOS_BINGHAM_VMS_H
#define	KRATOS_BINGHAM_VMS_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "vms.h"

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
    ///@name Kratos Elements
    ///@{

    /// An element for the simulation of Bingham plastics.
    /**
     * This class implements a stabilized formulation based on the
     * Variational Multiscale framework to simulate a non-Newtonian fluid using the
     * Bingham plasic model. The the subscales can be modeled
     * using either Algebraic Subgird Scales (ASGS) or Orthogonal Subscales (OSS).
     * In the case of OSS, the projection terms are treated explicitly (computed
     * using the results of the previous iteration) and the subscales are not
     * tracked in time. The choice of subscale model is made based on the ProcessInfo
     * variable OSS_SWITCH (OSS if 1, ASGS otherwise).
     * This class implements both the 2D and 3D versions of the element.
     *
     * Note that the terms in the stabilized momentum equation that depend on gradients
     * of viscosity are not taken into account, as in this model viscosity is constant
     * inside each element.
     *
     * This class requires at least the following variables:\n
     * On each Node, as solution step variables VELOCITY, PRESSURE, ACCELERATION, MESH_VELOCITY.\n
     * On ProcessInfo OSS_SWITCH, DYNAMIC_TAU, DELTA_TIME, YIELD_STRESS.\n
     * If OSS is used, the nodes also require NODAL_AREA, ADVPROJ and DIVPROJ as solution step variables.\n
     * To compute the effective viscosity, MU and TAU have to be defined on the elements.\n
     * Some additional variables can be used to print results on the element: TAUONE, TAUTWO, TAU, MU, VORTICITY.
     *
     * @see ResidualBasedEliminationBuilderAndSolver compatible monolithic solution strategy.
     * @see PressureSplittingBuilderAndSolver compatible segregated solution strategy.
     * @see TrilinosPressureSplittingBuilderAndSolver compatible mpi strategy.
     * @see ResidualBasedPredictorCorrectorVelocityBossakScheme time scheme that can use
     * OSS stabilization.
     */
    template< unsigned int TDim,
              unsigned int TNumNodes = TDim + 1 >
    class BinghamVMS : public VMS<TDim,TNumNodes>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of VMS
        KRATOS_CLASS_POINTER_DEFINITION(BinghamVMS);

        ///base type: an IndexedObject that automatically has a unique number
        typedef IndexedObject BaseType;

        ///definition of node type (default is: Node<3>)
        typedef Node < 3 > NodeType;

        /**
         * Properties are used to store any parameters
         * related to the constitutive law
         */
        typedef Properties PropertiesType;

        ///definition of the geometry type with given NodeType
        typedef Geometry<NodeType> GeometryType;

        ///definition of nodes container type, redefined from GeometryType
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        typedef Vector VectorType;

        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

        typedef VMS<TDim,TNumNodes> BaseElementType;

        ///@}
        ///@name Life Cycle
        ///@{

        //Constructors.

        /// Default constuctor.
        /**
         * @param NewId Index number of the new element (optional)
         */
        BinghamVMS(IndexType NewId = 0) :
            BaseElementType(NewId)
        {}

        ///Constructor using an array of nodes.
        /**
         * @param NewId Index of the new element
         * @param ThisNodes An array containing the nodes of the new element
         */
        BinghamVMS(IndexType NewId, const NodesArrayType& ThisNodes) :
            BaseElementType(NewId, ThisNodes)
        {}

        /// Constructor using a geometry object.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         */
        BinghamVMS(IndexType NewId, GeometryType::Pointer pGeometry) :
            BaseElementType(NewId, pGeometry)
        {}

        /// Constuctor using geometry and properties.
        /**
         * @param NewId Index of the new element
         * @param pGeometry Pointer to a geometry object
         * @param pProperties Pointer to the element's properties
         */
        BinghamVMS(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
            BaseElementType(NewId, pGeometry, pProperties)
        {}

        /// Destructor.
        virtual ~BinghamVMS()
        {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /// Create a new element of this type
        /**
         * Returns a pointer to a new BinghamVMS element, created using given input
         * @param NewId: the ID of the new element
         * @param ThisNodes: the nodes of the new element
         * @param pProperties: the properties assigned to the new element
         * @return a Pointer to the new element
         */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
                                PropertiesType::Pointer pProperties) const
        {
            return Element::Pointer(new BinghamVMS(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Obtain an array_1d<double,3> elemental variable, evaluated on gauss points.
        /**
         * If the variable is VORTICITY, computes the vorticity (rotational of the velocity)
         * based on the current velocity values. Otherwise, it assumes that the input
         * variable is an elemental value and retrieves it. Implemented for a
         * single gauss point only.
         * @param rVariable Kratos vector variable to get
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                                 std::vector<array_1d<double, 3 > >& rOutput,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {
            BaseElementType::GetValueOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
        }

        /// Obtain a double elemental variable, evaluated on gauss points.
        /**
         * If the variable is TAUONE or TAUTWO, calculates the corresponding stabilization
         * parameter for the element, based on rCurrentProcessInfo's DELTA_TIME and
         * DYNAMIC_TAU. Otherwise, it assumes that the input variable is an
         * elemental value and retrieves it.
         * For MU and TAU, returns the elemental values of viscosity and stress, as
         * given by the Bingham model. Implemented for a single gauss point only.
         * @param rVariable Kratos vector variable to compute
         * @param Output Will be filled with the values of the variable on integrartion points
         * @param rCurrentProcessInfo Process info instance
         */
        virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                 std::vector<double>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {
            if (rVariable == TAUONE || rVariable == TAUTWO || rVariable == MU || rVariable==TAU || rVariable==EQ_STRAIN_RATE)
            {
                double TauOne, TauTwo;
                double Area;
                array_1d<double, TNumNodes> N;
                boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim> DN_DX;
                GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

                array_1d<double, 3 > AdvVel;
                this->GetAdvectiveVel(AdvVel, N);

                double Density, KinViscosity;
                this->EvaluateInPoint(Density, DENSITY, N);
                this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

                double Viscosity;
                this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

                this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Viscosity, rCurrentProcessInfo);

                rValues.resize(1, false);
                if (rVariable == TAUONE)
                {
                    rValues[0] = TauOne;
                }
                else if (rVariable == TAUTWO)
                {
                    rValues[0] = TauTwo;
                }
                else if (rVariable == MU)
                {
                    rValues[0] = Viscosity;
                }
                else if (rVariable == TAU)
                {
                    double NormS = this->SymmetricGradientNorm(DN_DX);
                    rValues[0] = Density*Viscosity*NormS*1.414213562;
                }
                else if (rVariable == EQ_STRAIN_RATE)
                {
                    double NormS = this->SymmetricGradientNorm(DN_DX);
                    rValues[0] = NormS*1.414213562;
                }

            }
            else // Default behaviour (returns elemental data)
            {
                rValues.resize(1, false);
                /*
                 The cast is done to avoid modification of the element's data. Data modification
                 would happen if rVariable is not stored now (would initialize a pointer to &rVariable
                 with associated value of 0.0). This is catastrophic if the variable referenced
                 goes out of scope.
                 */
                const VMS<TDim, TNumNodes>* const_this = static_cast<const VMS<TDim, TNumNodes>*> (this);
                rValues[0] = const_this->GetValue(rVariable);
            }
        }

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                                 std::vector<array_1d<double, 6 > >& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                 std::vector<Vector>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

        /// Empty implementation of unused CalculateOnIntegrationPoints overloads to avoid compilation warning
        virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                                 std::vector<Matrix>& rValues,
                                                 const ProcessInfo& rCurrentProcessInfo)
        {}

        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "BinghamVMS #" << this->Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "BinghamVMS" << TDim << "D";
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


        /// Calculate the viscosity, as given by the Bingham model.
        /**
         * @param Density Fluid density evaluated at the integration point
         * @param MolecularViscosity Viscosity of the fluid, in kinematic units (m2/s)
         * @param rShapeFunc Elemental shape functions, evaluated on the integration point
         * @param rShapeDeriv Shape function derivatives, evaluated on the integration point
         * @param TotalViscosity Effective viscosity (output)
         * @param rCurrentProcessInfo ProcessInfo instance (Checked for YIELD_STRESS)
         */
        virtual void GetEffectiveViscosity(const double Density,
                                           const double MolecularViscosity,
                                           array_1d<double, TNumNodes>& rShapeFunc,
                                           const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TDim >& rShapeDeriv,
                                           double& TotalViscosity,
                                           const ProcessInfo& rCurrentProcessInfo)
        {
            const double yield = rCurrentProcessInfo.GetValue(YIELD_STRESS);

            TotalViscosity = MolecularViscosity;
            if ( yield != 0)
            {
                const double NormS = this->SymmetricGradientNorm(rShapeDeriv);

                //non newtonian case
                const double mcoef = rCurrentProcessInfo.GetValue(M);
                double aux_1 = 1.0 - exp(-(mcoef * 1.414213562 * NormS));
                if(1.414213562 * NormS < 1.0/(1000.0*mcoef))
                    TotalViscosity += yield*mcoef;
                else
                    TotalViscosity += (yield / (1.414213562 * NormS)) * aux_1;

                TotalViscosity/=Density;
            }
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

        virtual void save(Serializer& rSerializer) const;

        virtual void load(Serializer& rSerializer)
        {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElementType);
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
        BinghamVMS & operator=(BinghamVMS const& rOther);

        /// Copy constructor.
        BinghamVMS(BinghamVMS const& rOther);

        ///@}

    }; // Class BinghamVMS

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    template< unsigned int TDim,
              unsigned int TNumNodes >
    inline std::istream& operator >>(std::istream& rIStream,
                                     BinghamVMS<TDim, TNumNodes>& rThis)
    {
        return rIStream;
    }

    /// output stream function
    template< unsigned int TDim,
              unsigned int TNumNodes >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const BinghamVMS<TDim, TNumNodes>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}

    ///@} // Fluid Dynamics Application group

} // namespace Kratos.

#endif	/* KRATOS_BINGHAM_VMS_H */

