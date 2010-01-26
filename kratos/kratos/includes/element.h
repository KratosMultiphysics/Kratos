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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2007-10-12 10:28:14 $
//   Revision:            $Revision: 1.9 $
//
//


#if !defined(KRATOS_ELEMENT_H_INCLUDED )
#define  KRATOS_ELEMENT_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/process_info.h"
#include "utilities/indexed_object.h"

#include "containers/weak_pointer_vector.h"

namespace Kratos
{
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
    
    /// Base class for all Elements.
    /** 
     * This is the base class for all elements used in KRATOS 
     * Elements inherited from this class have to reimplement 
     * all public functions that are needed to perform their designated
     * tasks. Due to a dummy implementation of every function though,
     * not all of them have to be implemented if they are not needed for 
     * the actual problem
     */
    class Element : public IndexedObject
    {
        public:
            ///@name Type Definitions
            ///@{
            /// Pointer definition of Element
            KRATOS_CLASS_POINTER_DEFINITION(Element);
            
            ///base type: an IndexedObject that automatically has a unique number
            typedef IndexedObject BaseType;
            
            ///definition of node type (default is: Node<3>)
            typedef Node<3> NodeType;
            
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
            
            typedef PointerVectorSet<Dof<double> , IndexedObject> DofsArrayType;
            
            typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;
            
            ///Type definition for integration methods
            typedef GeometryData::IntegrationMethod IntegrationMethod;
            ///@}
            
            ///@name Life Cycle 
            ///@{ 
            /** 
             * Constructor.
             */
            Element(IndexType NewId = 0) : BaseType(NewId), 
                    mpGeometry(new GeometryType(NodesArrayType())), 
                               mpProperties(new PropertiesType) 
            {
            }
            
            /** 
             * Constructor using an array of nodes 
             */
            Element(IndexType NewId, const NodesArrayType& ThisNodes) : 
                    BaseType(NewId), 
                             mpGeometry(new GeometryType(ThisNodes)), 
                                        mpProperties(new PropertiesType) 
            {
            }
            
            /** 
             * Constructor using Geometry 
             */
            Element(IndexType NewId, GeometryType::Pointer pGeometry) : 
                    BaseType(NewId), mpGeometry(pGeometry), mpProperties(new PropertiesType) 
            {
            }
            
            /** 
             * Constructor using Properties
             */
            Element(IndexType NewId, GeometryType::Pointer pGeometry,  
                    PropertiesType::Pointer pProperties) : BaseType(NewId), 
                    mpGeometry(pGeometry), 
                               mpProperties(pProperties) 
            {
            }
            
            /// Copy constructor.
            Element(Element const& rOther) : BaseType(rOther), 
                    mpGeometry(rOther.mpGeometry), 
                               mpProperties(rOther.mpProperties) 
            {
            }
            
            /// Destructor.
            virtual ~Element(){}
            
            ///@}
            ///@name Operators 
            ///@{
            /// Assignment operator.
            Element& operator=(Element const& rOther)
            {
                BaseType::operator=(rOther);
                mpGeometry = rOther.mpGeometry;
                mpProperties = rOther.mpProperties;
                return *this;
            }
            
            ///@}
            ///@name Operations
            ///@{
            /**
             * creates a new element pointer
             * @param NewId: the ID of the new element
             * @param ThisNodes: the nodes of the new element
             * @param pProperties: the properties assigned to the new element 
             * @return a Pointer to the new element
             */
            virtual Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, 
                                   PropertiesType::Pointer pProperties) const
            {
                KRATOS_TRY
                return Element::Pointer(new Element(NewId, GetGeometry().Create(ThisNodes), pProperties));
                KRATOS_CATCH("");
            }
            
            /**
             * is called to initialize the element.
             * Must be called before any calculation is done!
             */
            virtual void Initialize(){}
            
            /**
             * this is called during the assembling process in order
             * to calculate all elemental contributions to the global system
             * matrix and the right hand side
             * @param rLeftHandSideMatrix: the elemental left hand side matrix
             * @param rRightHandSideVector: the elemental right hand side
             * @param rCurrentProcessInfo: the current process info instance
             */
            
            virtual void ResetConstitutiveLaw(){}
            
            virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                               VectorType& rRightHandSideVector, 
                                               ProcessInfo& rCurrentProcessInfo){}
            
            /**
             * this is called during the assembling process in order
             * to calculate all elemental contributions to the global system
             * matrix and the right hand sides for a matrix-formed RHS
             * @param rLeftHandSideMatrix: the elemental left hand side matrix
             * @param rDampMatrix: the elemental damping matrix
             * @param rRightHandSideVector1: the elemental right hand side (1)
             * @param rRightHandSideVector2: the elemental right hand side (1)
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                                               MatrixType& rDampMatrix, 
                                               VectorType& rRightHandSideVector1, 
                                               VectorType& rRightHandSideVector2, 
                                               ProcessInfo& rCurrentProcessInfo){}
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental left hand side matrix only
             * @param rLeftHandSideMatrix: the elemental left hand side matrix 
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                                                       ProcessInfo& rCurrentProcessInfo)
            {
                if(rLeftHandSideMatrix.size1()!=0)
                    rLeftHandSideMatrix.resize(0,0);
            }
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental right hand side vector only
             * @param rRightHandSideVector: the elemental right hand side vector
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void CalculateRightHandSide( VectorType& rRightHandSideVector, 
                    ProcessInfo& rCurrentProcessInfo)
            {
                if(rRightHandSideVector.size()!=0)
                    rRightHandSideVector.resize(0);
            }
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental left hand side vector only
             * @param rLeftHandSideVector: the elemental left hand side vector
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void CalculateLeftHandSide( VectorType& rLeftHandSideVector, 
                    ProcessInfo& rCurrentProcessInfo)
            {
                if(rLeftHandSideVector.size()!=0)
                    rLeftHandSideVector.resize(0);
            }
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental right hand side matrix only
             * @param rRightHandSideMatrix: the elemental right hand side matrix
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void CalculateRightHandSide( MatrixType& rRightHandSideMatrix, 
                    ProcessInfo& rCurrentProcessInfo)
            {
                if(rRightHandSideMatrix.size1()!=0)
                    rRightHandSideMatrix.resize(0,0);
            }
            
            /**
             * this determines the elemental equation ID vector for all elemental
             * DOFs
             * @param rResult: the elemental equation ID vector
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void EquationIdVector( EquationIdVectorType& rResult, 
                                           ProcessInfo& rCurrentProcessInfo)
            {
                if(rResult.size()!=0)
                    rResult.resize(0);
            }
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental mass matrix
             * @param rMassMatrix: the elemental mass matrix
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
            {
                if(rMassMatrix.size1()!=0)
                    rMassMatrix.resize(0,0);
            }
            
            /**
             * adds the mass matrix scaled by a given factor to the LHS
             * @param rLeftHandSideMatrix: the elemental LHS matrix
             * @param coeff: the given factor
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void AddMassMatrix( MatrixType& rLeftHandSideMatrix, 
                                        double coeff, ProcessInfo& rCurrentProcessInfo){}
            
            /**
             * this is called during the assembling process in order
             * to calculate the elemental damping matrix
             * @param rDampMatrix: the elemental damping matrix
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
            {
                if(rDampMatrix.size1()!=0)
                    rDampMatrix.resize(0,0);
            }
            
            /**
             * adds the inertia forces to the RHS --> performs residua = static_residua - coeff*M*acc
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void AddInertiaForces( VectorType& rRightHandSideVector, double coeff, 
                                           ProcessInfo& rCurrentProcessInfo){}
            /**
             * Calculate Damp matrix and add velocity contribution to RHS
            * @param rDampMatrix: the velocity-proportional "damping" matrix
            * @param rRightHandSideVector: the elemental right hand side matrix
             * @param rCurrentProcessInfo: the current process info instance
             */            
	    virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix,
							    VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
            {
                            if(rDampMatrix.size1()!=0)
                                 rDampMatrix.resize(0,0);

            }

            /**
             * determines the elemental list of DOFs
             * @param ElementalDofList: the list of DOFs
             * @param rCurrentProcessInfo: the current process info instance
             */
            virtual void GetDofList( DofsVectorType& ElementalDofList,
                                     ProcessInfo& CurrentProcessInfo)
            {
                if(ElementalDofList.size()!=0)
                    ElementalDofList.resize(0);
            }
            
            /**
             * this is called for non-linear analysis at the beginning of the iteration process
             */
            virtual void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo){}
            
            /**
             * this is called for non-linear analysis at the end of the iteration process
             */
            virtual void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo){}
            
            /**
             * this is called in the beginning of each solution step
             */
            virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo){}
            
            /**
             * this is called at the end of each solution step
             */
            virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo){}
            
            /**
             * deletes all obsolete data from memory
             */
            virtual void CleanMemory(){}
            
            virtual void GetValuesVector(Vector& values, int Step = 0){}
            virtual void GetFirstDerivativesVector(Vector& values, int Step = 0){}
            virtual void GetSecondDerivativesVector(Vector& values, int Step = 0){}
            
            //output On integration points
            //calculate on element
            virtual void Calculate( const Variable<double>& rVariable, 
                                    double& Output, 
                                    const ProcessInfo& rCurrentProcessInfo){}
            virtual void Calculate( const Variable<array_1d<double,3> >& rVariable, 
                                    array_1d<double,3>& Output, 
                                    const ProcessInfo& rCurrentProcessInfo){}
            virtual void Calculate( const Variable<Vector >& rVariable, 
                                    Vector& Output, 
                                    const ProcessInfo& rCurrentProcessInfo){}
            virtual void Calculate( const Variable<Matrix >& rVariable, 
                                    Matrix& Output, 
                                    const ProcessInfo& rCurrentProcessInfo){}
            
            //calculate on gauss points
            virtual void CalculateOnIntegrationPoints( const Variable<double>& rVariable, 
                                            Vector& Output, 
                                            const ProcessInfo& rCurrentProcessInfo){}
            virtual void CalculateOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                            std::vector< array_1d<double,3> >& Output, 
                                            const ProcessInfo& rCurrentProcessInfo){}
            virtual void CalculateOnIntegrationPoints( const Variable<Vector >& rVariable, 
                                            std::vector< Vector >& Output, 
                                            const ProcessInfo& rCurrentProcessInfo){}
            virtual void CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, 
                                            std::vector< Matrix >& Output, 
                                            const ProcessInfo& rCurrentProcessInfo){}
            
            /**
             * Access for variables on Integration points.
             * This gives access to variables stored in the constitutive law on each integration point.
             * Specialisations of element.h (e.g. the TotalLagrangian) must specify the actual 
             * interface to the constitutive law!
             * Note, that these functions expect a std::vector of values for the 
             * specified variable type that contains a value for each integration point!
             * SetValueOnIntegrationPoints: set the values for given Variable.
             * GetValueOnIntegrationPoints: get the values for given Variable.
             */
            virtual void SetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                           std::vector<double>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void SetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                           std::vector<array_1d<double,3> > rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
                                           std::vector<Vector>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                           std::vector<Matrix>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void GetValueOnIntegrationPoints( const Variable<double>& rVariable, 
                                           std::vector<double>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void GetValueOnIntegrationPoints( const Variable<array_1d<double,3> >& rVariable,
                                           std::vector<array_1d<double,3> > rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
                                           std::vector<Vector>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            virtual void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, 
                                           std::vector<Matrix>& rValues, 
                                           const ProcessInfo& rCurrentProcessInfo){}
            
            /**
             * returns the used integration method. In the general case this is the 
             * default integration method of the used geometry. I an other integration
             * method is used the method has to be overwritten within the element
             * @return default integration method of the used Geometry
             */
            virtual IntegrationMethod GetIntegrationMethod()
            {
                return mpGeometry->GetDefaultIntegrationMethod();
            }
            
            ///@}
            ///@name Elemental Data
            ///@{
            
            DataValueContainer& Data()
            {
                return mData;
            }
            
            template<class TVariableType> typename TVariableType::Type& GetValue(
                    const TVariableType& rThisVariable)
            {
                return mData.GetValue(rThisVariable);
            }
            template<class TVariableType> void SetValue( 
                    const TVariableType& rThisVariable, 
                    typename TVariableType::Type const& rValue)
            {
                mData.SetValue(rThisVariable, rValue);
            }
            
            template<class TVariableType> typename TVariableType::Type const& GetValue( 
                    const TVariableType& rThisVariable) const
            {
                return mData.GetValue(rThisVariable);
            }
            
            template<class TDataType> bool Has( const Variable<TDataType>& rThisVariable) const
            {
                return mData.Has(rThisVariable);
            }
            
            template<class TAdaptorType> bool Has( 
                    const VariableComponent<TAdaptorType>& rThisVariable) const
            {
                return mData.Has(rThisVariable);
            }
            
            ///@}
            ///@name Access
            ///@{
            
            GeometryType::Pointer pGetGeometry()
            {
                return mpGeometry;
            }
            
            const GeometryType::Pointer pGetGeometry() const
            {
                return mpGeometry;
            }
            
            GeometryType& GetGeometry()
            {
                return *mpGeometry;
            }
            
            GeometryType const& GetGeometry() const
            {
                return *mpGeometry;
            }
            
            PropertiesType::Pointer pGetProperties()
            {
                return mpProperties;
            }
            
            const PropertiesType::Pointer pGetProperties() const
            {
                return mpProperties;
            }
            
            PropertiesType& GetProperties()
            {
                return *mpProperties;
            }
            
            PropertiesType const& GetProperties() const
            {
                return *mpProperties;
            }
            
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
                buffer << "Element #" << Id();
                return buffer.str();
            }
            
            /// Print information about this object.
            virtual void PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "Element #" << Id();
            }
            
            /// Print object's data.
            virtual void PrintData(std::ostream& rOStream) const
            {
                mpGeometry->PrintData(rOStream);
            }
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
            /** A pointer to data related to this node. */
            DataValueContainer mData;
            
            /**
             * pointer to the element's geometry
             */
            GeometryType::Pointer mpGeometry;
            
            /**
             * pointer to the element's properties
             */
            Properties::Pointer mpProperties;
            
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
            ///@}    
    }; // Class Element 
    ///@} 
    ///@name Type Definitions       
    ///@{ 
    ///@} 
    ///@name Input and output 
    ///@{ 
    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, 
                                      Element& rThis);
    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const Element& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << " : " << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }
    ///@} 
    
    /**
     * definition of NEIGHBOUR_ELEMENTS variable
     */
    KRATOS_DEFINE_VARIABLE(WeakPointerVector< Element >, NEIGHBOUR_ELEMENTS)

}  // namespace Kratos.
#endif // KRATOS_ELEMENT_H_INCLUDED  defined 
