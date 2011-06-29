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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-04-24 10:30:22 $
//   Revision:            $Revision: 1.6 $
//
//


#if !defined(KRATOS_CONDITION_H_INCLUDED )
#define  KRATOS_CONDITION_H_INCLUDED



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

    /// Short class definition.

    /** Detail class definition.
     */
    class Condition : public IndexedObject
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of Condition
        KRATOS_CLASS_POINTER_DEFINITION(Condition);

        typedef IndexedObject BaseType;

        typedef Node < 3 > NodeType;

        typedef Properties PropertiesType;

        typedef Geometry<NodeType> GeometryType;

        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

        typedef Vector VectorType;

        typedef Matrix MatrixType;

        typedef std::size_t IndexType;

        typedef std::size_t SizeType;

        typedef std::vector<std::size_t> EquationIdVectorType;

        typedef std::vector< Dof<double>::Pointer > DofsVectorType;

        typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

        typedef VectorMap<IndexType, DataValueContainer> SolutionStepsConditionalDataContainerType;


        ///@}
        ///@name Life Cycle
        ///@{

        /** Constructor.
         */
        Condition(IndexType NewId = 0)
        : BaseType(NewId)
        , mpGeometry()
//        , mpGeometry(new GeometryType(NodesArrayType()))
        , mpProperties(new PropertiesType)
        {
        }

        /** Constructor using an array of nodes
         */
        Condition(IndexType NewId, const NodesArrayType& ThisNodes)
        : BaseType(NewId)
        , mpGeometry(new GeometryType(ThisNodes))
        , mpProperties(new PropertiesType)
        {
        }

        /** Constructor using Geometry
         */
        Condition(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId)
        , mpGeometry(pGeometry)
        , mpProperties(new PropertiesType)
        {
        }

        /** Constructor using Properties
         */
        Condition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId)
        , mpGeometry(pGeometry)
        , mpProperties(pProperties)
        {
        }

        /// Copy constructor.

        Condition(Condition const& rOther)
        : BaseType(rOther)
        , mpGeometry(rOther.mpGeometry)
        , mpProperties(rOther.mpProperties)
        {
        }


        /// Destructor.

        virtual ~Condition()
        {
        }


        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.

        Condition & operator=(Condition const& rOther)
        {
            BaseType::operator=(rOther);
            mpGeometry = rOther.mpGeometry;
            mpProperties = rOther.mpProperties;

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        virtual Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
        {
            return Condition::Pointer(new Condition(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            if (rLeftHandSideMatrix.size1() != 0)
                rLeftHandSideMatrix.resize(0, 0);
            if (rRightHandSideVector.size() != 0)
                rRightHandSideVector.resize(0);
        }

        virtual void CalculateLocalVelocityContribution(MatrixType& rDampMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            if (rDampMatrix.size1() != 0)
                rDampMatrix.resize(0, 0);
            /*if(rRightHandSideVector.size()!=0)
                    rRightHandSideVector.resize(0);*/
        }

        virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            if (rLeftHandSideMatrix.size1() != 0)
                rLeftHandSideMatrix.resize(0, 0);
        }

        virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            if (rRightHandSideVector.size() != 0)
                rRightHandSideVector.resize(0);
        }

        virtual void CalculateLeftHandSide(VectorType& rLeftHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
            if (rLeftHandSideVector.size() != 0)
                rLeftHandSideVector.resize(0);
        }

        virtual void CalculateRightHandSide(MatrixType& rRightHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            if (rRightHandSideMatrix.size1() != 0)
                rRightHandSideMatrix.resize(0, 0);
        }

        virtual void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
        {
            if (rResult.size() != 0)
                rResult.resize(0);
        }

        virtual void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo)
        {
            if (ConditionDofList.size() != 0)
                ConditionDofList.resize(0);
        }

        virtual void MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            if (rMassMatrix.size1() != 0)
                rMassMatrix.resize(0, 0);
        }

        //adds the mass matrix scaled by a given factor to the LHS

        virtual void AddMassMatrix(MatrixType& rLeftHandSideMatrix, double coeff, ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
        {
            if (rDampMatrix.size1() != 0)
                rDampMatrix.resize(0, 0);
        }
        //adds the inertia forces to the RHS --> performs residua = static_residua - coeff*M*acc

        virtual void AddInertiaForces(VectorType& rRightHandSideVector, double coeff, ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void Initialize()
        {
        }

        virtual void InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
        {
        }

        virtual void FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
        {
        }

        virtual void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
        {
        }

        virtual void FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
        {
        }

        virtual void CleanMemory()
        {
        }

        virtual void GetValuesVector(Vector& values, int Step = 0)
        {
        }

        virtual void GetFirstDerivativesVector(Vector& values, int Step = 0)
        {
        }

        virtual void GetSecondDerivativesVector(Vector& values, int Step = 0)
        {
        }

        virtual void Calculate(const Variable<double >& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void Calculate(const Variable< array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        //calculate on gauss points

        virtual void CalculateOnIntegrationPoints(const Variable<double>& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable, std::vector< array_1d<double, 3 > >& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable, std::vector< Vector >& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        /**
         * Access for variables on Integration points.
         * This gives access to variables stored in the constitutive law on each integration point.
         * Specialisations of condition must specify the actual interface to the integration points!
         * Note, that these functions expect a std::vector of values for the specified variable type that
         * contains a value for each integration point!
         * SetValueOnIntegrationPoints: set the values for given Variable.
         * GetValueOnIntegrationPoints: get the values for given Variable.
         */
        virtual void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void SetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable, std::vector<array_1d<double, 3 > > rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void SetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void SetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable, std::vector<array_1d<double, 3 > >& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable, std::vector<array_1d<double, 6 > >& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }

        virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
        {
        }
        ///@}
        ///@name Conditional Data
        ///@{

        DataValueContainer& Data()
        {
            return mData;
        }

        template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
        {
            return mData.GetValue(rThisVariable);
        }

        template<class TVariableType>
        void SetValue(const TVariableType& rThisVariable, typename TVariableType::Type const& rValue)
        {
            mData.SetValue(rThisVariable, rValue);
        }

        template<class TVariableType> typename TVariableType::Type const& GetValue(const TVariableType& rThisVariable) const
        {
            return mData.GetValue(rThisVariable);
        }

        template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
        {
            return mData.Has(rThisVariable);
        }

        template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
        {
            return mData.Has(rThisVariable);
        }

        /**
         * This function provides the place to perform checks on the completeness of the input.
         * It is designed to be called only once (or anyway, not often) typically at the beginning
         * of the calculations, so to verify that nothing is missing from the input
         * or that no common error is found.
         * @param rCurrentProcessInfo
         */
        virtual int Check(const ProcessInfo& rCurrentProcessInfo)
        {
            KRATOS_TRY
            if (this->Id() < 1)
            {
                KRATOS_ERROR(std::logic_error, "Element found with Id 0 or negative","")
            }
            if (this->GetGeometry().Area() < 0)
            {
                std::cout << "error on element -> " << this->Id() << std::endl;
                KRATOS_ERROR(std::logic_error, "Area can not be lesser than 0","")
            }
            return 0;

            KRATOS_CATCH("");
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
            buffer << "Condition #" << Id();
            return buffer.str();
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "Condition #" << Id();
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

        GeometryType::Pointer mpGeometry;

        Properties::Pointer mpProperties;


      ///@}
      ///@name Serialization
      ///@{

	friend class Serializer;

	virtual void save(Serializer& rSerializer) const
	{
 	  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
	  rSerializer.save("Data", mData);
          rSerializer.save("Geometry",mpGeometry);
	  rSerializer.save("Properties", mpProperties);
	}

	virtual void load(Serializer& rSerializer)
	{
 	  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
	  rSerializer.load("Data", mData);
          rSerializer.load("Geometry",mpGeometry);
	  rSerializer.load("Properties", mpProperties);
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



        ///@}

    }; // Class Condition 

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream & operator >>(std::istream& rIStream,
            Condition& rThis);

    /// output stream function

    inline std::ostream & operator <<(std::ostream& rOStream,
            const Condition& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << " : " << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}
    KRATOS_DEFINE_VARIABLE(WeakPointerVector< Condition >, NEIGHBOUR_CONDITIONS)

} // namespace Kratos.

#endif // KRATOS_CONDITION_H_INCLUDED  defined 


