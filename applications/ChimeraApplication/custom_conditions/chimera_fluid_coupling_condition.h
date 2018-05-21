//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Sonja Schneider
//                   Jordi Cotela
//

// From S. Schneider Implementation of a Chimera Technique TUM Master Thesis 2015

#ifndef KRATOS_CHIMERA_FLUID_COUPLING_CONDITION_H
#define KRATOS_CHIMERA_FLUID_COUPLING_CONDITION_H

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "geometries/geometry_data.h"
#include "includes/condition.h"
#include "includes/process_info.h"

// Application includes


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

/// A condition to impose continuity of tractions across fluid boundaries. To be used with VMS2D or VMS3D elements.
template< unsigned int TDim >
class ChimeraFluidCouplingCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ChimeraFluidCouplingCondition
    KRATOS_CLASS_POINTER_DEFINITION(ChimeraFluidCouplingCondition);

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

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    ChimeraFluidCouplingCondition(IndexType NewId = 0);

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    ChimeraFluidCouplingCondition(IndexType NewId,
                            const NodesArrayType& ThisNodes);

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    ChimeraFluidCouplingCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry);

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    ChimeraFluidCouplingCondition(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties);

    /// Copy constructor.
    ChimeraFluidCouplingCondition(ChimeraFluidCouplingCondition const& rOther);

    /// Destructor.
    virtual ~ChimeraFluidCouplingCondition();


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    ChimeraFluidCouplingCondition & operator=(ChimeraFluidCouplingCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new ChimeraFluidCouplingCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    virtual Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;


    /// Allocate memory for integration point data.
    virtual void Initialize() override;

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override;

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;


    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override;




    /// Check that all data required by this condition is available and reasonable
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) override;


    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo) override;


    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofsVectorType& ConditionDofList,
                            ProcessInfo& CurrentProcessInfo) override;

    virtual GeometryData::IntegrationMethod GetIntegrationMethod() override;


    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z,) PRESSURE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    virtual void GetValuesVector(Vector& Values,
                                 int Step = 0) override;


    virtual void Calculate(const Variable< array_1d<double,3> >& rVariable,
                           array_1d<double,3>& Output,
                           const ProcessInfo& rCurrentProcessInfo) override;


    /// Used to pass interpolated tensions on the integration point (with rVariable == TENSION)
    virtual void SetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                             std::vector<array_1d<double, 3 > > rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                             std::vector<array_1d<double, 3 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;



    virtual void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                             std::vector<double>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


    virtual void GetValueOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
                                             std::vector<array_1d<double, 6 > >& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    virtual void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
                                             std::vector<Vector>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;


    virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                             std::vector<Matrix>& rValues,
                                             const ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                              std::vector<double>& rOutput,
                                              const ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                                              std::vector< array_1d<double, 3 > >& Output,
                                              const ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
                                              std::vector< Vector >& Output,
                                              const ProcessInfo& rCurrentProcessInfo) override;

    virtual void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                                              std::vector< Matrix >& Output,
                                              const ProcessInfo& rCurrentProcessInfo) override;


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
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;


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

    void CalculateNormal(array_1d<double,3>& An );


    void ImposeTraction(MatrixType &rLocalMatrix, VectorType &rLocalVector);


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

    std::vector< array_1d<double,3> > mGPTraction;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;

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

}; // Class ChimeraFluidCouplingCondition


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< unsigned int TDim >
inline std::istream& operator >> (std::istream& rIStream,
                                  ChimeraFluidCouplingCondition<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template< unsigned int TDim >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ChimeraFluidCouplingCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_CHIMERA_FLUID_COUPLING_CONDITION_H
