//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#if !defined(CONSTRAINT_H)
#define CONSTRAINT_H
// System includes

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "containers/constraint_equation.h"
#include "containers/flags.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{
/** \brief Constraint
	* A class that implements the interface for different constraints to be applied on a system.
    */
class MasterSlaveConstraint :  public IndexedObject, public Flags
{

  public:
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveConstraint);

    typedef Dof<double> DofType;
    typedef Node<3> NodeType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef double ConstantType;
    typedef Matrix MatrixType;

    ///@name Life Cycle
    ///@{

    /**
	    Creates a MultipointConstraint object
	*/
    MasterSlaveConstraint() : Flags()
    {
    }
    /// Destructor.
    virtual ~MasterSlaveConstraint(){};        


    ///@}

    ///@name Access
    ///@{

    /**
	* Clears the maps contents
	*/
    void Clear()
    {
        mConstraintEquation.Clear();
    }


    /**
     * is called to initialize the constraint
     * if the constraint needs to perform any operation before any calculation is done
     */
    virtual void Initialize()
    {
    }

    /**
     * this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void EquationIdVector(EquationIdVectorType& rResult,
                                  ProcessInfo& rCurrentProcessInfo)
    {
        if (rResult.size() != 0)
            rResult.resize(0);
    }

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rTransformationMatrix the elemental left hand side matrix
     * @param rConstant the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLocalSystem(MatrixType& rTransformationMatrix,
                                      ConstantType& rConstant,
                                      ProcessInfo& rCurrentProcessInfo)
    {
      if (rTransformationMatrix.size1() != 0)
      {
    	rTransformationMatrix.resize(0, 0, false);
      }

	  rConstant = 0.00;
    }


    /**
	* Adds a master to the current master slave relation
	*/
    virtual void AddMaster(DofType const &rMasterDof, double Weight)
    {

    }

    /**
	* Adds a master to the current master slave relation
	*/
    virtual void AddSlave(DofType const &rSlaveDof, double Weight)
    {

    }    


    /**
	* Get the Total number of MasterDOFs for a given slave dof
	* @return Total number of MasterDOFs for a given slave dof
	*/
    virtual std::size_t GetNumberOfMasters()
    {
        return mConstraintEquation.GetNumberOfMasters();
    }

    /**
	* Returns the string containing a detailed description of this object.
	* @return the string with informations
	*/
    virtual void GetInfo() const
    {
    }

    ///@

    ///@name Static Operations
    ///
    //@{

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " Constraint base class !" << std::endl;
        mConstraintEquation.PrintInfo(rOStream);
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {
        this->mConstraintEquation.save(rSerializer);
    }

    virtual void load(Serializer &rSerializer) override
    {
        this->mConstraintEquation.load(rSerializer);
    }

  private:

    MasterSlaveRelation mConstraintEquation;
    ///@}

    ///@}
};



///@name Input/Output funcitons
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                MasterSlaveConstraint& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                const MasterSlaveConstraint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;

    return rOStream;
}
    

///@}



} // namespace Kratos

#endif // CONSTRAINT_H_INCLUDED