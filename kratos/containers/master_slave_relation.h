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

#if !defined(CONSTRAINT_EQUATION_H)
#define CONSTRAINT_EQUATION_H
// System includes
#include <vector>
#include <unordered_set>
#include <iostream>
#include <assert.h>

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/process_info.h"

namespace Kratos
{


/**
 * @class MasterSlaveRelation
 * @ingroup KratosCore
 * @brief This class stores the information regarding the MasterSlaveRelation equation. Naming convenction is defined like this. (each object of this class will store one equation in the given form
 * @details *   SlaveEquationId = w_1*MasterEquationId_1 + w_2*MasterEquationId_2 + ..... + w_n*MasterEquationId_n
 *
 *   each slaveDOF and each of its masterDOF have the following attributes
 *   a. dof ID
 *   b. dof KEY
 *   c. equation ID
 *
 *   one should be able to access the constraint equation both with (dofID , dofKEY) pair and/or equationID of the slave.
 *   This equation is imposed on the linear system of equations either element wise or on the global system.
 * @author Aditya Ghantasala
 */
class MasterSlaveRelation : public IndexedObject
{
  public:
    typedef std::size_t IndexType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef std::vector<IndexType> EquationIdVectorType;
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveRelation);

    // empty constructor and methods to add master and slave independently.
    MasterSlaveRelation() : IndexedObject(0)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }


    MasterSlaveRelation(IndexType const &rSlaveEquationId) : IndexedObject(rSlaveEquationId)
    {
        SetConstant(0.0);
        SetConstantUpdate(0.0);
    }

    void SetConstant(double Constant) { mConstant = Constant; }
    void SetConstantUpdate(double ConstantUpdate) { mConstantUpdate = ConstantUpdate; }
    double Constant() const { return mConstant; }
    double ConstantUpdate() const { return mConstantUpdate; }
    IndexType SlaveEquationId() const { return this->Id(); }



    // Get number of masters for this slave
    std::size_t GetNumberOfMasters() const
    {
        return mMasterDataSet.size();
    }


    /**
     * this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void EquationIdVector(EquationIdVectorType& rSlaveEquationIds,
                                  EquationIdVectorType& rMasterEquationIds,
                                  ProcessInfo& rCurrentProcessInfo)
    {
        if (rSlaveEquationIds.size() != 0)
            rSlaveEquationIds.resize(1);

        if (rMasterEquationIds.size() != 0)
            rMasterEquationIds.resize(this->GetNumberOfMasters());

        rSlaveEquationIds[0] = this->SlaveEquationId();

        auto master_it = mMasterDataSet.begin();
        for (IndexType i=0; i<this->GetNumberOfMasters(); i++)
            rMasterEquationIds[i] = (*(++master_it)).first;
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
                                      VectorType& rConstantVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
      if (rTransformationMatrix.size1() != 0)
      {
    	rTransformationMatrix.resize(1, this->GetNumberOfMasters(), false);
      }

      if (rConstantVector.size() != 0)
      {
    	rConstantVector.resize(1, false);
      }

        auto master_it = mMasterDataSet.begin();
        for (IndexType i=0; i<this->GetNumberOfMasters(); i++)
            rTransformationMatrix(0,i) = (*(++master_it)).second;

      rConstantVector(0) = Constant();
    }

    void PrintInfo(std::ostream& rOutput) const override
    {
        rOutput << "##############################" << std::endl;
        rOutput << "SlaveEquationId :: " << SlaveEquationId() << std::endl;
        rOutput << "Constant :: " << Constant() << std::endl;
        int index = 0;
        rOutput << "############################## :: Masters" << std::endl;
        for (auto &master : mMasterDataSet)
        {
            rOutput << index << " Master  equation id :: " << master.first << ", weight :: " << master.second << std::endl;
            index++;
        }
        rOutput << "##############################" << std::endl;
    }

    void Clear()
    {
        mMasterDataSet.clear();
    }

    // This is only for serializer. Not to be used outside
    void AddMaster(std::size_t MasterEquationId, double Weight)
    {
        auto res = mMasterDataSet.find(MasterEquationId);
        if (res != mMasterDataSet.end())
        {
            (*res).second += Weight;
        }
        else
        {
            mMasterDataSet[MasterEquationId] = Weight;
        }
    }

  private:
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {

    }

    virtual void load(Serializer &rSerializer) override
    {

    }

    ///@}


    std::unordered_map<IndexType, double> mMasterDataSet;

    double mConstant;
    double mConstantUpdate;

}; // End of ConstraintEquation class


} // namespace Kratos

#endif // CONSTRAINT_SLAVE_H_INCLUDED