#ifndef CONTROL_POINT_H
#define CONTROL_POINT_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ==============================================================================

namespace Kratos
{

//TODO: use kratos Nodes for the control points. you can pass the weight using the set value function of the node -- that is -- remove this class
class ControlPoint
{
public:

    // ==========================================================================
    // Type definitions for linear algebra including sparse systems
    // ==========================================================================
    //typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //typedef typename SparseSpaceType::MatrixType SparseMatrixType;
    //typedef typename SparseSpaceType::VectorType VectorType;
    //typedef std::vector<double> DoubleVector;
  typedef std::vector<double> DoubleVector;
    /// Pointer definition of ControlPoint
    //KRATOS_CLASS_POINTER_DEFINITION(ControlPoint);

    /// Default constructor.
    ControlPoint( double X, double Y, double Z, double W, unsigned int cp_id)
    {
      m_coordinates.push_back( X );
      m_coordinates.push_back( Y );
      m_coordinates.push_back( Z );
     //   m_displacements.push_back( 0.0 );
      //m_displacements.push_back( 0.0 );
      //m_displacements.push_back( 0.0 );
    m_w = W;
        m_cp_id = cp_id;
    }

    double getX()
    {
        //return (m_coordinates[0] + m_displacements[0]);
    return m_coordinates[0];
    }

    double getY()
    {
       //return (m_coordinates[1] + m_displacements[1]);
    return m_coordinates[1];
    }

    double getZ()
    {
      //return (m_coordinates[2] + m_displacements[2]);
    return m_coordinates[2];
    }

    //double getX0()
    //{
    //  return m_coordinates[0];
    //}
    //double getY0()
    //{
    //   return m_coordinates[1];
    //}
    //double getZ0()
    //{
    //  return m_coordinates[2];
    //}   
    //double getdX()
    //{
    //  return m_displacements[0];
    //}
    //double getdY()
    //{
    //   return m_displacements[1];
    //}
    //double getdZ()
    //{
    //  return m_displacements[2];
    //}  
    //void setdX(double dX)
    //{
    //    m_displacements[0] = dX;
    //}  
    //void setdY(double dY)
    //{
    //    m_displacements[1] = dY;
    //}  
    //void setdZ(double dZ)
    //{
    //    m_displacements[2] = dZ;
    //}          

    double getWeight()
    {
      return m_w;
    }

    int getCPId()
    {
      return m_cp_id;
    }

    //void SetMappingMatrixId(unsigned int id)
    //{
    //    m_mapping_matrix_id = id;
    //}    
    //int GetMappingMatrixId()
    //{
    //    if(m_mapping_matrix_id<0)
    //        KRATOS_THROW_ERROR(std::logic_error, "No mapping matrix ID specified for current control point", m_mapping_matrix_id);
    //  return m_mapping_matrix_id;
    //}
    //void SetRelevantForMapping()
    //{
    //  m_is_relevant_for_mapping = true;
    //}   
    //bool IsRelevantForMapping()
    //{
    //  return m_is_relevant_for_mapping;
    //}    

    /// Destructor.
    virtual ~ControlPoint()
    {
    }
    // ==============================================================================


//    evaluat

    // ==============================================================================


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ControlPoint";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ControlPoint";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


private:
    // ==============================================================================
    // General working arrays
    // ==============================================================================
    DoubleVector m_coordinates;
    //DoubleVector m_displacements;
    double m_w;
    unsigned int m_cp_id;
    //int m_mapping_matrix_id = -1;
    //bool m_is_relevant_for_mapping = false;


}; // Class ControlPoint

}  // namespace Kratos.

#endif // CONTROL_POINT_H
