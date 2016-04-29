//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-09-25 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//

#if !defined(KRATOS_LINECLUSTER3D_H_INCLUDED )
#define  KRATOS_LINECLUSTER3D_H_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <cmath>

// Project includes
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/node.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "utilities/indexed_object.h"
#include "containers/weak_pointer_vector.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_elements/cluster3D.h"

namespace Kratos
{
    class Element;
    class ProcessInfo;
    
    class LineCluster3D : public Cluster3D
    {
    public:

      /// Pointer definition of Cluster3D
        KRATOS_CLASS_POINTER_DEFINITION(LineCluster3D);

        LineCluster3D( );
        LineCluster3D( IndexType NewId, GeometryType::Pointer pGeometry );
        LineCluster3D( IndexType NewId, NodesArrayType const& ThisNodes);
        LineCluster3D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;      
        
      /// Destructor.
        virtual ~LineCluster3D();
   
        virtual void CustomInitialize(ProcessInfo& r_process_info);          

        virtual std::string Info() const
        {
	    std::stringstream buffer;
	    buffer << "Discrete Element #" << Id();
	    return buffer.str();
        }
      
        virtual void PrintInfo(std::ostream& rOStream) const
        {
	    rOStream << "Discrete Element #" << Id();
        }
      
        virtual void PrintData(std::ostream& rOStream) const
        {
	  //mpGeometry->PrintData(rOStream);
        }
      
    protected:
           
    private:
      

    }; // Class LineCluster3D

 
  /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, 
                    LineCluster3D& rThis){
        return rIStream;
    }

  /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, 
                    const LineCluster3D& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
}  // namespace Kratos.

#endif // KRATOS_LINECLUSTER3D_INCLUDED  defined
