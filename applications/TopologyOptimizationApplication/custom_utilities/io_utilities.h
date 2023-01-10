//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
//
// ==============================================================================

#if !defined(KRATOS_IO_UTILITIES_H_INCLUDED)
#define  KRATOS_IO_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes


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

/// Solution utility to filter results.
/** Detail class definition.

 */

class IOUtilities
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IOUtilities
    KRATOS_CLASS_POINTER_DEFINITION(IOUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IOUtilities( )
    {
    }

    /// Destructor.
    virtual ~IOUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- Save optimization results in a restart file (mdpa) --------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    void SaveOptimizationResults( const char* RestartInputFile, ModelPart& rModelPart, const char* RestartOutputFile )
    {
        KRATOS_TRY;

        KRATOS_INFO("[TopOpt]") << "::[Saving optimization results as restart file]::"<< std::endl;

        // Create an empty .mdpa restart file
        std::ofstream FileToBeCreated;
        FileToBeCreated.open(RestartOutputFile, std::ios::trunc);

        // Read the original .mdpa file
        std::ifstream FileToBeRead(RestartInputFile);
        if(!FileToBeRead.is_open())
            KRATOS_THROW_ERROR(std::invalid_argument, "Specified restart input file does not exist: ",RestartInputFile);

        // Write the given input file except the block covering X_PHYS (this is to be replaced by the current optimization results)
        bool DensityBlockActive = false;
        std::string LineString;
        while (std::getline(FileToBeRead, LineString))
        {
            if(LineString.find("Begin ElementalData X_PHYS")!=std::string::npos)
                DensityBlockActive = true;
            else if(LineString.find("End ElementalData")!=std::string::npos and DensityBlockActive)
            {
                DensityBlockActive = false;
                continue;
            }
            else if(DensityBlockActive == false)
                FileToBeCreated << LineString << "\n";
        }

        // Write the actual X_PHYS elemental data
        FileToBeCreated << "\nBegin ElementalData X_PHYS\n";
        for(ModelPart::ElementsContainerType::iterator elem_i = rModelPart.ElementsBegin(); elem_i!=rModelPart.ElementsEnd(); elem_i++)
        {
            FileToBeCreated << "    " << elem_i->Id() << " " << elem_i->GetValue(X_PHYS) << "\n";
        }
        FileToBeCreated << "End ElementalData\n";

        // Close files
        FileToBeCreated.close();
        FileToBeRead.close();

        KRATOS_INFO("[TopOpt]") <<"  Restart File succesfully generated under the name " << RestartOutputFile <<std::endl;

        KRATOS_CATCH("");
    }

    // ---------------------------------------------------------------------------------------------------------------------------------------------
    // --------------------------------- WRITE STL FILE FROM SURFACE MESH  -------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------

    /// Generates a .STL file from a a provided surface mesh
    void WriteSurfaceAsSTLFile(const char* file_name, ModelPart& rSurfaceModelPart)
    {
        KRATOS_TRY;

        KRATOS_INFO("[TopOpt]") <<"\n::[Generating STL file]::"<<std::endl;

        // Write stl of surface model part
        std::ofstream myfile;
        myfile.open (file_name);
        myfile << "solid Layer0\n";
        for ( ModelPart::ConditionIterator cond_i = rSurfaceModelPart.ConditionsBegin(); cond_i != rSurfaceModelPart.ConditionsEnd(); ++cond_i )
        {
            array_1d<double,3> area_normal = cond_i->GetValue(NORMAL);
            myfile << "  facet normal " << area_normal[0] <<" " << area_normal[1] << " " << area_normal[2] <<"\n";
            myfile << "    outer loop\n";

            Element::GeometryType& rNodes = cond_i->GetGeometry();
            if (rNodes.size()==3)
            {
                for(unsigned int i = 0; i<rNodes.size(); i++)
                    myfile << "      vertex "<< rNodes[i].X() <<" " << rNodes[i].Y() << " " << rNodes[i].Z() <<"\n";

                myfile << "    end loop\n";
                myfile << "  end facet\n";
            }
            else 
            {
                KRATOS_ERROR << "The number of nodes is: "<< rNodes.size() << ". This is not feasible for STL files!"<< std::endl;
            }

        }
        myfile << "endsolid Layer0\n";
        myfile.close();

        KRATOS_INFO("[TopOpt]") <<"  STL File succesfully generated under the name " << file_name <<std::endl;

        KRATOS_CATCH("");
    }

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
        return "IOUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "IOUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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
    //IOUtilities& operator=(IOUtilities const& rOther);

    /// Copy constructor.
    //IOUtilities(IOUtilities const& rOther);


    ///@}

}; // Class IOUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif	/* KRATOS_IO_UTILITIES_H_INCLUDED */
