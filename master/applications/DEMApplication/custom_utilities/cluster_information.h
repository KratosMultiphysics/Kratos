//
// File:   cluster_information.h
// Author: Miguel Angel Celigueta maceli@cimne.upc.edu
//


#if !defined(CLUSTER_INFORMATION_H)
#define  CLUSTER_INFORMATION_H

namespace Kratos
{

class ClusterInformation
{
    
public:
    ClusterInformation(){};
    std::string mName;
    double mSize;
    double mVolume;
    std::vector<double> mListOfRadii;
    std::vector<array_1d<double,3> > mListOfCoordinates;
    array_1d<double,3> mInertias;
    
    virtual ~ClusterInformation() {};
    
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        //rOStream << std::endl << this->mSize <<"  " << this->mVolume << std::endl;                                
    }
    
    virtual std::string Info() const {
        std::stringstream buffer;
        buffer << mName;
        return buffer.str();
    }  
    
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
         //TODO:   
    }

    void load(Serializer& rSerializer)
    {   
        //TODO:            
    }		

};

inline std::istream& operator >> (std::istream& rIStream, ClusterInformation& rThis);

inline std::ostream& operator << (std::ostream& rOStream, const ClusterInformation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : ";
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // CLUSTER_INFORMATION_H  defined 

