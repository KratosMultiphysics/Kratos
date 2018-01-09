//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//                    
//


#if !defined(KRATOS_CONTACT_PAIR_INCLUDED )
#define  KRATOS_CONTACT_PAIR_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <cmath>


template <class T>
class ContactPair {
public:
  T value[2];

  ContactPair() {}

  ContactPair(const T& First,const T& Second) {
    value[0]  = First;
    value[1]  = Second;
  }

  ~ContactPair() {}

  T& operator[](std::size_t index) {
    return value[index];
  }

  T const& operator[](std::size_t index) const {
    return value[index];
  }

  ContactPair& operator = (ContactPair& Pair) {
    value[0] = Pair[0];
    value[1] = Pair[1];

    return *this;
  }

  ContactPair& operator = (const ContactPair& Pair) {
    value[0] = Pair[0];
    value[1] = Pair[1];

    return *this;
  }


  inline bool operator == (const ContactPair& Pair){
    return  (value[0] == Pair[0]) && (value[1] == Pair[1]) ;
  }

  inline bool  operator == (const ContactPair& Pair) const {
    return ( (value[0] == Pair[0]) && (value[1] == Pair[1]) ) ;
  }

  inline std::ostream& operator << (std::ostream& OStream) {
    OStream << "Object 1 =  "<<  value[0] << "  " << "Object 2 =" << value[1] << std::endl;
    return OStream;
  }

  inline std::ostream& operator << (std::ostream& OStream) const {
    OStream << "Object 1 =  "<<  value[0] << "  " << "Object 2 =" << value[1] << std::endl;
    return OStream;
  }
};

#endif // KRATOS_CONTACT_PAIR_INCLUDED
