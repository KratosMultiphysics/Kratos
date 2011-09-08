/*
  class COMPLEX
  =============
  A minimal class used when passing complex arithmetic variables
  from C++ to FORTRAN 77.

  The template parameter is used for specification of precision:

  COMPLEX<float>  is equivalent to F77 COMPLEX
  COMPLEX<double> is equivalent to F77 DOUBLE COMPLEX

  Author: Carsten A. Arnholm,
  Updates:
      04-MAR-1996 initial, non-template version
      14-MAY-1996 Template version
      29-JUL-1996 Tested portability to SGI/Unix,
                  corrected operator=(const COMPLEX<T>& )
*/

#ifdef real
   // some people define real as a macro
   #undef real
   #pragma message(__FILE__" : warning: 'real' macro definition cancelled")
#endif

template<class T>
class COMPLEX {
public:
   COMPLEX();
   COMPLEX(const COMPLEX<T>& );
   COMPLEX(const T& re,const T& im);
   COMPLEX<T>& operator=(const COMPLEX<T>& );
  ~COMPLEX();
   const T& real();
   const T& imag();
private:
   T m_re;
   T m_im;
};

template<class T>
inline COMPLEX<T>::COMPLEX()
:m_re(T()),m_im(T())
{}

template<class T>
inline COMPLEX<T>::COMPLEX(const COMPLEX<T>& copy)
:m_re(copy.m_re),m_im(copy.m_im)
{}

template<class T>
inline COMPLEX<T>::COMPLEX(const T& re,const T& im)
:m_re(re),m_im(im)
{}

template<class T>
inline COMPLEX<T>& COMPLEX<T>::operator=(const COMPLEX<T>& copy)
{
   m_re = copy.m_re;
   m_im = copy.m_im;
   return *this;
}

template<class T>
inline COMPLEX<T>::~COMPLEX()
{}

template<class T>
inline const T& COMPLEX<T>::real()
{
   return m_re;
}

template<class T>
inline const T& COMPLEX<T>::imag()
{
   return m_im;
}
