

// Forward Declarations
template<typename... T> struct TypeTraits;
template<typename DERIVED_T> class BASE;
template<class T> class DERIVED;


// Recusive expression downcasting to derived type
/*
template<typename TYPE>
struct TypeTraits<TYPE>
{
    typedef typename TypeTraits<TYPE>::dtype    dtype;
};
*/


template<typename DERIVED_T>
class BASE
{
  public:  
    typedef typename TypeTraits<DERIVED_T>::dtype dtype;
};


template<class T>
class DERIVED : public BASE< DERIVED<T> >
{
  public:
    typedef T dtype;  
};





int main()
{
  return 0;
}