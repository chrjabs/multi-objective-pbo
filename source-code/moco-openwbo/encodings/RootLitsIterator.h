#ifndef ROOTLITSITERATOR_H
#define ROOTLITSITERATOR_H
#include "../partition.h"
#include "RootLitsTypes.h"
#include <iterator>
namespace rootLits
{
  enum class iterators {
    GeQ,
    Strat
  };

  template<class pointer_out = myType::iterator, 
	   class value_type_out = myType::iterator::value_type>  class GIterator{
  public:
  using iterator_category	= std::bidirectional_iterator_tag;
    using difference_type	= myType::iterator::difference_type;
    using reference		= myType::reference;
    using pointer = pointer_out;
    using value_type = value_type_out;
  };

  template<class pointer_out = myType::iterator, 
	   class value_type_out = myType::iterator::value_type>  class IteratorT: 
    GIterator<pointer_out, value_type_out>{
    
  public:
    using value_type = typename GIterator<pointer_out, value_type_out>::value_type;
    using pointer = typename GIterator<pointer_out, value_type_out>::pointer;
  
    IteratorT(const pointer& it) : m_ptr{it} {}
    value_type operator*(){return *m_ptr;};
    pointer operator->() { return m_ptr; }
    IteratorT& operator++(){++m_ptr; return *this;};
    IteratorT operator++(int){IteratorT tmp = *this; ++m_ptr; return tmp;};
    IteratorT& operator--(){--m_ptr; return *this;}
    friend bool operator== (const IteratorT& a, const IteratorT& b) { return a.m_ptr == b.m_ptr; };
    friend bool operator!= (const IteratorT& a, const IteratorT& b) { return a.m_ptr != b.m_ptr; };  

  private:
    pointer m_ptr;
  };

  using Iterator =  IteratorT<>;
  using CIterator = IteratorT<myType::const_iterator, myType::const_iterator::value_type>;
  using GeQIterator = Iterator;
  using StratIterator = Iterator;
}

#endif
