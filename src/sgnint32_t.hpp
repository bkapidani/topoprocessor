// file sgnint32_t<T>.hpp
#ifndef SGNINT32_T_HPP
#define SGNINT32_T_HPP

#include <cmath>
#include <iostream>
#include <stdlib.h>

template <typename T>
class sgnint32_t
{
   public:
      uint32_t val;
      T sign;
      sgnint32_t<T>(void) { this->val=0; this->sign=1; }
      sgnint32_t<T>(uint32_t v, T s) { this->val=v; this->sign=s; }
      uint32_t Val() const { return this->val; }
      T Sign()  const { return this->sign; } 
      sgnint32_t<T>& operator=(const sgnint32_t<T>& i);
      sgnint32_t<T>& operator-(const sgnint32_t<T>& i);
};

// sgnint32_t<T> operator-(const sgnint32_t<T>& i);
// bool operator==(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator!=(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator>(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator<(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator<(const int b, const sgnint32_t<T> a);
// bool operator<(const sgnint32_t<T> a, const int b);
// bool operator==(const sgnint32_t<T> a, const int b);
// bool operator!=(const sgnint32_t<T> a, const int b);
// bool operator<(const sgnint32_t<T> a, const double b);
// uint32_t abs(const sgnint32_t<T> a);
// bool signbit(const sgnint32_t<T>& a);
// int32_t Sign(const sgnint32_t<T>& a);

template <typename T>
sgnint32_t<T>& sgnint32_t<T> :: operator=(const sgnint32_t<T>& i)
{
   this->val=i.Val();
   this->sign=i.Sign();
   
   return *this;
}

template <typename T>
sgnint32_t<T> operator-(const sgnint32_t<T>& i)
{
   sgnint32_t<T> ret(i.Val(),-i.Sign());
   
   return ret;
}

template <typename T>
bool operator==(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator!=(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return false;
   else
      return true;
}

template <typename T>
bool operator>(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval > bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval < bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator==(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return true;
   else
      return false;
}

template <typename T>
bool operator!=(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return false;
   else
      return true;
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   auto asgn= a.Sign();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

template <typename T>
bool operator<(const int b, const sgnint32_t<T> a )
{
   auto aval= a.Val();
   auto asgn= a.Sign();
   
   if ( b == 0 && asgn < 0)
      return false;
   else if (asgn < 0)
      return (b<(-aval));
   else
      return (b<aval);
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const double b)
{
   auto aval= a.Val();
   auto asgn= a.Sign();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

template <typename T>
uint32_t abs(const sgnint32_t<T> a)
{   
   return a.Val(); 
}

template <typename T>
bool signbit(const sgnint32_t<T>& a)
{
   return a.Sign()<0; 
}

template <typename T>
T Sign(const sgnint32_t<T>& a)
{
   return a.Sign();
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const sgnint32_t<T>& s)
{
    T d = s.Sign()*(s.Val());
    os << d;
    
    return os;
}

#endif

