// file sgnint32_t.cpp
#include "sgnint32_t.hpp"

sgnint32_t& sgnint32_t :: operator=(const sgnint32_t& i)
{
   this->val=i.Val();
   this->sgn=i.Sgn();
   
   return *this;
}

sgnint32_t operator-(const sgnint32_t& i)
{
   sgnint32_t ret(i.Val(),-i.Sgn());
   
   return ret;
}

bool operator==(const sgnint32_t a, const sgnint32_t b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return true;
   else
      return false;
}

bool operator!=(const sgnint32_t a, const sgnint32_t b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return false;
   else
      return true;
}

bool operator>(const sgnint32_t a, const sgnint32_t b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval > bval )
      return true;
   else
      return false;
}

bool operator<(const sgnint32_t a, const sgnint32_t b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval < bval )
      return true;
   else
      return false;
}

bool operator==(const sgnint32_t a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return true;
   else
      return false;
}

bool operator!=(const sgnint32_t a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return false;
   else
      return true;
}

bool operator<(const sgnint32_t a, const int b)
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

bool operator<(const int b, const sgnint32_t a )
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return false;
   else if (asgn < 0)
      return (b<(-aval));
   else
      return (b<aval);
}

// bool operator<(const sgnint32_t a, int b)
// {
   // auto aval= a.Val();
   // auto asgn= a.Sgn();
   
   // if ( b == 0 && asgn < 0)
      // return true;
   // else if (asgn < 0)
      // return ((-aval)<b);
   // else
      // return (aval<b);
// }

// bool operator<(const sgnint32_t a, const unsigned b)
// {
   // auto aval= a.Val();
   // auto asgn= a.Sgn();
   
   // if ( b == 0 && asgn < 0)
      // return true;
   // else if (asgn < 0)
      // return ((-aval)<b);
   // else
      // return (aval<b);
// }

bool operator<(const sgnint32_t a, const double b)
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

uint32_t abs(const sgnint32_t a)
{   
   return a.Val(); 
}

bool signbit(const sgnint32_t& a)
{
   return a.Sgn()<0; 
}

int32_t sgn(const sgnint32_t& a)
{
   return a.Sgn();
}
