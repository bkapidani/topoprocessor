/*
 * This source file is part of Topoprocessor.
 *
 * Copyright (C) 2021, Bernard Kapidani bernard(dot)kapidani(at)gmail(dot)com
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * SEE THE GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS.
 * YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
 * ALONG WITH THIS PROGRAM. IF NOT, SEE <https://www.gnu.org/licenses/>.
 */
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
