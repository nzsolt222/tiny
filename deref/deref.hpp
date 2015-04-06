// The MIT License (MIT)
//
// Copyright (c) 2015 Nagy Zsolt
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef NZS_DEREF_HPP
#define NZS_DEREF_HPP

#include <type_traits>

namespace nzs
{
template <class T>
class is_dereferancable
{
    using Yes = char;
    using No = double;

    template <typename C>
    static No test(...);

    template <typename C, typename std::enable_if<
                              std::is_pointer<C>::value>::type* = nullptr>
    static Yes test(int);

    template <typename C,
              typename std::enable_if<std::is_class<C>::value>::type* = nullptr,
              typename std::enable_if<std::is_member_function_pointer<
                  decltype(&C::operator*)>::value>::type* = nullptr>
    static Yes test(int);

public:
    static bool const value = sizeof(test<T>(0)) == sizeof(Yes);
};

template <class T, class U = void>
struct dereference
{
    using type = T&;

    inline dereference(T& t) : t(t){};

    inline operator type() const noexcept { return t; }

private:
    type t;
};

template <class T>
struct dereference<T,
                   typename std::enable_if<is_dereferancable<T>::value>::type>
{
    using type = decltype(*std::declval<T>());

    inline dereference(T& t) : t(*t){};

    inline operator type() const noexcept { return t; }

private:
    type t;
};

template <class T>
typename dereference<T>::type deref(T& t)
{
    return dereference<T>(t);
}

} /* nzs */

#endif
