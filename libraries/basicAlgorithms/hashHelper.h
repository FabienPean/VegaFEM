// The following code is from Boost library.
// Boost has the following license text:
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)
#ifndef HASHHELPER_H
#define HASHHELPER_H
#include <cstddef>

inline std::size_t hashCombine(std::size_t seed, std::size_t newValue)
{
    seed ^= newValue + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
}


#endif
