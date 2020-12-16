
#include <stdexcept>
#include <string>

#pragma once

#ifndef _is_likely
  //asserts should always pass, help the jump predictions with compiler builtins
  #if defined(__GNUC__) || defined(__clang__)
    #define _is_likely(expr) ( __builtin_expect( (!!(expr)), 1) )
  #else // MSVC __assume is NOT the equivalent
    #define _is_likely(expr) (expr)
  #endif
#endif

class AssertFailed : public std::runtime_error {
public:
  explicit AssertFailed(const std::string& what)
    : std::runtime_error(what)
  {}
};

#define RIMACS_ASSERT(EXPR) do { if(! _is_likely(EXPR)) { throw AssertFailed(#EXPR) ; } } while(false)

