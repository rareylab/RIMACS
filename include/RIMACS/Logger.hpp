
#pragma once


#include <iostream>

#ifndef _is_unlikely
  //asserts should always pass, help the jump predictions with compiler builtins
  #if defined(__GNUC__) || defined(__clang__)
    #define _is_unlikely(expr) ( __builtin_expect( (!!(expr)), 0) )
    #define pass(expr) expr
  #else // MSVC __assume is NOT the equivalent
    #define _is_unlikely(expr) (expr)
    #define pass(expr) expr
  #endif
#endif

namespace Base {

  enum LogLevel {
    NoLog = 0,
    Warning, // Warning is written to error log
    Info,
    Debug
  };


class UserLog {
public:
  const LogLevel m_level = LogLevel::Warning;
  std::ostream& cout = std::cout;
  std::ostream& cerr = std::cerr;
};

#define LOGGING_IS_ENABLED_(INFO, LEVEL, EXPECT) (INFO && EXPECT(INFO->m_level <= LEVEL))
#define LOGGING_IS_ENABLED(INFO, LEVEL) LOGGING_IS_ENABLED_(INFO, LEVEL, pass)

#define LOG_DEBUG(INFO, LOG) if(LOGGING_IS_ENABLED_(INFO, Base::LogLevel::Debug, _is_unlikely)) {INFO->cout << LOG;}
#define LOG_WARNING(INFO, LOG) if(LOGGING_IS_ENABLED_(INFO, Base::LogLevel::Warning, _is_unlikely)) {INFO->cerr << LOG;}


} // namespace Base

