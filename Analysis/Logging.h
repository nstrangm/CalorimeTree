#ifndef GAMMACONV_Logging
#define GAMMACONV_Logging

#include <iostream>

using std::cerr; //  Preferably use cerr since cout is not always printed exactly where called
using std::cout; //  Now the std:: in std::cout can be omitted
using std::endl;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++ Logging +++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #define ENTER                                                                            \
    {                                                                                      \
      std::fprintf(stdout, "\033[37m[ ENTER ] %s -> %d\033[0m\n", __FUNCTION__, __LINE__); \
    }
  #define EXIT                                                                             \
    {                                                                                      \
      std::fprintf(stdout, "\033[37m[ EXIT  ] %s -> %d\033[0m\n", __FUNCTION__, __LINE__); \
    }
  #define INFO(s)                                               \
    {                                                           \
      std::fprintf(stdout, "\033[92m[ INFO  ] %s\033[0m\n", s); \
    }
  #define LOG(s)                                                \
    {                                                           \
      std::fprintf(stdout, "\033[96m[  LOG  ] %s\033[0m\n", s); \
    }
  #define WARN(s)                                               \
    {                                                           \
      std::fprintf(stderr, "\033[93m[ WARN  ] %s\033[0m\n", s); \
    }
  #define ERROR(s)                                                                    \
    {                                                                                 \
      std::fprintf(stderr, "\033[91m[ ERROR ] in line %d: %s\033[0m\n", __LINE__, s); \
    }
  #define ERRORETURN(s)                                                               \
    {                                                                                 \
      std::fprintf(stderr, "\033[91m[ ERROR ] in line %d: %s\033[0m\n", __LINE__, s); \
      return;                                                                         \
    }
  #define DEBUG                                                                                            \
    {                                                                                                      \
      std::fprintf(stderr, "\033[95m[ DEBUG ] %s -> %s -> %d\033[0m\n", __FILE__, __FUNCTION__, __LINE__); \
    }
  #define FATAL(s)                                                                                                    \
    {                                                                                                                 \
      std::fprintf(stderr, "\033[91m[ FATAL ] in %s -> %s -> %d:\n%s\033[0m\n", __FILE__, __FUNCTION__, __LINE__, s); \
      exit(0);                                                                                                        \
    }

#endif
