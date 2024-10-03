#ifndef GAMMACONV_Logging
#define GAMMACONV_Logging

#include <iostream>

using std::cerr; //  Preferably use cerr since cout is not always printed exactly where called
using std::cout; //  Now the std:: in std::cout can be omitted
using std::endl;

int debugLevel = 2; // 0: No output except for warnings, errors and fatals, 1: +Log statements, 2: +Info+ENTER/EXIT statementsm 3: All output including DEBUG

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++ Logging +++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define ENTER                                                                              \
  {                                                                                        \
    if (debugLevel > 1)                                                                    \
      std::fprintf(stdout, "\033[37m[ ENTER ] %s -> %d\033[0m\n", __FUNCTION__, __LINE__); \
  }
#define EXIT                                                                               \
  {                                                                                        \
    if (debugLevel > 1)                                                                    \
      std::fprintf(stdout, "\033[37m[ EXIT  ] %s -> %d\033[0m\n", __FUNCTION__, __LINE__); \
  }
#define INFO(s)                                                 \
  {                                                             \
    if (debugLevel > 1)                                         \
      std::fprintf(stdout, "\033[92m[ INFO  ] %s\033[0m\n", s); \
  }
#define LOG(s)                                                  \
  {                                                             \
    if (debugLevel > 0)                                         \
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
#define DEBUG                                                                                              \
  {                                                                                                        \
    if (debugLevel > 2)                                                                                    \
      std::fprintf(stderr, "\033[95m[ DEBUG ] %s -> %s -> %d\033[0m\n", __FILE__, __FUNCTION__, __LINE__); \
  }
#define FATAL(s)                                                                                                    \
  {                                                                                                                 \
    std::fprintf(stderr, "\033[91m[ FATAL ] in %s -> %s -> %d:\n%s\033[0m\n", __FILE__, __FUNCTION__, __LINE__, s); \
    exit(0);                                                                                                        \
  }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++ Extra functions ++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TString ReturnTimehms(int s)
{
  int h = s / 3600;
  s = s % 3600;
  int m = s / 60;
  s = s % 60;
  return Form("%d hours, %d minutes and %d seconds", h, m, s);
}

void PrintProgressBar(int i, int N, int Steps = 100)
{

  int barWidth = 100;

  static int Progress = 0;

  static clock_t t = clock();

  if (((double)i) / ((double)N) * Steps > Progress)
  {

    cout.flush();
    cout << " [" << Progress * 100 / Steps << "%]"
         << "[";
    int pos = barWidth * (((double)Progress) / Steps);
    for (int i = 0; i < barWidth; ++i)
    {
      if (i < pos)
        cout << "|";
      else
        cout << " ";
    }
    cout << "] - " << ReturnTimehms(((clock() - t) / CLOCKS_PER_SEC) * ((double)(N - i) / (double)i)) << " left.  "
         << "\r";

    Progress++;
  }
  if (i == N - 1)
  {
    cout.flush();
    cout << "[" << Progress << "%]"
         << "[";
    int pos = barWidth * (Progress / 100.);
    for (int i = 0; i < barWidth; ++i)
    {
      cout << "|";
    }
    cout << "] - "
         << "Processing finished in " << ReturnTimehms((clock() - t) / CLOCKS_PER_SEC) << endl;
  }
}

void PrintProgressNumber(int i, int N, int Steps = 100)
{
  static int Progress = 0;
  if (((double)i) / ((double)N) * Steps > Progress)
  {
    cout << "[" << Form("%.1f", Progress * 100. / Steps) << "%]" << endl;
    Progress++;
  }
}

void PrintProgress(int i, int N, int Steps = 100, bool isDebug = false)
{
  if (isDebug)
    PrintProgressBar(i, N);
  else
    PrintProgressNumber(i, N, 1000);
}


#endif
