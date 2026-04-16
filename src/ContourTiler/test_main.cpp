#include <ContourTiler/main.h>

// #include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/loglevel.h>
#include <log4cplus/configurator.h>
#include <iomanip>

#include <CppUnitLite2.h>
#include <TestResultStdErr.h>

#include <ContourTiler/test_common.h>

string data_dir;
string out_dir;
Number_type cl_delta;

using namespace log4cplus;

int main(int argc, char** argv)
{
  PropertyConfigurator::doConfigure("log4cplus.properties");

  log4cplus::Logger logger = log4cplus::Logger::getInstance("");
  //LOG4CPLUS_INFO(logger, "Configured logger");

  data_dir = "test_data";
  out_dir = ".";
  cl_delta = 0.01;
  string test;
  for (int i = 1; i < argc; ++i)
  {
    if ((string)argv[i] == "--data")
    {
      data_dir = argv[i+1];
      ++i;
    }
    else if ((string)argv[i] == "--out")
    {
      out_dir = argv[i+1];
      ++i;
    }
    else if ((string)argv[i] == "--delta")
    {
      cl_delta = boost::lexical_cast<Number_type>(argv[i+1]);
      ++i;
    }
    else if ((string)argv[i] == "--test")
    {
      test = argv[i+1];
      test = test + "Test";
      ++i;
    }
  }

  TestResultStdErr result;
  if (test != "")
    TestRegistry::Instance().Run(result, test.c_str());
  else
  {
    TestRegistry::Instance().Run(result);
  }
  return (result.FailureCount());
}

