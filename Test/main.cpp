#include <gtest/gtest.h>

std::string path;
char delimiter = '\0';

int main(int argc, char* argv[])
{
  std::string executable = *argv;
  path = executable.substr(0, executable.find_last_of("/\\") + 1);
  delimiter = path.back();

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}