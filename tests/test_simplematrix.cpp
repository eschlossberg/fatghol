#include <iostream>

#include "simple.hpp"

main(int argc, char **argv)
{
  SimpleMatrix M(100, 100);

  int i = 0;
  int j = 0;
  for (int x = 1; x < 10; x++)
    {
      M.entry(i, j, x);
      i += 13; i %= 100;
      j += 29; j %= 100;
    }

  std::cout << M.rank() << std::endl;
}
