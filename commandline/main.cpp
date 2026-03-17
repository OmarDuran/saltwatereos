
#include "SWEOSbash.h"
#include <iostream>
int main( int argc, char** argv )
{
  std::cout << "\033[36m"
            << "  Fork & TBB parallelization:\n"
            << "    Omar Duran, University of Bergen (UiB), 2026\n"
            << "  Original source:\n"
            << "    github.com/zguoch/salern\n"
            << "\033[0m" << std::endl;
  SWEOSbash::bash_run(argc, argv);
}
