#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char const* argv[])
{
   ofstream file;
   file.open("dualformetis.dat.epart.1");
   for (int i = 0; i < 5768; i++)
   {
      file << "0 \n";
   }

   return 0;
}
