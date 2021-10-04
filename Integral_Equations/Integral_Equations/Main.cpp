#ifndef ATTACH_H
#include "Attach.h"
#endif

// The aim of this project is to implement the methods used to numerically solve integral equations
// Code will focus on Fredholm and Volterra equations of the 2nd type
// R. Sheehan 4 - 10 - 2021

int main()
{

	testing::gaussj_test(); 

	testing::ludcmp_test(); 

	std::cout << "Press enter to close console\n";
	std::cin.get();

	return 0; 
}