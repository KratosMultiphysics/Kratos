#include <stdio.h>
#include <iostream>
using namespace std;

void showProgress(int step, int maxStep)
{
	//system("cls");
	std::cout << step <<" of "<<maxStep;
}

void printMessage(char* message)
{	 
  std::cout << message << "\n";	
};