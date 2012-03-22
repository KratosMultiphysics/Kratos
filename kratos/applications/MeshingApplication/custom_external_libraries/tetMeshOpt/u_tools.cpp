#include <stdio.h>
#include <iostream>
using namespace std;

void showProgress(int step, int maxStep)
{
	system("cls");
	cout << step <<" of "<<maxStep;
}

void printMessage(char* message)
{	 
  cout << message << "\n";	
};