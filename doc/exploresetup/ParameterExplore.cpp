/*
Cameron Kates
7/15/2013

ParameterExplore

ParameterExplore is a version of inputRadiusChanger modified and renamed for public distribution.

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits>

using namespace std;

void ReadMe()
{
	cout << "

int main()
{
	double startR; 
	double endR;
	double numR;
	double startRm; 
	double startE; 
	double endE;
	double fixedE;
	double lowRange;
	double highRange;
	double go;

	cout << "ParameterExplore 1.6 - 7/19/2013" << endl;
	cout << "To view ReadMe file enter -1 for minimum radius." << endl;
	
	cout << "Please enter the minimium radius: " << endl;
	cin >> startR;
	cout << "Please enter the maximum radius: " << endl;
	cin >> endR;
	
	//Calculates the total number of radii that will be explored.
	numR = ((endR-startR)+0.1)*10; 

	//Determines if the user wants to vary matching radius and acts accordingly.
	char rmChoice;
	cout << "Do you want to vary matching radius? y/n" << endl;
	cin >> rmChoice;
	int rmKey=1;
	while(rmChoice != 'n' && rmChoice != 'N' && rmChoice != 'y' && rmChoice != 'Y')
	{
		cout << "Input not recognized. Please try again." << endl;
		cin >>rmChoice;
	}
	if(rmChoice == 'y' || rmChoice == 'Y')
	{
		cout << "Please enter the minimium matching radius: " << endl;
		cin >> startRm;
		while(startRm > startR)
		{
			cout << "The matching radius must be less than or equal to the maximum radius. Please enter a new value." << 				endl;
			cin >> startRm;
		}
	}
	else if(rmChoice == 'n' || rmChoice == 'N')
	{
		rmKey = -1;
	}
	
	//Determines if the user wants to vary energy and acts accordingly.
	char eChoice;
	cout << "Do you want to vary energy? y/n" << endl;
	cin >> eChoice;
	int eKey=1;
	while(eChoice != 'n' && eChoice != 'N' && eChoice != 'y' && eChoice != 'Y')
	{
		cout << "Input not recognized. Please try again." << endl;
		cin >>eChoice;
	}
	if(eChoice == 'y' || eChoice == 'Y')
	{
		cout << "Please enter the minimium energy: " << endl;
		cin >> startE;
		cout << "Please enter the maximum energy: " << endl;
		cin >> endE;
	}
	else if(eChoice == 'n' || eChoice == 'N')
	{
		cout << "Please enter the fixed energy: " << endl;
		cin >> fixedE;		
		eKey = -1;
	}
	
	
///////////////////////////////////////////////////////////////////
	
	int caseCount=0;
	char userCheck;
	char warningCheck;

	//Case count. 0.001 added to compensate for possible roundoff error. This is essentially a null version of the work loops.
	for(double R = startR; R<=endR+0.001; R=R+0.1)
		{
			if(rmKey==-1)
			{
				startRm = R;
			}
			for(double Rm = startRm; Rm<=R+0.001; Rm=Rm+0.1)
			{
				double Rsmall = R-0.2;
				if(eKey==-1)
				{
					endE=fixedE;
					startE=endE;
				}
				for(double E = startE; E<=endE+0.001; E=E+0.1)
				{
					caseCount++;
				}
			}
		}

	cout << "Over what range would you like to calculate the log derivatives?" << endl;
	cin >> lowRange;
	cin >> highRange;
		
	cout << "This will generate " << caseCount << " cases. Do you want to proceed? y/n" << endl;	
	cin >> userCheck;
	if(userCheck == 'y' || userCheck == 'Y')
	{
		if(caseCount >= 4000)
		{
			cout << "WARNING: THIS WILL CREATE A POTENTIALLY DANGEROUS NUMBER OF CASES! Are you sure you wish to continue? y/n" << endl;
			cin >> warningCheck;
			userCheck = warningCheck;
		}
	}

	//Begin changing input file.	
	if(userCheck == 'y' || userCheck == 'Y')
	{
		//Creates the InputRecords file which keeps track of the users inputs for future reference.
		cout << "Creating records file..." << endl;
		time_t curr=time(0); //Finds current date and time.
		ofstream recordsfile;
		recordsfile.open ("InputRecords");
		recordsfile << "Created on: " << ctime(&curr);
		recordsfile << "Using ParameterExploreV1.6" << endl;			//MENTION OF VERSION
		recordsfile << " " << endl;
		recordsfile << "Rc: " << startR << " " << endR << endl;
		if(rmKey == 1)
		{
			recordsfile << "Rmatch: " << startRm << " " << endl;
		}
		if(rmKey == -1)
		{
			recordsfile << "Rmatch: N/A " << endl;
		}
		if(eKey == 1)
		{
			recordsfile << "Energy: " << startE << " " << endE << endl;
		}
		if(eKey == -1)
		{
			recordsfile << "Energy: N/A " << endl;
		}
		recordsfile.close();

		//Changes the 'in' file.
		cout << "Creating file..." << endl;
		ofstream myfile;
		myfile.open ("parameterExplore");					
		myfile << "echo '" <<caseCount<< " " << numR << " LOGDERIVERANGE " << lowRange << " " << highRange << "' >> in" << endl;
	
		for(double R = startR; R<=endR+0.001; R=R+0.1) //0.001 added to compensate for possible roundoff 			error.
		{
			if(rmKey==-1)
			{
				startRm = R;
			}
			for(double Rm = startRm; Rm<=R+0.001; Rm=Rm+0.1)
			{
				double Rsmall = R-0.2;
				if(eKey==-1)
				{
					endE=fixedE;
					startE=endE;
				}
				for(double E = startE; E<=endE+0.001; E=E+0.1)
				{
					myfile << "sed 's/En/"<<E<<"/' PseudoTemplate|sed 's/Rin/"<<R<<"/g' |sed 's/Rsmall/" 
					<<Rsmall<<"/' |sed 's/Rmatch/"<<Rm<<"/' >>in" << endl;	
				}
			}
		}
			
		myfile << "echo 'END' >> in" << endl;
		myfile.close();
		cout << "File created." << endl;
		cout << "Generating executable..." << endl;
		system("chmod 755 parameterExplore");	
		cout << "Executable file parameterExplore created." << endl;	
		cout << "Changing 'in' file values..." << endl;
		system("parameterExplore");
		cout << "Complete." << endl;
		return 0;
	}
	else if(userCheck == 'n' || userCheck == 'N')
	{
		cout << "Program terminated. Please try again." << endl;
		return 0;
	}
	else
	{
		cout << "Input not recognized. Program terminated. Please try again." << endl;
		return 0;
	}	
}
