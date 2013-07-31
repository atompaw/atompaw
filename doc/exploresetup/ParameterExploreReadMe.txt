/////////////////////////////////
/////////////////////////////////
Input Radius Changer, or ircv, is an executable written in C++ for the Atomic Structure Database for First Principles Modeling research project with the Natalie Holzwarth research group. It take an input of an atom specific 'in' file as well as a variable template 'u2' file. It will then give the user the option to vary the parameters radius Rc, matching radius Rm, and energy E. These parameters will then be used in combination with the u2 file to append the new EXPLORE cases to the end of the in file. The in file can then be submitted to the atomic modeling software atompaw.

Record Start:

Version Update - V5
-This version expands on V4 by adding the capability to ask the user if they would like to vary matching radii as well as energies or not. Additionally, the method to do this was changed from a multi-loop method in V3 to an equal-value loop termination method, where the loop is only allowed to run once if an input is not varied. This simplifies and streamlines the code.

6/11/2013 - Bug Fix
User JDrewery discovered a bug in the caseCount function. The program was reporting too few cases on the input file. Bug determined to be caused by a roundoff error. Recopied loop from actual file writing portion of the code and replaced file write line with caseCount++. Bug resolved.

6/18/2013
V5 Update - Feature addition, file organization change.
-Directory organization changed so the user never sees the version changes. 
-Adds the capability for the program to create an output file which records the parameters entered by the user and adds a time and   date stamp.
-ircv now outputs the total number of Rc radii it varies across in the same line as the total number of cases.  

6/19/2013
V5.1 - Minor Feature Addition
-Return more obvious warning to the user if they are going to create too many files.

6/20/2013
V5.1 - Bug Note
-If the user inputs a letter for the minimum radius the program gets stuck in an infinite loop. Issue is easily enough avoided that it doesn't warrant a fix. 

7/15/2013
Version Update - V6
-NAWH has made a change to atompaw which allows the user to specify the range of logderivative errors over which to calulate. This version adds the ability for the user to input these values. 
-Added a message to the user about accessing the ircvReadMe file. ircvReadMe file was previously known as ircvNotes. 

7/19/2013
Version Update - V1.7
-There is a desire to release IRCV online for the public. The name of IRCV has been changed to ParameterExplore and the u2 file has been renamed PseudoTemplate. Additionally, the version numbering standard has been changed. Now, the format is as follows:
V(MajorFunctionalityUpdate).(FunctionAddition/FeatureAddition.(BugFix/MinorChange within the current major and function numbers).


##2013##
##Cameron Kates##
##WFU##
/////////////////////////////////
/////////////////////////////////
