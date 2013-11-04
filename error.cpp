/********************************/
/* Error handling class			*/
/*	Written by Michael Marvin	*/
/********************************/


/************************************************************************************************************
    Error code levels:
<1 - SEVERE - AMDAT can not continue with this error and must exit
<10000 - MODERATE - AMDAT can not complete the current command with this error and will attempt to skip it
>=10000 - MINOR - AMDAT can still continue, but things may be incorrect

    Error code meanings, if you take a code please put it here!
0  - General SEVERE error
-1 - Input file not found
-2 - System file not found
-3 - End of if/loop not found
-4 - Inconsistent number of atoms
-5
-6 - Incorrect number of arguments given, kill execution version


1 - General MODERATE error
2 - Output file not found
3 - Constant not found
4 - Unbalanced parenthesis
5 - No atom set provided/found
6 - Incorrect number of arguments given, skip command version


10000 - General MINOR error
10001 - Number of processors invalid
10002 - Attempt to seek past the end of the input file

************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>


#include "error.h"
#include "control.h"

using namespace std;

/* Contructor method. Requires an error message and an error code */
Error::Error (string msg, int code, bool internal)
{
    if (code<1)
        throw_severe(msg, code, internal);
    else if (code<10000)
        throw_moderate(msg, code, internal);
    else
        throw_minor(msg, code, internal);

}

/* Contructor for external-only (input script) errors */
Error::Error (string msg, int code)
{
    Error (msg, code, false);
}

Error::Error()
{
    Error ("Improperly formatted error", 10000);
}

/*Copy constructor*/
Error::Error(const Error & copy)
{
}
  

/*Assignment operator*/
Error Error::operator =(const Error & copy)
{
  return *this;
}




void Error::throw_severe(string msg, int code, bool internal)
{   
    cerr << endl;
	if (!internal)
	    cerr << "Error! An error occured and AMDAT is not able to continue." << endl;
	else
	    cerr << "Error! An internal error occured and AMDAT is not able to continue." << endl;
    cerr << "Error message: ";
    print_msg_line(msg, internal);
    cerr << "Error Code: " << code << endl;
    exit(code);
}

void Error::throw_moderate(string msg, int code, bool internal)
{
    cerr << endl;
	if (!internal)
	    cerr << "Error! An error occured and AMDAT must skip execution of the current command." << endl;
	else
	    cerr << "Error! An internal error occured and AMDAT must skip execution of the current command." << endl;
    cerr << "Error message: ";
    print_msg_line(msg, internal);
    cerr << "Error Code: " << code << endl;
	int line_num=Control::get_line_number();
	if (line_num<Control::get_input_file_length())
	    Control::line_seek(line_num+1);
	else
	    Control::line_seek(Control::get_input_file_length());
}

void Error::throw_minor(string msg, int code, bool internal)
{
    cerr << endl;
	if (!internal)
	    cerr << "Error! An error occured but AMDAT will try to continue." << endl;
	else
	    cerr << "Error! An internal error occured but AMDAT will try to continue." << endl;
    cerr << "Error message: ";
    print_msg_line(msg, internal);
    cerr << "Error Code: " << code << endl;
}

void Error::print_msg_line(string msg, bool internal)
{
    if (!internal)
        cerr << "\"" << msg << "\" at command \"" << Control::get_raw_line(Control::get_line_number()-1) << "\"" << endl;
    else
        cerr << "\"" << msg << "\"" << endl;
}
