// Copyright (C) 2010-2025 David S. Simmons

#ifndef ERROR
#define ERROR

#include <string>
#include <stdio.h>


namespace std{
class Error{

void throw_severe(string, int, bool);
void throw_moderate(string, int, bool);
void throw_minor(string, int, bool);
void print_msg_line(string, bool);

public:
	Error (string, int, bool); //general error contructor
	Error (string, int); //external-only error constructor
	Error (); //null constructor
	Error (const Error &); //copy
	Error operator = (const Error &);	//assignment


};
}

#endif
