/*Declaration of function to take a line and convert it into substrings broken up by spaces and tabs*/
/*Written by David S. Simmons*/



#ifndef TOKENIZE
#define TOKENIZE
#include <string>

using namespace std;

int tokenize(string, string *);
bool in_string_array(string * tokens, int array_size, string target);
int find_in_string_array(string * tokens, int array_size, string target);

#endif