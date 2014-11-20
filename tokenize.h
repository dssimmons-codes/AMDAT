/*Declaration of function to take a line and convert it into substrings broken up by spaces and tabs*/
/*Written by David S. Simmons*/



#ifndef TOKENIZE
#define TOKENIZE
#include <string>
#include <vector>

#define ARGMAX 100

namespace std{

class Tokenize
{
    vector <string> tokens;
    bool not_empty;

public:
    Tokenize(){};			//default constructor does nothing
    Tokenize(string line);		//constructor to tokenize initial line
    Tokenize(const Tokenize & copy);	//copy constructor
    Tokenize operator=(const Tokenize & copy);	//equality operator
    //void operator()(string line);	//tokenize new line
    vector <string> operator()(string line);
    string operator()(int tokenindex){return tokens[tokenindex];};	//return token corresponding to tokenindex
    int operator()(string* tokenarray,int maxsize);	//returns number of tokens and sets tokenarray* to be an array of tokens up to size maxsize
    int operator()(string line, string* tokenarray, int maxsize=ARGMAX);	//tokenizes new line, returns number of tokens, and sets tokenarray to contain tokens
	
    
    //int operator(string line,string* tokenlist);
    int count(){return tokens.size();};		//returns number of tokens

    bool in_string_array(string target);	//check if a token is in the line
    int find_in_string_array(string target);	//return the position of a token if in the line
};

}
  
#endif


