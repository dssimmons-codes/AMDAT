/*Declaration of function to take a line and convert it into substrings broken up by spaces and tabs*/
/*Written by David S. Simmons*/



#ifndef TOKENIZE
#define TOKENIZE
#include <string>
#include <vector>
#include "vector_map.h"

#define ARGMAX 100

namespace std{

class Tokenize
{
    vector <string> tokens;
    Vector_Map <string,string> flags;
    bool not_empty;
    void parse_token(string line, size_t token_start);
    void parse_token(string line, size_t token_start, size_t token_end);
    string flagmarker;
    bool flagcheck;

public:
    Tokenize(){flagcheck=false;flagmarker="";};			//default constructor just sets flag checking to off
    Tokenize(string line);		//constructor to tokenize initial line
    Tokenize(const Tokenize & copy);	//copy constructor
    Tokenize operator=(const Tokenize & copy);	//equality operator
    //void operator()(string line);	//tokenize new line
    vector <string> operator()(string line);
    
    void setflagmarker(string marker){flagmarker=marker;flagcheck=true;};	//set the flagmarker and set to check for flags
    void setflagmarker(){flagmarker="";flagcheck=false;};	//clear the flagmarker and set to not check for flags
    
    string operator()(int tokenindex){return tokens[tokenindex];};	//return token corresponding to tokenindex
    string operator[](string flagid);
    int operator()(string* tokenarray,int maxsize);	//returns number of tokens and sets tokenarray* to be an array of tokens up to size maxsize
    int operator()(string line, string* tokenarray, int maxsize=ARGMAX);	//tokenizes new line, returns number of tokens, and sets tokenarray to contain tokens
	
    
    //int operator(string line,string* tokenlist);
    int count(){return tokens.size();};		//returns number of tokens
    bool isflagged(string flagid){return flags.count(flagid);};		//returns a 1 if flag given by flag id is in line; 0 otherwise.
    
    
    bool in_string_array(string target);	//check if a token is in the line
    int find_in_string_array(string target);	//return the position of a token if in the line
};

}
  
#endif


