// Copyright (C) 2010-2025 David S. Simmons

/*Function to take a line and convert it into substrings broken up by spaces and tabs*/
/*Written by David S. Simmons*/

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;

#include "tokenize.h"
#define ARGMAX 10000

Tokenize::Tokenize(string line)
{
  flagcheck=false;
  flagmarker="";
  (*this)(line);
};

Tokenize::Tokenize(const Tokenize & copy)
{
  tokens = copy.tokens;
  flags = copy.flags;
  flagmarker=copy.flagmarker;
  flagcheck=copy.flagcheck;
}

Tokenize Tokenize::operator=(const Tokenize & copy)
{
  if(this!=&copy)
  {
    tokens = copy.tokens;
    flags = copy.flags;
    flagmarker=copy.flagmarker;
    flagcheck=copy.flagcheck;
    return *this;
  }
}

vector <string> Tokenize::operator()(string line)
{
  tokens.clear();
  flags.clear();

  
  size_t token_start = 0;
  size_t next_space;
  size_t next_tab;
  size_t token_end;
  size_t token_size;
  
  tokens.clear();
  
  if(line=="")
  {
    not_empty = 0;
  }
  else
  {
    not_empty = 1;
  
  while(1)	//loop until break command, at end of line
  {
    next_space = line.find(' ',token_start);	//find next space
    next_tab = line.find("\t",token_start);	//find next tab
    
    if(next_space==line.npos && next_tab==line.npos)	//check if end of line has been reached
    {
      if(token_start!=line.length())		//check if there is any token left to copy
      {
        parse_token(line,token_start);
	//tokens.push_back(line.substr(token_start));	//if so, pull token from line
      }
      break;
    }
    
    else if(next_space == line.npos)		//if there are no more spaces, use tab as next delimiter
    {token_end = next_tab;}
    
    else if(next_tab == line.npos)		//if there are no more tabs, use space as next delimiter
    {token_end = next_space;}
    
    else
    {
      if(int(next_space)<int(next_tab)){token_end=next_space;}	//if space comes first, use it for the next tokenbreak
      else{token_end = next_tab;}	//if tab comes first, use it for the next tokenbreak
    }
    
    token_size = token_end-token_start;	//determine token size
    if(int(token_size) != 0)		//check if token size is zero (as in the case of consecutive tabs/spaces).  If not, use these delimiters to get the next token
    {
      parse_token(line,token_start,token_end);
      //tokens.push_back(line.substr(token_start,token_size));	//pull token from line
    } 
    
    token_start = token_end + 1;				//move up read start
  }
  }
  return tokens;
}




void Tokenize::parse_token(string line, size_t token_start)
{
  string firstchar = line.substr(token_start,1);
  string flagid, flag;
  
  if(flagcheck&&firstchar==flagmarker)
  {
    flagid = line.substr(token_start+1,1);
    flag = line.substr(token_start+2);
    flags.insert(flagid,flag);
  }
  else
  {
    tokens.push_back(line.substr(token_start));
  }
}



void Tokenize::parse_token(string line, size_t token_start, size_t token_end)
{
  string firstchar = line.substr(token_start,1);
  string flagid, flag;
  
  if(flagcheck&&firstchar==flagmarker)
  {
    flagid = line.substr(token_start+1,1);
    flag = line.substr(token_start+2,token_end-token_start-2);
    flags.insert(flagid,flag);
  }
  else
  {
    tokens.push_back(line.substr(token_start,token_end-token_start));
  }
}




bool Tokenize::in_string_array(string target)
{
  for(int stringii=0; stringii<tokens.size(); stringii++)
  {
    if(target==tokens[stringii])
    {return 1;}
  }
  return 0;
}


  

int Tokenize::find_in_string_array(string target)
{
  for(int stringii=0; stringii<tokens.size(); stringii++)
  {
    if(target==tokens[stringii])
    {return stringii;}
  }
  return -1;
}


int Tokenize::operator()(string*tokenarray, int maxsize)
{
  if(tokens.size()>maxsize)
  {
    cout<<"\nError: attempting to convert token vector into unacceptably small token array\n";
    exit(0);
  }
  
  for(int tokenii=0;tokenii<tokens.size();tokenii++)
  {
    tokenarray[tokenii]=tokens[tokenii];
  }
  
  return tokens.size();
}


string Tokenize::operator[](string flagid)
{
  int location;
  location = flags.find(flagid);
  if(location==-1)
  {
    return "";
  }
  else
  {
    return flags.at(flagid);
  }
}


int Tokenize::operator()(string line, string*tokenarray, int maxsize)
{
  (*this)(line);
  
  if(tokens.size()>maxsize)
  {
    cout<<"\nError: attempting to convert token vector into unacceptably small token array\n";
    exit(0);
  }
  
  for(int tokenii=0;tokenii<tokens.size();tokenii++)
  {
    tokenarray[tokenii]=tokens[tokenii];
  }
  
  return tokens.size();
}


#ifdef NEVER

int Tokenize::operator()(string line, string* tokenlist)
{
  (*this)(line);
  
  



int tokenize(string line, string * tokens)
{
  int n_tokens = 0;
  size_t token_start = 0;
  size_t next_space;
  size_t next_tab;
  size_t token_end;
  size_t token_size;
  
  while(1)	//loop until break command, at end of line
  {
    next_space = line.find(' ',token_start);	//find next space
    next_tab = line.find("\t",token_start);	//find next tab
    
    if(next_space==line.npos && next_tab==line.npos)	//check if end of line has been reached
    {
      if(token_start!=line.length())	//check if there is any token left to copy
      {
        tokens[n_tokens] = line.substr(token_start);	//if so, pull token from line
        n_tokens++;	//increment count of tokens
      }
      return n_tokens;	//since end of string has been reached, terminate method
    }
    
    else if(next_space == line.npos)	//if there are no more spaces, use tab as next delimiter
    {token_end = next_tab;}
    
    else if(next_tab == line.npos)	//if there are no more tabls, use space as next delimiter
    {token_end = next_space;}
    
    else
    {
      if(int(next_space)<int(next_tab)){token_end=next_space;}	//if space comes first, use it for the next tokenbreak
      else{token_end = next_tab;}	//if tab comes first, use it for the next tokenbreak
    }
    
    token_size = token_end-token_start;	//determine token size
    if(int(token_size) != 0)	//check if token size is zero (as in the case of consecutive tabs/spaces).  If not, use these delimiters to get the next token
    {
      tokens[n_tokens] = line.substr(token_start,token_size);	//pull token from line
      n_tokens++;	//increment count of tokens
    } 
    
    token_start = token_end + 1;	//move up read start
  }
}

#endif