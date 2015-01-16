/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Control class methods*/
/*Written by David S. Simmons*/

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <unistd.h>
#include <iterator>
#include <stdexcept>

#include "mean_square_displacement.h"
#include "van_hove_self.h"
#include "van_hove_distinct.h"
#include "control.h"
#include "wave_density.h"
#include "intermediate_scattering_function.h"
#include "incoherent_scattering_function.h"
#include "debyewaller_dist.h"
#include "stiffness_dist.h"
#include "displacement_distribution.h"
#include "non_gaussian_parameter.h"
#include "gaussian_comparison.h"
#include "fast_particles.h"
#include "radial_debye_waller.h"
#include "mean_square_displacement_2d.h"
#include "velocity_autocorrelation.h"
#include "strings.h"
#include "static_trajectory_list.h"
#include "rgtensor_stats.h"
#include "displacement_map.h"
#include "composition.h"
#include "n_fold_order_parameter.h"
#include "structure_factor.h"
#include "clustered_list.h"
#include "vector_autocorrelation.h"
#include "trajectory_list_decay.h"
#include "mean_displacement.h"
#include "edge_detector_timedependent.h"
#include "multibody_set.h"
#include "gyration_radius.h"
#include "multibody_list.h"
#include "multibody_analysis.h"
#include "version.h"
#include "tokenize.h"
#include "error.h"

using namespace std;



//TODO: Convert all argcheck() calls to bool_argcheck() for error handling
//      Convert all errors to use Error class
//      Add <<endl; to nearly all cout references

/*--------------------------------------------------------------------------------*/

string * Control::constants = new string[LISTSIZE];
string * Control::constant_names = new string[LISTSIZE];
int Control::n_constants = 0;
int Control::current_line = 0;
vector<string> Control::inputFileVector;

/*Constructor method - gets command input from file and calls appropriate method*/
Control::Control(char * filename_input, string * consts, string * consts_names, int n_consts,string tempfile)
{
  constants=consts;
  constant_names=consts_names;
  n_constants=n_consts;
  MAXTHREADS=omp_get_max_threads();
  int numLines = read_input_file(filename_input);

  if (!execute_commands(0, numLines))
    cout << "\nProgram terminated prior to completion! An error was encountered." << endl;
  else
    cout << "\nProgram completed successfully." << endl;
}



/*---------------------------------------------------------------------------------*/
int Control::read_input_file(char * filename_input)
{
  /** Reads the input file into a vector of strings called inputFileVector and filters out comments
  * @author Michael Marvin
  * @date 7/22/2013
  **/
  int argii;
  string line;
  string command;

  initialize_lists();

  input.open(filename_input);
  //if (!input.is_open()){cout << "Error opening file " << filename_input << ".\n"; exit(1);}
  if (!input.is_open()){Error (string("Error opening file ").append(filename_input), -1, true);}
  int numLines=0;
  while(!input.eof())
  {
    line = "";
    command = "";
    for(argii=0;argii<ARGMAX;argii++)
    {
      args[argii]="";
    }

    getline(input,line);

    int commentStart = line.find("#",0);

    if (commentStart != string::npos)
    {
        line=line.substr(0, commentStart);
        if (line == "")
        { continue; }
    }

    n_args = tokenize(line, args);
    if(n_args==0)
    { continue; }
    else
    { command = args[0]; }

    if (command == "skip")
    {
        skip();
        continue;
    }

    inputFileVector.push_back(line);
    numLines++;

  }
  input.close();

//  cout << "\n";
/*  for(int i=0; i<numLines; i++) // Display the contents of the vector, aka, display the file minus any comments. TODO: Make this an option. MDM 10/21/13
  {
    cout << inputFileVector[i] << endl;
  }*/
  cout << "Running with " << omp_get_max_threads() << " threads.\n" << endl;
  return numLines;
}

int Control::execute_commands(int iIndex, int fIndex)
{
  /** Executes the commands in the input file between the two supplied indices
  * Error codes: 0: error - stop script , 1: success - continue, 2: break loop
  **/
  string command;
  string line;
  line_seek(iIndex);
  while (current_line < fIndex)
  {
    command = "";
    line = "";
    for(int argii=0;argii<ARGMAX;argii++)
    {
      args[argii]="";
    }

    line=read_line();

    n_args = tokenize(line, args);
    if(n_args==0)
    {
      command = "none";
      continue;
    }
    else
    {
      command = args[0];
      cout << line << endl; // TODO: Need to make this an option. Either compiler or command line. Should it be on by default? MDM 10/21/13
    }
    if (command.find("#",0) != string::npos)
    {
      command = "none";
      continue;
    }
    else if(command == "system" || command =="system_nv")
    {
      analyte=new System(inputFileVector);
      //System system(inputFileVector);							//instantiate system object
      //analyte = &system;
    }
    else if(command == "system_np")
    {
      analyte = new System(inputFileVector,1);
      //System system(inputFileVector,1);	//instantiate system object in non-constant-volume ensemble
      //analyte = &system;
    }
    else if (command == "create_list")
    {
      create_list();
    }
    else if (command == "create_multibodies")
    {
      create_multibodies();
    }
    else if (command == "trajectories_from_multibodies")
    {
      trajectories_from_multibodies();
    }
    else if (command == "combine_trajectories")
    {
        combine_trajectories();
    }
    else if (command == "msd")
    {
      msd();
    }
    else if (command == "msd_2d")
    {
      msd_2d();
    }
    else if (command == "mean_displacement")
    {
      mean_displacement();
    }
    else if (command == "vac_function")
    {
      vacf();
    }
    #ifndef TACC
    else if (command == "vac_fourier")
    {
      vacf_fourier();
    }
    #endif
    else if (command == "vhs")
    {
      calc_vhs();
    }
    else if (command == "vhd")
    {
      calc_vhd();
    }
    else if (command == "vht")
    {
      calc_vht();
    }
    else if (command == "limit")
    {
      limit();
    }
    else if (command == "isfs")
    {isfs();}
    #ifndef TACC
    else if (command == "isfd")
    {vhd_fourier();}
    #endif
    else if (command == "structure_factor")
    {structure_factor();}
    else if (command == "isf")
    {isf();}
    else if (command == "isf_list")
    {isf_list();}
    else if(command == "u2dist")
    {u2dist();}
    else if(command == "stiffness_dist")
    {stiffness_dist();}
    else if(command == "displacement_dist")
    {displacement_dist();}
    else if(command == "ngp")
    {ngp();}
    else if(command == "compare_gaussian")
    {compare_gaussian();}
    else if(command == "find_fast")
    {find_fast();}
    else if(command == "radial_debye_waller")
    {radial_debye_waller();}
    else if(command == "strings")
    {strings();}
    else if(command == "rgtensor_stats")
    {rgtensor_stats();}
    else if (command == "displacement_map")
    {displacement_map();}
    else if (command == "write_starr")
    {analyte->write_starr();}
    else if (command == "write_single_particle")
    {write_single_particle();}
    else if (command == "write_list_trajectory")
    {write_list_trajectory();}
    else if (command == "write_list_trajectory_full")
    {write_list_trajectory_full();}
    else if (command == "create_bin_list")
    {create_bin_list();}
    else if (command == "remove_bin_list")
    {remove_bin_list();}
    else if (command == "write_bin_xyz")
    {write_bin_xyz();}
    else if (command == "thresholded_list")
    {thresholded_list();}
    else if (command == "composition")
    {composition();}
    else if (command == "nfold")
    {nfold();}
    else if (command == "vector_autocorrelation_function")
    {vector_autocorrelation_function();}
    else if (command == "crosscorrelate_value_lists")
    {crosscorrelate_value_lists();}
    else if (command == "autocorrelate_value_list")
    {autocorrelate_value_list();}
    else if (command == "clustered_list")
    {clustered_list();}
    else if (command == "invert_list")
    {invert_list();}
    else if (command == "trajectory_list_decay")
    {trajectory_list_decay();}
    else if (command == "gyration_radius")
    {gyration_radius();}
    else if (command == "find_edge")
    {find_edge();}
    else if (command == "skip")
    {skip();}
    else if (command == "exit")
    {cout << "\nExit: Execution terminated at user request.\n";exit(1);}
    else if (command == "break")
    {cout << "\nBreaking out of loop.\n";return 2;}
    else if (command == "print")
    {print();}
    else if (command == "wait")
    {
        cout << "\nWaiting " << args[1] << " seconds" << endl;
        sleep(atoi(args[1].c_str()));
    }
    else if (command == "constant")
    {set_constant(args[1], args[2]);}
    else if (command == "for")
    {if (!do_for_loop()) return 0;}
    else if (command == "while")
    {}//Nothing yet
    else if (command == "if")
    { int error = do_if_statement();
        if (error != 1)
            return error; }
    else if (command == "end" || command == "else")
    {}
    else if (command == "eval" || command == "evaluate")
    {evaluate_expression();}
    else if (command == "round")
    {round_const();}
    else if (command == "floor")
    {floor_const();}
    else if (command == "ceil" || command == "ceiling")
    {ceil_const();}
    else if (command == "user_input")
    {get_user_input(true);}
    else if (command == "processors")
    {change_processors();}
	else if (command == "shell")
	{shell_command();}
    else if (command == "none")
    {}
    else
    {
      Error(string("Command '").append(command)+"' not recognized.", 1);
//      cout << "Command '" << command << "' not recognized.\n";
//      return 0;
    }
  }
    return 1;
}

/*********************************************************************************************************/
/* These functions are public functions for reading parts of the input file vector or seeking through it */
/*********************************************************************************************************/

string Control::read_line()
{
   /** Returns the current line in the input vector and increments the current line number
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    current_line++;
    return get_line(current_line-1);
}

string Control::read_next_line()
{
   /** Returns the line after the current one in the input vector and increments the current line number
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    current_line++;
    return get_line(current_line);
}

string Control::get_line(int lineNum)
{
   /** Returns the line in inputFileVector with lineNum as its index
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    return replace_constants(inputFileVector[lineNum]);
}

string Control::get_raw_line(int lineNum)
{
   /** Returns the line in inputFileVector with lineNum as its index, without replacing any constants
   * @author Michael Marvin
   * @date 7/23/2013
   **/
	if (lineNum < get_input_file_length())
    	return inputFileVector[lineNum];
	else
		return "";
}

int Control::line_seek(int lineNum)
{
   /** Change the current line number to lineNum, returns the old line number
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    int old_line=current_line;
    if (lineNum > get_input_file_length())
    {
        Error("A function tried to seek beyond the end of the input file.", 10002, true);
        current_line=get_input_file_length();
        return old_line;
    }
//	else if (lineNum == get_input_file_length())
//		--lineNum;
    current_line=lineNum;
    return old_line;
}

int Control::get_line_number()
{
   /** Returns the current position of the "cursor" in the inputFileVector
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    return current_line;
}

int Control::get_input_file_length()
{
   /** Returns the size of the inputFileVector
   * @author Michael Marvin
   * @date 8/15/2013
   **/
    return inputFileVector.size();

}

/************************************************************************************/
/* These functions are for mathematical AMDAT commands (such as round and evaluate) */
/************************************************************************************/

void Control::round_const()
{
   /** Input script command for rounding a constant to an int
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    string constant = args[1];
    string const_val = get_constant(constant);
    if (const_val.empty())
    {
        Error("Cannot round constant! Constant "+constant+" not found.", 3);
		return;
    }
    int out = round_float(atof(const_val.c_str()));
    stringstream ss;
    ss << out;
    set_constant(constant, ss.str());
}

int Control::round_float(float f)
{
   /** Rounds a float to the closest integer
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    return (f >= 0) ? (int)(f + 0.5) : (int)(f - 0.5);
}

void Control::floor_const()
{
   /** Input script command for rounding a constant down to an int
   * @author Michael Marvin
   * @date 8/15/2013
   **/
    string constant = args[1];
    string const_val = get_constant(constant);
    if (const_val.empty())
    {
        Error("Cannot floor constant! Constant "+constant+" not found.", 3);
		return;
    }
	float f = atof(const_val.c_str());
	int out = floor(f);
    stringstream ss;
    ss << out;
    set_constant(constant, ss.str());
}

void Control::ceil_const()
{
   /** Input script command for rounding a constant up to an int
   * @author Michael Marvin
   * @date 8/15/2013
   **/
    string constant = args[1];
    string const_val = get_constant(constant);
    if (const_val.empty())
    {
        Error("Cannot ceiling constant! Constant "+constant+" not found.", 3);
		return;
    }
	float f = atof(const_val.c_str());
	int out = ceil(f);
    stringstream ss;
    ss << out;
    set_constant(constant, ss.str());
}

void Control::evaluate_expression()
{
   /** Evaluates a mathematical expression and sets the value inside a constant. Performs calculations linearly, DOES NOT FOLLOW ORDER OF OPERATIONS!
   * Parenthesis support added in version 0.430
   * @author Michael Marvin
   * @date 7/23/2013
   **/

//    string opArray[]={"+", "-", "*", "/", "^", "%", "(", ")"};  // Add new math operators here
//    vector<string> opVect(opArray, opArray+sizeof(opArray)/sizeof(string));

    string out_const = "";
    string line="";
    for (int i=1; i<ARGMAX; i++) //Rebuild the line as a whole string
        line=line+args[i];
    if (line.find("=",0)!=string::npos)
    {
        int pos=line.find("=",0); //The constant is the content before the "=", the expression is after
        out_const=line.substr(0, pos);
        line=line.substr(pos+1);
    }
    else
    {
        Error("Evaluate command requires the \"=\" symbol!", 1);
        return;
    }
    out_const.erase(remove_if(out_const.begin(), out_const.end(), ::isspace), out_const.end()); //remove whitespace from the constant name

	float value=process_expression(line);

    stringstream ss;
    ss << value;
    set_constant(out_const, ss.str());
}

float Control::process_expression(string exp)
{
   /** Evaluates the specific string passed as a mathematical expression. Needed for recursion.
   * @author Michael Marvin
   * @date 10/21/2013
   **/
    string opArray[]={"+", "-", "*", "/", "^", "%"};  // Add new math operators here
    vector<string> opVect(opArray, opArray+sizeof(opArray)/sizeof(string));
    string stack="";
    float stackVal=0.0;
    string nextOp="+";
    for(int i=0; i<exp.length(); i++)
    {
        stringstream s;
        s << exp[i];
        string ch = s.str();
        if (ch == "(")
        {
            int depth=0;
            int len=-1;
            for(int j=i; j<exp.length(); j++)
            {
                stringstream T;
                T << exp[j];
                string tmp = T.str();
                if (tmp=="(")
                {
                    depth++;
                }
                else if (tmp==")")
                {
                    --depth;
                }
                if (depth==0)
                {
                    len=j-i-1;
                    break;
                }
            }
            if (len<0)
            {
                Error("Unbalanced parenthesis in expression!", 4);
				return 0;
            }
            stackVal=eval_terms(nextOp, stackVal, process_expression(exp.substr(i+1,len)));
            stack="";
            i=i+len+1;
        }
        else if (find(opVect.begin(), opVect.end(), ch)!=opVect.end())
        {
            stackVal=eval_terms(nextOp, stackVal, atof(stack.c_str()));

            stack="";
            nextOp=ch;
        }
        else if (ch==" ")
            continue;
        else
        {
            stack=stack+ch;
        }
    }
    if (stack != "")
        stackVal=eval_terms(nextOp, stackVal, atof(stack.c_str())); //Perform the last operation on the remaining term

	return stackVal;
}

float Control::eval_terms(string oper, float a, float b)
{
   /** Performs a mathematical operation on two floats, used primarily for the evaluate command.
   * @author Michael Marvin
   * @date 10/21/2013
   **/
   // cout << "a=" << a << " b=" << b << " op=" << oper << endl;
    if (oper == "+")
        return a+b;
    else if (oper == "-")
        return a-b;
    else if (oper == "*")
        return a*b;
    else if (oper == "/")
        return a/b;
    else if (oper == "^")
        return pow(a,b);
    else if (oper == "%")
        return float(int(a) % int(b));
}

/*************************************************************************/
/* These functions are for logic and loops (if statements and for loops) */
/*************************************************************************/

int Control::do_if_statement()
{
   /** Does an if statement to compare the two given arguments
   * @author Michael Marvin
   * @date 7/23/2013
   **/
    bool isStringCompare=false;
    bool result = false;
    float test1 = atof(args[1].c_str());
    string compare = args[2];
    float test2 = atof(args[3].c_str());
    if (test1 == 0.0 && test2 == 0.0 && !((args[1] == "0" || args[1] == "0.0") && (args[3] == "0" || args[3] == "0.0")))
        isStringCompare=true;
    int initialPos = current_line;
    int endPos = locate_loop_end(initialPos-1);
   // cout << "current: " << current_line << " end: " << endPos << endl;
    //initialPos++;
    if (endPos == -1)
    {
        cout << "\nUnbalanced if statement! Could not find matching \"end\" command!" << endl;
        return 0;
    }

    if (!isStringCompare)
    {
        if (compare == "==")
            result=(test1==test2);
        else if (compare == "<")
            result=(test1<test2);
        else if (compare == ">")
            result=(test1>test2);
        else if (compare == "<=")
            result=(test1<=test2);
        else if (compare == ">=")
            result=(test1>=test2);
        else if (compare == "!=")
            result=(test1!=test2);
    }
    else
    {
        string test1=args[1];
        string test2=args[3];
        if (compare == "==")
            result=(test1==test2);
        else if (compare == "<")
            result=(test1<test2);
        else if (compare == ">")
            result=(test1>test2);
        else if (compare == "<=")
            result=(test1<=test2);
        else if (compare == ">=")
            result=(test1>=test2);
        else if (compare == "!=")
            result=(test1!=test2);
    }
    int error=1;
    if (result) // If the comparison is true
    {
        //cout << "executing if" << endl;
        int end = locate_if_end(initialPos-1, true); // Execute between the "if" and either "else" or "end"
        error = execute_commands(initialPos, end);
    }
    else if (endPos != locate_if_end(initialPos-1, true)) // If there is an "else" statement, this will be true
    {
       // cout << "executing else" << endl;
        int start = locate_if_end(initialPos-1, true); // Execute between "else" and "end"
        int end = locate_if_end(start, false);
        error = execute_commands(start+1, end);
    }
    line_seek(endPos);
    return error;
}


bool Control::do_for_loop()
{
   /** Does a for loop until the appropriate conditions are met. fValue(third argument) is exclusive!
   * @author Michael Marvin
   * @date 7/22/2013
   **/
    string constant = args[1];
    int iValue = atoi(args[2].c_str());
    int fValue = atoi(args[3].c_str());
    int step;
    if (args[4] != "") // If no step argument provided, assume step == 1
        step = atoi(args[4].c_str());
    else
        step = 1;
    int initialPos = current_line;
    int endPos = locate_loop_end(initialPos-1);
    if (endPos == -1)
    {
        cout << "\nUnbalanced for loop! Could not find matching \"end\" command!" << endl;
        return false;
    }
    //initialPos++; //This increment causes the execution to begin one line after the for loop declaration
    if (iValue<=fValue)
    {
        for (int indexii=iValue;indexii<fValue;indexii=indexii+step) // Loop while incrementing indexii to fValue from below by step
        {
            stringstream ss;
            ss << indexii;
            set_constant(constant, ss.str());
            int error = execute_commands(initialPos, endPos);
            if (error == 0)
                return false;
            else if (error == 2)
                break;
        }
    }
    else
    {
        for (int indexii=iValue;indexii>fValue;indexii=indexii-step) // Loop while decrementing indexii to fValue from above by step
        { //Key changes from above: 2nd term has < sign switched to >, third term subtracts by step instead of adds
            stringstream ss;
            ss << indexii;
            set_constant(constant, ss.str());
            int error = execute_commands(initialPos, endPos);
            if (error == 0)
                return false;
            else if (error == 2)
                break;
        }
    }
    line_seek(endPos);
    return true;
}

int Control::locate_if_end(int iIndex, bool isIF)
{
    int depth = 0;
    string line;
    for(int argii=0;argii<ARGMAX;argii++)
    {
      args[argii]="";
    }
    for (int i=iIndex; i<inputFileVector.size(); i++)
    {
        line = get_raw_line(i);
        tokenize(line, args);
        string command = args[0];
        if (isIF)
        {
            if (command == "for" || command == "while" || command == "if")
                depth++;
            else if (command == "end" || command == "else")
            {
                depth--;
                if (depth == 0)
                    return i;
            }
        }
        else
        {
            depth = 1;
            if (command == "for" || command == "while" || command == "if")
                depth++;
            else if (command == "end")
            {
                depth--;
                if (depth == 0)
                    return i;
            }
        }
    }
    return -1;
}

int Control::locate_loop_end(int iIndex)
{
   /** Finds the line number of the "end" command that belongs to the loop at the location passed
   * @author Michael Marvin
   * @date 7/22/2013
   **/
    int depth = 0;
    string line;
    for(int argii=0;argii<ARGMAX;argii++)
    {
      args[argii]="";
    }
    for (int i=iIndex; i<inputFileVector.size(); i++)
    {
        line = get_raw_line(i);
        tokenize(line, args);
        string command = args[0];
        if (command == "for" || command == "while" || command == "if")
            depth++;
        else if (command == "end")
        {
            depth--;
            if (depth == 0)
                return i;
        }
    }
    return -1;
}

/**********************************************/
/* These functions are more general functions */
/**********************************************/

void Control::shell_command()
{
    string command = "";
    for (int i=1; i<ARGMAX; i++)
    {
		//strcat(command, args[i].c_str());
        command=command+" "+args[i];
    }
    //cout << command << endl;
    int status=std::system(command.c_str());
	cout << "Command exited with exit code " << status << endl;

}


void Control::change_processors()
{
   /** Changes the number of processors that AMDAT can use, up to the maximum allotted at the beginning. Useful to ensure that non-threaded analysis don't try to use threads and give incorrect data
   * @author Michael Marvin
   * @date 8/14/2013
   **/
    int threads = atoi(args[1].c_str());
    if (threads <= MAXTHREADS && threads > 0)
    {
        cout << "Setting number of processors to " << threads << endl;
        omp_set_num_threads(threads);
    }
    else if (threads > MAXTHREADS)
    {
        stringstream ss;
        ss<<"Number of processors given (" << threads << ") is higher than the maximum. Instead setting to (" << MAXTHREADS << ")";
        Error(ss.str(), 10001);
        omp_set_num_threads(MAXTHREADS);
    }
    else
    {
        Error(string("Cannot set number of processors to ")+args[1]+". Instead setting to 1.", 10001);
        omp_set_num_threads(1);
    }
}

void Control::get_user_input(bool show_tips)
{
    string input = "";
    int initialEnd = inputFileVector.size();
    int prevCursorPos = get_line_number();
    int linesAdded = 0;
    bool cancelled = false;
    if (show_tips)
    {
        cout << "\nAwaiting user input! Type valid AMDAT commands here, line by line, pressing \"Enter\" at the end of each line. When you are finished, type \"done\" to execute the commands, or \"cancel\" to cancel execution." << endl;
        cout << "(As of now there is little to no error checking so be careful when entering input)" << endl; // TODO: When this is no longer true, remove this warning!
    }
    while (input != "done") //Loop to allow an indefinite amount of commands until "done" is given
    {
        //cin >> input;
        getline(cin, input);
        if (input == "cancel") //Allow the user to cancel input without executing commands
        {
            cancelled = true;
            break;
        }
        if (input != "")
        {
            inputFileVector.push_back(input); //Add the lines to the end of the input file and increment the number of lines added
            linesAdded++;
        }
        //cout << input << endl;
    }
    if (linesAdded == 1 && input == "done")
    {
        inputFileVector.pop_back();
        cout << "Returning to previous execution." << endl;
        //line_seek(prevCursorPos);
        return;
    }
	if (!cancelled)
    {
        cout << "Beginning execution of user input." << endl;
        if (execute_commands(initialEnd, inputFileVector.size()-1))
            cout << "\nUser input successfully executed." << endl;
        else
            cout << "\nAn error occured while executing user input!" << endl;
    }
    else
    {
        cout << "User cancelled execution of input." << endl;
    }

    for (int i=0; i<linesAdded; i++)  // Remove the lines that were added to the script
        inputFileVector.pop_back();
	if (prevCursorPos < get_input_file_length())
    	line_seek(prevCursorPos);

    cout << "Waiting for additional input... (type \"done\" with no other input to return to previous execution)" << endl;

    get_user_input(false);
	if (prevCursorPos < get_input_file_length())
    	line_seek(prevCursorPos);
    return;

}

/*---------------------------------------------------------------------------------*/

 void Control::initialize_lists()
 {
   //list of particle lists
   n_gaussian_comparisons = 0;
   vhs_defined = 0;
   n_trajectory_list_bins = 0;
   n_value_lists=0;

 }


/*Gives error if number of arguments does not match that expected*/
 void Control::argcheck(int expected)
 {
 if(n_args!=expected)
  {
    stringstream ss;
    ss << "Incorrect number of arguments for command "<< args[0] <<".\n"<< n_args-1 << " arguments given, "<< expected-1 << " expected.";
    Error(ss.str(), -6);
  }
}

/*Same as previous function but returns true if the number is correct, and false if incorrect, for error handling. */
bool Control::bool_argcheck(int expected)
{
    if(n_args!=expected)
    {
        stringstream ss;
        ss << "Incorrect number of arguments for command "<< args[0] <<".\n"<< n_args-1 << " arguments given, "<< expected-1 << " expected.";
        Error(ss.str(), 6);
        return false;
    }
    return true;
}


/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/

string Control::replace_constants(string line) {
    string constant_name;
    int constant_num;
    size_t constant_size;
    size_t constant_name_size;
    size_t constant_start=0;
    size_t constant_end=0;
    size_t constant_name_start=0;
    size_t constant_name_end=0;

    while(1) //loop until break command, at end of line
    {
        constant_start=line.find("${", 0); //find start of variable
        constant_end=line.find("}",0); //find end of varible

        constant_size = constant_end-constant_start+1;

        constant_name_start = constant_start + 2;
        constant_name_end = constant_end;

        constant_name_size = constant_name_end - constant_name_start;


        if (constant_start<line.npos)
        {
            constant_name = line.substr(constant_name_start,constant_name_size);
            constant_num = find_constant(constant_name);
            if (constant_num < 0)
            {
//            	cout<< "constant " << constant_name << " is not defined."<< endl; cout.flush();
//	            exit(1);
                Error("Constant \""+constant_name+"\" is not defined.", 3);
				return "";
            }

            line.replace(constant_start, constant_size, constants[constant_num]);

        }
        else
        {
        break;
        }

    }return line;
    }

 /*--------------------------------------------------------------------------------*/



void Control::set_constant(string constant_name, string constant_value) {
  int constant_num;

  constant_num = find_constant(constant_name);
  if(constant_num==-1)
  {
    constant_num = n_constants;
    n_constants++;
  }
    constant_names[constant_num]=constant_name;
    constants[constant_num] = constant_value;
    //cout << "Setting constant \"" << constant_name << "\" to value \"" << constant_value << "\" " << endl;
}


/*--------------------------------------------------------------------------------*/


int Control::find_constant(string constant_name) {

	for(int constantii=0;constantii<n_constants;constantii++)
	{
		if(constant_name==constant_names[constantii])
		    return constantii;
	}
	//cout<< "constant " << constant_name << " is not defined."<< endl; cout.flush();
	//exit(1);
	return -1;
}


/*--------------------------------------------------------------------------------*/
/* Function to get the value of a constant from the constant's name. */
string Control::get_constant(string constant_name) {
    int index=find_constant(constant_name);
    if (index>=0)
        return constants[index];
    else
        return "";
}

/*--------------------------------------------------------------------------------*/


// void Control::new_constant()
// {
   /** adds a variable to be used in input files
   * @author Daniel Hunsicker
   * @date 8/27/2012
   **/

//   string constant;
//   string constant_name;

//   constant_name=args[1];
//   constant_value = args[2];

//   add_constant(constant_name,constant);

// }



/*--------------------------------------------------------------------------------*/



/*Method to run a daughter class of Analysis, passing user instructions for which particles to loop over*/
void Control::run_analysis(Analysis* analyzer, string setline)
{
  string setargs[ARGMAX];		//array of arguments in runline
  int n_setargs;			//number of arguments in runline
  string command;			//command specifying type of set to loop over
  int expected;
  string listname;
  int listnum;

  n_setargs = tokenize(setline, setargs);
  if(n_setargs==0){Error("No atom set command found.", 5);}
  else {command = setargs[0];}

  /*
  if(command == "list")
  {
    expected = 2;
    setargcheck(expected, n_setargs, command);
    listname = setargs[1];
    listnum = find_atomlist(listname);

    if(listnum!=-1)
    {
      analyzer->atomlist(atomlists[listnum]);
    }
  }
  */
  if(command == "list")
  {
    expected = 2;
    setargcheck(expected, n_setargs, command);
    listname = setargs[1];

    if(find_trajectorylist(listname)!=0)
    {
	    analyzer->analyze(find_trajectorylist(listname));
    }
    else
    {
	  cout << "\nTrajectory list '"<<listname<<"' not found.";
    }
  }
  else
  {
    analyzer->analyze(setline);
  }

}


/*--------------------------------------------------------------------------------*/



/*Gives error if number of arguments does not match that expected by run_analysis*/
 void Control::setargcheck(int expected, int n_runargs, string command)
 {
 if(n_runargs!=expected)
  {
    cout << "\nError: Incorrect number of arguments for command "<< command <<".\n"<< n_runargs-1 << " arguments given, "<< expected-1 << " expected.\n";
    exit(1);
  }
 }


 /*--------------------------------------------------------------------------------*/



  /*finds trajectorylist object by custom name*/
Trajectory_List* Control::find_trajectorylist(string listname, bool allow_nofind)const
{
    Trajectory_List * trajectory_list;

  try
  {
    trajectory_list = trajectories.at(listname);
  }
  catch(out_of_range & sa)
  {
    if(allow_nofind)
    {
      trajectory_list=0;
    }
    else
    {
      cout << "\nError: trajectory_list " << listname << " does not exist.\n";
      exit(0);
    }
  }

  return trajectory_list;

}



 /*--------------------------------------------------------------------------------*/



void Control::add_trajectorylist(Trajectory_List * t_list, string listname)
{
 bool result;

  result=(trajectories.insert({listname,t_list})).second;

  if(!result)
  {
    cout << "\nWarning: trajectory_list "<< listname<<" not created because a trajectory_list with this name already exists. Replacement of a trajectory_list requires that you first delete the existing list with the same name.\n";
  }


}





/*--------------------------------------------------------------------------------*/


void Control::combine_trajectories()
{
    int argii;
    string newlistname;

    Trajectory_List * new_trajectory_list;
    new_trajectory_list = new Trajectory_List;

    newlistname=args[1];

    new_trajectory_list=find_trajectorylist(args[2]);

    for(argii=3;argii<n_args;argii++)
    {
        (*new_trajectory_list)=(*new_trajectory_list)||(*find_trajectorylist(args[argii]));
    }

    add_trajectorylist(new_trajectory_list,newlistname);

}



/*--------------------------------------------------------------------------------*/

Multibody_List* Control::find_multibody_list(string listname,bool allow_nofind)const
{
  Multibody_List * multibody_list;

  try
  {
    multibody_list = multibody_lists.at(listname);
  }
  catch(out_of_range & sa)
  {
    if(allow_nofind)
    {
      multibody_list=0;
    }
    else
    {
      cout << "\nError: multibody_list " << listname << " does not exist.\n";
      exit(0);
    }
  }

  return multibody_list;
}




/*--------------------------------------------------------------------------------*/


void Control::add_multibody_list(Multibody_List* multibody_list,string multibody_list_name)
{
  bool result;

  result=(multibody_lists.insert({multibody_list_name,multibody_list})).second;

  if(!result)
  {
    cout << "\nWarning:multibody_list "<< multibody_list_name<<" not created because a multibody_list with this name already exists. Replacement of a multibody_list requires that you first delete the existing list with the same name.\n";
  }
}




/*--------------------------------------------------------------------------------*/

void Control::delete_multibody_list(string listname)
{
  Multibody_List * multibody_list;

  multibody_list = find_multibody_list(listname,1);
  if(multibody_list=0)
  {
    cout << "\nWarning: multibody_list " << listname << " does not exist and therefore cannot be deleted.";
  }
  else
  {
    multibody_lists.erase(listname);
    delete [] multibody_list;
  }
}


/*--------------------------------------------------------------------------------*/



  /*finds trajectorylist object by custom name*/
int Control::find_value_list(string listname)
{
	int listii;

	for(listii=0;listii<n_value_lists;listii++)
	{
		if(listname==value_list_names[listii]) return listii;
	}
	return -1;
}



 /*--------------------------------------------------------------------------------*/



void Control::add_value_list(Value_List<float> * av_list, string listname)
{
  int listnum;

  listnum = find_value_list(listname);
  if(listnum==-1)
  {
    listnum = n_value_lists;
    n_value_lists++;
  }
  value_lists[listnum] = av_list;
  value_list_names[listnum]=listname;

}



 /*--------------------------------------------------------------------------------*/


/*Method to set maximum number of displacement time loops to be used by any analysis to follow*/
 void Control::limit()
{
  int looplimit;
  looplimit = atoi(args[1].c_str());
  analyte->set_limit(looplimit);
}



/*--------------------------------------------------------------------------------*/



void Control::create_list()
{
  string listname;
  string runline;
  int expected=2;
  int n_setargs;
  string setargs[ARGMAX];
  argcheck(expected);

  Static_Trajectory_List* trajectory;
  Trajectory_List * trajpointer;
  //Trajectory_List_Bins * bins;

  listname = args[1];		//user-input name of list

//  getline(input,runline);

  runline = read_line();

  n_setargs=tokenize(runline,setargs);
//   if (setargs[0]=="bin_list")
//   {
//     if(n_setargs!=5)
//     {
//       cout<<"Wrong bin_list syntax, consult manual."<<endl;
//       exit(1);
//     }
//
//     trajpointer = new Trajectory_List;
//     bins = binned_trajectories[find_trajectorylist_bins(setargs[1])];
//     trajpointer = &(*bins)(atoi(setargs[2].c_str())-1,atoi(setargs[3].c_str())-1,atoi(setargs[4].c_str())-1);
//     add_trajectorylist(trajpointer, listname);
//     cout<<"\nTrajectory list "<<listname<<" created from bin "<< atoi(setargs[2].c_str())-1<<atoi(setargs[3].c_str())-1<<atoi(setargs[4].c_str())-1<<"with "<<trajpointer->show_n_trajectories(0)<< " trajectories."<<endl;
//   }
//   else
//   {
    trajectory = new Static_Trajectory_List;
    trajectory->reset(analyte);
    trajpointer=(Trajectory_List*)trajectory;
    run_analysis(trajectory,runline);	//create list from system loops
    add_trajectorylist(trajpointer, listname);	//add trajectory list to array
    cout<<"\nTrajectory list "<<listname<<" created with "<<trajpointer->show_n_trajectories(0)<< " trajectories."<<endl;
//   }
}



/*--------------------------------------------------------------------------------*/




void Control::trajectories_from_multibodies()
{
  string trajectory_list_name, multibody_set_name, centertypename;
  bool centertype;
  Trajectory_Set * trajectory_set_pointer;
  Static_Trajectory_List * new_trajectory_list;

  new_trajectory_list = new Static_Trajectory_List;

  trajectory_list_name = args[1];
  multibody_set_name = args[2];
  centertypename = args[3];

  if(centertypename == "centroid")
  {
    centertype = 0;
  }
  else if(centertypename == "com")
  {
    centertype = 1;
  }
  else
  {
    cout << "\n Type of multibody center '" << centertypename << "' not recognized. Allowable options are 'centroid' and 'com'";
    exit (0);
  }

  trajectory_set_pointer = analyte->create_trajectory_set(trajectory_list_name,multibody_set_name,centertype);

  new_trajectory_list->set(analyte,trajectory_set_pointer);


  add_trajectorylist(new_trajectory_list, trajectory_list_name);

}




/*--------------------------------------------------------------------------------*/
void Control::create_multibodies()
{
  string multibody_list_name;
  Multibody_Set* multibody_set_pointer;
  Multibody_List* new_multibody_list;
  
  new_multibody_list=new Multibody_List;


  multibody_list_name = args[1];
 
  
  multibody_set_pointer = analyte->create_multibody_set (multibody_list_name, n_args, args);    //create multibody set with name that is the same as the multibody list. This is where the multibodies are created.



  new_multibody_list->set(analyte,multibody_set_pointer);
  add_multibody_list(new_multibody_list,multibody_list_name);
  
}




/*--------------------------------------------------------------------------------*/


void Control::write_list_trajectory()
{
  string trajname;
  string listname;
  int expected = 3;
  argcheck(expected);

  listname = args[1];
  trajname = args[2];

  find_trajectorylist(listname)->write_xyz(trajname);
}


/*--------------------------------------------------------------------------------*/


void Control::write_list_trajectory_full()
{
  string trajname;
  string listname;
  int expected = 3;
  argcheck(expected);

  listname = args[1];
  trajname = args[2];

  find_trajectorylist(listname)->write_full_xyz(trajname);
}


/*--------------------------------------------------------------------------------*/

/*Calculate and write to file mean square displacement as requested by user*/
void Control::msd()
{
  string filename;
  string runline;
  int expected=2;
  argcheck(expected);
  dynamic = 1;

  filename = args[1];

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Square_Displacement msd(analyte);
  cout << "\nCalculating mean square displacement.\n";cout.flush();
  start = time(NULL);
  run_analysis <Mean_Square_Displacement> (msd,runline,filename); // pass run_analysis template the analysis type 'Mean_Square_Displacement'

  finish = time(NULL);
  cout << "\nCalculated mean square displacement in " << finish-start<<" seconds."<<endl;
}




/*--------------------------------------------------------------------------------*/


/*Calculate and write to file 2-d mean square displacement as requested by user*/
void Control::msd_2d()
{
  string filename;
  string runline;
  string plane;
  int expected=3;
  argcheck(expected);
  dynamic = 1;

  filename = args[1];
  plane = args[2];

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Square_Displacement_2D msd(analyte,plane);
  cout << "\nCalculating mean square displacement.";
  start = time(NULL);
  run_analysis<Mean_Square_Displacement_2D>(msd, runline, filename);
  finish = time(NULL);
  cout << "\nCalculated mean square displacement in " << finish-start<<" seconds."<<endl;

}


/*--------------------------------------------------------------------------------*/

/*Calculate and write to file mean displacement as requested by user*/
void Control::mean_displacement()
{
  string filename;
  string runline;
  int expected=2;
  argcheck(expected);
  dynamic = 1;

  filename = args[1];

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Displacement md(analyte);
  cout << "\nCalculating mean displacement.\n";cout.flush();
  start = time(NULL);
  run_analysis <Mean_Displacement> (md,runline,filename); // pass run_analysis template the analysis type 'Mean_Square_Displacement'

  finish = time(NULL);
  cout << "\nCalculated mean displacement in " << finish-start<<" seconds.";
}


/*--------------------------------------------------------------------------------*/

/*Calculate and write to file velocity autocorrelation function as requested by user*/
void Control::vacf()
{
  string filename;
  string runline;
  int argcount;
  int n_timesteps;
  argcount = n_args;
  filename = args[1];
  if (argcount==3)
  {
    n_timesteps=atoi(args[2].c_str());
  }


//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Square_Displacement msd(analyte);
  cout << "\nCalculating mean square displacement.";
  start = time(NULL);
  run_analysis(&msd, runline);
  cout << "\nCalculating velocity autocorrelation function.";

  Velocity_Autocorrelation vaf;

  if(argcount==2)
  {
    vaf.initialize(&msd);
  }
  else if(argcount==3)
  {
    vaf.initialize(&msd,n_timesteps);
  }
  else
  {
    cout << "Error: Incorrect number of arguments (" << argcount << ") for command vac_function. This command requires 1 to 2 arguments.\n";
    exit(1);
  }

  finish = time(NULL);
  cout << "\nCalculated velocity autocorrelation function in " << finish-start<<" seconds."<<endl;
  vaf.write(filename);
}



/*--------------------------------------------------------------------------------*/


#ifndef TACC
/*Calculate and write to file velocity autocorrelation function as requested by user*/
void Control::vacf_fourier()
{
  string filename;
  string runline;
  int argcount;
  int n_timesteps;
  argcount = n_args;
  filename = args[1];
  if (argcount==3)
  {
    n_timesteps=atoi(args[2].c_str());
  }


//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Square_Displacement msd(analyte);
  cout << "\nCalculating mean square displacement.";
  start = time(NULL);
  run_analysis(&msd, runline);
  cout << "\nCalculating velocity autocorrelation function.";

  Velocity_Autocorrelation vaf;

  if(argcount==2)
  {
    vaf.initialize(&msd);
  }
  else if(argcount==3)
  {
    vaf.initialize(&msd,n_timesteps);
  }
  else
  {
    cout << "Error: Incorrect number of arguments (" << argcount << ") for command vac_function. This command requires 1 to 2 arguments.\n";
    exit(1);
  }

  cout << "\nPerforming Fourier transform of velocity autocorrelation function.";
  vaf.fourier_transform();

  finish = time(NULL);
  cout << "\nCalculated velocity autocorrelation function in " << finish-start<<" seconds.";
  vaf.write_fourier(filename);
}
#endif


/*--------------------------------------------------------------------------------*/

/*calculate and write to file self van hove as requested by user*/
void Control::calc_vhs()
{
  string filename;
  string runline;

  int n_bins;
  float max_range;
  int expected = 4;
  argcheck(expected);

  filename = args[1];
  max_range = float(atof(args[2].c_str()));
  n_bins = atoi(args[3].c_str());

   //getline(input,inputline);
   //   istringstream iss (inputline,istringstream::in);

  //   iss >> filename >> max_range >> n_bins >> command;

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;
  //analyte->unwrap();		//should already be unwrapped

  vhs.set(analyte, n_bins, max_range);
  cout << "\nCalculating self part of Van Hove correlation function.";
  start = time(NULL);
  run_analysis(&vhs,runline);
  finish = time(NULL);
  vhs_defined=1;

  cout << "\nCalculated self Van Hove in " << finish-start<<" seconds.";
  vhs.write(filename);

}


/*--------------------------------------------------------------------------------*/



/*Calculate and write to file distinct van hove as requested by user*/
void Control::calc_vhd()
{

  string runline;
  float max_range;
  int n_bins;
  string filename;

  int expected = 3;
  argcheck(expected);

  filename = args[1];
  max_range = float(atof(args[3].c_str()));
  n_bins = atoi(args[4].c_str());

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;
  analyte->boxify();

  vhd.set(analyte, n_bins, max_range);
  cout << "\nCalculating distinct Van Hove correlation function.";
  start = time(NULL);
  run_analysis(&vhd, runline);
  finish = time(NULL);
  cout << "\nCalculated distinct Van Hove in " << finish-start<<" seconds.\n";
  cout << "Writing distinct Van Hove to file. ";
  vhd.write(filename);
}


/*--------------------------------------------------------------------------------*/



/*Calculated and write to file total van hove as requested by user*/
void Control::calc_vht()
{
  string filename;
  int expected = 2;

  argcheck(expected);

  cout << "\nCalculating complete Van Hove correlation function.\n";
  filename = args[2];

  vht =  vhs + vhd;

  vht.write(filename);
}



/*--------------------------------------------------------------------------------*/



#ifndef TACC
/*OUTDATED and likely incorrect method to take fourier transform of self van hove*/
void Control::vhs_fourier()
{
  string filename;
  int expected = 2;

  argcheck(expected);

  filename = args[1];
  cout << "\nCalculating radial fourier transform of self van hove.";
  cout << "\nWarning: this method is outdated and likely incorrect.";
  vhs.write_spatial_inverse(filename);
}
#endif


/*--------------------------------------------------------------------------------*/


#ifndef TACC
/*OUTDATED and likely incorrect method to take fourier transform of distinct van hove*/
void Control::vhd_fourier()
{
  string filename;
  int expected = 2;

  argcheck(expected);

  filename = args[1];
  cout << "\nCalculating radial fourier transform of distinct van hove.";
  cout << "\nWarning: this method is outdated and likely incorrect.";
  vhd.write_spatial_inverse(filename);
}
#endif


/*--------------------------------------------------------------------------------*/

void Control::structure_factor()
{
  string filename, command, wavevector_root, symmetry, runline1, runline2, plane, listname1, listname2;
  int expected = 6;
  int listnum1, listnum2;
  bool fullblock=0;
  float max_length_scale = 0;
  dynamic = 0;

if (n_args !=4)
{
     argcheck(expected);
}


  filename = args[1];			//name of file to which to save calculated data
  symmetry = args[2];			//determine if atom sets are the same or different
  plane = args[3];
  fullblock = bool(atoi(args[5].c_str()));
  max_length_scale = atof(args[4].c_str());

  /*Reading wave vectors from file*/
  cout << "\nReading wave vectors from file.\n";cout.flush();
  Wave_Vectors wavevectors(analyte,plane,max_length_scale);

//  getline(input,runline1);
  runline1 = read_line();
  cout <<"\n"<< runline1;


  n_args = tokenize(runline1, args);
  if (args[0]=="list" )
  {
    Structure_Factor struc_fac(analyte,&wavevectors,fullblock);
    listname1 = args[1];
    //listnum1 = find_trajectorylist(listname1);
    cout <<"\n"<< analyte->show_n_trajectories()<<"\t"<<wavevectors.show_n_wavenumbers()<<"\t"<<find_trajectorylist(listname1)->show_n_trajectories(0)<<"\t"<<fullblock<<"\t";cout.flush();

    if (symmetry=="symmetric")
    {
      cout << "\nCalculating structure factor.\n";cout.flush();
      start = time(NULL);
      struc_fac.analyze(trajectories[listname1]);
      finish = time(NULL);
      cout << "\nCalculated structure factor in " << finish-start<<" seconds.\n";
    }
    else if(symmetry=="asymmetric")
    {
//      getline(input,runline2);
      runline2 = read_line();
      cout <<"\n"<< runline2;
      n_args = tokenize(runline2, args);
      listname2 = args[1];
      //listnum2 = find_trajectorylist(listname2);
      cout << "\nCalculating structure factor.\n";cout.flush();
      start = time(NULL);
      //calls bins
      struc_fac.analyze(find_trajectorylist(listname1),find_trajectorylist(listname2));
      finish = time(NULL);
      cout << "\nCalculated structure factor in " << finish-start<<" seconds.\n";
    }
    else
    {
      cout<<"Error: command for structure factor calculation atom set types not understood.\n";
      exit(1);
    }
    struc_fac.write(filename);
  }
  else if (args[0]=="bin_list")
  {
    if (symmetry=="symmetric")
    {
      cout << "\nCalculating structure factor.";cout.flush();
      start = time(NULL);
      Structure_Factor struc_fac(analyte,&wavevectors,fullblock);
      run_analysis <Structure_Factor> (struc_fac,runline1,filename);
      finish = time(NULL);
      cout << "\nCalculated structure factor in " << finish-start<<" seconds.";
    }
    else if (symmetry=="asymmetric")
    {
      cout<<"Error: Binning currently only functions with symmetric S(k).\n";
      exit(1);
    }
    else
    {
      cout<<"Error: command for structure factor calculation atom set types not understood.\n";
      exit(1);
    }
  }
  else
  {
    cout << "\nStructure Factor requires target to be input via trajectory list for symmetric or asymmetric, and if binning is used symmetric only.";
    exit(0);
  }

}



/*--------------------------------------------------------------------------------*/


/*Calculates intermediate scattering function for particles of set 1, considering interactions with particles of set 2 only.  If both set 1 and set 2 are all particles, this gives the usual coherent intermediate scattering function.  A clean division of intermediate scattering functions that can properly be added together to yield the coherent one is given by using all for set 2 and selecting mutually exclusive subsets for various calculations of set 1.*/
void Control::isf()
{
  string filename, runline, symmetry, method, structfac_filename, plane;
  int inner, outer;		//minimum and maximum wave vector indices to calculate
  bool block_parallel;		//determines whether multiple equally space pairs are used between each pair of blocks (1) or if only spacings between the initial element of each block is used (0)

  filename = args[1];			//name of file to which to save calculated data
  method = args[2]; 			//determine method of approach: either specify wavenumber range automatically, or begin by the calculating structure factor and then enter the range manually by user, allowing them to examine structure factor first.
  symmetry = args[3]; 			//determine if atom sets are the same are different
  plane = args[4];
  block_parallel = bool(atoi(args[5].c_str()));


  cout << "\nDetermining wave vectors.";cout.flush();
  Wave_Vectors wavevectors(analyte,plane);

  Wave_Density wavedensity1(analyte,&wavevectors);
  Wave_Density wavedensity2(analyte,&wavevectors);

  if(method == "auto") //if method is auto, simply use inner and outer wavenumber limits from file
  {
    inner = atoi(args[6].c_str());
    outer = atoi(args[7].c_str());

    wavedensity1.set(analyte,&wavevectors,inner,outer);
    wavedensity2.set(analyte,&wavevectors,inner,outer);

    cout << "\nCalculating fourier transform of density for first set of particles.";cout.flush();

//    getline(input,runline);
    runline = read_line();
    cout <<"\n"<< runline;
    run_analysis(&wavedensity1,runline);	//calculate wave density for first set of particles

    if (symmetry=="symmetric")
    {
      wavedensity2 = wavedensity1;
    }
    else if(symmetry=="asymmetric")
    {
      cout << "\nCalculating fourier transform of density for second set of particles.";cout.flush();
//      getline(input,runline);
      runline = read_line();
      cout <<"\n"<< runline;
      run_analysis(&wavedensity2,runline);  //calculate wave density for second set of particles
    }
  }
  else if(method=="manual") //if method is manual, calculate intermediate scattering function and save to file.  Then allow user to input desired wavenumber range after they can view this file.
  {
    cout << "\nCalculating fourier transform of density for first set of particles.";cout.flush();
    structfac_filename = args[7];
//    getline(input,runline);
    runline = read_line();
    run_analysis(&wavedensity1,runline);	//calculate wave density for first set of particles
    if (symmetry=="symmetric")
    {
      wavedensity2 = wavedensity1;
    }
    else if(symmetry=="asymmetric")
    {
      cout << "\nCalculating fourier transform of density for second set of particles.";cout.flush();
//      getline(input,runline);
      runline = read_line();
      cout <<"\n"<< runline;
      run_analysis(&wavedensity2,runline);  //calculate wave density for second set of particles
    }
    else
    {
      cout<<"Error: command for structure factor calculation atom set types not understood.\n";
      exit(1);
    }

    cout << "\nCalculating structure factor.";cout.flush();
    Intermediate_Scattering_Function struct_fac(analyte, 0, &wavedensity1, &wavedensity2);

    struct_fac.write(structfac_filename);

    cout << "Please input the lower and upper wavenumber index bounds, separated by spaces.\n";
    cin >> inner >> outer;
  }
  /*Calculate intermediate scattering function*/
  cout << "\nCalculating intermediate scattering function.";cout.flush();
  Intermediate_Scattering_Function is_fun(analyte, &wavedensity1, &wavedensity2,inner,outer,block_parallel);
  is_fun.write(filename);		//write intermediate scattering function to file

}


/*--------------------------------------------------------------------------------*/



void Control::isf_list()
{
  #ifdef NEVER
  string filename, wavevector_root, runline, method, structfac_filename, listname, plane;
  int inner, outer, listindex;		//minimum and maximum wave vector indices to calculate

  filename = args[1];			//name of file to which to save calculated data
  listname = args[2];
  plane = args[3];


  cout << "\nDetermining wave vectors.";cout.flush();
  Wave_Vectors wavevectors(analyte,plane);

  Wave_Density wavedensity1(analyte,&wavevectors);

  inner = atoi(args[4].c_str());
  outer = atoi(args[5].c_str());

  wavedensity1.set(analyte,&wavevectors,inner,outer);

  cout << "\nCalculating fourier transform of density for first set of particles.";cout.flush();
//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;
  run_analysis(&wavedensity1,runline);	//calculate wave density for first set of particles

  listindex = find_atomlist(listname);
  cout << "\nCalculating intermediate scattering function.";cout.flush();
  Intermediate_Scattering_Function is_fun(analyte, atomlists[listindex], &wavedensity1, inner, outer);
  is_fun.write(filename);		//write intermediate scattering function to file
#endif
}





void Control::isfs()
{
  int inner, outer;
  string filename, wavevector_root, runline, plane;
  bool fblock=0;
  float max_length_scale=0;
  dynamic = 1;

  if(n_args!=6&&n_args!=7)
  {
    cout<<"\nIncorrect number of arguments for isfs. 6,7, or 8 arguments expected.\n";
    exit(0);
  }

  filename = args[1];
  inner = atoi(args[2].c_str());
  outer = atoi(args[3].c_str());
  plane = args[4];
  max_length_scale = atof(args[5].c_str());



  if(n_args==7)
  {
    fblock = bool(atoi(args[6].c_str()));
  }

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;
  start = time(NULL);
  cout << "\nDetermining wave vectors.";cout.flush();
  Wave_Vectors wavevectors(analyte,plane,max_length_scale);

  Incoherent_Scattering_Function isfs(analyte,&wavevectors, inner, outer, fblock);
  //run_analysis(&isfs,runline);
  run_analysis <Incoherent_Scattering_Function> (isfs,runline,filename);
  finish = time(NULL);
  cout << "\nSelf Intermediate Scattering Function calculated in " << finish-start<<" seconds.\n";cout.flush();
}



/*--------------------------------------------------------------------------------*/



/*control method to calculate distribution of u^2 values at a given time separation.*/
void Control::u2dist()
{
  int n_bins, t;
  float maxvalue;
  string runline, filename;
  int expected = 5;

  argcheck(expected);

  filename = args[1];
  n_bins = atoi(args[2].c_str());
  maxvalue = float(atof(args[3].c_str()));
  t = atoi(args[4].c_str());

//  getline(input,runline);		//set which atoms to calculate this for
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped

  DebyeWaller_Dist dw_dist(analyte, n_bins, maxvalue, t);
  run_analysis(&dw_dist,runline);

  dw_dist.write(filename);
}



/*--------------------------------------------------------------------------------*/



/*control method to calculate distribution of 1/u^2 (stiffness) values at a given time separation.*/
void Control::stiffness_dist()
{
  int n_bins, t;
  float maxvalue;
  string runline, filename;
  int expected = 5;

  argcheck(expected);

  filename = args[1];
  n_bins = atoi(args[2].c_str());
  maxvalue = float(atof(args[3].c_str()));
  t = atoi(args[4].c_str());

//  getline(input,runline);		//set which atoms to calculate this for
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped

  Stiffness_Dist stiff_dist(analyte, n_bins, maxvalue, t);
  run_analysis(&stiff_dist,runline);

  stiff_dist.write(filename);
}


/*--------------------------------------------------------------------------------*/



/*control method to calculate distribution of u^N values at a given time separation.*/
void Control::displacement_dist()
{
  int n_bins, t;
  float maxvalue, power;
  string runline, filename;
  int expected = 6;

  argcheck(expected);

  filename = args[1];
  power = atof(args[2].c_str());
  n_bins = atoi(args[3].c_str());
  maxvalue = float(atof(args[4].c_str()));
  t = atoi(args[5].c_str());

//  getline(input,runline);		//set which atoms to calculate this for
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped

  Displacement_Distribution disp_dist(analyte, power, n_bins, maxvalue, t);
  run_analysis(&disp_dist,runline);

  disp_dist.write(filename);
}



/*--------------------------------------------------------------------------------*/



void Control::ngp()
{
  string runline;
  string filename;
  int expected = 2;

  dynamic = 1;

  argcheck(expected);

  filename = args[1];

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();		//should already be unwrapped
  Mean_Square_Displacement msd(analyte);
  run_analysis(&msd, runline);
  Non_Gaussian_Parameter ngpar(analyte, &msd);
  //run_analysis(&ngpar, runline);
  run_analysis <Non_Gaussian_Parameter> (ngpar,runline,filename); // pass run_analysis template the analysis type
  //ngpar.write(filename);
  cout << "\n Peak time index of non-Gaussian parameter is " << ngpar.max() << ".";

}



/*--------------------------------------------------------------------------------*/



void Control::compare_gaussian()
{
  string filename, runline;
  int expected = 2;
  int n_runargs;
  if(!vhs_defined){cout<<"Error: self van hove not yet defined.  Self van hove must be calculated prior to compare_gaussian calculation.\n";exit(1);}

  if(n_args!=expected&&n_args!=3)
  {
	  cout << "\nError: Incorrect number of arguments for gaussian comparison method.\n"<< n_args-1 << " arguments given, 2 or 3 expected.\n";
	  exit(1);
  }


  filename = args[1];
  n_runargs=n_args;

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  //analyte->unwrap();	//should already be unwrapped
  Mean_Square_Displacement msd(analyte);
  run_analysis(&msd, runline);
  cout << "\nCalculating non-Gaussian parameter";
  Non_Gaussian_Parameter ngpar(analyte, &msd);
  run_analysis(&ngpar, runline);
  cout << "\nComparing self Van Hove with Gaussian approximation.";
  gaussian_comparison = new Gaussian_Comparison;

  if(n_runargs==2)
  {
    gaussian_comparison->set(analyte, &ngpar, &vhs, &msd);
  }
  else if(n_runargs==3)
  {
	  gaussian_comparison->set(analyte, atoi(args[2].c_str()), &vhs, &msd);
  }

  gaussian_comparison->write(filename);

  n_gaussian_comparisons=1;


}



/*--------------------------------------------------------------------------------*/



void Control::find_fast()
{
  string filename, runline;
  Fast_Particles * fast_particles;
  string listname;

  Trajectory_List * trajpointer;

  fast_particles = new Fast_Particles;

  trajpointer=(Trajectory_List*)fast_particles;

  listname=args[1];
  filename = args[2];

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  if(n_gaussian_comparisons==0){cout<<"\nError: A guassian comparison analysis myst be defined prior to calculating fast particle list.\n";exit(1);}

  cout << "\nIdentifying 'fast' particles.";
  fast_particles->set(analyte,gaussian_comparison);
  run_analysis(fast_particles, runline);
  fast_particles->write_count(filename);

  add_trajectorylist(trajpointer, listname);	//add trajectory list to array
}



/*--------------------------------------------------------------------------------*/



void Control::radial_debye_waller()
{
  string filename, runline;
  int timeii, n_bins, maxrange;
  float xcenter, ycenter, zcenter;

  filename = args[1];
  timeii = atoi(args[2].c_str());
  n_bins = atoi(args[3].c_str());
  maxrange = atof(args[4].c_str());
  xcenter = atof(args[5].c_str());
  ycenter = atof(args[6].c_str());
  zcenter = atof(args[7].c_str());


//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<< runline;

  Coordinate coordinate(xcenter, ycenter, zcenter);


  Radial_Debye_Waller rdw(analyte, timeii, n_bins, maxrange, coordinate);
  run_analysis(&rdw, runline);
  rdw.write(filename);
}



/*Calculate and write to file distinct van hove as requested by user*/
void Control::strings()
{
	string runline,t_listname;
	int timegap1, timegap2;
	string filename, imgfilename, sigmatrixfilename;
	int expected = 6;
	float threshold;
	argcheck(expected);

	filename = args[1];
	timegap1=atoi(args[2].c_str());
	//timegap2=atoi(args[3].c_str());
	threshold=atof(args[3].c_str());
	sigmatrixfilename=args[4];
	t_listname=args[5];

//	getline(input,runline);
    runline = read_line();
	cout <<"\n"<< runline;

	cout << "\nFinding strings.";
	start = time(NULL);
	Strings stringlist(analyte,timegap1,threshold,sigmatrixfilename);
	run_analysis(&stringlist, runline);
	finish = time(NULL);
	cout << "\nFound strings in " << finish-start<<" seconds.\n";
	cout << "Writing string data to file. ";
	stringlist.write(filename);

	imgfilename=filename;
	imgfilename+="_img.xyz";
	stringlist.write_jump(imgfilename,0);

	cout << "\nGenerating string trajectory list" << endl;

        Trajectory_List * trajpointer;
        trajpointer = new Trajectory_List();

	stringlist.construct_trajectory_list(trajpointer);

        add_trajectorylist(trajpointer, t_listname);	//add trajectory list to array

        cout<<"\nTrajectory list "<<t_listname<<" created.";
}



void Control::rgtensor_stats()
{
  string runline;
  string filename, dist_filename,gyr_rad_filename;
  int nbins1,nbins2;
  float maxrg;
  int expected = 7;
  argcheck(expected);

  filename = args[1];
  dist_filename = args[2];
  nbins1=atoi(args[3].c_str());
  gyr_rad_filename = args[4];
  nbins2=atoi(args[5].c_str());
  maxrg=atof(args[6].c_str());

  //analyte->unwrap();		//should already be unwrapped

//  getline(input,runline);
  runline = read_line();
  cout <<"\n"<<runline;
  RgTensor_Stats rgtensorstats(analyte);
  run_analysis(&rgtensorstats,runline);
  rgtensorstats.write(filename);

  rgtensorstats.calc_rel_asphericity_dist(nbins1);
  rgtensorstats.write_rel_asphericity_dist(dist_filename);
  rgtensorstats.calc_gyration_rad_dist(nbins2,maxrg);
  rgtensorstats.write_gyration_rad_dist(gyr_rad_filename);
}



void Control::displacement_map()
{
  string runline, filename;
  int timegap_index;
  int firstblock, lastblock;
  float max_displacement;

  filename = args[1];
  timegap_index=atoi(args[2].c_str());
  firstblock=atoi(args[3].c_str());
  lastblock=atoi(args[4].c_str());

  if(n_args==5){max_displacement=0;}
  else if(n_args==6){max_displacement=atof(args[5].c_str());}
  else {argcheck(5);}

//  getline(input,runline);
  runline = read_line();
  cout << "\n" << runline;

  Displacement_Map disp_map(analyte, timegap_index, firstblock, lastblock,max_displacement);

  cout << "\nGenerating displacement map.";
  start = time(NULL);
  run_analysis(&disp_map, runline);
  finish = time(NULL);
  cout << "\nGenerated displacement map(s) in " <<finish - start << " seconds.";

  disp_map.write(filename);

}



void Control::write_single_particle()
{
	string filename;
	int trajii;
	int expected = 3;
	argcheck(expected);

	filename = args[1];
	trajii=atoi(args[2].c_str());

	analyte->write_single_particle(trajii,filename);
}



/*Methods to Create and find binned trajectory lists*/
void Control::create_bin_list()
{
  /** Creates new Binned Trajectory List Object
  * @author Mark Mackura
  * @date 4/11/2012
  **/

  cout << "\nCreating Binned Trajectory List";

  string listname;
  string setline;
  int expected=2, xbins, ybins, zbins;
  float xlow,xhigh,ylow,yhigh,zlow,zhigh;
  argcheck(expected); //checks first line of input section for format of create_bin_list <name>
  string setargs[ARGMAX];
  Trajectory_List_Bins * traj_list_bins;
  Trajectory_List_Bins * trajectory_list_bins_pointer;
  string type;
  analyte->boxify();


  listname = args[1];		//user-input name of list
//  getline(input,setline);
  setline = read_line();

     int n_setargs = tokenize(setline, setargs);	//number of arguments in setline

     if (n_setargs==0)
     {
	  cout << "Error: No set type found.";
	  exit(1);
     }
     else
     {
	  type = setargs[0];
     }

     stringstream ss[n_setargs]; //allocate memory to convert input arguments into primitive data types

     if (type == "all")
     {
	  expected = 4;
	  setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <n_xbins> <n_ybins> <n_zbins>
	  for(int setargsii=0;setargsii<n_setargs;setargsii++)
	  {
	    //creates array of setargs strings as stringstream objects
	    ss[setargsii] << setargs[setargsii];
	  }
	  ss[1] >> xbins;
	  ss[2] >> ybins;
	  ss[3] >> zbins;
	  traj_list_bins = new Trajectory_List_Bins(analyte, xbins, ybins, zbins); //creates Trajectory_List_Bins object
     }
     else if (type == "region")
     {
	  expected = 10;
	  setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <n_xbins> <n_ybins> <n_zbins> <xlo> <xhi> <ylo> <yhi> <zlo> <zhi>
	  for(int setargsii=0;setargsii<n_setargs;setargsii++)
	  {
	    //creates array of setargs strings as stringstream objects
	    ss[setargsii] << setargs[setargsii];
	  }
	  ss[1] >> xbins;
	  ss[2] >> ybins;
	  ss[3] >> zbins;
	  ss[4] >> xlow;
	  ss[5] >> xhigh;
	  ss[6] >> ylow;
	  ss[7] >> yhigh;
	  ss[8] >> zlow;
	  ss[9] >> zhigh;
	  //creates Trajectory_List_Bins object for region bounded by dimensions given, if any are == 0 then box boundary is used
	  traj_list_bins = new Trajectory_List_Bins(analyte, xbins, ybins, zbins, xlow, xhigh, ylow, yhigh, zlow, zhigh);
     }
     else if (type == "distance")
     {
         string subtype;
        subtype = setargs[1];

        if (subtype == "trajectory")
        {
            expected = 6;
            float thickness;
            int n_bins;
            Trajectory_List * list_to_bin;
            Trajectory_List * clust_list;
            setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <list to bin> <list to take distace from> <bin thickness>
            for(int setargsii=0;setargsii<n_setargs;setargsii++)
            {
                //creates array of setargs strings as stringstream objects
                ss[setargsii] << setargs[setargsii];
            }
            ss[4] >> thickness;
            ss[5] >> n_bins;
            list_to_bin = find_trajectorylist(setargs[2]);
            clust_list = find_trajectorylist(setargs[3]);
            traj_list_bins = new Trajectory_List_Bins(analyte, thickness, n_bins, list_to_bin, clust_list); //creates Trajectory_List_Bins object
        }
        else if (subtype == "point")
        {
            expected = 8;
            float thickness,x,y,z;
            int n_bins;
            Trajectory_List * list_to_bin;

            setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <list to bin> <list to take distace from> <bin thickness>
            for(int setargsii=0;setargsii<n_setargs;setargsii++)
            {
                //creates array of setargs strings as stringstream objects
                ss[setargsii] << setargs[setargsii];
            }
            ss[6] >> thickness;
            ss[7] >> n_bins;
            ss[3] >> x;
            ss[4] >> y;
            ss[5] >> z;
            Coordinate point(x,y,z);
            list_to_bin = find_trajectorylist(setargs[2]);
            traj_list_bins = new Trajectory_List_Bins(analyte, thickness, n_bins, list_to_bin, point); //creates Trajectory_List_Bins object
        }
        else if (subtype == "plane")
        {

            string direction;
            direction = setargs[2];
            expected = 8;
            float thickness,position;
            int n_bins;
            string plane;
            Trajectory_List * list_to_bin;

            setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <list to bin> <list to take distace from> <bin thickness>
            for(int setargsii=0;setargsii<n_setargs;setargsii++)
            {
                //creates array of setargs strings as stringstream objects
                ss[setargsii] << setargs[setargsii];
            }
            ss[6] >> thickness;
            ss[7] >> n_bins;
            ss[5] >> position;
            list_to_bin = find_trajectorylist(setargs[3]);
            plane = setargs[4];
            traj_list_bins = new Trajectory_List_Bins(analyte, thickness, n_bins, list_to_bin, plane, position, direction); //creates Trajectory_List_Bins object
        }
        else
        {
            cout << "distance binning type "<< subtype <<" not recognized, use trajectory, point or plane."<<endl; cout.flush();
        }
    }
     else
     {
	cout << "Error: bin_list of type " << setargs[0] << " invalid";
	exit(1);
     }
     trajectory_list_bins_pointer = (Trajectory_List_Bins*)traj_list_bins;
     add_trajectorylist_bins(trajectory_list_bins_pointer, listname);	//add trajectory list to array

   cout<<"\nBinned Trajectory list "<<listname<<" created.";
}





void Control::add_trajectorylist_bins(Trajectory_List_Bins * traj_list_bins, string listname)
{
  /** Adds new Binned Trajectory List object to stored values
  * @author Mark Mackura
  * @date 4/11/2012
  **/
  int trajnum;

  trajnum = find_trajectorylist_bins(listname);
  if(trajnum==-1)
  {
    trajnum = n_trajectory_list_bins;
    binned_trajectories.push_back(traj_list_bins);
    n_trajectory_list_bins++;
  }
  else
  {
    binned_trajectories[trajnum] = traj_list_bins;
  }
  trajectory_list_bin_names[trajnum] = listname;
}


int Control::find_trajectorylist_bins(string listname)
{
  /** Finds trajectory_list_bins object by custom name
  * @author Mark Mackura
  * @date 4/11/2012
  **/
    int listii;
    for(listii=0;listii<n_trajectory_list_bins;listii++)
    {
      if(listname==trajectory_list_bin_names[listii])
      {
	return listii;
      }
    }
    return -1;
}


void Control::remove_bin_list()
{
  /** Removes a Binned Trajectory List Object
  * @author Michael Marvin
  * @date 7/17/2013
  **/

  int expected=2;
  argcheck(expected); //checks first line of input section for format of remove_bin_list <name>

  string listname = args[1];		//user-input name of list to remove

  cout << "\nRemoving binned trajectory list " << listname << endl;

  int listii = find_trajectorylist_bins(listname);
  if (listii<0)
  {
      cout << "\nBinned list " << listname << " not found! Cannot remove it!" << endl;
  }
  else
  {
      delete binned_trajectories.at(listii);
      binned_trajectories.erase(binned_trajectories.begin() + listii);
      for (int i=listii; i<n_trajectory_list_bins; i++)
      {
          trajectory_list_bin_names[i]=trajectory_list_bin_names[i+1];
      }
      n_trajectory_list_bins--;
  }
}


void Control::write_bin_xyz()
{
  /** Writes xyz file for all bins
  * @author Mark Mackura
  * @date 5/22/2012
  **/

  cout << "\nCreating xyz files for bins" << endl;

  string listname;
  string setline;
  string file;
  int expected=3;
  argcheck(expected); 				//checks first line of input section for format of write_bin_xyz <bin_list_ID> <filename>
  string setargs[ARGMAX];
  string type;
  listname = args[1];				//user-input name of bin_list
  file = args[2];
//  getline(input,setline);
  setline = read_line();
  int n_setargs = tokenize(setline, setargs);	//number of arguments in setline
  if (n_setargs==0)
  {
      cout << "Error: No set type found.";
      exit(1);
  }
  else
  {
      type = setargs[0];
  }
  if (type == "all")
  {
      expected = 1;
      setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type>
      binned_trajectories[find_trajectorylist_bins(listname)]->write_bins_xyz(file);
  }
  else if (type == "single")
  {
      int xii,yii,zii;
      expected = 4;
      setargcheck(expected, n_setargs, args[0]); // checks second line of input for size corresponding to format of <type> <xii> <yii> <zii>
      stringstream ss[n_setargs]; //allocate memory to convert input arguments into primitive data types
      for(int setargsii=0;setargsii<n_setargs;setargsii++)
      {
	//creates array of setargs strings as stringstream objects
	ss[setargsii] << setargs[setargsii];
      }
      ss[1] >> xii;
      ss[2] >> yii;
      ss[3] >> zii;
      binned_trajectories[find_trajectorylist_bins(listname)]->write_single_bin_xyz(file,xii,yii,zii);
  }
  else
  {
    cout << "Error: Type " << setargs[0] << " invalid";
    exit(1);
  }
}


void Control::nfold()
{
  /** calculates n-fold orientational order parameter and saves to an analysis value list
  * @author Daniel Hunsicker
  * @date 6/13/2012
  **/
  string filename;
  string runline;
  string sigma_file;
  string f_stem;
  string orient;
  int ord;
  int s_time;
  int e_time;
  float cutoff;
  int expected=10;
  string nfold_listname="n_fold";
  if (n_args==expected)
{
  filename = args[1];
  ord = atof(args[2].c_str());
  sigma_file = args[3];
  orient = args[4];
  cutoff = atof(args[5].c_str());
  nfold_listname = args[6];
  f_stem = args[7];
  s_time = atoi(args[8].c_str());
  e_time = atoi(args[9].c_str());
}
else if (n_args == 7)
{
    filename = args[1];
    ord = atof(args[2].c_str());
    sigma_file = args[3];
    orient = args[4];
    cutoff = atof(args[5].c_str());
    nfold_listname = args[6];
    s_time = 0;
    e_time = -1;
}
else
{
    argcheck(expected);
}
//  getline(input,runline);
  runline = read_line();
  //analyte->unwrap();	//should already be unwrapped

  N_Fold_Order_Parameter* nfold;

  nfold = new N_Fold_Order_Parameter(analyte, ord, sigma_file, orient, cutoff, f_stem, s_time, e_time);



  cout << "\nCalculating n_fold order parameter.\n";
  start = time(NULL);
  run_analysis(nfold, runline);
  add_value_list(nfold,nfold_listname);



  finish = time(NULL);



  cout << "\nCalculated n_fold order parameter in " << finish-start<<" seconds.";
  nfold->write(filename);

}



void Control::thresholded_list()
{
  /** creates a thresholded trajectory list from an analysis value list
  * @author Daniel Hunsicker
  * @date 6/13/2012
  **/
  string av_listname;
  string t_listname;
  //string runline;
  int av_listnum;
  float threshold1;
  int expected=5;
  float threshold2;
  string thresh_command;



  av_listname = args[1];
  t_listname = args[2];//user-input name of list
  thresh_command = args[3];
  threshold1 = atof(args[4].c_str());
  if (n_args==5){threshold2=0;}
  else if(n_args==6){threshold2=atof(args[5].c_str());}
  else {argcheck(expected);}

  //getline(input,runline);



  cout << "\nGenerating thresholded trajectory list" << endl;

  Trajectory_List * trajpointer;
   trajpointer = new Trajectory_List();
  av_listnum = find_value_list(av_listname);

if(thresh_command=="greater"){value_lists[av_listnum]->construct_t_list(bool(1),threshold1,trajpointer);}
else if(thresh_command=="less"){value_lists[av_listnum]->construct_t_list(bool(0),threshold1,trajpointer);}
else if (thresh_command=="between"){value_lists[av_listnum]->construct_t_list(threshold1,threshold2,trajpointer);}
else {cout<< "thresholding commmand unrecognized. command can only be greater, less or between\n";}

  add_trajectorylist(trajpointer, t_listname);	//add trajectory list to array

   cout<<"\nTrajectory list "<<t_listname<<" created with "<<trajpointer->show_n_trajectories(0)<< " trajectories.";

}


void Control::composition()
{
  /** returns the average and time dependent composition
  * @author Daniel Hunsicker
  * @date 6/13/2012
  **/
    string runline;
  string filename;
  int expected=2;

  argcheck(expected);
  dynamic=0;
  filename=args[1];

  int num_xbins, num_ybins, num_zbins;
  float lx, ly, lz;

  string setargs[ARGMAX];		//array of arguments in runline
  int n_setargs;			//number of arguments in runline
  string command;			//command specifying type of set to loop over
//  int bin_expected;
  string listname;

  string bin_listname;
  int bin_listnum;

//  getline(input,runline);
  runline = read_line();
  n_setargs = tokenize(runline, setargs);
   if ( n_setargs==0 )
     {
	  cout << "Error: No atom set command found.";
	  exit(1);
     }
   else
     {
	  command = setargs[0];
     }

  cout << endl << command << endl; cout.flush();
  if ( command=="bin_list")
  {

     if(n_setargs<3)
	  {
	    cout << "Error: Insufficient number of arguments in bin_list target line.\n";
	    exit(0);
	  }
    //bin_expected = 3; // <command> <bin_list_ID> <list_ID>
     //setargcheck(bin_expected, n_setargs, command);
     bin_listname = setargs[1];
     bin_listnum = find_trajectorylist_bins(bin_listname);
     if(bin_listnum!=-1)
      {
	num_xbins = binned_trajectories[bin_listnum]->show_n_xbins();
	num_ybins = binned_trajectories[bin_listnum]->show_n_ybins();
	num_zbins = binned_trajectories[bin_listnum]->show_n_zbins();
	lx = binned_trajectories[bin_listnum]->show_lx();
	ly = binned_trajectories[bin_listnum]->show_ly();
	lz = binned_trajectories[bin_listnum]->show_lz();
      }
      else
      {
	  cout << "\nBinned Trajectory list '"<<bin_listname<<"' not found.";
	  exit(1);
      }
  }
  else
  {
     num_xbins = 1;
     num_ybins = 1;
     num_zbins = 1;
     lx = (analyte->size()).show_x();
     ly = (analyte->size()).show_y();
     lz = (analyte->size()).show_z();
  }
     Composition comp(analyte,num_xbins,num_ybins,num_zbins,lx,ly,lz);
     cout << "\nCalculating composition.\n"; cout.flush();
     start = time(NULL);
     run_analysis <Composition> (comp, runline, filename);

     finish = time(NULL);

     cout << "\nCalculated composition in " << finish-start<<" seconds.";
}


void Control::find_edge()
{
  string filename;
  string runline;
  int expected=5;
  argcheck(expected);
  dynamic = 0;
  float vectorx,vectory,vectorz;
  Coordinate vector;
  
  filename = args[1];
  vectorx = atof(args[2].c_str());
  vectory = atof(args[3].c_str());
  vectorz = atof(args[4].c_str());
  vector.set(vectorx,vectory,vectorz);
  
  runline = read_line();
  cout <<"\n"<< runline;
  
  Edge_Detector_Timedependent edge(analyte,vector);
  cout << "\nFinding edge.\n";cout.flush();
  start = time(NULL);
  run_analysis <Edge_Detector_Timedependent> (edge,runline,filename); // pass run_analysis template the analysis type 'Mean_Square_Displacement'

  finish = time(NULL);
  cout << "\nFound edges in " << finish-start<<" seconds."<<endl;
}



void Control::clustered_list()
{
  /** creates a clusted trajectory list from a trajectory list
  * @author Daniel Hunsicker
  * @date 11/14/2012
  **/

  string traj_listname;
  string t_listname;
  int traj_listnum;
  int expected=7;
  string sig_file;
  string plane;
  int primary;
  int secondary;

  argcheck(expected);

  traj_listname = args[1];
  t_listname = args[2];//user-input name of list
  sig_file = args[3];
  plane = args [4];
  primary = atoi(args[5].c_str());
  secondary = atoi(args[6].c_str());

  //traj_listnum = find_trajectorylist(traj_listname);

  cout << "\nGenerating clustered trajectory list" << endl;

  Trajectory_List * trajpointer;
  trajpointer = new Trajectory_List();



  Clustered_List * clusters;
  clusters = new Clustered_List(*find_trajectorylist(traj_listname));


  start = time(NULL);
  clusters->construct_clust_list(trajpointer,sig_file,plane,primary,secondary);
  finish = time(NULL);
  cout << "\nGenerated clustered trajectory list in " << finish-start<<" seconds.";



  add_trajectorylist(trajpointer, t_listname);	//add trajectory list to array

   cout<<"\nTrajectory list "<<t_listname<<" created with "<<trajpointer->show_n_trajectories(0)<< " trajectories.";

}


void Control::invert_list()
{
  /** creates a trajectory list that is the inversion of a trajectory list times another trajectory list.
  * @author Daniel Hunsicker
  * @date 11/14/2012
  **/
  string traj_listname;
  string t_listname;
  string original_trajectory;
  int traj_listnum;
  int original_traj_listnum;
  int expected=4;



   argcheck(expected);

   traj_listname = args[1];
   original_trajectory= args[2];
   t_listname = args[3];  //user-input name of list

   Trajectory_List * trajpointer;
  trajpointer = new Trajectory_List();

   //traj_listnum = find_trajectorylist(traj_listname);
   //original_traj_listnum = find_trajectorylist(original_trajectory);

   //trajectories[traj_listnum]->inversion(trajpointer, trajectories[original_traj_listnum]);
    trajectories[traj_listname]->inversion(trajpointer, find_trajectorylist(original_trajectory));

   add_trajectorylist(trajpointer, t_listname);	//add trajectory list to array
cout<<"\nTrajectory list "<<t_listname<<" created with "<<trajpointer->show_n_trajectories(0)<< " trajectories.";
}





void Control::vector_autocorrelation_function()	//calculates autocorrelation function of vector orientations
{
  int expected = 3;
  string filename;
  string vectorlist_filename;
  argcheck(expected);

  filename=args[1];
  vectorlist_filename=args[2];

  cout << "\nCalculating vector autocorrelation function.\n";

  Vector_Autocorrelation vector_autocorrelation(analyte,vectorlist_filename);

  vector_autocorrelation.write(filename);

}

void Control::crosscorrelate_value_lists()
{
  int expected = 5;
  string filename;
  string listname1, listname2;
  string correlation_type;
  argcheck(expected);
  int list1, list2;
  float correlation=0;

  correlation_type = args[1];
  filename = args[2];
  listname1 = args[3];
  listname2 = args[4];

  list1 = find_value_list(listname1);
  list2 = find_value_list(listname2);

  ofstream output;

  if(correlation_type=="static")
  {
    correlation = value_lists[list1]->static_crosscorrelation(*value_lists[list2]);
    output.open(filename.c_str());
    output << "Correlation between value lists "<< listname1<<" and "<< listname2<<" calculated by AMDAT v." << VERSION << "\n";
    output << "Relative correlation (<AB>/<A><B>) is " << correlation;
    output.close();
  }
  else if(correlation_type=="dynamic")
  {
    output.open(filename.c_str());
    output << "Dynamic cross-correlation between value lists "<< listname1<<" and "<< listname2<<" calculated by AMDAT v." << VERSION << "\n";
    output.close();
    value_lists[list1]->dynamic_crosscorrelation_function(filename,*value_lists[list2]);
  }
  else
  {
    cout << "Correlation type not understood. Command ignored.\n";
  }

}


void Control::autocorrelate_value_list()
{
  int expected = 3;
  string filename;
  string listname;

  argcheck(expected);

  int list;
  //float correlation=0;

  filename = args[1];
  listname = args[2];

  list = find_value_list(listname);

  ofstream output(filename.c_str());
  output << "Dynamic auto-correlation of value list "<< listname<<" calculated by AMDAT v." << VERSION << "\n";
  output.close();
  value_lists[list]->dynamic_autocorrelation_function(filename);
}


void Control::skip()
{
    string runline;
    getline(input, runline);

    while (runline != "endskip")
    {
        getline(input, runline);
    }
}

void Control::print()
{
  /** print the supplied text to the terminal
  * @author Michael Marvin
  * @date 7/18/2013
  **/
	//cout << "\n";
	for (int i=1;i<n_args;i++)
	{
		cout << args[i] << " ";
	}
	cout <<  endl;
}

void Control::trajectory_list_decay()
{
    string filename,t_list_name,runline;
    dynamic=1;

    int expected=2;
    if (n_args==expected)
    {
        filename = args[1];
    }

  cout << "\nCalculating trajectory list decay" <<endl;
//    getline(input, runline);
    runline = read_line();

    Trajectory_List_Decay t_list_decay(analyte);

    start = time(NULL);
    run_analysis <Trajectory_List_Decay> (t_list_decay, runline, filename);
    finish = time(NULL);

    cout << "\nCalculated trajectory list decay in " << finish-start<<" seconds.";
}



void Control::gyration_radius()
{
  string filename, multibody_list_name;
  Multibody_List * multibodylist;

  int expected=3;
  argcheck(expected);

  filename = args[1];
  multibody_list_name=args[2];

  multibodylist = find_multibody_list(multibody_list_name);

  Gyration_Radius gyrrad(analyte);
  cout << "\nCalculating gyration radius.\n";cout.flush();
  start = time(NULL);
  gyrrad.analyze(multibodylist); // pass run_analysis template the analysis type 'Mean_Square_Displacement'
  finish = time(NULL);
  cout << "\nCalculated gyration radius in " << finish-start<<" seconds."<<endl;
  gyrrad.write(filename);

}
