/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Control class*/
/*Written by David S. Simmons*/
/*Contributions by Mark Mackura and Daniel Hunsicker*/


#ifndef CONTROL
#define CONTROL

#include <vector>
#include <unordered_map>
#include <time.h>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "system.h"
#include "analysis.h"
#include "van_hove_self.h"
#include "van_hove_distinct.h"
#include "gaussian_comparison.h"
#include "trajectory_list_bins.h"
#include "tokenize.h"
#include "bin_dynamics_analysis.h"
#include "bin_static_analysis.h"
#include "value_list.h"
#include "multibody_list.h"



#define ARGMAX 100

#define LISTSIZE 100

namespace std{

class Control
{
    /*Members involved with input file read-in*/
    string args [ARGMAX];		//current arguments
    Tokenize tokenize;		//object to help parse input
    int n_args;				//current number of arguments
    static vector<string> inputFileVector; // Vector that holds the lines in the input file
    int read_input_file(char *); // Read the input file into the vector
    int execute_commands(int, int); // Execute all commands in the input file vector with indices between the two ints
    static int current_line; // A holder for where the "cursor"  is in the input file vector
    ifstream input;
    static string * constants;                //Array of environment constants
    static string * constant_names;      //Array of custom environment constant names
    static int n_constants;
    static int find_constant(string);		//return index of constant with given custom name
    void set_constant(string, string);      // Sets a constant to a given value
    static string replace_constants(string); // A function that replaces the constants in the given string with their values
    string get_constant(string);             //Gets the value of a constant by name
    
    /*Methods enabling advanced scripting in input file*/
    /*Michael Marvin*/
    void shell_command(); 			// executes a Linux command
    void print();                               // "prints" the arguments to screen
    bool do_for_loop();                         // executes a for loop from the current line to the appropriate "end" command
    int locate_loop_end(int);                   // locates the line number of the appropriate "end" command for a loop
    int do_if_statement();                      // executes an if statement from the current line to the appropriate "end" command
    int locate_if_end(int, bool);               // locates the line number of the appropriate "end" command for if statements. The bool is whether it should look from an "if" (true) or "else" (false)
    int round_float(float);                     // rounds a float to the closest integer
    void round_const();                               // Input call to round a constant to an integer
    void floor_const();                               // Input call to floor a constant to an integer
    void ceil_const();                               // Input call to ceiling a constant to an integer
    void evaluate_expression();                 // Performs basic mathematical operations and saves it into a constant
    float process_expression(string);
	float eval_terms(string, float, float);
    void get_user_input(bool);               // Pauses execution and allows the user to input commands on-the-fly
    
    
    /*Members dealing with parallelization*/
    int MAXTHREADS;
    void change_processors();


    System * analyte;			//system being analyzed
    bool vhs_defined;			//stores whether the self van hove has been calculated
    Van_Hove_Self vhs;			//self van hove for system
    Van_Hove_Distinct vhd;		//distinct van hove for system
    Space_Time_Correlation_Function vht;
    Gaussian_Comparison * gaussian_comparison;		//array of Gaussian comparison objects
    int n_gaussian_comparisons;
    time_t start;			//timer start
    time_t finish;			//timer stop


    /*Arrays to store analysis results with a name given by the user, for later recall and use in other analysis techniques*/
    //To be added
    //Space-Time_Correlation_Function * space-time_correlations [LISTSIZE];
    //string space-time_correlationnames [LISTSIZE];
    //int n_space-time_correlations;
    //int find_space-time_correlations();
    //void add_space-time_correlations();

    //void * analyses [LISTSIZE];			//declare
    //string analysis_names [LISTSIZE];
    //int n_analyses;


    /*Members to store and access trajectory_list objects*/
    unordered_map <string, Trajectory_List*> trajectories;
    Trajectory_List* find_trajectorylist(string, bool allow_nofind=0)const;
    void add_trajectorylist(Trajectory_List*, string);
    void trajectories_from_multibodies();
    void combine_trajectories();

    /*Members to store and access multibody_list objects*/
    unordered_map<string,Multibody_List*> multibody_lists;
    Multibody_List* find_multibody_list(string,bool allow_nofind=0)const;
    void add_multibody_list(Multibody_List*,string);
    void delete_multibody_list(string);
    
    /*Some general methods*/
    void system();			//create system object with input file data
    void argcheck(int);			//check if the number of arguments to a analysis method is correct
    bool bool_argcheck(int);    //check if the number of arguments to a analysis method is correct and returns true if it is
    void limit();			//set limit of how many time spacings to use per timegap
    void create_list();
    void new_constant();
    void write_list_trajectory();
    void write_list_trajectory_full();


    void create_multibodies();		//method to create set of user-defined multibodies from molecule, together with a list of multibodies enabling access to this set.

    /*Method to parse input for analysis class atom loop types*/
    void run_analysis(Analysis*, string);
    void setargcheck(int, int, string);		//check number of arguments for above method


    /*Daniel Hunsicker*/
    void initialize_lists();		//initialize some arrays

    /*Ryan Lang*/
    bool dynamic;		//Create boolean to choose dynamic or static binning method.  Value of 1 denotes dynamic, value of 0 denotes static.

    /*Mark Mackura*/
    /* Methods and variables to handle bin lists */
    void create_bin_list();
    template <class Analysis_type> void run_analysis(Analysis_type, string, string);
    void add_trajectorylist_bins(Trajectory_List_Bins *,string);
    int find_trajectorylist_bins(string);
    string trajectory_list_bin_names[LISTSIZE];				//custom name of binned trajectory list
    //Trajectory_List_Bins * binned_trajectories [LISTSIZE];		//array of binned trajectory list objects
    vector<Trajectory_List_Bins *> binned_trajectories;
    int n_trajectory_list_bins;						//number of binned trajectory list objects stored
    void remove_bin_list();                     // Removes a bin list from memory
    void write_bin_xyz();						//writes xyz file for binned_trajectories (all or single)
    void skip();

    /*Methods to handle value lists*/
    /*Daniel Hunsicker*/
    Value_List<float> * value_lists [LISTSIZE];		//array of value list objects
    string value_list_names [LISTSIZE];		//custom name of value list
    int n_value_lists;

    int find_value_list(string);		//return index of value list with given custom name
    void add_value_list(Value_List<float>*, string);

    /*Analysis method calls*/
    void msd();			//calculate mean square displacement
    void msd_2d();		//calculate 2d mean square displacement
    void mean_displacement();	//calculate mean displacement - primary purpose is identification of momentum buildup
    void vacf();		//calculate velocity autocorrelation function
    void vacf_fourier();	//calculate fourier transform of velocity autocorrelation function
    void calc_vhs();		//calculate self van hove
    void calc_vhd();		//calculate distinct van hove
    void calc_vht();		//calculate total van hove
    void vhs_fourier();		//calculate fourier transform of self van hove
    void vhd_fourier();		//calculate fourier transform of distinct van hove
    void isf();			//calculate coherent (total) intermediate scattering function
    void isf_list();		//same as above but works with particle list instead of atom set
    void isfs();		//calculate self intermediate scattering function
    void structure_factor();	//calculate structure factor or intermediate scattering function as function of wavenumber at a given displacement time.
    void u2dist();		//calculate distribution of debye-waller factor (or displacement)
    void stiffness_dist();	//calculate distribbution of inverse debye-waller factor
    void displacement_dist();	//calculated distribution of particle displacement magnitudes to user-specified power
    void ngp();			//calculate non-gaussian parameter
    void compare_gaussian();	//compare self van hove to gaussian approximation for self van hove to find fast particles
    void find_fast();		//find fast particles
    void radial_debye_waller();	//calculate Debye-Waller factor as a function of distance from the origin
    void strings();		//find strings
    void rgtensor_stats();	//calculate rg tensor statistics
    void nfold();		//calculate n-fold orientational order parameter
    void composition();     	//writes a file detailing the compostion of the system
    void thresholded_list();	//generates a list thresholded from an analysis value list
    void clustered_list();	//generates a list clustered from a trajectoory list
    void invert_list();     //inverts a trajectory_list
    void displacement_map();
    void trajectory_list_decay();
    void vector_autocorrelation_function();	//calculates autocorrelation function of vector orientations
    void crosscorrelate_value_lists();
    void autocorrelate_value_list();
    void write_single_particle();	//write single trajectory to file as simple list of coordinates
    void find_edge();
    void unsteady_velocity();
    /*Multibody analysis method calls*/
    void gyration_radius();		//calculate mean gyration radius of multibody list
    void incremental_mean_displacement();	//calculates mean displacement between adjacent frames

  public:
    Control(char *,string*,string*,int,string);		//constructor
    /*  Input file manipulation funtions    */
    static string read_line();  // read the current line and advance the cursor
    static string read_next_line(); // advance the cursor and read the new line
    static string get_line(int); // get the line at the passed index
    static string get_raw_line(int); // get the line at the passed index but without converting constants
    static int line_seek(int); // change the cursor position in the file vector, returns old position
    static int get_line_number(); // get the current cursor position in the file vector
	static int get_input_file_length(); // gets the length of the input file, in lines

};


template <class Analysis_type>
void Control::run_analysis(Analysis_type analyzer, string setline, string filename)
{
     /** Method to run analysis_type as requested by user input file with bin support writing out results as calculations finish
     * @param Analysis_type - analysis type (Analysis child) to be used when duplicating object for bins
     * @param analyzer - analysis object 'template' to be copied for bins
     * @param setline - input data necessary for analysis tool
     * @param filename - output filename stem
     * @author Mark Mackura
     * @date 6/06/2012
     **/

     int n_setargs;			//number of arguments in runline
     string command;			//command specifying type of set to loop over
     int expected;
     string listname, listname2;
     bool use_persistence=0;

     string bin_listname;
     int bin_listnum;
     Trajectory_List_Bins * binned_trajectory_list_to_analyze;		//declare binned trajectory list to analyze

     Tokenize tokenize(setline);
     n_setargs = tokenize.count();
     //n_setargs = tokenize(setline, setargs);
     if ( n_setargs==0 )
     {
	  cout << "Error: No atom set command found.";
	  exit(1);
     }
     else
     {
	  command = tokenize(0);
	  cout <<command<< endl;cout.flush();
     }
     if(command == "list")
     {
       if(n_setargs==2)
       {
	  listname = tokenize(1);
	  //listnum = find_trajectorylist(listname);

	  if(trajectories.count(listname))
	  {
		  analyzer.analyze(trajectories[listname]);
		  analyzer.write(filename);
	  }
	  else
	  {
	  	cout << "\nTrajectory list '"<<listname<<"' not found.";
		exit(1);
	  }
       }
       else if(n_setargs==3)
       {
	  listname = tokenize(1);
	  //listnum = find_trajectorylist(listname);
	  listname2 = tokenize(2);
	  //listnum2 = find_trajectorylist(listname2);

	  if(trajectories.count(listname)&&trajectories.count(listname2))
	  {
		  analyzer.analyze(trajectories[listname],trajectories[listname2]);
		  analyzer.write(filename);
	  }
	  else
	  {
	  	cout << "\nTrajectory list '"<<listname<<"' not found.";
		exit(1);
	  }
       }
     }
     else if (command == "bin_list")
     {
	  //expected = 3; // <command> <bin_list_ID> <list_ID>
	  if(n_setargs<3)
	  {
	    cout << "Error: Insufficient number of arguments in bin_list target line.\n";
	    exit(0);
	  }
	  //setargcheck(expected, n_setargs, command);
	  bin_listname = tokenize(1);
	  bin_listnum = find_trajectorylist_bins(bin_listname);
	  listname = tokenize(2);
	  if (n_setargs==4)
	  {
	    use_persistence=bool(atoi(tokenize(3).c_str()));
	  }
	  //listnum = find_trajectorylist(listname);
	  if (trajectories.count(listname))
	  {

		 binned_trajectory_list_to_analyze = binned_trajectories[bin_listnum];
		 if(dynamic&&use_persistence)
		 {
		    Bin_Dynamics_Analysis <Analysis_type> bin_analyzer(binned_trajectory_list_to_analyze,analyzer,analyte);	// creates analysis array of desired analysis type using template declared above
		    bin_analyzer.analyze(trajectories[listname]); 							// performs analysis on intersection of bins given by <bin_list_ID> and trajectories given by <list_ID>
		    bin_analyzer.write(filename);									// writes out data for each bin individually
		 }
		 else
		 {
		    Bin_Static_Analysis <Analysis_type> static_bin_analyzer(binned_trajectory_list_to_analyze,analyzer,analyte);	// creates analysis array of desired analysis type using template declared above
		    static_bin_analyzer.analyze(trajectories[listname]); 							// performs analysis on intersection of bins given by <bin_list_ID> and trajectories given by <list_ID>
		    static_bin_analyzer.write(filename);									// writes out data for each bin individually
		 }
	   }
	   else
	   {
		    cout << "\nBinned Trajectory list '"<<bin_listname<<"' not found.";
		    exit(1);
	   }
	}	 
     else
     {
     cout << "\nError: Command " << command << " not understood. Acceptable options are list or bin_list." << endl;
     exit(10);
     //analyzer.analyze(setline);
     //analyzer.write(filename);
     }

}


}

#endif
