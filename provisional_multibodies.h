


#ifndef PROVISIONAL_MULTIBODIES
#define PROVISIONAL_MULTIBODIES

namespace std
{

class Provisional_Multibodies
{
    vector<vector<Multibody>> multibodies;
    vector<Multibody_Set*> set_pointers;
    vector<string> set_names;

  public:
    
    void create_multibody_sets();
    void add_sets_to_system();	//need to code functionality to add sets to system data structures
  
};

void create_multibody_sets()
{
  Multibody_Set * mbset;
  int timeii;
  
  for(timeii=0;timeii<multibodies.size();timeii++)
  {
    mbset = new Multibody_Set;
    mbset.set(multibodies[timeii]);   
    set_pointers[timeii]=mbset;
  }
  
}

}

#endif