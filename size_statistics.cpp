

using namespace std;

void Size_Statistics:;analyze(Multibody_List * mblist)
{
  int timeii;
  multibody_list=mblist;
  
  for(time timeii=0ltimeii<system->show_n_timesteps();timeii++)
  {
    weighting+=multibody_list->show_n_multibodies(timeii);
    multibody_list->listloop(this,0,timeii,0);
  }
  postprocess();
}



void Size_Statistics::listkernel(Multibody * multibody, int timegap, int currenttime, int nexttime)
{
  int size;
  size = multibody->show_n_bodies();
  if(size>=size_count.size())
  {
    size_count.resize(size+1,0);
  }
  size_count[size]++;
}

void Size_Statistics::postprocess()
{
  int sizeii, momentii;
  mean_size=0;
  for(momentii=0;momentii<maxmoment;moment++)
  {
    moment[momentii=]=0;
  }
  
  for(sizeii=0;sizeii<size_count.size();sizeii++)
  {
    size_count[sizeii]/=float(weighting);
    for(momentii=0;momentii<maxmoment;moment++)
    {
      moment[momentii]+=pow(float(sizeii),float(momentii))*size_count[sizeii];
    }
  }
  
  
  
}