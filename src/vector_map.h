// Copyright (C) 2010-2025 David S. Simmons

/*Amorphous Molecular Dynamics Analysis Toolkit (AMDAT)*/
/*Vector_Map class*/
/*Written to replace functionality in unordered_map to allow use compilation of AMDAT in currently common gcc compiler versions*/
/*Written by David S. Simmons*/


#ifndef VECTOR_MAP
#define VECTOR_MAP

#include <vector>
#include <utility>
#include <iterator>
#include <stdexcept>

namespace std{
 
template <class Keyclass, class Dataclass>
class Vector_Map
{
    vector <Keyclass> keyvector;
    vector <Dataclass> datavector;
      
  public:
    Vector_Map(){};
    Dataclass& at(const Keyclass& key);
    const Dataclass& at(const Keyclass& key) const;
    void erase(const Keyclass & key); 
    int find (const Keyclass& key);
    int find (const Keyclass& key) const;
    bool insert (const Keyclass& key ,const Dataclass & data); 
    Dataclass& operator[](const Keyclass & key);
    int count(const Keyclass & key) const;
    void clear();
    
};

template <class Keyclass, class Dataclass>
void Vector_Map< Keyclass,  Dataclass>::erase(const Keyclass & key)
{
  int keyii;
  keyii=find(key);
  keyvector.erase(keyvector.begin()+keyii);
  datavector.erase(datavector.begin()+keyii);
}

template <class Keyclass, class Dataclass>
int Vector_Map< Keyclass,  Dataclass>::find(const Keyclass& key)
{
  int keyii;
  for(keyii=0;keyii<keyvector.size();keyii++)
  {
    if(key==keyvector[keyii])
    {
      return keyii;
    }
  }
  return -1;
}


template <class Keyclass, class Dataclass>
int Vector_Map< Keyclass,  Dataclass>::find(const Keyclass& key) const
{
  int keyii;
  for(keyii=0;keyii<keyvector.size();keyii++)
  {
    if(key==keyvector[keyii])
    {
      return keyii;
    }
  }
  return -1;
}


template <class Keyclass, class Dataclass>
Dataclass& Vector_Map<Keyclass, Dataclass>::at (const Keyclass& key)
{
  int keyii;
  keyii=find(key);
  if(keyii==-1)
  {
    throw std::out_of_range("key not found");
  }
  else
  {
    return datavector[keyii];
  }
}


template <class Keyclass, class Dataclass>
const Dataclass& Vector_Map<Keyclass, Dataclass>::at (const Keyclass& key) const
{
  int keyii;
  keyii=find(key);
  if(keyii==-1)
  {
    throw std::out_of_range("key not found");
  }
  else
  {
    return datavector[keyii];
  }
}


template <class Keyclass, class Dataclass>
bool Vector_Map<Keyclass, Dataclass>::insert (const Keyclass& key ,const Dataclass & data)
{
  int keyii;
  keyii=find(key);
  if(keyii==-1)
  {
    keyvector.push_back(key);
    datavector.push_back(data);
    return 1;
  }
  else
  {
    return 0;
  }
}

template <class Keyclass, class Dataclass>
Dataclass& Vector_Map< Keyclass,  Dataclass>::operator[] (const Keyclass& key)
{
  int keyii;
  keyii=find(key);
  if(keyii==-1)
  {
    keyvector.push_back(key);
    datavector.emplace_back();
  }
  else
  {
    return datavector[keyii];
  }
}

template <class Keyclass, class Dataclass>
int Vector_Map< Keyclass, Dataclass>::count(const Keyclass & key) const
{
  int keyii;
  keyii=find(key);
  if(keyii==-1)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}


template <class Keyclass, class Dataclass>
void Vector_Map< Keyclass, Dataclass>::clear() 
{
  keyvector.clear();
  datavector.clear();
}


}

#endif