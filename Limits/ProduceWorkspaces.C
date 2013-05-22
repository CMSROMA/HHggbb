#include "HH2ggbbFitter_mgg.cc"

void ProduceWorkspaces(){

  int mass = 300;
  cout << "mass = " << mass << endl; 
  // runfits(mass, true);
  runfits(mass, false);
  
  /*
  mass = 500;
  cout << "mass = " << mass << endl; 
  runfits(mass, true);
  runfits(mass, false);

  mass = 700;
  cout << "mass = " << mass << endl; 
  runfits(mass, true);
  runfits(mass, false);

  mass = 1000;
  cout << "mass = " << mass << endl; 
  runfits(mass, true);
  runfits(mass, false);
  */
}
