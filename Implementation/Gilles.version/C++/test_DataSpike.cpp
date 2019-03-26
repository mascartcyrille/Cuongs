#include "DataSpike.hpp"
int main(void)
{
  vec T = {0.1, 0.15, 0.34};
  mat DN = {{2, 0.1, 0.27}, {1, 0.14, 0}, {1, 0.46, 0}}; 
  uvec n = {1, 2, 1};
  uvec n1 = n.subvec(0,1);
  DataSpike DS(T, n);
  cout << DS._T << " " << DS._neur << endl;
  DataSpike DS1(T, n1);
  cout << DS1._T << " " << DS1._neur << endl;  
  cout << "fin" << endl;
  DataSpike DS2(DN);
  cout << DS2._T << " " << DS2._neur << endl;    
  cout << DS.toDataNeur() << endl;
  cout << DS2.toDataNeur() << endl;

  DataSpike DS3("../data/M=1000/N", 100, 500, 600);
  DS3._T.save("/tmp/W.txt", raw_ascii);
  DS3._neur.save("/tmp/W1.txt", raw_ascii);

}
