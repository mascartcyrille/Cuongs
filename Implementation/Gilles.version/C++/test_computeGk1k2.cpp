#include "misc.hpp"
#include "DataSpike.hpp"
#include "computeGk1k2_V2.hpp"
#include "computeG_V2.hpp"
#include "computeb_V2.hpp"

int main(void){
  vec spike({0.1, 0.4, 0.45, 0.5, 0.6, 0.66});
  cout << "spike = " << spike << endl;
  double delta = 0.2;
  vec spike2({0.02, 0.04, 0.15, 0.46, 0.57, 0.6, 0.76, 0.77, 0.95}); 
  cout << "spike2 = " << spike2 << endl;

  uvec ones6(6, fill::ones);
  vec v1 = vec({0.1,0.4, 0.45, 0.5, 0.6, 0.66});
  DataSpike DS(v1, ones6);
  uvec low, cnt;
  umat b = computeb_V2(1, 5, 0, 1., DS, 0.2, low, cnt);
  mat Gv2;
  Gv2 =  computeG_V2(1, 5, 0, 1., DS, 0.2, low);
  //  cout << " Gv2" << endl;
  //  cout << Gv2 << endl;
  //  cout << " avant computeGk1k2_V2 " << endl;
  mat Gk;
  Gk = computeGk1k2_V2(1, 5, 0, 1., DS, 0.2, low);
  //  cout << " Gk avec computeGk1k2_V2 " << endl;
  cout << Gk << endl;
}
