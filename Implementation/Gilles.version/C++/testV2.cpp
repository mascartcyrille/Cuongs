#include "misc.hpp"
#include "DataSpike.hpp"
#include "computeb_V2.hpp"
#include "computeG_V2.hpp"

int main(void)
{
  int M=1;
  vec DN0, DN1(1), DN2(1), DN3(1), DN4(1), DN8(2), DN9(3), DN10(2);
  uvec neur1(1, fill::ones);
  uvec neur2(1), neur8(2), neur9(3), neur10(2);
  double Tmin = 0.;
  double Tmax = 1.;
  double delta = 0.1;  
  int K = 1;
  mat G0, G1, G2, G3, G4, G5, G6, G7, G8, G9, G10;
  umat b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
  uvec low, cnt;
  DataSpike DS;
  /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  Tmax =3;
  K = 3;
  cout << "Test 0: One neuron, zero point - (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b0 = computeb_V2(M, K, Tmin, Tmax, DN0, neur1, delta, low, cnt);
  //  G0 = computeG_V2(M, K, Tmin, Tmax, DN0, neur1, delta, low);
  b0 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G0 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b0 <<endl;
  cout << G0 <<endl;
  /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  K = 1;
  Tmax = 1;
  DN1 << .5;
  DS  = DataSpike(DN1, neur1);
  cout << "Test 1: One neuron, one point - spike at 0.5 - (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b1 = computeb_V2(M, K, Tmin, Tmax, DN1, neur1, delta, low, cnt);
  //  G1 = computeG_V2(M, K, Tmin, Tmax, DN1, neur1, delta, low);
  b1 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G1 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  cout << b1 <<endl;
  cout << G1 <<endl;
  /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  DN2 << 0.;
  delta = 1.;
  DS  = DataSpike(DN2, neur1);
  cout << "Test 2: One neuron, one point - spike at 0, (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b2 = computeb_V2(M, K, Tmin, Tmax, DN2, neur1, delta, low, cnt);
  //  G2 = computeG_V2(M, K, Tmin, Tmax, DN2, neur1, delta, low);
  b2 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G2 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  cout << b2 <<endl;
  cout << G2 <<endl;
  /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  DN3 << 1.;
  DS  = DataSpike(DN3, neur1);
  cout << "Test 3: One neuron, one point - spike at 1, (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b3 = computeb_V2(M, K, Tmin, Tmax, DN3, neur1, delta, low, cnt);
  //  G3 = computeG_V2(M, K, Tmin, Tmax, DN3, neur1, delta, low);
  b3 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G3 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b3 <<endl;
  cout << G3 <<endl;
  /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  DN4 << 2.;
  DS  = DataSpike(DN4, neur1);
  cout << "Test 4: One neuron, one point - spike at 2, (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b4 = computeb_V2(M, K, Tmin, Tmax, DN4, neur1, delta, low, cnt);
  //  G4 = computeG_V2(M, K, Tmin, Tmax, DN4, neur1, delta, low);
  b4 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G4 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b4 <<endl;
  cout << G4 <<endl;
    /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  M = 2;
  delta = 0.1;
  neur2 << 1;
  DS  = DataSpike(DN1, neur2);
  cout << "Test 5: Two neurons, one point - spike at .5, (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b5 = computeb_V2(M, K, Tmin, Tmax, DN1, neur2, delta, low, cnt);
  //  G5 = computeG_V2(M, K, Tmin, Tmax, DN1, neur2, delta, low);
  b5 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G5 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  cout << b5 <<endl;
  cout << G5 <<endl;
      /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  M = 1; K = 2;
  delta = 0.2;
  DN8 << 0.5 <<  0.6;
  neur8 = uvec(2, fill::ones);
  DS  = DataSpike(DN8, neur8);
  cout << "Test 8: One neuron, two points at the wrong place (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
//  b8 = computeb_V2(M, K, Tmin, Tmax, DN8, neur8, delta, low, cnt);
//  G8 = computeG_V2(M, K, Tmin, Tmax, DN8, neur8, delta, low);
  b8 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G8 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b8 <<endl;
  cout << G8 <<endl;

        /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  M = 1; K = 2;
  delta = 0.2;
  DN9 << 0.5<< 0.6<< 0.7;
  neur9 = uvec(3, fill::ones);
  DS  = DataSpike(DN9, neur9);
  cout << "Test 9: With three points (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b9 = computeb_V2(M, K, Tmin, Tmax, DN9, neur9, delta, low, cnt);
  //  G9 = computeG_V2(M, K, Tmin, Tmax, DN9, neur9, delta, low);
  b9 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G9 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b9 <<endl;
  cout << G9 <<endl;

          /// -----------------------------------------------------------
  /// -----------------------------------------------------------
  M = 2; K = 1;
  delta = 0.1;
  double x = .2;
  DN10 << 0.5<< 0.5+x;
  neur10 <<1<< 2;
  DS  = DataSpike(DN10, neur10);
  cout << "Test 10: Avec deux neurones et correspondance croisee (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
//  b10 = computeb_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low, cnt);
//  G10 = computeG_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low);
  b10 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G10 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  cout << b10 <<endl;
  cout << G10 <<endl;

  x = 0;
  DN10 << 0.5<< 0.5+x;  
  DS  = DataSpike(DN10, neur10);
  cout << "Test 10: Avec deux neurones et correspondance croisee (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
//  b10 = computeb_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low, cnt);
//  G10 = computeG_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low);
  b10 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G10 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);
  
  cout << b10 <<endl;
  cout << G10 <<endl;

  x = 0.1;
  DN10 << 0.5 << 0.5+x;  
  DS  = DataSpike(DN10, neur10);
  cout << "Test 10: Avec deux neurones et correspondance croisee (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b10 = computeb_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low, cnt);
  //  G10 = computeG_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low);
  b10 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G10 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  
  cout << b10 <<endl;
  cout << G10 <<endl;

  K = 2; x = 0.2; 
  DN10 << 0.5 << 0.5+x;  
  DS  = DataSpike(DN10, neur10);
  cout << "Test 10: Avec deux neurones et correspondance croisee (Tmin, Tmax]=(" << Tmin<<", "<<Tmax<<"]" << endl;
  //  b10 = computeb_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low, cnt);
  //  G10 = computeG_V2(M, K, Tmin, Tmax, DN10, neur10, delta, low);
  b10 = computeb_V2(M, K, Tmin, Tmax, DS, delta, low, cnt);
  G10 = computeG_V2(M, K, Tmin, Tmax, DS, delta, low);

  cout << b10 <<endl;
  cout << G10 <<endl;
}

//int main(void)
//{
//  int M = 1;
//  int K = 5;
//  double delta = 0.2;
//  double Tmin = 0.;
//  double Tmax = 2.;
//  vec DN(6);
//  DN << .1, .4, .45, .5, .6, .66;
//  uvec neur = uvec::Constant(6, 1);
//  mat b;
//  uvec low, cnt;
//  b = computeb_V2(M, K, Tmin, Tmax, DN, neur, delta, low, cnt);
//  cout << b << endl;
//}
