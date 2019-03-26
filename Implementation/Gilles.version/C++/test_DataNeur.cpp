#include "DataNeur.hpp"

int main(void){
  for(size_t i=0;i<7;i++)
    cout << pow(10,i) << " " << get_digit(pow(10,i)) << endl;
  cout << 123456 << " " << get_digit(123456)  << endl;
  cout << setprecision(8) << pow(10,6)+1 << " " << get_digit(pow(10,6)+1)  << endl;

   mat A9 = DataNeur("../../data/M=1000/N", 100, 500, 600); 
   A9.save("/tmp/A9.txt", raw_ascii);
   mat A10 = DataNeur("../../data/M=1000/N", 100, 501, 600); 
   A10.save("/tmp/A10.txt", raw_ascii);
}
