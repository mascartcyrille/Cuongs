#include <iostream>
#include <cmath>
#include <cstring>
#include <armadillo>
#include "misc.hpp"

using namespace::std;
using namespace::arma;

// On crée un DataNeur
mat DataNeur(string mystr, int M, size_t debneur, size_t nmax)
///
/// A DataNeur is simply a rectangular matrix containing spikes of size M-by-nmax where M is the number of neurons. nmax is the maximum number of spikes  +1
/// The 1st column of the DataNeur contains the number of spikes for each neuron
/// The other columns contain spike values. nan if no spike

{
  mat D(M, nmax, fill::zeros);
  size_t N = get_digit(M +debneur-1); // there are M neurons from debneur to (M+debneur-1)
  size_t N0 = get_digit(debneur);
  int irow = 0;
  int ideb , ifin;
  mat Dstart(1,1);
  string s1;
  ifstream ifs;
  string buf;
  bool okfile;
  //  cout << "N=" << N << endl;
  //  cout << "mystr " << mystr << " M" << M << " debneur " << debneur << " nmax " << nmax << endl;
  //  cout << debneur << " " <<"N0=" << N0 << endl;
  for(size_t ip=0;ip<=N-N0;ip++)
    {
      if(ip==0) {
	ideb = debneur;
	ifin = std::min<size_t>(debneur+M-1, pow(10, N0)-1);} // si (502, 569)  et  M=68 ?
      else {
	ideb = pow(10, ip);
	ifin = std::min<size_t>( M+debneur-1, pow(10, ip+1)-1);
      }
      //      cout << N-N0-ip << endl;
      s1 = string(N-N0-ip, '0'); //(N-ip) times '0'
      //      cout << s1 << endl;

      // s = mystr + s1 + i + ".txt"
      for(size_t i=ideb;i<=ifin;i++)
	{
	  //	  mat tmp;
	  vec tmp;
	  string s;
	  s.append(mystr);

	  if(s1.size())
	    {
	      //	      cout << "DEBUG !!!!! " << N-ip-1 << " " << s1 << endl;
	      s.append(s1);
	    }
	  s.append(to_string(i));
	  s.append(".txt");
	  //	  cout << i << " size " << strlen(s.c_str()) << " " <<s << endl;
	  //	  cout << s.c_str() << endl;
	  // test if the file is empty or not - okfile = true => file is not empty, otherwise it is  empty
	  ifs.open(s.c_str());
	  ifs >> buf;
	  if(!ifs.eof())
	    okfile = true;
	  else okfile = false;
	  ifs.close();
	  if(okfile) {
	  tmp.load(s, raw_ascii);
	  Dstart(0,0) = tmp.size();
	  //	  cout << "irow " << irow << " " <<tmp.size() << endl;
	  D.submat(irow, 0, irow, tmp.size()) = join_rows(Dstart, tmp.st());}
	  else
	    D.submat(irow, 0, irow, 0) = 0; // no spike
	  irow++;
	}
    }
  //   cout << "Fin loop " << endl;
  nmax = max(D.col(0))+1;
  //  cout << "nmax vaut " << nmax << endl;
  D = D.submat(0, 0, irow-1, nmax);
  //  D.save("/tmp/D.txt", raw_ascii);
  return D;
}


// On crée un DataNeur
mat DataNeur(string mycsvfile, int M)
///
/// A DataNeur is simply a rectangular matrix containing spikes of size M-by-nmax where M is the number of neurons. nmax is the maximum number of spikes  +1
/// The 1st column of the DataNeur contains the number of spikes for each neuron
/// The other columns contain spike values. nan if no spike

{
  mat A;
  A.load(mycsvfile.c_str(), csv_ascii);
}
