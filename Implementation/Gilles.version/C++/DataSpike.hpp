#ifndef __DATASPIKE_HPP__
#define __DATASPIKE_HPP__
#include <armadillo>
#include "DataNeur.hpp"
using namespace::arma;
using namespace::std;
class DataSpike{
public:
  vec _T;
  uvec _neur;
  DataSpike(vec & T, uvec & neur): _T(T), _neur(neur) {}; 
  //  DataSpike(const vec T, const uvec neur): _T(T), _neur(neur) {}; 
  DataSpike(): _T(vec()), _neur(uvec()) {}; 
  //  DataSpike(ifstream ifs);
  //  DataSpike(string* fn);
  mat toDataNeur(int nmax=600) const {
    uvec ind = unique(_neur);
    mat DN(ind.size(), nmax);
    //    cout << "ind = " << endl;
    //    cout << ind << endl;
    int imax= -1;
    for(int i=0;i<ind.size();i++)
      {
	//	cout << "in loop toDN" << endl;
	//	cout << i << endl;
	//	cout << _neur << endl;
	uvec qi = find(_neur==(i+1));
	//cout << qi << endl;
	//cout << qi.size() << endl;
	DN(i,0) = qi.size();
	imax = std::max(imax, int(qi.size()));
	for(int j=0;j<qi.size();j++)
	  DN(i,j+1) = _T(qi(j));
      }
    return DN.submat(0, 0, ind.size()-1, imax);
  };

void  FromDataNeur(const mat & DN){
    size_t M = DN.n_rows;
    size_t nmax = DN.n_cols;
    vec flat(M*nmax, fill::ones);
    uvec tn(M*nmax, fill::ones);
    int deb = 0;
    for(int i=0;i<M;i++)
      {
	if(DN(i,0)>0.)
	  {
	    //	cout << DN.n_rows << " " << i << " " << DN.n_cols << " " << DN(i,0) << " " << deb<< " " << deb+int(DN(i,0)-1) << flush <<endl;
		//		cout << flat.subvec(deb, deb+int(DN(i,0)-1)) << " separation " <<  DN(i, span(1, int(DN(i,0)))) << endl; 
	flat.subvec(deb, deb+int(DN(i,0)-1)) = DN(i, span(1, int(DN(i,0)))).st(); 
        tn.subvec(deb, deb+int(DN(i,0)-1)) = (i+1)*uvec(int(DN(i,0)), fill::ones);
	deb += int(DN(i,0));
	  }
      }
    uvec  indb = sort_index(flat.subvec(0, deb-1));
    _T = flat(indb);
    _neur = tn(indb);
  };

  //  DataSpike(const DataSpike & DS) :_T(DS._T), _neur(DS._neur) {};

  DataSpike(string mystr, int M, size_t debneur, size_t nmax)
{
  mat DN = DataNeur(mystr, M, debneur, nmax);
  //  cout << DN << endl;
  FromDataNeur(DN);
};

  DataSpike(mat DN)
{
  FromDataNeur(DN);
};
};

#endif
