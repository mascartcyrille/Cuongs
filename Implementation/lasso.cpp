#include <iostream>
#include <vector>
#include <cmath>

/* Parameters of the computation */
const long double	_spike [] = { 0.1, 0.4, 0.45, 0.5, 0.6, 0.66 },	// The spike time stamps, in increasing order
					_del = 0.2,													// The size of the bins
					_T_max = 2,												// The maximum time
					_gamma = 3.0,												// An arbitrary number
					_epsilon = 0.1;											// An arbitrary precision bound
const long int		_N = 6,														// Number of spikes
					_K = floor( _T_max / _del );								// Number of bins

/* Arrays and matrices to be computed for the LASSO */
long int	**	_A;
long double *	_a,		// Array[1..]
			*	_b,
			*	_d,
			**	_G;

/**
 * Counts pairs of spikes in bins
 * @param  spike [description]
 * @param  N     [description]
 * @param  del   [description]
 * @param  k     [description]
 * @param  A     [description]
 * @return       [description]
 */
long int count_pair( const long double * spike, const long int N, const long double del, const long int k, long int ** A );
/**
 * [computeG description]
 * @param G     [description]
 * @param A     [description]
 * @param del   [description]
 * @param spike [description]
 * @param N     [description]
 * @param K     [description]
 * @param T_max [description]
 */
void computeG( long double ** G, long int ** A, const long double del, const long double * spike, const long int N, const long double K, const long double T_max );
/**
 * [computed description]
 * @param  A     [description]
 * @param  del   [description]
 * @param  spike [description]
 * @param  N     [description]
 * @param  T_max [description]
 * @param  k     [description]
 * @param  gamma [description]
 * @return       [description]
 */
long double computed( long int ** A, const long double del, const long double * spike, const long int N, const long double T_max, long int k, long double gamma );
/**
 * [count_conse description]
 * @param  v [description]
 * @return   [description]
 */
long double count_conse( std::vector<std::pair<long int, long int>>& v );
/**
 * [count_spike description]
 * @param  spike [description]
 * @param  del   [description]
 * @param  k     [description]
 * @param  t     [description]
 * @return       [description]
 */
long int count_spike( const long double * spike, long int N, long double del, long int k, long double t );
/**
 * Applies the LASSO reconstruction method to a series of spikes in order to retrieve the Hawkes point process parameters that created the spikes.
 * @param  argc nothing
 * @param  argv nothing
 * @return      0
 */
int main( int argc, char const **argv ) {
	/* Initialization of the A counting matrix */
	_A = new long int*[ _N ];
	for( long int i = 0; i < _N; ++i ) {
		_A[ i ] = new long int[ _N ];
		for( long int j = 0; j < _N; ++j ) {
			_A[ i ][ j ] = (long int) -2;
		}
	}

	/* Initialization of the a, b, d and G matrices */
	_a = new long double[ _K + 1 ];
	_b = new long double[ _K + 1 ];
	_d = new long double[ _K + 1 ];
	_G = new long double*[ _K + 1 ];
	for( long int i = 0; i <= _K; ++i ) {
		_a[ i ] = _b [ i ] = _d[ i ] = -1;
		_G[ i ] = new long double[ _K + 1 ];
		for( long int j = 0; j <= _K; ++j ) {
			_G[ i ][ j ] = 0;
		}
	}

	/* Compute b */
	_b[ 0 ] = _N;
	for( long int i = 0; i < _K; ++i ) {
		_b[ i + 1 ] = count_pair( _spike, _N, _del, i, _A );
	}

	/* Compute G */
	computeG( _G, _A, _del, _spike, _N, _K, _T_max );
	for( long int j = 0; j <= _K; ++j ) {
		for( long int i = 0; i < j; ++i ) {
			_G[ i ][ j ] = _G[ j ][ i ];
		}
	}

	/* Compute d */
	_d[ 0 ] = std::sqrt( 2 * _gamma * std::log( _T_max ) * _N ) + _gamma * std::log( _T_max ) / 3;
	for( int k = 0; k < _K; ++k ) {
		_d[ k + 1 ] = computed( _A, _del, _spike, _N, _T_max, k, _gamma );
	}

	/* Printing the vectors and matrices */
	std::cout.precision(5);
	std::cout << "B: (";
	for( long int i = 0; i < _K; ++i ) {
		std::cout << std::fixed << _b[ i ] << ", ";
	}
	std::cout << std::fixed << _b[ _K ] << ")" << std::endl;

	std::cout << "A: (";
	for( long int i = 0; i < _N - 1; ++i ) {
		std::cout << "(";
		for( long int j = 0; j < _N - 1; ++j ) {
			std::cout << std::fixed << _A[ i ][ j ] << ",\t";
		}
		std::cout << std::fixed << _A[ i ][ _N - 1 ] << ")" << std::endl << "\t";
	}
	std::cout << "(";
	for( long int j = 0; j < _N - 1; ++j ) {
		std::cout << std::fixed << _A[ _N - 1 ][ j ] << ",\t";
	}
	std::cout << std::fixed << _A[ _N - 1 ][ _N - 1 ] << "))" << std::endl;

	std::cout << "G: (";
	for( long int i = 0; i < _K; ++i ) {
		std::cout << "(";
		for( long int j = 0; j < _K; ++j ) {
			std::cout << std::fixed << _G[ i ][ j ] << ",\t";
		}
		std::cout << std::fixed << _G[ i ][ _K ] << ")" << std::endl << "\t";
	}
	std::cout << "(";
	for( long int j = 0; j < _K; ++j ) {
		std::cout << std::fixed << _G[ _K ][ j ] << ",\t";
	}
	std::cout << std::fixed << _G[ _K ][ _K ] << "))" << std::endl;

	std::cout << "D: (";
	for( long int i = 0; i < _K; ++i ) {
		std::cout << std::fixed << _d[ i ] << ", ";
	}
	std::cout << std::fixed << _d[ _K ] << ")" << std::endl;

	/* Clean up the memory */
	delete [] _a;
	delete [] _b;
	delete [] _d;
	for( long int i = 0; i < _N; ++i ) {
		delete [] _A[ i ];
	}
	delete [] _A;
	for( long int i = 0; i < _K; ++i ) {
		delete [] _G[ i ];
	}
	delete [] _G;

	/* Tell the system everything went well */
	return 0;
}

long int count_pair( const long double * spike, const long int N, const long double del, const long int k, long int ** A ) {
	long double	T_low = 0.0,
				T_up = 0.0;
	long int	count = 0,
				j_start = 0,
				i = 0;

	while( j_start < N && spike[ j_start ] <= 0 ) ++j_start;

	for( long int j = j_start; j < N; ++j ) {
		T_low = spike[ j ] - (k + 1) * del;
		T_up = spike[ j ] - k * del;
		i = j - 1;
		while( i >= 0 && i < N ) {
			if( T_low  <= spike[ i ] && spike[ i ] < T_up ) {
				count = count + 1;
				A[ i ][ j ] = k;
			} else if( T_low > spike[ i ] ) {
				break;
			}
			--i;
		}
	}

	return count;
}

void computeG( long double ** G, long int ** A, const long double del, const long double * spike, const long int N, const long double K, const long double T_max ) {
	long double	rescase = 0.0,
				rescase1 = 0.0,
				rescase2 = 0.0,
				rescase3 = 0.0,
				length_inter = 0.0;
	std::vector<std::pair<long int, long int>> v;
	for( long int l = 0; l < K; ++l ) {
		for( long int k = l + 1; k < K; ++k ) {
			rescase1 = rescase2 = 0.0;
			v.clear();
			for( long int i = 0; i < N; ++i ) {
				for( long int j = 0; j < N; ++j ) {
					if( A[ i ][ j ] == k - l ) {
						v.push_back( std::make_pair( i, j ) );
					}
				}
			}
			for( long int idex = 0; idex < v.size(); ++idex ) {
				length_inter = std::min( spike[ v[ idex ].first ] + (k + 1) * del, T_max ) - std::max( spike[ v[ idex ].second ] + l * del, (long double) 0.0 );
				if( length_inter > 0 ) {
					rescase1+=length_inter;
				}
			}
			v.clear();

			for( long int i = 0; i < N; ++i ) {
				for( long int j = 0; j < N; ++j ) {
					if( A[ i ][ j ] == k - l - 1 ) {
						v.push_back( std::make_pair( i, j ) );
					}
				}
			}
			for( long int idex = 0; idex < v.size(); ++idex ) {
				length_inter = std::min( spike[ v[ idex ].second ] + (l + 1) * del, T_max ) - std::max( spike[ v[ idex ].first ] + k * del, (long double) 0.0 );
				if( length_inter > 0 ) {
					rescase2 += length_inter;
				}
			}

			G[ k + 1 ][ l + 1 ] = rescase1 + rescase2;
		}
	}
	for( long int k = 0; k < K; ++k ) {
		rescase1 = rescase3 = 0.0;
		v.clear();
		for( long int i = 0; i < N; ++i ) {
			for( long int j = 0; j < N; ++j ) {
				if( A[ i ][ j ] == 0 ) {
					v.push_back( std::make_pair( i, j ) );
				}
			}
		}
		for( long int idex = 0; idex < v.size(); ++idex ) {
			length_inter = std::min( spike[ v[ idex ].first ] + ( k + 1 ) * del, T_max ) - std::max( spike[ v[ idex ].second ] + k * del, (long double) 0.0 );
			if( length_inter > 0.0 ) {
				rescase1 += length_inter;
			}
		}
		v.clear();

		for( long int i = 0; i < N; ++i ) {
			length_inter = std::min( spike[ i ] + (k + 1) * del, T_max ) - std::max( spike[ i ] + k * del, (long double) 0.0 );
			if( length_inter > 0.0 ) {
				rescase3 += length_inter;
			}
		}
		G[ k + 1 ][ k + 1 ] = rescase1 * 2 + rescase3;
	}
	_G[ 0 ][ 0 ] = T_max;
	for( long int k = 0; k < K; ++k ) {
		rescase = 0.0;
		for( long int i = 0; i < N; ++i ) {
			length_inter = std::min( spike[ i ] + ( k + 1 ) * del, T_max ) - std::max( spike[ i ] + k * del, (long double) 0.0 );
			if(length_inter > 0.0) {
				rescase += length_inter;
			}
		}
		G[ k + 1 ][ 0 ] = rescase;
	}
}

long double computed( long int ** A, const long double del, const long double * spike, const long int N, const long double T_max, long int k, long double gamma ) {
	long int z [N] = {0},
				max = 0.0;
	std::vector<std::pair<long int, long int>> v;
	for( int i = 0; i < N; ++i ) {
		if( ((spike[ i ] + (k + 1) * del) <= T_max) && ((spike[ i ] + (k + 1) * del) >= 0) ) {
			z[ i ] = count_spike( spike, N, del, k, spike[ i ] + (k + 1) * del );
		}
		if( ((spike[ i ] + k * del) >= 0) && ((spike[ i ] + k * del) < T_max) ) {
			z[ i ] = std::max( z[ i ], (long int) 1 );
		}
	}
	for( long int i = 0; i < N; ++i ) {
		for( long int j = 0; j < N; ++j ) {
			if( A[ i ][ j ] == k ) {
				v.push_back( std::make_pair( i, j ) );
			}
		}
	}

	for( int i = 0; i < N; ++i ) {
		if( z[ i ] > max ) {
			max = z[ i ];
		}
	}

	return gamma * std::log( T_max ) / 3 * max + std::sqrt( 2 * gamma * std::log( T_max ) * count_conse( v ) );
}

long double count_conse( std::vector<std::pair<long int, long int>>& v ) {
	long int	count = 0,
				conse = 1,
				compare,
				n = v.size();

	if( v.empty() ) {
		return 0;
	} else {
		compare = v[ 0 ].second;
		v.push_back( std::make_pair( (long int) -1, (long int) -1 ) );
		for( long int i = 0; i < n; ++i ) {
			if( v[ i + 1 ].second == compare ) {
				++conse;
			} else {
				count += conse * (conse - 1);
				conse = 1;
				compare = v[ i + 1 ].second;
			}
			++i;
		}

		return n + count;
	}
}

long int count_spike( const long double * spike, long int N, long double del, long int k, long double t ) {
	long int count = 0;
	for( int i = 0; i < N; ++i ) {
		if( (spike[ i ] >= (t - (k + 1) * del)) && (spike[ i ] < (t - k * del)) ) {
			++count;
		}
	}
	return count;
}