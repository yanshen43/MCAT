#include "hoNullModel.h"

float   hoNullModel::_alpha	 = 2;
int		hoNullModel::_asize	 = 4;
float   hoNullModel::_beta	 = 0;
int*	hoNullModel::_coeffs = NULL;
float*	hoNullModel::_conds	 = NULL;
float*	hoNullModel::_counts = NULL;
float	hoNullModel::_factor = 8; /* 4 * 2 as _asize and _alpha default to 4 and 2 resp. */
int		hoNullModel::_fields = 156; /* 5^0 + 5^1 + 5^2 + 5^3 (1 + mono-/di-/trinucleotides) as _order defaults to 2 */
float*	hoNullModel::_freqs	 = NULL;
bool	hoNullModel::_gaps	 = true;
int		hoNullModel::_order	 = 2;
float*	hoNullModel::_probs	 = NULL;

int		hoNullModel::_bits	  = 3;
float*	hoNullModel::_rconds  = NULL;
int		hoNullModel::_rfields = 512;
float*	hoNullModel::_rprobs  = NULL;

void hoNullModel::init( ss_type sequences, float alpha, float beta, int order, bool gaps, float* freqs ){

	if( alpha < 0 ){
		fprintf( stderr, "Negative pseudocounts factor (alpha): %f\n", alpha );
		exit(0);
	}
	_alpha = alpha;

	if( beta < 0 ){
		fprintf( stderr, "Negative counts offset (beta): %f\n", beta );
		exit(0);
	}
	_beta = beta;

	if( order < 0 ){
		fprintf( stderr, "Negative order: %d\n", order );
		exit(0);
	}
	_order = order;

	_gaps = gaps;

	_freqs = freqs;

	init( sequences );
}

void hoNullModel::init( ss_type sequences ){

	int i, k;
	int l, L;
	int n, N;

	int base, order;
	float ncounts;

	int* offsets;
	unsigned char* s;

	uint8_t *kmer = new uint8_t[_order+1];
	uint8_t *new_kmer = new uint8_t[_order+1];

	_asize = nAlpha( Global::A );
	_factor = static_cast<float>(_asize) * _alpha;

	base = _asize + _gaps;

  // _fields[k] == sum_{i=0}^k base^i
  _fields = 1 + base;
  for (int i = 2, powBase = base; i < _order + 2; i++) {
    powBase *= base;
    _fields += powBase;
  }

	_counts = ( float* )calloc( _fields, sizeof( float ) );

	_conds = ( float* )calloc( _fields, sizeof( float ) );

	_probs = ( float* )calloc( _fields, sizeof( float ) );

  // _coeffs[k] == pow( base, k );
  _coeffs = ( int* )calloc( _order+2, sizeof( int ) ); /* coeffs */
  for (int k = 0, powBase = 1; k < _order + 2; k++) {
    _coeffs[k] = powBase;
    powBase *= base;
  }

	offsets = ( int* )calloc( _order+2, sizeof( int ) ); /* offsets */
	for( k=0; k < _order+2; k++ )
		offsets[k] = k ? offsets[k-1]+_coeffs[k] : _coeffs[k];

	/* Restructuring parameters */

	_bits = ( int )ceil( log( _asize+1 )/log(2) );

	_rfields = ( int )pow( 2, ( _order+1 )*_bits );

	_rconds = ( float* )calloc( _rfields, sizeof( float ) );
	_rprobs = ( float* )calloc( _rfields, sizeof( float ) );

	/* Calculate counts */

	N = sequences->nent;
	for( n=1; n <= N; n++ ){ /* across sequences */

		L = sequences->entity[n]->n;

		if( _order >= L ){
			fprintf( stderr, "Null model order >= (length of) background sequence no. %d\n", n );
			exit(-1);
		}

		s = sequences->entity[n]->S[0];

		for( l=1; l <= L; l++ ){ /* across positions */

			for( k=0; k <= _order && l+k <= L && s[l+k]; k++ ) /* s[l+k]? = char <n>? */
				_counts[ sub2ind( s, l, k ) ]++;
		}
	}

	/* Substract counts offset ( beta ) */

	if( _beta > 0 )
		for( i=1; i < _fields; i++ )
			_counts[i] = std::max( 0.0f, ( float )( _counts[i]-_beta ) );

	/* Calculate total counts for 0th order */

	ncounts = 0;
	for( i=1; i <= _asize; i++ )
		ncounts += _counts[i];

	/* Calculate probs */

	order = 0; /* 1-mers */

	if( order <= _order ){
		if( _freqs != NULL ){
			for( n=1; n <= _asize; n++ ){
				_conds[n] = _probs[n] = ( _counts[n] + _factor*_freqs[n] ) / ( ncounts + _factor );
			}
		}
		else{
			_freqs = ( float* )calloc( _asize+1, sizeof( float ) );
			for( n=1; n <= _asize; n++ ){ /* no pseudocounts in 0th order */
				_conds[n] = _probs[n] = _freqs[n] = _counts[n] / ncounts;
			}
		}
	}

	/* 1-mers++ */

	for( order++; order <= _order; order++ )
		calculateKmerProbs( kmer, 0, order-1 );

	/* Calculate probs (for gapped kmers) */

	if( _gaps ){

		order = 0; /* 1-mers */
		if( order <= _order ){
			_conds[_asize+1] = _probs[_asize+1] = 1;
		}

		/* 1-mers++ */
		for( order++; order <= _order; order++ )
			calculateGappedKmerProbs( kmer, 0, order, -1 );
	}

	/* Printouts */

	if( Global::verbose ){
		printf( "___________\n"
				"           |\n"
				"NULL MODEL |\n"
				"___________|\n\n" );
		printHoNullModel( _counts, _fields-1, offsets );
		printHoNullModel( _conds, _fields-1, offsets );
		printHoNullModel( _probs, _fields-1, offsets );
	}

	/* Copy (conditional) probabilities to Holgers index structure */

	for( order=_order; order >= 0; order-- )
		copyProbs( kmer, new_kmer, 0, order );

	/* Frees */

	free( offsets );
  delete[] kmer;
  delete[] new_kmer;
}

void hoNullModel::save(){

	std::stringstream str;

	str.str( "" );
	str << Global::outputDirectory << "/" << Global::shortFileName << ".freqsBg";
	FILE* fptrFreqs = fopen( str.str().c_str(), "w" );

	str.str( "" );
	str << Global::outputDirectory << "/" << Global::shortFileName << ".condsBg";
	FILE* fptrConds = fopen( str.str().c_str(), "w" );

	str.str( "" );
	str << Global::outputDirectory << "/" << Global::shortFileName << ".probsBg";
	FILE* fptrProbs = fopen( str.str().c_str(), "w" );

	int i, k;

	int* offsets;
	offsets = ( int* )calloc( _order+2, sizeof( int ) );
	for( k=0; k < _order+2; k++ )
		offsets[k] = k ? offsets[k-1]+_coeffs[k] : _coeffs[k];

	for( i=1; i <= nAlpha( Global::A ); ++i )
		fprintf( fptrFreqs, "%f ", _freqs[i] );
	fprintf( fptrFreqs, "\n" );
	fclose( fptrFreqs );

	for( k=1, i=1; i < _fields; ++i ){
		if( i==offsets[k] ){

			fprintf( fptrConds, "\n" );
			fprintf( fptrProbs, "\n" );
			++k;
		}
		fprintf( fptrConds, "%f ", _conds[i] );
		fprintf( fptrProbs, "%f ", _probs[i] );
	}
	fprintf( fptrConds, "\n\n" );
	fprintf( fptrProbs, "\n\n" );

	free( offsets );
}

/* Call: calculateGappedKmerProbs( kmer[order+1], 0, order, -1 )*/
void hoNullModel::calculateGappedKmerProbs( unsigned char* kmer, int pos, int order, int lastGap ){

	if( pos > order ){
		calculateGappedKmerProbs( kmer, order, lastGap );
	}
	else{
		unsigned char k;

		for( k=1; k <= _asize; k++ ){

			kmer[pos] = k;
			calculateGappedKmerProbs( kmer, pos+1, order, lastGap );
		}
		kmer[pos] = k;
		calculateGappedKmerProbs( kmer, pos+1, order, pos );
	}
}

void hoNullModel::calculateGappedKmerProbs( unsigned char* kmer, int order, int lastGap ){

	int i, ii, k;
	unsigned char N, X;

	N = static_cast<unsigned char>(_asize+1);
	i = sub2ind( kmer, order );

	if( lastGap > -1 )
		_conds[i] = _probs[i] / _probs[sub2ind(kmer,order-1)];

	for( k=lastGap+1; k <= order; k++ ){

		X = kmer[k]; kmer[k] = N;
		ii = sub2ind( kmer, order );

		_probs[ii] += _probs[i];

		kmer[k] = X;
	}
}

/* Call: calculateKmerProbs( kmer[order+1], 0, order-1 )*/
void hoNullModel::calculateKmerProbs( unsigned char* kmer, int pos, int order_minus_1 ){

	for( unsigned char k=1; k <= _asize; k++ ){

		kmer[pos] = k;

		if( pos == order_minus_1 )
			calculateKmerProbs( kmer, order_minus_1 );
		else
			calculateKmerProbs( kmer, pos+1, order_minus_1 );
	}
}

void hoNullModel::calculateKmerProbs( unsigned char* kmer, int order_minus_1 ){

	int i, ii, p;
	int order;
	float conds_ii, S;

	int *index = new int[_asize+1];
	int *iindex = new int[_asize+1];

	S = 0;
	order = order_minus_1 + 1;

	for( unsigned char k=1; k <= _asize; k++ ){
		kmer[order] = k;

		ii = sub2ind( kmer, order );
		iindex[k] = ii;

		i = sub2ind( kmer, order-1 );
		index[k] = i;

		p = sub2ind( kmer+1, order-1 );

		conds_ii = ( _counts[ii] + _factor*_conds[p] ) / ( _counts[i] + _factor );
		S += conds_ii;

		_conds[ii] = conds_ii;
	}
	for( unsigned char k=1; k <= _asize; k++ ){

		i = index[k];
		ii = iindex[k];

		_conds[ii] /= S;
		_probs[ii] = _conds[ii] * _probs[i];
	}

	delete[] index;
	delete[] iindex;
}

void hoNullModel::copyProbs( unsigned char* kmer, unsigned char* new_kmer, int pos, int order ){

	if( pos > order ){
		copyProbs( kmer, new_kmer, order );
	}
	else{
		unsigned char k;

		for( k=1; k <= _asize; k++ ){

			kmer[pos] = k;
			new_kmer[pos] = k;

			copyProbs( kmer, new_kmer, pos+1, order );
		}
		kmer[pos] = k;
		new_kmer[pos] = 0;

		copyProbs( kmer, new_kmer, pos+1, order );
	}
}

void hoNullModel::copyProbs( unsigned char* kmer, unsigned char* new_kmer, int order ){

	int i = sub2ind( kmer, order );
	int new_i = sub2ind2( new_kmer, order );

	_rprobs[new_i] = _probs[i];
	_rconds[new_i] = _conds[i];
}
