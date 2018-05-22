#ifndef HOUTILS_H
#define HOUTILS_H

#include "../refinementPhase/MotifContainer.h"

using std::cout;
using std::endl;

template< typename T >
void printHoModel( const T a, const int fr, const int to, const int fields, const int* breaks );

template< typename T >
void printHoNullModel( const T a, const int fields, const int* breaks );

void saveHoModel( const MotifContainer &motifs );

template< typename T >
inline void printHoModel( const T a, const int fr, const int to, const int fields, const int* breaks ){

	int i, k, pos;

	for( pos=fr; pos <= to; ++pos ){
		for( k=1, i=1; i <= fields; ++i ){
			if( i==breaks[k] ){

				cout << endl; ++k;
			}
			cout << a[pos][i] << " ";
		}
		cout << endl << endl;
	}
}

template< typename T >
inline void printHoNullModel( const T a, const int fields, const int* breaks ){

	int i, k;

	for( k=1, i=1; i <= fields; ++i ){
		if( i==breaks[k] ){

			cout << endl; ++k;
		}
		cout << a[i] << " ";
	}
	cout << endl << endl;
}

inline void saveHoModel( MotifContainer &models ){

	std::stringstream str;
	FILE *fptrFreqs, *fptrConds, *fptrProbs;

	float* freqs;
	float** conds;
	double** probs;

	int fields;
	int* offsets;

	Motif* motif;
	list<Motif*>::const_iterator iter;

	bool number;
	if( models.getMotifs().size() > 1 )
		number = true;
	else
		number = false;

	int i, k, pos;
	int identifier = 1;
	for( iter=models.getMotifs().begin(); iter != models.getMotifs().end(); iter++, identifier++ ){

		str.str( "" );
		str << Global::outputDirectory << "/" << Global::shortFileName;
		if( number )
			str << "-" << identifier;
		str << ".freqs";
		fptrFreqs = fopen( str.str().c_str(), "w" );

		str.str( "" );
		str << Global::outputDirectory << "/" << Global::shortFileName;
		if( number )
			str << "-" << identifier;
		str << ".conds";
		fptrConds = fopen( str.str().c_str(), "w" );

		str.str( "" );
		str << Global::outputDirectory << "/" << Global::shortFileName;
		if( number )
			str << "-" << identifier;
		str << ".probs";
		fptrProbs = fopen( str.str().c_str(), "w" );

		motif = *iter;
		freqs = motif->getFreqs();
		conds = motif->getConds();
		probs = motif->getPWM();

		offsets = motif->getOffsets();
		fields = offsets[motif->getOrder()+1]-1;

		for( i=1; i <= nAlpha( Global::A ); ++i )
			fprintf( fptrFreqs, "%f ", freqs[i] );
		fprintf( fptrFreqs, "\n" );
		fclose( fptrFreqs );

		for( pos=motif->getFirstMotifColumn(); pos <= motif->getLastMotifColumn(); ++pos ){
			for( k=1, i=1; i <= fields; ++i ){
				if( i==offsets[k] ){

					fprintf( fptrConds, "\n" );
					fprintf( fptrProbs, "\n" );
					++k;
				}
				fprintf( fptrConds, "%f ", conds[pos][i] );
				fprintf( fptrProbs, "%f ", probs[pos][i] );
			}
			fprintf( fptrConds, "\n\n" );
			fprintf( fptrProbs, "\n\n" );
		}

		fclose( fptrConds );
		fclose( fptrProbs );
	}
}

#endif /* HOUTILS_H */
