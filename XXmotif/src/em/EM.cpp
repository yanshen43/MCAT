#include "EM.h"
#include "../output.h"
#include "../refinementPhase/Motif.h"

#include "hoNullModel.h"

#include <cstdlib>
#include <ctime>

void EM::cv( int folds ){

	/* Bisher keine Reinitialisierung */

	int N = Global::posSet->nent;

	int i, k, n;
	int TESTnr;
	int TRAININGnr;

	int maxTESTnr;
	int maxTRAININGnr;

	if( N % folds == 0 ){

		maxTESTnr = N / folds;
		maxTRAININGnr = N-maxTESTnr;
	}
	else{
		maxTESTnr = static_cast<int>(ceil( static_cast<float>(N) / static_cast<float>(folds) ));
		maxTRAININGnr = N-maxTESTnr+1;
	}

	int* TEST = ( int* )malloc( ( maxTESTnr+1 )*sizeof( int ) );
	int* TRAINING = ( int* )malloc( ( maxTRAININGnr+1 )*sizeof( int ) );

	for( i=1; i <= folds; i++ ){

		k = i;
		TESTnr = 0;
		TRAININGnr = 0;

		for( n=1; n <= N; n++ ){

			if( n == k ){
				TEST[++TESTnr] = n;
				k += folds;
			} else{
				TRAINING[++TRAININGnr] = n;
			}
		}

		go( TRAININGnr, TRAINING );
	}

	free( TEST );
	free( TRAINING );
}

void EM::go(){

	int N = Global::posSet->nent;

	int* sequenceNr = ( int* )malloc( ( N+1 )*sizeof( int ) );

	for( int n=1; n <= N; n++ )
		sequenceNr[n] = n;

	go( N, sequenceNr );

	free( sequenceNr );
}

void EM::go( int N, int sequenceNr[] ){

	if( Global::verbose )
		printf( "___\n   |\nEM |\n___|\n\n" );

	int b, i, j, k, n, nr;
	int L, LW1, W, w1;
	int Wbgo, Wo;

	int round = 0;

	bool iterate;
	double l, ll, llo, normFactor;

	float **counts, **conds;
	double **probs;

	int* offsets;
	unsigned char* s;

	float q = Global::q;
	float qIq = q / ( 1-q );
	float Iqq = ( 1-q ) / q;

//	const float eps = 0.0001f; // 10^-4
//	const float eps = 0.00001f; // 10^-5
	const float eps = 0.000001f; // 10^-6

	ss_type sequences = Global::posSet;

	if( Global::verbose )
		printf( "N: %d\n", N );

	int order;
	int bgo = Global::modelOrderBg;

	float* bgConds = hoNullModel::getConds();
	float* bgProbs = hoNullModel::getProbs();

	double** bg = ( double** )malloc( (N+1)*sizeof( double* ) );
	double** posterior = ( double** )malloc( (N+1)*sizeof( double* ) );

	Motif* motif;
	list<Motif*>::const_iterator iter;

	int Wmin = INT_MAX, Wmax = INT_MIN;
	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

		motif = *iter;
		W = motif->getMotifLength();

		if( W < Wmin )
			Wmin = W;

		if( W > Wmax )
			Wmax = W;
	}

	for( n=1; n <= N; n++ ){ /* across sequences */

		nr = sequenceNr[n];

		L = sequences->entity[nr]->n;
		s = sequences->entity[nr]->S[0];

		for( k=1; k <= L; k++ ){ /* Check for N's in sequence */
			if( !( s[k] ) ){
				fprintf( stderr, "N's in sequence %d!\n", nr );
				exit(-1);
			}
		}

		if( Wmax > L ){ /* Check for sequences shorter than longest motif(s) */
			fprintf( stderr, "Motif length (%d) > sequence (no. %d) length (%d)!\n", Wmax, nr, L );
			exit(-1);
		}

		bg[n] = ( double* )calloc( L-Wmin+2, sizeof( double ) );
		posterior[n] = ( double* )calloc( L-Wmin+2, sizeof( double ) );
	}

	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

		if( Global::verbose )
			printf( "___________\n           |\nNEXT MODEL |\n___________|\n\n" );

		motif = *iter;

		motif->setFactor( Global::eta );

		W = motif->getMotifLength();
		order = motif->getOrder();
		b = motif->getFirstMotifColumn();

		if( motif->isLog() )
			motif->exp1merProbs();

		conds = motif->getConds();
		probs = motif->getPWM();

		counts = motif->getCounts();

		w1 = W-1;

		Wo = W-order;
		Wbgo = W-bgo;

		llo = 0; /* old ll */
		iterate = true;

		/* Calculate null model probabilities */

		if( w1 > bgo ){ /* use (conditional) probabilities */

			for( n=1; n <= N; n++ ){ /* across sequences */

				nr = sequenceNr[n];

				L = sequences->entity[nr]->n;
				s = sequences->entity[nr]->S[0];

				LW1 = L-W+1;

				for( k=1; k <= LW1; k++ ){ /* across positions */

					l = bgProbs[ hoNullModel::sub2ind( s, k, bgo ) ];

					for( i=k+1; (i-k) < Wbgo; i++ )
						l *= bgConds[ hoNullModel::sub2ind( s, i, bgo ) ];

					bg[n][k] = l; /* Possible to speed up reusing (conditional) probabilities */
				}
			}
		}
		else{ /* use probabilities */

			for( n=1; n <= N; n++ ){ /* across sequences */

				nr = sequenceNr[n];

				L = sequences->entity[nr]->n;
				s = sequences->entity[nr]->S[0];

				LW1 = L-W+1;

				for( k=1; k <= LW1; k++ ){ /* across positions */

					bg[n][k] = bgProbs[ hoNullModel::sub2ind( s, k, w1 ) ];
				}
			}
		}

		do{
			/* E step */

			if( Global::verbose ){
				offsets = motif->getOffsets();
				printf( "_______\n"
						"       |\n"
						"COUNTS |\n"
						"_______|\n\n" );
				printHoModel( counts, b, b+W-1, offsets[order+1]-1, offsets );
			}

			ll = 0;

			for( n=1; n <= N; n++ ){ /* across sequences */

				nr = sequenceNr[n];

				L = sequences->entity[nr]->n;
				s = sequences->entity[nr]->S[0];

				LW1 = L-W+1;

				normFactor = 0;

				/* k=1...L-W+1 */

				for( k=1; k <= LW1; k++ ){ /* across positions */

					l = probs[ b+order ][ motif->sub2ind( s, k, order ) ];
					for( i=k+1; (i-k) < Wo; i++ )
						l *= conds[ b+(i-k)+order ][ motif->sub2ind( s, i, order ) ];

					l /= bg[n][k];

//					posterior[n][k] = l * q / LW1;
					posterior[n][k] = l * qIq / LW1; /* after introducing eta */

					normFactor += posterior[n][k];

					ll += l;
				}

				/* k=0 */

//				posterior[n][0] = 1-q;
				posterior[n][0] = 1; /* after introducing eta */

				normFactor += posterior[n][0];

				ll++;

				for( k=0; k <= LW1; k++ ) /* posterior normalizations */
					posterior[n][k] /= normFactor;
			}

			ll = fast_log( static_cast<float>(ll) );

			if( Global::verbose )
				printf( "LL: %g\n", ll );

			/* Check convergence */

			if( fabs( ll-llo ) < eps ){

				iterate = false;
			}
			else{

				llo = ll;

				/* M step */

				motif->resetCounts();

				for( n=1; n <= N; n++ ){ /* across sequences */

					nr = sequenceNr[n];

					L = sequences->entity[nr]->n;
					s = sequences->entity[nr]->S[0];

					LW1 = L-W+1;

					for( k=1; k <= LW1; k++ ){ /* across positions in sequence */

						for( j=0; j < W; j++ ){ /* across positions in motif */

							for( i=0; i <= order && j+i < W; i++ ){ /* across orders */

								counts[ b+j+i ][ motif->sub2ind( s, k+j, i ) ] += static_cast<float>(posterior[n][k]);
							}
						}
					}

//					Johannes' Beschleunigung
//
//					unsigned char* seqcur = s;
//					for( k=1; k <= LW1; k++ ){ /* across positions in sequence */
//
//						for( j=0; j < W; j++ ){ /* across positions in motif */
//
//							int index=0; // index for single nucleotide
//							int i1 = std::min(j,order);
//							for( i=0; i <= i1; i++ ){ /* across orders */
//								index += (*(s[k+j-i]) << (_bits*i))
//								counts[ b+j][ index ] += posterior[n][k];
//							}
//						}
//						seqcur += ;  //advance on position in sequence
//					}
				}

				motif->multipleCounts( Iqq ); /* after introducing eta */

				motif->substractCountsOffset();
				motif->calculateHoProbs();

				round++;
			}

		} while( iterate );

		if( Global::verbose ){
			offsets = motif->getOffsets();
			printf( "_______\n"
					"       |\n"
					"COUNTS |\n"
					"_______|\n\n" );
			printHoModel( counts, b, b+W-1, offsets[order+1]-1, offsets );
			printf( "__________________________\n"
					"                          |\n"
					"CONDITIONAL PROBABILITIES |\n"
					"__________________________|\n\n" );
			printHoModel( conds, b, b+W-1, offsets[order+1]-1, offsets );
			printf( "____________________\n"
					"                    |\n"
					"TOTAL PROBABILITIES |\n"
					"____________________|\n\n" );
			printHoModel( probs, b, b+W-1, offsets[order+1]-1, offsets );
			printf( "(after %d iterations)\n", round );
		}

		motif->resetFactor(); // _factor = _asize * _alpha
	}

	if( Global::save ){

		// Prints PWMs (Holger-like)
		//
		// For further processing use create_instances_from_pwm.pl

		for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

			motif = *iter;
			motif->log1merProbs();
		}
		Output::print_pwm_folder( _models );

		// Prints HO-PWMs
		//
		// For further processing use plotOrderContr.R

		for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

			motif = *iter;
			if( motif->isLog() )
				motif->exp1merProbs();
		}
		saveHoModel( _models );
	}

	/* Frees */

	for( n=1; n <= N; n++ ){
		free( bg[n] );
		free( posterior[n] );
	}
	free( bg );
	free( posterior );
}

void EM::score(){

	if( Global::verbose )
		printf( "________\n        |\nSCORING |\n________|\n\n" );

	if( !( Global::cv ) || _models.getMotifs().size() > 1 ){
		fprintf( stderr, "What the hell are you doing here?" );
		exit(-1);
	}

	std::stringstream str;
	str << Global::outputDirectory << "/" << Global::shortTestFileName << "-" << Global::modelOrder << "-" << Global::modelOrderBg << ".logOdds";
	FILE* fptr = fopen( str.str().c_str(), "w" );

	int b, i, k, n;
	int L, LW1, W, W1;
	int Wbgo, Wo;

	double f, l;

	float **conds;
	double **probs;

	unsigned char* s;

	ss_type sequences = Global::testSet;
	int N = sequences->nent;

	if( Global::verbose )
		printf( "N: %d\n", N );

	int order;
	int bgo = Global::modelOrderBg;

	float* bgConds = hoNullModel::getConds();
	float* bgProbs = hoNullModel::getProbs();

	Motif* motif;
	list<Motif*>::const_iterator iter;

	int Wmin = INT_MAX, Wmax = INT_MIN;
	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

		motif = *iter;
		W = motif->getMotifLength();

		if( W < Wmin )
			Wmin = W;

		if( W > Wmax )
			Wmax = W;
	}

	for( n=1; n <= N; n++ ){ /* across sequences */

		L = sequences->entity[n]->n;
		s = sequences->entity[n]->S[0];

		for( k=1; k <= L; k++ ){ /* Check for N's in sequence */
			if( !( s[k] ) ){
				fprintf( stderr, "N's in sequence %d!\n", n );
				exit(-1);
			}
		}

		if( Wmax > L ){ /* Check for sequences shorter than longest motif(s) */
			fprintf( stderr, "Motif length (%d) > sequence (no. %d) length (%d)!\n", Wmax, n, L );
			exit(-1);
		}
	}

	for( iter=_models.getMotifs().begin(); iter != _models.getMotifs().end(); iter++ ){

		motif = *iter;
		W = motif->getMotifLength();
		order = motif->getOrder();
		b = motif->getFirstMotifColumn();

		if( motif->isLog() )
			motif->exp1merProbs();

		conds = motif->getConds();
		probs = motif->getPWM();

		W1 = W-1;
		Wo = W-order;
		Wbgo = W-bgo;

		/* calculate log odds */

		for( n=1; n <= N; n++ ){ /* across sequences */

			L = sequences->entity[n]->n;
			s = sequences->entity[n]->S[0];

			LW1 = L-W+1;

			for( k=1; k <= LW1; k++ ){ /* across positions */

				/* calculate background probabilities */

				if( W1 > bgo ){ /* speed up */

					f = bgProbs[ hoNullModel::sub2ind( s, k, bgo ) ];

					for( i=k+1; (i-k) < Wbgo; i++ )
						f *= bgConds[ hoNullModel::sub2ind( s, i, bgo ) ];
				}
				else{

					f = bgProbs[ hoNullModel::sub2ind( s, k, W1 ) ];
				}

				/* calculate foreground probabilities */

				l = probs[ b+order ][ motif->sub2ind( s, k, order ) ];

				for( i=k+1; (i-k) < Wo; i++ )
					l *= conds[ b+(i-k)+order ][ motif->sub2ind( s, i, order ) ];

				/* write log odds */

				fprintf( fptr, "%f ", fast_log( static_cast<float>(l/f) ) );
			}

			fprintf( fptr, "\n" );
		}
	}

	fclose( fptr );
}
