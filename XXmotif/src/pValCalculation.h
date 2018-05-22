#ifndef _PVALCALC_H_
#define _PVALCALC_H_

#include "Globals.h"
#include "LogTable.h"
#include "NullModel.h"
#include "refinementPhase/Motif.h"

//#ifdef __SUNPRO_C
//#include <sunmedia_intrin.h>
//#else
//#include <emmintrin.h>
//#include <pmmintrin.h>  // SSE3
//#include <tmmintrin.h>  // SSSE3
//#include <smmintrin.h>  // SSE4
//#endif

class PVal_Calculator{

public:
	static PVal_Calculator &getInstance();

	template<class KGEN> double calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops);
	template<class KGEN> double __calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops);

	double calculatePvalCons(int seq, int pos, int firstMotifColumn, const motif_columns_type& sortedColumns);
	double calculateFinalPval(int len, int delta_0, double p);
	static double getCombinedPval_log(double log_prod_pcons, int motifNb);
	double getWeightedPval(double p1, double p2, double w);
	static double getWeightedPval_log(double p1, double p2, double w);

	double getConsCorrection(int length) { return pConsCorrection[length]; }

	motif_columns_type getSortedColumns(double** pwm, motif_columns_type& columns);
	void initPositionalProbs(int motifLength, int firstMotifColumn, region r, int offset); /* precalculate positional probabilities */
	//void initPositionalProbs_newLength(int motifLength, int offset);

private:
	PVal_Calculator();
	~PVal_Calculator();
	PVal_Calculator(const PVal_Calculator &);             // intentionally undefined
	PVal_Calculator & operator=(const PVal_Calculator &); // intentionally undefined

	int sse_count_mismatches_gaps(const uint8_t* query_profile, const uint8_t* db_sequence,
												   const int dbseq_length, int query_length, int pos);
	double* positionalProb;
	int** tmpArray;
	int** tmpArray2;
	double* sum_i;
	double* pConsCorrection;
	double* pOverrepCorrection;

	uint8_t* result_ali_free;

public:
	double getOverrepCorrection(const int num) { return pOverrepCorrection[num]; }

};

template<class KGEN>
inline double PVal_Calculator::__calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
		int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
		double& pValScore, bool mops){

	double finalPval, pValPos;
  pValScore = kmerHash.get_kmer_probability_all_positions(kmer);

//  fprintf(stderr, "seq: %d, pos: %d, pValue: %f, size: %d\t", seq, pos, pValScore, length);
//	for(int m=0; m<1; m++){
//	    for(int b= 0; b < length; b++)fprintf(stderr, "%c", AlphaChar(kmer[b], Global::A));
//		fprintf(stderr, "\n");
//	}

	if(pValScore == 1 || pValScore == -1) return pValScore;

	pValScore *= pOverrepCorrection[static_cast<int>(cols.size())];

	//if(pValScore < 1e-7)fprintf(stderr, "%4d/%4d\tpValScore: %e * correction %f = %e\n", seq, pos, pValScore, pow(1.1, cols.size()) , pValScore * pow(1.1, cols.size()));
	if(Global::usePositionalProbs && !Global::aa) {
		pValPos = positionalProb[pos]; 						/* position probability from precalculated values */
		if(mops) {
			/* combine Score and Pos pValue: weight pValPos by 0.5
			 * formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
			finalPval = (pValScore*sqrt(pValPos) - pValScore*pValScore*pValPos*0.5)*2;
			//double p = pValScore*pValPos;
			//finalPval = p*(1-log(p));
		}else{
			finalPval = calculateFinalPval(len, delta_0, pValScore*pValPos);
		}
		//if(finalPval < 0.05) {
		//	fprintf(stderr, "len: %d, delta_0: %d, pValScore*pValPos: %e\n", len, delta_0, pValScore*pValPos);
		//	fprintf(stderr, "%d/%d, pValPos: %e, pValScore: %e, finalPval: %e\n", seq, pos, pValPos, pValScore, finalPval);
		//}
	} else {
		if(mops){
			finalPval = pValScore;
		}else{
			if (Global::fixedPosition) {
				finalPval = pValScore / NullModel::getProbability("$", 1);
			} else {
				finalPval = 1-pow(1-pValScore,avgLength-length+1);
			}
		}
//		if(finalPval < 0.001){
//			for(int i=0; i< size; i++) cerr << AlphaChar(kmer[i], Global::A);
//			fprintf(stderr, "\tSeq: %d, Pos: %d, pValScore: %f, finalPval: %f\n", seq, pos, pValScore, finalPval);
//		}
	}

//	if(Global::posSet->max_MultSeq > 1){
//		//for(int i=0; i< size; i++) cerr << AlphaChar(kmer[i], Global::A);
//		//fprintf(stderr, "\tSeq: %d, Pos: %d, pValScore: %f, finalPval: %f\n", seq, pos, pValScore, finalPval);
//
//		double pValCons = calculatePvalCons(seq, pos, cols.front(), cols);
//		pValCons *= pConsCorrection[cols.size()];
//		if(pValCons > 1) pValCons = 1;
//		finalPval = getWeightedPval(finalPval, pValCons, Global::consPvalWeight);
//		//fprintf(stderr, "pCons: %f, finalPval: %f\n", pValCons, finalPval);
//	}

	return finalPval;
}


template<class KGEN>
inline double PVal_Calculator::calculatePval(KGEN& kmerHash, unsigned char* kmer, const motif_columns_type &cols, \
			int seq, int pos, int len, int delta_0, int size, int length, double avgLength,\
			double& pValScore, bool mops){

    memcpy(kmer, Global::posSet->entity[seq]->S[0] + pos, size);
    return __calculatePval(kmerHash, kmer, cols, seq, pos, len, delta_0, size, length, avgLength,
    			pValScore, mops);
}

inline double PVal_Calculator::calculateFinalPval(int len, int delta_0, double p){
	double plateau = pow(1-p*len/delta_0, delta_0);
	double hill = exp(-p*len*sum_i[delta_0]);
	return 1-plateau*hill;
}

inline double PVal_Calculator::getCombinedPval_log(double log_prod_pcons, int motifNb){

	double result = log_prod_pcons;

	double quotient = 1;
	double sum = 1;

	double overflowConst = 1e100;
	const double log_overflowConst = 230.258509299; // log(overflowConst);
	for(int i=1; i<motifNb; i++){
		quotient *= ( (-log_prod_pcons) / i );
		sum += quotient;

		if(sum > overflowConst){
			sum /= overflowConst;
			quotient /= overflowConst;
			result += log_overflowConst;
		}
		//fprintf(stderr, "i: %d, quotient: %e, sum: %e\n", i, quotient, sum);
	}
	//if(sum <= 0) fprintf(stderr, "log_prod_pcons: %f (%.2e), motifNb: %d, result: %.2e, sum: %.2e, final: %.2e\n", log_prod_pcons, exp(log_prod_pcons), motifNb, result, sum, result+log(sum));
	assert(sum > 0);
	//if(motifNb == 9025)
	//	fprintf(stderr, "result: %e + log(sum): %e = %e (%e)\n", result, log(sum), result+log(sum), exp(result+log(sum)));

	return result + log(sum); // log( p * sum from 0 to n-1 { (-ln p)^i / i! } )
}

//inline double PVal_Calculator::getWeightedPval_log(double p1, double p2, double w){
//	p1 = exp(p1);
//	p2 = exp(p2);
//
//	/* formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
//	double pComb = (p1*pow(p2,w) - pow(p1,1/w)*p2*w) / (1-w);
//	return log(pComb);
//}

inline double PVal_Calculator::getWeightedPval_log(double p1, double p2, double w){
	/* formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
	double pComb = p1 + w*p2 + log(1.0/(1-w) * (1-w*pow(exp(p1+w*p2),1.0/w-1)));
	return pComb;
}

inline double PVal_Calculator::getWeightedPval(double p1, double p2, double w){
	/* formula: (p1*p2^w - p1^(1/w)*p2*w) / (1-w) */
	double pComb = (p1*pow(p2,w) - pow(p1,1/w)*p2*w) / (1-w);
	fprintf(stderr, "PV: comb(%g, %g) = %g\n", p1, p2, pComb);
	return pComb;
}

inline int PVal_Calculator::sse_count_mismatches_gaps(const uint8_t* query_profile,
											   const uint8_t* db_sequence,
											   const int dbseq_length,
											   int query_length, int pos
											   //unsigned char* results,
											   )
{
//  pos += query_length;
//  __m128i Mismatches = _mm_set1_epi8(16);
//  __m128i *qji;               // query profile score in row j (for residue x_j)
//  __m128i *query_profile_it = (__m128i *) query_profile;
//
//  int bestRes=16;
//  for (int j=0; j<dbseq_length; ++j) // loop over db sequence positions
//  {
//      // Get address of query scores for row j
//      qji = query_profile_it + db_sequence[j];
//
//      // Shift Mismatches by 1 and shifting in a zero
//      Mismatches = _mm_slli_si128(Mismatches, 1);
//      Mismatches = _mm_adds_epu8(Mismatches, *qji);
//
//      //results[(int)*(((char*) &Mismatches) + query_length)]++;
//      //results[j] = (int)*(((char*) &Mismatches) + query_length);
//      if((int)*(((char*) &Mismatches) + query_length) < bestRes){
//          //if(pos != j)
//    	     bestRes = (int)*(((char*) &Mismatches) + query_length);
//    	  if(bestRes == 0) return 0;
//      }
//      //if(j==pos){
//      //   fprintf(stderr, "pos: %d, aliBased score: %d\n", pos, (int)*(((char*) &Mismatches) + query_length));
//      //}
//  }
//
//  return bestRes;
	return 0;

  //printf("bestRes: %d, bestPos: %d\n", bestRes, bestPos);
  //printf("dbSeq_length: %d\n", dbseq_length);

  /*
  //
  // build minimum with sse2 instructins ... same cpu time
  //
  int iter = dbseq_length / 16 + 1;
  int rounds = iter / 127 + 1;
  for( int round = 0; round < rounds; round++){
          int roundOffset = round*2032;

	  __m128i MismatchesMin = _mm_set1_epi8(16);
	  __m128i bestPosArray = _mm_set1_epi8(0);
	  __m128i MismatchBlock;
	  __m128i Tmp, Tmp2, Pos;
	  __m128i *p = (__m128i*)(results + roundOffset);

	  int iterations = std::min(iter - 127*round , 127);
	  for (int i = 0; i < iterations; ++i)
	  {
		Pos = _mm_set1_epi8(i);
		MismatchBlock = _mm_load_si128(p++);
		Tmp = _mm_cmplt_epi8(MismatchBlock, MismatchesMin);
		MismatchesMin = _mm_min_epu8(MismatchBlock, MismatchesMin);
		Tmp2 = _mm_and_si128(Tmp, Pos);
		bestPosArray = _mm_max_epu8(bestPosArray, Tmp2);

		//for (int b = 0; b<16; b++) // loop over bytes in XMM registers
		//    printf("%d %d   Pos: %3i   MismatchBlock: %3i   Tmp: %3i   MismatchesMin: %3i   Tmp2: %3i  bestPosArray: %3i\n",
		//        i, b, (int)*(((char*) &Pos) + b)*16+b + 127*round*16, (int)*(((char*) &MismatchBlock) + b),
		//	  (int)*(((char*) &Tmp) + b), (int)*(((char*) &MismatchesMin) + b), (int)*(((char*) &Tmp2) + b), (int)*(((char*) &bestPosArray) + b));

	  }
	  int bestResRound = 16;
          int bestPosRound = 0;
	  for(int i = 0; i < 16; ++i){
		//printf("i: %d, mismatch: %d, pos: %d\n", i, (int)*(((char*) &MismatchesMin) + i), (int)*(((char*) &bestPosArray) + i) );
		int mismatch = (int)*(((char*) &MismatchesMin) + i);
		if(mismatch < bestResRound){
			bestResRound = mismatch;
			bestPosRound = roundOffset + (int)*(((char*) &bestPosArray) + i) * 16 + i;
		}
	  }
	  //printf("round %d, bestRes: %d, bestPos: %d\n", round, bestResRound, bestPosRound);

  	  if(bestResRound < bestRes){
		bestRes = bestResRound;
 		bestPos = bestPosRound;
	  }
  }
  */
}

#endif
