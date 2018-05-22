#ifndef __CONSERVATION_H__
#define __CONSERVATION_H__

#include "../Globals.h"
#include <stdint.h>

/********************************* PUBLIC *******************************/
void setConservationProbs();

/********************************* PRIVATE *******************************/
float*** allocateArray(ss_type set);
void printIndex(uint32_t index);

inline float*** allocateArray(ss_type set){
	float*** probs = (float***)calloc(((int)pow(2,16)+1), sizeof(float**));
	/* initialize Array that can index at most A + C + G + T = MAX_CONSERVATION_LENGTH kmers
	 * e.g. A=10, C=0, G=0, T=0 */

	for(int A=0; A<=MAX_CONSERVATION_LENGTH; A++){
		for(int C=0; C<=MAX_CONSERVATION_LENGTH; C++){
			for(int G=0; G<=MAX_CONSERVATION_LENGTH; G++){
				for(int T=0; T<=MAX_CONSERVATION_LENGTH; T++){
					int length = A+C+G+T;
					if(length>MAX_CONSERVATION_LENGTH) continue;
					uint32_t index = A*mutIndex[0] + C*mutIndex[1] + G*mutIndex[2] + T*mutIndex[3];
					probs[index] = (float**)malloc(set->max_MultSeq * sizeof(float*));
					length = std::min(length, MAX_CONSERVATION_LENGTH);
					for(int m=1; m<set->max_MultSeq; m++){
						probs[index][m] = (float*)calloc(std::min((int)(length*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS) + 1, sizeof(float));
					}
				}
			}
		}
	}
	return probs;
}

inline int sse_count_mismatches(const uint8_t* query_profile,
								const uint8_t* db_sequence,
								const int      dbseq_length,
								char* results, int pos)
{
//int j; // position in db sequence (0,..,dbseq_length-1)
//__m128i Mismatches = _mm_set1_epi8(16);
//__m128i MismatchesMin = _mm_set1_epi8(16);
//__m128i *qji;               // query profile score in row j (for residue x_j)
//__m128i *query_profile_it = (__m128i *) query_profile;
//
////fprintf(stderr, "pos: %d\n", pos);
//for (j=0; j<pos; ++j) // loop over db sequence positions
//{
//     // Get address of query scores for row j
//      qji = query_profile_it + db_sequence[j];
//      // Shift Mismatches by 1 and shifting in a zero
//      Mismatches = _mm_slli_si128(Mismatches, 1);
//      Mismatches = _mm_adds_epu8(Mismatches, *qji);
//      MismatchesMin = _mm_min_epu8(MismatchesMin, Mismatches);
//
//      //fprintf(stderr, "j: %d\t", j);
//     // for(int b=0; b<16; ++b)
//    //	  fprintf(stderr, "%d: %d\t", b+1, *(((char*) &MismatchesMin) + b));
//     // fprintf(stderr, "\n");
//
//
//}
//
////Mismatches = _mm_set1_epi8(16);
////fprintf(stderr, "==================\n");
//
//for (j=pos; j<dbseq_length; ++j) // loop over db sequence positions
//{
//     // Get address of query scores for row j
//     qji = query_profile_it + db_sequence[j];
//     // Shift Mismatches by 1 and shifting in a zero
//     Mismatches = _mm_slli_si128(Mismatches, 1);
//     Mismatches = _mm_adds_epu8(Mismatches, *qji);
//     MismatchesMin = _mm_min_epu8(MismatchesMin, Mismatches);
//
//     //fprintf(stderr, "j: %d\t", j);
//     //for(int b=0; b<16; ++b)
//     //   fprintf(stderr, "%d: %d\t", b+1, *(((char*) &MismatchesMin) + b));
//     //fprintf(stderr, "\n");
//     //printf("\nDB-pos: %i\n", j+1);
//     //for (int b = 0; b<16; b++) // loop over bytes in XMM registers
//     //  printf(" %3i   %3i   %3i \n",(int)*(((char*) qji) + b), (int)*(((char*) &Mismatches) + b), (int)*(((char*) &MismatchesMin) + b));
//}
//
//for (int b = 0; b<16; b++) // loop over bytes in XMM registers
//    results[b] = *(((char*) &MismatchesMin) + b);

return 0;
}

#endif
