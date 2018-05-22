#include "conservationScores.h"
#include "../ProgressBar.h"

void printIndex(uint32_t index){
	int A = index % 16;	index /= 16;
	int C = index % 16;	index /= 16;
	int G = index % 16;	index /= 16;
	int T = index % 16;
	fprintf(stderr, "A: %d, C: %d, G: %d, T: %d\n", A, C, G, T);
}

void setConservationProbs(){
	ss_type set = Global::posSet;
	if(Global::negSet != NULL) set = Global::negSet;

	float*** consProbs = allocateArray(set);
	//float*** aliFree = allocateArray(set);

	int** tmpArray = (int**)malloc((MAX_CONSERVATION_LENGTH+1)*sizeof(int*));
	for(int i=0; i<=MAX_CONSERVATION_LENGTH; i++) tmpArray[i] = (int*)calloc(set->max_MultSeq, sizeof(int));

	int** tmpArray2 = (int**)malloc((MAX_CONSERVATION_LENGTH+1)*sizeof(int*));
	for(int i=0; i<=MAX_CONSERVATION_LENGTH; i++) tmpArray2[i] = (int*)calloc(set->max_MultSeq, sizeof(int));

	char result_ali_free[16];
	char result_ali_free_rev[16];
	unsigned char q[80];

	//clock_t start = clock();

	std::cout << std::endl;

	for(int i=1; i<= set->nent; i++){
		if(!Global::batch_mode){
			if(i%50 == 0 || i == set->nent)
				std::cout << "\rcalibrate alignment pValues:   "<< progressBar(i, set->nent, 40) << std::flush;
		}

		int seq_length = set->entity[i]->n;

		int multSeqNb = set->entity[i]->mseq;
		unsigned char** alignment = set->entity[i]->S;

		for(int j=1; j<= seq_length; j++){
			/* max length for start postition j */
			int max_conservation_length = std::min(seq_length-j+1, MAX_CONSERVATION_LENGTH);

			/* create query profile */
			for (int k=0; k < max_conservation_length*5; k++) q[k] = 1;
			for (int k=0; k < max_conservation_length; k++){
				//fprintf(stderr, "%c", AlphaChar(alignment[0][j+k], Global::A));
				q[alignment[0][j+k]*16 + k] = 0;
			}
			//fprintf(stderr, "\tseq: %d, pos: %d\n", i, j);

			/****
			 * count alignment free matches
			 ****/
			if(Global::useAliFree){
				int region = 50;

				int ali_free_start = std::max(1, j-region);
				int ali_free_end = std::min(seq_length, j+region);
				int ali_free_length = ali_free_end-ali_free_start+1;
				int ali_free_pos = j-region > 0 ? region : j;

				for(int m=1; m<multSeqNb; m++){
					uint8_t* db_seq  = &alignment[m][ali_free_start];
					//for(int b=0; b<ali_free_length; b++)fprintf(stderr, "%c", AlphaChar(db_seq[b], Global::A));
					//fprintf(stderr, "\n");

					sse_count_mismatches(q, db_seq, ali_free_length, result_ali_free, ali_free_pos);

					if(Global::revcomp){
						int ali_free_start_rev = std::max(1, seq_length-j-region);
						int ali_free_end_rev = std::min(seq_length, seq_length-j+region);
						int ali_free_length_rev = ali_free_end_rev - ali_free_start_rev + 1;
						int ali_free_pos_rev = seq_length-j-region > 0 ? region : seq_length-j;
						db_seq  = &alignment[m][ali_free_start_rev];
						//for(int b=0; b<ali_free_length_rev; b++)fprintf(stderr, "%c", AlphaChar(db_seq[b], Global::A));
						//fprintf(stderr, "\n");

						sse_count_mismatches(q, db_seq, ali_free_length_rev, result_ali_free_rev, ali_free_pos_rev);
						for(int pos=0; pos < 16; pos++){
							//fprintf(stderr, "m: %d, length: %d, +: %d, -: %d\n", m, pos+1, result_ali_free[pos], result_ali_free_rev[pos]);
							result_ali_free[pos] = std::min(result_ali_free[pos], result_ali_free_rev[pos]);
						}
					}
					//fprintf(stderr, "\n");
					//  exit(-1);

					for(int k=1; k<=max_conservation_length; k++){
						tmpArray2[k][m] = result_ali_free[k-1] + tmpArray2[k][m-1];
					}
				}
			}

			/****
			 * count alignment based matches
			 ****/
			uint32_t index = 0;
			int realMultSeqNb = multSeqNb;
			for(int k=1; k<=max_conservation_length; k++){
				int base = alignment[0][j+k-1];
				if(base == 0) break;

				index += mutIndex[base-1];
				for(int m=1; m<realMultSeqNb; m++){
					int multBase = alignment[m][j+k-1];
					// don't count mutations if part of sequence is not alignable
					if(multBase == 0){
						realMultSeqNb = m;
						break;
					}

					int mutations = tmpArray[k-1][m] - tmpArray[k-1][m-1] + tmpArray[k][m-1];
					if(multBase != base) mutations++;
					tmpArray[k][m] = mutations;

					/* store results in array */
					mutations = std::min(std::min((int)(k*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), tmpArray[k][m]);
					consProbs[index][m][mutations]++;
					//consProbs[k][m][mutations]++;

					//mutations = std::min(std::min((int)(k*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), tmpArray2[k][m]);
					//aliFree[index][m][mutations]++;
				}

				/* if not all sequences are alignable, treat them as completely mutated */
				/*for(int m = multSeqNb; m < set->max_MultSeq; m++){
					//fprintf(stderr, "\n seq: %d, info: %s, multSeq: %d, maxSeq: %d\n", i, Global::negSet->entity[i]->info[0], (int)Global::negSet->entity[i]->mseq, (int)Global::negSet->max_MultSeq);

					tmpArray[k][m] = tmpArray[k-1][m] - tmpArray[k-1][m-1] + tmpArray[k][m-1] + 1;
					int mutations = std::min(std::min((int)(k*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), tmpArray[k][m]);
					consProbs[index][m][mutations]++;

					tmpArray2[k][m] = tmpArray2[k-1][m] - tmpArray2[k-1][m-1] + tmpArray2[k][m-1] + 1;
					mutations = std::min(std::min((int)(k*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS), tmpArray2[k][m]);
					aliFree[index][m][mutations]++;
				}*/
			}
			/*
			fprintf(stderr, "\nalignment-based\n");
			for(int m=0; m<multSeqNb; m++){
				for(int k=1; k<=max_conservation_length; k++){
					if(m==0){
						int base = alignment[0][j+k-1];
						fprintf(stderr, "%c\t", AlphaChar(base, Global::A));
					}else{	fprintf(stderr, "%d\t", tmpArray[k][m]); }
				} fprintf(stderr, "\n");
			}
			fprintf(stderr, "\ntmpArray2 after alignment-free\n");
			for(int m=0; m<set->max_MultSeq; m++){
				for(int k=1; k<=max_conservation_length; k++){
					if(m==0){
						int base = alignment[0][j+k-1];
						fprintf(stderr, "%c\t", AlphaChar(base, Global::A));
					}else{ fprintf(stderr, "%d\t", tmpArray2[k][m]); }
				} fprintf(stderr, "\n");
			}
			exit(-1);
			*/
		}
	}

	// printf("\ntime:%ld ticks, min: %f\n",clock()- start, (clock()-start) / (CLOCKS_PER_SEC * 60.0));

	for(int i=0; i<=MAX_CONSERVATION_LENGTH; i++){ free(tmpArray[i]); free(tmpArray2[i]); }
	free(tmpArray);	free(tmpArray2);


	float alpha = 1;
	/* create cumulative list and normalize counts */
 	/*for(int length=1; length<=MAX_CONSERVATION_LENGTH; length++){
		for(int m=1; m<set->max_MultSeq; m++){
			int max_cons_mutations = std::min((int)(length*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS);
			for(int k=1; k<=max_cons_mutations; k++){
				consProbs[length][m][k] += consProbs[length][m][k-1];
			}
			float norm_ali_based = consProbs[length][m][max_cons_mutations];

			for(int k=0; k<=max_cons_mutations; k++){
				consProbs[length][m][k] /= norm_ali_based;
			}
		}
	}*/
	
 	for(int A=0; A<=MAX_CONSERVATION_LENGTH; A++){
		for(int C=0; C<=MAX_CONSERVATION_LENGTH-A; C++){
			for(int G=0; G<=MAX_CONSERVATION_LENGTH-A-C; G++){
				for(int T=0; T<=MAX_CONSERVATION_LENGTH-A-C-G; T++){
					int length = A+C+G+T;

					uint32_t index = A*mutIndex[0] + C*mutIndex[1] + G*mutIndex[2] + T*mutIndex[3];
					length = std::min(length, MAX_CONSERVATION_LENGTH);

					for(int m=1; m<set->max_MultSeq; m++){
						int max_cons_mutations = std::min((int)(length*m*MAX_MUTATED_FRACTION), MAX_CONS_MUTATIONS);

						//printIndex(index); fprintf(stderr, ", m: %d, %f + %d*%e = %f\n",
						//		m, consProbs[index][m][0], numberOfPseudo[length], pVal, consProbs[index][m][0] + numberOfPseudo[length]*pVal);
						consProbs[index][m][0] += alpha;
						//aliFree[index][m][0]   += alpha;
						for(int k=1; k<=max_cons_mutations; k++){
							consProbs[index][m][k]+=consProbs[index][m][k-1];
							//aliFree[index][m][k]+=aliFree[index][m][k-1];
						}
						float norm_ali_based = consProbs[index][m][max_cons_mutations];
						//float norm_ali_free =  aliFree[index][m][max_cons_mutations];
						for(int k=0; k<=max_cons_mutations; k++){
							consProbs[index][m][k] /= norm_ali_based;
							//aliFree[index][m][k] /= norm_ali_free;
						}
					}
				}
			}
		}
	}
	
	Global::conservationProbs = consProbs;
	//Global::alignmentFreeProbs = aliFree;

	std::cout << std::endl;

	/*
	uint32_t index = 4*mutIndex[0] + 4*mutIndex[1]  + 4*mutIndex[2] + 4*mutIndex[3];
	std::cerr << std::endl << std::endl << "index: " << index << "\t";
	printIndex(index);
	for(int i=1; i<set->max_MultSeq; i++){
		std::cerr << "mult Seq: " << i << "\t";
		for(int j=0; j<=MIN(int, 16*i, MAX_CONS_MUTATIONS); j++){
			fprintf(stderr, "%d(%.2f%%),", j, consProbs[index][i][j]*100 );
			if(j==17){
				fprintf(stderr, "..."); break;
			}
		}
		std::cerr << "\n";
	}

	fprintf(stderr, "\n\n");
	for(int i=1; i<set->max_MultSeq; i++){
		std::cerr << "mult Seq: " << i << "\t";
		for(int j=0; j<=MIN(int, 16*i, MAX_CONS_MUTATIONS); j++){
			fprintf(stderr, "%d(%.2f%%),", j, aliFree[index][i][j]*100 );
			if(j==17){
				fprintf(stderr, "..."); break;
			}
		}
		std::cerr << "\n";
	}*/

//	exit(-1);
}
