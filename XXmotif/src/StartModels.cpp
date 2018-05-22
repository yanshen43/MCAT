#include "StartModels.h"
#include "output.h"
#include "utils.h"
#include "ProgressBar.h"
#include "elongationPhase/ExtensionTable.h"

#include <list>
#include <iostream>
#include <string>

/**************************************
 * find overrepresented IUPAC patterns in positive set compared to the negative set
 **************************************/
 void StartModels::findInitialMotifs(MotifContainer& startModels, double pValThreshold){
	//long time1=time(NULL);

	if(Global::aa){
		elongList el;
		MadonaPro::findInitialMotifs(el);
		addElongatedToStartModels(el, startModels);
		Output::printMergedMotifs(startModels, 10);

		Pool_alloc<std::_List_node<Match> >::reset_Pool();
	}else{
		if(!Global::batch_mode) cout << endl;
		cout << "\tELONGATE IUPAC SEED PATTERNS" << endl << endl;

#ifdef SIXMER
		if(Global::GAPS == 0) Global::GAPS = 1;
		else if(Global::GAPS == 1) Global::GAPS = 6;
		else if(Global::GAPS == 2)	Global::GAPS = 21;
		else Global::GAPS = 56;
#else
		if(Global::GAPS == 0) Global::GAPS = 1;
		else if(Global::GAPS == 1) Global::GAPS = 5;
		else if(Global::GAPS == 2)	Global::GAPS = 15;
		else Global::GAPS = 35;
#endif
		if(Global::type == ALL || Global::type == FIVEMERS || Global::type == NOPALINDROME || Global::type == NOTANDEM){
			/* iterate through all allowed gap combinations */
			for (int k = 0; k < Global::GAPS; k++) {
				/* create elongation object with gapIndex k*/
				ElongCandidates ec(k, FIVEMERS, pValThreshold);

				/* fill list of all kmers which should be elongated */
				elongList el;
				ec.getElongationCandidates(el);
				Pool_alloc<std::_List_node<Match> >::reset_Pool();

				addElongatedToStartModels(el, startModels);
			}
			if(!Global::batch_mode) cout << endl;
		}
		if(Global::type == ALL || Global::type == PALINDROME || Global::type == NOTANDEM){
			/* iterate through all palindromic start motifs */
			Global::GAPS = 16;
			for (int k=0; k < Global::GAPS; k++){
				if(k + 6 > Global::maxMatchPositions) continue;
				/* create elongation object with gapIndex k*/
				ElongCandidates ec(k, PALINDROME, pValThreshold);

				/* fill list of all kmers which should be elongated */
				elongList el;
				ec.getElongationCandidates(el);
				Pool_alloc<std::_List_node<Match> >::reset_Pool();

				addElongatedToStartModels(el, startModels);
			}
			if(!Global::batch_mode) cout << endl;
		}
		if(Global::type == ALL || Global::type == TANDEM || Global::type == NOPALINDROME){
			/* iterate through all palindromic start motifs */
			Global::GAPS = 16;
			for (int k=0; k < Global::GAPS; k++){
				if(k + 6 > Global::maxMatchPositions) continue;
				/* create elongation object with gapIndex k*/
				ElongCandidates ec(k, TANDEM, pValThreshold);

				/* fill list of all kmers which should be elongated */
				elongList el;
				ec.getElongationCandidates(el);
				Pool_alloc<std::_List_node<Match> >::reset_Pool();

				addElongatedToStartModels(el, startModels);
			}
		}
		ExtensionTable::freeTable();
		ElongCandidates::freeMemory();

		//cerr << "]" << endl;
		//cerr << endl << "totalCounts: " << startModels.getMotifNb() << endl;
	}

	Pool_alloc<std::_List_node<UngappedKmer> >::reset_Pool();
	Pool_alloc<std::_List_node<std::shared_ptr<Kmer> > >::reset_Pool();

	if(startModels.getMotifNb() == 0) {
		std::cerr << std::endl << "WARNING: no start models, so no motifs found." << std::endl;
	}
    //fprintf(stderr,"\n\n------------time: %ld seconds (%0.2f minutes)---------------\n",
    //       time(NULL)-time1,(float)(time(NULL)-time1)/60.0);

    if( Global::em ){
    	startModels.sort();
    	return; /* no merging so far */
    }
    else if( Global::cv ){
    	startModels.sort();
    	Output::writeBestKmerSites( startModels, 100 );
    	return;
    }else if (Global::noRefinementPhase){
    	startModels.sort();
    	Output::printMergedMotifs(startModels, 10);
    	return;
    }else{
    	/* Calculate new p-values with the branch & bound algorithm and sort all motifs by this p-value */
    	startModels.update_Pval_and_sort( Global::pseudo );
    }

	startModels.removeRedundantMotifs();

	if(Global::maxIterations >= 0){
		startModels.filter(1000, std::numeric_limits<int>::max(), log(1e3));
	}

	Output::printMergedMotifs(startModels, 10);

	if (Global::maxIterations >= 0) {
		startModels.merge(false);
		//startModels.filter(50, LogTable::LOG_i[10]);
	}

	//Output::printTopN(startModels, 20);
	Output::printMergedMotifs(startModels, 10);
	//exit(-1);
}

void StartModels::initStartMotif(MotifContainer& startModels){
	Motif* m;
	if( Global::bindingSiteFile != NULL ){
		m = new Motif( Global::bindingSiteLength );
		m->initMotifWithBindingSites( Global::bindingSiteFile );
	}
	else{
		m = new Motif( PWM_LENGTH );
		if( Global::profFile != NULL ){
			m->InitProfMotif( Global::profFile );
		}
		else if( Global::startMotif != NULL ){
			m->InitStartMotif( Global::startMotif, Global::pseudo );
		}
	}
	if( Global::DEBUG )
		m->setTracked();
	startModels.add( m );
}

void StartModels::addElongatedToStartModels(elongList& el, MotifContainer& startModels){

	/* sort elongated kmers by P-value */
	el.sort(ElongCandidates::compareKmerResults);

	/* remove motifs having exactly the same positive instances */
	//printOverrepPatternsIUPAC(el, 100);
	removeRedundantMotifs(el);
	//printOverrepPatternsIUPAC(el, 100);
	//exit(-1);

	elongList::iterator it = el.begin();
	while(it != el.end()){
		//cerr << "Pval_pos: " << (**it).p_pos << endl;

		/* create new Motif instance with sorted match positions */
		Motif* m = new Motif(*it, PWM_LENGTH);
		if(m->getTotalSites() < (Global::aa ? 2 : 3)){
			delete m;
			it = el.erase(it);
			continue;
		}
		/* if there might be an enriched region, filter to it */
		//if(!Global::aa && Global::usePositionalProbs)
		//		m->filter_region();

		if( Global::em ){
			m->initMotifWithStartPos( Global::multipleOccurrence );
		}
		else{
			if( Global::multipleOccurrence )
				m->updatePWM_OOL_initial( Global::pseudo );
			else
				m->filter_oops_updatePWM(); /* filter to one occurrence per sequence */
		}

		//cout << *m;

		it = el.erase(it);
		startModels.add(m);
	}
	//startModels.clear();
}

int StartModels::removeRedundantMotifs(elongList& el){
	int totCounts = 0;

	elongList::iterator next = el.begin();
	if(next++ != el.end()) totCounts++;

	if(Global::aa){
		for(elongList::iterator it = el.begin(); next != el.end();){
			int identical = 0;
			/* check whether two motifs are completely identical */
			if((*it)->seeds.size() == (*next)->seeds.size()){
				identical = 0;

				MatchContainer& sp1 = (*it)->seeds;
				sp1.sort(Match::cmpMatches);

				MatchContainer& sp2 = (*next)->seeds;
				sp2.sort(Match::cmpMatches);

				MatchContainer::iterator it1 = sp1.begin();
				MatchContainer::iterator it2 = sp2.begin();

				while(*it1 == *it2){
					it1++;
					it2++;
					if(it1 == sp1.end() && it2 == sp2.end()){
						identical = 1;
						break;
					}
				}
			}
			if(identical){
				next = el.erase(next);
			}else{
				totCounts++;
				it++;
				next++;
			}
		}
	}else{
		for(elongList::iterator it = el.begin(); next != el.end();){
			/* check whether two motifs are completely identical */
			if((*next)->getKmer()->bestNucString().compare((*it)->getKmer()->bestNucString()) == 0){
        for (Global::tracked_t::iterator s = Global::trackedMotifs.begin(); s != Global::trackedMotifs.end(); ++s) {
          if((*next)->getKmer()->bestNucString().compare(*s) == 0){
            s = Global::trackedMotifs.erase(s);
            s = Global::trackedMotifs.insert(s, (*it)->getKmer()->bestNucString());
          }
        }
				//cerr << (*next)->getKmer()->bestNucString() << " vs " << (*it)->getKmer()->bestNucString() << endl;
				next = el.erase(next);
			}else{
				totCounts++;
				it++;
				next++;
			}
		}
	}
	return totCounts;
}

void StartModels::printOverrepPatternsIUPAC(elongList& el, int N){
	printf("\n\nTOP %d Motifs (first instance):\n", N);
	int offset = Global::posSet->max_leng-Global::downstream;

	std::iostream::fmtflags flags = cout.flags();
	std::streamsize prec = cout.precision();

	int printCount = 0;
	for(elongList::iterator it = el.begin(); it != el.end() && printCount < N; printCount++, it++){
		Kmer& skr = **it;
		region& r = skr.enrichment;
		if(!Global::aa){
			cout << " " << std::setw(3) << printCount +1 << ": " << skr.getKmer()->toString(0, ElongCandidates::IUPAC_CHARS) << " ";
		}else{
			cout << " " << std::setw(3) << printCount +1 << ": " << skr.getKmer()->toString(0) << " ";
		}
		cout << "posCount: " << std::setw(3) << skr.seeds.size() << "  prob: " << std::scientific << std::setprecision(2) << skr.p_pos;
		if(!Global::aa){
			printf(" my: %3d\tRegion: %3d/%3d", r.max-offset, r.startRegion-offset, r.endRegion-offset);
		}
		cout << "\t\tlog(pVal): " << skr.p_set << "\tE-Value: " << exp(skr.p_set) << endl;
	}	

	cout.flags(flags);
	cout << std::setprecision(static_cast<int>(prec));
}
