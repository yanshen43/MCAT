#ifndef MOTIF
#define MOTIF

#include "Sorted_Sites.h"
#include "../Globals.h"
#include "../LogTable.h"
#include "../ThresholdChecker.h"
#include "../elongationPhase/Match.h"
#include "../elongationPhase/Kmer.h"
#include "../elongationPhase/StartPos.h"
#include "../nucleotides/motifRegion.h"
#include "../memoryPool/pool_alloc.h"
#include <iomanip>
#include <string>
#include <list>
#include <map>

using std::ostream;
using std::string;
using std::endl;
using std::list;
using std::map;
using std::pair;
using std::cout;
using std::cerr;

class Motif{
public:
	/* constructors */
	Motif(int length);
	Motif(const std::shared_ptr<Kmer>& skr, int length);
	Motif(const Motif& other);

	void InitProfMotif(char* profFile);
	void InitStartMotif(char* startMotif, double pseudo);

	/* destructor */
	~Motif();

	Motif* getPalindrome();

	/* getter methods */
	int getLength() const{ return _length; }
	int getUsedLength() const { return static_cast<int>(_motifColumns.size()); }
	int getMotifLength() const { return _motifColumns.back() - _motifColumns.front() + 1; }
	double getWorstScore(){return _worstScore;}

	const motif_columns_type& getMotifColumns() const { return _motifColumns; }
	motif_columns_type& getMotifColumns() { return _motifColumns; }

	int getFirstMotifColumn() const { return _motifColumns.front(); }
	int getLastMotifColumn() const { return _motifColumns.back(); }
	int getTotalSites() const { return static_cast<int>(_sites.size()); }
	int getTotalBindingSites() const { return static_cast<int>(_bindingSites.size()); }
	int32_t getPosSetSize() const { return _posSetSize; }
	void setPosSetSize(int32_t size) { _posSetSize = size; }

	double** getPWM() const { return _pwm; }

	StartPosContainer& getStartPosList(){ return _sites; }
	MatchContainer& getBindingSites(){ return _bindingSites; }
	const StartPosContainer& getStartPosList() const { return _sites; }
	double getPval() const { return _pVal; }
	region getEnrichment(){ return _enrichment; }
	bool isTracked() { return _tracked; }

	bool isConverged() const { return _converged; }
	void setConverged() { _converged = true; }
	void resetConverged() { _converged = false; }

	int* getCountStartPos_counts() { return _countStartPos.counts; }
	int& getCountStartPos_size() { return _countStartPos.size; }

	const string getIUPACString(char* IUPAC_long = NULL) const;
	void offset(int offset);

	/* setter methods */
	void setPval(double pVal){ _pVal = pVal; }
	void setEnrichment(region other){ _enrichment = other; }
	void setMotifColumns(motif_columns_type& columns){ _motifColumns = columns; }
	void setTracked(){ _tracked = true; }
	void setWorstScore(double score){ _worstScore = score;}

	/* editing methods */
	void buildStatePWM(double *const *const _pwm);
	void buildPWM(double** pwm, double alpha, double* bg);
	void updatePWM_OOL(sorted_sites* sortedSites, double** pwm,
			motif_columns_type& motifColumns, int K, double pseudo);
	void updatePWM_OOL_initial(double pseudo);
	void updatePWM(double pseudo, bool maximizeMotifLength=false);
	void setBindingSites(const ThresholdChecker &pValThreshold);

  void makeReverseComplement();

	void updateEnrichedRegion(sorted_sites* bestStartPos);
	void updateEnrichedRegion();

	void filter_region();
	void filter_fitToSequence();
	void filter_oops_updatePWM();
	void updatePval(bool Debug = false);
	void updatePval2(bool Debug);
	void updateEmpiricalPval();

	/* calculations */
	double calculateMatrixScore(int seq, int startPos, int multSeq, const ss_type seqset=Global::posSet, const motif_columns_type *mcols = 0);
	double calculateMatrixScore_logOdds(int seq, int startPos);
	double getMaxMatrixScore();


	void printPWM(ostream &os, double** pwm) const;
	void printLogPWM(ostream &os, double** pwm, motif_columns_type& motifColumns, bool log) const;
	void printFullPWM(ostream &os) const;
	friend ostream& operator<< (ostream &os, const Motif& M);

	static bool cmpPtr(const Motif* m1, const Motif* m2) {
		if(m1->_pVal != m2->_pVal){
			return m1->_pVal < m2->_pVal;
		}else{
			const StartPosContainer& sp1 = m1->getStartPosList();
			const StartPosContainer& sp2 = m2->getStartPosList();

			if(sp1.size() != sp2.size()){
				return sp1.size() < sp2.size();
			}

			StartPosContainer::const_iterator it1 = sp1.begin();
			StartPosContainer::const_iterator it2 = sp2.begin();

			while(*it1 == *it2){
				it1++;
				it2++;
				if(it1 == sp1.end() && it2 == sp2.end()){
					return true;
				}
			}
			if(it1->seq != it2->seq) return it1->seq < it2->seq;
			return it1->pos < it2->pos;
		}
	}

	/* HO PWMs */

	void calculateHoProbs();

	void exp1merProbs();

	float** getConds(){
		return _conds;
	}

	float** getCounts(){
		return _counts;
	}

	int* getOffsets(){
		return _offsets;
	}

	int getOrder(){
		return _order;
	}

	void initMotifWithBindingSites( char* bindingSiteFile );

	void initMotifWithStartPos( bool multipleOccurrence );

	bool isLog(){
		return _log;
	}

	void log1merProbs();

	void multipleCounts( float factor );
	void resetCounts();

	void resetFactor(){
		_factor = static_cast<float>(_asize) * _alpha;
	}

	void setFactor( float factor ){
		_factor = factor;
	}

	void substractCountsOffset();

	static float* getFreqs(){
		return _freqs;
	}

	static void initHoMotif();

	/* Calculates array indices for sequence k-mers */
	static int sub2ind( unsigned char* sequence, int pos, int order, bool byrow=true );

private:
	int	_length;	    	/* length of motif model */
	region _enrichment;		/* set, my, startRegion and endRegion */

	bool _converged; /* has merging step changed converged matrix => iterate again */

	double _pVal; /* current pValue */
	double _worstScore; /* worst Score of Motif that has been chosen by order statistics */

	StartPosContainer _sites; /* stored start positions */
	MatchContainer _bindingSites; /* binding sites over a cerain threshold */

	double** _pwm; /* PWM for modelling motifs */

	motif_columns_type _motifColumns; /* list with important positions in PWM  */
	struct countStartPos_type{
		int* counts; /* counts of how many instances begin at which position in the sequence */
		int size;
	};
	countStartPos_type _countStartPos;

	int32_t _posSetSize;

	bool _tracked;

	void updateOptimizedPWM(double** pwm, StartPosContainer& sites, motif_type type, double alpha, double* bg);

	void updateCountsPWM(bool maximizeMotifLength=false);
	void updateCounts(int seq, int pos, double** pwm);
	void normalizeLOG(double pseudo, bool maximizeMotifLength = false);

	void fill_startPos_with_IUPAC_matches(const std::shared_ptr<Kmer>& skr);
	void fill_startPos_with_seeds(const std::shared_ptr<Kmer>& skr);
	void resetCountStartPos();

	/* HO PWMs */

	/* Calculates array indices for k-mers */
	static int sub2ind( unsigned char* sequence, int order, bool byrow=true );

	/* Calculates (conditional) probabilities for kmers */
	void calculateKmerProbs( unsigned char* kmer, int pos, int order_minus_1 );
	void calculateKmerProbs( unsigned char* kmer, int order_minus_1 );

	void updateHoCounts();
	void updateHoCounts( int sequenceNr, int startPos );

	/* Pseudocounts factor */
	static float _alpha;

	/* Alphabet size */
	static int _asize;

	/* Counts offset */
	static float _beta;

	/* Address coefficients */
	static int* _coeffs;

	/* Pseuodcounts factor */
	static float _factor;

	/* Field number in _conds/_counts/_probs arrays */
	static int _fields;

	/* Monomer background frequencies */
	static float* _freqs;

	/* Address offsets */
	static int* _offsets;

	/* Max. model order */
	static int _maxorder;

	/* HO kmer conditionals */
	float** _conds;

	/* HO kmer counts */
	float** _counts;

	/* fast_log (or linear) probs in _pwm for MONOnucleotides? */
	bool _log;

	/* Current model order */
	int _order;
};

inline void Motif::resetCounts(){

	for( int i=0; i <= _length; i++ ){
		for( int k=0; k < _offsets[_order+1]; k++ ){
			_counts[i][k] = 0.0f;
		}
	}
}

inline void Motif::multipleCounts( float factor ){

	for( int i=0; i <= _length; i++ ){
		for( int k=0; k < _offsets[_order+1]; k++ ){
			_counts[i][k] *= factor;
		}
	}
}

inline int Motif::sub2ind( unsigned char* sequence, int order, bool byrow ){

	int i, k, l;
	i = 0;

	if( byrow ){
		for( k=0, l=order; k<=order; ++k, --l ){
			i += _coeffs[l] * sequence[k];
		}
	}
	else{
		for( k=0; k<=order; ++k ){
			i += _coeffs[k] * sequence[k];
		}
	}

	return i;
}

inline int Motif::sub2ind( unsigned char* sequence, int pos, int order, bool byrow ){

	int i, k, l;
	i = 0;

	if( byrow ){
		for( k=0, l=order; k<=order; ++k, --l ){
			i += _coeffs[l] * sequence[pos+k];
		}
	}
	else{
		for( k=0; k<=order; ++k ){
			i += _coeffs[k] * sequence[pos+k];
		}
	}

	return i;
}

inline void Motif::resetCountStartPos(){
	memset(_countStartPos.counts, 0, sizeof(int)*Global::posSet->max_leng+1);
	_countStartPos.size = 0;
}

inline void Motif::normalizeLOG(double pseudo, bool maximizeMotifLength){
	/* normalize pwm and calculate log values */
	pseudo = getPseudocounts(pseudo, getTotalSites());

	int start = getFirstMotifColumn();
	int end = getLastMotifColumn();
	if(maximizeMotifLength){
		start = 1;
		end = _length;
	}

	double* bg = Global::negSet != NULL ? Global::negBg : Global::posBg;

	for(int i=start; i<=end; i++){
		double allPseudo = _pwm[i][0] + pseudo;
		for(int j=1; j<= nAlpha(Global::A); j++)
			_pwm[i][j] = fast_log(static_cast<float>((_pwm[i][j] + pseudo*bg[j]) / allPseudo));
	}
}

inline double Motif::calculateMatrixScore_logOdds(int seq, int startPos){
	double smin = 0;
	int firstColumn = getFirstMotifColumn();

	uint8_t* S = Global::posSet->entity[seq]->S[0];
	double* const &bg_log = Global::negSet != NULL ? Global::negBg_log : Global::posBg_log;
	for(motif_columns_type::const_iterator it = _motifColumns.begin(); it != _motifColumns.end(); it++){
		int base = S[startPos + *it - firstColumn];
		if(base==0){
			smin += (LogTable::LOG1_1000 + LogTable::LOG_i[2]);
			continue;
		}
		smin += _pwm[*it][base]-bg_log[base];
	}
	return smin;
}

inline double Motif::getMaxMatrixScore(){
	double maxMatrixScore = 0;
	double* const &bg_log = Global::negSet != NULL ? Global::negBg_log : Global::posBg_log;
	for(motif_columns_type::const_iterator it = _motifColumns.begin(); it != _motifColumns.end(); it++){
		double maxScore = 0;
		for(int base=1; base <= nAlpha(Global::A); base++){
			maxScore = std::max(maxScore, _pwm[*it][base]-bg_log[base]);
		}
		maxMatrixScore += maxScore;
	}
	return maxMatrixScore;
}

inline double Motif::calculateMatrixScore(int seq, int startPos, int multSeq, const ss_type seqset, const motif_columns_type *mcols){
	const motif_columns_type &columns = (mcols==NULL ? _motifColumns : *mcols);
	double smin = 0;
	const int firstColumn = columns.front();
	double* const &bg_log = Global::negSet != NULL ? Global::negBg_log : Global::posBg_log;
	if(multSeq == 1){
		uint8_t* S = seqset->entity[seq]->S[0];
		for(motif_columns_type::const_iterator it = columns.begin(); it != columns.end(); it++){
			int base = S[startPos + *it - firstColumn];
			//cerr << AlphaChar(base, Global::A);
			if(base==0){
				smin += (LogTable::LOG1_1000 + LogTable::LOG_i[2]);
				continue;
			}
			smin += (_pwm[*it][base]-bg_log[base]);
		}
		//cerr << endl;
	}else{
		int minmultSeq = std::min(multSeq, (int)seqset->entity[seq]->mseq);
		uint8_t** S = seqset->entity[seq]->S;
		for(motif_columns_type::const_iterator it = columns.begin(); it != columns.end(); it++){
			for(int l=0; l< minmultSeq; l++){
				int base = S[l][startPos + *it - firstColumn];
				if(base==0){
					smin += (LogTable::LOG1_1000 + LogTable::LOG_i[2]) / minmultSeq;
					continue;
				}
				smin += (_pwm[*it][base]-bg_log[base]) / minmultSeq;
			}
		}
	}
	return smin;
}

inline void Motif::offset(int offset){
	_enrichment.max -= offset;
	_enrichment.startRegion -= offset;
	_enrichment.endRegion -= offset;

	if(offset > 0){
		motif_columns_type::reverse_iterator it_column = _motifColumns.rbegin();
		for( ; it_column != _motifColumns.rend(); it_column++){
			for(int i=1; i<=nAlpha(Global::A); i++){
				_pwm[*it_column+offset][i] = _pwm[*it_column][i];
			}
		}
	}else{
		motif_columns_type::iterator it_column = _motifColumns.begin();
		for( ; it_column != _motifColumns.end(); it_column++){
			for(int i=1; i<=nAlpha(Global::A); i++){
				_pwm[*it_column+offset][i] = _pwm[*it_column][i];
			}
		}
	}

	motif_columns_type::iterator it_column = _motifColumns.begin();
	for( ;it_column != _motifColumns.end(); it_column++){
		*it_column += offset;
	}
}

ostream& operator<< (ostream &os, const Motif& M);

#endif /* MOTIF_H_ */
