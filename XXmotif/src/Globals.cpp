#include "Globals.h"
#include "backgroundDistribution.h"
#include "getopt_pp/getopt_pp.h"
#include "NullModel.h"
#include "aminoacids/NullModelCompositional.h"
#include "aminoacids/AbstractSupplementaryInformationProvider.h"
#include <limits.h>
#include <fstream>
#include <math.h>


#include "em/hoNullModel.h"

using GetOpt::GetOpt_pp;
using GetOpt::Option;
using GetOpt::OptionPresent;

a_type 		Global::A = NULL;
ss_type		Global::posSet = NULL;			/** positive sequence set **/
ss_type 	Global::negSet = NULL;			/** negative sequence set **/
bool		Global::usePositionalProbs = false;
bool    Global::useCompositionProbs = false;
int 		Global::backgroundOrder = 8;

int 		Global::GAPS = 0;				/** number of gap combinations allowed in the start motifs */
merge_type	Global::mergeMode; /* default is set below for proteins/dna */
double 		Global::gapOpening = 1;
double 		Global::gapExtension = 1;
int 		Global::maxMultipleSequences = 100;
int			Global::maxPosSetSize = INT_MAX;
ThresholdChecker		Global::instanceThreshold;

double 		Global::overrepCorrection=1.05;
double 		Global::consCorrection=1.00;
double 		Global::consPvalWeight=1.0/3;

int			Global::maxMotifsPerSequence = 1000;
int			Global::maxSeqCount = 1000;

bool 		Global::useAliFree = false;
bool 		Global::useRankPvalues = false;

bool   		Global::multipleOccurrence = false;	/** multiple occurences should be found in one sequence **/
bool   		Global::oneOccurrence = false;	/** multiple occurences should be found in one sequence **/
bool		Global::zeroOrOneOccurrence = false;
bool		Global::revcomp= false;			/** search on reverse complement of sequences **/
bool 		Global::repeatFiltering = false;
bool 		Global::lowComplexityFilter = false;
bool		Global::noRefinementPhase = false;

char* 		Global::startMotif = NULL;		/** start motif for motif discovery **/
char*		Global::profFile = NULL;		/** file with a the start profile for motif discovery **/
int			Global::startRegion = 0;		/** start of enriched region **/
int			Global::endRegion = 1;			/** end of enriched region **/
motif_type	Global::type = ALL;			/* type of motif: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM*/
seq_format	Global::seqFormat = FASTA;      /* format of input sequences: FASTA, CLUSTALW */

float***	Global::conservationProbs = NULL;
float***	Global::alignmentFreeProbs = NULL;
bool        Global::removeHomology = false;

double* 	Global::posBg_log = NULL;		/* logarithm of distribution of positive set **/
double* 	Global::posBg = NULL;			/* background of positive set **/
double* 	Global::negBg_log = NULL; 		/* logarithm of distribution of negative set **/
double*		Global::negBg = NULL;			/* background of negative set **/

double		Global::pseudo = 0.1;			/** pseudocounts for pwm for finding new instances**/
double		Global::plusFrac = 0;
int			Global::neff_pwm;				/** effictive number of different bases in one pwm column **/
int 		Global::neff_discrete;			/** effictive number of different bases before pwm phase **/

int			Global::downstream = 0;			/** distance between the alignment point and the end of the input sequences **/

char* 		Global::outputDirectory = NULL; /** output Directory **/
char*		Global::name = NULL;			/** input file name **/
char*		Global::shortFileName = NULL;	/** input file name without path and .**/
char*		Global::negFile = NULL;			/** negative input file name **/
char* 		Global::benchmarkFolder = NULL; /** folder in which the benchmark results are stored **/
char*		Global::pwmFolder = NULL;
int			Global::maxMotifLevel = 3;		/** max number of extensions tried per level */
double		Global::minCoverage;

int			Global::minMatchPositions;
int			Global::maxMatchPositions;

bool		Global::aa = false;
bool		Global::maximizeMotifLength = true;

StateType 	Global::type_of_states;

/* (begin) options for EM-learning IMMs (author: siebert) */
bool		Global::em = false;					/* IMM EM mode */
char*		Global::bindingSiteFile = NULL;		/* (input) file with binding sites for initializing HO model */
int			Global::bindingSiteLength = 30;		/* length of motif model necessary to model given binding sites */

int			Global::modelOrder = 2;				/* (max) order of HO model */

float		Global::q = 0.1f;					/* weak binding parameter */
float		Global::eta = 18.0f;				/* substitutes alpha parameter */

float		Global::alpha = eta * q / ( 1-q );	/* pseudocounts factor of HO model (no command line parameter) */
float		Global::beta = 0.0f;				/* counts offset of HO model (no command line parameter) */

int			Global::modelOrderBg = 2;			/* order of HO background model */
float		Global::alphaBg = 10.0f;			/* pseudocounts factor of HO background model */
float		Global::betaBg = 0.0f;				/* counts offset of HO background model (no command line parameter) */

bool		Global::cv = false;					/* IMM EM CV mode */
char*		Global::testFile = NULL;			/* input file name for test sequences to score model */

char*		Global::shortTestFileName = NULL;	/* input file name without file path and file ending (no command line parameter) */
std::string	Global::initFile;					/* output file with binding sites for THE ULTIMATE k-mer (no command line parameter) */
ss_type 	Global::testSet = NULL;				/* test sequences to score model (no command line parameter) */

float*		Global::freqs = NULL;				/* monomer background frequencies (no command line parameter) */
bool		Global::gaps = true;				/* calculate background (conditional) probabilities for kmers with gaps (no command line parameter) */
bool		Global::verbose = false;			/* printout models */
bool		Global::save = false;				/* save models */

/* (end) options for EM-learning IMMs (author: siebert) */

/* HO null model (used by Holger & Eckhart) */
int			Global::order = 2;					/* order of HO null model */
float		Global::pseudocountsFactor = 10.0f;	/* pseudocounts factor of HO null model */
float		Global::countsOffset = 0.0f;		/* counts offset of HO null model */

/* cs blast */
int			Global::cswlen;
std::string	Global::csprofiles;
int 		Global::csbest;

SuppInfMode Global::suppInfMode = SUPP_NO;
double      Global::disoconsWeight = 1.0/3;
float		Global::dollarScaleFactor = 3.0f;
TerminusMode Global::termMode = BOTH;
std::string Global::aaMtfFile;
int			Global::maxIterations = std::numeric_limits<int>::max();
Global::tracked_t Global::trackedMotifs;
bool 		Global::trackedElongation = false;
bool		Global::trackedOnly = false;
double		Global::extensionECut;
int			Global::extensionMinCut;
int			Global::extensionMaxCut;
double		Global::aaStateSigThresh = 2.0;
double		Global::aaSeqFreqThresh = 0.75;
bool		Global::batch_mode = false;

bool 		Global::fixedPosition = false;
double 		Global::finalFilterThreshold = 1e2;

AbstractSupplementaryInformationProvider const* Global::supplInf;
std::string Global::nnetFilename;

bool 		Global::DEBUG = false;
std::string Global::argv0;

Global::Global(int argc, char *argv[]){
	argv0 = std::string(argv[0]);

	if(argc < 2) printHelpOutput();
	seq_format format = FASTA;
	for(long i = 2; i < argc; i++){
		if(strcmp(argv[i], "--aa") == 0){
			A = MkAlpha("XACDEFGHIKLMNPQRSTVWY$");
			aa = true;
		}else if(strcmp(argv[i], "--terminus-mode") == 0){
			if(i+1 < argc && strcmp(argv[i+1], "NONE") == 0) termMode = NONE;
			else if(i+1 < argc && strcmp(argv[i+1], "POS") == 0) termMode = POS;
			else if(i+1 < argc && strcmp(argv[i+1], "NEG") == 0) termMode = NEG;
			else if(i+1 < argc && strcmp(argv[i+1], "BOTH") == 0) termMode = BOTH;
			else { fprintf(stderr, "Illegal terminus mode\n"); exit(1); }
		}else if(strcmp(argv[i], "--maxPosSetSize") == 0){
			if(i+1 == argc){ fprintf(stderr, "\n\nERROR: no value given for option --maxPosSetSize\n\n"); exit(-1); }
			maxPosSetSize = atoi(argv[i+1]);
		}else if(strcmp(argv[i], "--format") == 0){
			if(i+1 == argc){ fprintf(stderr, "\n\nERROR: no value given for option --format\n\n"); exit(-1); }
			if(i+1 < argc && strcmp(argv[i+1], "CLUSTALW") == 0) format = CLUSTALW;
			if(i+1 < argc && strcmp(argv[i+1], "MFASTA") == 0) format = MFASTA;
			if(i+1 < argc && strcmp(argv[i+1], "CUSTOM") == 0) format = CUSTOM;
		}else if(strcmp(argv[i], "--negSet") == 0){
			negFile = argv[i+1];
		}else if(strcmp(argv[i], "--lcf") == 0){
			lowComplexityFilter = true;
		}
	}

	struct stat sts;
	if (((stat (argv[2], &sts)) == -1) || S_ISDIR(sts.st_mode)){
		fprintf(stderr, "\n !!! SEQFILE %s does not exist !!! \n\n", argv[2]);
		exit(-1);
	}
	if (negFile != NULL && (((stat (negFile, &sts)) == -1) || S_ISDIR(sts.st_mode))){
		fprintf(stderr, "\n !!! NEGFILE %s does not exist !!! \n\n", negFile);
		exit(-1);
	}

	if (!aa) termMode=NONE;

	if(A==NULL) A = MkAlpha("NACGT");

	posSet = readSeqSet(argv[2], A, format, maxPosSetSize, termMode==BOTH || termMode==POS);
	calculateSequenceFeatures(posSet, A);

    posBg_log = (double*)calloc(nAlpha(A)+1, sizeof(double));
    posBg     = (double*)calloc(nAlpha(A)+1, sizeof(double));
    negBg_log = (double*)calloc(nAlpha(A)+1, sizeof(double));
    negBg 	  = (double*)calloc(nAlpha(A)+1, sizeof(double));

	startRegion = 1;
	endRegion = Global::posSet->max_leng;

	if(!readCommandLineOptions(argc, argv))	printHelpOutput();
	if(aa){
		if(negFile == NULL){ fprintf(stderr, "no negative set !!! \n"); exit(-1);}
		negSet = readSeqSet(negFile, A, format, INT_MAX, termMode==BOTH || termMode==NEG);
		calculateSequenceFeatures(negSet, A);
	}else{
		if(negFile != NULL){
			negSet = readSeqSet(negFile,A,format, INT_MAX);
		}
		if(negSet != NULL ){
			if(repeatFiltering){ filter_repeats(negSet, A); }
			if(lowComplexityFilter){ cerr << "negSet: "; filter_lowComplexity(negSet, A); }
			if(revcomp){ createRevcomp(negSet, A); }
			fillGapsInMultipleAlignment(negSet, A);
			filterMaxMultipleSequences(negSet, maxMultipleSequences);
			calculateSequenceFeatures(negSet, A);
		}

		if(repeatFiltering){ filter_repeats(posSet, A); }
		if(lowComplexityFilter){ cerr << "posSet: "; filter_lowComplexity(posSet, A); }
		if(revcomp){ createRevcomp(posSet, A); }
		fillGapsInMultipleAlignment(posSet, A);
		filterMaxMultipleSequences(posSet, maxMultipleSequences);

	}

//    for(int i=1; i<=posSet->nent; i++){
//       e_type ali = posSet->entity[i];
//       for(int m=0; m<ali->mseq; m++){
//          fprintf(stderr, "%d: ", m);
//          for(int pos = 1; pos <= ali->n; pos++){
//             /* overwrite all gap positions in not target specie with nucleotides from more distantly related species */
//             fprintf(stderr, "%c", AlphaChar(ali->S[m][pos], A));
//          }fprintf(stderr, "\n");
//       }fprintf(stderr, "\n");
//    }
//    for(int i=1; i<=negSet->nent; i++){
//       e_type ali = negSet->entity[i];
//       for(int m=0; m<ali->mseq; m++){
//          fprintf(stderr, "%d: ", m);
//          for(int pos = 1; pos <= ali->n; pos++){
//             fprintf(stderr, "%c", AlphaChar(ali->S[m][pos], A));
//          }fprintf(stderr, "\n");
//       }fprintf(stderr, "\n");
//    }
	//exit(-1);

	checkSequenceSet(); // do some tests whether the negative set is a good choice
	setBackgroundDistribution(); // get trinuc distribution and conditional trinuc distribution of the background set
	if (useCompositionProbs) NullModelCompositional::initialize();
}

void Global::printHelpOutput(){
	bool developerMode = false;

	printf("\n==============================================================================================================================\n");
	printf("== XXmotif version 1.6");
	printf("\n==============================================================================================================================\n");
	printf("\nUsage: XXmotif OUTDIR SEQFILE [options] \n\n");
	printf("\tOUTDIR:  output directory for all results\n");
	printf("\tSEQFILE: file name with sequences from positive set in FASTA format\n");
	printf("\n");
	printf("Options:\n");
	printf("\t--negSet <FILE>\t\t\t\tsequence set which has to be used as a reference set\n");
	printf("\t--zoops\t\t\t\t\tuse zero-or-one occurrence per sequence model (DEFAULT)\n");
	printf("\t--mops\t\t\t\t\tuse multiple occurrence per sequence model\n");
	printf("\t--oops\t\t\t\t\tuse one occurrence per sequence model\n");
	printf("\t--revcomp\t\t\t\tsearch in reverse complement of sequences as well (DEFAULT: NO)\n");
	printf("\t--background-model-order <NUMBER>\torder of background distribution (DEFAULT: 2, 8(--negset) )\n");//, 4(--aa) )\n");
	if(developerMode)printf("\t--write-pwm-file FOLDER\t\t\twrite PWMs into this folder (DEFAULT OUTDIR)\n");
	printf("\t--pseudo <NUMBER>\t\t\tpercentage of pseudocounts used (DEFAULT: 10)\n");
	printf("\t-g|--gaps <NUMBER>\t\t\tmaximum number of gaps used for start seeds [0-3] (DEFAULT: 0)\n");
	printf("\t--type <TYPE>\t\t\t\tdefines what kind of start seeds are used (DEFAULT: ALL)\n");
	printf("\t\t\t\t\t\t - possible types: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM\n");
	printf("\t--merge-motif-threshold <MODE>\t\t\tdefines the similarity threshold for merging motifs (DEFAULT: HIGH)\n");
	printf("\t\t\t\t\t\t - possible modes: LOW, MEDIUM, HIGH\n");
	if(developerMode)printf("\n");
	printf("\t--no-pwm-length-optimization\t\tdo not optimize length during iterations (runtime advantages)\n");
	if(developerMode)printf("\t--min-match-positions <INT>\t\tmin number of non-wildcard positions per motif (DEFAULT: 4 (NA), 2 (AA))\n");
	printf("\t--max-match-positions <INT>\t\tmax number of positions per motif (DEFAULT: 17, higher values will lead to very long runtimes)\n");
	printf("\n");
	printf("\t--batch\t\t\t\t\tsuppress progress bars (reduce output size for batch jobs)\n");
	if(developerMode)printf("\t--maxSeqOcc\t\t\tmaximum number of motif occurrences per sequence)\n");
	if(developerMode)printf("\t--debug\t\t\t\t\tshow matrix for every iteration step\n");
	printf("\t--maxPosSetSize <NUMBER>\t\tmaximum number of sequences from the positive set used [DEFAULT: all]\n");
	printf("\t--help\t\t\t\t\tprint this help page\n");
	printf("\t--trackedMotif <SEED>\t\t\tinspect extensions and refinement of a given seed (DEFAULT: not used)\n");
	if(developerMode)printf("\t--neff-states <NUMBER>\t\t\teffective number of different states in one IUPAC extension (DEFAULt: 6)\n");
	if(developerMode)printf("\t--neff-pwm <NUMBER>\t\t\teffective number of different states in one PWM column (DEFAULt: 10 (NA), 63(AA))\n");
	if(developerMode)printf("\t--gapOpening <NUMBER>\t\t\tbit penalty for every opened gap\n");
	if(developerMode)printf("\t--gapExtension <NUMBER>\t\t\tbit penalty for every extended gap position\n");
	printf("\n");
	printf("Using conservation information\n");
	printf("\t--format FASTA|MFASTA\t\t\tdefines what kind of format the input sequences have (DEFAULT: FASTA)\n");
	printf("\t--maxMultipleSequences <NUMBER>\t\tmaximum number of sequences used in an alignment [DEFAULT: all]\n");
	if(developerMode)printf("\t--cons-length <NUMBER>\t\tused nucleotides for conservation pVal calculation [DEFAULT: 8]\n");
	printf("\n");
	printf("Using localization information\n");
	printf("\t--localization\t\t\t\tuse localization information to calculate combined P-values \n\
			\t\t\t(sequences should have all the same length)\n");
	printf("\t--downstream <NUMBER>\t\t\tnumber of residues in positive set downstream of anchor point (DEFAULT: 0)\n");
	printf("\n");
	printf("Start with self defined motif:\n");
	printf("\t-m|--startMotif <MOTIF>\t\t\tStart motif (IUPAC characters)\n");
	printf("\t-p|--profileFile <FILE>\t\t\tprofile file\n");
	printf("\t--startRegion <NUMBER>\t\t\texpected start position for motif occurrences relative to anchor point (--localization)\n");
	printf("\t--endRegion <NUMBER>\t\t\texpected end position for motif occurrences relative to anchor point (--localization)\n");
	if( developerMode ) printf( "\n" );
	if( developerMode ) printf( "HO null model:\n" );
	if( developerMode ) printf( "\t--counts-offset <FLOAT>\t\t\tpseudocounts factor of HO null model\n" );
	if( developerMode ) printf( "\t--pseudocounts-factor <FLOAT>\t\tcounts offset of HO null model\n" );
	printf( "\n" );
	if(developerMode)printf("\n");
	if(developerMode)printf("\t--empirical-recalibration\t\trecalibrate pValues with negative set\n");
	if(developerMode)printf("\t--min-coverage <FLOAT>\t\t\tminimum fraction of sequences a motif has to be found in (DEFAULT: 0.0)\n");
	/*printf("\n");
	printf("Proteins only:\n");
	printf("\t--aa\t\t\t\t\tuse amino acids\n");
	printf("\t--aaMtfFile <STRING>\t\t\toutput file for AA motifs in MTF format (DEFAULT: directory/seqfile.mtf)\n");
	if(developerMode)printf("\t--extensionMinCut <INT>\t\t\tmin number of initial seeds extended (DEFAULT: 1000)\n");
	if(developerMode)printf("\t--extensionMaxCut <INT>\t\t\tmax number of initial seeds extended (DEFAULT: infty)\n");
	if(developerMode)printf("\t--extensionECut <FLOAT>\t\t\tE-value threshold for seed extension (DEFAULT: 1.0)\n");
	if(developerMode)printf("\t--aaStateSigThresh <FLOAT>\t\todds-threshold for state significance (DEFAULT: 2.0)\n");
	if(developerMode)printf("\t--aaSeqFreqThresh <FLOAT>\t\tconservation threshold for extension positions (DEFAULT: 0.75)\n");
	printf("\t--supplementary-information [NO|DISOCONS|NNET]\t\t\t\tuse supplementary information (DEFAULT: NO)\n");
	if(developerMode)printf("\t--disoconsWeight\t\t\tweight of disorder/conservation P-value (DEFAULT: 0.5)\n");
	if(developerMode)printf("\t--termFreqScale <NUMBER>\t\tscale frequencies of terminal character by NUMBER (DEFAULT: 3.0)\n");
	if(developerMode)printf("\t--trackedMotif <STRING>\t\t\toutput info on tracked motif for debugging (DEFAULT: NONE)\n");
	if(developerMode)printf("\t--trackedOnly\t\t\t\tdrop all but the tracked motifs (DEFAULT: no)\n");
	*/
	printf("\n==============================================================================================================================\n");
	printf("\n");
	exit(-1);
}

bool Global::readCommandLineOptions(int argc, char *argv[]){
	GetOpt_pp ops(argc,argv);
	if (ops >> OptionPresent('h', std::string("help"))) return false;

	/* process non-option arguments (dir and seqfile) */
	std::vector<std::string> args;
	ops >> Option(GetOpt_pp::EMPTY_OPTION, args);
	if (!(args.size()==2)) return false;
	int i;
	outputDirectory = String(args[0].c_str());
	createDirectory(outputDirectory);

	char* tmpFolder = (char*)calloc(1024, sizeof(char));
	sprintf(tmpFolder, "%s/tmp", outputDirectory);
	createDirectory(tmpFolder);
	free(tmpFolder);

	pwmFolder = String(outputDirectory);
	name = String(args[1].c_str());
	//printf("name: %s\t", name);
	/* extract filename from posSeqFile */
	i = 0;
	int start = 0, end = 0;
	while (name[++i] != '\0') {
		if (name[i] == '.') end = i - 1;
	}
	while (--i != 0 && Global::name[i] != '/');
	if (i == 0)	start = 0;
		else start = i + 1;
	char* fileName = (char*)malloc((end-start+2)*sizeof(char));
	//NEW(fileName, end-start+2, char);
	for (i = start; i <= end; i++) {
		fileName[i - start] = name[i];
	}
	fileName[i - start] = '\0';
	shortFileName = fileName;
	//printf("shortFileName: %s\n", shortFileName);

	/* flags */
	if (ops >> OptionPresent(' ', std::string("debug"))) DEBUG = 1;
	ops >> OptionPresent(' ', std::string("mops"), multipleOccurrence);
	ops >> OptionPresent(' ', std::string("oops"), oneOccurrence);
	ops >> OptionPresent(' ', std::string("zoops"), zeroOrOneOccurrence);
	if((oneOccurrence && multipleOccurrence) || (oneOccurrence && zeroOrOneOccurrence) || (zeroOrOneOccurrence && multipleOccurrence)){
		fprintf(stderr, "\n\n !!! one-occurrence per sequence model (--oops) cannot be used ");
		fprintf(stderr, "similtaniously with the multiple occurrence per sequence model (--mops) or zeror-or-one occurrence per sequence model (--zops) !!!\n\n");
		exit(-1);
	}

	ops >> OptionPresent(' ', std::string("localization"), usePositionalProbs);
	if(usePositionalProbs && (Global::posSet->max_leng != Global::posSet->min_leng)){
		fprintf(stderr, "\n\n !!! localization information can only be used if all input sequences have the same length !!!\n\n");
		exit(-1);
	}
	ops >> OptionPresent(' ', std::string("aa"), aa);
	ops >> OptionPresent(' ', std::string("revcomp"), revcomp);
	ops >> OptionPresent(' ', std::string("noRefinementPhase"), noRefinementPhase);
	ops >> OptionPresent(' ', std::string("ranks"), useRankPvalues);
	ops >> OptionPresent(' ', std::string("ali-free"), useAliFree);
	ops >> OptionPresent(' ', std::string("lcf"), lowComplexityFilter);
	if (ops >> OptionPresent(' ', std::string("filtering"))) repeatFiltering = true;

	if (aa) {
		std::string suppl_opt;
		ops >> Option(' ', std::string("supplementary-information"), suppl_opt, "");
    ops >> Option(' ', std::string("nnet-file"), nnetFilename, "");
		if (suppl_opt == "" || suppl_opt == "SUPP_NO") {
			suppInfMode = SUPP_NO;
		} else if (suppl_opt == "SUPP_DISOCONS") {
			suppInfMode = SUPP_DISOCONS;
		} else if (suppl_opt == "SUPP_NNET") {
			suppInfMode = SUPP_NNET;
			if (nnetFilename == "") {
			  std::cerr << "No neural net specified for supplementary information, use option --nnet-file" << std::endl;
			  exit(1);
			}
		} else {
			std::cerr << "Unknown supplementary information type: " << suppl_opt << std::endl;
			exit(1);
		}
		if (DEBUG) {
			std::cerr << "Supplementary information: " << suppInfMode << std::endl;
		}
		ops >> OptionPresent(' ', "trackedOnly", trackedOnly);
	}else{
	  suppInfMode = SUPP_NO;
	}

	/* options with argument */
	if (ops >> Option('m', std::string("startMotif"), startMotif)) printf("Start motif : %s\n",startMotif);
	if (ops >> Option('p', std::string("profileFile"), profFile)) printf("Profile file: %s\n", profFile);
	ops >> Option(' ', std::string("downstream"), downstream);
 	int offset = Global::posSet->max_leng - Global::downstream;
	if( ops >> Option(' ', std::string("startRegion"), startRegion)) startRegion += offset;
	else startRegion = 0;
	if( ops >> Option(' ', std::string("endRegion"), endRegion)) Global::endRegion += offset;
	else endRegion = posSet->max_leng;

	ops >> Option(' ', std::string("type"), type); if(type == NO_VALID_MOTIF_TYPE) return false;
	ops >> Option(' ', std::string("format"), seqFormat); if(seqFormat == NO_VALID_SEQ_FORMAT) return false;
	mergeMode = aa ? MEDIUM : HIGH;
	if (DEBUG) std::cout << "Merge motif mode default: " << mergeMode << std::endl;
	ops >> Option(' ', std::string("merge-motif-threshold"), mergeMode); if(mergeMode == NO_VALID_MERGE_MODE) return false;

	ops >> OptionPresent(' ', std::string("extra-homology-filter"), removeHomology);
	if (DEBUG) printf("Extra homology filter: %sabled\n", removeHomology ? "en" : "dis");

	if (ops >> Option(' ', std::string("pseudo"), pseudo)) {
		printf("pseudocounts: %.1f%%\n", pseudo);
		pseudo /= 100;
	}
	if (ops >> Option(' ', std::string("plusFrac"), plusFrac)) {
		printf("plus fraction: %.1f%%\n", plusFrac);
		plusFrac /= 100;
	}

	if (ops >> Option(' ', std::string("negSet"), negFile)) printf("negFileName: %s\n", negFile);
	if (ops >> Option(' ', std::string("benchmarkFolder"), benchmarkFolder)) printf("benchmark Folder: %s\n", benchmarkFolder);
	if (ops >> Option(' ', std::string("write-pwm-file"), pwmFolder)) printf("pwm folder: %s\n", pwmFolder);
	std::string dummyMode;
	if (ops >> Option(' ', std::string("terminus-mode"), dummyMode)) {
		/* already handled above before reading posSet */
	}
	if(aa && DEBUG) {
		printf("termini handling: ");
		switch (termMode) {
			case NONE: printf("no $ added\n"); break;
			case POS: printf("$ added to posset only\n"); break;
			case NEG: printf("$ added to negset only\n"); break;
			case BOTH: printf("$ added to posset and negset\n"); break;
		}
	}
	if (ops >> Option(' ', std::string("termFreqScale"), dollarScaleFactor) || DEBUG) printf("scaling termini frequencies by: %f\n", dollarScaleFactor);
	if (ops >> Option('g', std::string("gaps"), GAPS)) printf("maximum number of initial gaps: %d\n", GAPS);
	if (ops >> Option(' ', std::string("gapOpening"), gapOpening)) printf("gapOpening penalty: %f\n", gapOpening);
	if (ops >> Option(' ', std::string("gapExtension"), gapExtension)) printf("gapExtension penalty: %f\n", gapExtension);
	if (ops >> Option(' ', std::string("maxMultipleSequences"), maxMultipleSequences)) printf("maximal used sequences: %d\n", maxMultipleSequences);
	if (ops >> Option(' ', std::string("maxPosSetSize"), maxPosSetSize)) printf("maximal used sequences: %d\n", maxPosSetSize);
	if (ops >> Option(' ', std::string("maxSeqOcc"), maxMotifsPerSequence)) printf("maximal motifs per sequence: %d\n", maxMotifsPerSequence);

	if (ops >> Option(' ', std::string("maxMotifLevel"), maxMotifLevel)) printf("max motif level: %d\n", maxMotifLevel);

	/* (begin) options for EM-learning IMMs (author: siebert) */

	if( ops >> OptionPresent( ' ', std::string( "em" ) ) ){
		em = true;
		if( DEBUG ){
			printf( "-------------------------------------------------------\n" );
			printf( "IMM EM mode options\n" );
			printf( "-------------------------------------------------------\n" );
		}
	}

	if( em ){

		if( ops >> Option( ' ', std::string( "bindingSiteFile" ), bindingSiteFile ) ){
			if( DEBUG ){
				printf( "bindingSiteFile: %s\n", bindingSiteFile );
			}
			if( ops >> Option( ' ', std::string( "bindingSiteLength" ), bindingSiteLength ) ){
				if( DEBUG ){
					printf( "bindingSiteLength: %d\n", bindingSiteLength );
				}
			} else{
				fprintf( stderr, "Please specify the length of binding site sequences (--bindingSiteLength). Sequence lengths must not differ. \n");
				/* maybe better to read in 1st sequence and determine length automatically
				 * check later whether all sequences correspond to this length
				 * need to be changed in case the model can be elongated during the EM scheme */
				exit( -1 );
			}
		} else if( DEBUG ){
			printf( "bindingSiteFile: no" );
			printf( "bindingSiteLength: %d\n", bindingSiteLength );
			/* use instances of (subsets of) XXmotif's IUPAC patterns to initialize IMMs */
		}
		/* TODO: option to read IMM directly from file */

		if( ops >> Option( ' ', std::string( "modelOrder" ), modelOrder ) && DEBUG ){
			printf( "modelOrder: %d\n", modelOrder );
		} else if( DEBUG ){
			printf( "modelOrder: %d (default)\n", modelOrder );
		}
		if( ops >> Option( ' ', std::string( "q" ), q ) && DEBUG ){
			printf( "q: %.2f\n", q );
		} else if( DEBUG ){
			printf( "q: %.2f (default)\n", q );
		}
		if( ops >> Option( ' ', std::string( "eta" ), eta ) && DEBUG ){
			printf( "eta: %.2f\n", eta );
		} else if( DEBUG ){
			printf( "eta: %.2f (default)\n", eta );
		}

		printf( "alpha: %.2f (eta*(q/(1-q)))\n", alpha );

//		if( ops >> Option( ' ', std::string( "alpha" ), alpha ) && DEBUG ){
//			printf( "alpha: %.2f\n", alpha );
//		} else if( DEBUG ){
//			printf( "alpha: %.2f (default)\n", alpha );
//		}
//		if( ops >> Option( ' ', std::string( "beta" ), beta ) && DEBUG ){
//			printf( "beta: %.2f\n", beta );
//		} else if( DEBUG ){
//			printf( "beta: %.2f (default)\n", beta );
//		}

		if( ops >> Option( ' ', std::string( "modelOrderBg" ), modelOrderBg ) && DEBUG ){
			printf( "modelOrderBg: %d\n", modelOrderBg );
		} else if( DEBUG ){
			printf( "modelOrderBg: %d (default)\n", modelOrderBg );
		}
		if( ops >> Option( ' ', std::string( "alphaBg" ), alphaBg ) && DEBUG ){
			printf( "alphaBg: %.2f\n", alphaBg );
		} else if( DEBUG ){
			printf( "alphaBg: %.2f (default)\n", alphaBg );
		}

//		if( ops >> Option( ' ', std::string( "betaBg" ), betaBg ) && DEBUG ){
//			printf( "betaBg: %.2f\n", betaBg );
//		} else if( DEBUG ){
//			printf( "betaBg: %.2f (default)\n", betaBg );
//		}
	}

	if( ops >> OptionPresent( ' ', std::string( "cv" ) ) ){
		cv = true;
		if( DEBUG ){
			printf( "-------------------------------------------------------\n" );
			printf( "IMM EM CV mode options\n" );
			printf( "-------------------------------------------------------\n" );
		}
	}

	if( cv ){

		std::stringstream s;
		s << Global::outputDirectory << "/" << Global::shortFileName << "-CV.ini";
		initFile = s.str();

		if( ops >> Option( ' ', std::string( "testSet" ), testFile ) ){

			if( DEBUG ){
				printf( "testSet: %s\n", testFile );
				printf( "initFile: %s\n", initFile.c_str() );
			}
			testSet = readSeqSet( testFile, Global::A, FASTA, INT_MAX );

			int pos = 0, startPos = 0, lastPos = 0;

			while( testFile[++pos] != '\0' )
				if( testFile[pos] == '.' )
					lastPos = pos-1;
			while( --pos != 0 && testFile[pos] != '/' )
				;
			if( pos == 0 )
				startPos = 0;
			else
				startPos = pos+1;

			shortTestFileName = ( char* )malloc( ( lastPos-startPos+2 )*sizeof( char ) );

			for( pos=startPos; pos <= lastPos; pos++ )
				shortTestFileName[ pos-startPos ] = testFile[pos];
			shortTestFileName[ pos-startPos ] = '\0';
		} else if( DEBUG ){
			fprintf( stderr, "Please provide a FASTA sequence file (--testSet).\n");
			exit( -1 );
		}
	}

	if( ops >> OptionPresent( ' ', std::string( "verbose" ) ) ){
		verbose = true;
		if( DEBUG ){
			printf( "verbose: yes\n" );
		}
	} else if ( DEBUG ){
		printf( "verbose: no\n" );
	}

	if( ops >> OptionPresent( ' ', std::string( "save" ) ) ){
		save = true;
		if( DEBUG ){
			printf( "save: yes\n" );
		}
	} else if ( DEBUG ){
		printf( "save: no\n" );
	}

	if( DEBUG ){
		printf( "-------------------------------------------------------\n" );
	}

	/* (end) options for EM-learning IMMs (author: siebert) */

	/* options for cs version */
	ops >> Option(' ', std::string("cswlen"), cswlen, 13);
	ops >> Option(' ', std::string("csbest"), csbest, 0);
	ops >> Option(' ', std::string("csprofiles"), csprofiles);
	if(csbest != 0){
		printf("cs len: %d\n", cswlen);
		printf("cs best: %d\n", csbest);
		if(csprofiles.length() == 0){ printf("ERROR: no cs profiles output file given.\n\t option --csprofiles\n"); return false; }
		printf("cs profiles output file: %s\n", csprofiles.c_str());
	}

	double min_cov;
	if (ops >> Option(' ', std::string("min-coverage"), min_cov)) {
		if (min_cov <= 1.0) {
			minCoverage = min_cov * Global::posSet->nent;
		} else {
			minCoverage = min_cov;
		}
	} else {
		minCoverage = 0;
	}
	if (aa && DEBUG) printf("minimum coverage: %.0f sequences (%.1f%%)\n", minCoverage, minCoverage / Global::posSet->nent * 100);
	if (minCoverage > Global::posSet->nent) {
		printf("Minimum coverage is greater than the number of sequences. Do you think this makes sense?\n");
		exit(1);
	}

	if (ops >> Option(' ', std::string("min-match-positions"), minMatchPositions)) {
	} else { minMatchPositions = aa ? 2 : 5;	}
	if (DEBUG) printf("Min match positions: %d\n", minMatchPositions);

	if (ops >> Option(' ', std::string("max-match-positions"), maxMatchPositions)) {
		if(maxMatchPositions > 26){
			printf("maxMatchPositions is only not allowed to be greater than 26 !!!\n");
			exit(1);
		}
	} else { maxMatchPositions = aa ? 10 : 17;	}
	if (DEBUG) printf("Max match positions: %d\n", maxMatchPositions);

	backgroundOrder = (negFile == NULL) ? 2 : 8;
	if(aa) backgroundOrder = 4;

	if(ops >> Option(' ', std::string("background-model-order"), backgroundOrder) || DEBUG){
		printf("order of background probability model: %d\n", backgroundOrder);
	}

	/* HO null model */

	order = ( negFile == NULL )? 2 : 8;
	if( aa ) order = 4;

	bool setOffset = false;
	if( ops >> Option(' ', std::string( "counts-offset" ), countsOffset ) )
		setOffset = true;
	if( ops >> Option(' ', std::string( "pseudocounts-factor" ), pseudocountsFactor ) )
		setOffset = true;
	if( setOffset || DEBUG){
		printf( "counts offset of HO null model: %f\n", countsOffset );
		printf( "pseudocounts factor of HO null model: %f\n", pseudocountsFactor );
	}

	int comp_opt;
	bool useComp = (aa ? 1 : 0);
	if (ops >> Option(' ', std::string("use-composition"), comp_opt)) {
		useComp = (comp_opt!=0);
	}
	if (aa && DEBUG) printf("using sequence composition for probability calculations: %s\n", useComp ? "yes" : "no");

	std::string thresh_init;
	if (ops >> Option(' ', std::string("instance-threshold"), thresh_init)) {
		instanceThreshold = ThresholdChecker(thresh_init);
	} else {
		instanceThreshold = ThresholdChecker("<1.0");
	}

	if(aa && (ops >> Option(' ', std::string("fixed-position"), fixedPosition) || DEBUG)){
		printf("fixed position %sabled\n", fixedPosition ? "en" : "dis");
	}

	if(ops >> Option(' ', std::string("final-filter"), finalFilterThreshold) || DEBUG){
		finalFilterThreshold = log(finalFilterThreshold);
		printf("E-value threshold for final filtering: %g\n", exp(finalFilterThreshold));
	}else{
		finalFilterThreshold = log(finalFilterThreshold);
	}
	/* options for protein version */
	ops >> Option(' ', std::string("maxIterations"), maxIterations);

	std::stringstream s;
	s << Global::outputDirectory << "/" << Global::shortFileName << ".mtf";
	aaMtfFile = s.str();
	if( ops >> Option(' ', std::string("mtf-file"), aaMtfFile) || DEBUG){
		printf("mtf output file: %s\n", aaMtfFile.c_str());
	}

	std::string states;
	if (aa && ops >> Option(' ', std::string("states-type"), states, "BLOSUM62_21")) {
		type_of_states = StateType::fromString(states);
		if (type_of_states.type == UNDEF) {
			std::cerr << "Invalid states type " << states << std::endl;
			exit(1);
		}
		std::cout << "State Alphabet: " << StateType::toString(type_of_states) << std::endl;
	}

	if (aa) {
		ops >> Option(' ', std::string("use-composition"), useCompositionProbs, false);
		if( DEBUG ) fprintf(stderr, "Compositional probabilities: %sabled\n", useCompositionProbs ? "en" : "dis");
	}

	/* N_eff can be set by option --neff
	 * or seperately for discrete and pwm phase */
	neff_discrete = (aa ? StateType::sizeOfAlphabet(type_of_states) : 6 );
	neff_pwm = (aa ? neff_discrete : 10);

	int n_eff = -1;
	bool setNeff = false;
	if (ops >> Option(' ', std::string("neff"), n_eff)) {
		neff_pwm = n_eff;
		neff_discrete = n_eff;
		setNeff = true;
	}
	if (ops >> Option(' ', std::string("neff-states"), neff_discrete)) 	setNeff = true;
	if (ops >> Option(' ', std::string("neff-pwm"), neff_pwm)) setNeff = true;
	if(setNeff || DEBUG) printf("N_eff discrete: %d\t\tN_eff pwm: %d\n", neff_discrete, neff_pwm);

	if (aa && (ops >> Option(' ', std::string("disoconsWeight"), disoconsWeight) || DEBUG)) {
		printf("disocons weight for P value combination: %f\n", disoconsWeight);
		if (disoconsWeight<=0 || disoconsWeight >=1) {
			std::cerr << "ERROR: P-value combination weight must be in ]0,1[." << std::endl;
			exit(1);
		}
	}

	std::string tr;
	if (ops >> Option(' ', std::string("trackedMotif"), tr)) {
	  printf("tracking all submotifs of %s\n", tr.c_str());
	  trackedMotifs.insert(tr);
	}
	if (aa && (ops >> Option(' ', std::string("extensionECut"), extensionECut, 1.0) || DEBUG)) printf("E cutoff for seed extension: %f\n", extensionECut);
	(aa && ops >> Option(' ', std::string("extensionMinCut"), extensionMinCut, 1000));
	if (aa && DEBUG) printf("min number of seeds extended: %d\n", extensionMinCut);
	if (aa && (ops >> Option(' ', std::string("extensionMaxCut"), extensionMaxCut, std::numeric_limits<int>::max()) || DEBUG)) {
		printf("max number of seeds extended: %d\n", extensionMaxCut);
		if (extensionMinCut>extensionMaxCut) {
			std::cerr << "ERROR: Minimum number of extensions is greater than maximum number." << std::endl;
			exit(1);
		}
	}
	if (aa && ops >> Option(' ', std::string("aaStateSigThresh"), aaStateSigThresh)) {
		if (aaStateSigThresh<=0) {
			std::cerr << "ERROR: State significance threshold must be positive." << std::endl;
			exit(1);
		}
	}
	aaStateSigThresh = log(aaStateSigThresh);
	if (aa && DEBUG) printf("thresh for sig states: %f = ln(%.2f)\n", log(aaStateSigThresh), aaStateSigThresh);

	if (aa && (ops >> Option(' ', std::string("aaSeqFreqThresh"), aaSeqFreqThresh) || DEBUG)) {
		printf("min fraction of conserved positions for extension: %f\n", aaSeqFreqThresh);
		if (aaSeqFreqThresh<=0 || aaSeqFreqThresh>1) {
			std::cerr << "ERROR: fraction of conserved positions must be in ]0,1]." << std::endl;
			exit(1);
		}
	}

	if (ops >> OptionPresent(' ', std::string("batch"), batch_mode)){
		if (batch_mode && DEBUG) printf("running in batch mode (no progress bars)\n");
	}
	bool no_pwm_opt = false;
	if (ops >> OptionPresent(' ', std::string("no-pwm-length-optimization"), no_pwm_opt)) {
		if (no_pwm_opt) {
			maximizeMotifLength = false;
			if (DEBUG) printf("PWM length optimization disabled\n");
		}
	}

	/* fail if any options remain unprocessed */
	if (ops.options_remain()) {
		std::list<std::string> unp = ops.remaining_options();
		for (std::list<std::string>::const_iterator it = unp.begin(); it!=unp.end(); ++it) {
			std::cerr << "ERROR: Unknown option " << *it << std::endl;
		}
		return false;
	}

	return true;
}

Global::~Global(){
   	if(name != NULL)  free(name);
   	if(outputDirectory != NULL) free(outputDirectory);
   	if(shortFileName != NULL) free(shortFileName);
   	if(negFile != NULL) free(negFile);
   	if(startMotif != NULL) free(startMotif);
   	if(profFile != NULL) free(profFile);
   	if(benchmarkFolder != NULL) free(benchmarkFolder);
   	if(pwmFolder != NULL) free(pwmFolder);

	ss_type set = posSet;
	if(negSet != NULL) set = negSet;

	double POW_2_16 = pow(2,16);
	if(conservationProbs != NULL){
		for(int i=0; i<= POW_2_16; i++){
			if(conservationProbs[i] == NULL) continue;
			for(int j=1; j< set->max_MultSeq; j++){
				free(conservationProbs[i][j]);
			}
			free(conservationProbs[i]);
		}
		free(conservationProbs);
	}

	if(alignmentFreeProbs != NULL){
		for(int i=0; i<= POW_2_16; i++){
			if(alignmentFreeProbs[i] == NULL) continue;
			for(int j=1; j< set->max_MultSeq; j++){
				free(alignmentFreeProbs[i][j]);
			}
			free(alignmentFreeProbs[i]);
		}
		free(alignmentFreeProbs);
	}

	free(posBg_log);
   	free(negBg_log);
   	free(posBg);
   	free(negBg);

	NullModel::destruct();

   	if (useCompositionProbs) NullModelCompositional::destruct();

   	if( em )
   		hoNullModel::destruct();

   	if( cv )
   		if( shortTestFileName != NULL )
   			free( shortTestFileName );

   	if(negSet != NULL) NilSeqSet(negSet);

   	NilSeqSet(posSet);
   	NilAlpha(A);
}

std::ostream& operator<<(std::ostream &os, const merge_type &m) {
  switch(m) {
    case LOW: os << "LOW"; break;
    case MEDIUM: os << "MEDIUM"; break;
    case HIGH: os << "HIGH"; break;
    case NO_VALID_MERGE_MODE: os << "NO_VALID_MERGE_MODE"; break;
    default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit(1);
  }
  return os;
}

std::ostream& operator<<(std::ostream &os, const SuppInfMode &v) {
	switch (v) {
	case SUPP_NO: os << "NO"; break;
	case SUPP_DISOCONS: os << "DISOCONS"; break;
	case SUPP_NNET: os << "NNET"; break;
	default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit(1);
	}
	return os;
}

