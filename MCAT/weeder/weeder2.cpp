/* 
***************************************************************
This program is free software. It comes without any warranty, to
the extent permitted by applicable law. You can redistribute it
and/or modify it under the terms of the 

GNU GENERAL PUBLIC LICENSE (GPL) Version 3 

that you can find at:

http://159.149.160.51/modtools/downloads/weederlicense.txt


Authors: 

Giulio Pavesi - giulio.pavesi at unimi.it
Federico Zambelli - federico.zambelli at unimi.it
***************************************************************
*/

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <cstdlib>

//#define DEBUG1
//#define DEBUG3
//#define DEBUG4
//#define DEBUG5
//#define DEBUG6

using namespace std;

const char ALPHABET[4] = {'A','C','G','T'}, RC_ALPHABET[4] = {'T','G','C','A'}; //CONSTANTS DEFINITIONS
const short int ALPHABET_SIZE = 4, OLIGO_TYPES = 3, MAX_AVAILABLE_FREQS = 10;
const short int OLIGO_LENGTHS[OLIGO_TYPES] = {6,8,10}, OLIGO_MAX_MIS[OLIGO_TYPES] = {1,2,3};
const int GC_DIV = 1;
enum {SYNTAX_ERROR, MISSING_FASTA_FILE, BAD_FASTA_SEQUENCE, BAD_SCANNING_OLIGO, NO_FREQ_FILE, TOPN_NOT_POSITIVE, CH_NOT_ALLOWED, /*EM_NEGATIVE,*/ BAD_SIM_THRESHOLD};
enum {Seq, Pos};

enum {TTEST,HYPER,ZSCORE,BERNSTEIN};
const short int SCORE_TYPE = 4, DEFAULT_MAX_SEQ = 100;
unsigned int MAX_SEQ = -1;
bool CHIPSEQ = false;

class oli
{
	private:
	friend class sequence;
	friend class matrix;
	bool max_lev;
//	vector<vector <map<string,bool>::iterator> > o_nhr;
//	vector<vector<vector <char> > >u_pos;
	vector<float> *ofreq;
	short int Type;
	int mis_level;
//	void raise_one_level();
//	void give_level(int);
//	double get_freq();
//	double get_freq_without_freqs();
//	bool pal_check(string*);
	public:
	oli(string, int);
	#ifdef DEBUG4
	string see_level(int);
	#endif
};


class sequence 
{
        private:
	friend void compute_scores(vector<sequence>*, vector< map<string,vector<float> >  >*, map<string,oli>*,int,bool); 
	bool good_flag;
        string NAME;
	string OR_SEQ;
	string SEQ;
	float GC_COUNT;
        unsigned short int GC;
	string mask_seq;
	int used_map, used_length;
	map<string,oli> *OLI_MAP;
	vector< map<string,double> > scores;
	vector<bool> used;
        bool seq_cleaner();
	void set_used_map(int);
	void oligo_scan(string);
	double get_score(oli*,int,int);
	void masker(string*,int);
	bool oligo_check(string*,int,int);
        public:
	sequence(string, map<string,oli>*);
//	void o_Clean(vector<sequence>*);
	void mask_sequence(vector<string>*,int,map<string,double> *);
	void used_reset();
 //     bool assign(string, map<string,oli>*);
	void scan(vector<sequence>*,int,int);
//	void calculate_scores(map<string,vector<float> >*,map<string,oli>*,int);
	void mat_oligo_scan(string,vector<string>*);
        string name();
        string seq();
	string or_seq();
	bool good();
	unsigned short int gc();
	#ifdef DEBUG5
	string see_scores(int);
	#endif
	#ifdef DEBUG1
	string check_1(); //OUTPUT SEQUENCE OLIGO COUNTS
	#endif
};

class matrix
{
	friend void write_output(vector<sequence>*,vector<matrix>*, int, char**);
	private:
	unsigned int mpos;
	bool not_too_similar;  
	vector<vector <double> > Z;
	vector<vector <bool> > zrc;
	int **best_occ;
	double *best_occ_score;
	string sstart;
	double NSITES;
	double sscore;
	double min,max;
	int x_s, y_s, Type;
	double **s_matrix, **n_matrix;
	bool good_matrix;
	vector<sequence> *s_ptr;
	void normalize();
	void get_z();
	void get_p();
	void squash();
	void set_min_max();
//	void generate_pfm(double**);
	void scan_Z();
	bool similarity_check(vector<matrix>*);
	public:
	double relative_score(string*);
	double get_n_element(unsigned int, unsigned int);
	double get_s_element(unsigned int, unsigned int);
	matrix(string, vector<sequence>*, vector<matrix>*, unsigned int);
	string get_sstart();
//	void write_tompa_format();
	bool sim_check();
	unsigned int motif_pos();
	unsigned int size();
	double nsites();
	double get_min();
	double get_max();
//	void see_matrix(double**);
};

string inputfile = "", ORGANISM = "HS",ofile, omfile/*, mlist = "mlist.txt", pscan_background_file*/; //GLOBAL VAR DECLARATIONS
map <string,bool> *mis_map;
map<string,double> MAIN_OLIGO_MAP;
vector< map<string, vector<float> > > *OF;
map<string,double> *SC_PTR;
unsigned short int TOPN = 50, EM_CYCLES = 1, /*M_NUM = 1*/ MAX_MATRIX = 25;
bool no_freq = false, DOUBLE_STRAND = true, raw_matrix = false/*, parse_only = false*/;
float gcdiv[GC_DIV+1], SIM_THRESHOLD = 0.95;

void command_line_parser(int, char**); // FUNCTION DECLARATIONS
vector<string> fasta_reader(const char*);
void error_handler(int, string);
void main_scan(vector<sequence>*,int);
void compute_scores(vector<sequence>*, vector< map<string,vector<float> > >*, map<string,oli>*,int,bool);
void get_freq(int, map<string, vector<float> >*);
void sort(map<string,double>, vector<string>*);
//float infer_freq(const string*, map<string,float>*);
short int get_type_from_l(short int);
string string_rev_comp(const string*);
double avg_seq_length(vector<sequence>*);
int ati(char,string*);
//void run_pscan();
unsigned int rc_index(unsigned int);
//string pscan_parser(short int);
void write_output(vector<sequence>*,vector<matrix>*, int, char**);
int choose_matrix(string,vector<matrix>*);
char bool_to_strand(bool);
//void parse_to_tompa();
void dislay_help();

int main (int argc, char **argv)
{
	vector<string> seqs;
	vector<string> Good_motifs;
	vector<sequence> SEQUENCE;
	vector<matrix> MATRIX;
	vector< map<string, vector<float> > > oligo_freqs (OLIGO_TYPES, map<string,vector<float> >());
	map<string, double> BEST_SCORES;
	map<string,oli> oli_map;
	map<string,bool> Mis_map;

	//**********PROV GC DINUCLEOTIDE THRESHOLDS FOR HUMAN
//	gcdiv[0] = 0;
//	gcdiv[1] = 0.070;
//	gcdiv[2] = 0.140;
//	gcdiv[GC_DIV] = 1;
	//*************************************

	//**********PROV GC THRESHOLDS FOR HUMAN
        gcdiv[0] = 0;
    //    gcdiv[1] = 0.455;
      //  gcdiv[2] = 0.595;
        gcdiv[GC_DIV] = 1;
        //*************************************

	mis_map = &Mis_map;
	OF = &oligo_freqs;
	SC_PTR = &BEST_SCORES;

	command_line_parser(argc, argv);

    /** EDITED file paths by Jeff Robertson to work for CS4884 pipeline 3/04/16 (last updated 09-23-16)**/
	ofile = inputfile;
	omfile = inputfile;

    int iDirLen = inputfile.length() - 1;

    while (inputfile[iDirLen] != '/') {
        iDirLen--;
    }

    ofile += ".w2";
    ofile = "results" + ofile.substr(iDirLen);
	omfile += ".matrix.w2";
    omfile = "results" + omfile.substr(iDirLen);

//	if(parse_only)
//		parse_to_tompa();

	seqs = fasta_reader(inputfile.c_str());

	ofstream out(ofile.c_str());
	out.close();

	for(int t = 0; t < (int)seqs.size(); t++)
        {
                sequence one_more_sequence(seqs[t],&oli_map);

                if(!one_more_sequence.good())
                        error_handler(BAD_FASTA_SEQUENCE,seqs[t]);
                else
                        SEQUENCE.push_back(one_more_sequence);
        }

	seqs.clear();

	for(int i = 0; i < OLIGO_TYPES; i++)
		if(OLIGO_LENGTHS[i] <= MAX_AVAILABLE_FREQS)
	//		oligo_freqs.push_back(get_freq(OLIGO_LENGTHS[i]));
			get_freq(i, &oligo_freqs[i]);

	for(int i = 0; i < OLIGO_TYPES; i++)
	{
		oli_map.clear();
		Mis_map.clear();		
		MAIN_OLIGO_MAP.clear();

		cerr << "\nScanning " << OLIGO_LENGTHS[i] << "MERS...";
		main_scan(&SEQUENCE,i);
		cerr << "\tdone" << endl;

/*		if(i >= oligo_freqs.size())
		{
			map<string,vector<float> > nf;	
			oligo_freqs.push_back(nf);
			no_freq = true;

		}*/


		cerr << "Scoring " << OLIGO_LENGTHS[i] << "MERS...";
		compute_scores(&SEQUENCE, &oligo_freqs, &oli_map, i, no_freq);
		cerr << "\tdone";
	}

	cerr << "\nSorting motifs and building matrices...";
	oli_map.clear();
	sort(BEST_SCORES,&Good_motifs);

//	for(int i = 0; i < Good_motifs.size(); i++) //&& i < MAX_MATRIX; i++)
		
	unsigned int gmcount = 0;	

	for(int i = 0; i < Good_motifs.size() && gmcount < MAX_MATRIX; i++)
	{
		matrix nm = matrix(Good_motifs.at(i), &SEQUENCE, &MATRIX, i);

		if(nm.sim_check())
		{
			MATRIX.push_back(nm);
			gmcount++;
		}

//		MATRIX.push_back(matrix(Good_motifs.at(i), &SEQUENCE, &MATRIX));	
	}

	cerr << "\tdone";

	cerr << "\nWriting output in " << ofile << " and " << omfile;
	write_output(&SEQUENCE,&MATRIX, argc, argv);
	cerr << "\tdone" << endl;


/*	if(pscan_background_file.size() != 0)
	{
		string best_pscan;

		cerr << endl << "Running pscan..." << endl;
		run_pscan();

		best_pscan = pscan_parser(BERNSTEIN);

//		cerr << "\nbest_pscan = " << best_pscan << endl;

//		MATRIX[choose_matrix(best_pscan, &MATRIX)].write_tompa_format();
	}
*/
	#ifdef DEBUG5
	for(int u = 0; u < OLIGO_TYPES; u++)
		for(int i = 0; i < SEQUENCE.size(); i++)
			cout << SEQUENCE[i].name()  << endl << SEQUENCE[i].see_scores(u) << endl;
	#endif

	#ifdef DEBUG1
	for(int i = 0; i < SEQUENCE.size(); i++)
		cout << SEQUENCE[i].check_1();
	#endif
	
	return EXIT_SUCCESS;
}

unsigned int rc_index(unsigned int i)
{
	unsigned int i_r;

	switch(i)
	{
		case 0:
		i_r = 3;
		break;

		case 1:
		i_r = 2;
                break;

		case 2:
		i_r = 1;
                break;

		case 3:
		i_r = 0;
                break;	
	}

	return i_r;
}

void main_scan(vector<sequence> *SEQUENCE, int i)
{
	cerr << endl;

	for(int j = 0; j < SEQUENCE->size() && j < MAX_SEQ; j++)
		SEQUENCE->at(j).scan(SEQUENCE,OLIGO_LENGTHS[i], j);	
	

	for(int j = 0; j < SEQUENCE->size() && j < MAX_SEQ; j++)
		SEQUENCE->at(j).used_reset();

	return;
}

double avg_seq_length(vector<sequence> *SEQUENCE)
{
	double tot_length = 0;

	for(int i = 0; i < SEQUENCE->size(); i++)
		tot_length += SEQUENCE->at(i).seq().size();

	return tot_length/SEQUENCE->size();
}

void compute_scores(vector<sequence> *SEQUENCE, vector< map<string,vector<float> > >*freqs, map<string,oli>* oli_map, int u, bool no_freq)
{
//	map<string, double> SCORES = MAIN_OLIGO_MAP;
	vector<string> s;

	double avg_seq_l = avg_seq_length(SEQUENCE);

	for(map<string,double>::iterator mi = MAIN_OLIGO_MAP.begin(); mi != MAIN_OLIGO_MAP.end(); mi++)
	{
		mi->second -= (freqs->at(u).find(mi->first)->second[0] * avg_seq_l);
		mi->second /= ((SEQUENCE->size()) - 1);	
	}
	
	sort(MAIN_OLIGO_MAP, &s);

	#ifdef DEBUG6
	cout << endl << endl << OLIGO_LENGTHS[u] << "MERS" << endl;
	#endif

	for(int i = 0; i < SEQUENCE->size(); i++)
	{
		SEQUENCE->at(i).mask_sequence(&s, u, &MAIN_OLIGO_MAP);
		#ifdef DEBUG6
		cout << ">" << i << endl;
		cout << SEQUENCE->at(i).seq() << endl;
		#endif
	}

	for(int j = 0; j < s.size() && j < TOPN; j++)
		SC_PTR->insert(make_pair(s[j],MAIN_OLIGO_MAP[s[j]]));

	return;
}

void error_handler(int error, string s_err)
{
        cerr << endl;

        switch(error)
        {
                case SYNTAX_ERROR:
                {
                        cerr << "\nSyntax: weeder2 -f filename [OPTIONS]\n" << endl;
                        break;
                }
                case MISSING_FASTA_FILE:
                {
                        cerr << "\nMmmh.. i was sure " << s_err << " was here just a moment ago...\n";
                        break;
                }
                case BAD_FASTA_SEQUENCE:
                {
                        cerr << "\nI think this sequence has some problem, i'm skipping it...\n" << s_err << endl;
                        return;
                        break;
                }
		case BAD_SCANNING_OLIGO:
                {
                        cerr << "\nFatal error, scanning with an oligo of wrong size... \n";
                        break;
                }
		case NO_FREQ_FILE:
		{
			cerr << "\nMissing frequency file: " << s_err << endl;
			break;
		}
/*		case TOPN_NOT_POSITIVE:
		{
			cerr << "\n\"Used top scoring oligos\" value must be greater than 0, using default..." << endl;
			return;
			break;
		}*/

		case CH_NOT_ALLOWED:
		{
			cerr << "\nFatal error: character not allowed found: " << s_err << endl;
			break;
		}

/*		case EM_NEGATIVE:
		{
			cerr << "\nEM cycles value must be >= 0. Using default..." << endl;
			return;
			break;
		}*/
		case BAD_SIM_THRESHOLD:
                {
                        cerr << "\nSimilarity threshold must be 0<= sim <= 1. Using default..." << endl;
                        return;
                        break;
                }
                default:
                {
                        cerr << "\nSome weird error occurred...\n";
                        break;
                }
        }

        cerr << endl;

        exit(error);
}

vector<string> fasta_reader(const char *f_name)
{
        ifstream in(f_name);
        string line,buf;
        vector<string> seqs;
        bool flag = false;

        if(!in)
        {
                string tmp = f_name;
                error_handler(MISSING_FASTA_FILE, tmp);
        }

        while(getline(in,line))
        {
                if(line[0] == '>' && flag)
                {
                        seqs.push_back(buf);
                        buf.clear();
                        buf += line;
                        buf += '\n';
                }

                else if(line[0] == '>' && !flag)
                {
                        buf += line;
                        buf += '\n';
                        flag = true;
                }

                else
                {
                        buf += line;
                        buf += '\n';
                }
        }

        seqs.push_back(buf);
        buf.clear();

	
        in.close();

        return seqs;
}

void display_help()
{
	cerr 	<< "\nSYNTAX\n\nweeder2 -f input_file [-O frequency_file_organism_code] [options]" << endl << endl
		<< "input_file must be in Multi-FASTA format." << endl 
		<< "\nWhen no organism code for oligo frequency files is provided it is assumed to be HS (Homo sapiens)." << endl
		<< "\nOligo frequency files for the following organisms are available in the standard Weeder 2.0 package: " << endl
		<< "\nHomo sapiens - Code: HS" << endl
		<< "Mus musculus - Code: MM" << endl
		<< "Drosophila melanogaster - Code: DM" << endl
		<< "Saccharomyces cerevisiae - Code: SC" << endl
		<< "Arabidopsis thaliana - Code: AT" << endl
		<< "\nOther frequency files  may be added to the FreqFiles directory by using the \"Frequency maker\" program" << endl
		<< "available at http://www.beaconlab.it/modtools" << endl << endl
		<< "OPTIONS" << endl << endl
		<< "-chipseq" << endl
		<< "This flag activates the ChIP-Seq heuristic speeding-up the computation." << endl << endl
		<< "-top <num> (DEFAULT: " << DEFAULT_MAX_SEQ << ")" << endl 
		<< "If the -chipseq parameter is used Weeder 2.0 scans all the input sequences for occurrences of the oligos contained in the top <num> input sequences." << endl
		<< "Increase this value when your input has many more than <num> sequences to improve the chance of finding motifs enriched only " << endl
		<< "in a subset of your input sequences." << endl << endl
		<< "-b <num> (DEFAULT: " << TOPN << ")" << endl
		<< "Weeder 2.0 builds occurrences matrix profiles and outputs (if other conditions are met) only the top <num> scoring motifs" << endl
		<< "for each motif length. Increase this value to have more (lowest scoring) motifs in the output (see also -maxm)." << endl << endl
		<< "-maxm <num> (DEFAULT: " << MAX_MATRIX << ")" << endl
		<< "To limit the output length, Weeder 2.0 reports only the top <num> scoring motifs with their associated occurrences" << endl
		<< "matrix and occurrences list. Increase <num> to have longer outputs with more lowest scoring motifs." << endl << endl
		<< "-ss" << endl
		<< "Single strand mode." << endl << endl
		<< "ADVANCED OPTIONS" << endl << endl
		<< "-sim <num> (DEFAULT: " << SIM_THRESHOLD << " MIN: 0 MAX: 1)" << endl
		<< "Similarity threshold for the redundancy filter. This filter removes from the output those motifs that are too similar to other motifs" << endl
		<< "already reported. Values close to 0 mean a stricter filter and vice versa values close to 1 impose a looser filter." << endl
		<< "Set <num> to 1 to disable the filter altogether. Set it to 0 to have in the output only the top scoring oligo for each one of" << endl
		<< "the possible oligo lengths (6, 8 and 10)." << endl << endl 
		<< "-em <num> (DEFAULT: " << EM_CYCLES << " MIN: 0 MAX: 100)" << endl	
		<< "Weeder 2.0 has a built-in expectation maximization (EM) matrix profiles refinement step." << endl
		<< "<num> defines the number of EM cycles to be performed by Weeder 2.0." << endl
		<< "One (default) or few EM cycles should be sufficient to \"clean\" matrix profiles without overfitting them." << endl << endl;

	return;
}

void command_line_parser(int argc, char **argv)
{
	if(argc == 1)
	{
		display_help();
		exit(EXIT_SUCCESS);
	}

        if(argc < 2)
                error_handler(SYNTAX_ERROR,"");

        for(int i = 1; i < argc; i++)
        {
                string buf = argv[i];

                if(buf == "-f")
                {
                        if(i + 1 < argc)
                                inputfile = argv[++i];
                        continue;
                }
		else if(buf == "-top")
		{
			if(i + 1 < argc)
				MAX_SEQ = atoi(argv[++i]);
			continue;
		}
		else if(buf == "-maxm")
		{
			if(i + 1 < argc)
				MAX_MATRIX = atoi(argv[++i]);
			continue;
		}
		else if(buf == "-b")
		{
			if(i + 1 < argc)
			{
				string s = argv[++i];
				istringstream sr(s);	

				sr >> TOPN;

				if(TOPN <= 0)
				{
					error_handler(TOPN_NOT_POSITIVE,"");
					TOPN = 50;
				}
			}
		}

		else if(buf == "-em")
                {
                        if(i + 1 < argc)
                        {
                                string s = argv[++i];
                                istringstream sr(s);

                                sr >> EM_CYCLES;

                                if(EM_CYCLES > 100)
					EM_CYCLES = 100;
                        }
                }
		else if(buf == "-sim")
                {
                        if(i + 1 < argc)
                        {
                                string s = argv[++i];
                                istringstream sr(s);

                                sr >> SIM_THRESHOLD;

                                if(SIM_THRESHOLD < 0 || SIM_THRESHOLD > 1)
                                {
                                        error_handler(BAD_SIM_THRESHOLD,sr.str());
                                        SIM_THRESHOLD = 0.95;
                                }
                        }
                }

/*		else if(buf == "-pscan")
                {
                        if(i < argc)
                        {
                                string s = argv[++i];
                                istringstream sr(s);

                                sr >> pscan_background_file;
                        }
                }*/
		else if(buf == "-O")
		{
                        if(i + 1 < argc)
                                ORGANISM = argv[++i];
                        continue;
                }

		else if(buf == "-ss")
		{
			DOUBLE_STRAND = false;
			continue;
		}

/*		else if(buf == "-parse")
                {
                        parse_only = true;
                        continue;
                }
*/		
		else if(buf == "-raw")
                {
                        raw_matrix = true;
                        continue;
                }
		else if(buf == "-chipseq")
		{
			CHIPSEQ = true;
			continue;
		}
		else if(buf == "-h" || buf == "--help")
		{
			display_help();
			exit(EXIT_SUCCESS);
		}

                else
                {
                        cerr << "\nBad argument: " << buf << endl << endl;
                        exit(EXIT_FAILURE);
                }
        }

	if(inputfile.empty())
		error_handler(SYNTAX_ERROR,"");

	if(!CHIPSEQ)
		MAX_SEQ = -1;
	else
	{
		if(MAX_SEQ == -1)
			MAX_SEQ = DEFAULT_MAX_SEQ;
	}

        return;
}

void sort(map<string,double> scores, vector<string> *OLIGO)
{
	double INF = -1000000000;
	string oligo;

	for(map<string,double>::iterator mii = scores.begin(); mii != scores.end(); mii++)
	{
		double inf = INF;
		for(map<string,double>::iterator mi = scores.begin(); mi != scores.end(); mi++)
		{
			if(mi->second > inf)
			{
				inf = mi->second;
				oligo = mi->first;
			}

		}		

		if(scores[oligo] > 0)
			OLIGO->push_back(oligo);
	
		scores[oligo] = INF;
	}

	return;
}
		
void get_freq(int t, map<string, vector<float> > *freqs)
{
	ostringstream filename;
//	map<string,float> freqs;
	string line;
    
    /** EDITED file paths by Jeff Robertson to work for CS4884 pipeline 3/04/16 **/
	if(DOUBLE_STRAND)
		filename << "weeder/FreqFiles/" << ORGANISM << "_ds" << "." << OLIGO_LENGTHS[t] << ".freq";
	else
		filename << "weeder/FreqFiles/" << ORGANISM << "_ss" << "." << OLIGO_LENGTHS[t] << ".freq";
		
	ifstream in(filename.str().c_str());

	cerr << "\nopening " << filename.str();

	if(!in) {
        filename.str("");
        if(DOUBLE_STRAND)
            filename << "./FreqFiles/" << ORGANISM << "_ds" << "." << OLIGO_LENGTHS[t] << ".freq";
        else
            filename << "./FreqFiles/" << ORGANISM << "_ss" << "." << OLIGO_LENGTHS[t] << ".freq";

        ifstream in(filename.str().c_str());

        cerr << "\nopening " << filename.str();
	    if(!in)
		    error_handler(NO_FREQ_FILE,filename.str());
    }

	while(getline(in,line))
	{
		string motif;
		float occ;

		if(line.size())
		{
			istringstream s1(line);

			s1 >> motif;

//			freqs->insert(make_pair(motif,vector<float>()));
//			freqs->insert(make_pair(motif,vector<float>((OLIGO_MAX_MIS[t] + 1)*GC_DIV,0)));

			vector<float> vf((OLIGO_MAX_MIS[t] + 1)*GC_DIV,0);

			for(unsigned short int a = 0; a < (OLIGO_MAX_MIS[t] + 1) * GC_DIV; a++)
			{
				s1 >> occ;

				vf[a] = occ;
//				freqs->find(motif)->second.push_back(occ);
				
			}

			freqs->insert(make_pair(motif,vf));
		}
	}

	in.close();

	cerr << "\nclosing " << filename.str();

	for(unsigned int q = 0; q < GC_DIV; q++)
	{

		float sum = 0;

		for(map<string,vector<float> >::iterator mi = freqs->begin(); mi != freqs->end(); mi++)
			sum += mi->second[(OLIGO_MAX_MIS[t] + 1) * q];

		for(map<string,vector<float> >::iterator mi = freqs->begin(); mi != freqs->end(); mi++)
			for(unsigned short int a = 0; a < OLIGO_MAX_MIS[t] + 1; a++)
				mi->second[a + ((OLIGO_MAX_MIS[t] + 1) * q)] /= sum;	
	}

	#ifdef DEBUG3
	for(map<string,vector<float> >::iterator mi = freqs->begin(); mi != freqs->end(); mi++)
	{
		cout << mi->first;

		for(unsigned short int a = 0; a < OLIGO_MAX_MIS[t] + 1; a++)
			cout  << '\t' << mi->second[a];

		cout << endl;
	}
	#endif

	
	return;
}

short int get_type_from_l(short int l)
{
	for(short int i = 0; i < OLIGO_TYPES; i++)
		if(OLIGO_LENGTHS[i] == l)
			return i;
}

/*float infer_freq(const string *mi, map<string,float> *freq)
{
	int osize = freq->begin()->first.size(), 
            diff = (int)mi->size() - osize;	
	float f = 0;

	f = freq->find(mi->substr(0,osize))->second; 

	for(int i = 1; i <= diff; i++)
	{
		float f1, div = 0;

		string ov = mi->substr(i,osize);

		f1 = freq->find(ov)->second;	

		string::iterator si = ov.end();
		si--;

		for(int j = 0; j < ALPHABET_SIZE; j++)
		{
			*si = ALPHABET[j];	
			div += freq->find(ov)->second;
		}

		f1 /= div;

		f *= f1;
	}
			
	return f;
}
*/
string string_rev_comp(const string *seq)
{
	int size = (int)seq->size();

	string rc_seq(size,'N');

	for(int i = 0; i < size; i++)
	{
		for(int j = 0; j < ALPHABET_SIZE; j++)
		{
			if(seq->at(i) == ALPHABET[j])
			{
				rc_seq[size - 1 - i] = RC_ALPHABET[j];
				break;
			}
		}
	}

	return rc_seq;
}

int ati(char c, string *err)
{
        for(int i = 0; i < ALPHABET_SIZE; i++)
                if(ALPHABET[i] == c)
                        return i;

//	string err;

//	err.push_back(c);

	error_handler(CH_NOT_ALLOWED, *err);
}

/*void run_pscan()
{
	ostringstream com;

	com << "pscan -q " << inputfile  << " -p " << pscan_background_file << " -w2 " << ofile;

	if(!DOUBLE_STRAND)
		com << " -ss";

	cerr << endl << com.str() << endl;

	system(com.str().c_str());

	return;
}*/

/*string pscan_parser(short int method)
{
	ostringstream resfile;
	string line;
	int pos;
	vector<vector<double> > scores(SCORE_TYPE);
	vector<string> name;
	double max = 1000000000;

	resfile << inputfile << ".res";

	ifstream in(resfile.str().c_str());

	while(getline(in,line))
	{
		double score;
		istringstream ss(line);
		string trash;

		ss >> trash;

		name.push_back(trash);

		ss  >> trash >> score;

		scores[TTEST].push_back(score);

		ss >> score;

		scores[HYPER].push_back(score);

		ss >> trash >> trash >> score;

		scores[ZSCORE].push_back(score);

		ss >> score;

		scores[BERNSTEIN].push_back(score);
	}
	
	for(int i = 0; i < scores[method].size(); i++)
	{
		if(scores[method][i] < max && name[i].size() == OLIGO_LENGTHS[OLIGO_TYPES - 1]) 
		{
			max = scores[method][i];
			pos = i;
		}
	}

	return name[pos];
}*/

int choose_matrix(string motif, vector<matrix> *MATRIX)
{
	for(int i = 0; i < MATRIX->size(); i++)
		if(motif == MATRIX->at(i).get_sstart())
			return i;
}
		
void write_output(vector<sequence> *SEQUENCE, vector<matrix> *MATRIX, int argc, char **argv)
{
	ofstream out(ofile.c_str());
	ofstream outm(omfile.c_str());
	out.precision(4);
	outm.precision(4);

	out << "COMMAND LINE:\n\n";

	for(unsigned int i = 0; i < argc; i++)
		out << argv[i] << ' ';

	out << endl << endl;

	out << "MOTIFS SUMMARY:" << endl << endl;

	for(int i = 0; i < MATRIX->size(); i++)
	{
		out << i+1 << ')' << '\t' << MATRIX->at(i).sstart;	

		if(DOUBLE_STRAND)
			out << '\t' << '(' << string_rev_comp(&MATRIX->at(i).sstart) << ')';

		out << '\t' << MATRIX->at(i).sscore << endl;
	}
		
	out << endl <<  endl << "DETAILED RESULTS:" << endl << endl;

	for(int i = 0; i < MATRIX->size(); i++)
	{
		multimap<double, string> MSORTER;

//		if(MATRIX->at(i).sstart.size() < 10) //PRINT ONLY DECAMERS
//			continue; //PRINT ONLY DECAMERS

		out << i+1 << ')' << '\t' << MATRIX->at(i).sstart;

		if(DOUBLE_STRAND)
			out << '\t' << '(' << string_rev_comp(&MATRIX->at(i).sstart) << ')';

		out << '\t' << MATRIX->at(i).sscore << endl << endl;	

		if(!MATRIX->at(i).good_matrix)
		{
			out << "BAD MATRIX" << endl << "***********\n" << endl << endl;
                        continue;
		}
		
		out << "Matrix: MAT" << i+1 << '\t' << MATRIX->at(i).sstart  /*<< '\t' << "Max = " << MATRIX->at(i).get_max() << '\t' << "Min = " << MATRIX->at(i).get_min()*/ << endl;
		outm << ">MAT" << i+1 << '\t' << MATRIX->at(i).sstart << endl;

		for(int k = 0; k < MATRIX->at(i).y_s; k++)
		{
			out << ALPHABET[k] << '\t';
			outm << ALPHABET[k] << '\t';

			if(!raw_matrix)
				for(int l = 0; l < MATRIX->at(i).x_s; l++)
				{
					out << MATRIX->at(i).get_n_element(l,k) << '\t';
					outm << MATRIX->at(i).get_n_element(l,k) << '\t';
				}
			else
				for(int l = 0; l < MATRIX->at(i).x_s; l++)
				{
					out << MATRIX->at(i).get_s_element(l,k) << '\t';
					outm << MATRIX->at(i).get_s_element(l,k) << '\t';
				}
			
			out << endl;
			outm << endl;
		} 

		out << "\nOCCURRENCES:" << endl;

		for(int j = 0; j < SEQUENCE->size(); j++)
		{
			ostringstream outbuf;

			string s = SEQUENCE->at(MATRIX->at(i).best_occ[j][Seq]).seq().substr(MATRIX->at(i).best_occ[j][Pos],OLIGO_LENGTHS[MATRIX->at(i).Type]);

			if(MATRIX->at(i).zrc[MATRIX->at(i).best_occ[j][Seq]][MATRIX->at(i).best_occ[j][Pos]])
				s = string_rev_comp(&s);

			double rscore = MATRIX->at(i).relative_score(&s);

			outbuf << SEQUENCE->at(MATRIX->at(i).best_occ[j][Seq]).name().substr(0,20) << '\t' << MATRIX->at(i).best_occ[j][Seq] + 1 << '\t' 
			    << s  /*<< '\t' << MATRIX->at(i).best_occ_score[j]*/  << '\t' << rscore
                            << '\t' << /*-(int)SEQUENCE->at(MATRIX->at(i).best_occ[j][Seq]).seq().size() +*/  MATRIX->at(i).best_occ[j][Pos] + 1
                            << '\t' << bool_to_strand(MATRIX->at(i).zrc[MATRIX->at(i).best_occ[j][Seq]][MATRIX->at(i).best_occ[j][Pos]]) << endl;

			MSORTER.insert(make_pair(rscore, outbuf.str()));	
		}

		multimap<double, string>::reverse_iterator mi = MSORTER.rbegin();

		while(mi != MSORTER.rend())
		{
			out << mi->second;
			mi++;
		}

		out << "**********\n" << endl;
	}

	out.close();
	outm.close();

	return;
			
}

char bool_to_strand(bool s)
{
	if(s)
		return '-';
		
	return '+';
}

/*void parse_to_tompa()
{
	ifstream in(ofile.c_str());

	string line;

	cout << ">dataset\n";

        cout << inputfile.substr(0,inputfile.find(".")) << endl;

        cout << ">instances\n";

	while(getline(in,line))
	{
		istringstream ss(line);
		string buf;

		ss >> buf;

		if(buf.size() == OLIGO_LENGTHS[OLIGO_TYPES - 1])
		{	
			getline(in,line);

			while(getline(in,line) && line.size())
			{
				string seq,trash,motif,pos;
			
				istringstream ss2(line);

				ss2 >> seq >> motif >> trash >> pos >> trash;

				if(trash[0] == '-')
					motif = string_rev_comp(&motif);

				cout << seq << ',' << pos << ',' << motif << endl;
			}

			break;
		}
	}

	exit(EXIT_SUCCESS);
}
*/
sequence::sequence(string s, map<string, oli>* oli_map)
{
        istringstream str(s);
        string line;

        getline(str,line);

        NAME = line;

	OLI_MAP = oli_map;

        while(getline(str,line))
                SEQ += line;

	for(int i = 0; i < OLIGO_TYPES; i++)
	{
		map<string,double> m;
		scores.push_back(m);
	}

	good_flag = seq_cleaner();

	float gcpc = GC_COUNT / (float)SEQ.size();
        
        for(unsigned short int i = 0; i < GC_DIV; i++)
        {
                if(gcpc >= gcdiv[i] && gcpc < gcdiv[i+1])
                {
                        GC = i;
                        break;
                }
        }

        return;
}

bool sequence::good()
{
	return good_flag;
}

void sequence::masker(string *s, int mmis)
{
	string::iterator si = mask_seq.begin();

	for(int i = 0; i < SEQ.size() - s->size(); i++)
	{
		if(oligo_check(s, i, mmis))
		{
			string::iterator nsi = si;

			if(si != mask_seq.begin())
				nsi--;

			for(int j = -1; j <= (int)s->size(); j++)
			{
				if(i+j < 0 || i+j >= (int)SEQ.size())
					continue;

				*nsi = OR_SEQ[i+j];
				nsi++;
			}
		}

		si++;
	}

	return;
}	

bool sequence::oligo_check(string *s, int off, int mmis)
{
	int count = 0, pos = 0;
	bool flag = true;

	for(int i = off; i < off + s->size(); i++)
	{
		if(s->at(pos) != OR_SEQ[i])
			count++;
	
		if(count > mmis)
		{
			flag = false;
			break;
		}

		pos++;
	}

	if(flag || !DOUBLE_STRAND)
		return flag;

	else
	{
		string rc = string_rev_comp(s);
		count = 0;
		pos = 0;

		for(int i = off; i < off + rc.size(); i++)
		{
			if(rc.at(pos) != OR_SEQ[i])
				count++;	

			if(count > mmis)
				return false;

			pos++;
		}

	}

	return true;
}

void sequence::mask_sequence(vector<string> *oligo,int mmis, map<string,double> *SCORES)
{
	mmis = OLIGO_MAX_MIS[mmis];

	for(int i = 0; i < oligo->size() && i < TOPN; i++)
	{
		if(SCORES->find(oligo->at(i))->second > 0 || !i)
			masker(&oligo->at(i), mmis);	
	}

	SEQ = mask_seq;
	
	mask_seq.clear();

	mask_seq.insert(mask_seq.begin(), SEQ.size(), 'N');

	return;
}

bool sequence::seq_cleaner()
{
        if(!SEQ.size())
                return false;

	GC_COUNT = 0;

        for(int i = 0; i < SEQ.size(); i++)
        {
                bool flag = false;

                SEQ[i] = toupper(SEQ[i]);

		mask_seq.push_back('N'); //masked!!

                for(int a = 0; a < ALPHABET_SIZE; a++)
                {
                        if(SEQ[i] == ALPHABET[a])
			{
                                flag = true;

		//**********NUCLEOTIDE GC PERCENT******************************
                              if(a == 1 || a == 2)   //GC PERCENT
                                      GC_COUNT++;
                //***************************************************************

		//**********DINUCLEOTIDE GC PERCENT******************************
                //                if(i < SEQ.size() - 1)
                //                {
                //                        if((a == 2 && (toupper(SEQ[i+1]) == ALPHABET[1])) || (a == 1 && (toupper(SEQ[i+1]) == ALPHABET[2])))
                //                                GC_COUNT++;
                //                }
                //***************************************************************

				break;
			}
                }

                if(!flag && SEQ[i] != 'N')
                        return false;
        }

	used_reset();

	OR_SEQ = SEQ;

        return true;
}

void sequence::used_reset()
{

	used.clear();

	for(int i = 0; i < SEQ.size(); i++)
		used.push_back(false);

	return;
}

string sequence::name()
{
        return NAME;
}

string sequence::seq()
{
	return SEQ;
}

string sequence::or_seq()
{
	return OR_SEQ;
}

void sequence::set_used_map(int l)
{
	for(int i = 0; i < OLIGO_TYPES; i++)
	{
		if(l == OLIGO_LENGTHS[i])
		{
			used_map = i;
			used_length = l;
			return;
		}
	}

	error_handler(BAD_SCANNING_OLIGO,"");
}

void sequence::oligo_scan(string oligo)
{
	string rc_oligo = string_rev_comp(&oligo);
	int size = (int)oligo.size(), oligo_max_mis = OLIGO_MAX_MIS[used_map];
	vector<bool> uo(SEQ.size(), false);

	unsigned int o = 0;


	for(int i = 0; i < SEQ.size() - size; i++)
	{
/*		unsigned int dimiss = 0;

		for(int j = 0; j < size; j += 2)
		{
			if(((int)oligo[j] - (int)oligo[j+1]) == ((int)SEQ[i+j] - (int)SEQ[i+j+1]))
				dimiss++;
			
		}

		if(dimiss < 2)
			continue;	
*/ //BAH!
		int count = 0;

		for(int j = 0; j < size && count <= oligo_max_mis; j++)
		{
			if(oligo[j] != SEQ[i+j])
				count++;

			if(count > oligo_max_mis)
				break;
		}

		if(count > oligo_max_mis)
			continue;
		else if(count == oligo_max_mis)
		{
			o++;
			uo[i] = true;
		}
		else
		{
			o = 1;
			oligo_max_mis = count;
			uo.resize(SEQ.size(), false);
			uo[i] = true;
		}

		if(!count)
			used[i] = true;
	}

	if(rc_oligo != oligo && DOUBLE_STRAND)
	{
		for(int i = 0; i < SEQ.size() - size; i++)
		{
			if(uo[i])
				continue;

			int count = 0;

			for(int j = 0; j < size && count <= oligo_max_mis; j++)
			{
				if(rc_oligo[j] != SEQ[i+j])
					count++;

				if(count > oligo_max_mis)
					break;
			}

			if(count > oligo_max_mis)
                        	continue;

			else if(count == oligo_max_mis) 
				o++;
			else
                	{
                        	o = 1;
                        	oligo_max_mis = count;
                	}

			if(!count)
				used[i] = true;
		}
	}

	if(MAIN_OLIGO_MAP.find(oligo) == MAIN_OLIGO_MAP.end())
        	MAIN_OLIGO_MAP[oligo] = 0;

	map<string,oli>::iterator mmi = OLI_MAP->find(oligo);

	if(mmi == OLI_MAP->end())
        	mmi = OLI_MAP->insert(make_pair(oligo,oli(oligo,oligo_max_mis))).first;

	if(o)
		MAIN_OLIGO_MAP[oligo] += log(o / (get_score(&mmi->second, oligo_max_mis, used_map) * SEQ.size()));

	return;
}		

void sequence::mat_oligo_scan(string oligo, vector<string> *occ)
{
	string rc_oligo;
	int size = (int)oligo.size(), oligo_max_mis = OLIGO_MAX_MIS[get_type_from_l(size)];
	vector<bool> uo(SEQ.size(), false);

	if(DOUBLE_STRAND)
		rc_oligo = string_rev_comp(&oligo);

	for(int i = 0; i < SEQ.size() - size; i++)
	{
		int count = 0;

		for(int j = 0; j < size && count <= oligo_max_mis; j++)
		{
			if(oligo[j] != SEQ[i+j])
				count++;
   //G                     if(count > oligo_max_mis)
     //G                    break;
		}

		if(count <= oligo_max_mis)
		{
			string o = SEQ.substr(i,size);			

			if(o.find("N") == string::npos)
			{
				occ->push_back(o);
				uo[i] = true;
			}
		}
	}

	if(rc_oligo != oligo && DOUBLE_STRAND)
	{
		for(int i = 0; i < SEQ.size() - size; i++)
		{
			if(uo[i])
                                continue;

			int count = 0;

			for(int j = 0; j < size && count <= oligo_max_mis; j++)
			{
				if(rc_oligo[j] != SEQ[i+j])
					count++;
       //G                         if(count > oligo_max_mis)
         //G                       break;
			}

			if(count <= oligo_max_mis)
			{
				string o = SEQ.substr(i,size);

				if(o.find("N") == string::npos)
				{
					occ->push_back(string_rev_comp(&o));
					uo[i] = true;
				}
			}

		}
	}

	return;
}		
		
void sequence::scan(vector<sequence> *VS, int length, int spos)
{
	cerr << (char)13 << NAME.substr(0,30) << '\t' << spos+1 << "\\" << VS->size(); 
	
	if(length != used_length)
	{
		for(int i = 0; i < VS->size(); i++)
			VS->at(i).set_used_map(length);
	}

	for(int i = 0; i < SEQ.size() - length; i++)
	{
		string oligo = SEQ.substr(i,length);

		if(!used[i]) 
			if(oligo.find("N") == string::npos)
			{
				for(int j = 0; j < VS->size(); j++) 
					VS->at(j).oligo_scan(oligo);
			}
	}

	return;
}

/*void sequence::o_Clean(vector<sequence> *VS)
{
	for(int i = 0; i < oligo_occ_maps[used_map].size(); i++)
	{
		for(map<string,unsigned short int>::iterator mi = oligo_occ_maps[used_map][i].begin(); mi != oligo_occ_maps[used_map][i].end(); mi++)
		{
			unsigned short int count = 0;
			string mis = mi->first;

			for(int j = 0; j < VS->size(); j++)
				for(int k = 0; k < VS->at(j).oligo_occ_maps[used_map].size(); k++)
				{
					map<string,unsigned short int>::iterator mmi = VS->at(j).oligo_occ_maps[used_map][k].find(mis);

					if(mmi->second)
					{
						count++;
						break;
					}
				}
	

			if(count < ceil((double)VS->size()/2.0))	
			{
				for(int g = 0; g < VS->size(); g++) 
					for(int q = 0; q < VS->at(g).oligo_occ_maps[used_map].size(); q++)
	//					VS->at(g).oligo_occ_maps[used_map][q].erase(mis); //DO NOT WORK ON MACOSX
					       	if(VS->at(g).oligo_occ_maps[used_map][q].find(mis) != VS->at(g).oligo_occ_maps[used_map][q].end())
							VS->at(g).oligo_occ_maps[used_map][q].find(mis)->second = 0; 
			}
		}
	}

	return;
}	
*/					

double sequence::get_score(oli *oligo, int lev, int type)
{
	return oligo->ofreq->at(lev + ((OLIGO_MAX_MIS[type] + 1) * GC));
}	

/*void sequence::occ_clear()
{
	oligo_occ_maps.clear();
	return;
}*/


#ifdef DEBUG5 
string sequence::see_scores(int type)
{
	ostringstream out;
	for(map<string,double>::iterator mi = scores[type].begin(); mi != scores[type].end(); mi++)
		out << mi->first << '\t' << mi->second << endl;

	return out.str();
}
#endif

#ifdef DEBUG4
string oli::see_level(int lev)
{
	ostringstream str;

	for(int i = 0; i < o_nhr[lev].size(); i++)
	{
		str << o_nhr[lev][i] << '\t';
		for (int j = 0; j < u_pos[lev][i].size(); j++)
			str << u_pos[lev][i][j] << '\t';
		str << endl;
	}

	return str.str();
}
#endif

/* void oli::give_level(int lev)
{
	if(lev <= mis_level)
		return; 

	else
		while(mis_level < lev)
			raise_one_level();
	return; 
}*/

/*void oli::raise_one_level()
{
	vector<map<string,bool>::iterator > N_OLIGO;
	vector<vector<char> > POS;

	u_pos.push_back(POS);

	o_nhr.push_back(N_OLIGO);

	int o_size = (int)o_nhr[0][0]->first.size();
	int oos = (int)u_pos[mis_level].size();

	for(int l = 0; l < oos; l++)
	{
		char pos = 0;

		int uml = (int)u_pos[mis_level][l].size();

		while(pos < o_size)
		{
			int os = (int)u_pos[mis_level + 1].size();

			bool flag = false;
			
			for(int t = 0; t < uml; t++)
			{
				if(pos <= u_pos[mis_level][l][t])
					flag = true;
			}

			if(flag)
			{
				pos++;
				continue;
			}


			for(int b = 0; b < ALPHABET_SIZE - 1; b++)
				u_pos[mis_level + 1].push_back(u_pos[mis_level][l]);

			vector<string> n_oligo(ALPHABET_SIZE - 1);

			for(int q = 0; q < ALPHABET_SIZE - 1; q++)
				n_oligo[q].assign(o_size,'N');

			for(int j = 0; j < o_size; j++)
			{
				if(j == pos)
				{
					int count = 0;

					for(int k = 0; k < ALPHABET_SIZE; k++)
					{
						if(ALPHABET[k] != o_nhr[mis_level][l]->first.at(j))
						{
							n_oligo[count][j] = (ALPHABET[k]);
							count++;
						}	
					}
				}

				else
					for(int h = 0; h < n_oligo.size(); h++)
						n_oligo[h][j] = o_nhr[mis_level][l]->first.at(j);
			}


			for(int h = 0; h < n_oligo.size(); h++)
			{
				o_nhr[mis_level + 1].push_back(mis_map->insert(mis_map->begin(),make_pair(n_oligo[h],true)));
				u_pos[mis_level + 1][h+os].push_back(pos);
			}

			pos++;
		}
	}

	mis_level++;

	if(!no_freq)
		ofreq.push_back(get_freq());
	else
		ofreq.push_back(get_freq_without_freqs());
	
	if(mis_level == OLIGO_MAX_MIS[Type])
	{
		for(int i = 1; i <= OLIGO_MAX_MIS[Type]; i++)
			o_nhr[i].clear();			

		u_pos.clear();
		max_lev = true;
	}

	return;
}*/

unsigned short int sequence::gc()
{
        return GC;
}
					
#ifdef DEBUG1
string sequence::check_1()
{
	ostringstream str;

	str << NAME << endl;

	for(int i = 0; i < oligo_occ_maps.size(); i++)
	{
		for(int j = 0; j < oligo_occ_maps[i].size(); j++)
		{
			str << "\nmismatch = " << j << endl;
			for(map<string,int>::iterator mi = oligo_occ_maps[i][j].begin(); mi != oligo_occ_maps[i][j].end(); mi++)
				str << mi->first << '\t' << mi->second << endl;

		}
	}

	return str.str();
}
#endif

oli::oli(string s, int level)
{
//	vector<map<string,bool>::iterator > vs;
  //      vector<vector<char> > pv(1);
	Type = get_type_from_l((short int)s.size());

//	max_lev = false;

  //      pv[0].push_back(-1);

//	vs.push_back(mis_map->insert(mis_map->begin(),make_pair(s,true)));

//	o_nhr.push_back(vs);
//        u_pos.push_back(pv);

       mis_level = 0;
	
//	if(OF->at(Type).find(s) != OF->at(Type).end())
	
		ofreq = &OF->at(Type).find(s)->second;			
//	else
//		ofreq.push_back(infer_freq(&s, &OF->at(Type-1)));
/*
	if(DOUBLE_STRAND)
	{
		string rc = string_rev_comp(&s);

		if(rc != s)	
		{
	//		if(OF->at(Type).find(rc) != OF->at(Type).end())
	//		{
			for(short int i = 0; i < ofreq.size(); i++)
				ofreq[i] += OF->at(Type).find(rc)->second.at(i);
	//		}
//			else
//				ofreq[0] += infer_freq(&rc, &OF->at(Type-1));
		}
	}*/
/*
	while(mis_level < level && !max_lev)
		raise_one_level(); */

	return;
}

/*oli::~oli()
{
	o_nhr.clear();
	u_pos.clear();
	ofreq.clear();

	return;
}*/
		

/*double oli::get_freq()
{
        bool flag = false;
	vector<map<string,bool>::iterator> *olist = &o_nhr[mis_level];


        double sum = 0;

        for(int i = 0; i < olist->size(); i++)
		  sum += (double)OF->at(Type).find(olist->at(i)->first)->second;

	if(DOUBLE_STRAND)
	{
		for(int i = 0; i < olist->size(); i++)
		{
			string rc = string_rev_comp(&olist->at(i)->first);

			if(pal_check(&rc))
				sum += (double)OF->at(Type).find(rc)->second;
		}
	}

        return sum;
}*/

/*double oli::get_freq_without_freqs()
{
        bool flag = false;
	vector<map<string,bool>::iterator> *olist = &o_nhr[mis_level];

        double sum = 0;

        int nf = (int)OF->size() - 1;
        int of = (int)OF->size() - 2;

        for(int i = 0; i < olist->size(); i++)
        {
                if(OF->at(nf).find(olist->at(i)->first) != OF->at(nf).end())
                        sum += (double)OF->at(nf).find(olist->at(i)->first)->second;
                else
                {
                        double newfreq = infer_freq(&olist->at(i)->first, &OF->at(of));

                        sum += newfreq;

                        OF->at(nf).insert(make_pair(olist->at(i)->first, newfreq));
                }
        }

	if(DOUBLE_STRAND)
	{
		for(int i = 0; i < olist->size(); i++)
		{
			string rc = string_rev_comp(&olist->at(i)->first);

			if(pal_check(&rc))
			{
				if(OF->at(nf).find(rc) != OF->at(nf).end())
					sum += (double)OF->at(nf).find(rc)->second;
				else
				{
					double newfreq = infer_freq(&rc, &OF->at(of));
	
					sum += newfreq;
			
					O->at(nf).insert(make_pair(rc, newfreq));
				}	
			}
		}
	}

        return sum;
}*/

/*bool oli::pal_check(string *rc)
{
	int count = 0;

	for(int i = 0; i < o_nhr[0][0]->first.size(); i++)
	{
		if(rc->at(i) != o_nhr[0][0]->first.at(i))
			count++;		
		
		if(count > mis_level)
			return true;
	}

	return false;
}*/	

matrix::matrix(string oligo, vector<sequence> *SEQUENCE, vector<matrix> *MATRIX, unsigned int matpos)
{
	sstart = oligo;
	mpos = matpos;
	double a = 0;

	s_ptr = SEQUENCE;

	NSITES = ceil(double((double)s_ptr->size() / 2.0));

	y_s = ALPHABET_SIZE;

	x_s = (int)sstart.size();

	Type = get_type_from_l(x_s);

	sscore = SC_PTR->find(sstart)->second;

	vector<string> occur;

	best_occ = new int*[(int)SEQUENCE->size()];
	best_occ_score = new double[(int)SEQUENCE->size()];

	for(int i = 0; i < (int)SEQUENCE->size(); i++)
	{
		SEQUENCE->at(i).mat_oligo_scan(sstart, &occur);

		vector<bool> vb((int)SEQUENCE->at(i).seq().size(),false);
		vector<double> vd((int)SEQUENCE->at(i).seq().size() - x_s, 0);

		zrc.push_back(vb);
		Z.push_back(vd);
		best_occ[i] = new int[2];
	}

	s_matrix = new double*[x_s];	
	n_matrix = new double*[x_s];

	for(int i = 0; i < x_s; i++)
	{
		s_matrix[i] = new double[y_s];
		n_matrix[i] = new double[y_s];
	}


	for(int i = 0; i < x_s; i++)
		for(int j = 0; j < y_s; j++)
			s_matrix[i][j] = 0;

	for(int i = 0; i < occur.size(); i++)
		for(int j = 0; j < x_s; j++)
		{
//			cerr << occur[i] << endl;
//			cerr.flush();
			s_matrix[j][ati(occur[i][j],&occur[i])]++;
		}


	for(int j = 0; j < y_s; j++)
			a += s_matrix[0][j];

	if(a < ceil(double((double)SEQUENCE->size() / 2.0)))
		good_matrix = false;
	else
		good_matrix = true;

	if(good_matrix)
	{
		normalize();

		not_too_similar = similarity_check(MATRIX);	

		for(int i = 0; i < EM_CYCLES && not_too_similar; i++)
		{
			get_z();
			squash();
			get_p();

			if(i < EM_CYCLES - 1)
			{
				for(int j = 0; j < zrc.size(); j++)
					for(int k = 0; k < zrc[j].size(); k++)
						zrc[j][k] = false;
			}
		}

		if(not_too_similar)
		{
	//		generate_pfm(n_matrix);
			set_min_max();
			get_z();  //Motif selection//
			scan_Z(); //////////////////
		}
	}


	return;
}

bool matrix::similarity_check(vector<matrix> *MATRIX)
{
//	cerr << "\nMAT: " << mpos  << " simcheck. Start oligo = " << sstart << endl << endl;

	double max_dist = 100;
	for(unsigned int i = 0; i < MATRIX->size(); i++)
	{
		//	cerr << "MAT: " << i << "(" << MATRIX->at(i).get_sstart() << ") = ";

			if(x_s != MATRIX->at(i).size())
			{
		//		cerr << -1 << endl;
				continue;
			}


			for(unsigned int z = 0; z <= OLIGO_MAX_MIS[get_type_from_l(x_s)]; z++) 
			{
				for(unsigned int v = 0; v <= z; v++)
				{
					double dist = 0;

					for(unsigned int a = 0 + v; a < x_s - (z+v); a++)
					{
						double buf = 0;

						for(unsigned int b = 0; b < y_s; b++)
							buf += pow(n_matrix[a][b] - MATRIX->at(i).get_n_element(a-v,b),2);

						buf *= (1.0/sqrt(2.0));

						dist += buf;
					}

					dist *= 1.0/(double)(x_s - z);

					if(dist < max_dist)
						max_dist = dist;

	//				cerr << endl << "z = " << z << '\t' << "v = " << v << '\t' << "maxdist = " << 1.0 - max_dist << '\t' << "dist = " << 1.0 - dist;

					if(1.0 - max_dist >= SIM_THRESHOLD)
						break;
					//	return false;
				}

				if(1.0 - max_dist >= SIM_THRESHOLD)
					break;
                                              //  return false;

				for(unsigned int v = 1; v <= z; v++)
                               	{
                                       	double dist = 0;

                                       	for(unsigned int a = 0 + v; a < x_s - (z+v); a++)
                                       	{
                                               	double buf = 0;

                                               	for(unsigned int b = 0; b < y_s; b++)
                                                       	buf += pow(n_matrix[a - v][b] - MATRIX->at(i).get_n_element(a,b),2);

                                               	buf *= (1.0/sqrt(2.0));

                                               	dist += buf;
                                       	}

                                       	dist *= 1.0/(double)(x_s - z);

                                       	if(dist < max_dist)
                                               	max_dist = dist;

                                  //     	cerr << endl << "z1 = " << z << '\t' << "v1 = " << v << '\t' << "maxdist = " << 1.0 - max_dist << '\t' << "dist = " << 1.0 - dist;

                                       	if(1.0 - max_dist >= SIM_THRESHOLD)
						break;
                                               //	return false;
                               	}


				if(1.0 - max_dist >= SIM_THRESHOLD)
					break;
                                             //   return false;

				if(DOUBLE_STRAND)
				{
					for(unsigned int v = 0; v <= z; v++)
					{	
						double dist_r = 0;

						for(unsigned int a = 0 + v; a < x_s - (z+v); a++)
                        			{
                                			double buf = 0;

                                			for(unsigned int b = 0; b < y_s; b++)
                                        			buf += pow(n_matrix[x_s - a - 1][rc_index(b)] - MATRIX->at(i).get_n_element(a-v,b),2);

                                			buf *= (1.0/sqrt(2.0));

                                			dist_r += buf;
                        			}

						dist_r *= 1.0/(double)(x_s - z);

						if(dist_r < max_dist)
							max_dist = dist_r;

		//				cerr << endl << "-z = " << z << '\t' << "v = " << v << '\t' << "maxdist = " << 1.0 - max_dist << '\t' << "dist = " << 1.0 - dist_r;

						if(1.0 - max_dist >= SIM_THRESHOLD)
							break;
                                       // 		return false;
					}

					if(1.0 - max_dist >= SIM_THRESHOLD)
						break;
                                                     //   return false;

					for(unsigned int v = 1; v <= z; v++)
                                        {
                                                double dist_r = 0;

                                                for(unsigned int a = 0 + v; a < x_s - (z+v); a++)
                                                {
                                                        double buf = 0;

                                                        for(unsigned int b = 0; b < y_s; b++)
                                                                buf += pow(n_matrix[x_s - (a - v) - 1][rc_index(b)] - MATRIX->at(i).get_n_element(a,b),2);

                                                        buf *= (1.0/sqrt(2.0));

                                                        dist_r += buf;
                                                }

                                                dist_r *= 1.0/(double)(x_s - z);

                                                if(dist_r < max_dist)
                                                        max_dist = dist_r;

                                            //    cerr << endl << "-z1 = " << z << '\t' << "v1 = " << v << '\t' << "maxdist = " << 1.0 - max_dist << '\t' << "dist = " << 1.0 - dist_r;

                                                if(1.0 - max_dist >= SIM_THRESHOLD)
							break;
                                                       // return false;
                                        }

					if(1.0 - max_dist >= SIM_THRESHOLD)
						break;
                                                    //    return false;
				}
			}


//			cerr << endl << 1.0 - max_dist << endl;

			if(1.0 - max_dist >= SIM_THRESHOLD)
                		return false;
	}	

//	if(1.0 - max_dist >= SIM_THRESHOLD)
//		return false;

	return true;
}

void matrix::scan_Z()
{
 	int count = 0;


	while(count < s_ptr->size())
	{
		double min = 0;

		for(int i = 0; i < (int)Z.size(); i++)
			for(int j = 0; j < (int)Z[i].size(); j++)
				if(Z[i][j] > min)
				{
					min = Z[i][j];
					best_occ[count][Seq] = i;
					best_occ[count][Pos] = j;
				}

//		cerr << "best_occ[count[Seq]] = " << best_occ[count][Seq] << '\t' << "best_occ[count[Pos]] = " << best_occ[count][Pos] << endl;
		best_occ_score[count] = Z[best_occ[count][Seq]][best_occ[count][Pos]] / NSITES;
		Z[best_occ[count][Seq]][best_occ[count][Pos]] = 0;
		count++;
	}

//	cerr << endl << endl;

	return;
}

/*void matrix::see_matrix(double **matrix)
{
	cout << sstart << endl;

	if(good_matrix)
		cerr << "GOOD" << endl;
	else
		cerr << "BAD" << endl;

	for(int i = 0; i < y_s; i++)
	{
		cout << ALPHABET[i] << '\t';

		for(int j = 0; j < x_s; j++)
			cout << matrix[j][i] << '\t';

		cout << endl;
	}

	return;
}*/

/*void matrix::generate_pfm(double **matrix)
{
	ofstream mout(mlist.c_str(), ios::app);

	ostringstream mfile;

	mfile << "MAT" << M_NUM << "_" << sstart << ".pfm";

	mout << mfile.str() << endl;

	ofstream out(mfile.str().c_str());

        for(int i = 0; i < y_s; i++)
        {
                for(int j = 0; j < x_s; j++)
                        out << matrix[j][i] << ' ';

                out << endl;
        }

	M_NUM++;

        return;
}*/

void matrix::normalize()
{
	double sum = 0, corr; 

	for(int j = 0; j < y_s; j++)
                       sum += s_matrix[0][j];		

	corr = sqrt(sum);

	for(int x = 0; x < x_s; x++)
                for(int y = 0; y < y_s; y++)
			n_matrix[x][y] = s_matrix[x][y] + corr;

	sum = 0;

	for(int j = 0; j < y_s; j++)
                       sum += n_matrix[0][j];

	for(int x = 0; x < x_s; x++)
		for(int y = 0; y < y_s; y++)
			n_matrix[x][y] = (n_matrix[x][y]/sum); 	

//	set_min_max();

	return;
}

void matrix::set_min_max()
{
	min = 0;
	max = 0;

	for(int x = 0; x < x_s; x++)
	{
		double row_min = 1, row_max = 0;
			
                for(int y = 0; y < y_s; y++)
		{
			if(n_matrix[x][y] < row_min)
				row_min = n_matrix[x][y];
			
			if(n_matrix[x][y] > row_max)
				row_max = n_matrix[x][y];
		}

		min += row_min;
		max += row_max;
	}

	return;
}

double matrix::relative_score(string *oligo)
{
	double rscore = 0;	

	for(unsigned int i = 0; i < oligo->size(); i++)
		rscore += n_matrix[i][ati(oligo->at(i), oligo)];		

	rscore =  1 - ((max - rscore) / (max - min));

	return rscore;
}

void matrix::get_z()
{
//	double NSITES = ceil(double((double)s_ptr->size() / 2.0)), SUM = 0;	
	double SUM = 0;

	for(int i = 0; i < (int)s_ptr->size(); i++)
	{
		string seq = s_ptr->at(i).seq();

		for(int j = 0; j < (int)seq.size() - x_s; j++)
		{
			double z = 1, rz = 1;
			string ss = seq.substr(j,x_s);
			string rs = string_rev_comp(&ss);

			if(ss.find('N') != string::npos)
			{
		//		Z[i].push_back(0);
				Z[i][j] = 0;
				continue;
			}


			for(int k = 0; k < x_s; k++)
				z *= n_matrix[k][ati(ss[k],&ss)];

			if(DOUBLE_STRAND)
			{
				for(int k = 0; k < x_s; k++)
					rz *= n_matrix[k][ati(rs[k],&rs)];

				if(rz > z)
				{
					z = rz;
					ss = rs;
					zrc[i][j] = true;
				}
			}

//			if(OF->at(Type).find(ss) != OF->at(Type).end())
                		z /= (OF->at(Type).find(ss)->second[0]);
//			else
/*			{
				double f;
				f = infer_freq(&ss, &OF->at(Type - 1));
				OF->at(Type).insert(make_pair(ss, f));
				z /= f;
			}*/

			Z[i][j] = z;
		}
	}

	for(int i = 0; i < Z.size(); i++)
		for(int j = 0; j < Z[i].size(); j++)
			SUM += Z[i][j];	

	for(int i = 0; i < Z.size(); i++)
                for(int j = 0; j < Z[i].size(); j++)
		{
			Z[i][j] /= SUM;
			Z[i][j] *= NSITES;
		}

	return;
}

void matrix::get_p()
{
	for(int x = 0; x < x_s; x++)
		for(int y = 0; y < y_s; y++)
			n_matrix[x][y] = 0;

	for(int i = 0; i < Z.size(); i++)
	{
		for(int j = 0; j < Z[i].size(); j++)
		{
			if(Z[i][j] != 0)
			{
				string s;

				s = s_ptr->at(i).seq().substr(j,x_s);

				if(zrc[i][j])
					s = string_rev_comp(&s);

			//	cerr << s << endl;

				for(int k = 0; k < x_s; k++)
					n_matrix[k][ati(s[k],&s)] += Z[i][j];	
			}
		}
	}

	for(int x = 0; x < x_s; x++)
	{
		double sum = 0;

                for(int y = 0; y < y_s; y++)	
			sum += n_matrix[x][y];

		for(int x = 0; x < x_s; x++)
                	for(int y = 0; y < y_s; y++)
				n_matrix[x][y] /= sum;
	}


	return;
}

void matrix::squash()
{
        double total = 0;
	double Nsites = ceil(double((double)s_ptr->size() / 2.0));

        bool renorm = true;

	for(int i = 0; i < Z.size(); i++)
		for(int j = 0; j < Z[i].size(); j++)
			total += Z[i][j];

        while(renorm)
        {
                renorm = false;

                double norm = total/Nsites;
                total = 0;

                for(int i1 = 0; i1 < Z.size(); i1++)
                        for(int i2 = 0; i2 < Z[i1].size(); i2++)
                        {
                                double p = Z[i1][i2];

                                if(p < 1)
                                        p /= norm; 

                                if(p > 1)
                                {
                                        p = 1;
                                        Nsites--;
                                        renorm = true;
                                }

                                Z[i1][i2] = p;

                                if(p <= 1)
                                        total += p;
                        }
        }

        return;
}

string matrix::get_sstart()
{
	return sstart;
}

/*void matrix::write_tompa_format()
{
	cout << ">dataset\n";

	cout << inputfile.substr(0,inputfile.find(".")) << endl;

	cout << ">instances\n";

	for(int j = 0; j < s_ptr->size(); j++)	
	{
		string s = s_ptr->at(best_occ[j][Seq]).seq().substr(best_occ[j][Pos],OLIGO_LENGTHS[Type]);	

		cout << best_occ[j][Seq] << ','
		    << -(int)s_ptr->at(best_occ[j][Seq]).seq().size() + best_occ[j][Pos] << ','
		    << s << endl;
	}

	return;
}
*/
bool matrix::sim_check()
{
	return not_too_similar;
}

unsigned int matrix::motif_pos()
{
	return mpos;
}

double matrix::get_n_element(unsigned int x, unsigned int y)
{
	return n_matrix[x][y];
}

double matrix::get_s_element(unsigned int x, unsigned int y)
{
        return s_matrix[x][y];
}

unsigned int matrix::size()
{
	return x_s;
}

double matrix::nsites()
{
	return NSITES;
}

double matrix::get_min()
{
	return min;
}

double matrix::get_max()
{
	return max;
}
