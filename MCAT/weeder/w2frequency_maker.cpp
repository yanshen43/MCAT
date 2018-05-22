#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <math.h>
#include <cstdlib>

using namespace std;

//typedef unsigned long int my_int;
const char ALPHABET[4] = {'A','C','G','T'};
const char R_ALPHABET_TABLE[20] = {'T','N','G','N','N','N','C','N','N','N','N','N','N','N','N','N','N','N','N','A'};
const short int ALPHABET_SIZE = 4, OLIGO_TYPES = 3;
const short int OLIGO_LENGTHS[OLIGO_TYPES] = {6,8,10}, OLIGO_MAX_MIS[OLIGO_TYPES] = {1,2,3};

enum {SYNTAX_ERROR, MISSING_FASTA_FILE, BAD_FASTA_SEQUENCE};

class sequence
{
    private:
        bool GOOD;
        string NAME;
        string SEQ;
        bool seq_cleaner();
    public:
        sequence(string);
        bool good();
        string name();
        string seq();
};

string inputfile, scode;
bool DOUBLE_STRAND = false;

void fasta_reader(const char*, vector<sequence> *);
void error_handler(int, string);
void count_ss_oligo(string*,vector<sequence> *, short int);
void count_ds_oligo(string*,vector<sequence> *, short int);
void oligo_gen(short int, short int, map<string, vector<unsigned int> >*, string); 
void near_oligo_gen(short int, short int, string *, map<string, vector<unsigned int> > *);
string string_rev_comp(const string*);

int main(int argc, char **argv)
{
    if(argc != 4)
        error_handler(SYNTAX_ERROR,"");

    vector<sequence> SEQUENCE;
    inputfile = argv[1];
    scode = argv[2];
    string mode = argv[3];

    if(mode != "ss" && mode != "ds")
        error_handler(SYNTAX_ERROR,"");

    cerr << "\nReading sequences file (avoiding duplicates)...\t";

    fasta_reader(inputfile.c_str(), &SEQUENCE);

    if(mode == "ss")
        for(short int i = 0; i < OLIGO_TYPES; i++)
            count_ss_oligo(&inputfile,&SEQUENCE,i);
    else
    {
        DOUBLE_STRAND = true;

        for(short int i = 0; i < OLIGO_TYPES; i++)
            count_ds_oligo(&inputfile,&SEQUENCE,i);
    }

    return EXIT_SUCCESS;
}


void fasta_reader(const char *f_name, vector<sequence> *SEQUENCE)
{
    ifstream in(f_name);
    string line,buf;
    vector<string> seqs;
    bool flag = false;
    unsigned int dupcount = 0, totcount = 0;

    if(!in)
    {
        string tmp = f_name;
        error_handler(MISSING_FASTA_FILE, tmp);
    }

    while(getline(in,line))
    {
        if(line[0] == '>' && flag)
        {
            sequence new_seq(buf);
            bool nflag = false;

            if(new_seq.good())
                for(int i = 0; i < SEQUENCE->size(); i++)
                    if(new_seq.seq() == SEQUENCE->at(i).seq())	
                    {	
                        nflag = true;
                        dupcount++;
                        break;
                    }

            if(!nflag && new_seq.good())
                SEQUENCE->push_back(new_seq);	

            buf.clear();
            buf += line;
            buf += '\n';
            totcount++;
        }

        else if(line[0] == '>' && !flag)
        {
            buf += line;
            buf += '\n';
            flag = true;
            totcount++;
        }

        else
        {
            buf += line;
            buf += '\n';
        }
    }

    sequence fnew_seq(buf);
    bool fflag = false;

    for(int i = 0; i < SEQUENCE->size(); i++)
        if(fnew_seq.seq() == SEQUENCE->at(i).seq())
        {
            fflag = true;
            dupcount++;
            break;
        }

    if(!fflag && fnew_seq.good())
        SEQUENCE->push_back(fnew_seq);

    buf.clear();

    in.close();

    cerr << "done - Read sequence(s): " << totcount << " Good Sequences: " << SEQUENCE->size() << " Duplicated sequence(s): " << dupcount << endl;

    return;
}

void error_handler(int error, string s_err)
{
    cerr << endl;

    switch(error)
    {
        case SYNTAX_ERROR:
            {
                cerr << "\nSyntax: ofreq filename species_code ss|ds\n" << endl;
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
        default:
            {
                cerr << "\nSome weird error occurred...\n";
                break;
            }
    }

    cerr << endl;

    exit(error);
}

void count_ss_oligo(string *inputfile, vector<sequence> *SEQUENCE, short int a)
{
    ostringstream outfile;
    map<string, vector<unsigned int> >  m_freq; // (1 + OLIGO_MAX_MIS[a], map<string, unsigned int> ());
    string oligo(OLIGO_LENGTHS[a], 'N');

    //	outfile << *inputfile << "_" << OLIGO_LENGTHS[a] << "_ss.freq";
    outfile << scode << "_ss." << OLIGO_LENGTHS[a] << ".freq";

    ofstream out(outfile.str().c_str());

    cerr << "Building " << OLIGO_LENGTHS[a] << "mers list...";
    oligo_gen(a, 0, &m_freq, oligo);
    cerr << "\t\tdone" << endl;

    map<string, unsigned int> omap;

    cerr << "Finding exact occurrences...";

    for(unsigned int i = 0; i < SEQUENCE->size(); i++)
    {
        for(unsigned int j = 0; j < SEQUENCE->at(i).seq().size() - OLIGO_LENGTHS[a]; j++)
        {
            for(unsigned int k = 0; k < OLIGO_LENGTHS[a]; k++)
                oligo[k] = SEQUENCE->at(i).seq().at(j+k);	

            if(oligo.find("N") != string::npos)
                continue;

            if(omap.find(oligo) == omap.end())	
                omap[oligo] = 1;
            else
                omap[oligo]++;
        }
    }

    cerr << "\tdone" << endl;


    map<string, unsigned int>::iterator omi = omap.begin();
    unsigned int counter = 0;

    while(omi != omap.end())
        //	for(unsigned int i = 0; i < SEQUENCE->size(); i++)
    {
        cerr << (char)13 << "Counting " << OLIGO_LENGTHS[a] << "mers occurrences..." << counter+1 << '/' << omap.size();

        //		for(unsigned int j = 0; j < SEQUENCE->at(i).seq().size() - OLIGO_LENGTHS[a]; j++)
        //	{
        //			for(unsigned int k = 0; k < OLIGO_LENGTHS[a]; k++)
        //				oligo[k] = SEQUENCE->at(i).seq().at(j+k);

        //			if(oligo.find("N") != string::npos)
        //				continue;

        oligo = omi->first;

        map <string, vector<unsigned int> > oligos;

        oligos[oligo] = vector<unsigned int> ();

        //			cout << "OLIGO: " << oligo << endl << "**********" << endl;

        near_oligo_gen(a,1,&oligo,&oligos);

        map<string, vector<unsigned int> >::iterator mi = oligos.begin();

        /*			while(mi != oligos.end())
                    {
                    cout << mi->first << '\t' << mi->second.size() << endl;
                    mi++;
                    }*/

        while(mi != oligos.end())
        {
            //				if(omi->first == "AAATTTTT")
            //                              	cout << mi->first << '\t' << mi->second.size() << endl;
            //			if(mi->first == "TAATTT" || mi->first == "AAATTA")
            //                      	cout << endl << omi->first << '\t' << omi->second  << '\t' << mi->second.size() << '\t' << mi->first << endl;

            m_freq[mi->first][mi->second.size()] += omi->second;	
            mi++;
        }

        /*			map<string, vector<unsigned int> >::iterator mi = m_freq.begin(); //SLOOOW

                    while(mi != m_freq.end())
                    {
                    unsigned int mcount = 0;

                    for(unsigned int l = 0; l < OLIGO_LENGTHS[a]; l++)
                    {
                    if(oligo[l] != mi->first[l])
                    {
                    mcount++;

                    if(mcount > OLIGO_MAX_MIS[a])
                    break;
                    }
                    }			

                    if(mcount <= OLIGO_MAX_MIS[a])
                    mi->second[mcount]++;

                    mi++;
                    } */
        //	}	

        omi++;
        counter++;
    }

    cerr << "\tdone" << endl;

    cerr << "Writing freq file: " << outfile.str();

    map<string, vector<unsigned int> >::iterator mi = m_freq.begin();

    //	cerr << "Oligo List:\n";	

    while(mi != m_freq.end())
    {
        out << mi->first;

        for(unsigned int t = 0; t <= OLIGO_MAX_MIS[a]; t++)
            out << '\t' << mi->second[t];	

        out << endl;

        mi++;
    }

    out.close();

    cerr << "\tdone" << endl << endl;

    return;
}

void count_ds_oligo(string *inputfile, vector<sequence> *SEQUENCE, short int a)
{
    ostringstream outfile;
    map<string, vector<unsigned int> >  m_freq; // (1 + OLIGO_MAX_MIS[a], map<string, unsigned int> ());
    string oligo(OLIGO_LENGTHS[a], 'N');

    //	outfile << *inputfile << "_" << OLIGO_LENGTHS[a] << "_ds.freq";
    outfile << scode << "_ds." << OLIGO_LENGTHS[a] << ".freq";

    ofstream out(outfile.str().c_str());

    cerr << "Building " << OLIGO_LENGTHS[a] << "mers list...";
    oligo_gen(a, 0, &m_freq, oligo);
    cerr << "\t\tdone" << endl;

    map<string, unsigned int> omap;

    cerr << "Finding exact occurrences...";

    for(unsigned int i = 0; i < SEQUENCE->size(); i++)
    {
        for(unsigned int j = 0; j < SEQUENCE->at(i).seq().size() - OLIGO_LENGTHS[a]; j++)
        {
            for(unsigned int k = 0; k < OLIGO_LENGTHS[a]; k++)
                oligo[k] = SEQUENCE->at(i).seq().at(j+k);	

            if(oligo.find("N") != string::npos)
                continue;

            string rc_oligo = string_rev_comp(&oligo);

            //	cout << oligo << endl << rc_oligo << endl << endl;

            if(rc_oligo == oligo)
            {
                if(omap.find(oligo) == omap.end())	
                    omap[oligo] = 1;
                else
                    omap[oligo] += 1;
            }

            else
                //			if(oligo != rc_oligo)
            {
                if(omap.find(oligo) == omap.end() && omap.find(rc_oligo) == omap.end())
                    omap[oligo] = 1;
                else if(omap.find(oligo) == omap.end())
                    omap[rc_oligo]++;
                else if(omap.find(rc_oligo) == omap.end())
                    omap[oligo]++;
            }			
        }
    }

    //	cerr << endl << "omap.find(AAAAATTT)->second = " << omap.find("AAAAATTT")->second << '\t' << omap.end()->second << endl;
    //	cerr << endl << "omap.find(AAATTTTT)->second = " << omap.find("AAATTTTT")->second << '\t' << omap.end()->second << endl;

    cerr << "\tdone" << endl;


    map<string, unsigned int>::iterator omi = omap.begin();
    unsigned int counter = 0;

    while(omi != omap.end())
    {
        cerr << (char)13 << "Counting " << OLIGO_LENGTHS[a] << "mers occurrences..." << counter+1 << '/' << omap.size();

        oligo = omi->first;

        map <string, vector<unsigned int> > oligos;

        oligos[oligo] = vector<unsigned int> ();

        near_oligo_gen(a,1,&oligo,&oligos);

        map<string, vector<unsigned int> >::iterator mi = oligos.begin(), mie = oligos.end();

        while(mi != mie)
        {
            //			if(omi->first == "AAATTTTT")
            //				cout << mi->first << '\t' << mi->second.size() << endl;

            oligos[string_rev_comp(&mi->first)] = mi->second;

            mi++;
        }

        /*		if(omi->first == "AAATTTTT")
                {
                mi = oligos.begin();

                while(mi != oligos.end())
                {
                cout << mi->first << '\t' << mi->second.size() << endl;
                mi++;
                }
                }*/

        mi = oligos.begin();

        while(mi != oligos.end())
        {
            //			if(mi->first == "TAATTT" || mi->first == "AAATTA")
            //				cout << endl << omi->first << '\t' << omi->second  << '\t' << mi->second.size() << '\t' << mi->first << endl;
            if(mi->first != string_rev_comp(&mi->first))	//QUI
                m_freq[mi->first][mi->second.size()] += omi->second;	//ORIGINAL
            else    //NEW
                m_freq[mi->first][mi->second.size()] += (omi->second * 2); //NEW

            mi++;
        }

        omi++;
        counter++;
    }

    cerr << "\tdone" << endl;

    cerr << "Writing freq file: " << outfile.str();

    map<string, vector<unsigned int> >::iterator mi = m_freq.begin();

    //	cerr << "Oligo List:\n";	

    while(mi != m_freq.end())
    {
        out << mi->first;

        for(unsigned int t = 0; t <= OLIGO_MAX_MIS[a]; t++)
            out << '\t' << mi->second[t];	

        out << endl;

        mi++;
    }

    out.close();

    cerr << "\tdone" << endl << endl;

    return;
}

void oligo_gen(short int a, short int lev, map<string, vector<unsigned int> > *freq, string oligo)
{
    for(unsigned int b = 0; b < ALPHABET_SIZE; b++)
    {
        oligo[lev] = ALPHABET[b];

        if(lev == OLIGO_LENGTHS[a] - 1)
            freq->insert(make_pair(oligo,vector<unsigned int> (OLIGO_MAX_MIS[a]+1,1)));
        else
            oligo_gen(a, lev+1, freq, oligo);	
    }

    return ;
}

void near_oligo_gen(short int a, short int lev,string *oligo, map<string, vector<unsigned int> > *oligos)
{
    map<string,vector<unsigned int> >::iterator mi = oligos->begin(), mie = oligos->end();

    while(lev <= OLIGO_MAX_MIS[a])
    {
        while(mi != mie)
        {
            if(mi->second.size() !=  lev-1)
            {
                mi++;
                continue;
            }
            for(unsigned int i = 0; i < OLIGO_LENGTHS[a]; i++)
            {
                bool flag = true;

                for(unsigned int k = 0; k < mi->second.size(); k++)
                {
                    if(mi->second.at(k) == i)
                    {
                        flag = false;
                        break;
                    }
                }

                if(!flag)
                    continue;

                for(unsigned int j = 0; j < ALPHABET_SIZE; j++)
                {
                    if(ALPHABET[j] == oligo->at(i))
                        continue;

                    string noligo = mi->first;

                    noligo[i] = ALPHABET[j];

                    vector<unsigned int> vbuf = mi->second;		

                    vbuf.push_back(i);

                    if(!DOUBLE_STRAND)
                        oligos->insert(make_pair(noligo, vbuf));
                    else
                    {
                        if(oligos->find(string_rev_comp(&noligo)) == oligos->end())
                            oligos->insert(make_pair(noligo, vbuf));
                    }
                }
            }

            mi++;
        }

        lev++;

        near_oligo_gen(a, lev, oligo, oligos);
    }

    return;

}

string string_rev_comp(const string *seq)
{
    int size = (int)seq->size();

    string rc_seq(size,'N');

    for(int i = 0; i < size; i++)
        rc_seq[size - 1 - i] = R_ALPHABET_TABLE[seq->at(i) - 'A'];

    return rc_seq;
}

sequence::sequence(string s)
{       
    istringstream str(s);
    string line;

    getline(str,line);

    NAME = line;

    while(getline(str,line))
        SEQ += line;

    GOOD = seq_cleaner();

    return;
} 

bool sequence::seq_cleaner()
{
    //        if(!SEQ.size())
    //              return false;
    for(int i = 0; i < OLIGO_TYPES; i++)
    {
        if(SEQ.size() < OLIGO_LENGTHS[i])
            return false;
    }

    for(int i = 0; i < SEQ.size(); i++)
    {
        bool flag = false;

        SEQ[i] = toupper(SEQ[i]);

        for(int a = 0; a < ALPHABET_SIZE; a++)
        {
            if(SEQ[i] == ALPHABET[a])
                flag = true;
        }

        if(!flag && SEQ[i] != 'N')
            SEQ[i] = 'N';
        //         return false;
    }

    return true;
}


string sequence::seq()
{
    return SEQ;
}

bool sequence::good()
{
    return GOOD;
}

