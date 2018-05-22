  // Alphabet.cpp: implementation of the CAlphabet class.
//
//////////////////////////////////////////////////////////////////////

#include "Alphabet.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CAlphabet::CAlphabet()
{
	/*
		A=0,//Adenine
		C,	//Cytosine
		G,	//Guanine
		T,	//Thymine
	 *	R,	//A or G
		W,  //A or T
		S,  //C or G
		Y,	//C or T
		N,	//spacer
	 */
	relation[CAlphabet::A][0] = CAlphabet::A;
	relation[CAlphabet::C][0] = CAlphabet::C;
	relation[CAlphabet::G][0] = CAlphabet::G;
	relation[CAlphabet::T][0] = CAlphabet::T;

	relation[CAlphabet::R][0] = CAlphabet::A;
	relation[CAlphabet::R][1] = CAlphabet::G;

	relation[CAlphabet::W][0] = CAlphabet::A;
	relation[CAlphabet::W][1] = CAlphabet::T;

	relation[CAlphabet::S][0] = CAlphabet::C;
	relation[CAlphabet::S][1] = CAlphabet::G;

	relation[CAlphabet::Y][0] = CAlphabet::C;
	relation[CAlphabet::Y][1] = CAlphabet::T;
	
	relation[CAlphabet::N][0] = CAlphabet::A;
	relation[CAlphabet::N][1] = CAlphabet::C;
	relation[CAlphabet::N][2] = CAlphabet::G;
	relation[CAlphabet::N][3] = CAlphabet::T;

	extendnum[CAlphabet::A] = 1;
	extendnum[CAlphabet::C] = 1;
	extendnum[CAlphabet::G] = 1;
	extendnum[CAlphabet::T] = 1;
	extendnum[CAlphabet::R] = 2;
	extendnum[CAlphabet::W] = 2;
	extendnum[CAlphabet::S] = 2;
	extendnum[CAlphabet::Y] = 2;
	extendnum[CAlphabet::N] = 4;

	compatible[CAlphabet::A][CAlphabet::A] = true;
	compatible[CAlphabet::A][CAlphabet::C] = false;
	compatible[CAlphabet::A][CAlphabet::G] = false;
	compatible[CAlphabet::A][CAlphabet::T] = false;
	compatibleindex[CAlphabet::A][CAlphabet::A] = 0;

	compatible[CAlphabet::C][CAlphabet::A] = false;
	compatible[CAlphabet::C][CAlphabet::C] = true;
	compatible[CAlphabet::C][CAlphabet::G] = false;
	compatible[CAlphabet::C][CAlphabet::T] = false;
	compatibleindex[CAlphabet::C][CAlphabet::C] = 0;

	compatible[CAlphabet::G][CAlphabet::A] = false;
	compatible[CAlphabet::G][CAlphabet::C] = false;
	compatible[CAlphabet::G][CAlphabet::G] = true;
	compatible[CAlphabet::G][CAlphabet::T] = false;
	compatibleindex[CAlphabet::G][CAlphabet::G] = 0;

	compatible[CAlphabet::T][CAlphabet::A] = false;
	compatible[CAlphabet::T][CAlphabet::C] = false;
	compatible[CAlphabet::T][CAlphabet::G] = false;
	compatible[CAlphabet::T][CAlphabet::T] = true;
	compatibleindex[CAlphabet::T][CAlphabet::T] = 0;

	compatible[CAlphabet::R][CAlphabet::A] = true;
	compatible[CAlphabet::R][CAlphabet::C] = false;
	compatible[CAlphabet::R][CAlphabet::G] = true;
	compatible[CAlphabet::R][CAlphabet::T] = false;
	compatibleindex[CAlphabet::R][CAlphabet::A] = 0;
	compatibleindex[CAlphabet::R][CAlphabet::G] = 1;

	compatible[CAlphabet::W][CAlphabet::A] = true;
	compatible[CAlphabet::W][CAlphabet::C] = false;
	compatible[CAlphabet::W][CAlphabet::G] = false;
	compatible[CAlphabet::W][CAlphabet::T] = true;
	compatibleindex[CAlphabet::W][CAlphabet::A] = 0;
	compatibleindex[CAlphabet::W][CAlphabet::T] = 1;

	compatible[CAlphabet::S][CAlphabet::A] = false;
	compatible[CAlphabet::S][CAlphabet::C] = true;
	compatible[CAlphabet::S][CAlphabet::G] = true;
	compatible[CAlphabet::S][CAlphabet::T] = false;
	compatibleindex[CAlphabet::S][CAlphabet::C] = 0;
	compatibleindex[CAlphabet::S][CAlphabet::G] = 1;

	compatible[CAlphabet::Y][CAlphabet::A] = false;
	compatible[CAlphabet::Y][CAlphabet::C] = true;
	compatible[CAlphabet::Y][CAlphabet::G] = false;
	compatible[CAlphabet::Y][CAlphabet::T] = true;
	compatibleindex[CAlphabet::Y][CAlphabet::C] = 0;
	compatibleindex[CAlphabet::Y][CAlphabet::T] = 1;

	compatible[CAlphabet::N][CAlphabet::A] = true;
	compatible[CAlphabet::N][CAlphabet::C] = true;
	compatible[CAlphabet::N][CAlphabet::G] = true;
	compatible[CAlphabet::N][CAlphabet::T] = true;
	compatibleindex[CAlphabet::N][CAlphabet::A] = 0;
	compatibleindex[CAlphabet::N][CAlphabet::C] = 1;
	compatibleindex[CAlphabet::N][CAlphabet::G] = 2;
	compatibleindex[CAlphabet::N][CAlphabet::T] = 3;
}

CAlphabet::~CAlphabet()
{

}

void CAlphabet::printalp(CAlphabet::ALPHABET a)
{
	switch(a) {
	case 0:
		printf("A");
		break;
	case 1:
		printf("C");
		break;
	case 2:
		printf("G");
		break;
	case 3:
		printf("T");
		break;
	case 4:
		printf("R");
		break;
	case 5:
		printf("W");
		break;
	case 6:
		printf("S");
		break;
	case 7:
		printf("Y");
		break;
	case 8:
		printf("N");
		break;
	default:
		printf("?");
	}
}

bool CAlphabet::checkValid(int type, char one)
{
	switch(one) {
	case 'A':
	case 'C':
	case 'G':
	case 'T':
	case 'a':
	case 'c':
	case 'g':
	case 't':
		return true;
	case 'R':
	case 'W':
	case 'S':
	case 'Y':
	case 'N':
	case 'r':
	case 'w':
	case 's':
	case 'y':
	case 'n':
		if (type==2) {
			return true;
		} else
			return false;
	default:
		return false;
	}
}

int CAlphabet::getIndex(char one)
{
	switch(one) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	default:
		return -1;
	}
}
