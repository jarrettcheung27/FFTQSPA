#if !defined(AFX_BCJRQSPA_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
#define AFX_BCJRQSPA_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include "stdafx.h"
#include "util.h"
#include "QaryLDPC.h"
#include "Random.h"
#include "Mapper.h"


class BCJRQSPA
{
public:
	BCJRQSPA();
	BCJRQSPA(string filename , int max_iteration, string mapping_name);
	~BCJRQSPA();

	int get_block_length() const;
	int get_total_length() const;
	void no_encoder(int *uu , int * cc);
	void re_encoder(int* uu, int* cc);
	void encoder4BiBo(int *uu , int * cc);
	void sym_encoder(int *uu);
	bool is_codeword(int *b_code);
	bool is_codeword_sym(int *b_code);
	int count_unsat(int *b_code);
	int count_unsat_sym(int *b_code);
	int m_bit_len_info;
	int m_bit_len_cc;
	

private:
	int parity_row;   //row number of the parity check matrix
	int parity_column;//column number of the parity check matrix
	int rank;         //rank of the parity check matrix
	                  //determine the check bit of the encoding process
	int *sym;
	int **parity_matrix;
	int **enc_parity_matrix;
	int d;

	int *inds;
	int *sym_inds;

	int max_iteration;

	int * row_weight_parity;
	int **row_location_parity;
	int * row_weight_enc;
	int **row_location_enc;

	//for interleaving
	double **_int_temp;
	int *pai;
	int *bit_pai;
	void InterleavingMessage(double **in, double **out, int col_len, int row_len);
	void DeInterleavingMessage(double **in, double **out, int col_len, int row_len);

public:	
	void InterleavingSymbol(int * , int * , int);
	void Init_Interleaver(int is_int);
	int q_ary;
	double power; 
	void Malloc(string filename, string mapping_name);
	CModem m_modem;


	Tanner_Graph tanner;	
	int m_degree;//GF(2^m)�е�m
	int bit_length_cc;//������BIT���ĳ���
	int m_len_xx;

	int FFTQSPA_FB(double *rr, int *uu, double sigma);
	int FFTQSPA(double* rr, int* uu, double sigma);
	int FFTQSPA4BiBo(double *rr, int *uu);
	int FFTQSPABIA(double *rr, int *uu, double sigma);
	int FFTQSPA_SYM(double *rr, int *uu, double sigma);
	int FFTQSPATrunc(double *rr, int *uu, double sigma);
	int HardDecision(double *rr, int *uu, double sigma);


	double *total_mul;
	double *temp_qary;
//for encoding and decoding
	int *tempP;
 	int **m_encH;
	int **m_decH;

	int sym_unsat_num;


//for message exchange
	char *bd;

	double **channel_for_spa;
	void Sys();
};

#endif // _MSC_VER > 1000
