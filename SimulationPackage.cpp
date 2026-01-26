//Writer : Terrence Zhao
//Time   :2009.04.08
//Function:Simulation
#include "stdafx.h"
#include "Random.h"
#include "Sourcesink.h"
#include "util.h"
#include "BCJRQSPA.h"
#include "QaryLDPC.h"
#include "Qary_Gauss_E.h"
using namespace std;

CLCRandNum rndGen0;
CWHRandNum rndGen1;
int len_pad_sym = 5;
extern CFiniteField2 GF;

void simulation_qary();
void simulation_qary_bibo();
void simulation_qary_tb();
void simulation_qary_piggypacking();
void simulation_qary_bpsk();
void simulation_qary_sym();
void set_flag(int &flag , int number);
/*
void simulation_qary_isi();
void simulation_qary_isi_T();
void simulation_qary_isi_inteleaver();
*/
int _tmain(int argc, _TCHAR* argv[])
{
	simulation_qary_bibo();
	//simulation_qary_sym();
	//simulation_qary_tb();
	//simulation_qary_piggypacking();
	//simulation_qary();
	//simulation_qary_isi();
}

void set_flag(int &flag , int number)
{
	std::cout<<"Please input the seed for the "<<number<<"  random number generator : "<<std::endl;
	std::cout<<"   Notice : default(-1) ; random(0) ; setseeds(1)"<<std::endl<<"   Flag = : ";
	std::cin>>flag;
	return;
}


void simulation_qary_tb()
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma;

	FILE *fp_0, *fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_", filename_orth = "";
	double **orth;

    flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
    //set_flag(flag , 2);
	rndGen1.SetSeed(flag);
	
	std::cout<<"Please input the number of the input file  :  "<<std::endl<<">:";
	file_num = "2013";
	//std::cin>>file_num;
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if( !inFile)
	{
		std::cerr<<"Can't open the parameter input file : "<<file_name<<std::endl;
	}
//Malloc
	char * temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	low_snr = atof(temp_value);

    inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	step = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename = filename + temp_string ;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_mapping = filename_mapping + temp_string;
	
	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_orth = filename_orth + temp_string;

	inFile.close();	delete []temp_string; delete []temp_value;
//Source and Sink
	CSourceSink *source_sink = new CSourceSink();
	CSourceSink *source_sink_ber = new CSourceSink();

	BCJRQSPA *ij = new BCJRQSPA(filename , max_iteration , filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0;

	int * uu = new int[total_length];//orginal symbol string
	int * uu_hat = new int[total_length];//symbol string after decoding
	int * cc = new int[total_length];//codeword
	int * int_cc = new int[total_length];//interleaved codeword

	for(int i = 0 ; i < total_length ; i ++)
	{
		cc[i] = 0;int_cc[i] = 0;
	}
	int m_len_xx = ij->m_len_xx, inform_len = total_length * ij->m_degree;
	int *bit = new int[inform_len];
	int *bit_hat = new int[inform_len];
	double * rr = new double[m_len_xx];//output of Mapper
	double * yy = new double[m_len_xx];//output of AWGN
	double * yy_t = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	orth = new double*[inform_len];
	for(int i = 0 ; i < inform_len ; i ++)
	{
		orth[i] = new double[inform_len];
	}
	std::ifstream oFile_orth("OrthMatrix.txt" , ios_base::in); 
	for(int i = 0 ; i < inform_len ; i ++)
	{
		for(int j = 0 ; j < inform_len ; j ++)
			oFile_orth>>orth[i][j];
	}


	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt" , ios_base::app); 
	std::ofstream oFile_ser("ser.txt" , ios_base::app); 
	std::ofstream iter("iter.txt" , ios_base::app); 
	std::ofstream oFile_fer("fer.txt" , ios_base::app);
	std::ofstream oFile_err("Error_Pattern.txt" , ios_base::app); 
	std::ofstream oFile_data("Error_Data.txt" , ios_base::app); 
	if((! oFile_ber) || (! oFile_fer))
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}
	oFile_ber<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ber<<"  SNR                         BER "<<std::endl;
	oFile_fer<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_fer<<"  SNR                         FER "<<std::endl;
	oFile_ser<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ser<<"  SNR                         SER "<<std::endl;

	std::cout<<endl<<"----------------------SIMULATION BEGIN----------------------"<<endl;

	double rate = (double)length / (double)total_length;
	for(double snr = low_snr ; snr < max_snr ; snr = snr + step)
	{
		//snr = 2.0;
		oFile_err<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		oFile_data<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		var = pow(10.0, 0.1 * snr); 
		sigma = sqrt(1 / (2.0 * var * rate));//Cal channel parameters
		//sigma = sqrt(1 / (2.0 * var * ij->m_modem.m_len_signal));//Cal channel parameters
		source_sink->ClrCnt();source_sink_ber->ClrCnt();BER = 0.0;FER = 0.0;total_iter = 0;//clear source and sink
		while(source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetSymStr(uu, q_ary, length);//source
			ij->re_encoder(uu,cc);//encoding
			ij->m_modem.Mapping(cc, rr, total_length);//mapping
			rndGen1.Normal(yy,m_len_xx);//Gaussian noise
			for(int i = 0 ; i < m_len_xx ; i ++)
			{
				yy_t[i] = 0.0;
				for(int j = 0 ; j < m_len_xx ; j ++)
				{
					yy_t[i] += yy[j] * orth[i][j];
				}
			}
			for(int i = 0 ; i < m_len_xx; i ++)
			{
				yy[i] = rr[i] + sigma * yy_t[i];//AWGN
			}
			total_iter += ij->FFTQSPA(yy ,uu_hat, sigma);//soft decoding
//destination
			for(int i = 0 ; i < total_length ; i ++)
			{
				for(int ii = 0 ; ii < ij->m_degree ; ii ++)
				{
					bit[i*ij->m_degree + ii] = (cc[i]>>ii)&1;
					bit_hat[i*ij->m_degree + ii] = (uu_hat[i]>>ii)&1;
				}
			}
			source_sink_ber->CntErr(bit, bit_hat, inform_len, 1);
			source_sink->CntErr_q(cc, uu_hat, total_length, 1);
//error pattern
			count = 0;
			for(int i = 0 ; i < total_length ; i ++)
			{
				uu_hat[i] = GF.Add(cc[i], uu_hat[i]);
				if(uu_hat[i])
					count += 1;
			}
			if(count)
			{
//output the codeword and the received signal
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<cc[i]<<" ";
				}
				oFile_data<<endl;
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<yy[i]<<" ";
				}
				oFile_data<<endl;
//output the error pattern
				oFile_err<<count<<"  "<<ij->count_unsat(uu_hat)<<"     ";
				for(int i = 0 ; i < total_length ; i ++)
				{
					if(uu_hat[i])
					{
						oFile_err<<i<<" "<<uu_hat[i]<<"   ";
					}
				}
				oFile_err<<std::endl;
			}
			if((int)source_sink->_tot_blk_num % 1 == 0)
			{		
				SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
				BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
				FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
				std::cout<<"Frame number: "<<source_sink->_tot_blk_num<<" Error frame number : "<<source_sink->_err_blk_num<<std::endl;
				std::cout<<"SNR          BER            SER          FER"<<std::endl;
				std::cout<<snr<<"          "<<BER<<"          "<<SER<<"           "<<FER<<std::endl;
				std::cout<<"Average iteration numbers is "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0);rndGen1.PrintState(fp_1);fclose(fp_0);fclose(fp_1);
			}

		}
		BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
		SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
		FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
		oFile_ber<<snr<<"                 "<<BER<<std::endl;
        oFile_fer<<snr<<"                 "<<FER<<std::endl;
        oFile_ser<<snr<<"                 "<<SER<<std::endl;
        iter<<snr<<"                 "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
	}
//clear
	oFile_ber.close();oFile_fer.close();oFile_ser.close();oFile_err.close();iter.close();oFile_data.close();
	delete source_sink;delete []uu;delete []cc;delete []int_cc;delete []rr;delete []yy;delete []yy_t;delete ij;delete []bit;delete []bit_hat;
	for(int i = 0 ; i < inform_len ; i ++)
	{
		delete []orth[i];
	}
	delete orth;
}


void simulation_qary_piggypacking()
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma;

	FILE *fp_0, *fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_", filename_orth = "";
	double **orth;

    flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
    //set_flag(flag , 2);
	rndGen1.SetSeed(flag);
	
	std::cout<<"Please input the number of the input file  :  "<<std::endl<<">:";
	//file_num = "2013";
	std::cin>>file_num;
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if( !inFile)
	{
		std::cerr<<"Can't open the parameter input file : "<<file_name<<std::endl;
	}
//Malloc
	char * temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	low_snr = atof(temp_value);

    inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	step = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename = filename + temp_string ;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_mapping = filename_mapping + temp_string;
	
	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_orth = filename_orth + temp_string;

	inFile.close();	delete []temp_string; delete []temp_value;
//Source and Sink
	CSourceSink *source_sink = new CSourceSink();
	CSourceSink *source_sink_ber = new CSourceSink();

	BCJRQSPA *ij = new BCJRQSPA(filename , max_iteration , filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0, max_error=0;

	int * uu = new int[total_length];//orginal symbol string
	int * uu_hat = new int[total_length];//symbol string after decoding
	int * cc = new int[total_length];//codeword
	int * int_cc = new int[total_length];//interleaved codeword

	for(int i = 0 ; i < total_length ; i ++)
	{
		cc[i] = 0;int_cc[i] = 0;
	}
	int m_len_xx = ij->m_len_xx, inform_len = total_length * ij->m_degree;
	int *bit = new int[inform_len];
	int *bit_hat = new int[inform_len];
	double * rr = new double[m_len_xx];//output of Mapper
	double * yy = new double[m_len_xx];//output of AWGN
	double * yy_t = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	orth = new double*[inform_len];
	for(int i = 0 ; i < inform_len ; i ++)
	{
		orth[i] = new double[inform_len];
	}
	/*std::ifstream oFile_orth("OrthMatrix.txt" , ios_base::in); 
	for(int i = 0 ; i < inform_len ; i ++)
	{
		for(int j = 0 ; j < inform_len ; j ++)
			oFile_orth>>orth[i][j];
	}*/


	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt" , ios_base::app); 
	std::ofstream oFile_ser("ser.txt" , ios_base::app); 
	std::ofstream iter("iter.txt" , ios_base::app); 
	std::ofstream oFile_fer("fer.txt" , ios_base::app);
	std::ofstream oFile_err("Error_Pattern.txt" , ios_base::app); 
	std::ofstream oFile_data("Error_Data.txt" , ios_base::app); 
	if((! oFile_ber) || (! oFile_fer))
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}
	oFile_ber<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ber<<"  SNR                         BER "<<std::endl;
	oFile_fer<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_fer<<"  SNR                         FER "<<std::endl;
	oFile_ser<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ser<<"  SNR                         SER "<<std::endl;

	std::cout<<endl<<"----------------------SIMULATION BEGIN----------------------"<<endl;

	double rate = (double)length / (double)total_length;
	for(double snr = low_snr ; snr < max_snr ; snr = snr + step)
	{
		//snr = 2.0;
		oFile_err<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		oFile_data<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		var = pow(10.0, 0.1 * snr); 
		sigma = sqrt(1 / (2.0 * var * rate));//Cal channel parameters
		//sigma = sqrt(1 / (2.0 * var * ij->m_modem.m_len_signal));//Cal channel parameters
		source_sink->ClrCnt();source_sink_ber->ClrCnt();BER = 0.0;FER = 0.0;total_iter = 0;max_error=100000;//clear source and sink
		while(source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetSymStr(uu, q_ary, length);//source
			ij->re_encoder(uu,cc);//encoding
			ij->m_modem.Mapping(cc, rr, total_length);//mapping
			rndGen1.Normal(yy,m_len_xx);//Gaussian noise

			for(int i = 0 ; i < m_len_xx; i ++)
			{
				yy[i] = rr[i] + sigma * yy[i];//AWGN
			}
			total_iter = ij->HardDecision(yy ,uu_hat, sigma);//soft decoding
			if (max_error > total_iter)
			{
				max_error = total_iter;
			}
//destination
			for(int i = 0 ; i < total_length ; i ++)
			{
				for(int ii = 0 ; ii < ij->m_degree ; ii ++)
				{
					bit[i*ij->m_degree + ii] = (cc[i]>>ii)&1;
					bit_hat[i*ij->m_degree + ii] = (uu_hat[i]>>ii)&1;
				}
			}
			source_sink_ber->CntErr(bit, bit_hat, inform_len, 1);
			source_sink->CntErr_q(cc, uu_hat, total_length, 1);
//error pattern
			count = 0;
			for(int i = 0 ; i < total_length ; i ++)
			{
				uu_hat[i] = GF.Add(cc[i], uu_hat[i]);
				if(uu_hat[i])
					count += 1;
			}
			if(count)
			{
//output the codeword and the received signal
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<cc[i]<<" ";
				}
				oFile_data<<endl;
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<yy[i]<<" ";
				}
				oFile_data<<endl;
//output the error pattern
				oFile_err<<count<<"  "<<ij->count_unsat(uu_hat)<<"     ";
				for(int i = 0 ; i < total_length ; i ++)
				{
					if(uu_hat[i])
					{
						oFile_err<<i<<" "<<uu_hat[i]<<"   ";
					}
				}
				oFile_err<<std::endl;
			}
			if((int)source_sink->_tot_blk_num % 1 == 0)
			{		
				SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
				BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
				FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
				std::cout<<"Frame number: "<<source_sink->_tot_blk_num<<" Error frame number : "<<source_sink->_err_blk_num<<std::endl;
				std::cout<<"SNR          BER            SER          FER"<<std::endl;
				std::cout<<snr<<"          "<<BER<<"          "<<SER<<"           "<<FER<<std::endl;
				std::cout<<"Average iteration numbers is "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
				std::cout<<"Average iteration numbers is "<<max_error<<std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0);rndGen1.PrintState(fp_1);fclose(fp_0);fclose(fp_1);
			}

		}
		BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
		SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
		FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
		oFile_ber<<snr<<"                 "<<BER<<std::endl;
        oFile_fer<<snr<<"                 "<<FER<<std::endl;
        oFile_ser<<snr<<"                 "<<SER<<std::endl;
        iter<<snr<<"                 "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
	}
//clear
	oFile_ber.close();oFile_fer.close();oFile_ser.close();oFile_err.close();iter.close();oFile_data.close();
	delete source_sink;delete []uu;delete []cc;delete []int_cc;delete []rr;delete []yy;delete []yy_t;delete ij;delete []bit;delete []bit_hat;
	for(int i = 0 ; i < inform_len ; i ++)
	{
		delete []orth[i];
	}
	delete orth;
}

/// <summary>
/// //////////////
/// </summary>
void simulation_qary_bpsk()
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma;

	FILE *fp_0, *fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_", filename_orth = "";
	double **orth;

    flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
    //set_flag(flag , 2);
	rndGen1.SetSeed(flag);
	
	std::cout<<"Please input the number of the input file  :  "<<std::endl<<">:";
	//file_num = "2013";
	std::cin>>file_num;
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if( !inFile)
	{
		std::cerr<<"Can't open the parameter input file : "<<file_name<<std::endl;
	}
//Malloc
	char * temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	low_snr = atof(temp_value);

    inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	step = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename = filename + temp_string ;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_mapping = filename_mapping + temp_string;
	
	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_orth = filename_orth + temp_string;

	inFile.close();	delete []temp_string; delete []temp_value;
//Source and Sink
	CSourceSink *source_sink = new CSourceSink();
	CSourceSink *source_sink_ber = new CSourceSink();

	BCJRQSPA *ij = new BCJRQSPA(filename , max_iteration , filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0;

	int * uu = new int[total_length];//orginal symbol string
	int * uu_hat = new int[total_length];//symbol string after decoding
	int * cc = new int[total_length];//codeword
	int * int_cc = new int[total_length];//interleaved codeword

	for(int i = 0 ; i < total_length ; i ++)
	{
		cc[i] = 0;int_cc[i] = 0;
	}
	int m_len_xx = ij->m_len_xx, inform_len = total_length * ij->m_degree;
	int *bit = new int[inform_len];
	int *bit_hat = new int[inform_len];
	double * rr = new double[m_len_xx];//output of Mapper
	double * yy = new double[m_len_xx];//output of AWGN
	double * yy_t = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	orth = new double*[inform_len];
	for(int i = 0 ; i < inform_len ; i ++)
	{
		orth[i] = new double[inform_len];
	}
	/*std::ifstream oFile_orth("OrthMatrix.txt" , ios_base::in); 
	for(int i = 0 ; i < inform_len ; i ++)
	{
		for(int j = 0 ; j < inform_len ; j ++)
			oFile_orth>>orth[i][j];
	}*/


	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt" , ios_base::app); 
	std::ofstream oFile_ser("ser.txt" , ios_base::app); 
	std::ofstream iter("iter.txt" , ios_base::app); 
	std::ofstream oFile_fer("fer.txt" , ios_base::app);
	std::ofstream oFile_err("Error_Pattern.txt" , ios_base::app); 
	std::ofstream oFile_data("Error_Data.txt" , ios_base::app); 
	if((! oFile_ber) || (! oFile_fer))
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}
	oFile_ber<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ber<<"  SNR                         BER "<<std::endl;
	oFile_fer<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_fer<<"  SNR                         FER "<<std::endl;
	oFile_ser<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ser<<"  SNR                         SER "<<std::endl;

	std::cout<<endl<<"----------------------SIMULATION BEGIN----------------------"<<endl;

	double rate = (double)length / (double)total_length;
	for(double snr = low_snr ; snr < max_snr ; snr = snr + step)
	{
		//snr = 2.0;
		oFile_err<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		oFile_data<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		var = pow(10.0, 0.1 * snr); 
		sigma = sqrt(1 / (2.0 * var * rate));//Cal channel parameters
		//sigma = sqrt(1 / (2.0 * var * ij->m_modem.m_len_signal));//Cal channel parameters
		source_sink->ClrCnt();source_sink_ber->ClrCnt();BER = 0.0;FER = 0.0;total_iter = 0;//clear source and sink
		while(source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetSymStr(uu, q_ary, length);//source
			ij->re_encoder(uu,cc);//encoding
			ij->m_modem.Mapping(cc, rr, total_length);//mapping
			rndGen1.Normal(yy,m_len_xx);//Gaussian noise
			/*for(int i = 0 ; i < m_len_xx ; i ++)
			{
				yy_t[i] = 0.0;
				for(int j = 0 ; j < m_len_xx ; j ++)
				{
					yy_t[i] += yy[j] * orth[i][j];
				}
			}*/
			for(int i = 0 ; i < m_len_xx; i ++)
			{
				yy[i] = rr[i] + sigma * yy[i];//AWGN
			}
			total_iter += ij->FFTQSPA(yy ,uu_hat, sigma);//soft decoding
//destination
			for(int i = 0 ; i < total_length ; i ++)
			{
				for(int ii = 0 ; ii < ij->m_degree ; ii ++)
				{
					bit[i*ij->m_degree + ii] = (cc[i]>>ii)&1;
					bit_hat[i*ij->m_degree + ii] = (uu_hat[i]>>ii)&1;
				}
			}
			source_sink_ber->CntErr(bit, bit_hat, inform_len, 1);
			source_sink->CntErr_q(cc, uu_hat, total_length, 1);
//error pattern
			count = 0;
			for(int i = 0 ; i < total_length ; i ++)
			{
				uu_hat[i] = GF.Add(cc[i], uu_hat[i]);
				if(uu_hat[i])
					count += 1;
			}
			if(count)
			{
//output the codeword and the received signal
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<cc[i]<<" ";
				}
				oFile_data<<endl;
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<yy[i]<<" ";
				}
				oFile_data<<endl;
//output the error pattern
				oFile_err<<count<<"  "<<ij->count_unsat(uu_hat)<<"     ";
				for(int i = 0 ; i < total_length ; i ++)
				{
					if(uu_hat[i])
					{
						oFile_err<<i<<" "<<uu_hat[i]<<"   ";
					}
				}
				oFile_err<<std::endl;
			}
			if((int)source_sink->_tot_blk_num % 1 == 0)
			{		
				SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
				BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
				FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
				std::cout<<"Frame number: "<<source_sink->_tot_blk_num<<" Error frame number : "<<source_sink->_err_blk_num<<std::endl;
				std::cout<<"SNR          BER            SER          FER"<<std::endl;
				std::cout<<snr<<"          "<<BER<<"          "<<SER<<"           "<<FER<<std::endl;
				std::cout<<"Average iteration numbers is "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0);rndGen1.PrintState(fp_1);fclose(fp_0);fclose(fp_1);
			}

		}
		BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
		SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
		FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
		oFile_ber<<snr<<"                 "<<BER<<std::endl;
        oFile_fer<<snr<<"                 "<<FER<<std::endl;
        oFile_ser<<snr<<"                 "<<SER<<std::endl;
        iter<<snr<<"                 "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
	}
//clear
	oFile_ber.close();oFile_fer.close();oFile_ser.close();oFile_err.close();iter.close();oFile_data.close();
	delete source_sink;delete []uu;delete []cc;delete []int_cc;delete []rr;delete []yy;delete []yy_t;delete ij;delete []bit;delete []bit_hat;
	for(int i = 0 ; i < inform_len ; i ++)
	{
		delete []orth[i];
	}
	delete orth;
}




void simulation_qary_sym()
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma;
	FILE *fp_0, *fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_";

    flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
    //set_flag(flag , 2);
	rndGen1.SetSeed(flag);
	
	std::cout<<"Please input the number of the input file  :  "<<std::endl<<">:";
	//std::cin>>file_num;
	file_num = "2013";
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if( !inFile)
	{
		std::cerr<<"Can't open the parameter input file : "<<file_name<<std::endl;
	}
//Malloc
	char * temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	low_snr = atof(temp_value);

    inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	step = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename = filename + temp_string ;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_mapping = filename_mapping + temp_string;
	inFile.close();	delete []temp_string; delete []temp_value;
//Source and Sink
	CSourceSink *source_sink = new CSourceSink();
	CSourceSink *source_sink_ber = new CSourceSink();

	BCJRQSPA *ij = new BCJRQSPA(filename , max_iteration , filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0;

	int * uu = new int[length];//orginal symbol string
	int * uu_hat = new int[total_length];//symbol string after decoding
	int * cc = new int[total_length];//codeword
	int * int_cc = new int[total_length];//interleaved codeword
	for(int i = 0 ; i < total_length ; i ++)
	{
		cc[i] = 0;int_cc[i] = 0;
	}
	int m_len_xx = ij->m_len_xx, inform_len = total_length * ij->m_degree;
	int *bit = new int[inform_len];
	int *bit_hat = new int[inform_len];
	double * rr = new double[m_len_xx];//output of Mapper
	double * yy = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt" , ios_base::app); 
	std::ofstream oFile_ser("ser.txt" , ios_base::app); 
	std::ofstream iter("iter.txt" , ios_base::app); 
	std::ofstream oFile_fer("fer.txt" , ios_base::app);
	std::ofstream oFile_err("Error_Pattern.txt" , ios_base::app); 
	std::ofstream oFile_data("Error_Data.txt" , ios_base::app); 
	if((! oFile_ber) || (! oFile_fer))
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}
	oFile_ber<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ber<<"  SNR                         BER "<<std::endl;
	oFile_fer<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_fer<<"  SNR                         FER "<<std::endl;
	oFile_ser<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ser<<"  SNR                         SER "<<std::endl;

	double rate = (double)length / (double)total_length;
	for(double snr = low_snr ; snr < max_snr ; snr = snr + step)
	{
		oFile_err<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		oFile_data<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		var = pow(10.0, 0.1 * snr); 
		sigma = sqrt(1 / (2.0 * var * rate * ij->m_modem.m_len_signal));//Cal channel parameters
		source_sink->ClrCnt();source_sink_ber->ClrCnt();BER = 0.0;FER = 0.0;total_iter = 0;//clear source and sink
		while(source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetSymStr(uu, q_ary, total_length);//source
			ij->sym_encoder(uu);//encoding
			ij->m_modem.Mapping(uu, rr, total_length);//mapping
			rndGen1.Normal(yy,m_len_xx);//Gaussian noise
			for(int i = 0 ; i < m_len_xx; i ++)
			{
				yy[i] = rr[i] + sigma * yy[i];//AWGN
			}
			total_iter += ij->FFTQSPA_SYM(yy ,uu_hat, sigma);//soft decoding
//destination
			for(int i = 0 ; i < total_length ; i ++)
			{
				for(int ii = 0 ; ii < ij->m_degree ; ii ++)
				{
					bit[i*ij->m_degree + ii] = (uu[i]>>ii)&1;
					bit_hat[i*ij->m_degree + ii] = (uu_hat[i]>>ii)&1;
				}
			}
			source_sink_ber->CntErr(bit, bit_hat, inform_len, 1);
			source_sink->CntErr_q(uu, uu_hat, total_length, 1);
//error pattern
			count = 0;
			for(int i = 0 ; i < total_length ; i ++)
			{
				uu_hat[i] = GF.Add(uu[i], uu_hat[i]);
				if(uu_hat[i])
					count += 1;
			}
			if(count)
			{
//output the codeword and the received signal
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<uu[i]<<" ";
				}
				oFile_data<<endl;
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<yy[i]<<" ";
				}
				oFile_data<<endl;
//output the error pattern
				oFile_err<<count<<"  "<<ij->count_unsat(uu_hat)<<"     ";
				for(int i = 0 ; i < total_length ; i ++)
				{
					if(uu_hat[i])
					{
						oFile_err<<i<<" "<<uu_hat[i]<<"   ";
					}
				}
				oFile_err<<std::endl;
			}
			if((int)source_sink->_tot_blk_num % 100 == 0)
			{		
				SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
				BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
				FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
				std::cout<<"Frame number: "<<source_sink->_tot_blk_num<<" Error frame number : "<<source_sink->_err_blk_num<<std::endl;
				std::cout<<"SNR          BER            SER          FER"<<std::endl;
				std::cout<<snr<<"          "<<BER<<"          "<<SER<<"           "<<FER<<std::endl;
				std::cout<<"Average iteration numbers is "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0);rndGen1.PrintState(fp_1);fclose(fp_0);fclose(fp_1);
			}

		}
		BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
		SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
		FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
		oFile_ber<<snr<<"                 "<<BER<<std::endl;
        oFile_fer<<snr<<"                 "<<FER<<std::endl;
        oFile_ser<<snr<<"                 "<<SER<<std::endl;
        iter<<snr<<"                 "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
	}
//clear
	oFile_ber.close();oFile_fer.close();oFile_ser.close();oFile_err.close();iter.close();oFile_data.close();
	delete source_sink;delete []uu;delete []cc;delete []int_cc;delete []rr;delete []yy;delete ij;delete []bit;delete []bit_hat;
}



void simulation_qary()
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma;

	FILE *fp_0, *fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_", filename_orth = "";
	double **orth;

    flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
    //set_flag(flag , 2);
	rndGen1.SetSeed(flag);
	
	std::cout<<"Please input the number of the input file  :  "<<std::endl<<">:";
	file_num = "2013";
	//std::cin>>file_num;
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if( !inFile)
	{
		std::cerr<<"Can't open the parameter input file : "<<file_name<<std::endl;
	}
//Malloc
	char * temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	low_snr = atof(temp_value);

    inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	step = atof(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_value,1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename = filename + temp_string ;

	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_mapping = filename_mapping + temp_string;
	
	inFile.getline(temp_string,1024);
	inFile.getline(temp_string,1024);
	filename_orth = filename_orth + temp_string;

	inFile.close();	delete []temp_string; delete []temp_value;
//Source and Sink
	CSourceSink *source_sink = new CSourceSink();
	CSourceSink *source_sink_ber = new CSourceSink();

	BCJRQSPA *ij = new BCJRQSPA(filename , max_iteration , filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0;

	int * uu = new int[total_length];//orginal symbol string
	int * uu_hat = new int[total_length];//symbol string after decoding
	int * cc = new int[total_length];//codeword
	int * int_cc = new int[total_length];//interleaved codeword

	for(int i = 0 ; i < total_length ; i ++)
	{
		cc[i] = 0;int_cc[i] = 0;
	}
	int m_len_xx = ij->m_len_xx, inform_len = total_length * ij->m_degree;
	int *bit = new int[inform_len];
	int *bit_hat = new int[inform_len];
	double * rr = new double[m_len_xx];//output of Mapper
	double * yy = new double[m_len_xx];//output of AWGN
	double * yy_t = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	orth = new double*[inform_len];
	for(int i = 0 ; i < inform_len ; i ++)
	{
		orth[i] = new double[inform_len];
	}
	std::ifstream oFile_orth("OrthMatrix.txt" , ios_base::in); 
	for(int i = 0 ; i < inform_len ; i ++)
	{
		for(int j = 0 ; j < inform_len ; j ++)
			oFile_orth>>orth[i][j];
	}


	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt" , ios_base::app); 
	std::ofstream oFile_ser("ser.txt" , ios_base::app); 
	std::ofstream iter("iter.txt" , ios_base::app); 
	std::ofstream oFile_fer("fer.txt" , ios_base::app);
	std::ofstream oFile_err("Error_Pattern.txt" , ios_base::app); 
	std::ofstream oFile_data("Error_Data.txt" , ios_base::app); 
	if((! oFile_ber) || (! oFile_fer))
	{
		std::cerr<<"Can't open the file ber.txt or the file fer.txt !Please check it !"<<std::endl;
		exit(-1);
	}
	oFile_ber<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ber<<"  SNR                         BER "<<std::endl;
	oFile_fer<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_fer<<"  SNR                         FER "<<std::endl;
	oFile_ser<<"The following result corresponding the input file : "<<file_name<<std::endl;
	oFile_ser<<"  SNR                         SER "<<std::endl;

	std::cout<<endl<<"----------------------SIMULATION BEGIN----------------------"<<endl;

	double rate = (double)length / (double)total_length;
	for(double snr = low_snr ; snr < max_snr ; snr = snr + step)
	{
		oFile_err<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		oFile_data<<"The following result corresponding the SNR point : "<<snr<<" dB."<<std::endl;
		var = pow(10.0, 0.1 * snr); 
		sigma = sqrt(1 / (2.0 * var * rate * ij->m_degree));//Cal channel parameters
		//sigma = sqrt(1 / (2.0 * var * ij->m_modem.m_len_signal));//Cal channel parameters
		source_sink->ClrCnt();source_sink_ber->ClrCnt();BER = 0.0;FER = 0.0;total_iter = 0;//clear source and sink
		while(source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetSymStr(uu, q_ary, length);//source
			for(int i = 0 ; i < length ; i ++)
				uu[i] = 0;
			ij->re_encoder(uu,cc);//encoding
			ij->m_modem.Mapping(cc, rr, total_length);//mapping
			rndGen1.Normal(yy,m_len_xx);//Gaussian noise
			yy[0] = 40 * sigma;
			for(int i = 0 ; i < m_len_xx ; i ++)
			{
				yy_t[i] = 0.0;
				for(int j = 0 ; j < m_len_xx ; j ++)
				{
					yy_t[i] += yy[j] * orth[i][j];
				}
			}
			for(int i = 0 ; i < m_len_xx; i ++)
			{
				yy[i] = rr[i] + sigma * yy_t[i];//AWGN
			}

			//for(int i = 0 ; i < m_len_xx; i ++)
			//{
			//	yy[i] = rr[i] + sigma * yy[i];//AWGN
			//}

			total_iter += ij->FFTQSPA(yy ,uu_hat, sigma);//soft decoding
//destination
			for(int i = 0 ; i < total_length ; i ++)
			{
				for(int ii = 0 ; ii < ij->m_degree ; ii ++)
				{
					bit[i*ij->m_degree + ii] = (cc[i]>>ii)&1;
					bit_hat[i*ij->m_degree + ii] = (uu_hat[i]>>ii)&1;
				}
			}
			source_sink_ber->CntErr(bit, bit_hat, inform_len, 1);
			source_sink->CntErr_q(cc, uu_hat, total_length, 1);
//error pattern
			count = 0;
			for(int i = 0 ; i < total_length ; i ++)
			{
				uu_hat[i] = GF.Add(cc[i], uu_hat[i]);
				if(uu_hat[i])
					count += 1;
			}
			if(count)
			{
//output the codeword and the received signal
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<cc[i]<<" ";
				}
				oFile_data<<endl;
				for(int i = 0 ; i < total_length ; i ++)
				{
					oFile_data<<yy[i]<<" ";
				}
				oFile_data<<endl;
//output the error pattern
				oFile_err<<count<<"  "<<ij->count_unsat(uu_hat)<<"     ";
				for(int i = 0 ; i < total_length ; i ++)
				{
					if(uu_hat[i])
					{
						oFile_err<<i<<" "<<uu_hat[i]<<"   ";
					}
				}
				oFile_err<<std::endl;
			}
			if((int)source_sink->_tot_blk_num % 1 == 0)
			{		
				SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
				BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
				FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
				std::cout<<"Frame number: "<<source_sink->_tot_blk_num<<" Error frame number : "<<source_sink->_err_blk_num<<std::endl;
				std::cout<<"SNR          BER            SER          FER"<<std::endl;
				std::cout<<snr<<"          "<<BER<<"          "<<SER<<"           "<<FER<<std::endl;
				std::cout<<"Average iteration numbers is "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0);rndGen1.PrintState(fp_1);fclose(fp_0);fclose(fp_1);
			}

		}
		BER = (source_sink_ber->_err_sym_num) /(source_sink_ber->_tot_sym_num);
		SER = (source_sink->_err_sym_num) /(source_sink->_tot_sym_num);
		FER = (source_sink->_err_blk_num) /(source_sink->_tot_blk_num);
		oFile_ber<<snr<<"                 "<<BER<<std::endl;
        oFile_fer<<snr<<"                 "<<FER<<std::endl;
        oFile_ser<<snr<<"                 "<<SER<<std::endl;
        iter<<snr<<"                 "<<(double)total_iter / (double)(source_sink->_tot_blk_num)<<std::endl;
	}
//clear
	oFile_ber.close();oFile_fer.close();oFile_ser.close();oFile_err.close();iter.close();oFile_data.close();
	delete source_sink;delete []uu;delete []cc;delete []int_cc;delete []rr;delete []yy;delete []yy_t;delete ij;delete []bit;delete []bit_hat;
	for(int i = 0 ; i < inform_len ; i ++)
	{
		delete []orth[i];
	}
	delete orth;
}


void simulation_qary_bibo()
// 函数的功能：对二进制输入/二进制输出（BiBo）的LDPC译码链路进行蒙特卡洛仿真，在不同SNR下通过AWGN信道传输并统计BER，同时记录平均迭代次数及随机数发生器状态。
// 输入数据的内容及格式：从文本参数文件“input_file_q_ary_2013.txt”中读取仿真参数（low_snr、max_snr、step、max_block_num、max_err_num、max_iteration、码/映射相关文件名等）；运行过程中由随机源产生长度为 m_infobit_len 的0/1比特序列作为信息比特。
// 输出数据的内容及格式：将每个SNR点对应的BER追加写入“ber.txt”，将每个SNR点对应的平均迭代次数追加写入“iter.txt”；并在控制台输出当前帧计数、误帧计数、SNR、BER与平均迭代次数；同时输出随机数发生器状态到“rnd_state_0.txt”和“rnd_state_1.txt”。
{
	int flag = 0, count, value = 0;
	double low_snr = 0.0, max_snr = 15.0, step = 0.2, var, sigma, temp0, temp1,tempsum;

	FILE* fp_0, * fp_1;
	int max_block_num = 1000000, max_err_num = 1000, q_ary = 16, max_iteration = 20, m_len_input = 5;
	string filename = "", filename_tanner = "", filename_mapping = "", file_num = "1", file_name = "input_file_q_ary_", filename_orth = "";

	flag = -1;
	//set_flag(flag , 1);
	rndGen0.SetSeed(flag);
	//set_flag(flag , 2);
	rndGen1.SetSeed(flag);

	std::cout << "Please input the number of the input file  :  " << std::endl << ">:";
	file_num = "2013";
	//std::cin>>file_num;
	file_name = file_name + file_num + ".txt";
	std::ifstream inFile;
	inFile.open(file_name.c_str());
	if (!inFile)
	{
		std::cerr << "Can't open the parameter input file : " << file_name << std::endl;
	}
	//Malloc
	char* temp_string = new char[1024], * temp_value = new char[1024];
	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	low_snr = atof(temp_value);

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	max_snr = atof(temp_value);

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	step = atof(temp_value);

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	max_block_num = atoi(temp_value);;

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	max_err_num = atoi(temp_value);

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_value, 1024);
	max_iteration = atoi(temp_value);

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_string, 1024);
	filename = filename + temp_string;

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_string, 1024);
	filename_mapping = filename_mapping + temp_string;

	inFile.getline(temp_string, 1024);
	inFile.getline(temp_string, 1024);
	filename_orth = filename_orth + temp_string;

	inFile.close();	delete[]temp_string; delete[]temp_value;
	//Source and Sink
	CSourceSink* source_sink = new CSourceSink();
	CSourceSink* source_sink_ber = new CSourceSink();

	BCJRQSPA* ij = new BCJRQSPA(filename, max_iteration, filename_mapping);

	int length = ij->get_block_length(), total_length = ij->get_total_length(), total_iter = 0;

	int* uu = new int[total_length];//orginal symbol string
	int* uu_hat = new int[total_length];//symbol string after decoding
	int* cc = new int[total_length];//codeword
	int* int_cc = new int[total_length];//interleaved codeword

	for (int i = 0; i < total_length; i++)
	{
		cc[i] = 0; int_cc[i] = 0;
	}
	int m_infobit_len = length * ij->m_degree;
	int m_len_xx = ij->m_len_xx, m_codebit_len = total_length * ij->m_degree;
	int* bit = new int[m_codebit_len];
	int* bit_cc = new int[m_codebit_len];
	int* bit_hat = new int[m_codebit_len];
	double* rr = new double[m_len_xx];//output of Mapper
	double* yy = new double[m_len_xx];//output of AWGN
	double* yy_t = new double[m_len_xx];//output of AWGN
	q_ary = ij->q_ary;

	double SER = 0.0, BER = 0.0, FER = 0.0;
	std::ofstream oFile_ber("ber.txt", ios_base::app);
	std::ofstream iter("iter.txt", ios_base::app);
	if ((!oFile_ber))
	{
		std::cerr << "Can't open the file ber.txt or the file fer.txt !Please check it !" << std::endl;
		exit(-1);
	}
	oFile_ber << "The following result corresponding the input file : " << file_name << std::endl;
	oFile_ber << "  SNR                         BER " << std::endl;

	std::cout << endl << "----------------------SIMULATION BEGIN----------------------" << endl;

	double rate = (double)length / (double)total_length;
	for (double snr = low_snr; snr < max_snr; snr = snr + step)
	{
		var = pow(10.0, 0.1 * snr);
		sigma = sqrt(1 / (2.0 * var * rate));//Cal channel parameters
		source_sink->ClrCnt(); source_sink_ber->ClrCnt(); BER = 0.0; FER = 0.0; total_iter = 0;//clear source and sink
		while (source_sink->_tot_blk_num < max_block_num && source_sink->_err_blk_num < max_err_num)
		{
			source_sink->GetBitStr(bit, m_infobit_len);//source
			ij->encoder4BiBo(bit, bit_cc);//encoding
			for (int i = 0; i < m_codebit_len; i++)
			{
				if (bit_cc[i] == 0)
				{
					rr[i] = -1.0;
				}
				else
				{
					rr[i] = 1.0;
				}
			}
			rndGen1.Normal(yy, m_codebit_len);//Gaussian noise

			for (int i = 0; i < m_codebit_len; i++)
			{
				yy[i] = rr[i] + sigma * yy[i];//AWGN
				temp0 = exp(-0.5 * (yy[i] + 1) * (yy[i] + 1) / sigma);
				temp1 = exp(-0.5 * (yy[i] - 1) * (yy[i] - 1)/ sigma);
				temp0 = temp0 / (temp0 + temp1);
				yy[i] = temp0;
			}
			total_iter += ij->FFTQSPA4BiBo(yy, bit_hat);//soft decoding

			source_sink_ber->CntErr(bit_cc, bit_hat, m_codebit_len, 1);
			if ((int)source_sink_ber->_tot_blk_num % 1 == 0)
			{
				BER = (source_sink_ber->_err_sym_num) / (source_sink_ber->_tot_sym_num);
				std::cout << "Frame number: " << source_sink_ber->_tot_blk_num << " Error frame number : " << source_sink_ber->_err_blk_num << std::endl;
				std::cout << "SNR          BER              " << std::endl;
				std::cout << snr << "          " << BER << std::endl;
				std::cout << "Average iteration numbers is " << (double)total_iter / (double)(source_sink_ber->_tot_blk_num) << std::endl;
				if ((fp_0 = fopen("rnd_state_0.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				if ((fp_1 = fopen("rnd_state_1.txt", "w+")) == NULL)
				{
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				rndGen0.PrintState(fp_0); rndGen1.PrintState(fp_1); fclose(fp_0); fclose(fp_1);
			}
		}
		BER = (source_sink_ber->_err_sym_num) / (source_sink_ber->_tot_sym_num);
		oFile_ber << snr << "                 " << BER << std::endl;
		iter << snr << "                 " << (double)total_iter / (double)(source_sink->_tot_blk_num) << std::endl;
	}
	//clear
	oFile_ber.close(); iter.close();
	delete source_sink; delete[]uu; delete[]cc; delete[]int_cc; delete[]rr; delete[]yy; delete[]yy_t; delete ij; delete[]bit; delete[]bit_hat;
	delete[]bit_cc;
}