#include "stdafx.h"
#include "BCJRQSPA.h"

#define MAX_D 100000.0
#define SMALLPROB 1.0e-12
#define SUPSMALLPROB 1.0e-30

int MMAAXX = 20;
extern CLCRandNum rndGen0;
CFiniteField2 GF;
const int len_pad_sym = 5;//number of [0] padded to terminate the trellis
using namespace std;

BCJRQSPA::BCJRQSPA(std::string filename, int max_iteration, string mapping_name)
{
	this->max_iteration = max_iteration;
	this->Malloc(filename, mapping_name);

}
void BCJRQSPA::Malloc(string filename, string mapping_name)
{
//开始读取参数
	FILE *fp;
	char temp_string[80];
	int temp_col, temp_row, temp_weight, temp_gf_number;

	if( (fopen_s(&fp,filename.c_str(),"r")) !=0 )
	{   
		fprintf(stderr,"\n cannot open %s",filename.c_str());
		exit(1);
	}
    fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->parity_row);
	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->parity_column);
	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->rank);
	fscanf_s(fp,"%s",temp_string,80);
	fscanf_s(fp,"%d",&this->q_ary);
	this->row_weight_enc = new int[this->rank];
	this->row_weight_parity = new int[this->parity_row];
	m_degree = 1;
	while((1<<m_degree) < this->q_ary)
		m_degree ++;
	GF.Malloc(m_degree);
	this->parity_matrix = new int*[this->parity_row];
	this->row_location_parity = new int*[this->parity_row];

	fscanf_s(fp,"%s",temp_string,80);  //*****1's_location/per_row****
	for (int i = 0;i < this->parity_row ; i ++)
	{   
		fscanf_s(fp,"%d",&temp_row);
		fscanf_s(fp,"%d",&temp_weight);
		this->parity_matrix[i] = new int[temp_weight];
		this->row_location_parity[i] = new int[temp_weight];
		this->row_weight_parity[i] = temp_weight;
		for (int j = 0 ; j < temp_weight ; j ++)
		{    
			fscanf_s(fp,"%d",&temp_col);
			fscanf_s(fp,"%d",&temp_gf_number);
			this->row_location_parity[temp_row][j] = temp_col;
			this->parity_matrix[temp_row][j] = int(temp_gf_number);
		}
	}
	
	m_encH = new int*[this->parity_row];
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		this->m_encH[i] = new int[this->parity_column];
	}
	m_decH = new int*[this->parity_row];
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		this->m_decH[i] = new int[this->parity_column];
	}
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			this->m_encH[i][j] = 0;
		}
	}
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		for(int j = 0 ; j < this->row_weight_parity[i] ; j ++)
		{
			this->m_encH[i][this->row_location_parity[i][j]] = this->parity_matrix[i][j];
		}
	}
	tempP = new int[this->parity_column];
	for (int j = 0; j < this->parity_column; j++)
	{
		tempP[j] = j;
	}
	sym = new int[parity_column];
//读入系统化后的编码矩阵
	this->enc_parity_matrix = new int*[this->rank];
	this->row_location_enc = new int*[this->rank];
	fscanf_s(fp,"%s",temp_string,80);  //*****1's_location/per_row****
	for (int i = 0;i < this->rank ; i ++)
	{   
		fscanf_s(fp,"%d",&temp_row);
		fscanf_s(fp,"%d",&temp_weight);
		this->row_weight_enc[i] = temp_weight;
		this->enc_parity_matrix[i] = new int[temp_weight];
		this->row_location_enc[i] = new int[temp_weight];
		for (int j = 0 ; j < temp_weight ; j ++)
		{    
			fscanf_s(fp,"%d",&temp_col);
			fscanf_s(fp,"%d",&temp_gf_number);
			this->row_location_enc[temp_row][j] = temp_col;
			this->enc_parity_matrix[temp_row][j] = int(temp_gf_number);
		}
	}
	fclose(fp);
	tanner.Construct_Tanner_Graph(filename);//构造Tanner图
	m_modem.Malloc(mapping_name.c_str());
	if(m_modem.m_num_signal != q_ary)
	{
		fprintf(stderr,"\nCardinality of signal set and order of finite field don't match.\n");
		system("pause");
	}
	m_len_xx = parity_column * m_modem.m_len_signal;

	channel_for_spa = new double*[this->get_total_length()];
	for(int i = 0 ; i < this->get_total_length() ; i ++){ 
		channel_for_spa[i] = new double[this->q_ary];}
	total_mul = new double[this->q_ary];
	pai = new int[parity_column];
	for(int i = 0 ; i < parity_column ; i ++){ 
		pai[i] = i;}
	temp_qary = new double[q_ary];
	inds = new int[parity_column];
	sym_inds = new int[parity_row];
	m_bit_len_cc = parity_column * m_degree;
	m_bit_len_info = (parity_column - rank) * m_degree;
}


void BCJRQSPA::InterleavingMessage(double **in, double **out, int col_len, int row_len)
{
	int i,j;
	for(i = 0 ; i < col_len ; i ++)
	{
		for(j = 0 ; j < row_len ; j ++)
		{
			out[i][j] = in[pai[i]][j];
		}
	}
}
/*
pai    = [2] [3] [0] [1] [4]
input  = [a] [b] [c] [d] [e]
output[i] = input[pai[i]]
output = [c] [d] [a] [b] [e]
*/
void BCJRQSPA::InterleavingSymbol(int *in_cc , int *out_cc , int input_len)
{
	int i;
	for(i = 0 ; i < input_len ; i ++)
	{
		out_cc[i] = in_cc[pai[i]];
	}
}

void BCJRQSPA::DeInterleavingMessage(double **in, double **out, int col_len, int row_len)
{
	int i,j;
	for(i = 0 ; i < col_len ; i ++)
	{
		for(j = 0 ; j < row_len ; j ++)
		{
			out[pai[i]][j] = in[i][j];
		}
	}
}


void BCJRQSPA::Init_Interleaver(int is_int)
{
	int t;
	double random_num;
	int temp;
	int position;

	for (t = 0; t < parity_column; t++)
	{
		pai[t] = t;
	}

	if(is_int)
	{
		for (t = parity_column - 1; t > 0; t--)
		{
			random_num = rndGen0.Uniform();
			position = (int) (random_num * t);
			temp = pai[position];
			pai[position] = pai[t];
			pai[t] = temp;
		}
	}
	return;
}


void BCJRQSPA::Sys()
{
	int flag;
	int i;
	int ii;
	int jj;
	int j;
	int n;
	int m;
	this->rank = 0;
	int temp = 0;
	int temp_index;

	int **tempH = new int*[this->parity_row];
	for(i = 0 ; i < this->parity_row ; i ++)
	{
		tempH[i] = new int[this->parity_column];
	}

	for(i = 0 ; i < this->parity_row ; i ++)
	{
		for(j = 0 ; j < this->parity_column ; j ++)
		{
			tempH[i][j] = 0;
		}
	}

	for(i = 0 ; i < this->parity_row ; i ++)
	{
		for(j = 0 ; j < this->row_weight_parity[i] ; j ++)
		{
			tempH[i][this->row_location_parity[i][j]] = this->parity_matrix[i][j];
		}
	}
	for(i = 0 ; i < this->parity_column ; i ++)
	{
		tempP[i] = i;
	}

	for (i = 0; i < this->parity_row; i ++)
	{
		//从第i行第i个位置开始，找第一个非0元素，从上往下，从左往往右。先找行，再找列。
		flag = 0;
		for (jj = i; jj < this->parity_column; jj ++)
		{
			for (ii = i; ii < this->parity_row; ii ++)
			{
				if (this->m_encH[ii][jj] != 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 1)
			{
				this->rank++;
				break;
			}
		}

		if (flag == 0)
			break;
		else
		{
//swap i and ii row
//换行
			if (ii != i)
			{
				for (n = 0; n < parity_column; n++)
				{
					temp = m_encH[i][n];
					m_encH[i][n] = m_encH[ii][n];
					m_encH[ii][n] = temp;
				}
			}
//swap i and jj col
//换列
			if (jj != i)
			{
				temp_index = tempP[i];
				tempP[i] = tempP[jj];
				tempP[jj] = temp_index;

				for (m = 0; m < parity_row; m++)
				{
					temp = m_encH[m][i];
					m_encH[m][i] = m_encH[m][jj];
					m_encH[m][jj] = temp;
				}
			}
//elimination
//消元
			for (m = 0; m < parity_row; m++)
			{
				if (m != i && m_encH[m][i] != 0)
				{
					temp = GF.Div(1, m_encH[i][i]);
					temp = GF.Mult(temp, m_encH[m][i]);
					for (n = 0; n < parity_column; n++)
					{
						m_encH[m][n] = GF.Add(m_encH[m][n], GF.Mult(temp, this->m_encH[i][n]));
					}
				}
			}
			temp = GF.Div(1, m_encH[i][i]);
			for (n = 0; n < parity_column; n++)
			{
				m_encH[i][n] = GF.Mult(temp, this->m_encH[i][n]);
			}
		}
	}
	for(i = 0 ; i <  this->parity_row ; i ++)
	{
		for(j=0 ; j < this->parity_column ; j ++)
		{
			m_decH[i][j] = tempH[i][tempP[j]];
		}
	}
	n = 0;

	std::ofstream mat("matrix.txt" , ios_base::app);
	mat<<"row_number"<<std::endl;
	mat<<parity_row<<std::endl;
	mat<<"col_number"<<std::endl;
	mat<<parity_column<<std::endl;
	mat<<"rank"<<std::endl;
	mat<<rank<<std::endl;
	mat<<"ary"<<std::endl;
	mat<<q_ary<<std::endl;
	mat<<"1's_location_per_row"<<endl;

	for(i = 0 ; i < parity_row ; i ++)
	{
		n = 0;
		for(j = 0 ; j < parity_column ; j ++)
		{
			if(m_decH[i][j])
				n ++;
		}
		mat<<i<<"   "<<n<<"   ";
		for(j = 0 ; j < parity_column ; j ++)
		{
			if(m_decH[i][j])
				mat<<j<<"   "<<m_decH[i][j]<<"   ";
		}
		mat<<endl;
	}

	mat<<"1's_location_per_row_of_encoding_matrix"<<endl;

	for(i = 0 ; i < rank ; i ++)
	{
		n = 0;
		for(j = 0 ; j < parity_column ; j ++)
		{
			if(m_encH[i][j])
				n ++;
		}
		mat<<i<<"   "<<n<<"   ";
		for(j = 0 ; j < parity_column ; j ++)
		{
			if(m_encH[i][j])
				mat<<j<<"   "<<m_encH[i][j]<<"   ";
		}
		mat<<endl;
	}

	mat.close();

	for(i = 0 ; i < this->parity_row ; i ++)
	{
		delete []tempH[i];
	}
	delete []tempH;
}


BCJRQSPA::~BCJRQSPA()
{
	int des;
	// std::cout<<"in destructor:";
	// std::cin>>des;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		delete []this->parity_matrix[i];
	}
	delete []this->parity_matrix;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		delete []this->m_encH[i];
	}
	delete []this->m_encH;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		delete []this->m_decH[i];
	}
	delete []this->m_decH;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		delete []this->row_location_parity[i];
	}
	delete []this->row_location_parity;	
	for(int i = 0 ; i < this->rank ; i ++)
	{
		delete []this->enc_parity_matrix[i];
	}
	delete []this->enc_parity_matrix;
	for(int i = 0 ; i < this->rank ; i ++)
	{
		delete []this->row_location_enc[i];
	}
	delete []this->row_location_enc;
	delete []this->row_weight_enc;
	delete []this->row_weight_parity;
	
	for(int i = 0 ; i < this->parity_column ; i ++){ 
		delete []channel_for_spa[i];}
	delete []channel_for_spa;
	
	delete []total_mul;

	delete []tempP;

	delete []sym;

	delete []temp_qary;

	delete []inds;

	delete []sym_inds;
}

int BCJRQSPA::get_block_length() const
{
	return this->parity_column - this->rank;
}
int BCJRQSPA::get_total_length() const
{
	return this->parity_column;
}

void BCJRQSPA::no_encoder(int *m_uu, int *b_cc)
{
	int parity_value = 0;
	for(int i = 0 ; i < this->parity_column - this->rank ; i ++)
	{
		b_cc[i + this->rank] = 0;
		m_uu[i] = 0;
	}
	for(int i = 0 ; i < this->rank ; i ++)
	{
		parity_value = 0;
		for(int j = 1 ; j < this->row_weight_enc[i] ; j ++)
		{
			parity_value = GF.Add(parity_value, GF.Mult(b_cc[this->row_location_enc[i][j]], this->enc_parity_matrix[i][j]));
		}
		b_cc[i] = parity_value;
	}
}




void BCJRQSPA::re_encoder(int *m_uu, int *b_cc)
{
	int parity_value = 0;
	for(int i = 0 ; i < this->parity_column - this->rank ; i ++)
	{
		b_cc[i + this->rank] = m_uu[i];
	}
	for(int i = 0 ; i < this->rank ; i ++)
	{
		parity_value = 0;
		for(int j = 1 ; j < this->row_weight_enc[i] ; j ++)
		{
			parity_value = GF.Add(parity_value, GF.Mult(b_cc[this->row_location_enc[i][j]], this->enc_parity_matrix[i][j]));
		}
		b_cc[i] = parity_value;
	}
	//if(!this->is_codeword(b_cc))
	//{
	//	std::cout<<"Encoding failure!!!!!input any number!!!!";
	//}
	//std::cout<<this->is_codeword(b_cc)<<"   ";
}


//二进制输入，二进制输出的多元LDPC编码器。
//*******************************多元码编码器*******************************
void BCJRQSPA::encoder4BiBo(int* m_uu, int* b_cc)
/**
 * 函数的功能：
 *   将输入的信息比特序列按 m_degree 位打包成 q-ary 符号（q_ary = 2^m_degree），
 *   依据系统化编码矩阵计算校验符号生成码字，并将最终的 q-ary 码字再展开为二进制比特序列输出。
 *
 * 输入数据的内容及格式：
 *   m_uu：
 *     指向 int 数组的指针，输入为二进制序列（元素取 0/1），
 *     长度应为 (parity_column - rank) * m_degree，对应信息部分比特。
 *   b_cc：
 *     指向 int 数组的指针，用于输出码字比特序列，
 *     需至少可容纳 parity_column * m_degree 个元素（函数内部先写入 q-ary 符号，再展开为比特覆盖写回）。
 *
 * 输出数据的内容及格式：
 *   b_cc（输出参数）：
 *     输出为二进制码字序列，长度为 parity_column * m_degree；
 *     第 j 个符号对应的比特为 b_cc[j*m_degree + ii] = (symbol >> ii) & 1。
 */
{
	int temp_Q, i;
	int parity_value = 0;
	//将m_uu由 二进制序列 转换为 多进制序列
	for (int i = 0; i < this->parity_column - this->rank; i++)
	{
		temp_Q = 0;
		for (int ii = 0; ii < m_degree; ii++)
		{
			temp_Q += m_uu[i * m_degree + ii] << ii;
		}
		b_cc[i + this->rank] = temp_Q;
	}
	for (int i = 0; i < this->rank; i++)
	{
		parity_value = 0;
		for (int j = 1; j < this->row_weight_enc[i]; j++)
		{
			parity_value = GF.Add(parity_value, GF.Mult(b_cc[this->row_location_enc[i][j]], this->enc_parity_matrix[i][j]));
		}
		b_cc[i] = parity_value;
	}
	//将b_cc由 多进制序列 转换为 二进制序列
	for (int j = parity_column - 1; j >= 0; j--)
	{
		b_cc[j * m_degree] = b_cc[j];
	}
	for (int j = 0; j < parity_column; j++)
	{
		i = b_cc[j * m_degree];
		for (int ii = 0; ii < m_degree; ii++)
		{
			b_cc[j * m_degree + ii] = (i >> ii) & 1;
		}
	}
}
//**************************************************************************

void BCJRQSPA::sym_encoder(int *uu)
{
	ldpc_edge *edge;
	int parity_value = 0;
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		edge = this->tanner.m_row_head[i].left;
		sym[i] = 0;
		while(edge->col_num != -1)
		{
			sym[i] = GF.Add(sym[i], GF.Mult(edge->pai_element, uu[edge->col_num]));
			edge = edge->left;
		}
	}
}


bool BCJRQSPA::is_codeword(int *b_code)
{
	bool code_success = true;
	int parity_check = 0;
	sym_unsat_num = 0;

	ldpc_edge *edge;

	for(int i = 0 ; i < parity_column ; i ++)
	{
		inds[i] = 0;
	}
	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		sym_inds[i] = 0;
		edge = this->tanner.m_row_head[i].left;
		parity_check = 0;
		while(edge->col_num != -1)
		{
			parity_check = GF.Add(parity_check, GF.Mult(edge->pai_element, b_code[edge->col_num]));
			edge = edge->left;
		}
		if (parity_check != 0)
		{
			code_success = false;
			
			sym_unsat_num ++;

			sym_inds[i] = 1;
			
			edge = this->tanner.m_row_head[i].left;
			while(edge->col_num != -1)
			{
				inds[edge->col_num] ++;
				edge = edge->left;
			}
			//break;
		}
	}
	return code_success;
}


bool BCJRQSPA::is_codeword_sym(int *b_code)
{
	bool code_success = true;
	int parity_check = 0;

	ldpc_edge *edge;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		edge = this->tanner.m_row_head[i].left;
		parity_check = 0;
		while(edge->col_num != -1)
		{
			parity_check = GF.Add(parity_check, GF.Mult(edge->pai_element, b_code[edge->col_num]));
			edge = edge->left;
		}
		if (parity_check != sym[i])
		{
			code_success = false;
			break;
		}
	}
	return code_success;
}

int BCJRQSPA::FFTQSPA_SYM(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	ldpc_edge *temp_edge;

	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while(temp_edge->row_num != -1)
		{
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter = 0;
	for(iter = 0 ; iter < this->max_iteration ; iter ++)
	{
//begin update the information from variable node to check node
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}
				
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					curr_sum += total_mul[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] /= curr_sum;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					if(total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if(total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				if(curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while(temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if(curr_sum == 0.0)
				{
					std::cout<<"DIVIDE BY 0"<<std::endl;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if(temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if(temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
		}
//end update information from variable node to check node

		if(this->is_codeword_sym(uu))	
		{
			break;
		}
//begin information permutation
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
//end information permutation

//begin information update at the check node
		for(int j = 0 ; j < this->parity_row ; j ++)
		{
			//init
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = 1.0;
			}
			//cal the total multi
			temp_edge = this->tanner.m_row_head[j].left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *=temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//cal the output from check node to variable nodes
			temp_edge = temp_edge->left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_c2v[ii] = total_mul[ii] / temp_edge->m_v2c[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_qary[ii] = temp_edge->m_c2v[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_c2v[ii] = temp_qary[GF.Add(ii, sym[temp_edge->row_num])];
				}
				temp_edge = temp_edge->left;
			}
		}
//end information update at the check node
	}
	return iter;
}

int BCJRQSPA::count_unsat(int *b_code)
{
	int count = 0;
	int parity_check = 0;

	ldpc_edge *edge;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		edge = this->tanner.m_row_head[i].left;
		parity_check = 0;
		while(edge->col_num != -1)
		{
			parity_check = GF.Add(parity_check, GF.Mult(edge->pai_element, b_code[edge->col_num]));
			edge = edge->left;
		}
		if (parity_check != 0)
		{
			count ++;
		}
	}
	return count;
}



int BCJRQSPA::count_unsat_sym(int *b_code)
{
	int count = 0;
	int parity_check = 0;

	ldpc_edge *edge;

	for(int i = 0 ; i < this->parity_row ; i ++)
	{
		edge = this->tanner.m_row_head[i].left;
		parity_check = 0;
		while(edge->col_num != -1)
		{
			parity_check = GF.Add(parity_check, GF.Mult(edge->pai_element, b_code[edge->col_num]));
			edge = edge->left;
		}
		if (parity_check != 0)
		{
			count ++;
		}
	}
	return count;
}


int BCJRQSPA::FFTQSPABIA(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	int dex;
	int acc = 0;
	ldpc_edge *temp_edge;

	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while(temp_edge->row_num != -1)
		{
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter = 0;
	for(iter = 0 ; iter < this->max_iteration ; iter ++)
	{
//begin update the information from variable node to check node
		if(acc > 0)
		{
			acc --;
		}
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}
				
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					curr_sum += total_mul[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] /= curr_sum;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					if(total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if(total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				if(curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while(temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if(curr_sum == 0.0)
				{
					std::cout<<"DIVIDE BY 0"<<std::endl;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if(temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if(temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
		}
//end update information from variable node to check node

		if(this->is_codeword(uu))
		{
			break;
		}
		if((sym_unsat_num <= 6) && (iter > max_iteration - 60) && (acc <= 0))
		{
			acc = 6;
			for(int j = 0 ; j < this->parity_column ; j ++)
			{
				if(inds[j])
				{
					temp_edge = this->tanner.m_col_head[j].up;
					while(temp_edge->row_num != -1)
					{
						if(sym_inds[temp_edge->row_num])
						{
							curr_sum = -1.0;
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{
								if(curr_sum < temp_edge->m_v2c[ii])
								{
									dex = ii;
									curr_sum = temp_edge->m_v2c[ii];
								}
							}
							temp_edge->m_v2c[dex] *= 1.4;

							curr_sum = 0.0;
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{					
								curr_sum += temp_edge->m_v2c[ii];
							}
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{
								temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
							}
						}
						temp_edge = temp_edge->up;
				     }
				}
				else
				{
					temp_edge = this->tanner.m_col_head[j].up;
					while(temp_edge->row_num != -1)
					{
						if(sym_inds[temp_edge->row_num])
						{
							curr_sum = -1.0;
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{
								if(curr_sum < temp_edge->m_v2c[ii])
								{
									dex = ii;
									curr_sum = temp_edge->m_v2c[ii];
								}
							}
							temp_edge->m_v2c[dex] *= 0.8;

							curr_sum = 0.0;
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{					
								curr_sum += temp_edge->m_v2c[ii];
							}
							for(int ii = 0 ; ii < this->q_ary ; ii ++)
							{
								temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
							}
						}
						temp_edge = temp_edge->up;
				     }
				}
			}
		}
//begin information permutation
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
//end information permutation

//begin information update at the check node
		for(int j = 0 ; j < this->parity_row ; j ++)
		{
			//init
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = 1.0;
			}
			//cal the total multi
			temp_edge = this->tanner.m_row_head[j].left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *=temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//cal the output from check node to variable nodes
			temp_edge = temp_edge->left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_c2v[ii] = total_mul[ii] / temp_edge->m_v2c[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);
				temp_edge = temp_edge->left;
			}
		}
//end information update at the check node
	}
	return iter;
}


int BCJRQSPA::FFTQSPATrunc(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	ldpc_edge *temp_edge;

	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while(temp_edge->row_num != -1)
		{
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter;
	for(iter = 0 ; iter < this->max_iteration ; iter ++)
	{
//begin update the information from variable node to check node
		if(iter == max_iteration - 1)
			iter = iter;
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}
				
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					curr_sum += total_mul[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] /= curr_sum;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					if(total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if(total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				if(curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while(temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if(curr_sum == 0.0)
				{
					std::cout<<"DIVIDE BY 0"<<std::endl;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if(temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if(temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
/////////////////////////////////////////////////////////////
				FFTTruncLEFT(temp_edge->m_v2c, q_ary, MMAAXX);
/////////////////////////////////////////////////////////////
				temp_edge = temp_edge->up;
			}
		}
//end update information from variable node to check node

		if(this->is_codeword(uu))	
		{
			break;
		}
//begin information permutation
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
//end information permutation

//begin information update at the check node
		for(int j = 0 ; j < this->parity_row ; j ++)
		{
			///xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//////////
			//forward 
			temp_edge = this->tanner.m_row_head[j].left;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				temp_edge->m_alpha[ii] = 1.0;
			}
			while (temp_edge->col_num != -1) 
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->left->m_alpha[ii] = temp_edge->m_alpha[ii] * temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//backward
			temp_edge = this->tanner.m_row_head[j].right;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				temp_edge->m_beta[ii] = 1.0;
			}
			while (temp_edge->col_num != -1) 
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->right->m_beta[ii] = temp_edge->m_beta[ii] * temp_edge->m_v2c[ii];
					temp_edge->m_c2v[ii] = temp_edge->m_alpha[ii] * temp_edge->m_beta[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);
/////////////////////////////////////////////////////////////
				//std::cout<<std::endl;
				//std::cout<<std::endl;
				//double sssu = 0.0;
				//for(int ii = 0 ; ii < q_ary ; ii ++)
				//{
				//	std::cout<<temp_edge->m_c2v[ii]<<"  ";
				//	sssu += temp_edge->m_c2v[ii];
				//}
				//std::cout<<std::endl;
				FFTTruncRIGHT(temp_edge->m_c2v, q_ary, MMAAXX);
				//for(int ii = 0 ; ii < q_ary ; ii ++)
				//{
				//	std::cout<<temp_edge->m_c2v[ii]<<"  ";
				//}
				//std::cout<<std::endl;
/////////////////////////////////////////////////////////////
				temp_edge = temp_edge->right;
			}
		}
//end information update at the check node
	}
	return iter;
}


int BCJRQSPA::HardDecision(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	int num;
	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		curr_sum = -1;
		for (int j = 0; j < this->q_ary; j++)
		{
			if (curr_sum<channel_for_spa[i][j])
			{
				curr_sum = channel_for_spa[i][j];
				uu[i] = j;
			}
		}
		curr_sum = rndGen0.Uniform() *q_ary;
	    num = (int)curr_sum;
	    num = num % q_ary;
	    uu[i] = num;
	}
	return this->count_unsat(uu);
}


int BCJRQSPA::FFTQSPA_FB(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	ldpc_edge *temp_edge;

	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while(temp_edge->row_num != -1)
		{
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter = 0;
	for(iter = 0 ; iter < this->max_iteration ; iter ++)
	{
//begin update the information from variable node to check node
		if(iter == max_iteration - 1)
			iter = iter;
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}
				
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					curr_sum += total_mul[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] /= curr_sum;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					if(total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if(total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				if(curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while(temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if(curr_sum == 0.0)
				{
					std::cout<<"DIVIDE BY 0"<<std::endl;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if(temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if(temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
		}
//end update information from variable node to check node

		if(this->is_codeword(uu))	
		{
			break;
		}
//begin information permutation
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
//end information permutation

//begin information update at the check node
		for(int j = 0 ; j < this->parity_row ; j ++)
		{
			///xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//////////
			//forward 
			temp_edge = this->tanner.m_row_head[j].left;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				temp_edge->m_alpha[ii] = 1.0;
			}
			while (temp_edge->col_num != -1) 
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->left->m_alpha[ii] = temp_edge->m_alpha[ii] * temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//backward
			temp_edge = this->tanner.m_row_head[j].right;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				temp_edge->m_beta[ii] = 1.0;
			}
			while (temp_edge->col_num != -1) 
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->right->m_beta[ii] = temp_edge->m_beta[ii] * temp_edge->m_v2c[ii];
					temp_edge->m_c2v[ii] = temp_edge->m_alpha[ii] * temp_edge->m_beta[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);

				temp_edge = temp_edge->right;
			}
			///xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//////////
			////init
			//for(int ii = 0 ; ii < this->q_ary ; ii ++)
			//{
			//	total_mul[ii] = 1.0;
			//}
			////cal the total multi
			//temp_edge = this->tanner.m_row_head[j].left;
			//while(temp_edge->col_num != -1)
			//{
			//	for(int ii = 0 ; ii < this->q_ary ; ii ++)
			//	{
			//		total_mul[ii] *=temp_edge->m_v2c[ii];
			//	}
			//	temp_edge = temp_edge->left;
			//}

			////cal the output from check node to variable nodes
			//temp_edge = temp_edge->left;
			//while(temp_edge->col_num != -1)
			//{
			//	for(int ii = 0 ; ii < this->q_ary ; ii ++)
			//	{
			//		temp_edge->m_c2v[ii] = total_mul[ii] / temp_edge->m_v2c[ii];
			//	}
			//	temp_edge->Transform(1);
			//	temp_edge->Permutation(1);
			//	temp_edge->Norm(1);
			//	temp_edge = temp_edge->left;
			//}
		}
//end information update at the check node
	}
	return iter;
}


//fft-qspa译码器：输入rr_bits_prob为比特概率
//              ：输出uu_bit为译出的比特序列
//*******************************多元码译码器*******************************

int BCJRQSPA::FFTQSPA4BiBo(double* rr_bits_prob, int* uu_bit)
/**
 * 函数的功能：
 *   基于 FFT-QSPA（频域/Hadamard 变换的 q-ary SPA）对由二进制比特组成的接收概率信息进行迭代译码。
 *   先将每组 m_degree 个二进制比特的“取0概率”组合为一个 q-ary 符号的信道先验概率（q_ary = 2^m_degree），
 *   然后在 Tanner 图上进行变量结点/校验结点消息传递更新（含归一化、截断、置换、变换等），
 *   迭代直至满足校验（is_codeword）或达到最大迭代次数。最后将 q-ary 判决结果展开为二进制比特输出。
 *
 * 输入数据的内容及格式：
 *   rr_bits_prob：
 *     指向 double 数组的指针，长度应为 parity_column * m_degree。
 *     rr_bits_prob[i*m_degree + j] 表示第 i 组（第 i 个 q-ary 符号对应的）第 j 个二进制比特取值为 0 的概率 P(b=0)。
 *   uu_bit：
 *     指向 int 数组的指针，要求至少可容纳 parity_column * m_degree 个元素（函数会写入）。
 *     注意：在迭代过程中内部也会暂存每个符号的 q-ary 判决（uu_bit[j] = 符号值），结束前会展开覆盖为比特。
 *
 * 输出数据的内容及格式：
 *   返回值：
 *     实际执行的迭代次数 iter（若提前满足校验则为提前终止时的迭代计数，否则为 max_iteration）。
 *   uu_bit（输出参数）：
 *     写回译码后的二进制比特序列，长度为 parity_column * m_degree；
 *     对于每个符号 j，uu_bit[j*m_degree + ii] 为该符号展开后的第 ii 位（(symbol >> ii) & 1）。
 */
{
	double curr_sum = 0.0;
	int i, q, j, bit_temp;
	ldpc_edge* temp_edge;

//由二进制比特的概率向量rr_bits_prob转换为q-ary符号的概率向量
//其中rr_bits_prob为取0的概率
	for (i = 0; i < parity_column; i++)
	{
		for (q = 0; q < q_ary; q++)
		{
			curr_sum = 1.0;
			for (j = 0; j < m_degree; j++)
			{
				bit_temp = (q >> j) & 1;
				if (bit_temp == 0)
				{
					curr_sum *= rr_bits_prob[i * m_degree + j];
				}
				else
				{
					curr_sum *= (1 - rr_bits_prob[i * m_degree + j]);
				}
			}
			channel_for_spa[i][q] = curr_sum;
		}
		ProbClip(channel_for_spa[i], q_ary);
		curr_sum = 0.0;
		for (q = 0; q < q_ary; q++)
		{
			curr_sum += channel_for_spa[i][q];
		}
		for (q = 0; q < q_ary; q++)
		{
			channel_for_spa[i][q] = channel_for_spa[i][q] / curr_sum;
		}
	} //end i


	//init metric from variable node to check node
		 //metric from check node to variable node
	for (int i = 0; i < this->parity_column; i++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while (temp_edge->row_num != -1)
		{
			for (int j = 0; j < this->tanner.q_ary; j++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for (int j = 0; j < this->tanner.q_ary; j++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter = 0;
	for (iter = 0; iter < this->max_iteration; iter++)
	{
		//begin update the information from variable node to check node
		for (j = 0; j < this->parity_column; j++)
		{
			for (int ii = 0; ii < this->q_ary; ii++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while (temp_edge->row_num != -1)
			{
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}

				curr_sum = 0.0;
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					curr_sum += total_mul[ii];
				}
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					total_mul[ii] /= curr_sum;
				}
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					if (total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if (total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for (int ii = 0; ii < this->q_ary; ii++)
			{
				if (curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu_bit[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while (temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if (curr_sum == 0.0)
				{
					std::cout << "DIVIDE BY 0" << std::endl;
				}
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if (temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if (temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
		}
		//end update information from variable node to check node

		if (this->is_codeword(uu_bit))
		{
			break;
		}
		//begin information permutation
		for (j = 0; j < this->parity_column; j++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while (temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
		//end information permutation

		//begin information update at the check node
		for (j = 0; j < this->parity_row; j++)
		{
			//init
			for (int ii = 0; ii < this->q_ary; ii++)
			{
				total_mul[ii] = 1.0;
			}
			//cal the total multi
			temp_edge = this->tanner.m_row_head[j].left;
			while (temp_edge->col_num != -1)
			{
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					total_mul[ii] *= temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//cal the output from check node to variable nodes
			temp_edge = temp_edge->left;
			while (temp_edge->col_num != -1)
			{
				for (int ii = 0; ii < this->q_ary; ii++)
				{
					temp_edge->m_c2v[ii] = total_mul[ii] / temp_edge->m_v2c[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);
				temp_edge = temp_edge->left;
			}
		}
		//end information update at the check node
	}

	for (j = parity_column - 1; j >=0 ; j--)
	{
		uu_bit[j * m_degree] = uu_bit[j];
	}
	for (j = 0; j < parity_column; j++)
	{
		i = uu_bit[j * m_degree];
		for (int ii = 0; ii < m_degree; ii++)
		{
			uu_bit[j * m_degree + ii] = (i >> ii) & 1;
		}
	}
	return iter;
}

//*************************************************************************

int BCJRQSPA::FFTQSPA(double *rr, int *uu, double sigma)
{
	double curr_sum = 0.0;
	ldpc_edge *temp_edge;

	this->m_modem.Demapping(rr, channel_for_spa, sigma, parity_column);
//init metric from variable node to check node
     //metric from check node to variable node
	for(int i = 0 ; i < this->parity_column ; i ++)
	{
		temp_edge = this->tanner.m_col_head[i].up;
		while(temp_edge->row_num != -1)
		{
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_c2v[j] = 1.0;
			}
			for(int j = 0 ; j < this->tanner.q_ary ; j ++)
			{
				temp_edge->m_v2c[j] = 1.0;
			}
			temp_edge = temp_edge->up;
		}
	}
	int iter = 0;
	for(iter = 0 ; iter < this->max_iteration ; iter ++)
	{
//begin update the information from variable node to check node
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = channel_for_spa[j][ii];
			}
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *= temp_edge->m_c2v[ii];
				}
				
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					curr_sum += total_mul[ii];
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] /= curr_sum;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					if(total_mul[ii] < SMALLPROB)
						total_mul[ii] = SMALLPROB;
					else if(total_mul[ii] > 1 - SMALLPROB)
						total_mul[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
			curr_sum = -1.0;
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				if(curr_sum < total_mul[ii])
				{
					curr_sum = total_mul[ii];
					uu[j] = ii;
				}
			}
			curr_sum = 0.0;
			temp_edge = temp_edge->up;
			while(temp_edge->row_num != -1)
			{
				curr_sum = 0.0;
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = total_mul[ii] / temp_edge->m_c2v[ii];
					curr_sum += temp_edge->m_v2c[ii];
				}
				if(curr_sum == 0.0)
				{
					std::cout<<"DIVIDE BY 0"<<std::endl;
				}
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_v2c[ii] = temp_edge->m_v2c[ii] / curr_sum;
					if(temp_edge->m_v2c[ii] < SMALLPROB)
						temp_edge->m_v2c[ii] = SMALLPROB;
					else if(temp_edge->m_v2c[ii] > 1 - SMALLPROB)
						temp_edge->m_v2c[ii] = 1 - SMALLPROB;
				}
				temp_edge = temp_edge->up;
			}
		}
//end update information from variable node to check node

		if(this->is_codeword(uu))	
		{
			break;
		}
//begin information permutation
		for(int j = 0 ; j < this->parity_column ; j ++)
		{
			temp_edge = this->tanner.m_col_head[j].up;
			while(temp_edge->row_num != -1)
			{
				temp_edge->Norm(0);
				temp_edge->Permutation(0);
				//begin hadamard transform of information
				temp_edge->Transform(0);
				//end hadamard transform of information
				temp_edge = temp_edge->up;
			}
		}
//end information permutation

//begin information update at the check node
		for(int j = 0 ; j < this->parity_row ; j ++)
		{
			//init
			for(int ii = 0 ; ii < this->q_ary ; ii ++)
			{
				total_mul[ii] = 1.0;
			}
			//cal the total multi
			temp_edge = this->tanner.m_row_head[j].left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					total_mul[ii] *=temp_edge->m_v2c[ii];
				}
				temp_edge = temp_edge->left;
			}

			//cal the output from check node to variable nodes
			temp_edge = temp_edge->left;
			while(temp_edge->col_num != -1)
			{
				for(int ii = 0 ; ii < this->q_ary ; ii ++)
				{
					temp_edge->m_c2v[ii] = total_mul[ii] / temp_edge->m_v2c[ii];
				}
				temp_edge->Transform(1);
				temp_edge->Permutation(1);
				temp_edge->Norm(1);
				temp_edge = temp_edge->left;
			}
		}
//end information update at the check node
	}
	return iter;
}