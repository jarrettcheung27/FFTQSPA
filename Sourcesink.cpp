//Writer : Terrence Zhao
//Time   :2009.04.08
//Function:the source generator....


#include "stdafx.h"
#include "SourceSink.h"

extern CLCRandNum rndGen0;
extern CWHRandNum rndGen1;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CSourceSink::CSourceSink()
{

}

CSourceSink::~CSourceSink()
{

}

/////////////////////////////////////////////////////////////////////////////
void CSourceSink::GetSymStr(int *uu, int q_ary, int len)
{
	double s;
	int num;
	for (int t = 0; t < len; t++)
	{
		s = rndGen1.Uniform() *q_ary;
		num = (int)s;
		num = num % q_ary;
		uu[t] = num;
	}
}

/////////////////////////////////////////////////////////////////////////////
//generate the source bits string////////
void CSourceSink::GetBitStr(int *uu, int len)
{
	for (int t = 0; t < len; t++)
		uu[t] = (rndGen1.Uniform() < 0.5?0:1);
}

//////////////////////////////////////////////////////////////////////////////////
void CSourceSink::ClrCnt()
{
	_tot_blk_num = 0;
	_tot_sym_num = 0;
	_err_blk_num = 0;
	_err_sym_num = 0;
}

//////////////////////////////////////////////
//用于计算每次编译码后的错误数量
void CSourceSink::CntErr(int *uu, int *uu_hat, int len, int flag)
{
	_temp_err = 0;
	for (int t = 0; t < len; t++){
		if (uu_hat[t] != uu[t])
			_temp_err++;
	}
	if (flag == 1){
		if (_temp_err > 0){
			_err_sym_num += _temp_err;
			_err_blk_num += 1;
		}
		
		_tot_blk_num += 1.0;
		_tot_sym_num += len;

		_ser = _err_sym_num /_tot_sym_num;
		_fer = _err_blk_num / _tot_blk_num;
	}

	return;
}

void CSourceSink::CntErr_q(int *uu, int *uu_hat, int len, int flag)
{
	_temp_err = 0;
	for (int t = 0; t < len; t++){
		if (uu_hat[t] == uu[t])
		{
			continue;
		}
		else
		{
			_temp_err++;
		}
	}
	if (flag == 1){
		if (_temp_err > 0){
			_err_sym_num += _temp_err;
			_err_blk_num += 1;
		}
		
		_tot_blk_num += 1.0;
		_tot_sym_num += len;

		_ser = _err_sym_num /_tot_sym_num;
		_fer = _err_blk_num / _tot_blk_num;
	}

	return;
}

///////////////////////////////////////////////////
void CSourceSink::PrintResult(FILE *fp)
{
	fprintf(fp, "\n tot_blk = %d: err_blk = %d: err_bit = %d: ser = %12.10f", 
		(int)_tot_blk_num, _err_blk_num, _err_sym_num, _ser);	
	return;
}