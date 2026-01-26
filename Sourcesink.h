//Writer : Terrence Zhao
//Time   :2009.04.08
//lastupdate : 2009.04.08
//Function:the source generator....
#ifndef _SOURCESINK_H
#define _SOURCESINK_H



#include "stdafx.h"
#include "util.h"
#include "Random.h"

class CSourceSink  
{
public:
	CSourceSink();
	virtual ~CSourceSink();
	void GetSymStr(int *uu, int q_ary, int len);
	void GetBitStr(int *uu, int len);
	//void GetBitStr(int *uu, int len , bool type);
	void ClrCnt();
	void CntErr(int *uu, int *uu_hat, int len, int flag);
	void CntErr_q(int *uu, int *uu_hat, int len, int flag);
	void PrintResult(FILE *fp);
public:
	double _tot_blk_num;//所有的block 数量
	double _tot_sym_num;//所有的bit 数量
	int _err_blk_num;//错误的block 数量
	int _err_sym_num;//错误的bit数量
	int _temp_err;  //当前block的错误bit数量

	double _ser;
	double _fer;
};


#endif