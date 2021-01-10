/*
 * Read.cpp
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */

#include "Read.h"

void Read_DBG(struct DBG *graph,uint64_t kmer_len,uint8_t c)
{
	//1:in
	//2:in_c
	//4:out
	//8:out_c
	//16:uni
	//32:uni_c
	//64:br
	//128:br_c

	uint8_t x=1;
	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);
	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);
	struct para_getN tmp_inName;
	tmp_inName.kmerlen=kmer_len;
	tmp_inName.FileName=inputFileName;
	FILE* int_input;

	if(c&x==x)
	{
		tmp_inName.InOut_label=0;
		tmp_inName.isMiddle=0;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010in";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_inbranch=0;
		fseek(int_input,0,2);

		uint64_t x_in;
		total_hash_number_inbranch=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
		graph->p_inbranchNode=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number_inbranch*para.kmer64Len);
		fseek(int_input,0,0);
		x_in=fread(graph->p_inbranchNode,sizeof(uint64_t),total_hash_number_inbranch*para.kmer64Len,int_input);
		graph->p_inbranchNode_len=total_hash_number_inbranch;
		fclose(int_input);
	}
	if((c>>1)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=0;
		tmp_inName.isMiddle=2;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010inad";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_inad=0;
		fseek(int_input,0,2);
		uint8_t x_inad;
		total_hash_number_inad=ftell(int_input)/(sizeof(uint8_t)*para.kmer64Len);
		graph->p_inbranchNode_c=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number_inad*para.kmer64Len);
		fseek(int_input,0,0);
		x_inad=fread(graph->p_inbranchNode_c,sizeof(uint8_t),total_hash_number_inad*para.kmer64Len,int_input);
		fclose(int_input);

	}
	if((c>>2)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=1;
		tmp_inName.isMiddle=0;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010out";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_outbranch=0;
		fseek(int_input,0,2);
		uint64_t x_out;
		total_hash_number_outbranch=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
		graph->p_outbranchNode=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number_outbranch*para.kmer64Len);
		fseek(int_input,0,0);
		x_out=fread(graph->p_outbranchNode,sizeof(uint64_t),total_hash_number_outbranch*para.kmer64Len,int_input);
		graph->p_outbranchNode_len=total_hash_number_outbranch;
		fclose(int_input);
	}
	if((c>>3)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=1;
		tmp_inName.isMiddle=2;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010outad";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_outad=0;
		fseek(int_input,0,2);
		uint8_t x_outad;
		total_hash_number_outad=ftell(int_input)/(sizeof(uint8_t)*para.kmer64Len);
		graph->p_outbranchNode_c=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number_outad*para.kmer64Len);
		fseek(int_input,0,0);
		x_outad=fread(graph->p_outbranchNode_c,sizeof(uint8_t),total_hash_number_outad*para.kmer64Len,int_input);
		fclose(int_input);
	}
	if((c>>4)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=1;
		tmp_inName.isMiddle=0;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);
		inputFileName[7]='u';
		inputFileName[8]='n';
		inputFileName[9]='i';

	//	inputFileName="kmer010uni";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_nonbranch=0;
		fseek(int_input,0,2);
		uint64_t x_uni;
		total_hash_number_nonbranch=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
		graph->p_nonbranchNode=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number_nonbranch*para.kmer64Len);
		fseek(int_input,0,0);
		x_uni=fread(graph->p_nonbranchNode,sizeof(uint64_t),total_hash_number_nonbranch*para.kmer64Len,int_input);
		graph->p_nonbranchNode_len=total_hash_number_nonbranch;
		fclose(int_input);
	}
	if((c>>5)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=1;
		tmp_inName.isMiddle=2;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);
		inputFileName[7]='u';
		inputFileName[8]='n';
		inputFileName[9]='i';

	//	inputFileName="kmer010uniad";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_nonad=0;
		fseek(int_input,0,2);
		uint8_t x_uniad;
		total_hash_number_nonad=ftell(int_input)/(sizeof(uint8_t)*para.kmer64Len);
		graph->p_nonbranchNode_c=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number_nonad*para.kmer64Len);
		fseek(int_input,0,0);
		x_uniad=fread(graph->p_nonbranchNode_c,sizeof(uint8_t),total_hash_number_nonad*para.kmer64Len,int_input);
		fclose(int_input);
	}
	if((c>>6)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=0;
		tmp_inName.isMiddle=1;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_branch=0;
		fseek(int_input,0,2);
		uint64_t x_whole;
		total_hash_number_branch=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
		graph->p_branchNode=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number_branch*para.kmer64Len);
		fseek(int_input,0,0);
		x_whole=fread(graph->p_branchNode,sizeof(uint64_t),total_hash_number_branch*para.kmer64Len,int_input);
		graph->p_branchNode_len=total_hash_number_branch;
		fclose(int_input);
	}
	if((c>>6)&x==x)
	{
		tmp_inName.kmerlen=kmer_len;
		tmp_inName.FileName=inputFileName;
		tmp_inName.InOut_label=0;
		tmp_inName.isMiddle=3;
		tmp_inName.isposition=0;
		getFileName(tmp_inName);

	//	inputFileName="kmer010ad";
		int_input=fopen(inputFileName,"rb");
		uint64_t total_hash_number_branch_c=0;
		fseek(int_input,0,2);
		uint64_t x_whole_ad;
		total_hash_number_branch_c=ftell(int_input)/(sizeof(uint8_t)*para.kmer64Len);
		graph->p_branchNode_c=(uint8_t*)malloc(sizeof(uint8_t)*total_hash_number_branch_c*para.kmer64Len);
		fseek(int_input,0,0);
		x_whole_ad=fread(graph->p_branchNode_c,sizeof(uint8_t),total_hash_number_branch_c*para.kmer64Len,int_input);
		fclose(int_input);
	}
}
void Read_Pos_List(uint32_t** pl, uint32_t ** plen,uint32_t kmer_len)
{
	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);
	struct para_getN tmp_inName;
	tmp_inName.kmerlen=kmer_len;
	tmp_inName.FileName=inputFileName;
	tmp_inName.InOut_label=0;
	tmp_inName.isMiddle=1;
	tmp_inName.isposition=1000;
	getFileName(tmp_inName);
//	cout << inputFileName << endl;

	//read position list file
    FILE * p;
	p=fopen(inputFileName,"rb");
	free(inputFileName);
	uint64_t total_hash_number=0;
	fseek(p,0,2);
	uint32_t *h0;
	uint32_t x;
	total_hash_number=ftell(p)/(sizeof(uint32_t));
	h0=(uint32_t*)malloc(sizeof(uint32_t)*total_hash_number);
	fseek(p,0,0);
	x=fread(h0,sizeof(uint32_t),total_hash_number,p);

    uint32_t sum=0,sum_i=0,travel_i=0;
    while(travel_i<x)
    {
        sum+=h0[travel_i];
        sum_i++;
        travel_i=travel_i+h0[travel_i]+1;
    }
//    cout << sum_i << endl;
    uint32_t * pl_tmp;
	pl_tmp=(uint32_t*)malloc(sizeof(uint32_t)*sum);
    uint32_t * plen_tmp;
	plen_tmp=(uint32_t*)malloc(sizeof(uint32_t)*sum_i);

    uint32_t sum_ii=0;
    travel_i=0;
    uint32_t array_i=0;
    while(travel_i<x)
    {
        plen_tmp[sum_ii]=h0[travel_i];
        sum_ii++;
        for(uint32_t j=1;j<=h0[travel_i];j++)
        {
            pl_tmp[array_i]=h0[travel_i+j];
            array_i++;
        }
        travel_i=travel_i+h0[travel_i]+1;
    }

    for(uint32_t i=1;i<sum_i;i++)
    {
        plen_tmp[i]=plen_tmp[i]+plen_tmp[i-1];
    }

    free(h0);

    *pl=pl_tmp;
    *plen=plen_tmp;
}
void Read_Pos_List_only_length(uint32_t ** plen,uint32_t kmer_len)
{
	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);
	struct para_getN tmp_inName;
	tmp_inName.kmerlen=kmer_len;
	tmp_inName.FileName=inputFileName;
	tmp_inName.InOut_label=0;
	tmp_inName.isMiddle=1;
	tmp_inName.isposition=1000;
	getFileName(tmp_inName);
	cout << inputFileName << endl;

	//read position list file
    FILE * p;
	p=fopen(inputFileName,"rb");
	free(inputFileName);
	uint64_t total_hash_number=0;
	fseek(p,0,2);
	uint32_t *h0;
	uint32_t x;
	total_hash_number=ftell(p)/(sizeof(uint32_t));
	h0=(uint32_t*)malloc(sizeof(uint32_t)*total_hash_number);
	fseek(p,0,0);
	x=fread(h0,sizeof(uint32_t),total_hash_number,p);

    uint32_t sum=0,sum_i=0,travel_i=0;
    while(travel_i<x)
    {
        sum+=h0[travel_i];
        sum_i++;
        travel_i=travel_i+h0[travel_i]+1;
    }

    uint32_t * plen_tmp;
	plen_tmp=(uint32_t*)malloc(sizeof(uint32_t)*sum_i);

    uint32_t sum_ii=0;
    travel_i=0;
    uint32_t array_i=0;
    while(travel_i<x)
    {
        plen_tmp[sum_ii]=h0[travel_i];
        travel_i=travel_i+h0[travel_i]+1;
    }
    free(h0);

    *plen=plen_tmp;
}
void Free_DBG(struct DBG* g)
{
	if(g->p_branchNode!=NULL)
	{
		free(g->p_branchNode);
	}
	if(g->p_branchNode_c!=NULL)
	{
		free(g->p_branchNode_c);
	}
	if(g->p_inbranchNode!=NULL)
	{
		free(g->p_inbranchNode);
	}
	if(g->p_inbranchNode_c!=NULL)
	{
		free(g->p_inbranchNode_c);
	}
	if(g->p_outbranchNode!=NULL)
	{
		free(g->p_outbranchNode);
	}
	if(g->p_outbranchNode_c!=NULL)
	{
		free(g->p_outbranchNode_c);
	}
	if(g->p_nonbranchNode!=NULL)
	{
		free(g->p_nonbranchNode);
	}
	if(g->p_nonbranchNode_c!=NULL)
	{
		free(g->p_nonbranchNode_c);
	}

}
void Read_Mini_edge(uint32_t kmer_len,uint64_t ** p_e,uint64_t *pl_e)
{
	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);
	struct para_getN tmp_inName;
	tmp_inName.kmerlen=kmer_len;
	tmp_inName.FileName=inputFileName;
	tmp_inName.InOut_label=0;
	tmp_inName.isMiddle=1;
	tmp_inName.isposition=100;
	getFileName(tmp_inName);
//	cout << inputFileName << endl;

	//read position list file
    FILE * p;
	p=fopen(inputFileName,"rb");
	free(inputFileName);
	uint32_t total_hash_number=0;
	fseek(p,0,2);
	uint64_t *h0;
	uint64_t x;
	total_hash_number=ftell(p)/(sizeof(uint64_t));
	h0=(uint64_t*)malloc(sizeof(uint64_t)*total_hash_number);
	fseek(p,0,0);
	x=fread(h0,sizeof(uint64_t),total_hash_number,p);
	fclose(p);
	*p_e=h0;
	*pl_e=total_hash_number;
}
