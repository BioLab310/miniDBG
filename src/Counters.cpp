/*
 * Counters.cpp
 *
 *  Created on: Jan 2, 2021
 *      Author: bio
 */

#include "Counter.h"

uint64_t Get_total_edge_num(struct DBG graph)
{
	uint64_t n=0;
	uint8_t c_tmp=1;
	for(uint64_t i=0;i<graph.p_branchNode_len;i++)
	{
		for(int j=0;j<=3;j++)
		{
			if(graph.p_branchNode_c[i]&(c_tmp<<j))
			{
				n++;
			}
		}
	}
	n+=graph.p_nonbranchNode_len;
	cout<<"total edge number is:"<<n<<endl;
	return n;
}
uint64_t Get_total_vertex_num(struct DBG graph)
{
	uint64_t n=0;
	n=graph.p_branchNode_len+graph.p_nonbranchNode_len;
	cout<<"total vertex number is:"<<n<<endl;;
	return n;
}
uint64_t Get_inbranch_edge_num(struct DBG graph)
{
	uint64_t n=0;
	uint8_t c_tmp=1;
	for(uint64_t i=0;i<graph.p_inbranchNode_len;i++)
	{
		for(int j=4;j<=7;j++)
		{
			if((graph.p_inbranchNode_c[i])&(c_tmp<<j))
			{
				n++;
			}
		}
	}
	cout<<"in-branching edge number is:"<<n<<endl;
	return n;
}
uint64_t Get_outbranch_edge_num(struct DBG graph)
{
	uint64_t n=0;
	uint8_t c_tmp=1;
	for(uint64_t i=0;i<graph.p_outbranchNode_len;i++)
	{
		for(int j=0;j<=3;j++)
		{
			if((graph.p_outbranchNode_c[i])&(c_tmp<<j))
			{
				n++;
			}
		}
	}
	cout<<"out-branching edge number is:"<<n<<endl;
	return n;
}
uint64_t Get_mini_edge_num(uint32_t kmer_len)
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
	cout << inputFileName << endl;

	uint64_t x;
	x=Get_file_len(inputFileName);
	x=x/(sizeof(uint64_t));
	free(inputFileName);

	cout<<"minimal edge number is:"<<x<<endl;
	return x;
}
void Get_distribution_of_list_length(uint32_t bucket_num,int kmer_len)
{
	uint64_t mini_edge_num;
	uint32_t *p_list_len;

	Read_Pos_List_only_length(&p_list_len,kmer_len);
	uint32_t maxnum=1;
	for(uint64_t i=0;i<mini_edge_num;i++)
	{
		if(p_list_len[i]>maxnum)
		{
			maxnum=p_list_len[i];
		}
	}

	uint32_t *bucket;
	uint32_t bucket_size=1;
	for(bucket_size=1;;bucket_size++)
	{
		if(pow(2,bucket_size)>maxnum)
			break;
	}
	bucket=(uint32_t*)malloc(sizeof(uint32_t)*(bucket_size+1));
	for(uint32_t i=0;i<bucket_size+1;i++)
	{
		bucket[i]=0;
	}
	uint32_t j;
	for(uint32_t i=0;i<mini_edge_num;i++)
	{
		for(j=0;j<=bucket_size;j++)
		{
			if(pow(2,j)>p_list_len[i])
				break;
		}
		if(p_list_len[i]!=maxnum)
		{
			bucket[j]++;
		}
	}
	for(uint32_t i=0;i<=bucket_size;i++)
	{
		cout<<i<<' '<<bucket[i] << endl;
	}
	free(p_list_len);
	free(bucket);
}
uint64_t Get_file_len(char * p)
{
	//read position list file
	uint64_t n;
    FILE * fp;
	fp=fopen(p,"rb");
	n=ftell(fp);
	fclose(fp);
	return n;
}
void Counter(uint32_t bucket_num,uint32_t kmer_len)
{
	uint8_t c=255;
	struct DBG graph;
	Read_DBG(&graph,kmer_len,c);
	uint64_t ne=0;
	uint64_t nv=0;
	uint64_t nie=0;
	uint64_t noe=0;
	uint64_t nme=0;
	ne=Get_total_edge_num(graph);
	cout << "the total edge number is:" << ne <<endl;
	nv=Get_total_vertex_num(graph);
	cout << "the total vertex number is:" << nv <<endl;
	nie=Get_inbranch_edge_num(graph);
	cout << "the total in-branching edge number is:" << nie <<endl;
	noe=Get_outbranch_edge_num(graph);
	cout << "the total out-branching edge number is:" << noe <<endl;
	nme=Get_mini_edge_num(kmer_len);
	cout << "the total minimal edge number is:" << nme <<endl;
	Get_distribution_of_list_length(bucket_num,kmer_len);
	Free_DBG(&graph);
}
