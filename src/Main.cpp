/*
 * Main.cpp
 *
 *  Created on: Jan 2, 2021
 *      Author: bio
 */

#include "Basic.h"
#include "Builder.h"
#include "Locate.h"
#include "Counter.h"

int main(int argc, char** argv)
{
	struct timeval tvs,tve;
	uint32_t kmer_len;
	uint32_t memory_limit;
	uint32_t thread_num=1;
	uint32_t l;
	char *ref_path;
	char method;
	for(uint32_t i=1;i<argc;i++)
	{
		if(argv[i][0]=='-'&&argv[i][1]=='L')//Len_kmer
		{
			kmer_len=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='t')//thread number
		{
			thread_num=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='r')//reference paths
		{
			ref_path=argv[i+1];
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='l')//reference paths
		{
			l=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='m')//memory limits
		{
			memory_limit=atoi(argv[i+1]);
		}
		else if(argv[i][0]=='-'&&argv[i][1]=='M')//memory limits
		{
			method=argv[i+1][0];
		}
	}

	if(method=='E')
	{
		//生成最小边集
		gettimeofday(&tvs,NULL);
		cout <<"start..."<<endl;
		Generate_Minimal_edges(kmer_len,thread_num);
		cout << "end..."<< endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"time cost for generating minimal edges is: "<<span<<endl;
	}
	else if(method=='P')
	{
		//生成最小边集的位置列表
		gettimeofday(&tvs,NULL);
		cout <<"start..."<<endl;
		Generate_position_list_for_mini_edge(kmer_len,memory_limit,ref_path);
		cout << "end..."<< endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"time cost for generating position list of minimal edges is: "<<span<<endl;
	}
	else if(method=='S')
	{
		//计算图的统计信息
		gettimeofday(&tvs,NULL);
		cout <<"start..."<<endl;
		uint32_t bucket=10;
		Counter(bucket,kmer_len);
		cout << "end..."<< endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"time cost for generating position list of minimal edges is: "<<span<<endl;

	}
	else if(method=='G')
	{
		//定位边
		gettimeofday(&tvs,NULL);
		cout <<"start..."<<endl;
		struct super_edge_para tmp;
		struct mini_edge_index m;
		struct position_list p;
		uint8_t c=240;
		Read_DBG(&(tmp.graph),kmer_len,c);
		tmp.kmer_len=kmer_len;
		//tmp.edge=//待定位的边
		tmp.offset=0;
		tmp.p_index_br=Tgenerate_array(tmp.graph.p_branchNode_len);
		tmp.p_index_uni=Tgenerate_array(tmp.graph.p_nonbranchNode_len);
		tmp.plate=0;
		for(uint32_t i=0;i<kmer_len-1;i++)
		{
			tmp.plate=tmp.plate|3;
			tmp.plate=tmp.plate<<2;
		}
		tmp.plate=tmp.plate>>2;
		Read_Mini_edge(kmer_len,&(m.p_mini_edge),&(m.p_mini_edge_len));
		m.p_index_mini_edge=Tgenerate_array(m.p_mini_edge_len);
		Read_Pos_List(&(m.pl), &(m.plen),kmer_len);
		///////////////////////////////////////////////
		Locate_Edge_All(tmp,m,&p);
		///////////////////////////////////////////////
		Free_DBG(&(tmp.graph));
		free(m.p_mini_edge);
		free(m.pl);
		free(m.plen);
		Tfree_genarray(&(tmp.p_index_br));
		Tfree_genarray(&(tmp.p_index_uni));
		Tfree_genarray(&(m.p_index_mini_edge));
		cout << "end..."<< endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"time cost for locating an edge is: "<<span<<endl;
	}
	else if(method=='H')
	{
		//定位路径
		gettimeofday(&tvs,NULL);
		cout <<"start..."<<endl;
		struct Locate_Path_para pp;
		uint8_t cc=240;
		Read_DBG(&(pp.graph),kmer_len,cc);
		pp.kmer_len=kmer_len;
		pp.p_index_br=Tgenerate_array(pp.graph.p_branchNode_len);
		pp.p_index_uni=Tgenerate_array(pp.graph.p_nonbranchNode_len);
		pp.plate=0;
		for(uint32_t i=0;i<kmer_len-1;i++)
		{
			pp.plate=pp.plate|3;
			pp.plate=pp.plate<<2;
		}
		pp.plate1=pp.plate+3;
		pp.plate=pp.plate>>2;
		struct mini_edge_index mm;
		Read_Mini_edge(kmer_len,&(mm.p_mini_edge),&(mm.p_mini_edge_len));
//		cout << mm.p_mini_edge_len << endl;
		mm.p_index_mini_edge=Tgenerate_array(mm.p_mini_edge_len);
		Read_Pos_List(&(mm.pl), &(mm.plen),kmer_len);
		struct position_list rr;
		///////////////////////////////////////////////
		struct RefFilePath p_ref_path;
		getRefFilePathes(ref_path, &p_ref_path);

		char *seq;
		uint64_t seq_length;

		for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
		{
			ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
			for(uint32_t j=5290367;j<seq_length;j++)
			{
				pp.path=seq+j;
				pp.path_len=kmer_len+l;
				Locate_Path(pp,mm,&rr);
				free(rr.p_pos_list);
			}

		}
		///////////////////////////////////////////////
		Free_DBG(&(pp.graph));
		free(mm.p_mini_edge);
		free(mm.pl);
		free(mm.plen);
		Tfree_genarray(&(pp.p_index_br));
		Tfree_genarray(&(pp.p_index_uni));
		Tfree_genarray(&(mm.p_index_mini_edge));
		cout << "end..."<< endl;
		gettimeofday(&tve,NULL);
		double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
		cout <<"time cost for locating a path is: "<<span<<endl;

	}
}
