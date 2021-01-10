/*
 * Builder.cpp
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */
#include "Builder.h"

uint64_t get_next_node_bit(uint64_t p_branchNode,int nextchar,uint64_t kmer_len)
{
	uint64_t nextnode;
	uint64_t high=0;
	for(uint32_t i=0;i<kmer_len-1;i++){
		high=high|3;
		high=high<<2;
	}
	high=high>>2;
	nextnode=p_branchNode&high;
	nextnode=nextnode<<2;
	switch(nextchar)
	{
		case 0:
			break;
		case 1:
			nextnode=nextnode|1;
			break;
		case 2:
			nextnode=nextnode|2;
			break;
		case 3:
			nextnode=nextnode|3;
			break;
		default:
			break;
	}
	return nextnode;
}
uint64_t get_mini_edge_bit(uint64_t p_branchNode,int nextchar,uint64_t kmer_len)
{
	uint64_t miniedge;
	miniedge=p_branchNode<<2;
	switch(nextchar)
	{
		case 0:
			break;
		case 1:
			miniedge=miniedge|1;
			break;
		case 2:
			miniedge=miniedge|2;
			break;
		case 3:
			miniedge=miniedge|3 ;
			break;
		default:
			break;
	}
	return miniedge;
}
void *ParFun(void *parm)
{
	threadParm_t *p = (threadParm_t *) parm;
	int whole_threadnum = p->whole_threadnum;

	struct bit256KmerPara para;
	uint64_t ukmerfindret;

	uint64_t mini_edge_n=0;
	uint64_t sizeof_mini=10000;
	uint64_t *mini_edge;
	mini_edge=(uint64_t*)malloc(sizeof(uint64_t)*sizeof_mini);

	uint64_t **bkmer_ptr_in;
	uint64_t **bkmer_ptr_uni;
	bkmer_ptr_in = Tgenerate_array(p->graph.p_inbranchNode_len);
	bkmer_ptr_uni = Tgenerate_array(p->graph.p_nonbranchNode_len);

	uint8_t tmp_single;
	tmp_single=1;

	uint64_t start,end,unit;
	unit=(p->graph.p_outbranchNode_len)/whole_threadnum;
	start=(p->threadID)*unit;
	end=(p->threadID+1)*unit;

	if(p->threadID==(whole_threadnum-1))
	{
		end=p->graph.p_outbranchNode_len;
	}

	for(uint64_t i=start;i<end;i++)
	{
		for(int j=0;j<=3;j++)
		{
			if((p->graph.p_outbranchNode_c[i])&(tmp_single<<j))
			{
				uint64_t nextnode;
				uint64_t miniedge;
				nextnode=get_next_node_bit(p->graph.p_outbranchNode[i],j,p->kmer_len);
				miniedge=get_mini_edge_bit(p->graph.p_outbranchNode[i],j,p->kmer_len);
				ukmerfindret = Tfind_arrindexN(bkmer_ptr_uni, p->graph.p_nonbranchNode, &nextnode, 1);
				while(ukmerfindret!=ULLONG_MAX)
				{
					if(((p->graph.p_nonbranchNode_c[ukmerfindret])&15)==0)
					{
						break;
					}
					else
					{
						for(int k=0;k<=3;k++)
						{
							if((p->graph.p_nonbranchNode_c[ukmerfindret])&(tmp_single<<k))
							{
								nextnode=get_next_node_bit(p->graph.p_nonbranchNode[ukmerfindret],k,p->kmer_len);
							}
						}
						ukmerfindret = Tfind_arrindexN(bkmer_ptr_uni, p->graph.p_nonbranchNode, &nextnode, 1);
					}
				}
				ukmerfindret = Tfind_arrindexN(bkmer_ptr_in, p->graph.p_inbranchNode, &nextnode, 1);
				if(ukmerfindret!=ULLONG_MAX)
				{
					if(mini_edge_n==sizeof_mini)
					{
						sizeof_mini+=10000;
						mini_edge=(uint64_t*)realloc(mini_edge,sizeof(uint64_t)*sizeof_mini);
					}
					mini_edge[mini_edge_n++]=miniedge;
				}
			}
		}
	}
	p->p_set=mini_edge;
	p->mini_edge_num=mini_edge_n;

}
void gen_mini_edge_set(struct DBG graph,uint64_t kmer_len)
{
	struct timeval tvs_total,tve_total;
	gettimeofday(&tvs_total,NULL);

	uint64_t ukmerfindret;

	uint64_t mini_edge_n=0;
	uint64_t sizeof_mini=10000;
	uint64_t *mini_edge;
	mini_edge=(uint64_t*)malloc(sizeof(uint64_t)*sizeof_mini);

	uint64_t **bkmer_ptr_in;
	uint64_t **bkmer_ptr_uni;
	bkmer_ptr_in = Tgenerate_array(graph.p_inbranchNode_len);
	bkmer_ptr_uni = Tgenerate_array(graph.p_nonbranchNode_len);

	uint8_t tmp_single;
	tmp_single=1;
	for(uint64_t i=0;i<graph.p_outbranchNode_len;i++)
	{
		for(int j=0;j<=3;j++)
		{
			if((graph.p_outbranchNode_c[i])&(tmp_single<<j))
			{
				uint64_t nextnode;
				uint64_t miniedge;
				nextnode=get_next_node_bit(graph.p_outbranchNode[i],j,kmer_len);
				miniedge=get_mini_edge_bit(graph.p_outbranchNode[i],j,kmer_len);
				ukmerfindret = Tfind_arrindexN(bkmer_ptr_uni, graph.p_nonbranchNode, &nextnode, 1);
				while(ukmerfindret!=ULLONG_MAX)
				{
					if(((graph.p_nonbranchNode_c[ukmerfindret])&15)==0)
					{
						break;
					}
					else
					{
						for(int k=0;k<=3;k++)
						{
							if((graph.p_nonbranchNode_c[ukmerfindret])&(tmp_single<<k))
							{
								nextnode=get_next_node_bit(graph.p_nonbranchNode[ukmerfindret],k,kmer_len);
							}
						}
						ukmerfindret = Tfind_arrindexN(bkmer_ptr_uni, graph.p_nonbranchNode, &nextnode, 1);
					}
				}
				ukmerfindret = Tfind_arrindexN(bkmer_ptr_in, graph.p_inbranchNode, &nextnode, 1);
				if(ukmerfindret!=ULLONG_MAX)
				{
					if(mini_edge_n==sizeof_mini)
					{
						sizeof_mini+=10000;
						mini_edge=(uint64_t*)realloc(mini_edge,sizeof(uint64_t)*sizeof_mini);
					}
					mini_edge[mini_edge_n++]=miniedge;
				}
			}
		}
	}
	gettimeofday(&tve_total,NULL);
	double span_check_total = tve_total.tv_sec-tvs_total.tv_sec + (tve_total.tv_usec-tvs_total.tv_usec)/1000000.0;
	cout << "the total cost time is: "<<span_check_total<<endl;

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

	FILE* int_input;
	int_input=fopen(inputFileName,"wb+");
	fwrite(mini_edge,sizeof(uint64_t),mini_edge_n,int_input);
	free(inputFileName);
	free(mini_edge);
	fclose(int_input);
}
void gen_mini_edge_set_pth(struct DBG graph,uint64_t kmer_len,uint32_t threadnum)
{
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	pthread_t threads[threadnum];
	threadParm_t threadParm[threadnum];
	int tn;

    for(tn=0;tn<threadnum;tn++)
    {
    	threadParm[tn].whole_threadnum = threadnum;
    	threadParm[tn].threadID = tn;
    	threadParm[tn].graph = graph;
    	threadParm[tn].kmer_len = kmer_len;
    	threadParm[tn].mini_edge_num = 0;
    	threadParm[tn].p_set = 0;
        pthread_create(&threads[tn],NULL,ParFun,(void *)&threadParm[tn]);
    }
    for(tn = 0; tn < threadnum; tn++)
    {
    	pthread_join(threads[tn], NULL);
    }

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

	FILE* int_input;
	int_input=fopen(inputFileName,"wb+");
    for(tn = 0; tn < threadnum; tn++)
    {
    	fwrite(threadParm[tn].p_set,sizeof(uint64_t),threadParm[tn].mini_edge_num,int_input);
    }
	free(inputFileName);
	fclose(int_input);
    for(tn = 0; tn < threadnum; tn++)
    {
    	free(threadParm[tn].p_set);
    }

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
    cout<<"gen_mini_edge_set_use_"<<threadnum<<"_threads:"<<span<<"without_savefiles"<<endl;
}
void Generate_Minimal_edges(uint64_t kmer_len,uint32_t threadnum)
{
	struct DBG graph;
	uint8_t c=63;
	Read_DBG(&graph,kmer_len,c);
	if(threadnum==1)
	{
		gen_mini_edge_set(graph,kmer_len);
	}
	else
	{
		gen_mini_edge_set_pth(graph,kmer_len,threadnum);
	}
	Free_DBG(&graph);
}
void gen_position_list(struct gen_mini_list_para tmp)
{
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);

	uint64_t total_length;
	total_length=0;
	uint64_t kmer_len=tmp.kmer_len;
	char * p_ref=tmp.p_ref;
	char * p_edge=tmp.p_edge;
	uint64_t memory_limit=tmp.memory_limit;
	uint64_t task_size=memory_limit*1024*1024*1024/8;

	cout << memory_limit << endl;
	cout << task_size << endl;

	struct bit256KmerPara para;
	para.kmer1Len=kmer_len*2;
	para.remainer1to64=para.kmer1Len%64;
	para.kmer64Len=para.kmer1Len/64+(para.remainer1to64?1:0);

	uint64_t *p_mini_edge_set;
	uint64_t mini_edge_num;

	FILE* int_input;
	int_input=fopen(p_edge,"rb");
	fseek(int_input,0,2);
	uint64_t x;
	mini_edge_num=ftell(int_input)/(sizeof(uint64_t)*para.kmer64Len);
	p_mini_edge_set=(uint64_t*)malloc(sizeof(uint64_t)*mini_edge_num*para.kmer64Len);
	fseek(int_input,0,0);
	x=fread(p_mini_edge_set,sizeof(uint64_t),mini_edge_num*para.kmer64Len,int_input);
	fclose(int_input);

	cout << "total number of the file is:" << mini_edge_num << endl;

	struct RefFilePath p_ref_path;
	getRefFilePathes(p_ref, &p_ref_path);

	char *seq;
	uint64_t seq_length;

	char* inputFileName;
	inputFileName=(char*)malloc(sizeof(char)*13);
	struct para_getN tmp_inName;
	if(tmp.label_e_or_v==1)
	{
		tmp_inName.kmerlen=kmer_len-1;
	}
	else
	{
		tmp_inName.kmerlen=kmer_len;
	}
	tmp_inName.FileName=inputFileName;
	tmp_inName.InOut_label=0;
	tmp_inName.isMiddle=1;
	tmp_inName.isposition=1000;
	getFileName(tmp_inName);
	cout << inputFileName << endl;
	cout << task_size << endl;

	FILE* int_output;
	int_output=fopen(inputFileName,"wb+");
	uint32_t **pos_list;
	uint32_t *pos_list_length;
	uint32_t *pos_list_size;
	uint64_t ukmerfindret;
	uint64_t **bkmer_ptr_miniedge;

	uint32_t loop_num;
	loop_num=mini_edge_num/task_size+1;
	for(uint32_t loop_i=0;loop_i<loop_num;loop_i++)
	{
		uint64_t task_size_cur;
		task_size_cur=task_size;
		if(loop_i==loop_num-1)
		{
			task_size_cur=mini_edge_num%task_size;
		}
		cerr << task_size_cur << endl;
		bkmer_ptr_miniedge = Tgenerate_array(task_size_cur);

		pos_list=(uint32_t**)malloc(sizeof(uint32_t*)*task_size_cur);
		pos_list_length=(uint32_t*)malloc(sizeof(uint32_t)*task_size_cur);
		pos_list_size=(uint32_t*)malloc(sizeof(uint32_t)*task_size_cur);

		for(uint32_t i=0;i<task_size_cur;i++)
		{
			pos_list_length[i]=0;
		}
		for(uint32_t i=0;i<task_size_cur;i++)
		{
			pos_list[i]=(uint32_t*)malloc(sizeof(uint32_t)*2);
			pos_list_size[i]=2;
		}

		uint64_t current;
		for(uint32_t ref_i=0;ref_i<p_ref_path.NumberOfPathes;ref_i++)
		{
			ReadSeq(&seq,&seq_length,p_ref_path.pRefFilePath[ref_i]);
			current=cal_hash_value_directly(seq,kmer_len);
			ukmerfindret = Tfind_arrindexN(bkmer_ptr_miniedge, p_mini_edge_set+loop_i*task_size, &current, 1);
			if(ukmerfindret!=ULLONG_MAX)
			{
				if(pos_list_length[ukmerfindret]==pos_list_size[ukmerfindret])
				{
					pos_list_size[ukmerfindret]+=2;
					pos_list[ukmerfindret]=(uint32_t*)realloc(pos_list[ukmerfindret],sizeof(uint32_t)*pos_list_size[ukmerfindret]);
				}
//				pos_list[ukmerfindret][pos_list_length[ukmerfindret]]=(0<<8)|ref_i;
				pos_list[ukmerfindret][pos_list_length[ukmerfindret]]=0;
				pos_list_length[ukmerfindret]++;
			}
			for(uint64_t i=1;i<seq_length-(kmer_len)+1;i++)
			{
				current=cal_hash_value_indirectly(seq+i,current,kmer_len);
				ukmerfindret = Tfind_arrindexN(bkmer_ptr_miniedge, p_mini_edge_set+loop_i*task_size, &current, 1);
				if(ukmerfindret!=ULLONG_MAX)
				{
					if(pos_list_length[ukmerfindret]==pos_list_size[ukmerfindret])
					{
						pos_list_size[ukmerfindret]+=2;
						pos_list[ukmerfindret]=(uint32_t*)realloc(pos_list[ukmerfindret],sizeof(uint32_t)*pos_list_size[ukmerfindret]);
					}
//					pos_list[ukmerfindret][pos_list_length[ukmerfindret]]=(i<<8)|ref_i;
					pos_list[ukmerfindret][pos_list_length[ukmerfindret]]=i;
					pos_list_length[ukmerfindret]++;
				}
			}
			free(seq);
		}

		for(uint32_t i=0;i<task_size_cur;i++)
		{
			fwrite(pos_list_length+i,sizeof(uint32_t),1,int_output);
			fwrite(pos_list[i],sizeof(uint32_t),pos_list_length[i],int_output);
			total_length+=pos_list_length[i];
		}

		cout << "the total length of the position list is :" << total_length << endl;
		for(uint32_t i=0;i<task_size_cur;i++)
		{
			free(pos_list[i]);
		}

		free(pos_list);
		free(pos_list_length);
		free(pos_list_size);
		Tfree_genarray(&bkmer_ptr_miniedge);
	}

	free(inputFileName);
	fclose(int_output);
	free(p_mini_edge_set);

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
	cout<<"gen_miniedge_position_list_time:"<<span<<endl;
}
void Generate_position_list_for_mini_edge(uint32_t kmer_len,uint64_t memory_limit,char * p_ref)
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

	struct gen_mini_list_para tmp;
	tmp.kmer_len=kmer_len+1;
	tmp.label_e_or_v=1;
	tmp.p_edge=inputFileName;
	tmp.memory_limit=memory_limit;
	tmp.p_ref=p_ref;

	gen_position_list(tmp);
	free(inputFileName);
}

