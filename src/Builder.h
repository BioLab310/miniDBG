/*
 * Builder.h
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */

#ifndef BUILDER_H_
#define BUILDER_H_

#include "Basic.h"
#include "Binary_Search.h"
#include "Read.h"
uint64_t get_next_node_bit(uint64_t p_branchNode,int nextchar,uint64_t kmer_len);
void gen_mini_edge_set(uint64_t **p_set,uint64_t *mini_edge_num, struct DBG graph,uint64_t kmer_len);
uint64_t get_mini_edge_bit(uint64_t p_branchNode,int nextchar,uint64_t kmer_len);
void Generate_Minimal_edges(uint64_t kmer_len,uint32_t threadnum);
void Generate_position_list_for_mini_edge(uint32_t kmer_len,uint64_t memory_limit,char * p_ref);
struct DBG
{
	uint64_t * p_outbranchNode;
	uint8_t  * p_outbranchNode_c;
	uint64_t   p_outbranchNode_len;

	uint64_t * p_inbranchNode;
	uint8_t  * p_inbranchNode_c;
	uint64_t   p_inbranchNode_len;

	uint64_t * p_nonbranchNode;
	uint8_t  * p_nonbranchNode_c;
	uint64_t   p_nonbranchNode_len;

	uint64_t * p_branchNode;
	uint8_t * p_branchNode_c;
	uint64_t  p_branchNode_len;
};
struct gen_mini_list_para
{
	uint32_t kmer_len;
	uint64_t memory_limit;
	char * p_edge;
	char * p_ref;
	uint32_t label_e_or_v;
};
typedef struct
{
	int whole_threadnum;
	int threadID;
	uint64_t *p_set;
	uint64_t mini_edge_num;
	struct DBG graph;
	uint64_t kmer_len;
} threadParm_t;


#endif /* BUILDER_H_ */
