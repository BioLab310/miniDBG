/*
 * Locate.h
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */

#ifndef LOCATE_H_
#define LOCATE_H_

#include "Basic.h"
#include "Binary_Search.h"
#include "Builder.h"

struct p_miniedge
{
	uint8_t type;
	int32_t offset;
	uint64_t p_mini_edge;
};
struct super_edge_para
{
    struct DBG graph;
    uint64_t **p_index_br;
    uint64_t **p_index_uni;
    uint64_t edge;
    uint64_t kmer_len;
    int32_t offset;
    uint32_t label;
    uint64_t plate;
};
struct super_next
{
    uint64_t mini_edge;
    uint8_t mini_edge_c;
    int is_super_edge;
    int32_t offset;
};
struct position_list
{
    uint32_t * p_pos_list;
    uint32_t p_pos_list_len;
    int32_t offset;
};
struct mini_edge_index
{
    uint64_t * p_mini_edge;
    uint64_t p_mini_edge_len;
    uint64_t **p_index_mini_edge;
    uint32_t* pl;
    uint32_t* plen;
};
struct Locate_Path_para
{
    struct DBG graph;
    uint64_t **p_index_br;
    uint64_t **p_index_uni;
//    uint64_t **p_index_mini;
    uint64_t kmer_len;
    char * path;
    uint32_t path_len;
    uint64_t plate;
    uint64_t plate1;
};

bool cmp(struct position_list a,struct position_list b);
uint32_t in_or_out(uint8_t c);
uint8_t get_node_ad(struct DBG graph,uint64_t n,uint64_t **p_index_br,uint64_t **p_index_uni,uint8_t &l);
void Gen_super_edge(struct super_edge_para tmp,struct super_next* r);
void Locate_Edge_All(struct super_edge_para tmp,struct mini_edge_index m,struct position_list *p);
void Locate_Edge(struct super_edge_para tmp,vector<struct p_miniedge> &p_mini_edge_vector);
void combine_pos_list(struct position_list a,struct position_list b,struct position_list* p1);
void combine_pos_list_multi(vector <struct position_list> a, struct position_list* p);
void get_mini_pos_list(struct mini_edge_index m,uint64_t e,struct position_list* p);
void mergy_pos_list(struct position_list p1,struct position_list p2,struct position_list *p3);
void mergy_pos_list_multi(vector <struct position_list> L, struct position_list *r);
uint32_t Locate_Path(struct Locate_Path_para p,struct mini_edge_index m,struct position_list* r);
#endif /* LOCATE_H_ */
