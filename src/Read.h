/*
 * Read.h
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */

#ifndef READ_H_
#define READ_H_

#include "Basic.h"
#include "Builder.h"

void Read_DBG(struct DBG *graph,uint64_t kmer_len,uint8_t c);
void Read_Pos_List(uint32_t** pl, uint32_t ** plen,uint32_t kmer_len);
void Read_Pos_List_only_length(uint32_t ** plen,uint32_t kmer_len);
void Free_DBG(struct DBG* g);
void Read_Mini_edge(uint32_t kmer_len,uint64_t ** p_e,uint64_t *pl_e);
#endif /* READ_H_ */
