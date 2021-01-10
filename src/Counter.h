/*
 * Counter.h
 *
 *  Created on: Jan 2, 2021
 *      Author: bio
 */

#ifndef COUNTER_H_
#define COUNTER_H_
#include "Basic.h"
#include "Builder.h"

uint64_t Get_file_len(char * p);
void Counter(uint32_t bucket_num,uint32_t kmer_len);



#endif /* COUNTER_H_ */
