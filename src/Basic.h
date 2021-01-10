/*
 * Basic.h
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */

#ifndef BASIC_H_
#define BASIC_H_

#include<algorithm>
#include<vector>
#include<cstdio>
#include <map>
#include<cstdlib>
#include<cstring>
#include <ctime>
#include<cmath>
#include <fstream>
#include <iostream>
#include "stdint.h"
#include "stdio.h"
#include "stdlib.h"
#include <pthread.h>
#include <sys/time.h>
using namespace std;

struct para_getN
{
	uint64_t kmerlen;
	char* FileName;
	uint32_t InOut_label;
	uint32_t isMiddle;
	uint32_t isposition;
};
struct bit256KmerPara
{
	uint32_t kmer1Len;
	uint32_t kmer64Len;
	uint32_t remainer1to64;
	uint64_t codefor1;
};
struct RefFilePath
{
	char **pRefFilePath;
	uint32_t NumberOfPathes;
};

uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len);
void getFileName(struct para_getN tmp);
void getRefFilePathes(char* pathFile, struct RefFilePath* p);
void ReadSeq(char **seq1,uint64_t *seq_length,char* p_ref);
uint64_t cal_hash_value_directly(char *seq,uint32_t len);
uint64_t cal_hash_value_indirectly(char *seq,uint64_t current,uint32_t len);
#endif /* BASIC_H_ */
