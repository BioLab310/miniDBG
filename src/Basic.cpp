/*
 * Basic.cpp
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */
#include "Basic.h"

uint32_t cmp256BitKmer(uint64_t*a,uint64_t*b,uint32_t len)
{
	uint32_t r=2;
	for(uint32_t i=0;i<len;i++)
	{
		if(a[i]<b[i])
		{
			r=0;
			break;
		}
		else
		{
			if(a[i]>b[i])
			{
				r=1;
				break;
			}
		}
	}
	return r;
}
void getFileName(struct para_getN tmp)
{
	uint32_t kmerlen=tmp.kmerlen;
	char* FileName=tmp.FileName;
	uint32_t label=tmp.InOut_label;
	uint32_t isMiddle=tmp.isMiddle;
	uint32_t isposition=tmp.isposition;

	if(isposition==0)
	{
		char init[5]="kmer";
		strcpy(FileName,init);
	}
	else if(isposition==10000)
	{
		char init[5]="unip";
		strcpy(FileName,init);
	}
	else if(isposition==100)
	{
		char init[5]="mini";
		strcpy(FileName,init);
	}
	else if(isposition==1000)
	{
		char init[5]="plst";
		strcpy(FileName,init);
	}
	else if(isposition==100000)
	{
		char init[5]="plin";
		strcpy(FileName,init);
	}
	else
	{
		char init[5]="posi";
		strcpy(FileName,init);
	}
	char swap[4];
	sprintf(swap,"%d",kmerlen);

	char in_add[3]="in";
	char out_add[4]="out";
	char ad[3]="ad";

	if(kmerlen<10)
	{
		FileName[4]='0';
		FileName[5]='0';
		strcpy(FileName+6,swap);
		if(isMiddle==1)
		{
			FileName[7]='\0';
		}
		else if(isMiddle==3)
		{
			strcpy(FileName+7,ad);
			FileName[9]='\0';
		}
		else
		{
			if(label==0)
			{
				strcpy(FileName+7,in_add);
				if(isMiddle==2)
				{
					strcpy(FileName+9,ad);
					FileName[11]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[9]='0';
						FileName[10]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[12]='\0';
					}
					else if(isposition<101)
					{
						FileName[9]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[12]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+9,swap);
						FileName[12]='\0';
					}
				}
				else
				{
					FileName[9]='\0';
				}
			}
			else
			{
				strcpy(FileName+7,out_add);
				if(isMiddle==2)
				{
					strcpy(FileName+10,ad);
					FileName[12]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[10]='0';
						FileName[11]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+12,swap);
						FileName[13]='\0';
					}
					else if(isposition<101)
					{
						FileName[10]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[13]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[13]='\0';
					}
				}
				else
				{
					FileName[10]='\0';
				}
			}
		}
	}
	else if(kmerlen<100)
	{
		FileName[4]='0';
		strcpy(FileName+5,swap);
		if(isMiddle==1)
		{
			FileName[7]='\0';
		}
		else if(isMiddle==3)
		{
			strcpy(FileName+7,ad);
			FileName[9]='\0';
		}
		else
		{
			if(label==0)
			{
				strcpy(FileName+7,in_add);
				if(isMiddle==2)
				{
					strcpy(FileName+9,ad);
					FileName[11]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[9]='0';
						FileName[10]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[12]='\0';
					}
					else if(isposition<101)
					{
						FileName[9]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[12]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+9,swap);
						FileName[12]='\0';
					}
				}
				else
				{
					FileName[9]='\0';
				}
			}
			else
			{
				strcpy(FileName+7,out_add);
				if(isMiddle==2)
				{
					strcpy(FileName+10,ad);
					FileName[12]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[10]='0';
						FileName[11]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+12,swap);
						FileName[13]='\0';
					}
					else if(isposition<101)
					{
						FileName[10]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[13]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[13]='\0';
					}
				}
				else
				{
					FileName[10]='\0';
				}
			}
		}
	}
	else
	{
		strcpy(FileName+4,swap);
		if(isMiddle==1)
		{
			FileName[7]='\0';
		}
		else if(isMiddle==3)
		{
			strcpy(FileName+7,ad);
			FileName[9]='\0';
		}
		else
		{
			if(label==0)
			{
				strcpy(FileName+7,in_add);
				if(isMiddle==2)
				{
					strcpy(FileName+9,ad);
					FileName[11]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[9]='0';
						FileName[10]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[12]='\0';
					}
					else if(isposition<101)
					{
						FileName[9]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[12]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+9,swap);
						FileName[12]='\0';
					}
				}
				else
				{
					FileName[9]='\0';
				}
			}
			else
			{
				strcpy(FileName+7,out_add);
				if(isMiddle==2)
				{
					strcpy(FileName+10,ad);
					FileName[12]='\0';
				}
				else if(isposition!=0)
				{
					if(isposition<11)
					{
						FileName[10]='0';
						FileName[11]='0';
						char swappo[2];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+12,swap);
						FileName[13]='\0';
					}
					else if(isposition<101)
					{
						FileName[10]='0';
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+11,swap);
						FileName[13]='\0';
					}
					else
					{
						char swappo[3];
						sprintf(swap,"%d",isposition-1);
						strcpy(FileName+10,swap);
						FileName[13]='\0';
					}
				}
				else
				{
					FileName[10]='\0';
				}
			}
		}
	}
}
void getRefFilePathes(char* pathFile, struct RefFilePath* p)
{
	struct timeval tvs,tve;
	gettimeofday(&tvs,NULL);
	ifstream int_input;

	cout << pathFile << endl;
	int_input.open(pathFile,ios::in);
	uint64_t total_hash_number=0;
	char hv_tmp[64];
	while(int_input.getline(hv_tmp,63))
	{
		total_hash_number++;
	}
	int_input.close();
	cout << "total ref number is: " <<total_hash_number << endl;

	p->pRefFilePath=(char**)malloc(sizeof(char*)*total_hash_number);
	for(uint32_t i=0;i<total_hash_number;i++)
	{
		p->pRefFilePath[i]=(char*)malloc(sizeof(char)*256);
	}

	uint32_t line_ref_path=0;
	int_input.open(pathFile,ios::in);
	while(int_input.getline(hv_tmp,63))
	{
		strcpy(p->pRefFilePath[line_ref_path],hv_tmp);
		line_ref_path++;
	}
	int_input.close();
	p->NumberOfPathes=line_ref_path;

	gettimeofday(&tve,NULL);
	double span = tve.tv_sec-tvs.tv_sec + (tve.tv_usec-tvs.tv_usec)/1000000.0;
}
void ReadSeq(char **seq1,uint64_t *seq_length,char* p_ref)
{
	uint32_t buffer_size=256;
	char buffer_line[256];
	memset(buffer_line,0,buffer_size);

	FILE *fp;
	fp = fopen(p_ref,"r+");
	if(fp==NULL)
	{
		cout <<"file can not be open!" << endl;
		return;
	}

	uint64_t total_size=0;
	fseek(fp,0,2);
	total_size=ftell(fp);

	char *seq;
	seq=(char*) malloc (sizeof(char)*total_size);

	fseek(fp,0,0);

	uint64_t num_of_N=0;
	uint64_t len=0;
	while (fgets(buffer_line,buffer_size-1,fp)!=NULL)
	{
		if(buffer_line[0]=='>')
			continue;
		else
		{
			for(uint32_t i=0;i<buffer_size;i++)
			{
				if(buffer_line[i]=='\n'||buffer_line[i]=='\0')
				{
					break;
				}
				if(buffer_line[i]>='a')
				{
					buffer_line[i]-=32;
				}
				if(buffer_line[i]!='A'&&buffer_line[i]!='C'&&buffer_line[i]!='G'&&buffer_line[i]!='T')
				{
					num_of_N++;
				}
				seq[len]=buffer_line[i];
				len++;
			}
		}
		memset(buffer_line,0,buffer_size);
	}
	*seq_length=len;
	*seq1=seq;
	cout << "the length of seq is: " << len << endl;
	cout << "the num of 'N' is: " << num_of_N <<endl;
}
uint64_t cal_hash_value_directly(char *seq,uint32_t len)
{
	uint64_t current=0;
	char *k_mer_temp=seq;
	for(uint32_t i=0;i<len;i++)
	{

		switch(k_mer_temp[i])
		{
			case 'A':
				current=current<<2;
				break;
			case 'C':
				current=current|1;
				current=current<<2;
				break;
			case 'G':
				current=current|2;
				current=current<<2;
				break;
			case 'T':
				current=current|3;
				current=current<<2;
				break;
			default:
				current=current<<2;
				break;
		}
	}
	current=current>>2;
	return current;
}
uint64_t cal_hash_value_indirectly(char *seq,uint64_t current,uint32_t len)
{

	char *k_mer_temp=seq;
	uint64_t high=0;
	for(uint32_t i=0;i<len*2-3;i++){
		high=high|1;
		high=high<<1;
	}
	high=high|1;
//	high=high<<(HashSize-2);
	current=current&high;//0000 0000 0000 0011 1111 1111 1111 1111
	current=current<<2;
	switch(k_mer_temp[len-1])
	{
		case 'A':
			//current=current|0;
			break;
		case 'C':
			current=current|1;
			//current=current<<2;
			break;
		case 'G':
			current=current|2;
			//current=current<<2;
			break;
		case 'T':
			current=current|3 ;
			//current=current<<2;
			break;
		default:
			break;
	}
	return current;
}
