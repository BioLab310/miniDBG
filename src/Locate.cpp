/*
 * Locate.cpp
 *
 *  Created on: Jan 1, 2021
 *      Author: bio
 */
#include "Locate.h"
uint32_t in_or_out(uint8_t c)
{
    //0:不分叉
    //1:入分叉
    //2:出分叉
    //3:出入都分叉
    uint8_t a;

    a=c&15;
    if(a!=0&&a!=1&&a!=2&&a!=4&&a!=8)
    {
        a=c>>4;
        if(a!=0&&a!=1&&a!=2&&a!=4&&a!=8)
        {
            return 3;
        }
        else
        {
            return 2;
        }
    }
    else
    {
        a=c>>4;
        if(a!=0&&a!=1&&a!=2&&a!=4&&a!=8)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}
bool cmp(struct position_list a,struct position_list b)
{
    return a.p_pos_list_len<b.p_pos_list_len;
}
uint8_t get_node_ad(struct DBG graph,uint64_t n,uint64_t **p_index_br,uint64_t **p_index_uni,uint8_t &l)
{
    //l=2:没有出现在图上
    //l=1:分叉节点
    //l=0:不分叉节点
	uint64_t tmp;

	tmp = Tfind_arrindexN(p_index_br, graph.p_branchNode, &n, 1);
	if(tmp!=ULLONG_MAX)
    {
        l=1;
        return graph.p_branchNode_c[tmp];
    }
    else
    {
        tmp = Tfind_arrindexN(p_index_uni, graph.p_nonbranchNode, &n, 1);
        if(tmp!=ULLONG_MAX)
        {
            l=0;
            return graph.p_nonbranchNode_c[tmp];
        }
        else
        {
            l=2;
            return 0;
        }
    }
}
void Gen_super_edge(struct super_edge_para tmp,struct super_next* r)
{
    //如果是mini-super-edge上的边，返回mini-edge的第一条边和变异量
    //如果不是mini-super-edge上的边，返回一条边这条边是分叉节点的唯一入边或是唯一出边

    if(tmp.label==0)
    {
        uint64_t bits_next;
        bits_next=tmp.plate;
//        bits_next=tmp.plate<<2;
//        bits_next=bits_next|3;

        uint64_t next,next_edge;
        next_edge=tmp.edge;
        next=next_edge&bits_next;

        //一直向右扩展，直到分叉节点
        uint8_t tmp_l;
        uint8_t tmp_c;
        int32_t shift_next=tmp.offset;
        tmp_c=get_node_ad(tmp.graph,next,tmp.p_index_br,tmp.p_index_uni,tmp_l);
        while(tmp_l==0)
        {
            shift_next--;
            uint8_t tmp_cc=15;

            if((tmp_cc&tmp_c)==1)
            {
                next_edge=next<<2;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==2)
            {
                next_edge=(next<<2)|1;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==4)
            {
                next_edge=(next<<2)|2;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==8)
            {
                next_edge=(next<<2)|3;
                next=next_edge&bits_next;
            }
            tmp_c=get_node_ad(tmp.graph,next,tmp.p_index_br,tmp.p_index_uni,tmp_l);
        }

        //若不是入分叉
        uint32_t label_tmp;
        label_tmp=in_or_out(tmp_c);
        if(label_tmp==2)
        {
            r->is_super_edge=1;
            r->mini_edge=next_edge;
            r->mini_edge_c=tmp_c;
            r->offset=shift_next;
        }
        else//若是入分叉，从原始边向左扩展，直到分叉节点
        {
            uint64_t pre,pre_edge;
            pre_edge=tmp.edge;
            pre=pre_edge>>2;

            int32_t shift_pre=tmp.offset;
            tmp_c=get_node_ad(tmp.graph,pre,tmp.p_index_br,tmp.p_index_uni,tmp_l);
            while(tmp_l==0)
            {
                shift_pre++;
                uint64_t tmp_a;

                if((tmp_c>>4)==1)
                {
                    pre_edge=pre;
                    pre=pre_edge>>2;
                }
                else if((tmp_c>>4)==2)
                {
                    tmp_a=1;
                    pre_edge=pre|(tmp_a<<(2*(tmp.kmer_len)));
                    pre=pre_edge>>2;
                }
                else if((tmp_c>>4)==4)
                {
                    tmp_a=2;
                    pre_edge=pre|(tmp_a<<(2*(tmp.kmer_len)));
                    pre=pre_edge>>2;
                }
                else if((tmp_c>>4)==8)
                {
                    tmp_a=3;
                    pre_edge=pre|(tmp_a<<(2*(tmp.kmer_len)));
                    pre=pre_edge>>2;
                }
                tmp_c=get_node_ad(tmp.graph,pre,tmp.p_index_br,tmp.p_index_uni,tmp_l);
            }

            //若是不是出分叉
            label_tmp=in_or_out(tmp_c);
            if(label_tmp==1)
            {
                r->is_super_edge=-1;
                r->mini_edge=pre_edge;
                r->mini_edge_c=tmp_c;
                r->offset=shift_pre;
            }
            else//若是出分叉
            {
                r->is_super_edge=0;
                r->mini_edge=pre_edge;
                r->mini_edge_c=tmp_c;
                r->offset=shift_pre;
            }
        }
    }
    else if(tmp.label==1)
    {
        uint64_t bits_next=0;
        for(int i=0;i<tmp.kmer_len-1;i++)
        {
            bits_next=bits_next|3;
            bits_next=bits_next<<2;
        }
        bits_next=bits_next|3;

        uint64_t next,next_edge;
        next_edge=tmp.edge;
        next=next_edge&bits_next;

        //一直向右扩展，直到分叉节点
        uint8_t tmp_l;
        uint8_t tmp_c;
        int32_t shift_next=tmp.offset;
        tmp_c=get_node_ad(tmp.graph,next,tmp.p_index_br,tmp.p_index_uni,tmp_l);
        while(tmp_l==0)
        {
            shift_next--;
            uint8_t tmp_cc=15;

            if((tmp_cc&tmp_c)==1)
            {
                next_edge=next<<2;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==2)
            {
                next_edge=(next<<2)|1;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==4)
            {
                next_edge=(next<<2)|2;
                next=next_edge&bits_next;
            }
            else if((tmp_cc&tmp_c)==8)
            {
                next_edge=(next<<2)|3;
                next=next_edge&bits_next;
            }
            tmp_c=get_node_ad(tmp.graph,next,tmp.p_index_br,tmp.p_index_uni,tmp_l);
        }

        //若不是入分叉
        uint32_t label_tmp;
        label_tmp=in_or_out(tmp_c);
        if(label_tmp==2)
        {
            r->is_super_edge=1;
            r->mini_edge=next_edge;
            r->mini_edge_c=tmp_c;
            r->offset=shift_next;
        }
        else//若是入分叉，从原始边向左扩展，直到分叉节点
        {
            r->is_super_edge=0;
            r->mini_edge=tmp.edge;
            r->mini_edge_c=tmp_c;
            r->offset=tmp.offset;
        }
    }
    else if(tmp.label==2)
    {
    	uint8_t tmp_l;
        uint8_t tmp_c;
        uint64_t pre,pre_edge;
        pre_edge=tmp.edge;
        pre=pre_edge>>2;

        int32_t shift_pre=tmp.offset;
        tmp_c=get_node_ad(tmp.graph,pre,tmp.p_index_br,tmp.p_index_uni,tmp_l);
        while(tmp_l==0)
        {
            shift_pre++;
            uint64_t tmp_a;
            if((tmp_c>>4)==1)
            {
                pre_edge=pre;
                pre=pre_edge>>2;
            }
            else if((tmp_c>>4)==2)
            {
                tmp_a=1;
                pre_edge=pre|(tmp_a<<(2*tmp.kmer_len));
                pre=pre_edge>>2;
            }
            else if((tmp_c>>4)==4)
            {
                tmp_a=2;
                pre_edge=pre|(tmp_a<<(2*tmp.kmer_len));
                pre=pre_edge>>2;
            }
            else if((tmp_c>>4)==8)
            {
                tmp_a=3;
                pre_edge=pre|(tmp_a<<(2*tmp.kmer_len));
                pre=pre_edge>>2;
            }
            tmp_c=get_node_ad(tmp.graph,pre,tmp.p_index_br,tmp.p_index_uni,tmp_l);
         }

        //若是不是出分叉
        uint32_t label_tmp;
        label_tmp=in_or_out(tmp_c);
        if(label_tmp==1)
        {
            r->is_super_edge=-1;
            r->mini_edge=pre_edge;
            r->mini_edge_c=tmp_c;
            r->offset=shift_pre;
        }
        else//若是出分叉
        {
            r->is_super_edge=0;
            r->mini_edge=pre_edge;
            r->mini_edge_c=tmp_c;
            r->offset=shift_pre;
        }
    }
}
void Locate_Edge_All(struct super_edge_para tmp,struct mini_edge_index m,struct position_list *p)
{
	vector <struct p_miniedge> p_mini_edge_vector;
	Locate_Edge(tmp,p_mini_edge_vector);
    struct position_list p_tmp;
    vector <struct position_list> a;
    for(uint32_t i=0;i<p_mini_edge_vector.size();i++)
    {
        get_mini_pos_list(m,p_mini_edge_vector[i].p_mini_edge,&p_tmp);
        struct position_list p_tmpp;
        p_tmpp.offset=0;
        p_tmpp.p_pos_list_len=p_tmp.p_pos_list_len;
        p_tmpp.p_pos_list=(uint32_t*)malloc(sizeof(uint32_t)*p_tmp.p_pos_list_len);
        for(uint32_t j=0;j<p_tmp.p_pos_list_len;j++)
        {
            p_tmpp.p_pos_list[j]=p_tmp.p_pos_list[j]+p_mini_edge_vector[i].offset;
        }
        a.push_back(p_tmpp);
    }
    combine_pos_list_multi(a,p);
    for(uint32_t i=0;i<a.size();i++)
    {
    	free(a[i].p_pos_list);
    	a[i].p_pos_list=NULL;
    }
}
void Locate_Edge(struct super_edge_para tmp,vector<struct p_miniedge> &p_mini_edge_vector)
{
//	tmp.plate=(tmp.plate<<2)+3;
    struct p_miniedge p_mini_edge_tmp;
    //label=0,first call
    struct super_edge_para tmp_vec;
    if(tmp.label==0)
    {
        struct super_next tmp_next;
        Gen_super_edge(tmp,&tmp_next);
        if(tmp_next.is_super_edge==0)
        {
            p_mini_edge_tmp.offset=tmp_next.offset;
            p_mini_edge_tmp.p_mini_edge=tmp_next.mini_edge;
            p_mini_edge_vector.push_back(p_mini_edge_tmp);
            return;
        }
        else if(tmp_next.is_super_edge==1)
        {
            uint8_t tmp_c;
            uint64_t next=(tmp_next.mini_edge&tmp.plate)<<2;
            tmp_vec=tmp;
            tmp_vec.offset=tmp_next.offset-1;
            tmp_vec.label=1;
            tmp_c=1;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=2;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+1;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=4;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+2;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=8;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+3;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
        }
        else if(tmp_next.is_super_edge==-1)
        {
            uint8_t tmp_c;
            uint64_t next=tmp_next.mini_edge>>2;
            tmp_vec=tmp;
            tmp_vec.offset=tmp_next.offset+1;
            tmp_vec.label=2;
            tmp_c=16;
            uint64_t tmp_in;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=32;
			tmp_in=1;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=64;
			tmp_in=2;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=128;
			tmp_in=3;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
        }
    }
    else if(tmp.label==1)
    {
        //向右找分叉节点
        struct super_next tmp_next;
        Gen_super_edge(tmp,&tmp_next);
        if(tmp_next.is_super_edge==0)
        {
            p_mini_edge_tmp.offset=tmp_next.offset;
            p_mini_edge_tmp.p_mini_edge=tmp_next.mini_edge;
            p_mini_edge_vector.push_back(p_mini_edge_tmp);
            return;
        }
        else if(tmp_next.is_super_edge==1)
        {
            uint8_t tmp_c;
            uint64_t next=(tmp_next.mini_edge&tmp.plate)<<2;
            tmp_vec=tmp;
            tmp_vec.offset=tmp_next.offset-1;
            tmp_vec.label=1;
            tmp_c=1;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=2;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+1;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=4;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+2;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=8;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next+3;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
        }

    }
    else if(tmp.label==2)
    {
        //向左找分叉节点
        struct super_next tmp_next;
        Gen_super_edge(tmp,&tmp_next);
        if(tmp_next.is_super_edge==0)
        {
            p_mini_edge_tmp.offset=tmp_next.offset;
            p_mini_edge_tmp.p_mini_edge=tmp_next.mini_edge;
            p_mini_edge_vector.push_back(p_mini_edge_tmp);
            return;
        }
        else if(tmp_next.is_super_edge==-1)
        {
            uint8_t tmp_c;
            uint64_t next=tmp_next.mini_edge>>2;
            tmp_vec=tmp;
            tmp_vec.offset=tmp_next.offset+1;
            tmp_vec.label=2;
            tmp_c=16;
            uint64_t tmp_in;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next;
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=32;
			tmp_in=1;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=64;
			tmp_in=2;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
			tmp_c=128;
			tmp_in=3;
			if((tmp_c&tmp_next.mini_edge_c)==tmp_c)
			{
                tmp_vec.edge=next|(tmp_in<<(2*tmp.kmer_len));
				Locate_Edge(tmp_vec,p_mini_edge_vector);
			}
        }
    }
}
void combine_pos_list(struct position_list a,struct position_list b,struct position_list* p1)
{
    uint32_t l;
    l=a.p_pos_list_len+b.p_pos_list_len;
    uint32_t * p;
    p=(uint32_t*)malloc(sizeof(uint32_t)*l);

    uint32_t x=0,y=0,z=0;
    while(x<a.p_pos_list_len&&y<b.p_pos_list_len)
    {
        if(a.p_pos_list[x]+a.offset<b.p_pos_list[y]+b.offset)
        {
            p[z]=a.p_pos_list[x]+a.offset;
            x++;
            z++;
        }
        else if(a.p_pos_list[x]+a.offset>b.p_pos_list[y]+b.offset)
        {
            p[z]=b.p_pos_list[y]+b.offset;
            y++;
            z++;
        }
        else
        {
            p[z]=a.p_pos_list[x]+a.offset;
            x++;
            y++;
            z++;
        }
    }
    if(x<a.p_pos_list_len)
    {
        for(uint32_t i=x;i<a.p_pos_list_len;i++)
        {
            p[z]=a.p_pos_list[x]+a.offset;
            x++;
            z++;
        }
    }
    else if(y<b.p_pos_list_len)
    {
        for(uint32_t i=y;i<b.p_pos_list_len;i++)
        {
            p[z]=b.p_pos_list[y]+b.offset;
            y++;
            z++;
        }
    }

    p=(uint32_t*)realloc(p,sizeof(uint32_t)*z);

    p1->p_pos_list=p;
    p1->p_pos_list_len=z;
    p1->offset=0;
}
void combine_pos_list_multi(vector <struct position_list> a, struct position_list* p)
{
    struct position_list c;
    struct position_list c_tmp;

    if(a.size()<2)
    {
    	c.p_pos_list=(uint32_t*)malloc(sizeof(uint32_t)*a[0].p_pos_list_len);
    	c.p_pos_list_len=a[0].p_pos_list_len;
    	for(uint32_t i=0;i<a[0].p_pos_list_len;i++)
    	{
    		c.p_pos_list[i]=a[0].p_pos_list[i];
    	}
    	c.offset=a[0].offset;
    	p->offset=c.offset;
    	p->p_pos_list=c.p_pos_list;
    	p->p_pos_list_len=c.p_pos_list_len;
        return;
    }
    else
    {
        combine_pos_list(a[0],a[1],&c);
        for(uint32_t i=2;i<a.size();i++)
        {
            c_tmp=c;
            combine_pos_list(a[i],c_tmp,&c);
            free(c_tmp.p_pos_list);
        }
        p->p_pos_list=c.p_pos_list;
        p->p_pos_list_len=c.p_pos_list_len;
        p->offset=c.offset;
    }
}
void get_mini_pos_list(struct mini_edge_index m,uint64_t e,struct position_list* p)
{
//	cout << e << endl;
   	uint64_t tmp;
	tmp = Tfind_arrindexN(m.p_index_mini_edge, m.p_mini_edge, &e, 1);
//	cout << m.p_mini_edge[tmp]<< endl;

	if(tmp==ULLONG_MAX)
    {
        cout << "error: get minimal edge position list failed!" << endl;
    }
    else
    {
        if(tmp==0)
        {
            p->p_pos_list=m.pl;
            p->p_pos_list_len=m.plen[0];
        }
        else
        {
            p->p_pos_list=m.pl+m.plen[tmp-1];
            p->p_pos_list_len=m.plen[tmp]-m.plen[tmp-1];
        }
    }
}
void mergy_pos_list(struct position_list p1,struct position_list p2,struct position_list *p3)
{
    uint32_t* r;
    uint32_t* q1;
    uint32_t* q2;
    uint32_t s1,lr;
    uint32_t I1,I2;
    uint32_t O1,O2;

    if(p1.p_pos_list_len>=p2.p_pos_list_len)
    {
        q1=p1.p_pos_list;
        I1=p1.p_pos_list_len;
        O1=p1.offset;
        q2=p2.p_pos_list;
        I2=p2.p_pos_list_len;
        O2=p2.offset;
        s1=p1.p_pos_list_len/p2.p_pos_list_len;
    }
    else
    {
        q2=p1.p_pos_list;
        I2=p1.p_pos_list_len;
        O2=p1.offset;
        q1=p2.p_pos_list;
        I1=p2.p_pos_list_len;
        O1=p2.offset;
        s1=p2.p_pos_list_len/p1.p_pos_list_len;
    }
	r=(uint32_t*)malloc(sizeof(uint32_t)*I2);
    uint32_t i1,i2;
    i1=0;
    i2=0;
    lr=0;
    while(i1<I1&&i2<I2)
    {
        if(q1[i1]+O1==q2[i2]+O2)
        {
            r[lr]=q1[i1]+O1;
            lr++;
            i1++;
            i2++;
        }
        else if(q1[i1]+O1>q2[i2]+O2)
        {
            i2++;
        }
        else
        {
//            if(q1[i1+s1]-O1<q2[i2]-O2)//数组越界
//            {
//                i1+=s1;
//            }
//            else
//            {
			i1++;
//            }
        }
    }
    if(lr==0)
    {
    	p3->p_pos_list=NULL;
    	p3->p_pos_list_len=0;
    	return;
    }
    else if(lr!=I2)
    {
    	r=(uint32_t*)realloc(r,sizeof(uint32_t)*lr);
    }
	p3->p_pos_list=r;
	p3->p_pos_list_len=lr;
	p3->offset=0;
}
void mergy_pos_list_multi(vector <struct position_list> L, struct position_list *r)
{
	if(L.size()==1)
	{
	    struct position_list c;
    	c.p_pos_list=(uint32_t*)malloc(sizeof(uint32_t)*L[0].p_pos_list_len);
    	c.p_pos_list_len=L[0].p_pos_list_len;
    	for(uint32_t i=0;i<L[0].p_pos_list_len;i++)
    	{
    		c.p_pos_list[i]=L[0].p_pos_list[i];
    	}
    	c.offset=L[0].offset;
		r->offset=c.offset;
		r->p_pos_list=c.p_pos_list;
		r->p_pos_list_len=c.p_pos_list_len;
		return;
	}
    sort(L.begin(),L.end(),cmp);
    struct position_list r_tmp;
    mergy_pos_list(L[0],L[1],r);
    if(r->p_pos_list_len==0)
    {
    	r->p_pos_list_len=0;
        r->p_pos_list=NULL;
        return;
    }

    for(uint32_t i=2;i<L.size();i++)
    {
        if(r->p_pos_list_len==0)
        {
//        	cout << i << endl;
        	r->p_pos_list_len=0;
            r->p_pos_list=NULL;
            return;
        }
        else
        {
            r_tmp.offset=r->offset;
            r_tmp.p_pos_list=r->p_pos_list;
            r_tmp.p_pos_list_len=r->p_pos_list_len;
            mergy_pos_list(L[i],r_tmp,r);
            free(r_tmp.p_pos_list);
        }
    }
}
uint32_t Locate_Path(struct Locate_Path_para p,struct mini_edge_index m,struct position_list* r)
{
	if(p.path_len<=p.kmer_len)
	{
		cout << "error: the input path is too short to be located!" << endl;
		return 0;
	}
	if(p.path_len==p.kmer_len+1)
	{
		cout << "this path is just an edge!" << endl;
		struct super_edge_para tmp;
		tmp.graph=p.graph;
		tmp.kmer_len=p.kmer_len;
		tmp.label=0;
		tmp.offset=0;
		tmp.p_index_br=p.p_index_br;
		tmp.p_index_uni=p.p_index_uni;
		tmp.plate=p.plate1;
		tmp.edge=cal_hash_value_directly(p.path,p.kmer_len+1);
		Locate_Edge_All(tmp,m,r);
		return 1;
	}
    //return 0, if the path is a false path;
    vector<struct p_miniedge> v_mini_edge;

    uint64_t kmer_tmp;
    uint8_t tmp_c;
    uint8_t tmp_l;
    uint8_t *p_c;
    p_c=(uint8_t*)malloc(sizeof(uint8_t)*(p.path_len-p.kmer_len+1));
    uint64_t *p_e;
    p_e=(uint64_t*)malloc(sizeof(uint64_t)*(p.path_len-p.kmer_len));
    kmer_tmp=cal_hash_value_directly(p.path,p.kmer_len);
    for(uint32_t i=0;i<p.path_len-p.kmer_len;i++)
    {
        tmp_c=get_node_ad(p.graph,kmer_tmp,p.p_index_br,p.p_index_uni,tmp_l);
        if(tmp_l==2)
        {
            r->p_pos_list=NULL;
            r->p_pos_list_len=0;
            return 0;//这个顶点不在图中出现的情况
        }
        else
        {
            uint8_t tmp_cc=0;
            uint64_t tmp_ui=0;
            if(p.path[i+p.kmer_len]=='A')
            {
                tmp_cc=1;
                tmp_ui=0;
            }
            else if(p.path[i+p.kmer_len]=='C')
            {
                tmp_cc=2;
                tmp_ui=1;
            }
            else if(p.path[i+p.kmer_len]=='G')
            {
                tmp_cc=4;
                tmp_ui=2;
            }
            else if(p.path[i+p.kmer_len]=='T')
            {
                tmp_cc=8;
                tmp_ui=3;
            }
            if((tmp_c&tmp_cc)==0)
            {
                r->p_pos_list=NULL;
                r->p_pos_list_len=0;
                return 0;//i顶点开始的这条边在图上没有出现的情况
            }
            else
            {
                //顶点和边都出现了，保存顶点邻接信息，边，计算下个顶点
                p_c[i]=tmp_c;
                p_e[i]=(kmer_tmp<<2)+tmp_ui;
                kmer_tmp=(kmer_tmp&p.plate)<<2;
                kmer_tmp=kmer_tmp+tmp_ui;
            }
        }
    }
    tmp_c=get_node_ad(p.graph,kmer_tmp,p.p_index_br,p.p_index_uni,tmp_l);
    if(tmp_l==2)
    {
        r->p_pos_list=NULL;
        r->p_pos_list_len=0;
        return 0;//这个顶点不在图中出现的情况
    }
    else
    {
        p_c[p.path_len-p.kmer_len]=tmp_c;
    }

    int32_t start,end;
    int32_t label_start;
    label_start=-1;
    start=0;
    end=0;
    struct p_miniedge tmp_mini;
    for(uint32_t i=0;i<p.path_len-p.kmer_len+1;i++)
    {
        uint32_t label;
        label=in_or_out(p_c[i]);
        if(label==0)
        {
            continue;
        }
        if(label_start==-1)
        {
            if(i!=0&&(label==1||label==3))
            {
                tmp_mini.p_mini_edge=p_e[0];
                tmp_mini.offset=-start;
                tmp_mini.type=1;
                v_mini_edge.push_back(tmp_mini);
            }
        }
        else if(label_start==1)
        {
            if(label==1||label==3)
            {
                end=i;
                tmp_mini.p_mini_edge=p_e[start];
                tmp_mini.offset=-start;
                tmp_mini.type=0;
                v_mini_edge.push_back(tmp_mini);
                label_start=0;
            }
        }
        if((label==2||label==3)&&i!=p.path_len-p.kmer_len)
        {
            label_start=1;
            start=i;
        }
    }
    if((label_start==-1||label_start==1)&&end!=p.path_len-p.kmer_len&&start<p.path_len-p.kmer_len)
    {
        tmp_mini.p_mini_edge=p_e[start];
        tmp_mini.offset=-start;
        tmp_mini.type=1;
        v_mini_edge.push_back(tmp_mini);
    }

    vector <struct position_list> a;
    struct position_list a_tmp;
    for(uint32_t i=0;i<v_mini_edge.size();i++)
    {
        if(v_mini_edge[i].type!=1)
        {
            get_mini_pos_list(m,v_mini_edge[i].p_mini_edge,&a_tmp);
            if(a_tmp.p_pos_list_len!=0)
            {
            	a_tmp.offset=v_mini_edge[i].offset;
				a.push_back(a_tmp);
            }
        }
    }

    struct position_list pp;
    if(a.size()!=0)
    {
		mergy_pos_list_multi(a,&pp);
		if(pp.p_pos_list_len==0)
		{
			r->p_pos_list_len=0;
			r->p_pos_list=NULL;
			return 0;
		}
		else
		{
			struct position_list p_tmpp;
			p_tmpp.p_pos_list=NULL;
			for(uint32_t i=0;i<v_mini_edge.size();i++)
			{
				if(v_mini_edge[i].type==1)
				{
					struct position_list p_tmp;
					p_tmp.p_pos_list=NULL;
					struct super_edge_para tmp;
					tmp.edge=v_mini_edge[i].p_mini_edge;
					tmp.graph=p.graph;
					tmp.kmer_len=p.kmer_len;
					tmp.label=0;
					tmp.offset=v_mini_edge[i].offset;
					tmp.plate=p.plate1;
					tmp.p_index_br=p.p_index_br;
					tmp.p_index_uni=p.p_index_uni;

					Locate_Edge_All(tmp,m,&p_tmp);

					if(p_tmp.p_pos_list_len==0)
					{
						if(r->p_pos_list!=NULL)
						{
							free(r->p_pos_list);
						}
						r->p_pos_list_len=0;
						r->p_pos_list=NULL;
						return 0;
					}
					else
					{
						p_tmpp=pp;
						mergy_pos_list(p_tmp,p_tmpp,&pp);
						free(p_tmp.p_pos_list);
						free(p_tmpp.p_pos_list);
						p_tmp.p_pos_list=NULL;
						p_tmp.p_pos_list_len=0;
						p_tmpp.p_pos_list=NULL;
						p_tmpp.p_pos_list_len=0;
					}
					if(pp.p_pos_list_len==0)
					{
						if(pp.p_pos_list!=NULL)
						{
							free(pp.p_pos_list);
						}
						pp.p_pos_list=NULL;
						return 0;
					}
				}
			}
		}
    }
    else
    {
    	vector <struct position_list> pa;
		for(uint32_t i=0;i<v_mini_edge.size();i++)
		{
			if(v_mini_edge[i].type==1)
			{
				struct position_list p_tmp;
				p_tmp.p_pos_list=NULL;
				struct super_edge_para tmp;
				tmp.edge=v_mini_edge[i].p_mini_edge;
				tmp.graph=p.graph;
				tmp.kmer_len=p.kmer_len;
				tmp.label=0;
				tmp.offset=v_mini_edge[i].offset;
				tmp.plate=p.plate1;
				tmp.p_index_br=p.p_index_br;
				tmp.p_index_uni=p.p_index_uni;

				Locate_Edge_All(tmp,m,&p_tmp);
				pa.push_back(p_tmp);
			}
		}
		mergy_pos_list_multi(pa,&pp);
		for(uint32_t i=0;i<pa.size();i++)
		{
			free(pa[i].p_pos_list);
		}
    }
    r->p_pos_list=pp.p_pos_list;
    r->p_pos_list_len=pp.p_pos_list_len;
    r->offset=pp.offset;
	return 1;
}



