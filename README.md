# MiniDBG 81.21.0101

Release Date: 8th January, 2021

Author

	~Changyong Yu (Northeastern University in CHINA)
	~Chu Zhao (Northeastern University in CHINA)
	~Jianyu Jin(Northeastern University in CHINA)
1.Introduction

MiniDBG is a tool，it can construct the minimal position list of de bruijn graph (cDBG) of one genome or multiple genomes. The built position lists can be used for locating any edge and path on the graph. MiniDBG also provides the locating methods in the program.

	The k-mers of the genome sequences are the vertex of the graph, the length of the k-mer is an important parameter of the graph. StLiter limits the length of k-mer less than or equal to 32. In our experience, 32 is long enough for various kinds of aims of data analysis when using cDBG. But also, we can modify the codes to support longer k-mer vertex.

	The compressed de bruijn graph is composed of three parts of information, 1) the branching k-mers, 2) the uni-paths and 3) the postion lists of the edges. For large-scale, especially huge-scale genome or multiple genomes, we recommend MiniDBG model for indexing the position of the edges for reducing the cost of memory. 

	There is a shortcoming of MiniDBG. MiniDBG is designed base on the assumption that the alphbet is {A,C,G,T}. Therefore, StLiter can just recognize the four letters. It cannot deal with the other characters with recognizing them. If it encounter with other characters, it will treat them as letter 'A'. This surely will result in false positive results. The further verification on the genome sequence will correct and delete the false positive results. It is worth attention for the users of cDBG for this shortcoming. We think it is not vitally important compared with the key problem of how to filter out true negitive results in data analysis such as read mapping. Therefore, it is accecptable.

	Overall, MiniDBG can effectively construct position list of minimal edges of cDBG for multiple large-scale genomes for locating edges or paths. We believe that the tool will improve the performance of data analysis with the bottleneck of large-scale data.

2.Test Data

	MiniDBG constructs the position lists of the minimal edges of compressed de bruijn graph of genome sequences. The genome sequences are .fa format, and the graph is build by StLiter, a tool for constructing CDBG of genome which can be found using webset: https://github.com/BioLab-cz/MiniDBG. we tested MiniDBG with the genome sequences used in paper "TwoPaCo: An efficient algorithm to build the compressed de Bruijn graph from many complete genomes". we download the data using website: https://github.com/medvedevgroup/TwoPaCo.

3.Building Notes

Our directory contains only a content : 1.codes to build compressed de bruijn graph.

THe directory named #src# consisted of the code necessary to build the MiniDBG.
To build the MiniDBG, change the directory to src and type

	$ make

After that, you can type commands to run tests using MiniDBG. MiniDBG provides 5 methods for building,locating and analyzing the compressed de bruijn graph of large-scale genomes.

4.Usage Notes

1)E method. Generate minimal edges of the cDBG.

	$ ./MiniDBG -M <method> -L <kmer_len> -t <thread_num> 
	
	E method calcultes effectively the minimal edges of the graph. It output one file, take L=16 for example, it will output one file, named "mini016". It is the file saving the minimal edges as (k+1)-mers with k=16. The (k+1)-mers in the files are sorted with the ascending lexicographic order. 

Running test:

	$ ./StLiter -M E -L 16 -t 4

	The above commond takes "kmer016in", "kmer016inad", "kmer016out", "kmer016outad", "kmer016uni", "kmer016uniad","kmer016" and "kmer016ad" binary files as input. All of these files are calculating results of program StLiter. And it will output "mini016" binary file saving the minimal edges with each as an uint64_t integer.

2)P method. Generate position lists of minimal edges of the cDBG.

	$ ./MiniDBG -M <method> -L <kmer_len> -r <ref_path> -m <hashtable_memory> 

	P method calculates position lists of the minimal edges based on minimal edge set. 

Running test:

	$ ./MiniDBG -M P -L 16 -r dataset1-m 10
		
	Input files:
		mini016, dataset1, nomLeu2.fa(reference sequence whose path is saved in dataset1)

	Output files:
		plst016

	The above commond calculates position lists of the minimal edges based on the given minimal edge set with length 16. The memory limit is 10GB.

3)S method. Calculate the counters of the edges and vertex of the graph.

	$ ./MiniDBG -M <method> -L <kmer_len> 

	S method calculates numbers of inbranching edges, out-branching edges, minimal edges. It also output the total number of vertex and edges of the graph. Moreover, it output the distribution of the length of the position list of the minimal edges. These numbers can tell the users the size of the graph and infer the index size of the graph. 

Running test:

	$ ./MiniDBG -M S -L 16
		
	Input files:
		kmer016in, kmer016inad, kmer016out, kmer016outad, kmer016uni, kmer016uniad,kmer016, kmer016ad,mini016 and plst016

	Output files:
		the numbers and the distribution

	The above commonds calculate the numbers of the graph with k-mer length 16.

4)G method. Locate any (k+1)-mer on the reference sequence.

	 $ ./MiniDBG -M <method> -L <kmer_len> -e <edge>

	G method locates any (k+1)-mers on the reference sequence.

Running test:

	$ ./MiniDBG -M G -L 16 -e 625

	Input files:
		kmer016, kmer016uni, mini016 and plst016

	Output files:
		the positions of edge(625 as a (k+1)-mer) on the reference sequence

	The above commonds locate an edge on the reference, the edge can be on the graph or not.  

5)H method. Locate any path (any "ACGT" string) on the reference.

	$ ./MiniDBG -M <method> -L <kmer_len> -h <p_path>

	H method locates any path (any string longer than k+1 and composed of characters in {A,C,G,T}) on the reference sequence. 

Running test:

	$ ./MiniDBG -M H -L 16 -h "ACGTACGTACGTACGTACGT"

	Input files:
		kmer016, kmer016uni, mini016 and plst016

	Output files:
		the positions of the path on the reference sequence

	The above locates the path on the reference.

5.Parameter Settings

The format of a parameter of MiniDBG in the command line is a pair of strings, here we denote the pair as (-p, [q]) or (-p,<q>). String p is the name of the parameter. String q is the value of the parameter input in the command line. [q] represents that the parameter is a optional parameter. <q> represents that the parameter is a necessary parameter.

@parameter (-r,<ref_path>)

	Parameter 'r' gives the path of a text file which saves the filenames of the genome files. For example 'ref_path'="dataset1", than dataset1 is a text file which saves the filename "nomLeu2.fa".

@parameter (-L,<kmer_len>)

	Parameter 'L' set the length of the k-mer, the program will generate the name of file required such as "kmer016" and "mini016" automaticlly.

@parameter (-m, <hashtable_memory>)

	Parameter 'm' decides the maximal memory cost of the hash table used by E or F method in GB. Note that the memory cost of the hash table is not the total memory cost. Actually, the total memory cost of StLiter is hard to be estimated precisely with different genome sequences. In our experience, the parameter can be set as follows.

	Memory available	hashtable_memory	
	>=32GB			5
	>=80GB			10
	>=160GB			15
	other			20

@parameter (-t, <thread_num>)

	Parameter 't' decide the number of threads in the E method.

6.Format of cDBG

	-kmerxxxin
	-kmerxxxout
	-kmerxxx
		~binary files of branching k-mers
		~xxx is the k-mer length
		~each r Bytes as a k-mer,where r=8*([kmer_len/32]+sign(kmer_len/32)).

	-kmerxxxinad
	-kmerxxxoutad
	-kmerxxxad
		~binary files of adjacent information of branching k-mers
		~xxx is the k-mer length
		~the i-th Byte assigns the adjacent information of the i-th k-mer in the corresponding files. 

	-minixxx
		~binary file of bit vector of minimal edges
		~xxx is the k-mer length
	-plstxxx
		~binary file of bit vector of position list
		~xxx is the k-mer length

7.License

	See LICENSE.txt

8.Contacts

	Please e-mail your feedback at cyyneu@126.com.



