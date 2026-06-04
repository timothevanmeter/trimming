## FAST TRIMMING


---

Fast trimming is a script to trim NGS results to find a target portion of the reads. Additionally, there is multiple filters to insure that the targeted sequence to retrieve has:
1. the expected length
2. an overall quality (Qscore) above a user-defined threshold
3. contains no stop codons

---


### Installation and requirements


**WARNING**: The script expects a Linux OS to run!\
For example, the script relies heavily on **`bash`** and **`awk`**.


In order to run the script you need the following requirements to be installed:
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [vsearch](https://github.com/torognes/vsearch)
- [zlib](https://www.zlib.net/)\
this is only needed if you want to run the script directly on compressed files .fastq.gz
- [jq](https://jqlang.org/)\
this is only needed to produce a `.json` file saving the trimming metadata for each processed file. However, I *highly* recommend to keep it.



#### Installing all the requirements


```
#CUTADAPT
sudo apt install cutadapt

#SEQKIT
wget https://github.com/shenwei356/seqkit/releases/download/v2.13.0/seqkit_darwin_amd64.tar.gz
tar -xvzf seqkit_darwin_amd64.tar.gz
sudo cp seqkit /usr/local/bin/

#VSEARCH
wget https://github.com/torognes/vsearch/archive/v2.31.0.tar.gz
tar -xvzf v2.31.0.tar.gz
cd vsearch-2.30.6
./autogen.sh
./configure CFLAGS="-O2" CXXFLAGS="-O2"
make ARFLAGS="cr"
sudo make install

#ZLIB
wget https://zlib.net/current/zlib.tar.gz
tar -xvzf zlib.tar.gz
cd zlib
./configure
make
sudo make install

#JQ
sudo apt install jq
```

---


### Usage


This script processes a `.fastq` or `.fastq.gz` file and ouptus a `.counts` file listing, in decreasing order, the unique sequence variants and their counts.\
It can be called from the command line as follows:


Mandatory options:

	```./fast_trimming.sh -a1 <adapter_1> -a2 <adapter_2> -is <insert_size> -i <file_to_process.fastq>```

example:

	```./fast_trimming.sh -a1 CAACCAGAGC -a2 TCCACCACTG -is 21 -i B2687-1.fastq```


Options:

`-i --input`	fastq file to process.\
	`-h --help`	Show this information.\
	`-k --keep`	Keep all intermediary files generated during runtime.\
	`-t --threshold`	Define the expected error threshold above which sequences are discarded. By default it is set to 0.01\
	`-a1 --adapter1`	Define the 5` sequence around the sequence of interest \
	`-a2 --adapter2`	Define the 3` sequence around the sequence of interest \
	`-is --insertSize`       Define the size of the insert or sequence of interest that you expect: only sequences of this exact size will be kept\
	`-jn --jobNumber`	       Define the job number for this process; only for when using parallel processing\
	`-z --summary`	       Produce a summary file with the number of sequences kept at each filtering step\
	`-s --silent`	       Run silently without any messages



#### Parallel processing
Parallel processing can be done by using the bash command [parallel](https://www.gnu.org/software/parallel/sphinx.html)\

This command expects that all the files to be processed are in the current folder and all end with the `.fastq.gz` extension.\
```
parallel --eta -j 4 fast_trimming.sh -z -t 0.1 -is 21 -a1 <adapter1> -a2 <adapter2> -i {} -jn {#} ::: *.fastq.gz
```
The number of concurrent job numbers can be increased (here 4 as indicated by `-j 4`), but remember that cutadapt is also distributing its load across CPUs which is controlled by the `--cores=<number_of_cores_used>` option.