# Aim: run kGWAS from the start

### Publication of the scritps in Gaurav et al., 2020 (https://www.nature.com/articles/s41587-021-01058-4)
### The original repository is in https://github.com/wheatgenetics/owwc/tree/master/kGWAS

## A) build a k-mer matrix:
1. Create k-mers from each of the samples directly from raw reads. We used 51-mers.
	- script ```run_jellyfish.sh```

Example:
```sh
source jellyfish-2.1.4

sample_id='sample_name'

in_dir=/path/to/raw_reads/${sample_id}
out_dir=/path/out/kmers/${sample_id}
mkdir -p $out_dir

# (zcat $in_dir/*.gz) will assume raw reads per sample per folder
jellyfish count <(zcat $in_dir/*.gz) \
-t 32 \
-C \
-m 51 \
-s 3G \
-o $out_dir/${sample_id}.jf
```

2. Create dump.txt files for each of the samples k-mers created in step 1. Step 1 and two are combined in the same script.
	- script ```run_jellyfish.sh```

Example:
```sh
jellyfish dump -L 2 -ct $out_dir/${sample_id}.jf > $out_dir/${sample_id}.dump.txt
```
Example of a dump.txt file
```
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA     1558
AAAACCCTCCTACAAATCTATATCCTGAGATTAAAGAATACAAAGAAAGGT     4
AAACCAGCTTCTAGTGTTGATGCTGAAATTATAGGTGGATCCACATTAAGT     2
AGACTAGAGACGACAGACTAGAGACTAGAGACCAGAGACGTTAGACTAGAG     4
ACGTATTTAAAAACTATAGTAGGATTACCACTTAATTATTATTGTTTTGTA     3
```

3. Create a configuration file with the name of each sample and the path to the k-mer database created in step 1. Example of configuration file with two samples. The final file must have all the samples to be included in the k-mer matrix. 

Example ```kmers.cfg```
	
	```
	sample_name1	path/to/sample_name1.jf
	sample_name2	path/to/sample_name2.jf
	```

4. For each of the samples create an invividual presence/absence matrix.
	- script: ```run_create_matrix.sh```

This script will output ```presenceMatrix_sample_name1.txt``` file for each of the samples

Example:
 ```sh
 singularity exec ~/tmp/quirozc/python3.img python3 create_matrix.py \
-a sample_name \
-j sample_name.dump.txt \
-c jellies.cfg \
-o presenceMatrix_sample_name.txt
```

Example of the k-mer matrix output:

```
ACCATCCAATTTCGAAAACATCATAACAAGTGAAAAAGAATGATACTCGGG     11111111111111110110101
ACCATCCAATTTCGAAAACATCATAACAAGTGAAAAAGAATGATACTTGGG     11011111111111110110101
ACCATCCAATTTCGAAAACATCATAACAAGTGAAAATGAATGATACTTGGG     00000010000000000000100
ACCATCCAATTTCGAAAACATCATAACAAGTGAAATAGAATGATACTCGGG     00000000000000000000100
```


5. Concatenate each of the individual matrices into one and sort by the k-mer column.

- script ```concatenate_matrix.sh```

```sh
in_dir='/path/to/individual/matrices'
mtx_dir='/path/out'

cat $in_dir/presenceMatrix_sample_* > ${mtx_dir}/presenceMatrix_all_samples_merged.txt

cd ${mtx_dir}
sort \
-k 1 presenceMatrix_all_samples_merged.txt \
-T ${mtx_dir} \
-o presenceMatrix_all_samples_merged_sorted.txt

wait
rm presenceMatrix_samples_merged.txt
```

6. If the k-mer matrix is too big, split the matrix in parts. The final k-mer matrix must be gzip compressed.

```sh
split -n l/5 presenceMatrix_all_samples_merged_sorted -d all_samples_mat_ && gzip all_samples_mat_*
```

7. Create a header file with the name of the samples in a single column in the same order as in the configuration file created in step 3 column 1.

	Example ```accessions.txt```
	```
	sample_name1
	sample_name2
	sample_name...n
	```

8. Sample 100k k-mers from the matrix created in step 6. This will contain the name of the accesion in column one and each of the k-mer in the rest of the columns.

Example:
```sh
singularity exec ~/tmp/quirozc/python3.img python3 ../scripts/random_kmers.py \
-k <(shuf -n 100000 presenceMatrix_all_samples_merged_sorted.txt) \
-a accessions.txt \
-o kmers_100k.tsv
```

The output will be the ```kmers_100k.tsv``` which will be used to compute the PCA and account for population structure when runing kGWAS downstream analysis.

## B) Prepare genome reference:
1. Using the genome reference in FASTA fromat run the script ref_splitter.py. This script will split the sequence in the reference in 10,000 bp sequence blocks, but will keep a single file.

	- script ```split_reference.sh```

Example:
```sh
reference=reference_name
ref_id=${reference}.fa


ref_dir='/path/to/reference'
out_dir='/path/out'
mkdir -p $out_dir

singularity exec ~/tmp/quirozc/python3.img python3 ref_splitter.py \
-a ${ref_dir}/${ref_id} \
-s 10000 \
-o $out_dir/${reference}.fa
```

The output will be:

```
>chr1A:0_10000
CTTAGTTCTACGAATCTTCAGAAATAATCCTGAAGTGGCTCCGCTCTCCTTGAACGCAGGCTCGTAGCTCTTGAGCTGAGATAGTCTTGCGTAGGAAGCCTTGAAAGTTTTGAGTTAACTTGTCCCTTAGGTTTTCCCACGAGTATATTGTACCCGGATG...
>chr1A:10000_20000
CACCTCATGAACCTAGCTAGACCCCCAAAACCTGGCCCATGGAGACGACGCGCCCGCCGCCGCCCCCGCGTCGTCCCCCTCGGCTTCGTCGTTGGTGATGGCGACGTGCTCGCCATGCTCCGTGGTGGGCTCCTCCCCCGGGGCCGGGCCCCGACCTTTT...
etc...
```

2. Create a txt file with the chromosome names in column one and the length of each chromosome in column two. This can be obtained using the original FASTA reference file.

Example:
```sh
## get chr lengths
awk '/^>/ {if (sle){print sle}; print ;sle=0;next; } { sle = sle +length($0)}END{print sle}' \
${ref_dir}/${ref_id} >> $out_dir/chr_sizes_${reference}.txt
```

The output will be ```chromosome_lengths.txt```

```
chr1A	40621098
chr1B	35710944
chr2A	35425885
chr2B	34643735
```

## C) Prepare phenotypes:
1. Create a file with the name of the samples in column one and the scored phenotype average in column two. NOTE: the names of the samples must match the ones used in the header file in step 7 of  A).

Example ```my_phenotype.txt```

```
sample_name1	30.56
sample_name2	41.56
sample_name3	39.89
sample_name4	29.89
sample_namen	35.11

```
2. Create a "Usable file". This will have a simgle column with the name of the samples to use in the GWAS analysis. The names must match the names in the phenotype file and the header file. This file can have fewer samples than phenotype file but it cannot have more (this is a filter file).

Example ```usable.txt```

```
sample_name1
sample_name2
sample_name4
```

## D) run kmerGWAS:
1. Put all the scripts in a single folder and run the "RunAssociation_GLM.py" script. This will use each of the other scrip in the folder.
	- KmerProjection_GLM.py
	- Phenotype_GLM.py
	- RunAssociation_GLM.py

2. If you have the k-mer matrix split in parts, run this script foe each of the parts in parallel.

- script ```01_run_kmerGWAS.sh```

Example using an array:

```sh
i=$SLURM_ARRAY_TASK_ID

reference=tef
base_dir=/base/path/
mtx_dir=${base_dir}/matrix
ref_dir=${base_dir}/genome_parts

declare -a kmer_mtxs=(\
"all_samplesmat_00" \
"all_samplesmat_01" \
"all_samplesmat_02" \
"all_samplesmat_03" \
"all_samplesmat_04" \
)

mtx=${kmer_mtxs[$i]}

phenotype=panicle_len
mc=4

out_dir=${base_dir}/01_results/${phenotype}
mkdir -p $out_dir

singularity exec ~/tmp/quirozc/python3.img python3 ${base_dir}/scripts/RunAssociation_GLM.py \
-i ${base_dir}/matirx/${mtx}.gz \
-hd ${base_dir}/matirx/accessions.txt \
-a $ref_dir/${reference}.fa \
-p ${base_dir}/phenotypes/${phenotype}.txt \
-u ${base_dir}/matirx/usable.txt \
-s ${base_dir}/matirx/kmers_100k.tsv \
-o $out_dir/${reference}_${phenotype}_${mtx}_mc${mc}.txt \
-mc ${mc}
```

Example of the ouput file results:

Column1: chromosome
Column2: start
Column3: end
Column4: association
Column5: correlation of having a k-mer with and a phenotype
Column6: number of k-mer in the window (start-end) with the same association and the same correlation.

```
chr1A	1620000	1630000	6.25	0.32	1
chr1B	5800000	5810000	6.09	0.34	3
chr1B	7320000	7330000	6.12	0.32	1
chr1B	7330000	7340000	6.51	0.34	2
chr1B	7330000	7340000	6.34	0.33	2
```

## E) plot results:
1. If kmerGWAS was run in parts, merge the individual results with the "merge_output.py" script.
```sh
reference=reference_name
phenotype=my_phenotype
mc=4

base_dir=/base/path/
ref_dir=$base_dir/genome_parts
script_dir=$base_dir/scripts
gwas_dir=$base_dir/01_results/${phenotype}

# Merge for filter ploting
singularity exec ~/tmp/quirozc/python3.img python3 ${script_dir}/merge_output.py \
$gwas_dir $gwas_dir/${reference}_${phenotype}_mc${mc}_filtered.txt

```
2. Using the merged file, and the chromosome length created in step 2 of B.

```sh
singularity exec ~/tmp/quirozc/python3.img python3 ${script_dir}/gwas_plot.py \
$gwas_dir/${reference}_${phenotype}_mc${mc}_filtered.txt \
$ref_dir/chromosome_lengths_s.txt \
$gwas_dir/${phenotype}_${reference}_mc${mc}_filtered.jpg
```

DONE
