# Cancer Genomics Data Analysis Exercise 1
This pipeline presumes you are running an unix system and have installed python 3 with seaborn, numpy and matplotlib as well as BWA for alignment and Samtools for read depth analysis.

## Step 1 - fetching the data
To download the data we can simply use wget as follows:
```sh
wget https://gear.embl.de/data/.exercise/tu.r1.fq.gz
wget https://gear.embl.de/data/.exercise/tu.r2.fq.gz
wget https://gear.embl.de/data/.exercise/wt.r1.fq.gz
wget https://gear.embl.de/data/.exercise/wt.r2.fq.gz
```

Besides the data we also need to fetch the reference data:
```sh
wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz
```
While there are more recent human genomes assembled, we shall use grch37 as the data is restricted using coordinates relevant to this assembly `chrX:20000000-40000000`.

## Step 2 - alignment
Before aligning to our reference genome, we need to index it as follows:
```sh
bwa index Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz
```

The alignment itself is done using `bwa mem`, which returns sam format on its stdout. Thus we pipe it into samtools for BAM conversion and sorting. Both the tumor and the wildtype reads are processed independently. The `-t 8` argument of bwa can be changed to fit the number of threads you wish to dedicate to the alignment.
```sh
bwa mem -t 8 Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz tu.r1.fq.gz tu.r2.fq.gz | samtools view -O BAM - | samtools sort -T tu_temp -O bam -o tu.sorted.bam -
bwa mem -t 8 Homo_sapiens.GRCh37.dna.chromosome.X.fa.gz wt.r1.fq.gz wt.r2.fq.gz | samtools view -O BAM - | samtools sort -T wt_temp -O bam -o wt.sorted.bam -
```

## Step 3 - read depth
As stated above, our data is not covering the whole `chrX`, so we will first trim it. To do so, we also need to index our sorted data.
```sh
samtools index tu.sorted.bam
samtools index wt.sorted.bam
samtools view -b tu.sorted.bam X:20000000-40000000 > tu.sorted.trimmed.bam
samtools view -b wt.sorted.bam X:20000000-40000000 > wt.sorted.trimmed.bam
```

Then `samtools` can be used to get the read depth:
```sh
samtools index tu.sorted.trimmed.bam
samtools index wt.sorted.trimmed.bam
samtools depth -a -r X:20000000-40000000 tu.sorted.trimmed.bam > tu.depth
samtools depth -a -r X:20000000-40000000 wt.sorted.trimmed.bam > wt.depth
```

## Step 4 - plotting
In the final step we use a python script to generate plots of the read depth.
```sh
python3 readDepthPlot.py --chromosome=X --fromPos=20000000 --toPos=40000000 --window=10000 tu.depth wt.depth
```

## Results
The results are then found in three png files:
- `wt.depth.png` - the read depth plot of the wild type
- `tu.depth.png` - the read depth plot of the tumor
- `ratio.png` - the log2 ratio of tu/wt
