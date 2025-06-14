
![Spectre](./logo.png)
# Spectre - Long read CNV caller
Spectre is a long read copy number variation (CNV) caller. 
Spectre is designed to detect large CNVs (>100kb) in a couple of minutes depending on your hardware.

To calculate CNVs Spectre uses primarily the coverage (Read depth) data. 
However, it can also use SNV data to detect loss of heterozygosity (LoH) regions.
Additionally, Spectre can use the breakpoint (SNF) data from Sniffles to improve the CNV calling. However, it has to be converted to the SNFJ format using [snf2json](https://github.com/philippesanio/snf2json).


The CNV output of Spectre is stored in three files, VCF, BED and .SPC which can be used in the population mode.

Furthermore, Spectre offers a population mode, which can be used to search for CNV support in multiple samples. 
Compared to other tools, Spectre searches not only in the final CNVs but also in CNV candidates which did not qualify for the final output of Spectre.
## Install Spectre
Tested on **Python** version: **3.11** and **3.10**

Install Spectre with Pip:
```bash
pip install spectre-cnv
```
Get the latest changes by building Spectre from source on your own and install it locally in your conda environment.
```bash
pip install build
git clone https://github.com/fritzsedlazeck/Spectre.git
cd ./Spectre
python3 -m build
pip install dist/spectre_cnv-<VERSION>.tar.gz # replace <VERSION> with e.g. 0.2.0
```

Setup a conda environment for Spectre (copy and paste the following commands)
```bash
conda create -n spectre python=3.10 pysam==0.22.0 numpy==1.24.3 pandas==2.0.1 matplotlib==3.7.1 scipy==1.10.1 -y
conda activate spectre
```
Alternatively, you can use pip for installing the packages stored in the requirements txt

```bash
conda create -n spectre python=3.10 pip -y
conda activate spectre
pip install -r requirements.txt
```
or install everything manually (check for package version in the requirements.txt file)

|Program| Conda                                       |
|-------|---------------------------------------------|
| python3 | conda install python=3.10                   |
| pysam | conda install -c bioconda pysam=0.22.0      |
| pandas| conda install -c anaconda pandas==2.0.1     |
| numpy| conda install -c anaconda numpy==1.24.3     |
| scipy| conda install -c anaconda scipy==1.10.1     |
| matplotlib| conda install -c anaconda matplotlib==3.7.1 |


## How to run
Spectre need as input:

Prerequisites:
Extract the coverage data from a BAM using [Mosdepth](https://github.com/brentp/mosdepth).
Example command:
```bash
mosdepth -t 8 -x -b 1000 -Q 20 "${out_path}/${sample_id}" "${bam_path}"
```

>IMPORTANT: We recommend to run **Mosdepth** with a **bin size of 1kb** and a **mapping quality of at least 20** (-Q 20), as Spectre is optimized for that.

Coverage produced at this resolution forms the basis of Spectre's ploidy
calculation. Each 1 kb window is normalised by the genome wide median
coverage to yield an estimated ploidy value. When creating the genome wide
plot these values are smoothed by averaging in a window of approximately one
megabase derived from the median distance between consecutive coverage
entries. The resulting smoothed ploidy values are then
coloured for plotting using a five step palette ranging from ``#b2182b`` at
zero through ``#fbea6cf9`` at two to ``#313695`` at four.

- The region coverage file (mosdepth)
- SampleID e.g.
- Output directory 
- Reference genome (can be bgzip compressed)

Optional
- **MDR** file (if not already generated, Spectre will do that for you. You can also use the MDR file for every sample which has been aligned to the same reference genome)
- VCF file containing SNV
- SNF data from Sniffles (if parsed through [snf2json](https://github.com/philippesanio/snf2json))


## Run Spectre
### MDR file
MDR files hold information about N regions in the reference genome and restrict Spectre's use of data from those regions. 
We are providing sample MDR files for the reference genomes GRCh37 and GRCh38.

Spectre will generate an MDR file for you if not provided, which can take some time. 
Thus, we highly recommend generating an MDR file for your reference genome before running Spectre on multiple samples aligned to the same reference.


Providing an MDR file will save you a substantial amount of time, as Spectre will not have to calculate the N regions for every sample.

The generation of the MDR file can be done with either the `RemoveNs` or `CNVCaller` command. In the latter case, the MDR (metadata. MDR) file will be saved in the sample's output directory.
Running the `RemoveNs` will generate an MDR file under the name in the directory you specified.
```bash
spectre RemoveNs \
  --reference reference.fasta.gz \
  --output-dir output_directory_path/ \
  --output-file name_of_metadata_file.mdr
  
```
### Blacklists
The blacklist is a supplementary file to the MDR file. It contains regions which should be ignored by Spectre.
Those regions are based on gap data from USCS. 
During testing we found that the gap data is not totally sufficient masking high frequency coverage regions such as telomeric and centromeric regions.
Thus we have extended the especially those problematic regions in the blacklist file. (grch37_blacklist_spectre_refined.bed and grch38_blacklist_spectre.bed)

### Run Spectre with a single sample
```bash
spectre CNVCaller \
  --coverage mosdepth/sampleid/mosdepth.regions.bed.gz \
  --sample-id sampleid \
  --output-dir sampleid_output_directory_path/ \
  --reference reference.fasta.gz
```
### Run Spectre with multiple samples
Run Spectre with multiple samples:
>INFO: This will start the population mode automatically. All provided settings will be applied to all samples.

>NOTE: If population flag is not set, Spectre will run in single sample mode. Thus, calculating only CNVs for the first sample. 

```bash
spectre.py CNVCaller \
  --coverage mosdepth/sampleid-1/mosdepth.regions.bed.gz mosdepth/sampleid-2/mosdepth.regions.bed.gz \
  --sample-id sampleid-1 sampleid-2 \
  --output-dir sampleid_output_directory_path/ \
  --reference reference.fasta.gz
  --population
```

### Population mode
Run Spectre in population mode with two or more samples:
>INFO: Spectre produces an intermediate file (.spc) which contains all calculated CNVs from a given samples. They are 
> located in the output folder of given sample.

```bash
spectre population \
  --candidates /path/to/sample1.spc /path/to/sample2.spc \
  --sample-id output_name \
  --output-dir sampleid_output_directory_path/
```


### Help
```
vcf_utils <command> [<args>]
    Spectre:
        CNVCaller:
            [Required]
                --coverage     Path to the coverage file from Mosdepth output. Expects the following files:
                                   <prefix>.regions.bed.gz
                                   <prefix>.regions.bed.gz.csi
                               Can be one or more paths. However, providing multiple samples is only intended to
                               work with the --population flag. Example:
                                    --coverage /path/md1.regions.gz /path/md2.regions.gz
                --sample-id    Sample name/ID. Can be one or more ID. However, providing multiple sample ids is only
                               intended to work with the --population flag. Example:
                                    --sample-id id1 id2
                --output-dir   Output directory
                --reference    Reference sequence used for mapping (for N removal)
            [Optional, if missing it will be created]
                --metadata     Metadata file for Ns removal (this will speed up Spectre massively if provided)
                --n-size       Required amount of consecutive Ns to be considered an NRegion 
                               in the reference sequence (Default = 5)

            [Optional]
                --blacklist    Blacklist in bed format for sites that will be ignored (Default = "")
                --only-chr     Comma separated list of chromosomes to use (e.g. chr1,chr2,chr3)
                --ploidy       Set the ploidy for the analysis, useful for sex chromosomes (Default = 2)
                --ploidy-chr   Comma separated list of key:value-pairs for individual chromosome ploidy control
                               (e.g. chrX:2,chrY:1) If chromosome is not specified, the default ploidy will be used.
                --snfj         Breakpoints from from Sniffle which has been converted from the SNF to the SNFJ format.
                               SNFJ files can be generated using the program snf2json.
                --min-cnv-len  Minimum length of CNV (Default 100kb)
                --snv          VCF file containing the SNV for the same sample CNV want to be called
                --cancer       Set this flag if the sample is cancer (Default = False) This will disable some safety 
                               checks, when determining the DEL and DUP thresholds. 

            [Optional, Coverage]
                --sample-coverage-overwrite     Overwrites the calculated sample coverage, which is used to normalize
                                                the coverage. e.g. a value of 30 equals to 30X coverage.
                --disable-max-coverage          Disables the maximum coverage check. This will allow to call CNVs

            [Optional, LoH (requires --snv)]
                --loh-min-snv-perkb             Minimum number of SNVs per kilobase for an LoH region (default=5)
                --loh-min-snv-total             Minimum number of SNVs total for an LoH region (default=100)
                --loh-min-region-size           Minimum size of a region for a LoH region (default=100000)

                --population   Runs the population mode on all provided samples. It will apply all the provided
                               configurations as well as the default population mode values to all samples.
                --threads      Amount of threads (This will boost performance if multiple samples are provided)
        RemoveNs:
            [Required]
                --reference    Reference genome used for mapping
                --output-dir   Output dir
                --output-file  Output file for results
            [Optional]
                --n-size       Required amount of consecutive Ns to be considered an NRegion 
                               in the reference sequence (Default = 5)
                --save-only    Will only save the metadata file and not show the results on screen (Default = False)

        Population:
            [Required]
                --candidates   At least 2 .spc sample files which should be used in the population mode.
                               (e.g. sample1.spc sample2.spc)
                --sample-id    The name of the sample-id will be added accordingly at the output.
                               (e.g. population_mode_<sample-id>.vcf.gz)
                --output-dir   Path of the output directory
            [Optional]
                --reference    Reference sequence
                --reciprocal-overlap        Minimum reciprocal overlap for supporting CNVs [0.0 - 1.0] (Default = 0.8)
                --disable-quality-filter    Disables the quality filter for the population mode. Spectre will also
                                            search for supporting CNVs in the .SPC files, which have not been reported
                                            as final CNVs in the VCF and BED file.
        Version:
            version    Shows current version/build
```

## Genome CNV plot

Spectre provides a genome wide plot summarising coverage and CNV calls. The
coverage input stems from Mosdepth run with a 1 kb bin size. After normalising
each 1 kb window by the genome wide median coverage Spectre obtains ploidy
estimates along every chromosome. For visualisation these estimates are
smoothed by averaging in windows of roughly one megabase derived from the
median spacing of the coverage data. 

