
# SexDetermine.R

## Overview
`SexDetermine.R` is an R script designed for the genetic sex identification of human genomes. The script processes genomic data in BAM or samtools idxstats format and provides sex assignments utilizing several metrics:  Kx, Ky (Yuncu et al., 2024), as well as Rx (Mittnik et al., 2016) and Ry (Skoglund et al., 2013).


## Prerequisites
- R programming environment
- Samtools installed and accessible from your system's PATH
- Bedtools installed and accessible from your system's PATH

## Installation
No additional installation is required beyond R, Samtools, and Bedtools.

## Usage


### Command Line Arguments
- `-ref`: Specify the reference genome (default: 'hg19').
- `-threads`: Set the number of threads for processing (default: 1).
- `-bam`: Specify the BAM file location.
- `-idx`: Specify the samtools idxstats file location.
- `-samtools`: Specify the full path to the samtools executable (default: 'samtools').
- `-intersectBed`: Specify the full path to the bedtools intersectBed executable (default: 'intersectBed').
- `-help, -h`: Display help and exit.

### Running the Script
Run the script from the command line by passing the appropriate arguments. For example:
```bash
Rscript SexDetermine.R -ref hg38 -threads 4 -bam /path/to/file.bam -samtools /usr/local/samtools-1.18/samtools -intersectBed /usr/local/bedtools2/bin/intersectBed
```
Alternatively, using the default parameters:
```bash
Rscript SexDetermine.R -bam /path/to/file.bam
```
Or:
```bash
Rscript SexDetermine.R -idx /path/to/file.idxstats
```

### Important Note

Ensure that the required BED files are present in the same directory as the SexDetermine.R. If these files are missing, the script will stop running and display an error message. You can download the necessary BED files from this repository.

## Citation
Please cite this paper when using SexDetermine.R for your publications.

> Female Lineages and Changing Kinship Patterns in Neolithic Çatalhöyük </br>
> Yüncü, E., Doğu, A. K., Kaptan, D., Kılıç, M. S., Mazzucato, C., Güler, M. N., Eker, E., Katırcıoğlu, B., Chyleński, M., Vural, K. B., Sevkar, A., Atağ, G., Altınışık, N. E., Baloğlu, F. K., Bozkurt, D., Pearson, J., Milella, M., Karamurat, C., Aktürk, Ş., Sağlıcan, E., Yıldız, N., Koptekin, D., Yorulmaz, S., Kazancı, D. D., Aydoğan, A., Karabulut, N. B., Gürün, K., Schotsmans, E. M. J., Anvari, J., Rosenstock, E., Byrnes, J., Biehl, P. F., Orton, D., Lagerholm, V. K., Gemici, H. C., Vasic, M., Marciniak, A., Atakuman, Ç., Erdal, Y. S., Kırdök, E., Pilloud, M., Larsen, C. S., Haddow, S. D., Götherström, A., Knüsel, C. J., Özer, F., Hodder, I., Somel, M. </br>
> https://doi.org/10.1101/2024.06.23.600259


