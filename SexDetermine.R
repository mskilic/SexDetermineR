args <- commandArgs(TRUE)

print_help <- function() {
    cat("Usage: Rscript SexDetermine.R [options]\n")
    cat("Options:\n")
    cat("  -ref <reference genome>                Specify the reference genome (default: 'hg19')\n")
    cat("  -threads <number>                      Set the number of threads for processing (default: 1)\n")
    cat("  -bam <path to BAM file>                Specify the BAM file location\n")
    cat("  -idx <path to idxstats file>           Specify the samtools idxstats file location\n")
    cat("  -samtools <path to samtools>           Specify the full path to the samtools executable (default: 'samtools')\n")
    cat("  -intersectBed <path to intersectBed>   Specify the full path to the bedtools intersectBed executable (default: 'intersectBed')\n")
    cat("  -help, -h                              Display this help and exit\n")
    cat("\n")
    cat("\nExamples:\n")
    cat("  Using all options:\n")
    cat("    Rscript SexDetermine.R -ref hg38 -threads 4 -bam /path/to/file.bam -samtools /usr/local/samtools-1.18/samtools -intersectBed /usr/local/bedtools2/bin/intersectBed\n")
    cat("  Using only input file with default parameters:\n")
    cat("    Rscript SexDetermine.R -bam /path/to/file.bam\n")
    cat("    Rscript SexDetermine.R -idx /path/to/file.idxstats\n\n") 
    
    cat("Notes:\n")
    cat(" 1) Ensure that 'samtools' and 'intersectBed' are in your system's PATH. If not, you can provide their paths with the -samtools and -intersectBed options.\n")
    cat(" 2) If an idxstats file is provided, the 'Ky' statistic cannot be calculated and 'Kxy Sex' assignment cannot be done. Please provide a BAM file to enable it.\n")
}

# Check for help option at the beginning
if ('-help' %in% args || '-h' %in% args) {
  print_help()
  quit(save = "no", status = 0)
}

valid_options <- c("-ref", "-threads", "-bam", "-idx", "-samtools", "-intersectBed", "-help", "-h")

for (i in seq_along(args)) {
  if (startsWith(args[i], "-") && !(args[i] %in% valid_options)) {
    stop(paste("Invalid option", args[i]))
  }
}


Samtools_Location <- "samtools"
intersectBed_Location <- "intersectBed"
Reference <- "hg19"  # Default reference genome
numthreads <- 1
ApplyQF <- 1

# Function to parse command line arguments
parse_args <- function(args) {
  options <- list(ref = Reference, threads = numthreads, filetype = "", location = "", samtools = Samtools_Location, intersectBed = intersectBed_Location)
  for (i in seq_along(args)) {
    switch(args[i],
           '-ref' = options$ref <- args[i + 1],
           '-threads' = options$threads <- as.numeric(args[i + 1]),
           '-bam' = {options$filetype <- '-bam'; options$location <- args[i + 1]},
           '-idx' = {options$filetype <- '-idx'; options$location <- args[i + 1]},
           '-samtools' = options$samtools <- args[i + 1],
           '-intersectBed' = options$intersectBed <- args[i + 1]
    )
  }
  return(options)
}

# Parsing arguments
options <- parse_args(args)


if (!is.na(as.numeric(options$threads)) && as.numeric(options$threads) %% 1 == 0) {
  numthreads <- as.numeric(options$threads)
} else {
  stop("The value for -threads must be an integer.")
}

# Applying parsed arguments
Reference <- options$ref
InputFileType <- options$filetype
Input_Location <- options$location
Samtools_Location <- options$samtools
intersectBed_Location <- options$intersectBed

script_dir <- getwd()

bed_locations <- list(
  "hg38" = paste(script_dir, "/hg38_Y_95.bed", sep=""),
  "hg19" = paste(script_dir, "/hg19_Y_95.bed", sep=""),
  "GRCh37"  = paste(script_dir, "/GRCh37_Y_95.bed", sep="")
)

# Determine the paths based on the reference genome
if (Reference == "hg38" || Reference == "GRCh38") {
  bed_Location <- bed_locations[["hg38"]]
  
  if (!file.exists(bed_Location)) {
    stop("The required file 'hg38_Y_95.bed' was not found in the same directory as this script. Please ensure the file is present.\n       If it is missing, you can download it from the GitHub repository: https://github.com/mskilic/SexDeterminer.")
  }

  UniqueYLength <- 10139328
} else if (Reference == "hg19" || Reference == "hs37d5") {
  bed_Location <- bed_locations[["hg19"]]
  
  if (!file.exists(bed_Location)) {
    stop("The required file 'hg19_Y_95.bed' was not found in the same directory as this script. Please ensure the file is present.\n       If it is missing, you can download it from the GitHub repository: https://github.com/mskilic/SexDeterminer.")
  }  
  
  UniqueYLength <- 10139645
} else if (Reference == "GRCh37") {
  bed_Location <- bed_locations[["GRCh37"]]
  
  if (!file.exists(bed_Location)) {
    stop("The required file 'GRCh37_Y_95.bed' was not found in the same directory as this script. Please ensure the file is present.\n       If it is missing, you can download it from the GitHub repository: https://github.com/mskilic/SexDeterminer.")
  }  
  
  UniqueYLength <- 10139645
} else {
  stop("Unrecognized reference genome. Please ensure you are using one of the following references: 'hg19', 'hs37d5', 'GRCh37', 'hg38', 'GRCh38'.")
}

# Set up variables for the types of input files
if (InputFileType == "-bam") {
  BAM_Location <- Input_Location
  if (!file.exists(BAM_Location)) {
    stop("BAM file does not exist at the specified path.")
  }
} else if (InputFileType == "-idx") {
  IDX_Location <- Input_Location
  if (!file.exists(IDX_Location)) {
    stop("IDX file does not exist at the specified path.")
  }
} else {
  stop("Unsupported file type specified. Use '-bam' or '-idx'.")
}

CheckSamtools <- Sys.which(Samtools_Location)

if (CheckSamtools == "") {
  if (Samtools_Location != "samtools") {
    stop("samtools executable does not exist at the specified path.")
  } else {
    stop("samtools executable is not accessible from your system's PATH. Please specify the full path to the Samtools executable using -samtools option.")
  }
  
}


CheckIntersectBed <- Sys.which(intersectBed_Location)

if (CheckIntersectBed == "") {
  if (intersectBed_Location != "intersectBed") {
    stop("intersectBed executable does not exist at the specified path.")
  } else {
    stop("intersectBed executable is not accessible from your system's PATH. Please specify the full path to the bedtools intersectBed executable using -intersectBed option.")
  }
  
}


if (InputFileType == "-idx" && file.exists(IDX_Location)){
  
  SortedIDXfile <- paste(IDX_Location,"sorted",sep = ".")
  
  if(Reference == "hg38" || Reference == "GRCh38" || Reference == "GRCh37"){
    
    run_idx=paste("cat",IDX_Location,"| awk '$1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$|^chrY$/' | sort  -k1,1V >",SortedIDXfile,sep=' ')
    system(run_idx) 
  
  } else {
    
    run_idx=paste("cat",IDX_Location,"| awk '$1 ~ /^[1-9]$|^1[0-9]$|^2[0-2]$|^X$|^Y$/' | sort  -k1,1V >",SortedIDXfile,sep=' ')
    system(run_idx) 
    
  }
  idx <- read.table(SortedIDXfile,nrows=24,row.names = 1)
  templist <- strsplit(IDX_Location, split = "/")
  IDX <- sapply(templist, tail, 1)
  PREFIX <- strsplit(IDX, split =".idxstats")
  PREFIX <- as.character(PREFIX)
  
  templist2 <-  strsplit(PREFIX, split = "_")
  templist3 <- lapply(templist2, function(x) x[1])
  templist3 <- as.character(templist3)
  
  Sample <-  strsplit(templist3, split = "-")
  Sample <- lapply(templist3, function(x) x[1])
  Sample <- as.character(Sample)
  
  remover <- paste("rm",SortedIDXfile,sep = ' ')
  system(remover)
  
} else if (InputFileType == "-bam" && file.exists(BAM_Location)){
  templist <- strsplit(BAM_Location, split = "/")
  BAMFILE <- sapply(templist, tail, 1)
  PREFIX <- strsplit(BAMFILE, split =".bam")
  PREFIX <- as.character(PREFIX[[1]])
  
  templist2 <-  strsplit(PREFIX, split = "_")
  templist3 <- lapply(templist2, function(x) x[1])
  templist3 <- as.character(templist3)
  
  Sample <-  strsplit(templist3, split = "-")
  Sample <- lapply(templist3, function(x) x[1])
  Sample <- as.character(Sample)
  
  if (ApplyQF == 0) {
    IDXFILE <- paste(PREFIX[[1]],'idxstats',sep='.')
    if(Reference == "hg38" || Reference == "GRCh38" || Reference == "GRCh37"){
      run_idx=paste(Samtools_Location,"idxstats","-@",numthreads,BAM_Location,"| awk '$1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$|^chrY$/' | sort  -k1,1V >",IDXFILE,sep=' ')
      system(run_idx) 
    } else {
      run_idx=paste(Samtools_Location,"idxstats","-@",numthreads,BAM_Location,"| awk '$1 ~ /^[1-9]$|^1[0-9]$|^2[0-2]$|^X$|^Y$/' | sort  -k1,1V >",IDXFILE,sep=' ')
      system(run_idx)
    }
    idx <- read.table(IDXFILE,nrows=24,row.names = 1)
  } else {
    QSample <- paste(Sample,"Q",sep = '_')
    QBAMFILE <- paste(QSample,"bam",sep = '.')
    QIDXFILE <- paste(QSample,'idxstats',sep='.')
    QualityFilter = paste(Samtools_Location,"view","-@",numthreads,"-q 30 -b",BAM_Location,">",QBAMFILE,sep = ' ')
    system(QualityFilter)
    Indexer = paste(Samtools_Location,"index",QBAMFILE,sep = ' ')
    system(Indexer)
    run_idx=paste(Samtools_Location,"idxstats","-@",numthreads,QBAMFILE,">",QIDXFILE,sep=' ')
    system(run_idx)
    if(Reference == "hg38" || Reference == "GRCh38" || Reference == "GRCh37"){
      run_idx=paste(Samtools_Location,"idxstats","-@",numthreads,QBAMFILE,"| awk '$1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$|^chrY$/' | sort  -k1,1V >",QIDXFILE,sep=' ')
      system(run_idx) 
    } else {
      run_idx=paste(Samtools_Location,"idxstats","-@",numthreads,QBAMFILE,"| awk '$1 ~ /^[1-9]$|^1[0-9]$|^2[0-2]$|^X$|^Y$/' | sort  -k1,1V >",QIDXFILE,sep=' ')
      system(run_idx)
    }
    idx <- read.table(QIDXFILE,nrows=24,row.names = 1)
    remover <- paste("rm",QIDXFILE,sep = ' ')
    system(remover)
    IntersectedReads <- paste(Sample,"intersectedreads",sep = '.')
    run_intersect =  paste(intersectBed_Location,"-abam",QBAMFILE,"-b",bed_Location,"-bed | wc >",IntersectedReads,sep = ' ')
    system(run_intersect)
    remover <- paste("rm",QBAMFILE,sep = ' ')
    system(remover)
    remover <- paste("rm",paste(QBAMFILE,"bai",sep = '.'),sep = ' ')
    system(remover)
    UniqueY <- read.table(IntersectedReads)
    remover <- paste("rm",IntersectedReads,sep = ' ')
    system(remover)
  }
} else {
  stop("Input file not properly specified or does not exist.")
}

TotalRead <- sum(idx$V3)
if (Reference == "hg38" || Reference == "GRCh38" || Reference == "38") {
  idx$UngappedLength <- c(231223641,240863511,198255541,189962376,181358067,170078524,158970135,144768136,122084564,133263006,134634058,133137821,97983128,91660769,85089576,83378703,83481871,80089650,58440758,63944268,40088623,40181019,154893034,26452288)
} else {
  idx$UngappedLength <- c(225934550,238204522,194797140,188042934,177695260,167395067,155536559,142964911,120626573,131314747,131169619,130481395,95589878,88289540,81694769,78884753,78129607,74661510,56060841,59505520,35134224,34894566,151100560,25653566)
}
idx$ReadsOverUngappedLength <- idx$V3/idx$UngappedLength
idx$ReadsOverGappedLength <- idx$V3/idx$V2 # for Rx
temp1  <- idx[-24,]
temp1  <- temp1[-23,]
Xreads <- idx[23,"V3"]
Yreads <- idx[24,"V3"]

if (Xreads == 0){
  stop("Zero Reads on X chromosome")
}

ratio <- temp1$ReadsOverUngappedLength
sd <- sqrt(sum((ratio - mean(ratio)) ^ 2 / (length(ratio) - 1)))
SE <- sd/sqrt(22)
mean <- mean(ratio)

LCI <- mean - 2.08*SE
UCI <- mean + 2.08*SE
LowExpectedXX <- idx[23,"UngappedLength"]*LCI
ExpectedXX <- idx[23,"UngappedLength"]*mean
UpExpectedXX <- idx[23,"UngappedLength"]*UCI
LKx <- Xreads/UpExpectedXX
Kx  <- Xreads/ExpectedXX
UKx <- Xreads/LowExpectedXX

Ky <- "NA"

if (LKx > 0.8) {
  Kxsex = "XX"
} else if (UKx < 0.6 && Kx != 0) {
  Kxsex = "XY"
} else if (LKx > 0.6 && UKx > 0.8) {
  Kxsex = "XX?"
} else if (LKx < 0.6 && UKx < 0.8) {
  Kxsex = "XY?"
} else {
  Kxsex = "NA"
}

# Rx 
Xratio <- idx[23,"ReadsOverGappedLength"]
temp1$Rx <- Xratio/temp1$ReadsOverGappedLength

Rxs <- temp1$Rx
sdrx <- sqrt(sum((Rxs - mean(Rxs)) ^ 2 / (length(Rxs) - 1)))
SErx <- sdrx/sqrt(22) 
Rx <- mean(Rxs)
LRx <- Rx - 1.96*SErx
URx <- Rx + 1.96*SErx

if (Rx =="Inf" || Rx =="NaN"){
  Rxsex = "NA"
} else {
  if (LRx > 0.8) {
    Rxsex = "XX"
  } else if (URx < 0.6) {
    Rxsex = "XY"
  } else if (LRx > 0.6 && URx > 0.8) {
    Rxsex = "XX?"
  } else if (LRx < 0.6 && URx < 0.8) {
    Rxsex = "XY?"
  } else {
    Rxsex = "NA"
  }
}


#Ry
n=Yreads+Xreads
Ry=Yreads/n
SEry=sqrt((Ry*(1.0-Ry))/n)
LRy <- Ry - 1.96*SEry
URy <- Ry + 1.96*SEry

if(Ry =="NaN"){
  Rysex = "NA"
}else {
  if (Ry == 0) {
    Rysex = "XX??"
  } else if (URy < 0.016) {
    Rysex = "XX"
  } else if (LRy > 0.075) {
    Rysex = "XY"
  } else if (LRy > 0.016 && URy > 0.075) {
    Rysex = "XY?"
  } else if (LRy < 0.016 && URy < 0.075) {
    Rysex = "XX?"
  }else {
    Rysex = "NA"
  }
}


# Kxy
if(ApplyQF != 0 && InputFileType == "-bam") {
  
  UniqueYreads <- UniqueY$V1
  LowExpectedXY <- UniqueYLength*LCI/2
  ExpectedXY <- UniqueYLength*mean/2
  UpExpectedXY <- UniqueYLength*UCI/2
  LKy <- UniqueYreads/UpExpectedXY
  Ky  <- UniqueYreads/ExpectedXY
  UKy <- UniqueYreads/LowExpectedXY
  
  #Dxy
  
  if(Kx < 0.5) {
    DKx <- 0.5
  } else if(Kx > 1) {
    DKx <- 1
  } else {
    DKx <- Kx
  } 
  
  if(Ky > 1) {
    DKy <- 1
  } else {
    DKy <- Ky
  }
  
  DistanceToXX <- sqrt((1 - DKx)^2 + (DKy)^2)
  DistanceToXY <- sqrt((0.5 - DKx)^2 + (1 - DKy)^2)
  DistanceToXXY <- sqrt((1 - DKx)^2 + (1 -DKy)^2)
  DistanceToX0 <- sqrt((0.5 - DKx)^2 + (DKy)^2)
  
  if (LKx > 0.8 && UKy < 0.15  ) {
    Kxysex = "XX"
    Distance <- DistanceToXX
  } else if (LKx > 0.6 && UKx > 0.8 && UKy < 0.15 ) {
    Kxysex = "XX?"
    Distance <- DistanceToXX
  } else if (LKx > 0.8 && LKy < 0.15 && UKy < 0.70  ) {
    Kxysex = "XX?"
    Distance <- DistanceToXX
  } else if (LKx > 0.6 && UKx > 0.8 && LKy < 0.15 && UKy < 0.70  ) {
    Kxysex = "XX??"
    Distance <- DistanceToXX
  } else if (UKx < 0.6 && LKy > 0.70 ) {
    Kxysex = "XY"
    Distance <- DistanceToXY
  } else if (LKx < 0.6 && UKx < 0.8 && LKy > 0.70 ) {
    Kxysex = "XY?"
    Distance <- DistanceToXY
  } else if (UKx < 0.6 && LKy > 0.15 && UKy > 0.70 ) {
    Kxysex = "XY?"
    Distance <- DistanceToXY
  } else if (LKx < 0.6 && UKx < 0.8 && LKy > 0.15 && UKy > 0.70 ) {
    Kxysex = "XY??"
    Distance <- DistanceToXY
  } else if (LKx > 0.8 && LKy > 0.70 && UniqueYreads > 50 ) {
    Kxysex = "XXY"
    Distance <- DistanceToXXY
  } else if (LKx > 0.6 && UKx > 0.8 && UKy > 0.70  && UniqueYreads > 50 ) {
    Kxysex = "XXY?"
    Distance <- DistanceToXXY
  } else if (LKx > 0.8 && LKy > 0.10 && UKy > 0.70  && UniqueYreads > 50 ) {
    Kxysex = "XXY?"
    Distance <- DistanceToXXY
  } else if (LKx > 0.6 && UKx > 0.8 && LKy > 0.15 && UKy > 0.70 && UniqueYreads > 50 ) {
    Kxysex = "XXY??"
    Distance <- DistanceToXXY
  } else if (UKx < 0.6 && UKy < 0.15  && Xreads > 125 ) {
    Kxysex = "X0"
    Distance <- DistanceToX0
  } else if (LKx < 0.6 && UKx < 0.8 && UKy < 0.15 && Xreads > 125 ) {
    Kxysex = "X0?"
    Distance <- DistanceToX0
  } else if (UKx < 0.6 && LKy < 0.15 && UKy < 0.70 && Xreads > 125 ) {
    Kxysex = "X0?"
    Distance <- DistanceToX0
  } else if (LKx < 0.6 && UKx < 0.8 && LKy < 0.15 && UKy < 0.70 && Xreads > 125) {
    Kxysex = "X0??"
    Distance <- DistanceToX0
  }else {
    Kxysex = "NA"
    Distance <- NA
  } 
} else {
  Kxysex = "NA"
  Distance <- NA
}


cat("Sample: ",Sample," Total Reads: ",TotalRead," Xreads: ",Xreads," Yreads: ",Yreads," Kx: ",Kx," Kx 95% CI: ",LKx,"-",UKx," Kx Sex: ",Kxsex," Ky: ",Ky," Kxy Sex: ",Kxysex," Distance: ",Distance," Rx: ",Rx," Rx Sex: ",Rxsex," Ry: ",Ry," Ry Sex: ",Rysex,"\n")

#Kx-Ky Plot
if (ApplyQF != 0 && InputFileType == "-bam"){
  cat("------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n")
  
  cat("\t <<< Kx-Ky Plot >>> | Ky: Y Axis // Kx: X Axis | \n")
  
  Kxpos <- floor((Kx-0.3)*50)
  if ( Kxpos > 15) {
    Kxpos <- Kxpos -1
  }
  Kypos <- 21 - floor(Ky*17)
  
  if(Kypos > 20) {
    Kypos <- 20
  }
  if(Kypos < 1) {
    Kypos <- 1
  }
  if(Kxpos > 50) {
    Kypos <- 50
  }
  if(Kxpos < 1) {
    Kxpos <- 1
  }
  
  for (i in 1:20) {
    cat(0.060*(21-i))
    for(j in 1:50){
      if(j == 1 && i ==1){
        cat(" |")
      }else if(j == 1 && i ==11){
        cat(" |")
      }else if(j == 1){
        cat("|")
      } 
      if(i == 5 && j == 11) {
        cat("M")
      }
      if(i == 20 && j == 35) {
        cat("F")
      }
      if(i == 20 && j == 11) {
        cat("T")
      }
      if(i == 5 && j == 35) {
        cat("K")
      }
      if(i == Kypos && j == Kxpos){
        cat("o")
      } else {
        cat(" ")
      }
      
    }
    cat("\n")
  }
  cat("    -----------------------------------------------------------------")
  cat ("\n")
  for(k in 1:55){
    if(k %% 5 == 0 ){
      cat("|")
    } else {
      cat(" ")
    }
  }
  cat("\n")
  cat("  ","0.3","","0.4","","0.5","","0.6","","0.7","","0.8","","0.9"," ","1"," ","1.1","","1.2","","1.3")
  
  cat ("\n\n o for Sample // M for Male // F for Female // K for Klinefelter // T for Turner \n")
}
