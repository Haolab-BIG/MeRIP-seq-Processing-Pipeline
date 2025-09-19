#!/bin/bash

# MeRIP-seq Complete Analysis Pipeline
# 包含从fastq到peak calling和motif分析的完整流程

print_help() {
  echo "Usage: $0 --IpFq <treatment_fastq> --InputFq <control_fastq> --outputdir <output_directory> --referencedir <reference_index_directory> --gtf <gtf_file> --sif <singularity_environment> [options]"
  echo "--IpFq          : Path to the treatment fastq (required)"
  echo "--InputFq       : (optional) Path to the input control fastq"
  echo "--outputdir     : Path to the output directory (required)"
  echo "--referencedir  : Path to the directory where bowtie reference build with prefix (required)"
  echo "--gtf           : Path to the GTF annotation file (required for peak calling)"
  echo "--sif           : Path to the singularity environment file (required)"
  echo "--threads       : (optional) Number of threads to use, default 8"
  echo "--binSize       : (optional) Number of binsize to use, default 10"
  echo "--genome        : (optional) Genome version (hs for human, mm for mouse), default hs"
  echo "--skip_peak     : (optional) Skip peak calling step"
  echo "--skip_motif    : (optional) Skip motif analysis step"
}

if [[ "$#" -eq 0 || "$1" == "--help" ]]; then
  print_help
  exit 0
fi

# 参数初始化
IpFq=""
InputFq=""
outputdir=""
referencedir=""
gtf=""
sif=""
threads=8
binSize=10
genome="hs"
skip_peak=false
skip_motif=false

# 解析参数
while [[ $# -gt 0 ]]; do
  case $1 in
    --IpFq) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: fastq file $2 not found!"
        exit 1
      fi
      IpFq="$2"; shift 2 ;;
    --InputFq)
      if [[ -z "$2" ]]; then
        InputFq=""
        shift 1
      else
        if [[ ! -f "$2" ]]; then
          echo "Error: fastq file $2 not found!"
          exit 1
        fi
        InputFq="$2"
        shift 2
      fi;;
    --outputdir) 
      if [[ -z "$2" ]]; then
        echo "Error: output directory $2 not found!"
        exit 1
      fi
      outputdir="$2"; shift 2 ;;
    --referencedir) 
      if [[ -z "$2" ]]; then
        echo "Error: reference directory $2 not found!"
        exit 1
      fi
      referencedir="$2"; shift 2 ;;
    --gtf) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: GTF file $2 not found!"
        exit 1
      fi
      gtf="$2"; shift 2 ;;
    --sif) 
      if [[ -z "$2" || ! -f "$2" ]]; then
        echo "Error: singularity sif file $2 not found!"
        exit 1
      fi
      sif="$2"; shift 2 ;;
    --threads) 
      threads="$2"; shift 2 ;;
    --binSize) 
      binSize="$2"; shift 2 ;;
    --genome) 
      genome="$2"; shift 2 ;;
    --skip_peak)
      skip_peak=true; shift 1 ;;
    --skip_motif)
      skip_motif=true; shift 1 ;;
    *) 
      echo "Unknown option $1"; print_help; exit 1 ;;
  esac
done

# 检查必需参数
if [[ -z "$IpFq" || -z "$outputdir" || -z "$referencedir" || -z "$sif" ]]; then
  echo "Error: Missing required parameters"
  print_help
  exit 1
fi

# 如果进行peak calling，检查GTF文件
if [[ "$skip_peak" = false && -z "$gtf" ]]; then
  echo "Error: GTF file is required for peak calling. Use --skip_peak to skip peak calling or provide --gtf"
  exit 1
fi

########### 设置目录结构
mkdir -p "${outputdir}"
mkdir -p "${outputdir}/qc"
FastQCdir="${outputdir}/qc"
mkdir -p "${outputdir}/trim"
trimedfadir="${outputdir}/trim"
mkdir -p "${outputdir}/bam"
alignoutdir="${outputdir}/bam"
mkdir -p "${outputdir}/bw"
bwoutdir="${outputdir}/bw"
mkdir -p "${outputdir}/peak"
pcoutdir="${outputdir}/peak"
mkdir -p "${outputdir}/multiqc"
multiqcdir="${outputdir}/multiqc"
mkdir -p "${outputdir}/logs"

########### 设置Singularity命令
[[ -n "$IpFq" ]] && IpFq=$(readlink -f "$IpFq")
[[ -n "$InputFq" ]] && InputFq=$(readlink -f "$InputFq")
[[ -n "$outputdir" ]] && outputdir=$(readlink -f "$outputdir")
[[ -n "$referencedir" ]] && referencedir=$(readlink -f "$referencedir")
[[ -n "$gtf" ]] && gtf=$(readlink -f "$gtf")
[[ -n "$sif" ]] && sif=$(readlink -f "$sif")

bind_dirs=()
[[ -n "$IpFq" ]] && bind_dirs+=("$(dirname "$IpFq")")
[[ -n "$InputFq" ]] && bind_dirs+=("$(dirname "$InputFq")")
[[ -n "$outputdir" ]] && bind_dirs+=("$outputdir")
[[ -n "$referencedir" ]] && bind_dirs+=("$(dirname "$referencedir")")
[[ -n "$gtf" ]] && bind_dirs+=("$(dirname "$gtf")")

# 去重绑定目录
declare -A seen
bind_dirs_unique=()
for dir in "${bind_dirs[@]}"; do
    real_dir=$(readlink -f "$dir")
    if [[ -n "$real_dir" && -z "${seen[$real_dir]}" ]]; then
        bind_dirs_unique+=("$real_dir")
        seen["$real_dir"]=1
    fi
done

SING_EXEC="singularity exec --cleanenv"
for dir in "${bind_dirs_unique[@]}"; do
    SING_EXEC+=" -B $dir:$dir"
done
SING_EXEC+=" $sif"

########### gzip压缩（如果需要）
if [[ "$IpFq" != *.gz ]]; then
    echo "gzip... ${IpFq}"
    $SING_EXEC gzip -f "$IpFq"
    IpFq="${IpFq}.gz"
fi
if [[ -n "$InputFq" && "$InputFq" != *.gz ]]; then
    echo "gzip... ${InputFq}"
    $SING_EXEC gzip -f "$InputFq"
    InputFq="${InputFq}.gz"
fi

########### 处理fastq文件
files=()
[[ -n "$IpFq" ]] && files+=("$IpFq")
[[ -n "$InputFq" ]] && files+=("$InputFq")

for fqpath in "${files[@]}"; do
    fileName=$(basename "$fqpath" | sed -E 's/\.f(ast)?q\.gz$//')
    echo "Processing file: $fileName"
    
    # FastQC
    echo "Running FastQC..."
    $SING_EXEC fastqc -o "${FastQCdir}" -t ${threads} -q "$fqpath"
    
    # Trim Galore
    echo "Running Trim Galore..."
    $SING_EXEC trim_galore -q 25 --cores ${threads} --phred33 --length 36 -e 0.1 \
        --stringency 3 --fastqc --illumina --basename "trim_${fileName}" -o "${trimedfadir}" "$fqpath"
    
    # STAR比对
    echo "Running STAR alignment..."
    $SING_EXEC STAR --runThreadN ${threads} \
         --genomeDir "${referencedir}" \
         --readFilesIn "${trimedfadir}/trim_${fileName}_trimmed.fq.gz" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${alignoutdir}/${fileName}." \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 1 \
         --outFilterMismatchNmax 3 \
         --outFilterMismatchNoverLmax 0.1 \
         --outSAMattributes All >"${alignoutdir}/${fileName}.STAR.log" 2>&1
    
    # 去除重复reads
    echo "Removing duplicate reads..."
    $SING_EXEC samtools sort -n -@ ${threads} "${alignoutdir}/${fileName}.Aligned.sortedByCoord.out.bam" | \
        $SING_EXEC samtools view -q 255 -b | \
        $SING_EXEC samtools fixmate -m -r - - | \
        $SING_EXEC samtools sort -@ ${threads} - | \
        $SING_EXEC samtools markdup -r -s - "${alignoutdir}/${fileName}.DeDup.bam" 2> "${alignoutdir}/${fileName}.markdup.log"
    
    $SING_EXEC samtools flagstat "${alignoutdir}/${fileName}.DeDup.bam" > "${alignoutdir}/${fileName}.flagstat.txt"
    $SING_EXEC samtools index "${alignoutdir}/${fileName}.DeDup.bam"
    
    # 生成bigWig文件
    echo "Generating BigWig files..."
    $SING_EXEC bamCoverage -b "${alignoutdir}/${fileName}.DeDup.bam" \
        -o "${bwoutdir}/${fileName}.reverse.bw" \
        --filterRNAstrand forward --ignoreForNormalization chrM \
        -p ${threads} --binSize ${binSize} --normalizeUsing CPM
    
    $SING_EXEC bamCoverage -b "${alignoutdir}/${fileName}.DeDup.bam" \
        -o "${bwoutdir}/${fileName}.forward.bw" \
        --filterRNAstrand reverse --ignoreForNormalization chrM \
        -p ${threads} --binSize ${binSize} --normalizeUsing CPM
done

########### MultiQC
echo "Running MultiQC..."
$SING_EXEC multiqc "${FastQCdir}"/* "${trimedfadir}"/* "${alignoutdir}"/*.Log.final.out \
    "${alignoutdir}"/*flagstat.txt -o "${multiqcdir}" --force

########### Peak Calling (如果不需要跳过)
if [[ "$skip_peak" = false ]]; then
    echo "Starting peak calling analysis..."
    
    # 获取IP和Input的BAM文件
    IP_BAM_NAME=$(basename "$IpFq" | sed -E 's/\.f(ast)?q\.gz$//')
    IP_BAM="${alignoutdir}/${IP_BAM_NAME}.DeDup.bam"
    
    if [[ -n "$InputFq" ]]; then
        INPUT_BAM_NAME=$(basename "$InputFq" | sed -E 's/\.f(ast)?q\.gz$//')
        INPUT_BAM="${alignoutdir}/${INPUT_BAM_NAME}.DeDup.bam"
    else
        echo "Error: Input control is required for peak calling with exomePeak2"
        exit 1
    fi
    
    # 检查BAM文件是否存在
    if [[ ! -f "$IP_BAM" ]]; then
        echo "Error: IP BAM file $IP_BAM does not exist. Cannot run peak calling."
        exit 1
    fi
    
    if [[ ! -f "$INPUT_BAM" ]]; then
        echo "Error: Input BAM file $INPUT_BAM does not exist. Cannot run peak calling."
        exit 1
    fi
    
    # 创建peak calling的R脚本
    R_SCRIPT=$(cat <<EOF
# 设置工作目录
setwd("$pcoutdir")

# 安装并加载必要的包
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
    library(BiocManager)
}

if (!require("exomePeak2")) {
    BiocManager::install("exomePeak2")
}

# 设置基因组版本对应的BSgenome包
genome_version <- "$genome"
if (genome_version == "hs") {
    genome_version <- "hg38"
    bsgenome_pkg <- "BSgenome.Hsapiens.UCSC.hg38"
} else if (genome_version == "mm") {
    genome_version <- "mm10"
    bsgenome_pkg <- "BSgenome.Mmusculus.UCSC.mm10"
} else {
    stop("Unsupported genome version: ", genome_version)
}

# 检查并安装BSgenome包
if (!require(bsgenome_pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing required genome package:", bsgenome_pkg, "\\n")
    BiocManager::install(bsgenome_pkg, update = FALSE, ask = FALSE)
    if (!require(bsgenome_pkg, character.only = TRUE, quietly = TRUE)) {
        stop("Failed to install or load: ", bsgenome_pkg)
    }
}

library(exomePeak2)

# 设置文件路径
ip_bam <- "$IP_BAM"
input_bam <- "$INPUT_BAM"
GENE_ANNO_GTF <- "$gtf"

cat("Running exomePeak2 with parameters:\\n")
cat("IP BAM:", ip_bam, "\\n")
cat("Input BAM:", input_bam, "\\n")
cat("Genome:", genome_version, "\\n")
cat("GTF:", GENE_ANNO_GTF, "\\n")

# 运行exomePeak2 - 只保留有input的情况
result <- exomePeak2(
    bam_ip = ip_bam,
    bam_input = input_bam,
    gff = GENE_ANNO_GTF,
    genome = genome_version,
    parallel = $threads
)
result
# 保存结果
saveRDS(result, file = "exomePeak2_result.rds")
cat("Peak calling completed successfully!\\n")
EOF
)

    # 将R脚本写入文件并运行
    echo "$R_SCRIPT" > "${pcoutdir}/run_exomePeak2.R"
    
    echo "Running exomePeak2 peak calling..."
    $SING_EXEC Rscript "${pcoutdir}/run_exomePeak2.R" 2>&1 | tee "${outputdir}/logs/exomePeak2.log"
    
    if [[ $? -eq 0 ]]; then
        echo "exomePeak2 analysis completed successfully!"
        echo "Peak calling results are available in: ${pcoutdir}"
        
        ########### Motif Analysis (如果不需要跳过)
        if [[ "$skip_motif" = false ]]; then
            echo "Starting motif analysis..."
            
            # 设置参考基因组文件路径
            # 这里需要根据实际情况设置基因组fasta文件的路径
            # 假设基因组文件与参考目录在同一位置，且命名为GRCh38.p14.genome.fa
            genome_fasta="${referencedir}/../GRCh38.p14.genome.fa"
            
            if [[ ! -f "$genome_fasta" ]]; then
                echo "Warning: Genome fasta file not found at $genome_fasta"
                echo "Please provide the correct path to the genome fasta file for motif analysis"
                skip_motif=true
            fi
            
            if [[ "$skip_motif" = false ]]; then
                mkdir -p "${pcoutdir}/motif_analysis"
                motifdir="${pcoutdir}/motif_analysis"
                
                # 转换BED12到BED6格式
                echo "Converting BED12 to BED6 format..."
                $SING_EXEC bed12ToBed6 -i "${pcoutdir}/exomePeak2_output/peaks.bed" -n > "${motifdir}/bed6_peak.bed"
                
                # 索引基因组文件（如果尚未索引）
                echo "Indexing genome fasta file..."
                $SING_EXEC samtools faidx "$genome_fasta"
                
                # 从BED文件提取fasta序列
                echo "Extracting fasta sequences from peaks..."
                $SING_EXEC bed2fasta -s -both -o "${motifdir}/peak.fa" "${motifdir}/bed6_peak.bed" "$genome_fasta"
                
                # 运行STREME进行motif发现
                echo "Running STREME for motif discovery..."
                $SING_EXEC streme --oc "${motifdir}" --rna --minw 5 --maxw 5 --thresh 0.1 --align center --p "${motifdir}/peak.fa"
                
                if [[ $? -eq 0 ]]; then
                    echo "Motif analysis completed successfully!"
                    echo "Motif results are available in: ${motifdir}"
                else
                    echo "Warning: STREME motif analysis failed. Check the log for details."
                fi
            fi
        fi
    else
        echo "Error: exomePeak2 analysis failed. Check the log for details."
        exit 1
    fi
fi