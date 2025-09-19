# MeRIP-seq-Processing-Pipeline
The MeRIP-seq analysis pipeline processes raw FASTQ data through a series of steps including adapter trimming, quality control, genome mapping, peak calling, and motif analysis.

# Part I Introduction
## i. Workflow
Here stands an throughout workflow of MeRIP-seq data analysis.
<img width="1561" height="328" alt="{60F67D9B-A185-4AD0-A189-2B85B24551FC}" src="https://github.com/user-attachments/assets/9ec8a35c-4dd1-43d5-92fb-2997df464168" />



## ii. Features
This pipeline provides a fully containerized Singularity environment that bundles all required tools and dependencies. With a single command, the entire MeRIP-seq workflow—from raw FASTQ input through trimming, quality control, genome alignment, peak calling, and motif analysis—can be executed reproducibly on any compatible system.

# Part II Requirements
1.  **Recommended System Configuration**:

      * 8-core CPU
      * 24 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **Download Files**:

      * `run_MeRIP-seq.sh`
      * `MeRIPseq.sif` (The Singularity container)
	  
4.  **Reference Data**: A directory containing bowtie index (Below are the detailed steps for the human hg38 genome. For other reference genomes, please download the corresponding files and replace them as needed).
      ```bash
      mkdir basement_data
      cd basement_data
      # Download Genome FASTA and GTF
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.p14.genome.fa.gz
	  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz
      # Unzip the files
      gunzip GRCh38.p14.genome.fa.gz
      gunzip gencode.v46.annotation.gtf.gz
      # Remove scafford
      awk '/^>/ {p=0} /^>chr[0-9XYM]/ {p=1} p' GRCh38.primary_assembly.genome.fa > GRCh38.primary_assembly.genome.chr.fa
      # Build index
      mkdir hg38_chr_star_index 
	  singularity exec MeRIPseq.sif STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./hg38_chr_star_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v46.primary_assembly.annotation.gtf \
     --sjdbOverhang 149
      # Remove unnecessary files
      rm GRCh38.primary_assembly.genome.chr.fa
      rm GRCh38.primary_assembly.genome.fa
      ```
6.   **Required File Structure**
      ```bash
      basement_data/
      ├── MeRIPseq.sif
      ├── hg38_chr_star_index/
      └── run_MeRIP-seq.sh
      ```

# Part III Running

   * **Example code**

      ```bash
      bash ./final_all_merip.sh --IpFq ./MERIP-ip_test.fastq.gz \
                                     --InputFq ./MERIP-input_test.fastq.gz \
                                     --outputdir /mnt/guanli/lyc.MeRIP.singularity/newresult \
                                     --referencedir ./hg38_chr_star_index \
                                     --gtf /mnt/guanli/lyc.MeRIP.singularity/gencode.v46.annotation.gtf \
                                     --sif ./MeRIPseq.sif --threads 20 --binSize 10 --genome hs 
      ```
   * **Command Parameters**

      - `--IpFq          : Path to the treatment fastq (required)"
      - `--InputFq       : (optional) Path to the input control fastq"
      - `--outputdir     : Path to the output directory (required)"
      - `--referencedir  : Path to the directory where bowtie reference build with prefix (required)"
      - `--gtf           : Path to the GTF annotation file (required for peak calling)"
      - `--sif           : Path to the singularity environment file (required)"
      - `--threads       : (optional) Number of threads to use, default 8"
      - `--binSize       : (optional) Number of binsize to use, default 10"
      - `--genome        : (optional) Genome version (hs for human, mm for mouse), default hs"
      - `--skip_peak     : (optional) Skip peak calling step"

# Part IV Output

   * **Output Structure**
      ```bash
      result/
      ├── bam
            ├── MERIP-input_test.Aligned.sortedByCoord.out.bam
            ├── MERIP-input_test.DeDup.bam
            ├── MERIP-input_test.DeDup.bam.bai
            ├── MERIP-input_test.flagstat.txt
            ├── MERIP-input_test.Log.final.out
            ├── MERIP-input_test.Log.out
            ├── MERIP-input_test.Log.progress.out
            ├── MERIP-input_test.markdup.log
            ├── MERIP-input_test.SJ.out.tab
            └── MERIP-input_test.STAR.log
			├── MERIP-ip_test.Aligned.sortedByCoord.out.bam
			├── MERIP-ip_test.DeDup.bam
			├── MERIP-ip_test.DeDup.bam.bai
			├── MERIP-ip_test.flagstat.txt
			├── MERIP-ip_test.Log.final.out
			├── MERIP-ip_test.Log.out
			├── MERIP-ip_test.Log.progress.out
			├── MERIP-ip_test.markdup.log
			└── MERIP-ip_test.SJ.out.tab
			└── MERIP-ip_test.STAR.log
      ├── bw/
            ├── MERIP-input_test.forward.bw
            └── MERIP-input_test.reverse.bw
			├── MERIP-ip_test.forward.bw
			└── MERIP-ip_test.reverse.bw
      ├── multiqc/
            ├── multiqc_data/
            └── multiqc_report.html
      └── peak/
            ├── exomePeak2_output/
            ├── motif_analysis/
            ├── exomePeak2_result.rds
            └── run_exomePeak2.R
          
      ```
   * **Output Interpretation**

      - **`*.Log.final.out`**

        - **Content**: Contains STAR alignment summary statistics, including the total number of reads processed, reads aligned, reads discarded, and uniquely mapped reads. It provides an overview of mapping quality and efficiency for each FASTQ file.
        - **Application**: Used to assess alignment quality and sequencing library performance. These statistics help in troubleshooting mapping issues, evaluating experiment success, and can be parsed by downstream tools like MultiQC for visualization and comparison across samples.
          
          <img width="795" height="839" alt="{B7A76EDD-A701-4770-8CA1-617388366207}" src="https://github.com/user-attachments/assets/802e4f57-1028-4721-a7c8-4afa8479c92a" />


      - **`*.DeDup.bam`**

        - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bam.html.
        - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.bw`**

        - **Content**: A BigWig file that represents the End-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format. For more information please refer to: https://genome.ucsc.edu/goldenpath/help/bigWig.html
        - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes.

      - **`*.dedup.bam.bai`**

        - **Content**: This is the index file for the BAM file.
        - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

      - **`*.flagstat.txt`**

        - **Content**: Contains alignment statistics generated by samtools flagstat after removing duplicate reads. It reports the total number of reads, mapped reads, properly paired reads, singletons, and the number of duplicate reads removed, providing a summary of the final, deduplicated BAM file.
        - **Application**: Used to evaluate the quality of the deduplicated alignment, check library complexity, and ensure that downstream analyses (e.g., peak calling, coverage calculation) are based on high-quality, non-redundant reads.

         <img width="542" height="366" alt="{59D74160-7FE3-4DFD-9833-4D538CA679BE}" src="https://github.com/user-attachments/assets/0a992e27-de49-415a-bdd6-dfdddcd70720" />


      - **`*.markdup.log`**

        - **Content**: Log file generated by `samtools markdup`, summarizing read duplication. It includes READ (total number of input reads), WRITTEN (reads retained after removing duplicates), EXCLUDED, EXAMINED, counts of PAIRED and SINGLE reads, as well as DUPLICATE SINGLE/PAIR and DUPLICATE TOTAL.
        - **Application**: Used to evaluate library complexity and duplication rate. A high WRITTEN/READ ratio indicates low duplication and good library complexity, while a low ratio suggests high PCR duplication or low-complexity sequencing.
        
		  <img width="328" height="341" alt="{3BF121D4-0171-49B3-AAB7-3F4D471F2B4E}" src="https://github.com/user-attachments/assets/b58e3811-611f-4b90-bb44-744edee64b0c" />

      - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively.

        - **General Statistics**: A combined table summarizing important metrics for each sample:
	  
          <img width="2071" height="193" alt="{45F2B3B7-CD78-40D9-9BFC-68C45834C1E8}" src="https://github.com/user-attachments/assets/afd221e9-2e25-4041-b260-cc06d6e60560" />


        - **FastQC**: Quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', 'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content', 'Sequence Length Distribution', 'Sequence Duplication Levels', 'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content':
	      
          - Sequence Quality Histograms: The mean quality value across each base position in the read.
	  
            <img width="2070" height="722" alt="{A50E5F0E-B257-4B26-ABC1-4EE27D5D5A46}" src="https://github.com/user-attachments/assets/2cd29260-a086-4755-ad41-7ff5446664cc" />


          - Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
	  
            <img width="2108" height="639" alt="{C31B6B0A-6CEA-4504-8202-4B1B0B6E7380}" src="https://github.com/user-attachments/assets/b7d8c1f2-a163-4f3f-b3e3-6e2675d82240" />


        - **Samtools**: This module parses the output from samtools flagstat to report the percentage of total, mapped, and properly paired reads, providing a summary of alignment quality. Helps evaluate the effectiveness of deduplication and ensures that downstream analyses (e.g., peak calling, coverage profiling) are based on unique, non-redundant reads.
	  
          <img width="2158" height="957" alt="{03AC7145-11CE-4CAD-BBC5-45E931A5F49C}" src="https://github.com/user-attachments/assets/4fea9875-4017-4ffe-ae95-81836f42f91a" />


        - **STAR**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:
	  
          <img width="2044" height="368" alt="{AF510E4F-CF64-4118-AC3A-6C49E456D585}" src="https://github.com/user-attachments/assets/2d305060-2870-41e8-b04f-3db66b3517bd" />


        - **`exomePeak2_output/peaks.csv`**
		- **peaks.bed**: 12-column BED format
        - **peaks.csv**:
        - **Columns**:
        | Column | Description |
        |--------|-------------|
        | chr    | Chromosome name |
        | chromStart  | Start position of the broad peak (0-based) |
        | chromEnd    | End position of the broad peak (not inclusive) |
        | name   | Unique identifier for each peak |
        | strand  | Strand orientation of the peak |
        | blockCount | Number of exonic blocks comprising the peak |
        | blockSizes | Comma-separated sizes of each exonic block |
        | blockStarts | Comma-separated relative start positions of each block |
        | RPM.IP | Normalized read coverage in the Immunoprecipitation (IP) sample |
        | RPM.input | Normalized read coverage in the Input control sample |
        | geneID | Gene ID associated with the peak |
        | log2FC | Logarithm base 2 of the fold change between IP and Input |
		| pvalue | Statistical significance of peak enrichment |
		| fdr | Adjusted p-value controlling for multiple testing |
		| score | Composite score representing peak quality and significance |

          <img width="1016" height="423" alt="{C6D8B9EB-BEDF-4F72-9C72-3DAE97C72D0F}" src="https://github.com/user-attachments/assets/5e790d8b-41e5-4198-8af4-8a70b8b31f38" />
        - **gc_fit.pdf**: This file visualizes the GC-content bias correction model fitted by exomePeak2 to normalize sequencing coverage variations caused by GC composition differences.
	      <img width="590" height="364" alt="{5AD64503-2C51-4B3B-83A3-F879DE3D91D2}" src="https://github.com/user-attachments/assets/c7dcc331-d38e-409c-905c-1b8e96cd2030" />

	    - **`motif_analysis/streme.html`**
          <img width="1585" height="526" alt="{52337D6C-9973-4DE8-A663-1A76A6888940}" src="https://github.com/user-attachments/assets/baaf50ea-5b72-4c30-95bc-ea22f1ad4703" />


