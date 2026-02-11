"""
Quality Control Module
======================

QC rules for virome data preprocessing.

Input: Raw paired-end reads
Output: Cleaned paired-end reads in {OUTDIR}/clean_reads/
"""

# ================================================================================
# FastQC and Quality Assessment
# ================================================================================

# Generate FastQC output patterns based on config naming convention
auto_config = config.get("sample_auto_detection", {})
r1_pattern = auto_config.get("r1_pattern", "*_R1.fastq.gz")
r2_pattern = auto_config.get("r2_pattern", "*_R2.fastq.gz")

# Extract base names by removing extensions
# "*_R1_001.fastq.gz" becomes "{sample}_R1_001"
r1_base = r1_pattern.replace("*", "{sample}").replace(".fastq.gz", "").replace(".fq.gz", "")
r2_base = r2_pattern.replace("*", "{sample}").replace(".fastq.gz", "").replace(".fq.gz", "")

rule fastqc_raw:
    """
    Run FastQC on raw reads to assess initial quality and identify artifacts
    """
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        html_r1 = f"{OUTDIR}/fastqc/raw/{r1_base}_fastqc.html",
        html_r2 = f"{OUTDIR}/fastqc/raw/{r2_base}_fastqc.html",
        zip_r1 = f"{OUTDIR}/fastqc/raw/{r1_base}_fastqc.zip",
        zip_r2 = f"{OUTDIR}/fastqc/raw/{r2_base}_fastqc.zip"
    log:
        f"{OUTDIR}/logs/fastqc_raw/{{sample}}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o $(dirname {output.html_r1}) {input.r1} {input.r2} 2>&1 | tee {log}
        """

rule clumpify_optical_duplicates:
    """
    Remove optical duplicates using Clumpify (BBTools)
    These are Illumina sequencing artifacts from patterned flow cells

    Note: Memory requirement scales with read count (~0.5GB per million reads)
          Allocated 48GB to handle samples up to ~95 million reads
    """
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        r1 = temp(f"{OUTDIR}/clumpify/{{sample}}_R1.fastq.gz"),
        r2 = temp(f"{OUTDIR}/clumpify/{{sample}}_R2.fastq.gz")
    log:
        f"{OUTDIR}/logs/clumpify/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 48000
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        clumpify.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            dedupe optical \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule fastp_trim:
    """
    CRITICAL STEP: fastp with NovaSeq-specific settings

    1. PolyG tail removal (NovaSeq 2-channel chemistry artifact) - MUST BE FIRST
    2. Adapter trimming (auto-detect Illumina adapters)
    3. Quality filtering (Q20 threshold)
    4. Length filtering (100bp minimum due to library prep artifacts)
    5. Complexity filtering

    This handles most NovaSeq + Illumina library prep artifacts
    """
    input:
        r1 = f"{OUTDIR}/clumpify/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/clumpify/{{sample}}_R2.fastq.gz"
    output:
        r1 = temp(f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz"),
        r2 = temp(f"{OUTDIR}/fastp/{{sample}}_R2.fastq.gz"),
        html = f"{OUTDIR}/fastp/{{sample}}_fastp.html",
        json = f"{OUTDIR}/fastp/{{sample}}_fastp.json"
    log:
        f"{OUTDIR}/logs/fastp/{{sample}}.log"
    threads: 8
    params:
        min_length = config.get("min_read_length", 100),
        qual_threshold = config.get("quality_threshold", 20)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --trim_poly_g \
            --poly_g_min_len 10 \
            --detect_adapter_for_pe \
            --qualified_quality_phred {params.qual_threshold} \
            --length_required {params.min_length} \
            --low_complexity_filter \
            --complexity_threshold 30 \
            --correction \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            2>&1 | tee {log}
        """

rule fastqc_trimmed:
    """
    FastQC on trimmed reads to verify polyG removal and adapter trimming success
    """
    input:
        r1 = f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/fastp/{{sample}}_R2.fastq.gz"
    output:
        html_r1 = f"{OUTDIR}/fastqc/trimmed/{{sample}}_R1_fastqc.html",
        html_r2 = f"{OUTDIR}/fastqc/trimmed/{{sample}}_R2_fastqc.html",
        zip_r1 = f"{OUTDIR}/fastqc/trimmed/{{sample}}_R1_fastqc.zip",
        zip_r2 = f"{OUTDIR}/fastqc/trimmed/{{sample}}_R2_fastqc.zip"
    log:
        f"{OUTDIR}/logs/fastqc_trimmed/{{sample}}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o $(dirname {output.html_r1}) {input.r1} {input.r2} 2>&1 | tee {log}
        """

rule flag_phix:
    """
    Flag PhiX174 contamination without removing reads

    For VLP-enriched viromes, we flag contamination rather than remove it to avoid
    accidentally discarding legitimate viral sequences. This generates QC metrics
    to identify samples with abnormal PhiX levels (which may indicate sequencing issues).

    Uses BBDuk in detection-only mode to quantify PhiX contamination.
    """
    input:
        r1 = f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/fastp/{{sample}}_R2.fastq.gz"
    output:
        stats = f"{OUTDIR}/contamination_flagging/phix/{{sample}}_phix_stats.txt"
    log:
        f"{OUTDIR}/logs/flag_phix/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 8000
    params:
        phix_ref = REFERENCES["phix"]
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        # Run bbduk in stats-only mode (no read filtering)
        # Redirect stdout to /dev/null since we only want stats
        bbduk.sh \
            in1={input.r1} in2={input.r2} \
            ref={params.phix_ref} \
            k=31 hdist=1 \
            stats={output.stats} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule flag_univec:
    """
    Flag vector/plasmid contamination without removing reads

    Checks for common cloning vectors, plasmids, and synthetic sequences that may
    indicate laboratory contamination. Generates QC metrics to identify outlier samples.

    Uses comprehensive UniVec-derived database including common lab vectors.
    """
    input:
        r1 = f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/fastp/{{sample}}_R2.fastq.gz"
    output:
        stats = f"{OUTDIR}/contamination_flagging/univec/{{sample}}_univec_stats.txt"
    log:
        f"{OUTDIR}/logs/flag_univec/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 8000
    params:
        vector_ref = REFERENCES["vector_contaminants"]
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        # Run bbduk in stats-only mode to detect vector contamination
        bbduk.sh \
            in1={input.r1} in2={input.r2} \
            ref={params.vector_ref} \
            k=31 hdist=1 \
            stats={output.stats} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule remove_primer_b_forward:
    """
    Remove forward primer B sequences using left-side trimming

    Step 1 of 2-step primer B removal process. Removes forward primer B sequences
    that may contaminate the 5' end of reads. Each sample should have primarily one
    primer B variant; presence of multiple variants indicates cross-contamination.

    This step generates stats for downstream cross-contamination analysis.
    """
    input:
        r1 = f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/fastp/{{sample}}_R2.fastq.gz"
    output:
        r1 = temp(f"{OUTDIR}/primer_b/step1/{{sample}}_R1.fastq.gz"),
        r2 = temp(f"{OUTDIR}/primer_b/step1/{{sample}}_R2.fastq.gz"),
        stats = f"{OUTDIR}/primer_b/stats/step1/{{sample}}_pb_fwd_stats.txt"
    log:
        f"{OUTDIR}/logs/primer_b/{{sample}}_step1.log"
    threads: 4
    resources:
        mem_mb = 8000
    params:
        ref = REFERENCES["primer_b_forward"],
        ktrim = config["primer_b"]["bbduk_forward"]["ktrim"],
        k = config["primer_b"]["bbduk_forward"]["k"],
        mink = config["primer_b"]["bbduk_forward"]["mink"],
        rcomp = config["primer_b"]["bbduk_forward"]["rcomp"],
        ordered = config["primer_b"]["bbduk_forward"]["ordered"],
        ow = config["primer_b"]["bbduk_forward"]["ow"]
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        bbduk.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            ref={params.ref} \
            stats={output.stats} \
            ktrim={params.ktrim} \
            k={params.k} \
            mink={params.mink} \
            rcomp={params.rcomp} \
            ordered={params.ordered} \
            ow={params.ow} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule remove_primer_b_rc:
    """
    Remove reverse complement primer B sequences using right-side trimming

    Step 2 of 2-step primer B removal process. Removes reverse complement primer B
    sequences that may contaminate the 3' end of reads after step 1.

    This step generates stats for downstream cross-contamination analysis.
    """
    input:
        r1 = f"{OUTDIR}/primer_b/step1/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/primer_b/step1/{{sample}}_R2.fastq.gz"
    output:
        r1 = temp(f"{OUTDIR}/primer_b/step2/{{sample}}_R1.fastq.gz"),
        r2 = temp(f"{OUTDIR}/primer_b/step2/{{sample}}_R2.fastq.gz"),
        stats = f"{OUTDIR}/primer_b/stats/step2/{{sample}}_pb_rc_stats.txt"
    log:
        f"{OUTDIR}/logs/primer_b/{{sample}}_step2.log"
    threads: 4
    resources:
        mem_mb = 8000
    params:
        ref = REFERENCES["primer_b_rc"],
        ktrim = config["primer_b"]["bbduk_rc"]["ktrim"],
        k = config["primer_b"]["bbduk_rc"]["k"],
        mink = config["primer_b"]["bbduk_rc"]["mink"],
        rcomp = config["primer_b"]["bbduk_rc"]["rcomp"],
        ordered = config["primer_b"]["bbduk_rc"]["ordered"],
        ow = config["primer_b"]["bbduk_rc"]["ow"]
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        bbduk.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            ref={params.ref} \
            stats={output.stats} \
            ktrim={params.ktrim} \
            k={params.k} \
            mink={params.mink} \
            rcomp={params.rcomp} \
            ordered={params.ordered} \
            ow={params.ow} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule host_depletion:
    """
    Remove host genome sequences using minimap2

    For VLP-enriched samples, high host contamination (>10%) indicates VLP prep failure
    This step serves dual purpose:
    1. Remove host sequences from analysis
    2. QC metric for VLP prep efficiency

    Note: Now operates after primer B removal (PhiX/vector flagging doesn't filter reads)
    Note: Minimap2 loads host genome index (~12GB for human) plus read processing
          Allocated 24GB to handle large samples
    """
    input:
        r1 = f"{OUTDIR}/primer_b/step2/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/primer_b/step2/{{sample}}_R2.fastq.gz",
        host_ref = REFERENCES["host_genome"],
        # Ensure contamination flagging happens first (for workflow ordering)
        phix_stats = f"{OUTDIR}/contamination_flagging/phix/{{sample}}_phix_stats.txt",
        univec_stats = f"{OUTDIR}/contamination_flagging/univec/{{sample}}_univec_stats.txt"
    output:
        r1 = temp(f"{OUTDIR}/host_depleted/{{sample}}_R1.fastq.gz"),
        r2 = temp(f"{OUTDIR}/host_depleted/{{sample}}_R2.fastq.gz"),
        stats = f"{OUTDIR}/host_depleted/{{sample}}_host_stats.txt"
    log:
        f"{OUTDIR}/logs/host_depletion/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 24000
    conda:
        "../envs/mapping.yaml"
    shell:
        """
        # Map to host genome
        minimap2 -ax sr -t {threads} {input.host_ref} {input.r1} {input.r2} 2>> {log} | \
        samtools view -@ {threads} -bS - | \
        samtools sort -@ {threads} -n -o - | \
        samtools fastq -@ {threads} -f 12 -F 256 \
            -1 {output.r1} -2 {output.r2} - 2>> {log}

        # Calculate host mapping stats
        echo "Sample: {wildcards.sample}" > {output.stats}
        echo "Host reference: {input.host_ref}" >> {output.stats}
        minimap2 -ax sr -t {threads} {input.host_ref} {input.r1} {input.r2} 2>> {log} | \
        samtools view -@ {threads} -c >> {output.stats}
        """

rule remove_rrna:
    """
    Remove ribosomal RNA sequences using BBDuk

    Uses SILVA database for comprehensive rRNA removal
    Critical for RT-based protocols where all RNA gets converted to cDNA

    Note: Full SILVA SSU+LSU database loads 123M kmers (~5GB)
          Additional memory needed for read processing on large samples
          Allocated 32GB to handle database + up to 95M reads

    Note: Files not marked as temp() to avoid broken symlinks in QC-only mode.
          These may be intermediate (when PCR dedup enabled) or final clean reads.
          Symlinks in clean_reads/ always point to the true final outputs.
    """
    input:
        r1 = f"{OUTDIR}/host_depleted/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/host_depleted/{{sample}}_R2.fastq.gz"
    output:
        r1 = f"{OUTDIR}/rrna_removed/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/rrna_removed/{{sample}}_R2.fastq.gz",
        stats = f"{OUTDIR}/rrna_removed/{{sample}}_rrna_stats.txt"
    log:
        f"{OUTDIR}/logs/rrna_removal/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 32000
    params:
        rrna_ref = REFERENCES["rrna"]
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        bbduk.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            ref={params.rrna_ref} \
            k=31 \
            stats={output.stats} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}
        """

rule remove_pcr_duplicates:
    """
    Remove PCR duplicates using clumpify.sh from BBTools

    For Round A/B virome data with 40 cycles of PCR, this step removes
    PCR duplicates, getting closer to "unique molecules present" rather
    than "PCR-amplified copies."

    Placement rationale:
    - Last QC step (after all quality filtering and rRNA removal)
    - Quality filtering has already removed error-containing reads
    - More accurate duplicate detection on clean viral reads
    - Prevents keeping "wrong" copies with sequencing errors

    Optional: Controlled by config['deduplication']['pcr']['enabled']
    - Enable for high-cycle PCR data (Round A/B with 40 cycles)
    - Disable for non-PCR or low-cycle data

    Note: Uses exact matching for duplicate detection (dedupe optical=f).
          This is appropriate here because quality filtering has already
          removed error-containing reads, so remaining sequence variation
          likely represents real biological diversity (viral quasispecies)
          rather than PCR/sequencing errors.

    Note: Files not marked as temp() to avoid broken symlinks in QC-only mode.
          When enabled, these are the final clean reads accessed via symlinks.

    Note: Input depends on rRNA removal toggle - uses rrna_removed if enabled,
          host_depleted if disabled.
    """
    input:
        r1 = lambda wc: (
            f"{OUTDIR}/rrna_removed/{wc.sample}_R1.fastq.gz"
            if config.get("rrna_removal", {}).get("enabled", True)
            else f"{OUTDIR}/host_depleted/{wc.sample}_R1.fastq.gz"
        ),
        r2 = lambda wc: (
            f"{OUTDIR}/rrna_removed/{wc.sample}_R2.fastq.gz"
            if config.get("rrna_removal", {}).get("enabled", True)
            else f"{OUTDIR}/host_depleted/{wc.sample}_R2.fastq.gz"
        )
    output:
        r1 = f"{OUTDIR}/pcr_deduplicated/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/pcr_deduplicated/{{sample}}_R2.fastq.gz",
        stats = f"{OUTDIR}/pcr_deduplicated/{{sample}}_pcr_dup_stats.txt"
    log:
        f"{OUTDIR}/logs/pcr_deduplication/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 32000  # Similar to rRNA removal, scales with read count
    conda:
        "../envs/bbtools.yaml"
    shell:
        """
        clumpify.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            dedupe optical=f \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            2>&1 | tee {log}

        # Extract stats for reporting
        # clumpify.sh outputs stats to stderr/stdout, captured in log
        grep -E "^(Input|Duplicates|Result)" {log} > {output.stats} || true
        """

# Quality assessment relies on directly measurable metrics: host_percent, rrna_percent, final_reads

rule final_fastqc:
    """
    Final FastQC on clean reads

    Uses the final QC step output based on config toggles.
    """
    input:
        r1 = lambda wc: get_clean_reads_source(wc.sample, "R1"),
        r2 = lambda wc: get_clean_reads_source(wc.sample, "R2")
    output:
        html_r1 = f"{OUTDIR}/fastqc/final/{{sample}}_R1_fastqc.html",
        html_r2 = f"{OUTDIR}/fastqc/final/{{sample}}_R2_fastqc.html",
        zip_r1 = f"{OUTDIR}/fastqc/final/{{sample}}_R1_fastqc.zip",
        zip_r2 = f"{OUTDIR}/fastqc/final/{{sample}}_R2_fastqc.zip"
    log:
        f"{OUTDIR}/logs/fastqc_final/{{sample}}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc -t {threads} -o $(dirname {output.html_r1}) {input.r1} {input.r2} 2>&1 | tee {log}
        """

def get_clean_reads_source(sample, read):
    """
    Determine the source of clean reads based on config toggles.

    Priority order (latest step that is enabled):
    1. pcr_deduplicated (if PCR dedup enabled)
    2. rrna_removed (if rRNA removal enabled)
    3. host_depleted (fallback when both disabled)
    """
    pcr_dedup_enabled = config.get("deduplication", {}).get("pcr", {}).get("enabled", False)
    rrna_enabled = config.get("rrna_removal", {}).get("enabled", True)

    if pcr_dedup_enabled:
        return f"{OUTDIR}/pcr_deduplicated/{sample}_{read}.fastq.gz"
    elif rrna_enabled:
        return f"{OUTDIR}/rrna_removed/{sample}_{read}.fastq.gz"
    else:
        return f"{OUTDIR}/host_depleted/{sample}_{read}.fastq.gz"

rule symlink_clean_reads:
    """
    Create symlinks to final clean reads for easy access

    Links to the appropriate QC step based on config toggles:
    - pcr_deduplicated (if PCR dedup enabled)
    - rrna_removed (if rRNA removal enabled, default)
    - host_depleted (if both rRNA and PCR dedup disabled)
    """
    input:
        r1 = lambda wc: get_clean_reads_source(wc.sample, "R1"),
        r2 = lambda wc: get_clean_reads_source(wc.sample, "R2")
    output:
        r1 = f"{OUTDIR}/clean_reads/{{sample}}_R1.fastq.gz",
        r2 = f"{OUTDIR}/clean_reads/{{sample}}_R2.fastq.gz"
    shell:
        """
        ln -sf $(readlink -f {input.r1}) {output.r1}
        ln -sf $(readlink -f {input.r2}) {output.r2}
        """

# ================================================================================
# Reporting Rules
# ================================================================================

rule count_reads:
    """
    Count reads at each QC step for tracking

    Tracks reads through the full pipeline including primer B removal
    Note: phix_removed step removed since we now flag instead of filter
    Note: Conditionally includes rrna_removed step if enabled
    Note: Conditionally includes pcr_deduplicated step if enabled
    Note: Can take several hours for large samples (50M+ reads)
    """
    input:
        raw_r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        clumpify_r1 = f"{OUTDIR}/clumpify/{{sample}}_R1.fastq.gz",
        fastp_r1 = f"{OUTDIR}/fastp/{{sample}}_R1.fastq.gz",
        pb_step1_r1 = f"{OUTDIR}/primer_b/step1/{{sample}}_R1.fastq.gz",
        pb_step2_r1 = f"{OUTDIR}/primer_b/step2/{{sample}}_R1.fastq.gz",
        host_r1 = f"{OUTDIR}/host_depleted/{{sample}}_R1.fastq.gz",
        rrna_r1 = lambda wc: f"{OUTDIR}/rrna_removed/{wc.sample}_R1.fastq.gz" if config.get("rrna_removal", {}).get("enabled", True) else [],
        pcr_dedup_r1 = lambda wc: f"{OUTDIR}/pcr_deduplicated/{wc.sample}_R1.fastq.gz" if config.get("deduplication", {}).get("pcr", {}).get("enabled", False) else [],
        clean_r1 = lambda wc: get_clean_reads_source(wc.sample, "R1")
    output:
        temp(f"{OUTDIR}/reports/read_counts/{{sample}}.tsv")
    resources:
        time_min = 360
    params:
        rrna_enabled = config.get("rrna_removal", {}).get("enabled", True),
        pcr_dedup_enabled = config.get("deduplication", {}).get("pcr", {}).get("enabled", False)
    conda:
        "../envs/qc.yaml"
    shell:
        """
        echo -e "sample\tstep\treads" > {output}
        echo -e "{wildcards.sample}\traw\t$(zcat {input.raw_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        echo -e "{wildcards.sample}\tclumpify\t$(zcat {input.clumpify_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        echo -e "{wildcards.sample}\tfastp\t$(zcat {input.fastp_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        echo -e "{wildcards.sample}\tprimer_b_fwd\t$(zcat {input.pb_step1_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        echo -e "{wildcards.sample}\tprimer_b_rc\t$(zcat {input.pb_step2_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        echo -e "{wildcards.sample}\thost_depleted\t$(zcat {input.host_r1} | wc -l | awk '{{print $1/4}}')" >> {output}

        # Conditionally add rRNA removal stats
        if [ "{params.rrna_enabled}" = "True" ]; then
            echo -e "{wildcards.sample}\trrna_removed\t$(zcat {input.rrna_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        fi

        # Conditionally add PCR dedup stats
        if [ "{params.pcr_dedup_enabled}" = "True" ]; then
            echo -e "{wildcards.sample}\tpcr_deduplicated\t$(zcat {input.pcr_dedup_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        fi

        # Final clean reads (from whatever the last enabled step is)
        echo -e "{wildcards.sample}\tclean\t$(zcat {input.clean_r1} | wc -l | awk '{{print $1/4}}')" >> {output}
        """

rule merge_read_counts:
    """
    Merge read counts from all samples
    """
    input:
        expand(f"{OUTDIR}/reports/read_counts/{{sample}}.tsv", sample=SAMPLES)
    output:
        f"{OUTDIR}/reports/read_counts.tsv"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        """

rule aggregate_contamination_stats:
    """
    Aggregate PhiX and vector contamination statistics across all samples

    Parses BBDuk stats files to create a summary table showing contamination
    levels per sample. This helps identify outlier samples with abnormal
    contamination that may indicate sequencing or library prep issues.
    """
    input:
        phix_stats = expand(f"{OUTDIR}/contamination_flagging/phix/{{sample}}_phix_stats.txt", sample=SAMPLES),
        univec_stats = expand(f"{OUTDIR}/contamination_flagging/univec/{{sample}}_univec_stats.txt", sample=SAMPLES)
    output:
        f"{OUTDIR}/reports/contamination_summary.tsv"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/aggregate_contamination_stats.py"

rule plot_contamination:
    """
    Generate contamination visualization plots

    Creates multiple plot types to help identify outlier samples:
    - Bar plot: contamination levels per sample
    - Box plot: distribution of contamination levels
    - Scatter plot: correlation between PhiX and vector contamination
    - Heatmap: overview of contamination across all samples
    """
    input:
        f"{OUTDIR}/reports/contamination_summary.tsv"
    output:
        f"{OUTDIR}/reports/contamination_bars.png",
        f"{OUTDIR}/reports/contamination_heatmap.png"
    params:
        output_prefix = f"{OUTDIR}/reports/contamination"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/plot_contamination.py"

rule analyze_primer_b_contamination:
    """
    Analyze primer B cross-contamination across all samples

    Parses BBDuk stats from primer B removal steps to detect cross-contamination.
    Each sample should have primarily one primer B variant; multiple variants
    indicate cross-contamination during library prep or sequencing.

    Outputs:
    - Summary table with dominant primer B and contamination metrics
    - Distribution data for heatmap visualization
    """
    input:
        step1 = expand(f"{OUTDIR}/primer_b/stats/step1/{{sample}}_pb_fwd_stats.txt", sample=SAMPLES),
        step2 = expand(f"{OUTDIR}/primer_b/stats/step2/{{sample}}_pb_rc_stats.txt", sample=SAMPLES)
    output:
        summary = f"{OUTDIR}/reports/primer_b_contamination_summary.tsv",
        distribution = f"{OUTDIR}/reports/primer_b_distribution.tsv"
    params:
        sample_assignments = config["primer_b"].get("sample_assignments", None),
        contamination_threshold = config["primer_b"]["contamination_threshold"]
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/analyze_primer_b_contamination.py"

rule plot_primer_b_heatmap:
    """
    Generate primer B cross-contamination heatmap

    Creates a heatmap showing primer B variant distribution across samples.
    Helps identify cross-contamination patterns and verify sample identity.
    """
    input:
        distribution = f"{OUTDIR}/reports/primer_b_distribution.tsv",
        summary = f"{OUTDIR}/reports/primer_b_contamination_summary.tsv"
    output:
        heatmap = f"{OUTDIR}/reports/primer_b_heatmap.png"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/plot_primer_b_heatmap.py"

rule qc_metrics:
    """
    Generate QC quality metrics for each sample based on:
    - % host reads (VLP prep efficiency)
    - % rRNA reads (biological contamination)
    - Final read count after QC (sequencing depth)

    Provides quality categories rather than binary pass/fail assessments.
    """
    input:
        read_counts = f"{OUTDIR}/reports/read_counts.tsv"
    output:
        f"{OUTDIR}/reports/sample_qc_metrics.tsv"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/generate_qc_metrics.py"

rule multiqc:
    """
    Aggregate all QC reports into single MultiQC dashboard

    Now includes contamination flagging summary
    """
    input:
        # FastQC reports
        expand(f"{OUTDIR}/fastqc/raw/{r1_base}_fastqc.zip", sample=SAMPLES),
        expand(f"{OUTDIR}/fastqc/trimmed/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{OUTDIR}/fastqc/final/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        # fastp reports
        expand(f"{OUTDIR}/fastp/{{sample}}_fastp.json", sample=SAMPLES),
        # Read counts
        f"{OUTDIR}/reports/read_counts.tsv",
        # QC metrics
        f"{OUTDIR}/reports/sample_qc_metrics.tsv",
        # Contamination summary
        f"{OUTDIR}/reports/contamination_summary.tsv"
    output:
        f"{OUTDIR}/multiqc/multiqc_report.html"
    log:
        f"{OUTDIR}/logs/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc \
            {OUTDIR} \
            -o {OUTDIR}/multiqc \
            -n multiqc_report.html \
            --force \
            2>&1 | tee {log}
        """

rule virome_report:
    """
    Generate user-friendly virome QC reports (replaces MultiQC)

    Creates dual HTML outputs for different use cases:
    - virome_report.html: Lightweight version with external images (~2-3MB)
    - virome_report_standalone.html: Self-contained with embedded images (~16MB)
    Plus structured JSON data for virome-specific quality assessment.
    """
    input:
        # QC data sources
        read_counts=f"{OUTDIR}/reports/read_counts.tsv",
        contamination_summary=f"{OUTDIR}/reports/contamination_summary.tsv",
        qc_metrics=f"{OUTDIR}/reports/sample_qc_metrics.tsv",
        primer_b_summary=f"{OUTDIR}/reports/primer_b_contamination_summary.tsv",
        # Plot files for report
        contamination_plots=[
            f"{OUTDIR}/reports/contamination_bars.png",
            f"{OUTDIR}/reports/contamination_heatmap.png"
        ],
        primer_b_plot=f"{OUTDIR}/reports/primer_b_heatmap.png",
        # Optional: FastQC files for additional metrics
        fastqc_files=expand(f"{OUTDIR}/fastqc/final/{{sample}}_R1_fastqc.zip", sample=SAMPLES)
    output:
        html_report=f"{OUTDIR}/reports/virome_report.html",
        html_report_standalone=f"{OUTDIR}/reports/virome_report_standalone.html",
        json_data=f"{OUTDIR}/reports/virome_report_data.json"
    log:
        f"{OUTDIR}/logs/virome_report.log"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/generate_virome_report.py"
