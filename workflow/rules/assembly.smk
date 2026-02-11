"""
Assembly Module - Two-Stage Workflow
=====================================

Two-stage assembly workflow for virome contigs based on validated research
showing that global coassembly performs poorly. This workflow supports:
- by_sample (default): One MEGAHIT assembly per sample
- by_group: One MEGAHIT assembly per user-defined group

Workflow (same for both strategies):
1. Stage 1: MEGAHIT assemblies (per-sample or per-group)
2. Stage 2: Rename contigs with unique IDs (sample/group prefix)
3. Stage 3: Concatenate all renamed contigs
4. Stage 4: Flye meta-assembly for final polished assembly

NOTE: Global coassembly has been removed - research shows it performs worst
on all assembly metrics. Use 'by_sample' or 'by_group' instead.

Input: Cleaned paired-end reads (from QC module or user-provided)
Output: Final Flye meta-assembled contigs + assembly stats
"""

# ================================================================================
# Helper Functions for Entry Point Logic
# ================================================================================

def get_reads_for_assembly(wildcards):
    """
    Get cleaned reads for assembly based on pipeline entry point.

    Returns paths to R1/R2 files depending on where pipeline starts.
    Uses clean_reads symlinks which automatically point to the correct
    final QC step based on config toggles (rRNA removal, PCR dedup).
    """
    start_from = config["pipeline"].get("start_from", "raw_reads")

    if start_from == "raw_reads":
        # Use QC module outputs via clean_reads symlinks
        # These automatically point to the correct final step based on config
        return {
            "r1": f"{OUTDIR}/clean_reads/{wildcards.sample}_R1.fastq.gz",
            "r2": f"{OUTDIR}/clean_reads/{wildcards.sample}_R2.fastq.gz"
        }
    elif start_from == "cleaned_reads":
        # Use user-provided cleaned reads
        cleaned_dir = config["pipeline"].get("cleaned_reads_dir")
        if not cleaned_dir:
            raise ValueError(
                "Pipeline start_from='cleaned_reads' but 'cleaned_reads_dir' not specified in config!"
            )
        return {
            "r1": f"{cleaned_dir}/{wildcards.sample}_R1.fastq.gz",
            "r2": f"{cleaned_dir}/{wildcards.sample}_R2.fastq.gz"
        }
    else:
        raise ValueError(
            f"Invalid start_from value: '{start_from}'. "
            f"Must be 'raw_reads' or 'cleaned_reads'"
        )


def get_group_reads_for_assembly(wildcards):
    """
    Get reads for all samples in a group (for by_group strategy).

    Returns dict with lists of R1/R2 files for all samples in the group.
    Uses clean_reads symlinks which automatically point to the correct
    final QC step based on config toggles (rRNA removal, PCR dedup).
    """
    group_id = wildcards.group
    group_samples = GROUPS_TO_SAMPLES.get(group_id, [])

    start_from = config["pipeline"].get("start_from", "raw_reads")

    if start_from == "raw_reads":
        base_dir = f"{OUTDIR}/clean_reads"
    elif start_from == "cleaned_reads":
        base_dir = config["pipeline"].get("cleaned_reads_dir")
    else:
        raise ValueError(f"Invalid start_from value: '{start_from}'")

    return {
        "r1": [f"{base_dir}/{sample}_R1.fastq.gz" for sample in group_samples],
        "r2": [f"{base_dir}/{sample}_R2.fastq.gz" for sample in group_samples]
    }


def get_renamed_contigs_for_concatenation(wildcards):
    """
    Get all renamed contig files for concatenation based on assembly strategy.

    Args:
        wildcards: Snakemake wildcards object (unused but required for input functions)

    Returns list of renamed contig file paths.
    """
    strategy = config["pipeline"].get("assembly_strategy", "by_sample")

    if strategy == "by_sample":
        return expand(
            f"{OUTDIR}/assembly/per_sample/{{sample}}/renamed.contigs.fa",
            sample=SAMPLES
        )
    elif strategy == "by_group":
        return expand(
            f"{OUTDIR}/assembly/per_group/{{group}}/renamed.contigs.fa",
            group=GROUP_IDS
        )
    else:
        raise ValueError(f"Unknown assembly strategy: {strategy}")


# ================================================================================
# Stage 1: MEGAHIT Assemblies
# ================================================================================

rule megahit_individual:
    """
    Individual sample assembly using MEGAHIT (Stage 1 of two-stage workflow)

    Assembles each sample separately using paired-end reads directly.
    This is the default strategy and works well for most virome studies.
    The individual assemblies are then renamed and combined via Flye meta-assembly.
    """
    input:
        unpack(get_reads_for_assembly)
    output:
        contigs = f"{OUTDIR}/assembly/per_sample/{{sample}}/final.contigs.fa",
        done = touch(f"{OUTDIR}/assembly/per_sample/{{sample}}/.done")
    params:
        out_dir = f"{OUTDIR}/assembly/per_sample/{{sample}}",
        min_contig = config.get("assembly", {}).get("min_contig_length", 1000)
    log:
        f"{OUTDIR}/logs/assembly/megahit_{{sample}}.log"
    threads: 12
    resources:
        mem_mb = 50000  # 50GB per sample
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # Remove output directory if it exists
        rm -rf {params.out_dir}

        # Run megahit with paired-end reads
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &> {log}

        # Verify output
        [ -f {output.contigs} ] || {{ echo "Error: Assembly failed" >> {log}; exit 1; }}
        """


rule megahit_custom_groups:
    """
    Group-based assembly using MEGAHIT (Stage 1 of two-stage workflow)

    Assembles all samples within a user-defined group together using paired-end
    reads directly. Use this when samples are biologically related (e.g., same
    patient, same treatment).

    Groups are defined in a TSV file specified by pipeline.groups_file.
    """
    input:
        unpack(get_group_reads_for_assembly)
    output:
        contigs = f"{OUTDIR}/assembly/per_group/{{group}}/final.contigs.fa",
        done = touch(f"{OUTDIR}/assembly/per_group/{{group}}/.done")
    params:
        out_dir = f"{OUTDIR}/assembly/per_group/{{group}}",
        min_contig = config.get("assembly", {}).get("min_contig_length", 1000)
    log:
        f"{OUTDIR}/logs/assembly/megahit_group_{{group}}.log"
    threads: 24
    resources:
        mem_mb = 100000  # 100GB for group assemblies
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # Remove output directory if it exists
        rm -rf {params.out_dir}

        # Create comma-separated file lists for MEGAHIT
        R1=$(echo {input.r1} | tr ' ' ',')
        R2=$(echo {input.r2} | tr ' ' ',')

        # Run megahit with multiple paired-end input files
        megahit \
            -1 $R1 \
            -2 $R2 \
            -o {params.out_dir} \
            --min-contig-len {params.min_contig} \
            --k-list 21,29,39,59,79,99,119,141 \
            -t {threads} \
            &> {log}

        # Verify output
        [ -f {output.contigs} ] || {{ echo "Error: Assembly failed" >> {log}; exit 1; }}
        """


# ================================================================================
# Stage 2: Rename Contigs with Unique IDs
# ================================================================================

rule rename_contigs_individual:
    """
    Rename contigs with sample-specific prefixes (Stage 2 of two-stage workflow)

    Adds sample name prefix to all contig headers to ensure unique IDs when
    contigs from multiple samples are concatenated.

    Example: >contig_1 becomes >sample001contig_1

    Note: Empty assemblies are allowed (creates empty output with warning).
    This can happen with low-complexity or low-read-count samples.
    """
    input:
        contigs = f"{OUTDIR}/assembly/per_sample/{{sample}}/final.contigs.fa"
    output:
        renamed = f"{OUTDIR}/assembly/per_sample/{{sample}}/renamed.contigs.fa"
    log:
        f"{OUTDIR}/logs/assembly/rename_contigs_{{sample}}.log"
    threads: 2
    resources:
        mem_mb = 8000
    shell:
        """
        echo "Renaming contigs with prefix '{wildcards.sample}'" > {log}

        # Check if input has contigs
        ORIGINAL_COUNT=$(grep -c '^>' {input.contigs} || echo 0)
        echo "Original contigs: $ORIGINAL_COUNT" >> {log}

        if [ "$ORIGINAL_COUNT" -eq 0 ]; then
            echo "WARNING: Sample '{wildcards.sample}' has no contigs - creating empty output" >> {log}
            touch {output.renamed}
        else
            sed 's/>/>{wildcards.sample}/' {input.contigs} > {output.renamed} 2>> {log}
            RENAMED_COUNT=$(grep -c '^>' {output.renamed} || echo 0)
            echo "Renamed contigs: $RENAMED_COUNT" >> {log}
        fi
        """


rule rename_contigs_custom_groups:
    """
    Rename contigs with group-specific prefixes (Stage 2 of two-stage workflow)

    Adds group name prefix to all contig headers to ensure unique IDs when
    contigs from multiple groups are concatenated.

    Example: >contig_1 becomes >patient_Acontig_1

    Note: Empty assemblies are allowed (creates empty output with warning).
    """
    input:
        contigs = f"{OUTDIR}/assembly/per_group/{{group}}/final.contigs.fa"
    output:
        renamed = f"{OUTDIR}/assembly/per_group/{{group}}/renamed.contigs.fa"
    log:
        f"{OUTDIR}/logs/assembly/rename_contigs_group_{{group}}.log"
    threads: 2
    resources:
        mem_mb = 8000
    shell:
        """
        echo "Renaming contigs with prefix '{wildcards.group}'" > {log}

        # Check if input has contigs
        ORIGINAL_COUNT=$(grep -c '^>' {input.contigs} || echo 0)
        echo "Original contigs: $ORIGINAL_COUNT" >> {log}

        if [ "$ORIGINAL_COUNT" -eq 0 ]; then
            echo "WARNING: Group '{wildcards.group}' has no contigs - creating empty output" >> {log}
            touch {output.renamed}
        else
            sed 's/>/>{wildcards.group}/' {input.contigs} > {output.renamed} 2>> {log}
            RENAMED_COUNT=$(grep -c '^>' {output.renamed} || echo 0)
            echo "Renamed contigs: $RENAMED_COUNT" >> {log}
        fi
        """


# ================================================================================
# Stage 3: Concatenate Renamed Contigs
# ================================================================================

rule concatenate_all_contigs:
    """
    Concatenate all renamed contigs for Flye meta-assembly (Stage 3)

    Combines all renamed contig files (from individual samples or groups)
    into a single file for the Flye meta-assembly step.

    Output is strategy-specific to allow comparison between strategies.
    """
    input:
        get_renamed_contigs_for_concatenation
    output:
        concatenated = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/concatenated_contigs.fa"
    log:
        f"{OUTDIR}/logs/assembly/{ASSEMBLY_STRATEGY}/concatenate_contigs.log"
    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        echo "Concatenating renamed contigs from {input}" > {log}

        # Concatenate all input files
        cat {input} > {output.concatenated} 2>> {log}

        # Verify output
        [ -s {output.concatenated} ] || {{ echo "Error: Concatenated file is empty" >> {log}; exit 1; }}

        # Log statistics
        TOTAL_CONTIGS=$(grep -c '^>' {output.concatenated} || echo 0)
        TOTAL_SIZE=$(stat -f%z {output.concatenated} 2>/dev/null || stat --printf="%s" {output.concatenated})
        echo "Total concatenated contigs: $TOTAL_CONTIGS" >> {log}
        echo "Total file size: $TOTAL_SIZE bytes" >> {log}
        """


# ================================================================================
# Stage 4: Flye Meta-Assembly
# ================================================================================

rule flye_meta_assembly:
    """
    Final meta-assembly using Flye (Stage 4 of two-stage workflow)

    Uses Flye's --subassemblies mode to polish and combine the concatenated
    contigs from all samples/groups into a final high-quality assembly.

    Parameters are based on the validated hecatomb approach:
    - --subassemblies: Input mode for pre-assembled contigs
    - --plasmids: Include plasmid detection (important for viromes)
    - -g 1g: Estimated genome size (1 Gbp placeholder for metagenomes)

    Output is strategy-specific to allow comparison between strategies.

    Resources: 16-24 threads, 64-100GB memory, 24-hour time limit
    """
    input:
        contigs = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/concatenated_contigs.fa"
    output:
        assembly = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye/assembly.fasta",
        info = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye/assembly_info.txt"
    params:
        out_dir = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye"
    log:
        f"{OUTDIR}/logs/assembly/{ASSEMBLY_STRATEGY}/flye_meta.log"
    threads: 24
    resources:
        mem_mb = 100000,  # 100GB
        time_min = 1440   # 24 hours
    conda:
        "../envs/flye.yaml"
    shell:
        """
        echo "Starting Flye meta-assembly" > {log}
        echo "Input contigs: {input.contigs}" >> {log}
        echo "Output directory: {params.out_dir}" >> {log}

        # Remove output directory if it exists (Flye won't overwrite)
        rm -rf {params.out_dir}

        # Run Flye in subassemblies mode
        # Parameters from validated hecatomb approach
        flye --subassemblies {input.contigs} \
             -t {threads} \
             --plasmids \
             -o {params.out_dir} \
             -g 1g \
             2>&1 | tee -a {log}

        # Verify outputs
        [ -f {output.assembly} ] || {{ echo "Error: Flye assembly failed - no output" >> {log}; exit 1; }}
        [ -f {output.info} ] || {{ echo "Warning: Assembly info file not found" >> {log}; }}

        # Log statistics
        CONTIG_COUNT=$(grep -c '^>' {output.assembly} || echo 0)
        echo "Final assembly contig count: $CONTIG_COUNT" >> {log}
        """


# ================================================================================
# Create Convenience Symlinks
# ================================================================================

rule link_final_assembly:
    """
    Create convenience symlink for final Flye assembly output

    Makes the final assembly easily accessible at a standard location for
    downstream pipelines (e.g., phage-analysis pipeline).

    Output is strategy-specific to allow comparison between strategies.
    """
    input:
        f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye/assembly.fasta"
    output:
        f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/final.contigs.fa"
    shell:
        """
        ln -sf flye/assembly.fasta {output}
        """


# ================================================================================
# Assembly Statistics and QC
# ================================================================================

rule assembly_stats:
    """
    Calculate assembly statistics

    Metrics calculated for the final Flye assembly:
    - Number of contigs
    - Total assembly size
    - N50, L50
    - Longest contig
    - GC content

    Also includes per-sample/per-group MEGAHIT stats for comparison.

    Output is strategy-specific to allow comparison between strategies.
    """
    input:
        flye_assembly = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye/assembly.fasta",
        flye_info = f"{OUTDIR}/assembly/{ASSEMBLY_STRATEGY}/flye/assembly_info.txt",
        megahit_contigs = (
            expand(f"{OUTDIR}/assembly/per_sample/{{sample}}/final.contigs.fa", sample=SAMPLES)
            if config["pipeline"].get("assembly_strategy", "by_sample") == "by_sample"
            else expand(f"{OUTDIR}/assembly/per_group/{{group}}/final.contigs.fa", group=GROUP_IDS)
        )
    output:
        stats = f"{OUTDIR}/reports/assembly_stats_{ASSEMBLY_STRATEGY}.tsv"
    params:
        strategy = config["pipeline"].get("assembly_strategy", "individual")
    log:
        f"{OUTDIR}/logs/assembly/{ASSEMBLY_STRATEGY}/assembly_stats.log"
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/assembly_stats.py"
