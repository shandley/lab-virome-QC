# Learning Resources for Lab Virome QC Development

A curated list of resources to help lab members learn the tools and practices we use.

---

## Quick Start Guides

### Absolute Beginners

**Never used Git/GitHub before? Start here:**

1. **GitHub Learning Lab** (Interactive, ~2 hours)
   - https://lab.github.com/
   - Courses: "Introduction to GitHub", "Communicating using Markdown"

2. **Git Branching Game** (Learn by playing, ~1 hour)
   - https://learngitbranching.js.org/
   - Makes Git concepts visual and fun

3. **GitHub Flow** (15 minutes reading)
   - https://guides.github.com/introduction/flow/
   - Our workflow in a nutshell

### Command Line Basics

**Need command line refresher?**

- **Unix/Linux Tutorial** (2-3 hours)
  - https://swcarpentry.github.io/shell-novice/
  - From Software Carpentry (excellent for scientists)

---

## Git & GitHub

### Beginner Level

**Concepts to Master:**
- Commit, push, pull
- Branches
- Pull requests
- Merge conflicts

**Resources:**

| Resource | Time | Why It's Good |
|----------|------|---------------|
| [Git Handbook](https://guides.github.com/introduction/git-handbook/) | 30 min | Quick overview |
| [GitHub Skills](https://skills.github.com/) | 2-4 hours | Interactive courses |
| [Oh My Git!](https://ohmygit.org/) | 1-2 hours | Learn Git through a game |
| [Atlassian Git Tutorial](https://www.atlassian.com/git/tutorials) | Varies | Comprehensive reference |

### Intermediate Level

**Concepts to Master:**
- Rebasing
- Cherry-picking
- Resolving complex conflicts
- Git hooks

**Resources:**
- **Pro Git Book** (free)
  - https://git-scm.com/book/en/v2
  - Chapters 2, 3, 5 are most relevant

- **Git Katas** (exercises)
  - https://github.com/eficode-academy/git-katas
  - Practice specific Git skills

### Cheat Sheets

- **GitHub Git Cheat Sheet** (PDF)
  - https://education.github.com/git-cheat-sheet-education.pdf

- **Atlassian Git Cheat Sheet**
  - https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet

---

## Snakemake

### Beginner Level

**Start Here:**

1. **Official Snakemake Tutorial** (4-6 hours)
   - https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
   - Best introduction to Snakemake concepts

2. **Snakemake Short Tutorial** (HBC Training)
   - https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/04_snakemake.html
   - Bioinformatics-focused

**Concepts to Master:**
- Rules and DAGs
- Wildcards
- Input functions
- Conda integration

### Intermediate Level

**Advanced Topics:**

1. **Best Practices** (Official docs)
   - https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html

2. **Cluster Execution**
   - https://snakemake.readthedocs.io/en/stable/executing/cluster.html

3. **Configuration Files**
   - https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html

### Video Tutorials

- **Snakemake for Beginners** (YouTube, Johannes Köster)
  - Creator of Snakemake explains concepts

- **Bioinformatics Workflow Management with Snakemake**
  - Various conference talks on YouTube

### Example Workflows

**Learn by reading:**
- **Snakemake Wrappers**
  - https://snakemake-wrappers.readthedocs.io/
  - Pre-built rules for common tools

- **nf-core (Nextflow, but similar concepts)**
  - https://nf-co.re/
  - Well-documented pipelines

---

## Python for Bioinformatics

### Beginner Level

**If you know some Python but not for bioinformatics:**

1. **Python for Biologists**
   - https://pythonforbiologists.com/
   - Free tutorials specifically for biologists

2. **Rosalind** (Problem-solving platform)
   - https://rosalind.info/problems/list-view/
   - Learn by solving bioinformatics problems

### Our Specific Needs

**What you'll use Python for in this project:**
- Parsing output files (text processing)
- Calculating QC metrics
- Generating reports

**Key Libraries:**
- `pandas` - Data manipulation
- `pathlib` - File path handling
- `argparse` - Command-line arguments

**Quick Tutorials:**
- **Pandas 10 minutes**
  - https://pandas.pydata.org/docs/user_guide/10min.html

- **Pathlib Tutorial**
  - https://realpython.com/python-pathlib/

---

## YAML

**Configuration file format we use**

### Quick Learn (15 minutes)

- **Learn YAML in Y minutes**
  - https://learnxinyminutes.com/docs/yaml/

- **YAML Syntax**
  - https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html

### Tools

- **YAML Validator**
  - http://www.yamllint.com/
  - Check syntax online

---

## Conda & Bioconda

**Package management for bioinformatics**

### Beginner Level

1. **Conda User Guide**
   - https://docs.conda.io/projects/conda/en/latest/user-guide/index.html
   - Start with "Getting Started"

2. **Bioconda Tutorial**
   - https://bioconda.github.io/
   - Bioinformatics-specific packages

**Key Commands:**
```bash
# Create environment
conda create -n myenv python=3.10

# Activate
conda activate myenv

# Install packages
conda install -c bioconda snakemake

# Export environment
conda env export > environment.yaml
```

---

## Virome Analysis & QC

### Background Reading

**Understanding VLP and Viromes:**

1. **"Detecting contamination in viromes using ViromeQC"**
   - Zolfo et al. (2019) Nature Biotechnology
   - THE paper on virome QC

2. **"Evaluation of methods to purify VLPs"**
   - Conceição-Neto et al. (2015) BMC Genomics
   - VLP preparation methods

3. **"The long and short of it: benchmarking viromics"**
   - Warwick-Dugdale et al. (2024)
   - Sequencing technology comparison

### NovaSeq Artifacts

- **fastp documentation** (polyG trimming)
  - https://github.com/OpenGene/fastp

- **Illumina sequencing chemistry**
  - Official Illumina documentation

---

## Code Review

### How to Review Code

**Resources:**
- **Google's Code Review Guide**
  - https://google.github.io/eng-practices/review/
  - Industry standard practices

- **How to Make Your Code Reviewer Fall in Love**
  - https://mtlynch.io/code-review-love/
  - From the author's perspective

### What to Look For

**Checklist:**
- [ ] Does it work? (test it!)
- [ ] Is it readable?
- [ ] Are there comments for complex parts?
- [ ] Does it follow our style guide?
- [ ] Is documentation updated?
- [ ] Are there tests?

---

## Writing Good Documentation

**Resources:**

1. **Write the Docs**
   - https://www.writethedocs.org/guide/
   - Community for documentation

2. **Technical Writing Course** (Google, free)
   - https://developers.google.com/tech-writing
   - Learn technical writing basics

3. **Markdown Guide**
   - https://www.markdownguide.org/
   - We use Markdown for all docs

---

## HPC / Cluster Computing

**For running pipelines on compute clusters:**

1. **SLURM Basics**
   - https://slurm.schedmd.com/quickstart.html
   - If your HPC uses SLURM

2. **Snakemake Cluster Execution**
   - https://snakemake.readthedocs.io/en/stable/executing/cluster.html

3. **WashU RIS Documentation** (if applicable)
   - Internal documentation for your cluster

---

## Bioinformatics Tools Used

### BBTools
- **BBMap Guide**
  - https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/
  - Comprehensive guide to BBDuk, Clumpify, etc.

### FastQC & MultiQC
- **FastQC Tutorial**
  - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

- **MultiQC Docs**
  - https://multiqc.info/docs/

### ViromeQC
- **GitHub + Paper**
  - https://github.com/SegataLab/viromeqc

---

## Practice Datasets

**Where to get test data:**

1. **SRA (Sequence Read Archive)**
   - https://www.ncbi.nlm.nih.gov/sra
   - Search for "virome" datasets

2. **ENA (European Nucleotide Archive)**
   - https://www.ebi.ac.uk/ena/browser/home

3. **Simulated Data**
   - Create small test datasets for development

**Lab Test Data:**
- Located at: `[Add path to lab shared test data]`
- Use for testing pipeline changes

---

## Lab-Specific Resources

### Internal Documentation

- **Lab Wiki:** [Add link]
- **Lab Protocols:** [Add link]
- **Compute Resources:** [Add link]
- **Data Storage:** [Add link]

### Lab Meetings & Training

**Weekly Schedule:**
- **Monday 10am:** Lab meeting
- **Friday 2pm:** Git/GitHub office hours
- **Ad-hoc:** Pair programming sessions

**Monthly Workshops:**
- First Friday: Bioinformatics tools
- Third Friday: Data analysis & visualization

### Who to Ask

**Expertise Map:**

| Topic | Expert | Contact |
|-------|--------|---------|
| Git/GitHub | [Name] | [Email] |
| Snakemake | [Name] | [Email] |
| Python | [Name] | [Email] |
| Virome biology | [Name] | [Email] |
| HPC/Cluster | [Name] | [Email] |
| Statistics | [Name] | [Email] |

---

## Recommended Learning Path

### Week 1: Git & GitHub Basics
- [ ] Complete GitHub Learning Lab "Introduction to GitHub"
- [ ] Play Git branching game
- [ ] Make your first PR to this repo (even just fixing a typo!)
- [ ] Read CONTRIBUTING.md

### Week 2: Snakemake Basics
- [ ] Complete official Snakemake tutorial
- [ ] Read through our Snakefile with comments
- [ ] Run the pipeline locally
- [ ] Understand the DAG

### Week 3: Python for Scripts
- [ ] Review our Python scripts
- [ ] Understand pandas basics
- [ ] Modify a script (add a comment, improve output)
- [ ] Write a small test

### Week 4: Virome QC Concepts
- [ ] Read key papers (see above)
- [ ] Understand VLP prep
- [ ] Learn about NovaSeq artifacts
- [ ] Interpret MultiQC reports

### Week 5+: Contribute!
- [ ] Pick a "good first issue"
- [ ] Implement a small feature
- [ ] Go through full PR process
- [ ] Help review someone else's PR

---

## Study Groups

**Want to learn together?**

Consider forming study groups:
- **Git Study Group** (for Git beginners)
- **Snakemake Study Group** (work through tutorial together)
- **Paper Discussion** (read virome papers monthly)

Schedule via lab Slack!

---

## Staying Current

### Blogs & Newsletters

- **Snakemake Blog**
  - https://snakemake.github.io/blog/

- **GitHub Blog**
  - https://github.blog/

- **Bioinformatics Training**
  - https://bioinformaticstraining.org/

### Twitter/X Follows

- @johanneskoester (Snakemake creator)
- @github (GitHub updates)
- #bioinformatics hashtag

### Conferences & Workshops

- **BOSC** (Bioinformatics Open Source Conference)
- **ISMB** (Intelligent Systems for Molecular Biology)
- Local bioinformatics meetups

---

## Contributing to This Document

Found a great resource? **Add it!**

1. Create branch
2. Edit this file
3. Submit PR with description of resource

---

## Questions?

- Open an issue with `question` label
- Ask in lab Slack #learning channel
- Email: [Lab email]

**Remember:** Everyone was a beginner once. Don't hesitate to ask questions! 🚀
