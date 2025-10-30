# Contributing to Lab Virome QC

Welcome! This guide will help you contribute effectively to the Lab Virome QC pipeline. Whether you're fixing a bug, adding a feature, or improving documentation, we appreciate your help.

---

## Table of Contents

- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Coding Standards](#coding-standards)
- [Testing Guidelines](#testing-guidelines)
- [Documentation](#documentation)
- [Submitting Changes](#submitting-changes)
- [Communication](#communication)

---

## Getting Started

### 1. Fork and Clone

```bash
# Fork the repository on GitHub (click "Fork" button)

# Clone your fork
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
cd lab-virome-QC

# Add upstream remote
git remote add upstream https://github.com/shandley/lab-virome-QC.git

# Verify remotes
git remote -v
```

### 2. Set Up Development Environment

```bash
# Create conda environment with Snakemake
conda create -n virome-qc-dev -c conda-forge -c bioconda snakemake python=3.10
conda activate virome-qc-dev

# Install development dependencies
conda install -c conda-forge pre-commit black flake8 pytest
```

### 3. Create a Branch

**Branch Naming Convention:**
- `feature/description` - New features (e.g., `feature/add-primer-trimming`)
- `fix/description` - Bug fixes (e.g., `fix/viromeqc-parsing`)
- `docs/description` - Documentation updates (e.g., `docs/update-readme`)
- `refactor/description` - Code refactoring (e.g., `refactor/qc-rules`)

```bash
# Create and switch to new branch
git checkout -b feature/your-feature-name

# OR for bug fix
git checkout -b fix/bug-description
```

---

## Development Workflow

### The GitHub Flow (Simplified)

```
main (protected)
  ↓
  ├─ feature/new-qc-rule → Pull Request → Review → Merge
  ├─ fix/bug-in-fastp → Pull Request → Review → Merge
  └─ docs/update-guide → Pull Request → Review → Merge
```

### Daily Workflow

1. **Start your day: Sync with upstream**
```bash
git checkout main
git fetch upstream
git merge upstream/main
git push origin main
```

2. **Work on your branch**
```bash
git checkout feature/your-feature
# Make changes...
git add .
git commit -m "Clear description of changes"
```

3. **Keep your branch updated**
```bash
# Periodically merge main into your feature branch
git checkout main
git pull upstream main
git checkout feature/your-feature
git merge main
```

4. **Push your changes**
```bash
git push origin feature/your-feature
```

5. **Create Pull Request on GitHub**

---

## Coding Standards

### Snakemake Rules

**Good Rule Example:**
```python
rule fastp_trim:
    """
    Clear description of what this rule does

    Key parameters:
    - Quality threshold: Q20
    - Minimum length: 100bp
    """
    input:
        r1 = "path/to/input_R1.fastq.gz",
        r2 = "path/to/input_R2.fastq.gz"
    output:
        r1 = "path/to/output_R1.fastq.gz",
        r2 = "path/to/output_R2.fastq.gz"
    log:
        "logs/fastp/{sample}.log"
    threads: 8
    conda:
        "envs/qc.yaml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --thread {threads} \
            2>&1 | tee {log}
        """
```

**Best Practices:**
- ✅ Always include docstrings for rules
- ✅ Use `log:` directive for all rules
- ✅ Specify `threads:` and `resources:` appropriately
- ✅ Use conda environments (don't hardcode paths to tools)
- ✅ Use `temp()` for intermediate files when appropriate
- ✅ Include meaningful variable names

### Python Scripts

Follow PEP 8 style guide:

```python
# Good
def parse_viromeqc_output(input_file):
    """
    Parse ViromeQC output file to extract enrichment score.

    Args:
        input_file (str): Path to ViromeQC output file

    Returns:
        float: Enrichment score or None if not found
    """
    with open(input_file) as f:
        for line in f:
            if "enrichment" in line.lower():
                return float(line.split()[-1])
    return None
```

**Use type hints when possible:**
```python
from pathlib import Path
from typing import Optional

def parse_file(filepath: Path) -> Optional[float]:
    """Parse file and return score."""
    pass
```

### Configuration Files

**YAML Formatting:**
```yaml
# Good: Clear structure with comments
qc_thresholds:
  # ViromeQC enrichment score threshold
  min_enrichment_score: 10

  # Maximum acceptable host contamination (%)
  max_host_percent: 10
```

---

## Testing Guidelines

### Before Submitting a PR

**1. Test Your Changes Locally**

```bash
# Dry run to check syntax
snakemake --use-conda -n

# Run on test data (if available)
snakemake --use-conda --cores 4

# Check specific rule
snakemake --use-conda -n results/fastp/test_sample_R1.fastq.gz
```

**2. Test Different Scenarios**

- ✅ Single sample
- ✅ Multiple samples
- ✅ Different file paths
- ✅ Edge cases (small files, missing inputs)

**3. Lint Your Code**

```bash
# Python scripts
flake8 workflow/scripts/

# Snakemake
snakemake --lint
```

### Writing Tests

For Python scripts, add tests in `tests/`:

```python
# tests/test_qc_flags.py
def test_parse_viromeqc():
    """Test ViromeQC parsing."""
    result = parse_viromeqc("test_data/sample1_viromeqc.txt")
    assert result is not None
    assert result > 0
```

Run tests:
```bash
pytest tests/
```

---

## Documentation

### What to Document

**1. Code Changes**
- Add docstrings to new rules and functions
- Update inline comments for complex logic
- Include parameter descriptions

**2. User-Facing Changes**
- Update README.md if workflow changes
- Update config/config.yaml comments
- Add examples for new features

**3. Developer Notes**
- Update CONTRIBUTING.md for new processes
- Document dependencies in environment files
- Add troubleshooting tips

### Documentation Style

**README Updates:**
```markdown
## New Feature Name

Brief description of what it does.

### Usage

\`\`\`bash
# Example command
snakemake --use-conda --config new_param=value
\`\`\`

### Parameters

- `param_name`: Description of parameter
```

---

## Submitting Changes

### Pull Request Process

**1. Prepare Your PR**

```bash
# Make sure you're on your feature branch
git checkout feature/your-feature

# Ensure it's up to date with main
git fetch upstream
git merge upstream/main

# Push to your fork
git push origin feature/your-feature
```

**2. Create PR on GitHub**

- Go to https://github.com/shandley/lab-virome-QC
- Click "Pull requests" → "New pull request"
- Click "compare across forks"
- Select your fork and branch
- Fill out the PR template (see below)

**3. PR Title Format**

```
[Type] Brief description

Examples:
[Feature] Add primer A trimming rule
[Fix] Correct ViromeQC parsing logic
[Docs] Update installation instructions
[Refactor] Simplify read counting script
```

**4. PR Description Should Include:**

- **What** changed
- **Why** the change was needed
- **How** to test it
- **Related issues** (if any): Fixes #123

**5. Request Review**

Tag 1-2 lab members as reviewers

### Code Review Guidelines

**As a Reviewer:**

✅ **Look for:**
- Code works as intended
- Follows coding standards
- Includes appropriate tests
- Documentation is clear
- No breaking changes (unless intended)

✅ **Provide constructive feedback:**
```
Good: "Consider using `Path` instead of string for file paths.
       This handles OS differences better."

Bad: "This is wrong."
```

✅ **Approve or Request Changes**

**As an Author:**

✅ Respond to all comments
✅ Make requested changes
✅ Re-request review after updates
✅ Don't take feedback personally - we're all learning!

### Merging

**Who Can Merge:**
- Repository maintainers only
- Requires 1 approval
- All CI checks must pass

**How to Merge:**
1. Squash and merge (for feature branches)
2. Delete branch after merge
3. Pull updated main to your local

```bash
git checkout main
git pull upstream main
```

---

## Communication

### Where to Communicate

**GitHub Issues:**
- Bug reports
- Feature requests
- General questions
- Task tracking

**Pull Requests:**
- Code review discussions
- Implementation details

**Lab Meetings:**
- Major design decisions
- Roadmap planning
- Training sessions

### Issue Guidelines

**Before Creating an Issue:**
1. Search existing issues
2. Check if it's already fixed in `main`
3. Gather relevant information

**Good Issue Example:**

```markdown
**Bug Report: ViromeQC fails with error**

**Environment:**
- OS: Ubuntu 22.04
- Snakemake version: 7.32.4
- ViromeQC version: 1.0

**Steps to Reproduce:**
1. Run pipeline with config X
2. See error at ViromeQC step

**Expected behavior:**
ViromeQC should complete successfully

**Actual behavior:**
Error: "KeyError: enrichment"

**Logs:**
[Attach log file]

**Possible solution:**
Maybe need to update parsing logic?
```

### Commit Message Guidelines

**Format:**
```
<type>: <subject>

<body>

<footer>
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation only
- `style`: Formatting (no code change)
- `refactor`: Code change (no functional change)
- `test`: Adding tests
- `chore`: Maintenance

**Examples:**

```
feat: add random primer A trimming rule

- Implements hard trimming of 6bp from 5' end
- Adds configuration parameter for primer length
- Updates documentation

Closes #42
```

```
fix: correct ViromeQC enrichment score parsing

The regex was failing on scientific notation. Now handles
both decimal and scientific notation.

Fixes #38
```

---

## Lab-Specific Practices

### Weekly Sync

- **When:** Every Monday 10am
- **What:** Review PRs, discuss blockers, assign tasks
- **Where:** Lab meeting room / Zoom

### Pair Programming

Encouraged for:
- Complex features
- Learning new tools
- Debugging tricky issues

Schedule via Slack/email

### Learning Resources

**Snakemake:**
- Official tutorial: https://snakemake.readthedocs.io/
- Lab's tutorial videos: [Add link]

**Git/GitHub:**
- GitHub Learning Lab: https://lab.github.com/
- Git branching game: https://learngitbranching.js.org/

**Python:**
- Real Python: https://realpython.com/
- Lab's Python workshop notes: [Add link]

### Mentorship

**New to the project?** Pair with experienced member:
- Lab member A: Snakemake expert
- Lab member B: ViromeQC specialist
- Lab member C: Git/GitHub expert

---

## Getting Help

**Stuck? Here's what to do:**

1. **Check Documentation**
   - README.md
   - Issue discussions
   - Snakemake docs

2. **Ask on GitHub**
   - Open an issue with "Question" label
   - Tag relevant person

3. **Ask in Lab Slack**
   - #virome-qc channel
   - Quick questions

4. **Schedule 1-on-1**
   - For complex issues
   - Pair programming sessions

---

## Recognition

All contributors will be:
- Listed in README.md
- Acknowledged in publications using the pipeline
- Credited in releases

Thank you for contributing! 🎉

---

## Quick Reference

```bash
# Setup
git clone https://github.com/YOUR_USERNAME/lab-virome-QC.git
git remote add upstream https://github.com/shandley/lab-virome-QC.git

# Daily workflow
git checkout main
git pull upstream main
git checkout -b feature/my-feature
# ... make changes ...
git add .
git commit -m "feat: description"
git push origin feature/my-feature
# Create PR on GitHub

# Sync with main
git fetch upstream
git merge upstream/main

# Testing
snakemake --use-conda -n
snakemake --use-conda --cores 4

# Clean up after merge
git checkout main
git pull upstream main
git branch -d feature/my-feature
```

---

## Questions?

Open an issue or contact:
- Scott Handley (PI): scott.handley@wustl.edu
- [Lab member]: [email]
- [Lab member]: [email]
