# Lab Team Onboarding Guide

**Welcome to the Lab Virome QC project!**

This quick guide will get you up and running with the project and collaborative GitHub workflow.

---

## 🎯 Goals of This Project

1. **Build a robust virome QC pipeline** for our VLP-enriched sequencing data
2. **Learn collaborative software development** using Git and GitHub
3. **Create reusable, documented tools** for the lab and community
4. **Practice modern bioinformatics** best practices

---

## 📋 First Week Checklist

### Day 1: Account Setup (30 minutes)

- [ ] **Create GitHub account** (if you don't have one)
  - Go to: https://github.com/signup
  - Use your institutional email (e.g., @wustl.edu)

- [ ] **Accept repository invitation**
  - Check your email for invitation to `lab-virome-QC`
  - Accept the invitation

- [ ] **Set up two-factor authentication** (2FA)
  - GitHub Settings → Security → Two-factor authentication
  - Use an authenticator app (Google Authenticator, Authy)
  - **Important for lab security!**

- [ ] **Configure your profile**
  - Add a profile picture (helps team recognize you)
  - Add your name
  - Optional: Add bio, location

### Day 2: Git Setup (1 hour)

- [ ] **Install Git**
  ```bash
  # macOS (using Homebrew)
  brew install git

  # Linux (Ubuntu/Debian)
  sudo apt-get install git

  # Windows
  # Download from: https://git-scm.com/download/win
  ```

- [ ] **Configure Git**
  ```bash
  git config --global user.name "Your Name"
  git config --global user.email "your.email@wustl.edu"

  # Check configuration
  git config --list
  ```

- [ ] **Set up SSH keys** (recommended for easier pushing)
  ```bash
  # Generate SSH key
  ssh-keygen -t ed25519 -C "your.email@wustl.edu"

  # Copy public key
  cat ~/.ssh/id_ed25519.pub

  # Add to GitHub: Settings → SSH and GPG keys → New SSH key
  ```

### Day 3: Project Setup (1-2 hours)

- [ ] **Clone the repository**
  ```bash
  cd ~/Code  # or wherever you keep projects
  git clone https://github.com/shandley/lab-virome-QC.git
  cd lab-virome-QC
  ```

- [ ] **Set up development environment**
  ```bash
  # Install Conda/Mamba if not already installed
  # Then create Snakemake environment:
  conda create -n virome-qc-dev -c conda-forge -c bioconda snakemake python=3.10
  conda activate virome-qc-dev
  ```

- [ ] **Explore the codebase**
  ```bash
  # Look at project structure
  ls -la

  # Read the README
  cat README.md

  # Look at the Snakefile
  cat workflow/Snakefile
  ```

### Day 4: Learn the Workflow (2-3 hours)

- [ ] **Read documentation**
  - [ ] README.md (pipeline overview)
  - [ ] CONTRIBUTING.md (how to contribute)
  - [ ] CODE_OF_CONDUCT.md (community guidelines)

- [ ] **Complete Git tutorial**
  - Start with: https://lab.github.com/
  - Course: "Introduction to GitHub"
  - Takes ~30 minutes

- [ ] **Play Git branching game** (fun way to learn!)
  - https://learngitbranching.js.org/
  - Complete first 4 levels

### Day 5: Make Your First Contribution (1-2 hours)

- [ ] **Find a "good first issue"**
  - Go to: Issues tab → Filter by `good first issue` label
  - Or: Start with documentation improvements (fix typos, clarify instructions)

- [ ] **Create a branch**
  ```bash
  git checkout -b fix/my-first-fix
  ```

- [ ] **Make your change**
  - Even fixing a typo is great practice!
  - Add a comment to code
  - Improve documentation

- [ ] **Commit and push**
  ```bash
  git add .
  git commit -m "docs: fix typo in README"
  git push origin fix/my-first-fix
  ```

- [ ] **Create Pull Request**
  - Go to GitHub repository
  - Click "Pull requests" → "New pull request"
  - Fill out the template
  - Request review from a teammate

🎉 **Congratulations!** You've completed your first contribution!

---

## 🔄 Daily Workflow (Once You're Up to Speed)

### Morning Routine (5 minutes)

```bash
# 1. Start your day by syncing with main
git checkout main
git pull origin main

# 2. Check what you were working on
git branch
```

### Working on a Task

```bash
# 1. Create a branch for your work
git checkout -b feature/add-something

# 2. Make changes to files
# (edit files in your editor)

# 3. Test your changes
snakemake --use-conda -n  # dry run

# 4. Commit your work (frequently!)
git add .
git commit -m "feat: add new QC rule for X"

# 5. Push to your branch
git push origin feature/add-something
```

### Creating a Pull Request

1. Go to GitHub repository
2. Click "Pull requests" → "New pull request"
3. Fill out the template:
   - What did you change?
   - Why did you change it?
   - How did you test it?
4. Request review from 1-2 teammates
5. Address feedback
6. Once approved, it gets merged!

---

## 🤝 Team Communication

### Where to Communicate

| Topic | Where | When |
|-------|-------|------|
| **Bug reports** | GitHub Issues | Anytime |
| **Feature ideas** | GitHub Issues | Anytime |
| **Code questions** | GitHub Issues | Anytime |
| **Quick questions** | Lab Slack #virome-qc | Daily |
| **Design discussions** | Lab meetings | Monday 10am |
| **Git/GitHub help** | Office hours | Friday 2pm |

### How to Get Help

**Stuck on something?**

1. **Check documentation first**
   - README.md
   - LEARNING_RESOURCES.md
   - Existing issues

2. **Search Google/Stack Overflow**
   - Most Git/GitHub questions have been answered

3. **Ask on Slack**
   - #virome-qc channel
   - "@mention" person if urgent

4. **Create an issue**
   - Use "Question" template
   - Provide context and what you've tried

5. **Schedule pair programming**
   - Great for complex problems
   - Learn faster by working together

---

## 📚 Essential Reading

**Week 1:**
- [ ] README.md (this project)
- [ ] CONTRIBUTING.md (this project)
- [ ] GitHub Flow: https://guides.github.com/introduction/flow/

**Week 2:**
- [ ] Snakemake Tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
- [ ] Git Handbook: https://guides.github.com/introduction/git-handbook/

**Week 3:**
- [ ] ViromeQC paper: Zolfo et al. (2019) Nature Biotech
- [ ] Our key papers list: [Add link to lab papers]

**Ongoing:**
- Explore LEARNING_RESOURCES.md for tutorials on specific topics

---

## 🎯 Expectations & Best Practices

### What We Expect

**From Everyone:**
- Respectful communication (see CODE_OF_CONDUCT.md)
- Ask questions when stuck
- Document your work
- Test your code before submitting
- Respond to code review feedback
- Help onboard new members

**From Contributors:**
- One PR at a time (keeps things manageable)
- Commit messages follow format (see CONTRIBUTING.md)
- Code is tested locally before pushing
- Documentation updated for new features

### What You Can Expect

**From the Team:**
- Welcoming environment for beginners
- Patient code reviews with explanations
- Help when you're stuck
- Recognition for your contributions

**From Leadership:**
- Clear project direction
- Regular meetings and check-ins
- Resources for learning
- Authorship on publications

---

## 🚀 Skill Development Plan

### Beginner Level (Weeks 1-4)

**Goals:**
- Comfortable with Git basics (clone, branch, commit, push)
- Can navigate GitHub (issues, PRs)
- Understand pipeline structure
- Can make simple documentation changes

**Milestones:**
- [ ] Make 3 successful PRs
- [ ] Review someone else's PR
- [ ] Fix a simple bug

### Intermediate Level (Months 2-3)

**Goals:**
- Can write Snakemake rules
- Understand Python scripts in pipeline
- Can debug issues independently
- Comfortable with Git branching strategies

**Milestones:**
- [ ] Add a new QC rule
- [ ] Write a Python script
- [ ] Help onboard a new team member

### Advanced Level (Months 4+)

**Goals:**
- Lead feature implementation
- Mentor new contributors
- Make architectural decisions
- Present work in lab meetings

**Milestones:**
- [ ] Design and implement major feature
- [ ] Present pipeline at lab meeting
- [ ] Co-author on manuscript

---

## ⚠️ Common Pitfalls to Avoid

### For Git/GitHub Beginners

❌ **Don't:**
- Commit directly to `main` (it's protected anyway)
- Make very large PRs (hard to review)
- Forget to pull before starting work
- Be afraid to ask questions

✅ **Do:**
- Create branches for all changes
- Make small, focused PRs
- Sync with main regularly
- Ask for help early

### For Code Development

❌ **Don't:**
- Push untested code
- Hardcode file paths
- Skip documentation
- Ignore error messages

✅ **Do:**
- Test locally first (dry run, full run)
- Use config files for paths
- Document as you go
- Read error messages carefully

---

## 🏆 Success Metrics

**After 1 Month:**
- Made 5+ contributions (even small ones)
- Comfortable with basic Git workflow
- Understand pipeline structure
- Know who to ask for help

**After 3 Months:**
- Implemented a feature or fixed bugs
- Can review others' code
- Helping onboard new members
- Comfortable with most tools

**After 6 Months:**
- Leading feature development
- Mentoring newer members
- Contributing to project direction
- Ready for co-authorship

---

## 📞 Who to Contact

| Question | Contact |
|----------|---------|
| **Repository access** | Scott Handley (@shandley) |
| **Git/GitHub help** | [Git expert name] |
| **Snakemake questions** | [Snakemake expert name] |
| **Python help** | [Python expert name] |
| **Virome biology** | [Biology expert name] |
| **General questions** | Slack #virome-qc |

---

## 🎊 Welcome to the Team!

We're excited to have you working on this project. Remember:

- **Everyone was a beginner once** - don't be intimidated
- **Mistakes are learning opportunities** - Git makes them reversible!
- **Questions are encouraged** - help us improve documentation
- **Your contributions matter** - even small improvements help
- **Have fun!** - Science should be enjoyable

Looking forward to your contributions! 🧬🚀

---

## Quick Links

- **Repository:** https://github.com/shandley/lab-virome-QC
- **Project Board:** [Add link after creating]
- **Lab Slack:** #virome-qc
- **Meeting Notes:** [Add link to shared doc]
- **Lab Calendar:** [Add link]

---

**Questions about onboarding?**
- Open an issue with "Question" label
- Ask in Slack #virome-qc
- Email: scott.handley@wustl.edu

**Ready to get started? Pick up a "good first issue" and dive in!**
