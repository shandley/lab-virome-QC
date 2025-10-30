# GitHub Repository Setup Guide

This guide will help the lab set up proper GitHub settings, branch protection, and team management for collaborative development.

---

## Table of Contents

- [Repository Settings](#repository-settings)
- [Branch Protection Rules](#branch-protection-rules)
- [Team Management](#team-management)
- [GitHub Actions Setup](#github-actions-setup)
- [Project Boards](#project-boards)
- [Labels and Milestones](#labels-and-milestones)
- [Notifications](#notifications)

---

## Repository Settings

### General Settings

1. **Navigate to Settings** (must have admin access)
   - Go to: https://github.com/shandley/lab-virome-QC/settings

2. **Features to Enable:**
   ```
   ✅ Issues
   ✅ Projects
   ✅ Wiki (optional - for detailed documentation)
   ✅ Discussions (optional - for community Q&A)
   ❌ Sponsorships (probably not needed)
   ```

3. **Pull Requests Settings:**
   ```
   ✅ Allow squash merging (recommended - keeps history clean)
   ✅ Allow merge commits
   ❌ Allow rebase merging (can cause confusion for beginners)
   ✅ Automatically delete head branches (keeps branches clean)
   ✅ Allow auto-merge
   ```

4. **Default Branch:**
   - Ensure `main` is the default branch

---

## Branch Protection Rules

**Critical for preventing accidental changes to main!**

### Setting Up Branch Protection

1. **Go to Settings → Branches → Add rule**

2. **Branch name pattern:** `main`

3. **Enable these protections:**

```
✅ Require a pull request before merging
   ✅ Require approvals (set to 1)
   ✅ Dismiss stale pull request approvals when new commits are pushed
   ❌ Require review from Code Owners (optional, for larger teams)

✅ Require status checks to pass before merging
   ✅ Require branches to be up to date before merging
   Status checks to require:
      - linting
      - dry-run
      - validation

✅ Require conversation resolution before merging

✅ Require linear history (optional - enforces clean history)

❌ Require deployments to succeed (not needed for this project)

✅ Lock branch (only for releases, not main)

✅ Do not allow bypassing the above settings
   - Admins can bypass (useful for emergencies)
   - Enforce for administrators (recommended after initial setup)

✅ Restrict who can push to matching branches
   - Add: Lab members who should have write access
```

### Why These Protections Matter

| Protection | Why It Helps |
|------------|--------------|
| **Require PR** | Prevents direct commits to main; forces code review |
| **Require approvals** | Ensures someone else reviews your code |
| **Status checks** | Automated testing before merge (catches errors) |
| **Conversation resolution** | Ensures all feedback is addressed |
| **Up to date** | Prevents merge conflicts |

---

## Team Management

### Creating a GitHub Team

1. **Organization Required:**
   - Individual accounts can't create teams
   - Consider creating a GitHub Organization for the lab
   - Free for public repos: https://github.com/organizations/plan

2. **If Using Organization:**

   **Create team:**
   ```
   Organization Settings → Teams → New team

   Team name: virome-qc-team
   Description: Lab Virome QC Pipeline developers
   Visibility: Visible (lab members can see)
   ```

   **Add team members:**
   ```
   Members tab → Add member

   Roles:
   - Maintainer: PI, senior lab members (can manage team)
   - Member: Everyone else (can contribute)
   ```

3. **If Using Individual Repo (Current Setup):**

   **Add collaborators:**
   ```
   Settings → Collaborators → Add people

   Permissions:
   - Admin: PI only
   - Write: All lab members
   - Read: External collaborators (if any)
   ```

### Recommended Roles

| Person | Role | Why |
|--------|------|-----|
| **Scott Handley (PI)** | Admin | Final approval, settings management |
| **Senior lab member** | Admin or Write | Help with code review, mentorship |
| **Bioinformatician** | Write | Primary development, review PRs |
| **Lab members** | Write | Contribute features, fixes |
| **Rotation students** | Write | Contribute with oversight |
| **External collaborators** | Read | Can view, clone, but not push |

---

## GitHub Actions Setup

**Already configured in `.github/workflows/ci.yml`**

### What the CI Does

1. **Linting** - Checks code style
2. **Dry-run** - Validates Snakemake syntax
3. **Validation** - Checks YAML files, required files
4. **Documentation** - Checks markdown links

### Viewing CI Results

- Go to: Actions tab in GitHub
- See status of each workflow run
- Red ❌ = failed, Green ✅ = passed

### Required Status Checks

Configure in Branch Protection (see above):
- `linting`
- `dry-run`
- `validation`

This prevents merging if tests fail!

---

## Project Boards

**Use GitHub Projects for task management**

### Creating a Project Board

1. **Go to Projects tab → New project**

2. **Choose template:**
   - **Team backlog** (recommended)
   - **Feature** (for specific features)
   - **Bug triage** (for bug tracking)

3. **Recommended Setup:**

   **Board name:** `Lab Virome QC Development`

   **Columns:**
   ```
   📋 Backlog       - Ideas, future tasks
   🎯 To Do         - Prioritized tasks
   🏃 In Progress   - Currently working on
   👀 In Review     - PR submitted, needs review
   ✅ Done          - Completed
   ```

4. **Using the Board:**

   - Convert issues to project cards
   - Drag cards between columns
   - Assign to team members
   - Set due dates for milestones

### Example Workflow

```
User creates issue (#42: Add primer trimming)
   ↓
Moved to "To Do" column
   ↓
Lab member picks it up → "In Progress"
   ↓
Creates PR → "In Review"
   ↓
PR approved & merged → "Done"
   ↓
Issue automatically closed
```

---

## Labels and Milestones

### Essential Labels

**Create these labels** (Settings → Labels):

| Label | Color | Purpose |
|-------|-------|---------|
| `bug` | #d73a4a | Something isn't working |
| `enhancement` | #a2eeef | New feature or request |
| `documentation` | #0075ca | Documentation improvements |
| `good first issue` | #7057ff | Good for newcomers |
| `help wanted` | #008672 | Extra attention needed |
| `question` | #d876e3 | Further information requested |
| `wontfix` | #ffffff | This will not be worked on |
| `duplicate` | #cfd3d7 | Duplicate issue |
| `priority: high` | #ff0000 | Urgent |
| `priority: medium` | #ffaa00 | Normal priority |
| `priority: low` | #00ff00 | Can wait |
| `status: blocked` | #b60205 | Waiting on something |
| `status: in-progress` | #fbca04 | Currently being worked on |

### Creating Milestones

**For planning releases/phases:**

```
Settings → Milestones → New milestone

Examples:
- v1.0 - Initial Release (Target: 2025-02-01)
- Phase 1 - Basic QC Pipeline (Target: 2025-01-15)
- Phase 2 - Advanced Features (Target: 2025-03-01)
```

**Link issues to milestones:**
- Helps track progress toward goals
- Shows % completion
- Good for lab meeting updates

---

## Notifications

### Configuring Notifications

**Personal Settings:**
```
GitHub → Settings → Notifications

Recommended settings:
✅ Email notifications for: Participating and @mentions
✅ Web notifications for: Everything
✅ Watching: Automatic watching for repositories you contribute to
```

### Watch Settings for This Repo

```
Repository → Watch button (top right)

Options:
- All Activity (recommended for core team)
- Participating and @mentions (for occasional contributors)
- Ignore (not recommended)

Customize:
✅ Issues
✅ Pull Requests
✅ Releases
❌ Discussions (unless you want)
```

### Managing Notification Overload

**Tips:**
1. Use email filters (Gmail rule: `to:lab-virome-QC`)
2. Check notifications once or twice a day
3. Unsubscribe from threads you're not involved in
4. Use `@mentions` to get someone's attention
5. Set "Do Not Disturb" hours in GitHub settings

---

## Security Settings

### Required for Lab Compliance

1. **Enable Security Features:**
   ```
   Settings → Security → Code security and analysis

   ✅ Dependency graph
   ✅ Dependabot alerts
   ✅ Dependabot security updates
   ✅ Secret scanning (if available)
   ```

2. **What This Does:**
   - Alerts you to vulnerable dependencies
   - Automatically creates PRs to update dependencies
   - Prevents accidental commit of secrets (passwords, keys)

---

## Recommended Workflow for New Lab Members

### Day 1: Setup

```bash
# 1. Accept invitation to repository
# 2. Clone the repo
git clone https://github.com/shandley/lab-virome-QC.git
cd lab-virome-QC

# 3. Set up development environment
conda create -n virome-qc-dev snakemake
conda activate virome-qc-dev

# 4. Configure git
git config user.name "Your Name"
git config user.email "your.email@wustl.edu"

# 5. Read CONTRIBUTING.md
# 6. Pick a "good first issue" from the project board
```

### Weekly Routine

**Monday:**
- Sync with team (lab meeting or Slack)
- Review project board
- Pick tasks for the week

**During Week:**
- Create feature branches
- Make commits
- Push changes
- Create PRs
- Review others' PRs

**Friday:**
- Update project board
- Document blockers
- Plan next week

---

## Automation Tips

### Auto-assign Reviewers

**Create `.github/CODEOWNERS`:**
```
# Auto-request review from these users

# All files
* @shandley

# Snakemake files
workflow/* @bioinformatician

# Python scripts
*.py @python-expert

# Documentation
*.md @documentation-lead
```

### Issue Templates (Already Created!)

Located in `.github/ISSUE_TEMPLATE/`:
- `bug_report.md`
- `feature_request.md`
- `question.md`

These guide users to provide necessary information.

---

## Troubleshooting Common Issues

### "I can't push to main"

✅ **This is intentional!** Branch protection is working.

**Solution:** Create a branch and PR instead.

### "CI is failing but my code works locally"

**Check:**
1. Did you commit all files?
2. Are there YAML syntax errors?
3. Check the Actions tab for specific error

### "I accidentally committed to main"

**Don't panic!**
```bash
# If you haven't pushed:
git reset --soft HEAD~1

# If you have pushed (and have permission):
git revert HEAD
git push
```

**Then:** Create proper PR next time!

### "Merge conflicts!"

```bash
# Update your branch with main
git checkout your-branch
git fetch origin
git merge origin/main

# Resolve conflicts in files
# Then:
git add .
git commit -m "Resolve merge conflicts"
git push
```

---

## Monthly Maintenance Tasks

**PI or Maintainer should:**

- [ ] Review and close stale issues
- [ ] Update milestones
- [ ] Check security alerts
- [ ] Review team permissions
- [ ] Archive completed project boards
- [ ] Update documentation
- [ ] Review PR review turnaround time
- [ ] Check CI status and fix if needed

---

## Resources

### GitHub Guides
- [GitHub Flow](https://guides.github.com/introduction/flow/)
- [Understanding Projects](https://docs.github.com/en/issues/planning-and-tracking-with-projects)
- [About Protected Branches](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/defining-the-mergeability-of-pull-requests/about-protected-branches)

### Lab-Specific
- Slack channel: `#virome-qc`
- Lab meetings: Monday 10am
- Git/GitHub office hours: Friday 2pm

---

## Questions?

**For GitHub setup questions:**
- Open an issue with `question` label
- Ask in lab Slack
- Email: scott.handley@wustl.edu

**For permissions issues:**
- Contact: @shandley (PI)
- Must be logged into GitHub
- Check spam for invitation emails

---

## Checklist for Repository Setup

**PI/Admin should complete:**

- [ ] Enable branch protection on `main`
- [ ] Add all lab members as collaborators
- [ ] Create project board
- [ ] Set up labels
- [ ] Create first milestone
- [ ] Configure notifications
- [ ] Enable security features
- [ ] Review and approve CI workflow
- [ ] Set up CODEOWNERS (optional)
- [ ] Schedule regular maintenance

**Each lab member should:**

- [ ] Accept repository invitation
- [ ] Fork or clone repository
- [ ] Read CONTRIBUTING.md
- [ ] Set up development environment
- [ ] Configure git with name and email
- [ ] Make a test branch and PR
- [ ] Introduce yourself in an issue/PR
- [ ] Pick a "good first issue"

---

**Last updated:** 2025-01-30
**Next review:** 2025-02-28
