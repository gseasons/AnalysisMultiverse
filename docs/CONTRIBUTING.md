## About the Contribution Guidelines
The Analysis Multiverse software (AMS) and its documentation should assume no computer science or programming background.

The AMS is intended to process functional magnetic resonance imaging (fMRI) data. The anticipated userbase for the AMS will be researchers from the fields of psychology, neuropsychology, neuroscience, and biology. These researchers may include student researchers. With this in mind, code should be as self-explanatory as possible. The documentation should assume no previous experience
with running or installing software from GitHub, it should assume no previous experience with high performance computing (HPC), and should not assume any fluency with Python or any other programming languages. In essence, the documentation should be as accessible as possible, to facilitate the use of AMS by all potential users.

These guidelines, as well as the documentation itself, are living documents. This means that they are meant to be updated and amended as necessary, dictated by the needs of the AMS contributors as well as the needs of the AMS users.

## How to Contribute to the Repository
> [!IMPORTANT]
> This repository enforces branch-naming rules using [Conventional Branch](https://conventionalbranch.org), and commit message formatting using [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).
>
> Using naming conventions like these gives us consistent naming practices across users and across time, allowing anyone to refer to the conventions to easily interpret documented changes and the current status of the repository. 

Git is a version control system, while GitHub is an online repository that _uses_ Git. If you are not familiar with what Git is, or with using it locally on your computer, check out the resources available from [Git](https://git-scm.com/learn) or [Github](https://docs.github.com/en/get-started/git-basics). If you are new to GitHub, please check out some [tutorials](https://www.youtube.com/watch?v=r8jQ9hVA2qs) or read the [GitHub docs](https://docs.github.com/en/get-started/start-your-journey/about-github-and-git) before proceeding. Tips for specific tutorials will be provided at the beginning of each section.

### Issues[^1]
Issues are meant to draw attention to areas that need work or have potential for growth. Issues are a feature of GitHub.

> [!TIP]
> If you're new to GitHub, learn all about [issues](https://docs.github.com/en/issues).
>
> Specifically, check out [what an issue is](https://docs.github.com/en/issues/tracking-your-work-with-issues/learning-about-issues/about-issues) and [how to create an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/creating-an-issue).

This section lists the suggested use for issues with this repository, but please keep in mind that the codebase is still being developed and tested, and the documentation is a work in progress. The issue you highlight may already be over on our [TODO](docs/TODO.md) list, so please check that first.

Use issues to...

 - Document and report bugs in the code
 - Identify functions, processes, or code files that lack adequate documentation within the code, or on our documentation website
 - Suggest clarifications for a specific section of the documentation, or in general
 - Suggest content that should be included in the documentation
 - Suggest features to include in future versions

### Pull Requests
Pull requests (PRs) are a way to include changes you've made into the main repository. PRs are a feature of GitHub.

> [!TIP]
> If you are unfamiliar with the process of committing your changes with Git, please read the [Commits](#commits) section below ***first***.
> 
> If you're new to Git and GitHub, learn all about [forks](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks), [branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches), and [pull requests](https://docs.github.com/en/pull-requests).
> 
> Specifically, check out [how to fork a repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo), [how to clone a repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), [how to create a branch](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-and-deleting-branches-within-your-repository), [what a pull request is](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests), and [how to create a pull request from a forked repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).

The process for submitting a PR to this repository is:

1. Fork the AnalysisMultiverse repository to your own GitHub account.
2. Clone the forked repository to your local machine.
3. Create a branch using [Conventional Branch](https://conventionalbranch.org) formatting, adding your username after a dash to the end of the branch name. For example, if user "emazerol" was developing a documentation update, she would name her branch "chore/update-docs-emazerol". If user "k8redfern" was developing a fix for a bug in the code, she would name her branch "fix/bug-in-the-code-k8redfern".
4. Make the changes to the files and / or code.
5. Push the changes to the forked repository on your GitHub account.
6. Go to the original AnalysisMultiverse repository, and click on the **Compare & pull request** button.

A template for your PR will be provided when you go to submit the PR. Simply replace the text under each provided header. If a specific header is not applicable, replace the template text with "N/A".

This information will be used to create the commit message for your PR if it is accepted, so please be thorough in your documentation. Improper documentation of your changes will result in a request for changes until the documentation is adequate.

### Commits
A commit is essentially a "snapshot" of the Git repository that can show the changes have been made since the last commit. Commits are saved to your Git repository with a message that you write to describe the changes made. Commits are a feature of Git.

> [!TIP]
> If you're new to Git, learn all about [commits](https://github.com/git-guides/git-commit).

Why do good commit messages matter?

Commit messages make up a permanent log of changes to a Git repository. Depending on the repository, the timeframe, and the number of contributors, the commit history can vary widely. As said by [Chris Beams](https://chris.beams.io/git-commit):
> \[...\] a well-crafted Git commit message is the best way to communicate context about a change to fellow developers (and indeed to their future selves). A diff[^2] will tell you ***what*** changed, but only the commit message can properly tell you ***why***.

Beams expands further on what will be our main point here:

> A project’s long-term success rests (among other things) on its maintainability, and a maintainer has few tools more powerful than [their] project’s log. It’s worth taking the time to learn how to care for one properly. What may be a hassle at first soon becomes habit, and eventually a source of pride and productivity for all involved.

Commit messages are a form of data and as such, part of the data management plan for AMS includes the format of commit messages as a part of our *commit*ment (pun intended) to the [FAIR principles](https://www.go-fair.org/fair-principles/) of Accessibility and Interoperability. In this spirit, we will use the advice from Beams' article "[How to Write a Git Commit Message](https://chris.beams.io/git-commit)" with the formatting structure defined by [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).

The required format will be as follows, with further details below:

| Element | Format |
|-----------|-------|
| SUBJECT | \<type\>\[optional scope\]: \<description\> |
| EMPTY LINE | empty |
| CO-AUTHORS | \[optional\] |
| EMPTY LINE | empty |
| BODY | \[optional\] |
| EMPTY LINE | empty |
| FOOTER | \[optional\] |

#### In General

1. The underlying intent of the change takes priority.

When a change bridges multiple concepts, the underlying intent takes priority[^21]. For example, if a formatting change is made to a unittest but does not affect the function of the unittest, the type should be _style_, not _test_. This helps to make the overall change clear right in the subject line.

2. Make atomic commits.

Wallace Freitas wrote an article about [best practices for commits](https://dev.to/wallacefreitas/best-practices-to-make-a-good-commit-writing-clean-effective-commit-messages-5eg9), and we will be utilizing his recommendation to make atomic commits. _Atomic_ refers to making only one change per commit. This makes changes easier to review, track, and if it comes down to it, revert. 

3. Use the scope option whenever possible.

The subject line scope is optional, but if it is included, it should consist of a noun describing a section of the codebase or documentation surrounded by parentheses, for example: `Docs(README): fix typo in line 5`[^22]. If it provides enough context in conjunction with the rest of the subject, using a scope may remove the need to complete the body of the commit message.

#### Subject Line

The subject should be limited to 50 characters and should include the type of change being made. This is to prevent important information from being automatically line-wrapped and hidden from initial view. 

The commit type can be one of the following, and should utilize the short form in parentheses to save space[^3]:
- Feature (_feat_): Introduces a new feature to the codebase[^4].
- Fix (_fix_): Patches a bug in the codebase[^5].
- Performance (_perf_): Adjustments made to optimize resource usage without altering the program's primary function[^6].
- Style (_style_): Aesthetic or formatting alterations that do not influence execution[^7].
- Refactor (_ref_): Modifications that reorganize existing code to boost overall maintainability or human-readability, strictly excluding any changes that qualify merely as stylistic or performance-based[^8].
- Documentation (_docs_): Updates applied to textual resources, including comments within the code, README files, and typo corrections[^9].
- Test (_test_): The creation or revision or automated testing scripts[^10].
- Continuous Integration (_ci_): Adjustments to continuous integration and deployment pipelines or scripts (e.g., GitHub Actions setups)[^11].
- Build (_build_): Alterations to external dependencies, compilation configurations, and project build tools[^12].
- Chore (_chore_): Any remaining development tasks that do not align with the specific labels above[^14].

In the subject description, use the imperative mood. According to [Beams](https://chris.beams.io/git-commit):
> [i]mperative mood just means "spoken or written as if giving a command or instruction". A few examples:
> - Clean your room
> - Close the door
> - Take out the trash
> 
> [... ]The imperative can sound a little rude; that’s why we don’t often use it. But it’s perfect for Git commit subject lines. One reason for this is that [...] when you write your commit messages in the imperative, you’re following Git’s own built-in conventions.

Beams recommends a quick way to gauge whether your subject line description is written in the imperative mood:
> **A properly formed Git commit subject line should always be able to complete the following sentence**:
> - If applied, this commit will  _your subject line here_

If the completed sentence makes grammatical sense, you're on the right track!

#### Co-Authors

If there are any co-authors on the commit, these should be added first after an inserted empty line[^15]. Before you can add a co-author, you must know the appropriate email to use if you want them to get contribution credit on GitHub[^16].

The structure to add co-author(s) to your commit message is as follows[^17]:

```shell
$ git commit -m "Refactor usability tests.
>
> Co-authored-by: NAME <NAME@EXAMPLE.COM>
> Co-authored-by: ANOTHER-NAME <ANOTHER-NAME@EXAMPLE.COM>"
```

Add another empty line[^18] before continuing to the content of the body. If there are no co-authors, only one empty line[^19] is required between the subject line and the body. The empty lines will aid in readability.

#### Body

While not entirely necessary for smaller, self-explanatory changes, further explanation than the brief subject line may be required or beneficial in many cases. The information that can be found in the body can provide important context for implementations that may seem unnecessary or counterintuitive but were made for very specific reasons, such as streamlining a particular workflow, or to take advantage of specific architecture. If it is not immediately obvious from the subject line description, please include the following information:
- **Where**: If necessary, give a detailed description of where the change is being made. Include file names and / or line numbers.
- **Current Status**: If applicable, provide a description of how things worked before the change, and what problems arose from this. Include instructions on how to reproduce the error.
- **Change**: If applicable, provide a description of how things work now, given the change.
- **Why**: Describe why you decided to solve the problem, and your reasoning behind the methods you used.
- **Testing**: If applicable, include a detailed description of how the changes were tested. Include manual testing.
- **Additional Notes**: If necessary, include more context with additional notes, such as references.[^20]

Structurally, ensure the body text wraps at 72 characters. This follows the docstrings and comments line length restrictions followed by [Python's PEP 8](https://peps.python.org/pep-0008/#maximum-line-length), and allows for users on a variety of screen resolutions to read your message without needing to scroll.

#### Footer

If applicable, link [directly](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/linking-a-pull-request-to-an-issue) to any existing related issues or PRs in the footer. Use the keywords provided by GitHub to refer to these issues or PRs. This allows reviewers to easily find the related issues and PRs to close or amend them, as necessary.

### In Conclusion

These contributing guidelines have been created to facilitate the long-term maintenance and development of the AMS by a team that spans different locations, timezones, and timeframes. Having consistent methods of contributing and documentation will also aid the future publishing process.

###### References:

[^1]: Inspiration for the ***Issues*** section drawn from Dr. James Hughes' [cs101 repository](https://github.com/jameshughes89/cs101).
[^2]: [Diff](https://git-scm.com/docs/git-diff): "Show changes between the working tree and the index or a tree, changes between the index and a tree, changes between two trees, changes resulting from a merge, changes between two blob objects, or changes between two files on disk."
[^3]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^4]: Modified from [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).
[^5]: Modified from [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).
[^6]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^7]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^8]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^9]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^10]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^11]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^12]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^14]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^15]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^16]: Creating a commit with multiple authors \(from [_GitHub_](https://docs.github.com/en/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/creating-a-commit-with-multiple-authors)\).
[^17]: Creating a commit with multiple authors \(from [_GitHub_](https://docs.github.com/en/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/creating-a-commit-with-multiple-authors)\).
[^18]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^19]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^20]: Chris Beams, [_How to Write a Git Commit Message_](https://chris.beams.io/git-commit).
[^21]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^22]: Modified from [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).
