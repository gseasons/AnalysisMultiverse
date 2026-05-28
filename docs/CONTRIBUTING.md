# About the Contribution Guidelines

The Analysis Multiverse software (AMS) and its documentation should assume no computer science or programming background.

The AMS is intended to process functional magnetic resonance imaging (fMRI) data. The anticipated userbase for the AMS will be researchers from the fields of psychology, neuropsychology, neuroscience, and biology. These researchers may include student researchers. With this in mind, code should be as self-explanatory as possible. The documentation should assume no previous experience
with running or installing software from GitHub, it should assume no previous experience with high performance computing (HPC), and should not assume any fluency with Python or any other programming languages. In essence, the documentation should be as accessible as possible, to facilitate the use of AMS by all potential users.

These guidelines, as well as the documentation itself, are living documents. This means that they are meant to be updated and amended a necessary, dictated by the needs of the AMS contributors as well as the needs of the AMS users.

# How to Contribute to the Repository

## Issues[^1]

This section lists the suggested use for issues on GitHub, but please keep in mind that the documentation is still a work in progress, and the codebase is still being tested and developed. The issue you highlight may already be on our "To Do" list! 

Use issues to...

 - Document and report bugs in the code
 - Identify functions, processes, or code files that lack adequate documentation within the code, or on our documentation website
 - Suggest clarifications for a specific section of the documentation, or in general
 - Suggest content that should be included in the documentation

## Pull Requests
Fork the repository to submit pull requests (PRs) to contribute to the
documentation or the codebase. Use the template provided when you go to submit the PR, replacing the text under each provided header. If specific header is not applicable, replace the template text with "N/A".

Use an appropriate title based on the commit message subject formatting listed in the ***Commits*** section below.

This information will be used to create a commit message if your pull request is accepted, so please be thorough in your documentation. Improper documentation of your changes will result in a request for changes until the documentation is adequate.

## Commits
Why do good commit message matter?

Commit messages make up a permanent log of changes to a Git repository. Depending on the repository, the timeframe, and the number of contributors, the commit history can vary widely. As said by [Chris Beams](https://chris.beams.io/git-commit):
> \[...\] a well-crafted Git commit message is the best way to communicate context about a change to fellow developers (and indeed to their future selves). A diff[^2] will tell you ***what*** changed, but only the commit message can properly tell you ***why***.

Beams expands further on what will be our main point here: 

> A project’s long-term success rests (among other things) on its maintainability, and a maintainer has few tools more powerful than [their] project’s log. It’s worth taking the time to learn how to care for one properly. What may be a hassle at first soon becomes habit, and eventually a source of pride and productivity for all involved.

Commit messages are a form of data and as such, part of the data management plan for AMS includes the format of commit messages as a part of our *commit*ment (pun intended) to the [FAIR principles](https://www.go-fair.org/fair-principles/) of Accessibility and Interoperability. In this spirit, we will blend the advice from Beams' article "[How to Write a Git Commit Message](https://chris.beams.io/git-commit)" with the formatting structure defined by [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).

The required format will be as follows, with further details below:

| Element | Format |
|-----------|-------|
| SUBJECT | \<Type\>\[optional scope\]: \<description\> |
| EMPTY LINE | empty |
| CO-AUTHORS | \[optional\] |
| EMPTY LINE | empty |
| BODY | \[optional\] |
| EMPTY LINE | empty |
| FOOTER | \[optional\] |

### Subject Line

The subject should be limited to 50 characters, and must include a capitalized type. The commit type can be one of the following, and should utilize the short form in parentheses[^3]:
- Feature (_Feat_): Introduces a new feature to the codebase[^4].
- Fix (_Fix_): Patches a bug in the codebase[^5].
- Performance (_Perf_): Adjustments made to optimize resource usage without altering the program's primary function[^6].
- Style (_Style_): Aesthetic or formatting alterations that do not influence execution[^7].
- Refactor (_Ref_): Modifications that reorganize existing code to boost overall maintainability or human-readability, strictly excluding any changes that qualify merely as stylistic or performance-based[^8].
- Documentation (_Docs_): Updates applied to textual resources, including comments within the code, README files, and typo corrections[^9].
- Test (_Test_): The creation or revision or automated testing scripts[^10].
- Continuous Integration (_CI_): Adjustments to continuous integration and deployment pipelines or scripts (e.g., GitHub Actions setups)[^11].
- Build (_Build_): Alterations to external dependencies, compilation configurations, and project build tool[^12].
- Revert (_Revert_): A commit that is being reverted[^13].
- Chore (_Chore_): Any remaining development tasks that do not align with the specific labels above[^14].

Do not end the subject line with a period, as it makes unnecessary use of 2% of the subject line space!

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

### Co-Authors

If there are any co-authors on the commit, these should be added first after an inserted empty line[^15]. But before you can add a co-author, you must know the appropriate email to use if you want them to get contribution credit on GitHub[^16].

The structure to add co-author(s) to your commit message is as follows[^17]:

```shell
$ git commit -m "Refactor usability tests.
>
> Co-authored-by: NAME <NAME@EXAMPLE.COM>
> Co-authored-by: ANOTHER-NAME <ANOTHER-NAME@EXAMPLE.COM>"
```

Add another empty line[^18] before continuing to the content of the body. If there are no co-authors, only one empty line[^19] is required between the subject line and the body.

### Body

While not entirely necessary for smaller, self-explanatory changes, further explanation than the brief subject line may be required or beneficial in many cases. If it is not immediately obvious from the subject line description, include the following information:
- **Where**: If necessary, give a detailed description of where the change is being made. Include file names and line numbers.
- **Current Status**: If applicable, provide a description of how things worked before the change, and what problems arose from this. Include instructions on how to reproduce the error.
- **Change**: If applicable, provide a description of how things work now, given the change.
- **Why**: Describe why you decided to solve the problem, and your reasoning behind the methods you used.
- **Testing**: If applicable, include a detailed description of how the changes were tested. Include manual testing.
- **Additional Notes**: If necessary, include more context with additional notes, such as references.[^20]

Structurally, ensure the body text wraps at 72 characters. This follows the docstrings and comments line length restrictions followed by [Python's PEP 8](https://peps.python.org/pep-0008/#maximum-line-length).

### Footer

If applicable, link to any existing related issues or PRs. Link [directly](https://docs.github.com/en/issues/tracking-your-work-with-issues/using-issues/linking-a-pull-request-to-an-issue) if you can. Use the keywords provided by GitHub to refer to these issues or PRs.

### In General

> 1. The underlying intent of the change takes priority.

When a change bridges multiple concepts, the underlying intent takes priority[^21]. For example, if a formatting change is made to a unittest but does not affect the function of the unittest, the type should be _Style_, not _Test_.

> 2. Make atomic commits.

Wallace Freitas wrote an article about [best practices for commits](https://dev.to/wallacefreitas/best-practices-to-make-a-good-commit-writing-clean-effective-commit-messages-5eg9), and we will be utilizing his recommendation to make atomic commits. _Atomic_ refers to making only one change per commit. This makes changes easier to review, track and, if it comes down to it, revert. 

> 3. Use the scope option whenever possible.

The subject line scope is optional, but if it is included, it should consist of a noun describing a section of the codebase or documentation surrounded by parentheses, for example: `Docs(README): fix typo in line 5`[^22]. Using the scope may remove the need to complete the body of the commit message, if it provides enough information when combined with the description.

#### References[^23]:

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
[^13]: Modified from [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).
[^14]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^15]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^16]: Creating a commit with multiple authors \(from [_GitHub_](https://docs.github.com/en/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/creating-a-commit-with-multiple-authors)\).
[^17]: Creating a commit with multiple authors \(from [_GitHub_](https://docs.github.com/en/pull-requests/committing-changes-to-your-project/creating-and-editing-commits/creating-a-commit-with-multiple-authors)\).
[^18]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^19]: If using a text editor on command line to create your commit message, use two consecutive newlines.
[^20]: Chris Beams, [_How to Write a Git Commit Message_](https://chris.beams.io/git-commit).
[^21]: Modified from the [Wikipedia entry on the Conventional Commits specification](https://en.wikipedia.org/wiki/Conventional_Commits_Specification).
[^22]: Modified from [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).
[^23]: The Markdown of this document was edited with [StackEdit](https://stackedit.io/).