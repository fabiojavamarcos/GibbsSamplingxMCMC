{\rtf1\ansi\ansicpg1252\cocoartf1671\cocoasubrtf600
{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl280\partightenfactor0

\f0\fs24 \cf2 \expnd0\expndtw0\kerning0
# Contributing Guide\
\
After reading through this guide, check out some good first issues to contribute to by checking the labels.\
\
## Understanding the basics\
\
We welcome contributions and encourage you to follow the GitHub workflow specified below. If you are not familiar with this type of workflow, take a look at GitHub's excellent overview on the [GitHub flow](https://guides.github.com/introduction/flow/index.html) and the explanation of [Feature Branch Workflow](https://atlassian.com/git/tutorials/comparing-workflows#feature-branch-workflow) for the idea behind this kind of development.\
\
1. Get the code on your local machine.     \
	1. Fork the project into your GitHub account.\
    	2. Clone your forked repository on your local machine.\
2. **Create a new branch** (such as `fix-for-issue-121`). Be sure to create a **separate branch** for each improvement you implement.\
3. Do your work on the **new branch - not the master branch.** \
4. Create a pull request. For an overview of pull requests, take a look at GitHub's [pull request help documentation](https://help.github.com/articles/about-pull-requests/).\
5. In case your pull request is not yet complete or not yet ready for review, consider creating a [draft pull request](https://github.blog/2019-02-14-introducing-draft-pull-requests/) instead.\
\
## Formal requirements for a pull request\
\
The main goal of the formal requirements is to provide credit to you and to be able to understand the patch.\
\
### Add your change to `CHANGELOG.md`\
\
You should edit the [CHANGELOG.md](CHANGELOG.md) located in the root directory of the project source.\
Add a line with your changes in the appropriate section.\
\
If you did internal refactorings or improvements not visible to the user (e.g., UI, .bib file), then you don't need to put an entry there.\
\
#### Format of keyboard shortcuts\
\
Example: `<kbd>Ctrl</kbd> + <kbd>Enter</kbd>`\
\
In case you add keys to the changelog, please follow these rules:\
\
- `<kbd>` tag for each key\
- First letter of key capitalized\
- Combined keys separated by `+`\
- Spaces before and after separator `+`\
\
### Author credits\
\
You will be given credit in the [`AUTHORS`](AUTHORS) file in the root of the repository \
\
Please, **do not add yourself at JavaDoc's `@authors`**.\
The contribution information is tracked via the version control system.\
\
Your contribution is considered being made under [MIT license](https://tldrlegal.com/license/mit-license).\
\
### Write a good commit message\
\
See [good commit message] or [commit guidelines section of Pro Git].\
The first line of your commit message is automatically taken as the title for the pull-request.\
All other lines make up the body of the pull request. Add the words `fixes #xxx` to your PR to auto-close the corresponding issue.\
\
### Test your code\
\
We know that writing test cases takes a lot of time.\
Nevertheless, we rely on our test cases to ensure that a bug fix or a feature implementation doesn't break anything.\
### When adding a library\
\
Please try to use a version available at JCenter In any case, describe the library at [external-libraries.txt](external-libraries.txt).\
\
Also add a txt file stating the license in `libraries/`.\
\
\
### When making an architectural decision\
\
In case you add a library or do major code rewrites, we ask you to document your decision.\
Recommended reading: <https://adr.github.io/>.\
\
We simply ask to create a new markdown file in `docs/adr` following the template presented at <https://adr.github.io/madr/>.\
\
In case you want to directly add a comment to a class, simply use following template (based on [sustainable architectural decisions](https://www.infoq.com/articles/sustainable-architectural-design-decisions)):\
\
```text\
In the context of <use case/user story u>,\
facing <concern c>\
we decided for <option o>\
and neglected <other options>,\
to achieve <system qualities/desired consequences>,\
accepting <downside / undesired consequences>,\
because <additional rationale>.\
```\
\
\
## Create a pull request\
\
Create a pull request on GitHub following GitHub's guide "[Creating a pull request from a fork](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork)".\
For text inspirations, consider [How to write the perfect pull request](https://github.com/blog/1943-how-to-write-the-perfect-pull-request).\
\
If you want to indicate that a pull request is not yet complete **before** creating the pull request, you may consider creating a [draft pull request](https://github.blog/2019-02-14-introducing-draft-pull-requests/).\
Alternatively, once the PR has been created, you can add the prefix `[WIP]` (which stands for "Work in Progress") to indicate that the pull request is not yet complete, but you want to discuss something or inform about the current state of affairs.\
\
[commit guidelines section of Pro Git]: http://git-scm.com/book/en/Distributed-Git-Contributing-to-a-Project#Commit-Guidelines\
[good commit message]: https://github.com/joelparkerhenderson/git_commit_message\
\
Feel free to contact me: fabiojavamarcos@gmail.com \
\
Thanks!\
}