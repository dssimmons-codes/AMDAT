# Contributing

First off, thank you for considering contributing to AMDAT! Please take a moment to review this document **before submitting a pull request**.

Following these guidelines helps to communicate that you respect the time of the developers managing and developing this open source project. In return, they should reciprocate that respect in addressing your issue, assessing changes, and helping you finalize your pull requests.

AMDAT is an open source project and we love to receive contributions from our community — you! There are many ways to contribute, from writing tutorials or blog posts, improving the documentation, submitting bug reports and feature requests or writing code which can be incorporated into AMDAT itself.


## Pull requests

**Please ask first before starting work on any significant new features.**

It's never a fun experience to have your pull request declined after investing a lot of time and effort into a new feature.
 To avoid this from happening, we request that contributors create [a feature request](https://github.com/dssimmons-codes/AMDAT/discussions/new?categories=ideas) to first discuss any new ideas.
Your ideas and suggestions are welcome!

<!-- Please ensure that the tests are passing when submitting a pull request. -->
<!-- If you're adding new features to AMDAT, please include tests. -->

## Where do I go from here?

For any questions, support, or ideas, etc. [please create a GitHub discussion](https://github.com/dssimmons-codes/AMDAT/discussions/new).
If you've noticed a bug, [please submit an issue][new issue].

### Fork and create a branch

If this is something you think you can fix, then [fork AMDAT] and create
a branch with a descriptive name.

### Implement your fix or feature

Develop and test your fix on your fork.
At this point, you're ready to make your changes.
Feel free to ask for help.

### Create a Pull Request

At this point, if your changes look good and tests are passing, you are ready to create a pull request.

<!-- Github Actions will run our test suite against all supported AMDAT versions. It's possible that your changes pass tests in one AMDAT version but fail in another. In that case, you'll have to setup your development. -->

## Merging a PR (maintainers only)

A PR can only be merged into main by a maintainer if: CI is passing, approved by another maintainer and is up to date with the default branch.
Any maintainer is allowed to merge a PR if all of these conditions ae met.

## Shipping a release (maintainers only)

Maintainers need to do the following to push out a release:

* Create a feature branch from main and make sure it's up to date.
* Run standard test suite.
* Review and merge the PR.
* [Create a GitHub Release](https://github.com/dssimmons-codes/AMDAT/releases/new) by selecting the tag and generating the release notes.

[new issue]: https://github.com/dssimmons-codes/AMDAT/issues/new
