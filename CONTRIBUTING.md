# Contibuting to Palabos

Please read carefully before contributing to Palabos.

## Code contribution to the Palabos 

In order to contribute to the Palabos source code it is required that you sign off on the [Developer Certificate of Origin (DCO)](https://gitlab.com/unigespc/palabos/blob/dco/dco/README.md).

In line with our value that everyone can contribute to Palabos, we strive to make our process as developer-friendly and frictionless as possible. As such, we consider the act of contributing to the code by submitting a *Merge Request* as the "Sign off" or agreement to the certifications and terms of the DCO and AGPL3 license. No further action is required.

## Bug reports

When you think you have identified a bug, please open an issue and do not forget to give as many details as possible for us to be able to
reproduce the bug. Ideally the report should contain a minimal snippets of code where the bug happens, along with the output you are getting.
Otherwise it is close to impossible for us to try to figure out what is happening.

## Novel contributions

Contributions to Palabos are accepted by community members who prove sufficient programming skills, familiarity with the Palabos source code, and level of maturity with respect to the proposed numerical or physical models.

Before contributing to Palabos, make sure to be familiar with the programming 
guidelines listed in the User's Guide (see [here](https://palabos.unige.ch/palabos-documentation/)). Also, before jumping into a large development project, it is good practice to interact with the core developers and the community on the forum, to identify if the development project answers a need in the community and is aligned with the general development goals of Palabos.

In order to contribute to Palabos you must create a fork of the project (click on
the fork button, see [here](https://docs.gitlab.com/ee/workflow/forking_workflow.html#creating-a-fork) for more info).
This will create an exact copy of the Palabos repo in your namespace.

On your fork of the project you are free to make any modification you want to
the source code, but for your code to be reusable by others you should not modify existing classes
but create your own new classes.

When implementing a model from an existing paper please include the reference to the paper in the
comments. If the paper is not published yet please include a link to a preprint (ideally on arXiv).

You should illustrate how to use your contribution. To do so
add an example test case in `examples/codesByTopic`. See how the examples there are built
there and try to make something similar. If your contribution is used in a
physical context not present in `examples/showCases` you could instead create
a new show case.

When you think your contributions are ready to be added to Palabos, you should create
a `merge request` (see [here](https://docs.gitlab.com/ee/gitlab-basics/add-merge-request.html) for how to 
create one). In your `merge request` briefly summarize what your additions are and highlight the particularly important ones.

## Feature requests

Feature requests can be created as issues. Keep in mind that we are an open source community
and that it may not be possible for us to answer to any of the feature requests
since we have limited amount of time. But we may be able to provide guidance
for you to create the new features you would like to see in Palabos.
