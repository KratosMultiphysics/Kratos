# To contribute to the project, you need to do it through pull/merge request

First you need to fork the repository into your own account. You can do that simply by clicking the fork button on the gitlab interface.

https://gitlab.inria.fr/solverstack/morse_cmake/forks/new

Then, clone the repository on your laptop:
```
#!shell
git clone git@gitlab.inria.fr:username/forkname.git
```

Once this is done, you can setup the morse_cmake repository as the upstream of your clone to simplify the update of your fork repository.
```
#!shell
git remote add upstream git@gitlab.inria.fr:solverstack/morse_cmake.git
```

Now, you have your repository configured, and you want to create a new pull request. The first step is to create a branch from the HEAD of the your fork repository.
```
#!shell
git branch your_branch_name
git checkout your_branch_name
```

Apply your modifications in your branch. Then, you need to push this branch on your online repository
```
#!shell
git push -f origin your_branch_name
```

or without -f, if the branch already exists online, and you just want to update it.

Once your branch is online, on the gitlab interface, go to the branches webpage, select the branch you want to push as a merge request, and push the button !!!

***Be careful to check the 'close after merge' check box, and to push to the solverstack/morse_cmake repository*** By default the checkbox is not checked, and the default repository is your fork.


# To review locally a private pull request submitted by someone else

Get the patch from the pull request (Need to update that !!!! Coming from bitbucket)
```
#!shell
curl https://bitbucket.org/api/2.0/repositories/icldistcomp/parsec/pullrequests/#PR/patch > pr#PR.patch
```

Than apply the patch on your local copy

```
#!shell

git apply pr#PR.patch
```
