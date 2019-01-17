# To contribute to the project, you need to do it through pull/merge request

First you need to fork the repository into your own account. You can do that simply by clicking the fork button on the gitlab interface.

https://gitlab.inria.fr/solverstack/spm/forks/new

Then, clone the repository on your laptop:
```
#!shell
git clone git@gitlab.inria.fr:username/forkname.git
```

Once this is done, you can setup the SpM repository as the upstream of your clone to simplify the update of your fork repository.
```
#!shell
git remote add upstream git@bitbucket.org:solverstack/spm.git
```

Now, you have your repository configured, and you want to create a new pull request. The first step is to create a branch from the HEAD of the your fork repository.
```
#!shell
git branch your_branch_name
git checkout your_branch_name
```

Commit your modifications into your branch. Then, you need to push this new branch to your online repository
```
#!shell
git push origin your_branch_name
```

Once your branch is online, on the gitlab interface, go to the branches webpage, select the branch you want to push as a merge request, and push the button !!!

***Be careful to check the 'close after merge' check box, and to push to the solverstack/spm repository*** By default the checkbox is not checked, but the repository should be the correct one by default.


# To review locally a private pull request submitted by someone else

Follow the instructions describes by the "Check out branch" button associated to the merge request.
