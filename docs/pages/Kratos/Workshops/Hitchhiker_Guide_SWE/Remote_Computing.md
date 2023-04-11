---
title: Remote Computing
keywords: 
tags: [Remote_Computing.md]
sidebar: kratos_workshops
summary: 
---
# Remote Computing
Before continuing further with the remote computing, make sure the following has been done:
- Have the necessary tools for Remote Computing installed ([eduVPN](Installation_Guides.html#41-vpn){:target="_blank"} and [Git Bash](Installation_Guides.html#42-secure-shell-ssh-client){:target="_blank"}).
- Check if the GiD Model of the building as well as the wind characteristics are correct.
- Check if the simulation [parameters](Preprocessing.html#3-simulation-parameters){:target="_blank"} and [processes](Preprocessing.html#2-processes){:target="_blank"} of the ProjectParametersCustom.json and MainKratosCustom.py are adjusted.
- Test if the run begins (3-5 time steps, (with all possible necessary output turned on and additional processes)).

Now you are prepared for the full simulation on the chairs cluster.

___
## 1. Useful Linux Commands
As a start, it is useful to get to know a couple of linux commands. The angle brackets "< >" are used to indicate a placeholder value, which you should replace with your own input.

| Command | Function | 
| -------- | -------- |
| pwd | Shows the current working directory. | 
| cd `<path>` | Goes to the selected path. | 
| cd .. | Goes back to the parent folder. | 
| ls | Lists all the files and folders that are at the current location | 
| mv `<file-name>` `<new-path>` | Moves a file to a certain location. It is also used to change the name of a file. | 
| mv -r `<folder-name>` `<new-path>` | The "-r" option means "recursively". Moves the folder and all inside to a certain location. | 
| cp `<file-name>` `<new-name>` | Copies a file. | 
| cp -r `<file-name>` `<new-name>` | The "-r" option means "recursively". Copies a folder and all inside it. | 
| vim | Visualizes a text file. Pressing "I" you enter in "insert mode" to edit it. To exit press "Esc" and the ":w" to save it and ":q" to quit. If the name of the file does not exist it creates a new one. | 
| mkdir | Creates a directory | 
| htop | Checks the usage of processors and RAM | 
| tail - n `<num-lines>` `<file-name>` | Prints in the screen the last <num-lines> rows of a text file. Very useful to check if the simulation is running. | 
| sh | Execute an .sh file | 

___
## 2. Login to cluster
 
You will need to initially be connected to the university network to find the computers. For that, use EduVPN.

The Credentials for the CIP-Pool computer (normal computer) and the Cluster will be given during the course. You will need the following:

- User name
- IP-address
- Password

The login is done directly through the "ssh" command in the Git Bash terminal:

```shell
$ ssh -oHostKeyAlgorithms=+ssh-dss <user-name>@<IP-address>
```

**Note: Some operating systems might have problems connecting via ssh to the cluster, due to the old operating system at the cluster, that is why we use the parameter "-oHostKeyAlgorithms=+ssh-dss".*

___
## 3. Transfer files
  
The file structure in the cluster is divided in the following folders:
- Software:
  - Kratos &rarr; Folder where Kratos is compiled. **You should not touch it.**
  - setup_kratos.sh &rarr; File to add Kratos to the path. **You should not touch it.**
- Documents:
  - Tests &rarr; Folder with some examples and tests.
  - Templates &rarr; Here are saved some important files you will need to copy.
  - Group`<group-data>` &rarr; Multiple folders, each for one group. Here you can save all your documents and simulations.
  
The next step is to transfer your simulation file to your group folder in the cluster. To transfer files from your personal computer to the remote computers (or the other way around) we use the secure copy ("scp") command:
  
- To copy from your personal computer to the remote location, type in your terminal (mind the spaces):
  ```shell  
  $ scp -r -oHostKeyAlgorithms=+ssh-dss <folder-to-copy> <user-name>@<IP-address>:<path>
  ```

- To copy from the statik computer to the remote location, type in your terminal:
  ```shell
  $ scp -r -oHostKeyAlgorithms=+ssh-dss <user-name>@<IP-address>:<path-folder-to-copy> <paste-destination>
  ```
  
Remember to put the parameter `-r` so the folder and all the contents can be copied. If you only want a file, omit the parameter `-r`. 

___
## 4. Running in the cluster

First, here are some useful commands for the cluster:

| Command | Function | 
| -------- | -------- |
| qstat | Shows own submitted jobs. | 
| qdel `<job-id>` | Cancel a job | 
| qstat -f | Show the state of all the nodes in the computer, including where each job is running.| 

Let us imagine, that you want to run a script (e.g. MainKratosCustom.py) located in a certain folder of the cluster. These would be the steps to launch the job with the qsub command:

- Log in on the head of the Statik cluster (see instructions above).
- Navigate to the folder where you have the script you want to run:
```shell
$ cd <path-to-simulation-folder>
```
- Using the template you have in `../Documents/Templates/q_run.sh` create a **q_run.sh** file and save it in the same folder than your **MainKratosCustom.py**. This file is in charge of running the simulation. **Please specify your** ***\<job-name\>*** **(it must start with a "G" followed by the number of your group, and preferably something below 8 characters total)**. You can also change the error (after -e) and output (after -o) file names, which will serve to output any errors that might happen during the simulation and the standard prompt output respectively. By last, write the command that you want to execute. In your case: 
```shell
$ mpirun -x LD_LIBRARY_PATH --bind-to core --map-by socket python3 MainKratosCustom.py
```
The simulation will be using MPI parallel processing feature, which speeds up your simulation. Please make sure you change the `"parallel_type"` from `"OpenMp"` to `"MPI"` in "ProjectParametersCustom.json" as:
```json
{
    "problem_data"     : {
        "parallel_type" : "MPI"
    }
}
```

- Using the template you have in `../Documents/Templates/q_submit.sh` create a **q_submit.sh** file and save it in the same folder than your **MainKratosCustom.py**. This file helps you launch the job, calling the **q_run.sh** file. To do it, just execute it:
```shell
$ sh q_submit.sh
```

- As an alternative to the previous step, one can also run the command to launch the job directly in the terminal, so the **q_submit.sh** file is not needed:
```shell
$ qsub -pe impi_tight_fu <num-of-cores> -V q_run.sh
```

- Check that the simulation is running with the **qstat** or **qstat -f** commands. If it is running, go to one of the output files and use the following command to visualize the last lines and check that they are being written:
```shell
$ tail -n <number-of-lines-to-print> <output-file-name>
```   

When the simulation ends, the nodes should be liberated automatically. You can then [copy the results back to your computer](#3-transfer-files), the same way you copied the results from your computer to the cluster, just the other way around. The next step is the [postprocessing of the simulation results](Postprocessing.html){:target="_blank"}.


