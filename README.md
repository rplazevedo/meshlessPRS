# modifiedPRS

This set of files is used to calculate the evolution of a domain wall network on a 2D grid using the PRS equaations, both in their standard form and in a version where only grid points that correspond to walls or radiation are evolved.

## File description

input.ini - Where the network and simulation parameters are stores. You can also name the set of data files that are output

run_files.py - The main UI. If you're just using the code you should only care about this file

cleanup.py - Deletes all the data files generated by the main code. Can delete everything or only particular files

sprs_init.py - Python file used for running the initial standard PRS simulation

sprs_cont.py - Python file used for running the rest of the standard PRS simulation

meshless_prs.py - Python file used for running the continuation of the standard simulation using a meshless algorithm

merge_v.py - Python file used for merging the velocity data files

merge_phi.py - Python file used for merging the field data files

reduce_data.py - Reduces the original data organizes it for plotting

plot_data.py - Plots the data

## How to use

In general, to obtain the data all you have to do is run "run_files.py" in a Python 3 interpreter. The cleanup file is run first, then the simulation is performed. The UI asks which parts you want to run, with the default option highlighted.

The data can be plotted using "plot_data.py".
