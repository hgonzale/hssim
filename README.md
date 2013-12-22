# Hybrid System Simulator

Each folder contains the implementation of an example in Section V in
[1], as detailed below.

# Forced Linear Oscillator with Stop

The folder 'mass_wall' contains a Matlab implementation of the
algorithm from Fig. 9 in [1], and the forced oscillator in Section V-A
in [1].

To generate the Figures 10d-f in [1] just type:
>> plot_comp

To generate Figure 10b, first type:
>> fmw_iterate
then type:
>> plot_grnd( grnd_trjs, user, 1 )

To generate Figure 10c, first type:
>> fmw_zeno
then type:
>> plot_grnd( grnd_trjs, user, 0 )

To implement your own version of the forced linear oscillator with
stop, you will need to begin by setting your own simulation parameters
as is done on lines 6-11 of 'fmw_iterate.m'. After setting up these
parameters you can load the forced linear oscillator with stop model
as is done on line 14 of 'fmw_iterate.m'. 

Different parameters for the force linear oscillator with stop, as in
Table 1 in [1], can be setup using the code on lines 17-30 of
'fmw_iterate.m'. Our simulation algorithm is then called using the
'fwd_RK2.m' file as is done on line 49 of 'fmw_iterate.m'. If you want
to call the PS method, use the 'PS_method.m' file as is done on line
56 of 'fmw_iterate.m'.

This code was tested on Matlab R2011b and R2013a.


# Verification of a Navigation Hybrid System

The folder 'navigation' contains a Matlab implementation of the
algorithm from Fig. 9 in [1], and the navigation benchmark in Section
V-B in [1].

Figures 11a-c were generated using randomly sampled trajectories but
we have code to generate figures similar to them.

To generate Figures 11a type:
>> nav_example
and then run:
>> plot_navbox( user, 4 )

To generate Figures 11b type:
>> nav_infexample
and then run:
>> plot_navbox( user, 4 )

To generate Figures 11c type:
>> nav_zenoexample
and then run:
>> plot_navbox( user, 4 )

To implement your own version of the verification task, you will need
to begin by setting your own simulation parameters as is done on lines
4-9 of 'nav_example.m'. After setting up these parameters you can load
the navigation model as is done on line 12 of 'nav_example.m'.

Different parameters for the navigation model can be setup using the
code on lines 18-19 of 'nav_example.m'. Our simulation algorithm is
then called using the 'fwd_RK2.m' file as is done on line 47 of
'nav_example.m'.

The code was tested on Matlab R2011b and R2013a.


# Polyped Locomotion Hybrid System Model

This folder contains a Python implementation of the algorithm from
Fig. 9 in [1], and the multi-legged locomotion example from Section
V-C in [1].

The following commands will generate Figs. 12 & 13 from [1]:
$ python sch.py pronk.cfg
$ python seq.py pronk.cfg
    
The figures are saved in both .pdf and .eps format in the 'fig' and
'seq' subdirectories, respectively.

The numerical simulation algorithm from Fig. 9 in [1] is implemented
in the Euler function in the 'relax.py' module. Users who wish to
simulate their own hybrid system can subclass the "hybrid dynamical
system" HDS class from 'relax.py'.

A particularly simple example to start from is the bouncing ball 'BB'
class in 'relax.py'. To simulate this model and generate a plot
containing a phase portrait and energy exchange versus time, simply
run 'relax.py' from iPython:
$ run relax.py

To modify the initial conditions or parameters for the polyped
locomotion simulation simply change the 'x0' and 'p0' variables in
'pronk.cfg'.

This code was tested on Python 2.7, NumPy 1.7, SciPy 0.12 (GCC 4.7),
Matplotlib 1.2, and iPython 0.13.

# References

[1] S.A. Burden, H. Gonzalez, R. Vasudevan, R. Bajcsy, and
S.S. Sastry, "Metrization and Simulation of Controlled Hybrid
Systems".  arXiv:1302.4402
