# AtmosphericAna
leoperes@pos.if.ufrj.br
===============================================

Module of Atmospheric Neutrino Analysis in DUNE 

It started to be developed at v39_00_00 of dunetpc repository, but it should work fine with newest versions.

The goal of this module is to provide a series of reconstructed variables that can be used in a Boosted decision tree analysis in order to obtain a NC/CC classification of the events and study the directionality of the original atmospheric neutrino.

================================================

Instructions to run:

You should install a local version of dunetpc.

Then, add this repository in the source, a path like

/dune/app/users/$USER/$WORK_DIRECTORY/srcs/dunetpc/dune/AtmosphericAna

You should also include the new directory in the CMakeLists file.

Install this code just once with

>>> cd $MRB_BUILDDIR 

>>> mrb i

Now, you are able to run and make changes.
After a change to compile again you must

>>> cd $MRB_BUILDDIR/dunetpc/dune/AtmosphericAna

>>> make install

And just AtmosphericAna codes are reinstalled.

================================================

The samples can be found in 

https://wiki.dunescience.org/wiki/High-E_and_Non-accelerator_Physics



