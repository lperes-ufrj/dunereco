# dunereco with Pand+atm+bdm Codes 

## To install
[![iuricode](https://github-readme-stats.vercel.app/api/pin/?username=iuricode&repo=readme-template)](https://github.com/iuricode/)

$ source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
Then you need to setup larsoft, I used the version v09_53_02
$ setup larsoft v09_53_02 -q e20:prof
Create a env in a work folder with
$ mrb newDev
The source the local products 
$ source localProducts*/setup
Then we can fork in the source direcroty
$ cd srcs
$ git clone (this repo)
You also gonna need dunesim repo from DUNE (check the tags for the same larsift version)
$ git clone (dunesim)
Update cmake lists files
$ mrb uc
Set environment
$ mrbsetenv 
If everything is ok, you can install the packages 
$ cd $MRB_BUILDDIR
$ mrb i -j4

And voilá! You should probably be able to run the codes. 

