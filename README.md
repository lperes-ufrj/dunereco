# dunereco with Pand+atm+bdm Codes 

## To install

Source DUNE:
```md
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
```
Then you need to setup larsoft, I used the version v09_53_02
```md
setup larsoft v09_53_02 -q e20:prof
```
Create a env in a work folder with
```md
mrb newDev
```
The source the local products 
```md
source localProducts*/setup
```

Then you can fork in the source direcroty
```md
cd srcs
```
```md
git clone (this repo)
```
Update cmake lists files
```md
mrb uc
```

Set environment
```md
mrbsetenv 
```
If everything is ok, you can install the packages 
```md
cd $MRB_BUILDDIR
```
```md
mrb i -j4
```

And voilá! You should probably be able to run the codes. 

