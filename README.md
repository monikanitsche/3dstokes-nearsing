/src      contains all source F90 files
/src/common: all common modules and driver
/src/sphere: all modules specific to sphere
/src/ellipse: all modules specific to ellipse

/bld        here you build driver.out by typing make
/bld/sphere: driver.out flow past sphere at 45 degrees
/bld/ellipse: make different executables for different ellipses/uinf
              driver1.out flow past ellipse1 : 123 at 0 deg
              driver2.out flow past ellipse2 : ??
              driver3.out flow past ellipse2 : ??
              they use different init files


/runs            contains script for runs, all input, 
                 all output files, all plotting routines
/runs/sphere: flow past sphere at 45 degrees
/runs/sphere/datout4_crosssec: pb4 - velo in crossection with/without

/runs/ellip1: flow past 123 ellipse at 0 degrees
/runs/ellip1/datout4_crosssec: pb4 - velo in crossection with/without
