types.o : params.o
globalvars.o : types.o params.o

#old.o : params.o flags.h
mod_stokes.o : flags.h params.o types.o globalvars.o mod_products.o mod_SDLP.o mod_compcoeffdens.o mod_EHpqr.o 
mod_SDLP.o : types.o globalvars.o
mod_runge.o : params.o globalvars.o mod_geom.o
mod_products.o : params.o
mod_EHpqr.o : flags.h params.o types.o old.o
mod_compcoeffdens.o : params.o types.o globalvars.o

driver.o : params.o globalvars.o mod_init.o mod_geom.o mod_target.o 
mod_init.o : globalvars.o mod_geom.o
mod_geom.o : params.o types.o globalvars.o mod_products.o
mod_target.o : params.o globalvars.o mod_geom.o
testcorr.o : params.o types.o globalvars.o mod_geom.o mod_EHpqr.o mod_stokes.o

