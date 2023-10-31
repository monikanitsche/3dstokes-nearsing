types.o : params.o
globalvars.o : types.o params.o

mod_stokes.o : flags.h params.o types.o globalvars.o mod_products.o mod_SDLP.o mod_compcoeffdens.o mod_EHpqr.o
mod_mod_SDLP.o : types.o
mod_runge.o : params.o globalvars.o mod_geom.o
mod_products.o : params.o
mod_EHpqr.o : flags.h params.o types.o
mod_compcoeffdens.o : params.o types.o globalvars.o

driver.o : params.o globalvars.o mod_init.o mod_geom.o mod_target.o
mod_init.o : globalvars.o mod_geom.o
mod_geom.o : params.o types.o globalvars.o mod_products.o
mod_target.o : params.o

