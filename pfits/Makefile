cc = gcc
#objects = pfits_fv pfits_describe pfits_plot pfits_process pfits_process_old pfits_statistics pfits_draw
objects = pfits_fv pfits_describe pfits_plot pfits_statistics pfits_draw


all:$(objects)
pfits_fv:pfits_fv.c
	$(cc) -o pfits_fv pfits_fv.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

pfits_describe:pfits_describe.c
	$(cc) -o pfits_describe pfits_describe.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

pfits_plot:pfits_plot.c
	$(cc) -o pfits_plot pfits_plot.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

pfits_process:pfits_process.c pfits_setup.c
	$(cc) -o pfits_process pfits_process.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

pfits_process_old:pfits_process_old.c
	$(cc) -o pfits_process_old pfits_process_old.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

pfits_statistics:pfits_statistics.c
	$(cc) -o pfits_statistics pfits_statistics.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11  -lgfortran -lcpgplot  -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
	
pfits_draw:pfits_draw.c
	$(cc) -o pfits_draw pfits_draw.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot  -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio


#rm files compiled

.PHONY : clean
clean :
	rm $(objects)

