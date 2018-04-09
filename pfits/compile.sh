gcc -o pfits_fv pfits_fv.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_describe pfits_describe.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_plot pfits_plot.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_process pfits_process.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_process_old pfits_process_old.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_statistics pfits_statistics.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11  -lgfortran -lcpgplot  -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio
gcc -o pfits_draw pfits_draw.c pfits_loader.c pfits_setup.c -I/public/home/mcc/psrsoft/pgplot -L/public/home/mcc/psrsoft/pgplot -lpgplot -lX11 -lgfortran -lcpgplot  -L/public/home/mcc/psrsoft/cfitsio/lib -lcfitsio

