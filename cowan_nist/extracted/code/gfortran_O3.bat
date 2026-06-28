set LIBPATH="C:\MinGW\mingw64\lib"
gfortran.exe -fshort-enums -malign-double -ftracer -fno-backslash -O3 -o %1.exe %1.F 2>gfortran_out.txt
