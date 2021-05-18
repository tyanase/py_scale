# py_scale


## Motivation:
* call SCALE subroutine from python
* reference:
    * SCALE, by Team-SCALE@R-CCS : https://scale.riken.jp/
    * Fortran+Python, by Nakao-san@R-CCS: https://mnakao.net/data/2018/HPFPC.pdf


## How to
### Fortran-side:
* SCALE compiled with shared library options
    * add "-fPIC" option to sysdep/Makedef* file
* Interfacial subroutine program file
    * call SCALE subroutine
    * be called by python script
### Python-side: 
* use "ctypes" library
    * call interfacial subroutine


## Example
* set environmental parameters for SCALE compilation (scale 5.3.3, gfortran 4.8.5, netcdf 4.7.3, and openmpi 4.0.2)
```
$ export SCALE_SYS="Linux64-gnu-ompi"
$ export SCALE_NETCDF_LIBS="-L/usr/local/netcdf-c-4.7.3/lib -lnetcdff -lnetcdf"
$ export SCALE_NETCDF_INCLUDE="-I/usr/local/netcdf-c-4.7.3/include"
$ export SCALE_ENABLE_OPENMP="F"
```

* download SCALE
```
$ wget https://scale.riken.jp/archives/scale-5.3.3.tar.gz
$ tar zxvf scale-5.3.3.tar.gz
$ mv scale-5.3.3 py_scale-5.3.3
```

* modify sysdep/Makedef.Linux64-gnu-ompi
```
FFLAGS_CORE = -cpp -m64                                                                     \
              -std=f2003 -fall-intrinsics -pedantic-errors -fimplicit-none -fmodule-private \
              -fconvert=big-endian -frecord-marker=4 -ffree-form -ffree-line-length-none    \
              -fPIC
FFLAGS_DYN   = $(FFLAGS) -fPIC
CFLAGS_FAST  = -O3 -m64 -fPIC
```

* compile SCALE (refer to "SCALE User's Guide Chapter 5.3 How to use scale library")
```
$ cd py_scale-5.3.3/scalelib/src
$ make -j 2
```

* compile interfacial subroutine program (sample.f90)
```
$ mpif90 -shared -fPIC -o sample.so sample.f90 -J py_scale-5.3.3/include `nc-config --cflags` -L py_scale-5.3.3/lib -lscale `nc-config --libs`
```

* run python script (sample.py; python 3.7.4 w/ matplotlib, numpy)
    * diagnostic of surface momentum flux by calling fortran subroutine (sample.so)
```
$ python sample.py
```

* check whether SFLX_MOM*.png files are made
