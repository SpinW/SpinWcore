# SpinW_core
This is the core functionality of SpinW written in C++. At the moment it is
in development as a proof of concept and _DOES NOT WORK!_

## Requirements
In order to build SpinW_core the following is required:
- CMake >= 3.6 
- gcc > 5.3
- Armadillo 
- Intel MKL

And OpenMP development libraries if parallelization is required (recomended). 

## Testing - MATLAB
The library can be used in MATLAB by calling the library:

```matlab
[notfound, warnings] = loadlibrary(fullfile(pwd,'libSpinW.so'), 'Hspinw.h');
```

Then transferring the SpinW object (s):

```matlab
latt  = libstruct('lattice',setfield(s.lattice,'nSymOp',size(s.lattice.sym,3)));
unit_cell = libstruct('unit_cell',setfield(s.unitcell,'nAtom',size(s.unit_cell.r,2)));
twin = libstruct('twin',setfield(s.twin,'nTwin',size(s.twin,2)));
unit = libstruct('unit',s.unit);
mag_str = s.mag_str;
mag_str.F_real = real(mag_str.F);
mag_str.F_imag = imag(mag_str.F);
mag_str.nMagExt = size(mag_str.F,2);

mag_str.nK = size(mag_str.F,3);
mag_str = rmfield(mag_str,'F');
mag_str = libstruct('mag_str',mag_str);

c_sw = calllib('libSpinW','create_sw',latt, unit_cell, twin, mag_str, unit);
```

Options for the spinwave generation have to be given in a slightly different syntax:

```matlab
opt = libstruct('spinwave_opt');
opt.notwin = false;
opt.sortMode = true;
opt.optmem = 0;
opt.toll = 1e-4;
opt.hermit = true;
opt.omega_toll = 1e-5;
opt.formfact = false;
opt.gtensor = false;
opt.cmplxBase = false;
opt.tid = -1;
opt.fid = -1;
```

The HKL range is given by creating a C-pointer where the first 2 elements are
the number of HKL end points and number of elements in HKL. i.e from [0 0 0] to [0 0 1] is 2 and the number of elements in [0 0 1] is 3.

```matlab
libpointer('doublePtr',[2, 3, 0 0 0, 0 0 1, 10])
```

The Spectrum can be calculated by calling ```sw_spinwave``` and a double C-Pointer to the memory address is returned.

```matlab
sw_spec = calllib('libSpinW','sw_spinwave',c_sw,hkl,opt)
``` 

## Testing - Python
A C interface for the spinwave function is available but more development is needed.
