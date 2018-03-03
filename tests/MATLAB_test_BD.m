
% Try to load the library.
if not(libisloaded('libSpinW'))
    try
    [notfound, warnings] = loadlibrary(...
        fullfile(pwd,'..','libSpinW.so'),...
        fullfile(pwd,'..','include','Hspinw.h'));
    assert(isempty(notfound), 'Could not load test library')
    catch ME
        disp(ME.message)
        exit;
    end
end

% Create a base model
s = sw_model('triAF',1);

% Create the C strucures.
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

% Make the C SpinW object
c_sw = calllib('libSpinW','create_sw',latt, unit_cell, twin, mag_str, unit);

% This is a test to see if HKL generation works.
% hkl = calllib('libSpinW','sw_qscan',libpointer('doublePtr',[2, 3, 0 0 0, 0 0 1, 10]));

% Create the spinwave options C struct. 
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
opt.tid = -2;
opt.fid = -1;

% Call the spinwave function.
% HKL is [0 0 0] -> [0 0 1] and 10 points.
hkl = libpointer('doublePtr',[2, 3, 0 0 0, 0 0 1, 10]);
newSpec = calllib('libSpinW', 'sw_spinwave', c_sw, hkl, opt);

% Unload the library
clear variables 
unloadlibrary('libSpinW')
