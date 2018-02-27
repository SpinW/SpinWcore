if not(libisloaded('libSpinW'))
    [notfound, warnings] = loadlibrary(fullfile(pwd,'cmake-build-debug-gcc_7','libSpinW.so'), 'Hspinw.h');
%     [notfound, warnings] = loadlibrary(fullfile(pwd,'libSpinW.so'), 'Hspinw.h');
    assert(isempty(notfound), 'Could not load test library')
%     c = onCleanup(@() unloadlibrary('libSpinW'));
end

s = sw_model('triAF',1);

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

% c1 = calllib('libSpinW','sw_basisvector',c_sw,false);
% c1.setdatatype(c1.DataType,3,3)
% if all(reshape(s.basisvector,1,[]) == reshape(c1.Value,1,[]))
%     disp('Success! :-)')
% end
% 
% c2 = calllib('libSpinW','sw_basisvector',c_sw,true);
% c2.setdatatype(c2.DataType,3,3)
% if all(reshape(s.basisvector(true),1,[]) == reshape(c2.Value,1,[]))
%     disp('Norm Success! :-)')
% end


% a = calllib('libSpinW','sw_qscan',libpointer('doublePtr',[2, 3, 0 0 0, 0 0 1, 10]));


opt = libstruct('spinwave_opt');
% opt = calllib('libSpinW','spinwave_opt_ini',opt);
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
t = calllib('libSpinW','sw_spinwave',c_sw,libpointer('doublePtr',[2, 3, 0 0 0, 0 0 1, 10]),opt)

% calllib('libSpinW','destroy_sw',c_sw)
clear variables 
unloadlibrary('libSpinW')