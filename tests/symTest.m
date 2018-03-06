
% Try to load the library.
if not(libisloaded('libSpinW'))
    try
    [notfound, warnings] = loadlibrary(...
        fullfile(pwd,'..','cmake-build-debug-gcc_7','libSpinW.so'),...
        fullfile(pwd,'..','include','Hspinw.h'));
    assert(isempty(notfound), 'Could not load test library')
    catch ME
        disp(ME.message)
        exit;
    end
end

sym_path = '/home/ward_s/Development/SpinW/spinw/dat_files/';

sym = calllib('libSpinW','loadsym',sym_path);

exit()