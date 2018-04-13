
% Try to load the library.
if not(libisloaded('libSpinW'))
    try
    [notfound, warnings] = loadlibrary(...
        fullfile(pwd,'..','cmake-build-debug-gcc_7','libSpinW.so'),...
        fullfile(pwd,'..','include','Hspinw.h'));
    assert(isempty(notfound), 'Could not load test library')
    catch ME
        disp(ME.message)
        exit(1);
    end
end

sym_path = fullfile(sw_rootdir,'dat_files',filesep);

try
    sym = calllib('libSpinW','loadsym',sym_path);
    disp('Success')
    exit(0)
catch ME
    disp(ME.message)
    exit(1)
end