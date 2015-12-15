% script setup.m: Install and add to path.
p = pwd;
w = [p filesep 'DWT'];
disp(['Adding directories to path and saving path...']);
addpath(w);
addpath([p filesep 'Graphics']);
addpath([p filesep 'Misc']);
addpath([p filesep 'Polarization']);
savepath;

disp(['Checking for compiled .mex files...']);
if ismac
    disp('Apple computers lack a native C compiler.');
    disp(['To use wavelet subroutines, please compile .mex files in ' ...
        w ' manually.']);
    return;
end
if ispc
    suff = 'w';
else
    suff = 'g';
end
try
    f = dir([w filesep '*mex' suff '*']);
catch
    disp('DWT subdirectory not found; wavelet routines won''t work.');
    return;
end

if size(f,1) < 3
    err = 0;
    cd(w);
    disp(['Trying to compile .mex files for wavelet subroutines.']);
    try
        disp(['Building dwptpyr (DWT, DWPT)...']);
        mex dwptpyr.c;
    catch
        err = 1;
        disp(['dwptpyr.c failed to compile.']);
    end
    try
        disp(['Building modwptpy (MODWT, MODWPT)...']);
        mex modwptpy.c;
    catch
        err = 1;
        disp(['modwptpy.c failed to compile.']);
    end
    try
        disp(['Building imodwptp (inverse MODWT, inverse MODWPT)...']);
        mex imodwptp.c;
    catch
        err = 1;
        disp(['imodwptp.c failed to compile.']);
    end
    if err == 1
        disp('Wavelet routines error! Some or all wavelet routines will not work.');
        disp('Check that your system supports C compilation via the mex command.');
        cd(p);
        return
    end
else
    disp(['Compatible .mex files found, nothing to do.'])
end
cd(p);
disp('Setup complete.');
