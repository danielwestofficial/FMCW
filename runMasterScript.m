function runMasterScript()
% runMasterScript - Executes the BladeRF master batch script.
%
% Assumes runMaster.bat is in the MATLAB path or working directory.

    status = system('runMaster.bat');
    if status ~= 0
        error('Failed to execute runMaster.bat. Check file path or permissions.');
    end
end