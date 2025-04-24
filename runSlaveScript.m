function runSlaveScript()
% runSlaveScript - Executes the BladeRF slave batch script.
%
% Assumes runSlave.bat is in the MATLAB path or working directory.

    status = system('runSlave.bat');
    if status ~= 0
        error('Failed to execute runSlave.bat. Check file path or permissions.');
    end
end