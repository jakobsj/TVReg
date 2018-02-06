%INSTALL Script used to install TVReg on the Windows/MAC/Linux platform.
%
% Compiles the mex files for the TVReg package.
%
% For Linux, tested with the gcc compiler.
%
% See readme.txt for further instructions.


args = '';
com = computer;

if strcmp(com(1:5), 'PCWIN')
    % If you want to be able to break the execution of the programs, try to
    % set CTRLCBREAK = 1, which uses a non-documented MATLAB API.
    % If you do not have libut.lib in your Matlab dir, try CTRLCBREAK = 2.
    % Default.
    %
    CTRLCBREAK=0;
    
    if CTRLCBREAK==0
        args = [args ''];
    elseif CTRLCBREAK==1
        args = [args '-DLIBUT -L"' matlabroot 
                  '\extern\lib\win32\lcc" -llibut'];
    elseif CTRLCBREAK==2
        args = [args '-DLIBUT -Lexternlib -llibut'];
    else
        error('Not a valid option for CTRLCBREAK')
    end
else % insert any unix compile args here
    args = args;
end

if strcmp(com(end-1:end), '64')
    args = [args ' -largeArrayDims'];
else
    args = [args ' -compatibleArrayDims'];
end


any_error = false;

try
    cs = sprintf('mex %s c/tvreg_upn_c.c c/tools.c c/tv_core.c ', args);
    eval(cs)
catch ME
    any_error= true;
    rethrow(ME);
    
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use tvreg_upn because the compilation failed.')
    disp('You may not have a working compiler installed or mex configured correctly.')
    disp('Please check readme.txt for directions.')
end


try
    cs = sprintf('mex %s c/tvreg_gpbb_c.c c/tools.c c/tv_core.c ', args);
    eval(cs)
catch ME
    any_error= true;
    rethrow(ME);
    
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use tvreg_gpbb because the compilation failed.')
    disp('You may not have a working compiler installed or mex configured correctly.')
    disp('Please check readme.txt for directions.')
end

if any_error == false
    disp('Install completed successfully.')
else
    disp('Installation did not complete successfully.')
end
