% ROME_SETUP Script to set up the ROME Environment

function rome_setup()
    % Get Path
    str = which(mfilename);        % find full path of current file
    ind = findstr('/', str);       % get indices of / or \
    if(isempty(ind))
        ind = findstr('\', str);
    end
    
    rome_dir = str(1:ind(end)-1);    % get directory name of rome
    addpath(rome_dir, ...
        [rome_dir, '/utilityfuncs']);
        
    % display current ROME version
    rome_version;
    
    % display ok message
    disp('ROME successfully installed!');
end

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
