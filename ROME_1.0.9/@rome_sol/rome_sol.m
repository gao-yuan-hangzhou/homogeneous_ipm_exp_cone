% ROME_SOL Describes a solution object in ROME. 
% a rome_sol object is output AFTER optimization, and is used for
% instantiating with uncertianties.
% 
% Data components:
% 1. LDRAffineMap: LDR component
% 2. DeflectCoeff : Coefficients of the deflected components
% 3. DeflectVals  : Values of the deflected components
% 
% Note that we assume data to be stored in the ()^- negative part 
% 
% Modification History: 
% 1. Joel 
% 

classdef rome_sol
    properties
        LDRAffineMap = [];        % linear part of the solution in an affine map format
        DeflectCoeff  = [];       % Coefficients of the deflected components
        DeflectAffineMap   = [];  % Values of the deflected components
        Size = []  ;              % actual size of the solution
    end    
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
