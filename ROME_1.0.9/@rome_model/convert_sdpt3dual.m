% ROME_MODEL\convert_sdpt3dual Converts the model object into sdpt3 dual
% solvable.
%
% 1. Jingxi April 24 09

function [blk,At,C,b, add_obj]=convert_sdpt3dual(model_obj)


%initilization
blk_ind=1;
b=[];

if(isa(model_obj.ObjFn, 'rome_var'))
    % zero pad to include aux model variables.
    b = model_obj.ObjFn.BiAffineMap(1, 2:end).';
    b = [spalloc(model_obj.ObjFn.NumUnmappedVars, 1, 0); b];
    b = [b; spalloc(model_obj.NumVars - length(b), 1,0)];
    add_obj = model_obj.ObjFn.BiAffineMap(1, 1);
    if(model_obj.MinMaxFlag == rome_model.MINIMIZE)
        b = -b;
    end
end

%Aqt

if(~isempty(model_obj.QC))
    QC=model_obj.QC;
    blk{blk_ind,1}='q';
    blk{blk_ind,2}=[];
    C{blk_ind,1}=[];
    ind=[];
    for i=1:size(QC,1)
        sub_blk_size=1-QC{i}(1)+QC{i}(2);
        blk{blk_ind,2}=[blk{blk_ind,2}, sub_blk_size];
        C{blk_ind,1}=[C{blk_ind,1};spalloc(sub_blk_size,1,0)];
        ind=[ind,QC{i}(1):1:QC{i}(2)];
    end
    
    At{blk_ind}=sparse(1:size(C{blk_ind,1},1),ind,-1);
    num_append=model_obj.NumVars-size( At{blk_ind},2);
    if (num_append>=1)
        At{blk_ind}=[At{blk_ind} spalloc(size( At{blk_ind},1),num_append,0)];
    end    
    blk_ind=blk_ind+1;
end

if(~isempty(model_obj.LC))
    Al=-model_obj.LC.BiAffineMap;
    %resize Al
    if (model_obj.LC.NumUnmappedVars>=1)        
        Al=[Al(:,1),spalloc(size(Al,1),model_obj.LC.NumUnmappedVars,0),Al(:,2:end)];        
    end
    post_zeros=model_obj.NumVars-model_obj.LC.NumUnmappedVars-model_obj.LC.NumMappedVars;
    if(post_zeros>=1)
        Al=[Al,spalloc(size(Al,1),post_zeros,0)];
    end
    
    %Aut (equality constraints)
    if( ~isempty(model_obj.IndEq))        
        Au=Al(model_obj.IndEq,:);        
        blk{blk_ind,1}='u';
        blk{blk_ind,2}=size(Au,1);
        At{blk_ind}=Au(:,2:end);
        RHS=-Al(:,1);        
        C{blk_ind,1}=RHS(model_obj.IndEq);
        % update Al and c
        
        index=setdiff([1:size(Al,1)]',model_obj.IndEq);        
        Al=Al(index,:);        
        blk_ind=blk_ind+1;
    end
    
    %Alt
    C{blk_ind,1}=-Al(:,1);
    At{blk_ind}=Al(:,2:end);    
else
    At{blk_ind}=[];
    C{blk_ind,1}=[];
end

%first deal with ub and lb
LB = model_obj.LB;
UB = model_obj.UB;
Al_LB=[];
Al_UB=[];
blk{blk_ind,1}='l';


if(~isempty(LB))
    Al_LB=sparse(1:size(LB,1),1:size(LB,1),-1);
    num_append_LB = model_obj.NumVars - length(LB);
    
    if(num_append_LB > 0)
        Al_LB = [Al_LB, spalloc(size(LB,1),num_append_LB,0)];
    end
    index=find(~isinf(LB));
    LB=LB(index,:);
    Al_LB= Al_LB(index,:);
    
    At{blk_ind}=[At{blk_ind};Al_LB];
    C{blk_ind,1}=[C{blk_ind,1};-LB];    
end

if(~isempty(UB))
    Al_UB=sparse(1:size(UB,1),1:size(UB,1),1);
    num_append_UB = model_obj.NumVars - length(UB);
    
    if(num_append_UB > 0)
        Al_UB = [Al_UB, spalloc(size(UB,1),num_append_UB,0)];
    end
    index=find(~isinf(UB));
    UB=UB(index,:);
    Al_UB= Al_UB(index,:);
    
    At{blk_ind}=[At{blk_ind};Al_UB];
    C{blk_ind,1}=[C{blk_ind,1};UB];    
end

blk{blk_ind,2}=size(At{blk_ind},1);

% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
