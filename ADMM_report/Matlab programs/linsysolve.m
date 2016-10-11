%%**********************************************************************
    function q = linsysolve(L,r) 

    if isfield(L,'perm')
       if strcmp(L.matfct,'chol')
          q(L.perm,1) = mextriang(L.R, mextriang(L.R,r(L.perm),2) ,1);
       elseif strcmp(L.matfct,'spchol')
          q(L.perm,1) = mexbwsolve(L.Rt,mexfwsolve(L.R,r(L.perm)));
       end
    else
       if strcmp(L.matfct,'chol')
          q = mextriang(L.R, mextriang(L.R,r,2) ,1);
       elseif strcmp(L.matfct,'spchol')
          q = mexbwsolve(L.Rt,mexfwsolve(L.R,r));
       end        
    end
%%**********************************************************************