function[ par ] = updateF( t, D, par, priors, info )

expect = getExpectations( par );
%tnow = info.time == t;
tind = t+1;

%% update Fprec

if t == 0
    % if t is the first time point in the model
    ppart1 = priors.Fprec{tind};
else
    ppart1 = expect.sig;
end

if t == max(info.tuniq)
    % if t is the last time point in the model
    ppart2 = 0;
else
    ppart2 = par.Smean' * expect.sig * par.Smean;
    % cut that extra row and column out
    %ppart2 = temp(2:end,2:end);
end

if sum(t==info.tuniq)>0
    % if the time point t has data
    ppart3 = par.mship*diag(par.Kmean) * diag(expect.gam) * (par.mship*diag(par.Kmean))';
else
    ppart3 = 0;
end

% 4356
% ppart1
% ppart2
% ppart3
par.Fprec{tind} = ppart1 + ppart2 + ppart3;

%tind
%logdet(par.Fprec{tind})
% if(prod(diag(par.Fprec{tind}))==0)
%     ppart1
%     error('Determinant is zero.');
% end

%% update Fmean (uses current Fprec just calculated)
if t == 0
    % if t is the first time point in the model
    mpart1 = priors.Fprec{tind} * priors.Fmean(:,tind);
else
    %tind
    %blah = expect.S * par.Fmean(:,tind-1)
    %blah = (expect.S * par.Fmean(:,tind-1) + par.Sconstmean)
    mpart1 = expect.sig * (par.Smean * par.Fmean(:,tind-1) + par.Sconstmean);
end

if t == max(info.tuniq)
    % if t is the last time point in the model
    mpart2 = 0;
else
    mpart2 = par.Smean' * expect.sig * (par.Fmean(:,tind+1) - par.Sconstmean);
end

if sum(t==info.tuniq)>0
    % if the time point t has data
    mpart3 = par.mship*diag(par.Kmean) * diag(expect.gam) * par.mumean(tind,:)';
else
    mpart3 = 0;
end

% t
% ppart1\mpart1
% ppart2\mpart2
% ppart3\mpart3
% par.Fprec{tind} 
% (mpart1+mpart2+mpart3)

warning off
par.Fmean(:,tind) = par.Fprec{tind} \ (mpart1+mpart2+mpart3);
warning on

%par.Fmean

end