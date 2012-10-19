function[ par ] = updateH( t, D, par, priors, info )

% OBSOLETE!! NOT USED ANY MORE!!

% use this function for time points without data
% the parameters are still "F", but the function is "H"

expect = getExpectations( par );

%tnow = info.time == t;
tind = t+1;

par.Fprec{tind} = priors.Fprec{tind} + expect.S' * par.Fprec{tind+1} * expect.S;

if t == min(info.tuniq)
    temp1 = 0;
else
    temp1 = par.Fprec{tind} * par.Smean * par.Fmean(:,tind-1);
end

if t == max(info.tuniq)
    temp2 = 0;
else
    temp2 = par.Smean' * par.Fprec{tind+1} * par.Fmean(:,tind+1);
end


temp2 = par.Smean' * par.Fprec{tind+1} * par.Fmean(:,tind+1);

par.Fmean(:,tind) = par.Fprec{tind} \ (temp1+temp2);

end