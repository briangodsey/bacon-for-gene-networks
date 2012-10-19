function[ par ] = updateSmeanHYP( D, par, priors, info )

%expect = getExpectations( par );

SprecMAT = zeros(info.nclus,info.nclus);
for i = 1:info.nclus
    SprecMAT(i,:) = diag(par.Sprec{i})';
end
par.SmeanHYPprec = priors.SmeanHYPprec + sum(sum(SprecMAT));

temp = zeros(info.nclus,info.nclus);
for i = 1:info.nclus
%     % subtract 1 from the diagonal elements of par.Smean because we want
%     % our prior over the transitions matrix to reflect in some way the
%     % identity matrix
%     temp1 = zeros(1,info.nclus);
%     temp1(i) = 1;
    temp(i,:) = diag(par.Sprec{i})' .* par.Smean(i,:);
end
par.SmeanHYPmean = (priors.SmeanHYPprec*priors.SmeanHYPmean ...
    + sum(sum(temp)) ) / par.SmeanHYPprec;

end