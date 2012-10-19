function[ par ] = updateMU( t, exper, D, par, priors, info )

expect = getExpectations( par.exper{exper} );

%% one t at a time
tnow = info.time{exper} == t;
tind = info.Tdatind{exper}(info.Tdat{exper}==t);

par.exper{exper}.muprec(tind,:) = expect.gam(:,tind)'*par.mship + ...
    info.N{exper}(t+1) * expect.L;

temp1 = ( expect.gam(:,tind) .* par.exper{exper}.Fmean(:,t+1) )'*par.mship;
temp2 = expect.L * sum( D{exper}(tnow,:),1 );

par.exper{exper}.mumean(tind,:) = (temp1 + temp2) ./ ...
    par.exper{exper}.muprec(tind,:);

end