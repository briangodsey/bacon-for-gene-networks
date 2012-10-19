function[ par ] = updateSmeanHYPdoublePREC( D, par, priors, info )

par.SmeanHYPdoublePRECshape = priors.SmeanHYPdoublePRECshape ...
    + 0.5 + 0.5*info.nclus^2;

expSmeanMIDsq = sum(sum(par.SmeanMiddleGuymean.^2 + 1./par.SmeanMiddleGuyprec)); %matrix sum
expSmeanMIDxSmeanHYP = sum(sum(par.SmeanMiddleGuymean * par.SmeanHYPmean)); %matrix sum
expSmeanHYPsq = par.SmeanHYPmean^2 + 1/par.SmeanHYPprec; %scalar
expSmeanHYPxSmeanHYPprior = par.SmeanHYPmean * priors.SmeanHYPmean; %scalar
expSmeanHYPpriorSq = priors.SmeanHYPmean^2; %scalar

par.SmeanHYPdoublePRECscale = priors.SmeanHYPdoublePRECscale...
    + 0.5*( expSmeanHYPsq - 2*expSmeanHYPxSmeanHYPprior + expSmeanHYPpriorSq ) ...
    + 0.5*( expSmeanMIDsq - 2*expSmeanMIDxSmeanHYP + info.nclus^2*expSmeanHYPsq );

precLimit = 10^-3;
if( par.SmeanHYPdoublePRECshape / par.SmeanHYPdoublePRECscale < precLimit )
    par.SmeanHYPdoublePRECscale = par.SmeanHYPdoublePRECshape / precLimit;
end

end
