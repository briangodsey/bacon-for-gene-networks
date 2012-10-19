function [priors] = setpriors( D, info )

smallnum = 10^-10;
largenum = 10^10;

priors.mship = repmat(1/info.nclus, info.nclus, info.G);
%priors.mship(1,:) = 0;
%priors.mship = [ 1 1 1 0 0 0; 0 0 0 1 1 1];

priors.mumean = repmat(0,info.Tmax+1,info.G);
%priors.mumean = mus;
priors.muprec = repmat(smallnum,info.Tmax+1,info.G);

priors.deltamean = repmat( 1, 1, size(D,1) );
priors.deltaprec = repmat( smallnum, 1, size(D,1) );

priors.gamshape = repmat( smallnum, 1, info.G );
priors.gamscale = repmat( smallnum, 1, info.G );

priors.Lshape = smallnum;
priors.Lscale = smallnum;

%priors.Kmean = repmat( 1, 1, info.G );
%priors.Kprec = repmat( largenum, 1, info.G );
priors.Kmean = repmat( 1, 1, 1 );
priors.Kprec = repmat( largenum, 1, 1 );

priors.Smean = zeros( info.nclus );
%priors.Smean = Smat;
%priors.Smean(1,1) = 0;
priors.Sprec = {};
for i = 1:info.nclus
    priors.Sprec{i} = eye(info.nclus)*smallnum;
end

priors.Sconstmean = repmat(0,info.nclus,1);
priors.Sconstprec = repmat(smallnum,info.nclus,1);

% for Fmean, add T+1 for convenience; it stays zero
%priors.Fmean = repmat( [1:info.nclus]', 1, info.Tmax + 1) - info.nclus/2;
priors.Fmean = zeros(info.nclus, info.Tmax + 1);
priors.Fprec = {};
for t = 1:info.Tmax+1
    priors.Fprec{t} = eye( info.nclus )*smallnum;
end

priors.SIGshape = repmat(smallnum,1,info.nclus);
priors.SIGscale = repmat(smallnum,1,info.nclus);

%% old
% smallnum = 10^-5;
% 
% priors.mship = repmat(1/info.nclus, info.nclus, info.G);
% 
% priors.mumean = repmat(0,info.Tmax+1,info.G);
% priors.muprec = repmat(smallnum,info.Tmax+1,info.G);
% 
% priors.deltamean = repmat( 0, 1, size(D,1) );
% priors.deltaprec = repmat( 1, 1, size(D,1) );
% 
% priors.gamshape = repmat( smallnum, 1, info.G );
% priors.gamscale = repmat( smallnum, 1, info.G );
% 
% priors.Lshape = smallnum;
% priors.Lscale = smallnum;
% 
% priors.Kmean = repmat( 1, 1, info.G );
% priors.Kprec = repmat( smallnum, 1, info.G );
% 
% priors.Smean = eye( info.nclus+1 );
% for i = 1:info.nclus+1
%     priors.Sprec{i} = eye(info.nclus+1);
% end
% 
% % for Fmean, add T+1 for convenience; it stays zero
% priors.Fmean = repmat( [1:info.nclus]', 1, info.Tmax + 1) - info.nclus/2; 
% priors.Fmean(1,:) = 1;
% for t = 1:info.Tmax+1
%     priors.Fprec{t} = eye( info.nclus );
% end

end