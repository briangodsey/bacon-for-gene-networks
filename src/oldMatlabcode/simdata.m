
%% create data
nclus = 5;
numreps = 3;
genesperclus = 10;
timepoints = 8;

clusternumbers = nclus;

Smat = 1.5*(rand(nclus)-0.5);

Sconst = 10*normrnd(0,0.0001,nclus,1);
time0 = normrnd(0,10,nclus,1);

preD = 10*time0;
for i = 2:timepoints
    preD = [ preD Smat*preD(:,i-1)+Sconst ];
end

preD = preD + normrnd(0,0.01,size(preD));

% %%%% use special data set below %%%%%%%
% Smat = Smat2;
% preD = preD2;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mus = [];
for i = 1:size(preD,1)
    mus = [ mus repmat(preD(i,:)', 1 ,genesperclus) ];
end
mus = mus + normrnd( 0, 0.01, size(mus) );

% D = [];
% for i = 1:size(preD,1)
%     D = [ D repmat(preD(i,:)',numreps,genesperclus) ];
% end
D = repmat( mus, numreps, 1 );
D = D + normrnd( 0, 0.001, size(D) );

% simulate deltas
deltas = repmat(0,size(D));
%deltas(1:8,:) = 2;
%deltas(9:16,:) = -2;
D = D+deltas;

%% special data sets

% set 1, 10 clusters; the dynamics get worse through time (compare
% postD to par.Fmean or preD
Smat1 = ...
[    0.0351    0.4349    0.3670    0.4414    0.4784    0.4023   -0.2469    0.2642    0.0496   -0.3397; ...
   -0.0283   -0.4042    0.4471   -0.1038    0.3331    0.1454   -0.3899   -0.4802    0.3476   -0.1633; ...
    0.4600    0.3939    0.2304   -0.2684   -0.2716   -0.0703   -0.0256    0.1729   -0.2043   -0.1857; ...
    0.1586   -0.0410    0.1437   -0.0033    0.4319    0.1259   -0.2898   -0.3405   -0.4757   -0.2708; ... 
    0.1839    0.4006    0.3414    0.1811    0.1608    0.1569    0.2354    0.0424    0.4554    0.3517; ...
    0.2951    0.0165   -0.2112    0.3566    0.3605   -0.3232   -0.2707   -0.4715    0.4606    0.1098; ...
    0.3123    0.0633    0.4442    0.3238   -0.1878    0.3617   -0.4033   -0.2716   -0.1534    0.4278; ...
    0.2870   -0.4993    0.3234   -0.4952   -0.4164    0.1788   -0.0761    0.2740    0.3416    0.3799; ...
   -0.1843    0.3398   -0.1134    0.1942   -0.0179   -0.1921   -0.0205   -0.0747    0.4764   -0.0145; ...
    0.0874   -0.4898   -0.0222   -0.1127    0.4414   -0.3549   -0.2240   -0.4489    0.1480    0.4454 ];

preD1 = ...
[  110.3366 -129.6083   33.1487  -60.8949   14.7966  -16.4128   -0.2649   15.2718; ...
   18.0974   15.6739  -35.7059   34.2149   15.5344  -35.6780   65.2258  -62.4407; ...
  106.6276  142.0443  -14.6211    2.8393  -32.0747   -2.5387  -24.0821   13.2291; ...
  -82.3405  105.6472   -8.1256   33.6618   -9.9642  -27.1207   31.6663  -30.9777; ...
 -115.8992  -50.6839    9.7213  -48.8548   18.8294   26.6792  -32.2997   36.2133; ...
  -46.5387  -52.1200  -94.6764   64.2252    1.2589    2.0008   23.9624   -4.7284; ...
  102.5072   77.4401   29.9663  -58.0099   86.4291  -57.5902   -1.9651   36.6561; ...
 -208.5974   -0.1180  -77.5465  -44.4677   -2.7874  -19.9638   12.4444  -18.8571; ...
 -154.3653  -90.5468   -0.1454    5.9737   24.1458   13.7852   -5.4783   22.8856; ...
  -53.3376   -2.8441  -69.9099   56.4733  -11.3025  -16.1813   46.8961  -42.8189 ];
  
% set 2, 7 clusters; everything is right, but par.Sprec and expect.sig are
% low
Smat2 = ...
   [ 0.2255   -0.0764   -0.4031   -0.3559    0.0325   -0.1977   -0.4215; ...
    0.4436   -0.4885   -0.4198    0.4515    0.2075   -0.1496    0.1158; ...
   -0.0191    0.1709   -0.2189    0.2263   -0.0921   -0.2782   -0.4462; ...
   -0.1154   -0.0170   -0.1733    0.4462   -0.0559   -0.1915    0.3398; ...
    0.3584    0.2415    0.3423   -0.3604    0.0100   -0.0554    0.0406; ...
    0.2220    0.2557   -0.3874   -0.4308    0.4077   -0.3295    0.2935; ...
   -0.4257   -0.4933   -0.3517    0.0878   -0.0729   -0.2532   -0.0037 ];
   
preD2 = ...
 [  -83.2041  -33.8096  -15.5519  -37.0893  -11.9585  -30.1210  -11.2570  -23.1852; ...
 -130.9442   37.4523  -11.5495   14.1032   -6.1913    2.2624    0.5378   -0.4963; ...
   37.2700  -46.1962    1.6486  -10.5471   19.4160   -5.9633   11.0878   -3.3422; ...
   17.2915   -0.3348   45.7768   25.7660   28.5462   22.1443   21.5477   18.8721; ...
  114.5245  -55.8084  -14.2694  -25.1016  -20.7399   -8.1418  -18.8192   -6.7866; ...
   37.3488  -39.3209   31.0372  -35.2134   -7.6730  -22.4878   -7.1666  -14.7915; ...
    0.4509   70.6216   25.8988    8.8364   25.5098    7.1847   22.0233    5.6212 ];
    
%% create info parameter
% info.G = size(D,2);  % number of genes
% info.slides = size(D,1);
% info.nclus = nclus;
% 
% info.exper = repmat( 1, 1, size(D,1) );
% info.names = { 'simulation' };
% info.expnames = info.names( info.exper );
% info.timeunit = min( info.timeORIG( info.timeORIG ~= 0 ) );
% info.time = info.timeORIG / info.timeunit;
% info.tuniq = unique(info.time);
% info.Tmax = max(info.time);
% %info.Tdat = size( unique(info.timeORIG), 1 );
% %info.mutime(info.tuniq) = 1:length(info.tuniq);
% 
% if( round(info.time) ~= info.time )
%     error('Vector info.time: smallest element is not a divisor of all other elements.\n');
% end
% 
% info.N = repmat(0,1,info.Tmax);
% for i = 1:info.Tmax+1
%     info.N(i) = sum(info.time+1==i);
% end

%% now the priors
% 
% smallnum = 10^-8;
% largenum = 10^8;
% 
% priors.mship = repmat(1/info.nclus, info.nclus, info.G);
% %priors.mship(1,:) = 0;
% %priors.mship = [ 1 1 1 0 0 0; 0 0 0 1 1 1];
% 
% priors.mumean = repmat(0,info.Tmax+1,info.G);
% %priors.mumean = mus;
% priors.muprec = repmat(smallnum,info.Tmax+1,info.G);
% 
% priors.deltamean = repmat( 0, 1, size(D,1) );
% priors.deltaprec = repmat( smallnum, 1, size(D,1) );
% 
% priors.gamshape = repmat( smallnum, 1, info.G );
% priors.gamscale = repmat( smallnum, 1, info.G );
% 
% priors.Lshape = smallnum;
% priors.Lscale = smallnum;
% 
% priors.Kmean = repmat( 1, 1, info.G );
% priors.Kprec = repmat( largenum, 1, info.G );
% 
% priors.Smean = eye( info.nclus );
% %priors.Smean = Smat;
% %priors.Smean(1,1) = 0;
% priors.Sprec = {};
% for i = 1:info.nclus
%     priors.Sprec{i} = eye(info.nclus)*smallnum;
% end
% 
% priors.Sconstmean = repmat(0,info.nclus,1);
% priors.Sconstprec = repmat(smallnum,info.nclus,1);
% 
% % for Fmean, add T+1 for convenience; it stays zero
% priors.Fmean = repmat( [1:info.nclus]', 1, info.Tmax + 1) - info.nclus/2;
% %priors.Fmean(:,1) = [ 0.7; 1.5 ];
% %priors.Fmean = preD;
% %priors.Fmean = [ 0 0 0; 1 1 1 ];
% priors.Fprec = {};
% for t = 1:info.Tmax+1
%     priors.Fprec{t} = eye( info.nclus )*smallnum;
% end
% 
% priors.SIGshape = repmat(smallnum,1,info.nclus);
% priors.SIGscale = repmat(smallnum,1,info.nclus);

%% par
%