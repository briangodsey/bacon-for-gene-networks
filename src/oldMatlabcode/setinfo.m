function [ info ] = setinfo( D, timeORIG, exper, names )

info.G = size(D,2);  % number of genes
info.slides = size(D,1);

info.timeORIG = timeORIG;
info.exper = exper;
info.names = names;
info.expnames = info.names( info.exper );
info.timeunit = min( info.timeORIG( info.timeORIG ~= 0 ) );
info.time = info.timeORIG / info.timeunit;
info.tuniq = unique(info.time);
info.Tmax = max(info.time);
%info.Tdat = size( unique(info.timeORIG), 1 );
%info.mutime(info.tuniq) = 1:length(info.tuniq);

if( round(info.time) ~= info.time )
    error('Vector info.time: smallest element is not a divisor of all other elements.\n');
end

info.N = repmat(0,1,info.Tmax);
for i = 1:info.Tmax+1
    info.N(i) = sum(info.time+1==i);
end

end