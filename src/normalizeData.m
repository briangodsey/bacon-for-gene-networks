function nD = normalizeData( D, n )

nD = D;

if ( size(D,1)/n ~= round(size(D,1)/n) )
   disp('Warning: The number of rows is not divisble by the n parameter value'); 
end

%% if ( n == 0 )
%%     n = size(D,1);
%% end

for i = 1:floor( size(D,1)/n )
    Dsub = D( ((i-1)*n+1):(i*n), : );
    Dmeans = repmat( mean(Dsub,1), size(Dsub,1), 1 );
    residuals = Dsub - Dmeans;
    noise = repmat( std(residuals,0,1), size(Dsub,1), 1 );
    nD( ((i-1)*n+1):(i*n), : ) = residuals ./ noise;
end

end