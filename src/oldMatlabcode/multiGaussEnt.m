function [ ent ] = multiGaussEnt( prec, N )

%ent = log( ((2*pi*exp(1))^N *1/det(prec) )^0.5 );
ent = 0.5*( N*log(2*pi*exp(1)) - logdet(prec) );
end