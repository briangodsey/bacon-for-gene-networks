function[ intsigs Ssigs ] = getInteractionSignificance( Smat, Sprec, mship )

  n = size(mship,1);
  g = size(mship,2);

  if length(Sprec) == 1
    Ssigs =  Smat ~= 0;
  else
    Svarmat = zeros(n);
    for i = 1:n
      %%Svarmat(i,:) = diag(pinv(Sprec{i}));
      Svarmat(i,:) = diag(pinv(Sprec{i}));
    end
    %% Ssigs = 2*normcdf( abs(Smat)./(Svarmat.^0.5) )-1;
    Ssigs = 2*normcdf( Smat ./ (Svarmat.^0.5) )-1;

    %% abs(Smat)
    %% Svarmat
    %% abs(Smat)./(Svarmat.^-0.5)
    %% Ssigs

  end

  intsigs = zeros(g);
  for i = 1:g
    for j = 1:g
      %% intsigs(i,j) = mship(:,i)' * Ssigs * mship(:,j);
      intsigs(i,j) = ( mship(:,i)' * Smat * mship(:,j) ) ...
	  / ( mship(:,i)'.^2 * Svarmat * mship(:,j).^2 )^0.5;
    end
  end

end