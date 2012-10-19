function [ info ] = setinfo( D, timeORIG, names )

info.G = size(D{1},2);  % number of genes
%%info.slides = size(D{1},1);

info.timeORIG = timeORIG;
%%info.exper = exper;
info.names = names;
%%info.expnames = info.names( info.exper );
info.timeunit = min( info.timeORIG{1}( info.timeORIG{1} ~= 0 ) );

info.time = {};
info.Tmax = {};
info.Tdat = {};
info.Tdatind = {};
info.Tall = {};
info.timeD = {};
info.N = {};
for i = 1:size(info.timeORIG,2)
  info.time{i} = info.timeORIG{i} / info.timeunit;
  %% info.tuniq = unique(info.time);
  info.Tmax{i} = max(info.time{i});
  info.Tdat{i} = sort( unique(info.time{i}) );
  info.Tdatind{i} = 1:length(info.Tdat{i});
  info.Tall{i} = min(info.time{i}):max(info.time{i});

  for j = 1:length(info.time{i})
    info.timeD{i}(j) = find(info.Tdat{i}==info.time{i}(j));
  end
  
  if( round(info.time{i}) ~= info.time{i} )
    error('Vector info.time: smallest element is not a divisor of all other elements.\n');
  end

  info.N{i} = repmat(0,1,info.Tmax{i});
  for j = info.Tall{i}+1
    info.N{i}(j) = sum(info.time{i}+1==j);
  end
  
end


end