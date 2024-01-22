function GpK=KMeans(abd,Y,k) 
% The K-Means algorithm
b=pinv(x2fx(abd))*Y;
b(1)=[];
for i=length(k):-1:1
    GpK(i,:)=kmeans(b,k(i));
    % scoreK(i)=GpScore(expect',GpK(i,:)');
end
end