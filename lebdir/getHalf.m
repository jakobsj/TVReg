function P = getHalf(P)
%GETHALF Discard duplicate vectors differing only by the sign
% 
% P = getHalf(P) takes 2*N times 3 array with unit vectors in rows and
% returns N times 3 array, where half the vectors have been removed. The
% removed vectors are equivalent to the remaining vectors except for their
% sign, which is opposite. 
%
% getHalf is an auxiliary function for getLebedevDirections.

PPT   = P*P';
[I,J] = find( abs(PPT+1) < 1e-10 );
IJ    = [I,J];
sIJ   = sortrows(IJ,1);
P2    = length(P)/2;

list      = zeros(P2,2);
visited   = [];
listcount = 1;

for k= 1:2*P2
    % check if visited
    if any(visited==sIJ(k,1))
        % do nothing
    else
        list(listcount,:) = sIJ(k,:);
        listcount = listcount+1;
        visited = [visited,sIJ(k,:)];
    end
end

P = P(list(:,1),:);
