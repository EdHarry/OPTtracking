function paths = allPaths(wt, startnode, endnode, n, startN, revisitNodes) % from: http://stackoverflow.com/questions/15524623/all-possible-path-from-starting-node-to-ending-node, EHarry Dec 2014, will make this function use space matricies to optimise processing
wt = logical(sparse(wt));
lastpath = startnode; %We begin with the path containing just the startnode
%costs = 0; %The cost of this path is zero because we haven't yet crossed any edges
paths = []; %The set of solution paths is empty (I'm assuming startnode!=endnode)
N = size(wt,1); %Obtain the number of nodes in the graph
if nargin < 4
    n = N;
end
if nargin < 5
    startN = 2;
end
if nargin < 6
    revisitNodes = false;
end
assert(N==size(wt,2)); %Assert that the adjacency matrix is a square matrix
nDiag = spdiags((1:N)',0,N,N);
for i = startN : n
    sizeLastPath = size(lastpath);
    %Creates a matrix with a row for each path and a 1 in a column where there's a possible move from the last visited node in a path to this column
    nextmove = wt(lastpath(:, i - 1), :);
    
    %Zero out any nodes we've already visited
    d = spdiags((1:sizeLastPath(1))',0,sizeLastPath(1),sizeLastPath(1));
    if ~revisitNodes
        nrows = d * ones(sizeLastPath);
        inds = sub2ind(size(nextmove), reshape(nrows,[],1), reshape(lastpath,[],1));
        nextmove(inds) = false;
    end
    
    % If there are no more available moves we're done
    if ~any(nextmove)
        break;
    end
    
    %For each true entry in our nextmove matrix, create a new path from the old one together with the selected next move
    nextmoverow = d * nextmove;
    nextmovecol = nextmove * nDiag;
    rowlist = reshape(nonzeros(nextmoverow),[],1);
    collist = reshape(nonzeros(nextmovecol),[],1);
    nextpath = [lastpath(rowlist,:), collist];
    
    %     % Compute the costs of the new set of paths by adding the old ones to the cost of each newly traversed edge
    %     inds = sub2ind([N,N],nextpath(:, i-1),nextpath(:,i));
    %     costs = costs(rowlist) + wt(inds);
    
    % For any path finishing on the end node, add it to the return list (and it's corresponding cost)
    %     reachedend = nextpath(:,i) == endnode;
    reachedend = ismember(nextpath(:,i),endnode);
    if any(reachedend)
        %         paths = [paths; {nextpath(reachedend, :)},{costs(reachedend)}];%#ok<AGROW>
        paths = [paths; {nextpath(reachedend, :)}];%#ok<AGROW>
    end
    
    lastpath = nextpath;
    %     %Then remove it from the list of paths still being explored
    %     lastpath = nextpath(~reachedend, :);
    %     %     costs = costs(~reachedend);
    %
    %     %If there are no more paths, we're done
    %     if isempty(lastpath)
    %         break;
    %     end
end
end
