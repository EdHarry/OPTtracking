function trackPairs = trackPairing2D(tracks2D)
%TRACKPAIRING2D pair-up tracks 
%   EHarry Nov 2014

%% GET PAIRING PARAMETERS
pairingParam = trackPairing_parameters;

% unpack
maxZDifference = pairingParam.maxZDifference;

%% GET TRACK MATRICIES
trackFeats = deal(repmat(struct('mat',[]),2,1));

for i = 1:2
    feat = convStruct2MatIgnoreMS(tracks2D(i).trackInfo);
    if i == 2
        [m,nTimePoints_m] = size(feat);
        [n,nTimePoints_n] = size(trackFeats(1).mat);
        if nTimePoints_m > nTimePoints_n % pad matricies if nessesary to give both series the same number of timepoints
            trackFeats(1).mat = [trackFeats(1).mat NaN(n,(nTimePoints_m - nTimePoints_n))];
        elseif nTimePoints_m < nTimePoints_n
            feat = [feat NaN(m,(nTimePoints_n - nTimePoints_m))];%#ok<AGROW>
        end
    end
    trackFeats(i).mat = feat;
end

%% CREATE COST MAT
costMat = -ones(n,m);

%% MAIN LOOP
for iTrack = 1:n
    zN = getZ(1,iTrack);
    for jTrack = 1:m
        zM = getZ(2,jTrack);
        zDiff = abs(zN - zM);
        if nanmax(zDiff) < maxZDifference % if max Z diff is less than cutoff, consider for linking
            xCorr = crossCorr(diff(zN)',diff(zM)',0); % get cross correlation between rate-of-change of each series
            xCorr = xCorr(1);
            if xCorr > 0 % if correlation is positive, add inverse to cost matrix
                costMat(iTrack,jTrack) = 1 / xCorr;
            end
        end
    end
end

% end if no possible links
if all(costMat(:) == -1)
    disp('no pairs found');
    trackPairs = [];
    return
end

%% LINEAR ASSIGNMENT
links = lap(costMat,[],[],true); % solve lap
links = double(links(1:n));
goodLinks = links <= m;
links = [find(goodLinks) links(goodLinks)];

%% EXTRA LINKS
unlinked = cell(2,1);
unlinked{1} = find(~goodLinks); 
unlinked{2} = setdiff((1:m)',links(:,2));

for i = 1:2
    j = mod(i,2) + 1;
    for link = links(:,i)'
        switch i
            case 1
                costMatSub = costMat(link,:)';
            case 2
                costMatSub = costMat(:,link);
        end
        linksSub = find(costMatSub ~= -1);
        linksSub = linksSub(ismember(linksSub, unlinked{j}));
        if ~isempty(linksSub)
            switch i
                case 1
                    linksSub = [link*ones(size(linksSub)) linksSub];%#ok<AGROW>
                case 2
                    linksSub = [linksSub link*ones(size(linksSub))];%#ok<AGROW>
            end
            links = [links; linksSub];%#ok<AGROW>
        end
    end
end

%% RETURN
trackPairs = links;

%% SUBFUNCTIONS

    function z = getZ(iSeries,iTrack) % returns Z coords, the 2nd coordinate ("y") in both cases
        z = trackFeats(iSeries).mat(iTrack,2:8:end);
    end

end

