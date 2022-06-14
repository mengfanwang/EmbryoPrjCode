function [edgeCost, stdEstFrom2Cur] = distance2EdgeCost_cell(neiVox, curVox,  curStd, neiStd, curRatioSquare,g, timeGap)
if nargin<7
    timeGap = 1;
end
% use chi-square
df = 0;
%chiSquTestStat = 0;
if g.stdCombined
    stdEstFrom2Cur = sqrt(curRatioSquare(:,1)'.*curStd.^2 + curRatioSquare(:,2)'.*neiStd.^2);
else
    stdEstFrom2Cur = curStd;
end
if 0 % way 1 strcmp(g.varEstMethod,'independent')
    if ~isempty(uNei) && isempty(find(isnan(neiStd) | neiStd==0, 1))
        chiSquTestStat = sum((ovDistanceRegion(curVox, neiVox)./neiStd).^2) ;
        df = df+length(neiStd);
    end
    if ~isempty(curVox)
        chiSquTestStat = chiSquTestStat + (ovDistanceRegion(curVox, neiVox)./curStd).^2);
        df = df+length(neiStd);
    end
else
    muXYZ = ovDistanceRegion(curVox, neiVox);
    % problem: variance use which one
    if ~isempty(muXYZ)
        chiSquTestStat = (muXYZ./stdEstFrom2Cur).^2;
        df = df+length(stdEstFrom2Cur);
        pValue = 1-chi2cdf(chiSquTestStat, df);
    else
        error('No muXYZ exist!!!\n');
    end
end
if isfield(g,'jumpCost')
    if length(g.jumpCost)>1
        pValue = pValue*g.jumpCost(timeGap);
        pValue = sum(min([g.jumpCost; g.jumpCost*0+pValue])); % why?
        edgeCost = pvalue2edgeCost(pValue, g);
    elseif length(g.jumpCost)==1
        edgeCost = pvalue2edgeCost(pValue, g);
        edgeCost = edgeCost+(timeGap-1)*g.jumpCost;
    else
        edgeCost = pvalue2edgeCost(pValue, g);
    end
else
    edgeCost = pvalue2edgeCost(pValue, g);
end
end