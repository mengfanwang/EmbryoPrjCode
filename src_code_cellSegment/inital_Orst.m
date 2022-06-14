function OrSt = inital_Orst(varMap)
% parameters for fg detection and gap significance test

% % temporal variance
%[sigma2, xx1] = temporalVar(vid, vid2, 1);
% use noise variance
%[~, OrSt.stbVarMap, xx] = calVarianceStablizationBY(vid_stb, 0.8);
%[OrSt.NoStbVar, OrSt.NoStbVarMap, xx] = calVarianceStablizationBY(vid, 0.8);
OrSt.gapTestWay = 'localOrdStats';%'orderStats';%
if numel(varMap)==3
    OrSt.stbVarMap = cell(3,1);
    OrSt.NoStbVarMap = cell(3,1);
    OrSt.NoStbVar = nan(3,1);
    for i=1:3
        if ~isempty(varMap{i})
            OrSt.NoStbVarMap{i} = varMap{i}{1,1};
            OrSt.NoStbVar(i) = varMap{i}{2,1};
            OrSt.stbVarMap{i} = varMap{i}{1,2};
        end
    end
else
    OrSt.NoStbVarMap = varMap{1,1};
    OrSt.NoStbVar = varMap{2,1};
    OrSt.stbVarMap = varMap{1,2};
end
% figure;plot(xx1), hold on; plot(xx2);legend('t', 's');
%OrSt.noiseVar = calVarianceStablizationBY(vid_stb, 0.8);
OrSt.mu = [];
OrSt.sigma = [];
OrSt.p_thres = 0.01; % pvalue threshold
OrSt.imProcMethod = 'noStb';% 'stb' or 'noStb'
OrSt.fgVarSamplewiseEst = true;% or 'noStb'
OrSt.refine_cc = 'single_seed'; % use only the connected-component containing the testing seed
OrSt.fgTestWay = 'KSecApprox';%'ttest_varKnown';KSecApprox

if strcmp(OrSt.fgTestWay, 'lookupTable')
    muSigma = paraP3D(0.05, 0,0,255,4, 500, 0.1, 10);
    OrSt.NoStbMu = muSigma.mu;
    OrSt.NoStbSig = muSigma.sigma;
end