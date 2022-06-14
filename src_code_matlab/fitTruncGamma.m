function phat = fitTruncGamma(data)


%!!!: let's see if we directly use gamma
% phat = gamfit(data);
% return
disp_gap = 0.03;
% data = clsstSftOvrlpDist;
data(data==0) = [];
% figure; histogram(data,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);


phatOri = gamfit(data);
% temp = gamrnd(phatOri(1),phatOri(2),length(data),1);
% hold on; histogram(temp,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);
% title(num2str(phatOri));
% hold off;

phat0 = phatOri;
truncThr0 = max(data)*2;
truncThr = max(data);

tol = quantile(data,0.99)*0.01;

while (truncThr0 - truncThr) > tol
    truncThr0 = truncThr;
%     disp(truncThr0);
    
    truncThr = gaminv(1-0.05, phat0(1), phat0(2));
    truncedVec = data(data <= truncThr);
%     figure; histogram(truncedVec,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);
    
    pdf_truncgamma = @(x,a,b) gampdf(x,a,b) ./ gamcdf(truncThr,a,b);
    %[phat,~] = mle(truncedVec, 'pdf',pdf_truncgamma, 'start',phat0);
    phat = mleNewtonMethod(pdf_truncgamma, truncedVec, phat0);
%     temp = gamrnd(phat(1),phat(2),length(data),1);
%     hold on; histogram(temp,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);
%     hold off;
    phat0 = phat;
end
fprintf('----------Fitted Gamma: %f, %f-------------\n', phat(1), phat(2));
% figure; histogram(truncedVec,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);
% temp = gamrnd(phat(1),phat(2),length(data),1);
% hold on; histogram(temp,'BinWidth',disp_gap); xlim([min(data),quantile(data,0.99)]);
% title(num2str(phat));
% hold off;






