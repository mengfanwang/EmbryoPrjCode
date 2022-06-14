function [order_vec, alpha, beta] = order_info(fg, bg, invalid_rg)
% get the order info

if nargin == 2
    invalid_rg = [];
end
M = length(fg);
N = length(bg);
order_vec = zeros(M+N, 2);

n = M+N;
valueAndlabel = zeros(M+N,2);

valueAndlabel(1:N,1) = bg;
valueAndlabel(N+1:n,1) = fg;
valueAndlabel(1:N,2) = -1;
valueAndlabel(N+1:n,2) = 1;

[valueAndlabel, ~] = sortrows(valueAndlabel, 1); %ascending 

order_vec(:,1) = valueAndlabel(:,2);

order_vec(:,2) = valueAndlabel(:,2);
if isempty(invalid_rg)

    labels = valueAndlabel(:,2);
    ratio_v = (1:length(labels))'./length(labels);
    
    %avg_fg_ratio = 1-mean(ratio_v(labels==1));
    bias = 1/length(labels);
    avg_nei_ratio = mean(ratio_v(labels==-1));
    k1 = N/n;
    k2 = avg_nei_ratio;
    dl = 2*k2-k1-bias;
    if dl <= 0 
        beta = 0;
    else
        beta = sqrt(dl/(1-k1));
    end
%     b = -k1;
%     c = k1-2*k2;
%     delta = b*b-4*c;
%     beta = 0;
%     if delta >= 0
%         beta = 0.5*(-b+sqrt(delta));
%         if beta <= 0 || beta >= 1
%             beta = 0.5*(-b-sqrt(delta));
%         end
%         if beta <= 0 || beta >= 1
%             beta = 0;
%         end
%     end
    if beta>1
        beta=1;
    end
    alpha = k1*(1-beta);
    % gamma = 1-alpha-beta;
    invalid_rg = round(alpha*n):round((alpha+beta)*n);
    invalid_rg(invalid_rg<1 | invalid_rg > n) = [];
    order_vec(invalid_rg, 2) = nan;
else
    order_vec(invalid_rg, 2) = nan;
    beta = length(invalid_rg)/size(order_vec, 1);
    k1 = N/n;
    alpha = k1*(1-beta);
    
end
end