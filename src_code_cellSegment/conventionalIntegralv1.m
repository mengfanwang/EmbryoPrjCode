function [mu, sigma] = conventionalIntegralv1(interSectV, outerBorderV, order_vec)

%%%%%%%%%%%%
% interSectV: the values of the pixels in the gap
% outerBorderV: the values of the pixels in the neighbor of the gap
% groupid is the vector containing the group information for the ordered
% variables, with 1 for group 1 and -1 for group 2.
% M is the number of elements in group1 and N is the number of
% elements in the group 2
if nargin == 2
    order_vec = [];
end
if isempty(order_vec)
    M = length(outerBorderV);
    N = length(interSectV);
    valueAndlabel = zeros(M+N,2);
    n = M+N;
    valueAndlabel(1:N,1) = interSectV;
    valueAndlabel(N+1:n,1) = outerBorderV;
    valueAndlabel(1:N,2) = -1;
    valueAndlabel(N+1:n,2) = 1;

    valueAndlabel = sortrows(valueAndlabel, 1);
    weights = zeros(1,M+N);
    weights(valueAndlabel(:,2) <0) = -(M+N)/N;
    weights(valueAndlabel(:,2) >0) = (M+N)/M;
else % +1 fg, -1 bg, nan not used
    n = size(order_vec,1);
    weights = zeros(1,n);
    N = sum(order_vec(:,1) < 0); % real number of bg
    M = sum(order_vec(:,1) > 0); % real number of fg
    % not consider overlapped part
    weights(order_vec(:,2) <0) = -n/N;
    weights(order_vec(:,2) >0) = n/M;
end
% n=2000;
% weights = zeros(5000, 1);
% weights(rid(1:1000)) = -1/1000;
% weights(rid(1001:2000)) = 1/1000;

u = [1/((n+1)):1/((n+1)):n/((n+1))];
%u = [1/(10*(n+1)):1/(10*(n+1)):10*n/(10*(n+1))];

v = u;
J = weights;
Fmu = J.*norminv(u);
mu = 1/2*(1/((n+1)))*(2*sum(Fmu) - Fmu(1) - Fmu(end));
%mu = 1/2*(1/(10*(n+1)))*(2*sum(Fmu) - Fmu(1) - Fmu(end));


F1 = J.*u./normpdf(norminv(u)); % Here the normal distribution is used
F2 = J.*(1-v)./normpdf(norminv(v));
F_matrix = F1' * F2;
F_matrix_2 = triu(F_matrix,1);
S = 2*1/4*(1/((n+1)))^2*(F_matrix_2(1,1) + F_matrix_2(end,end) + F_matrix_2(1,end) + F_matrix_2(end,1)...
    + 2*sum(F_matrix_2(2:end-1,1)) +2*sum(F_matrix_2(2:end-1,end)) + 2*sum(F_matrix_2(1,2:end-1)) ...
    + 2*sum(F_matrix_2(end,2:end-1)) + 4*sum(sum(F_matrix_2(2:end-1,2:end-1))));
F_matrix_3 = diag(F_matrix);
S1 = 1/4*(1/((n+1)))^2*(F_matrix_3(1,1) + F_matrix_3(end,end) + F_matrix_3(1,end) + F_matrix_3(end,1)...
    + 2*sum(F_matrix_3(2:end-1,1)) +2*sum(F_matrix_3(2:end-1,end)) + 2*sum(F_matrix_3(1,2:end-1)) ...
    + 2*sum(F_matrix_3(end,2:end-1)) + 4*sum(sum(F_matrix_3(2:end-1,2:end-1))));

sigma = sqrt(S + S1)/sqrt(n);
end


% b = 0.3;
% a = 0.000001;
% out = b*norminv(b)-a*norminv(a)+normpdf(norminv(b))-normpdf(norminv(a))

