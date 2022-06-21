% function [his,variance,parameters] = histogramCount(data,kerSize,options)

    % This function is used to calculate variance of the input data with
    % the memory constriant. It can calculate each image seperately and 
    % then give an overall variance curve.
    % The input data should be 4D. Otherwise, there are better way to do.
    
clc;clear;close all;

for tt = 0:49
    tt_ind = num2str(100+tt);
    tt_ind = tt_ind(2:3);
    load(['H:\Embryo\TM0-49\multiview\' tt_ind '.mat']);
    
    kerSize = [3 3 3];
    histEdges = 0.5:254.5;
    binEdges = -0.5:255.5;
    
    realSignal = zeros(dataSize);
    realSignal(:,:,:,ii) = medfilt3(data(:,:,:,ii),kerSize);
    kerRadi = floor(kerSize/2);
    data = data(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
        kerRadi(3)+1:end-kerRadi(3));
    realSignal = realSignal(kerRadi(1)+1:end-kerRadi(1),kerRadi(2)+1:end-kerRadi(2),...
        kerRadi(3)+1:end-kerRadi(3));
end