function p = track_disp_para(numTracks)
% the parameters for drawing the trajectories

p.shiftScale = [50 50 10];
p.orgDataDim = 1.2; % make the data dimmer, so tracks will stick out
p.maxInt = 220;
p.particleSize = 1;
p.particleCl = [1 1 0.1]*p.maxInt;
p.otherParticleCl = [0 0 1]*p.maxInt;% color for all particles other than those in the track
p.stoppedParticleCl = [1 0 1]*p.maxInt;
p.lineWidth = 0.3;
p.RegLineCl = [0 1 0]*p.maxInt; % common line color
p.stoppeLineCl = [1 0 0]*p.maxInt; % when the line stopped

if nargin < 1
    numTracks = 100;
end
p.cmap = hsv(numTracks);
p.cmap = p.cmap(randperm(size(p.cmap,1)),:)*p.maxInt;

end