netname = 'TA';
staname = 'H65A';
twin = [1 30];
MTCfile = [RFDIR netname '_' staname '_RF.mat'];

RFmat = load(MTCfile, 'time','radRF','qRF','svRF','bin');
%RFmat = load('NM_HENM_filteredRF.mat');

R = RFmat.radRF';
t = RFmat.time;
rayP = RFmat.bin(:, 2);
epiDist = RFmat.bin(:, 1);

% Remove NaN traces or zero traces
validTraces = ~any(isnan(R), 2) & ~all(R == 0, 2);
R = R(validTraces, :);
rayP = rayP(validTraces);
epiDist = epiDist(validTraces);

scale = 0.5;
singlePlot = 1;
pws=0;

% parse velocity model
H  = vmodel(1, :);
Vp = vmodel(2, :);
Vs = vmodel(3, :);
% calculate arrivals
shft = 0;
[t1, t2, t3] = travelTimesAppx(Vp, Vs, H, rayP, 3, 2);
tPhase = [t1, t2, t3] - shft;

%tPhase = [];
h1=figure(2);clf;
RFWigglePlot(R, t, rayP, epiDist, [], twin, scale, singlePlot,pws);
sgtitle(sprintf('%s.%s MTC-RF', netname, staname));
