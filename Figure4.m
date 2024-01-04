% load the data
%%% Run a0_setupparameters first
netname = 'BK';
staname = 'JRSC';

dir_path = [localBaseDir 'Prj7_RadonT/Prj7_US_Earthscope/fista_results/' netname '_' staname '/'];
MTCfile = [RFDIR netname '_' staname '_RF.mat'];
RFmat = load(MTCfile, 'time','radRF','qRF','svRF','bin');

tx = RFmat.radRF';
t = RFmat.time;
rayP = RFmat.bin(:, 2)';
dist= RFmat.bin(:, 1);

dt = t(2) - t(1);
% cut data
tBegin = find(t > 0, 1);
tEnd   = find(t > 25, 1);
t = t(tBegin:tEnd);
tx = tx(:, tBegin:tEnd);

dq = (qmax-qmin)/(nq-1);  %value increament of qs
qs = qmin+dq*(0:1:nq-1);   % generate all q values from qmin to qmax

% load the .mat file
load([dir_path 'fista_on_' netname '_' staname  'nocrust.mat']);

% plot the resulting Radon images for every combination of the parameters
for i=1:10%size(results,1)
    fig = figure();
    m_out = cell2mat(results(i,3));
    lambda = cell2mat(results(i,1));
    maxiter = cell2mat(results(i,2));

    RadonPlot(t, qs, m_out, twin, qwin, [], []);
    ylabel('$q (s/km)^2$', 'interpreter', 'latex', 'Fontsize', 20 );

    %save the figure in the specified directory
    saveas(fig,[dir_path 'fista_results_noCrust' num2str(i)],'epsc');
    saveas(fig,[dir_path 'fista_resultsnoCrust' num2str(i)],'png');

    close(fig);
end

% Ask the user for the index of the best model image
best_index = input('Please enter the index of the best image: ');

% Save the m_out and lambda value of the best image
best_m_out = cell2mat(results(best_index,3));

% Display the lambda value of the best_m_out on the screen
best_lambda = cell2mat(results(best_index,1));
disp(['Lambda value of the best image: ' num2str(best_lambda)]);
%%
[Kmask, qps_q] = KdiagMask(best_m_out', t, qs, H, Vp, Vs, 4,'diagonal');%filter model

Kmask = Kmask';
mk = best_m_out.*Kmask;

[dhat] = forward_radon_freq_umd(mk,dt,rayP,qs,N,fmin,fmax);%transform model back to data

rawRF = tx';
filteredRF = dhat';
%%
% Determine the lengths of t and filteredRF
length_t = length(t);
length_filteredRF = size(filteredRF, 2);

% Check if length of t is greater than length of filteredRF, if so pad filteredRF with zeros
if length_t > length_filteredRF
    filteredRF = [filteredRF, zeros(size(filteredRF, 1), length_t - length_filteredRF)];
end
%%
filteredRF_Struct = struct('filteredRF', filteredRF, 'taus', t, 'rayP', rayP);

rawRF_stack = mean(rawRF, 1);
rawRF_stack = rawRF_stack / max(abs(rawRF_stack));
filteredRF_stack= mean(filteredRF, 1);
filteredRF_stack = filteredRF_stack / max(abs(filteredRF_stack));

% % Create original time vector
original_taus = linspace(min(taus), max(taus), length(filteredRF_stack)); 

% % Interpolate to new time vector
filteredRF_stack_resampled = interp1(original_taus, filteredRF_stack, taus, 'linear', 'extrap');
filteredRF_stack_resampled = filteredRF_stack_resampled/ max(abs(filteredRF_stack_resampled));

%filteredRF_stack = pws_stack(rawRF, filteredRF);
disp('done')
%%
figure(2);clf;
subplot(5,3,1:3)
jbfill(t, max(rawRF_stack, 0), zeros(1, length(t)), [0 0 1],'k', 1, 1.0);
jbfill(t, min(rawRF_stack, 0), zeros(1, length(t)), [1 0 0],'k', 1, 1.0);
xlim(twin);
ylim([-1 1]);
hold on;

subplot(5,3,[4:3:10,6:3:12])
RadonPlot(t, qs, best_m_out, [0 25], qwin, [], []);  
%xlim(twin);  % Set the x-axis limits to match the other subplots
hold on;
plot(t, qps_q, 'r')  
plot(t, qps_q + 40, 'r--')
plot(t, qps_q + 300, 'r--')

subplot(5,3,13:15)
jbfill(t, max(filteredRF_stack, 0), zeros(1, length(t)), [0 0 1],'k', 1, 1.0);
jbfill(t, min(filteredRF_stack, 0), zeros(1, length(t)), [1 0 0],'k', 1, 1.0);
xlim(twin);grid on
ylim([-1 1])
set(gca, 'xtick', xticks, 'xticklabel', new_labels);
hold on;
