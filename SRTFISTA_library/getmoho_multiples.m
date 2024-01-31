function mask = getmoho_multiples(best_m_out, Vp, Vs, rayP)
  % Input:
  % best_m_out: the radon transformed model matrix
  % Vp, Vs, rayP: parameters needed for travelTimesAppx function
  % Output:
  % mask: filtered matrix to isolate the second Moho multiple
  
  % Constants
  time_range = 3:7; % given time range to find moho phase
  
  % Find strongest moho amplitude in the positive q domain between 3 to 7 seconds
  pos_q_domain = best_m_out(1:size(best_m_out,1)/2, time_range);
  [~, max_idx] = max(pos_q_domain(:));
  [max_row, max_col] = ind2sub(size(pos_q_domain), max_idx);
  t = time_range(max_col);
  
  % Calculate the moho depth H
  H = t / (1/Vs - 1/Vp);
  
  % Calculate theoretical arrival times for moho and multiples
  [~, ~, tP1p2s] = travelTimesAppx(Vp, Vs, H, rayP, 1, 2);
  
  % Finding the corresponding indices in best_m_out
  tP1p2s_idx = round(tP1p2s - min(time_range) + 1);
  
  % Create mask to isolate second multiple
  mask = zeros(size(best_m_out));
  mask(size(best_m_out,1)/2 + 1:end, tP1p2s_idx) = 1;
  
  % Keep only negative values
  mask = mask & (best_m_out < 0);
  
  % Apply mask to best_m_out to isolate the second multiple
  best_m_out_filtered = best_m_out .* mask;
end

