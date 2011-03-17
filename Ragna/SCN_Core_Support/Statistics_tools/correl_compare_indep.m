function [r1, r2, rdiff, rdiff_p] = correl_compare_indep(y_g1, y_g2)
%
% Compute statistics on the difference between correlation coefficients for indepndent samples g1 and g2.
% correlations are computed over y (n observations x k variables) for each
% group.
%
% Tor Wager, Aug. 2008
%
% Source, Hubert Blalock, Social Statistics, NY: McGraw-Hill, 1972:  406-407. 
% From notes on Garson stats website.
%

k = size(y_g1, 2);

if k ~= size(y_g2, 2), error('Num. variables for each group must match.'); end

r1 = corrcoef(y_g1);
r2 = corrcoef(y_g2);

n1 = size(y_g1, 1);
n2 = size(y_g1, 1);

% Fisher's Z scores and stats on each group
[tmp, r1_sig, r1_z] = r2z(r1(:), n1);

r1_sig = reshape(r1_sig, k, k);
r1_z = reshape(r1_z, k, k);

[tmp, r2_sig, r2_z] = r2z(r2(:), n2);

r2_sig = reshape(r2_sig, k, k);
r2_z = reshape(r2_z, k, k);

disp('Group 1'); correlation_to_text( r1, r1_sig);
disp('Group 2'); correlation_to_text( r2, r2_sig);

% difference and stats

rdiff = r1 - r2;

rdiff_z = r1_z - r2_z;
rdiff_se =  sqrt(1./(n1 - 3) + 1./(n2 - 3));

rdiff_zscore = rdiff_z ./ rdiff_se;

rdiff_p = 2 * (1 - normcdf(abs(rdiff_zscore)));

disp('Difference'); 
print_matrix(rdiff);

disp('Difference p-values'); 
print_matrix(rdiff_p);

end


