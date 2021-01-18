close all
clear all

data = load('profiling.dat');

ncards     = 1;
gflops_vec = 2;
bwidth_vec = 3;
gflops_rdc = 4;
bwidth_rdc = 5;
gflops_stc = 6;
bwidth_stc = 7;
gflops_spm = 8;
bwidth_spm = 9;


vec_gflops = [];
vec_bwidth = [];
rdc_gflops = [];
rdc_bwidth = [];
stc_gflops = [];
stc_bwidth = [];
spm_gflops = [];
spm_bwidth = [];

max_cards = max(data(:,ncards));

for n = 1:max_cards
    idx = find(data(:,ncards) == n);
    avg = sum(data(idx,:)) / length(idx);

    vec_gflops = [vec_gflops avg(gflops_vec)];
    vec_bwidth = [vec_bwidth avg(bwidth_vec)];
    rdc_gflops = [rdc_gflops avg(gflops_rdc)];
    rdc_bwidth = [rdc_bwidth avg(bwidth_rdc)];
    stc_gflops = [stc_gflops avg(gflops_stc)];
    stc_bwidth = [stc_bwidth avg(bwidth_stc)];
    spm_gflops = [spm_gflops avg(gflops_spm)];
    spm_bwidth = [spm_bwidth avg(bwidth_spm)];
end

figure(1);
set(gca, 'FontSize', 18);

plot(1:max_cards, vec_gflops, 'ko-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');
hold on
plot(1:max_cards, rdc_gflops, 'ro-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

plot(1:max_cards, stc_gflops, 'go-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

plot(1:max_cards, spm_gflops, 'bo-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

xlim([0.75 max_cards + 0.25]);
set(gca, 'xtick', 1:max_cards);

xlabel('#devices (3 Tesla C2070 + Intel Core i7)')
ylabel('Effective GFLOPS')

legend('Vector arithmetic', 'Reduction', 'Stencil conv', 'SpMV', ...
		'location', 'northwest');
legend boxoff

print('-depsc', 'gflops.eps');

figure(2)

set(gca, 'FontSize', 18);

plot(1:max_cards, vec_bwidth, 'ko-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');
hold on
plot(1:max_cards, rdc_bwidth, 'ro-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

plot(1:max_cards, stc_bwidth, 'go-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

plot(1:max_cards, spm_bwidth, 'bo-', ...
		'linewidth', 2, 'markersize', 6, 'markerfacecolor', 'w');

xlim([0.75 max_cards + 0.25]);
set(gca, 'xtick', 1:max_cards);

xlabel('#devices (3 Tesla C2070 + Intel Core i7)')
ylabel('Effective bandwidth (GB/sec)')

legend('Vector arithmetic', 'Reduction', 'Stencil conv', 'SpMV', ...
		'location', 'northwest');
legend boxoff

print('-depsc', 'bwidth.eps');

