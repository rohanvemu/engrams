%% Code
load session3_noisy_a_46.mat
neuron = neuron_network_imaging'; 
fire_seq = zeros(size(neuron, 1), size(neuron, 2));
for i = 1:size(neuron,1)
    [pks,locs] = findpeaks(neuron(i, :)); 
    y = i.*ones(1, length(locs));
    x = (locs./100);
    fire_seq(i, locs(locs<100)) = 1;
    hold on 
    drawnow
    plot(x, y, '.k') 
    xlabel("Time")
    ylabel("Neuron ID")
    title("Rasterization")
    xlim([0 10])
    ylim([0 200])
end 