%A script to scale the statistical wavelet and ricker wavelet to the deterministic wavelet and
%to combine/weight the deterministic wavelets, creating a weighted
%deterministic wavelet that closely matched the scaled statistical wavelet.

%Author: T. Berridge
%Date: 14/07/2020

clc; % Clear Command Window
close all; %Close all figures
clear; %Erase all existing variables


%load wavelets with 4ms smapling interval and 256ms lnegth
%data in form [time(ms):amplitude]
load('wavelet_C1.txt') %deterministic at 20/06a-C1
load('wavelet_4.txt') %deterministic at 20/06-4
load('wavelet_B14Z.txt') %deterministic at 20/06a-B14Z
exstat1=load('wavelet_extendedstatistical1.txt'); %Created across my study area using parameters CNOOC used for whole Buzzard field
ricker=load('wavelet_ricker.txt');
%% - Scale statistical wavelet
stacked_wavelet = (wavelet_C1 + wavelet_4 + wavelet_B14Z)/3; %equally weighted

scale_factor1 = mean(stacked_wavelet(23:43,2)./exstat1(23:43,2)); % scale factor taken from aplitudes over 80ms window across the centre of the wavelet
exstat1_scaled = [exstat1(:,1),exstat1(:,2)*scale_factor1]; %scaled statistical wavelet

%% - Scale ricker wavelet
scale_factor2 = stacked_wavelet(33,2)/ricker(33,2); % scale factor taken from aplitudes over 80ms window across the centre of the wavelet
ricker_scaled = [ricker(:,1),ricker(:,2)*scale_factor2]; %scaled ricker wavelet

%% - calculate well weightings to 1 d.p.
x1=linspace(0.1,0.9,9); %min. well weighting is 0.1 and max. well weighting is 0.9
y1=linspace(0.1,0.9,9);
z1=linspace(0.1,0.9,9);

for i=1:length(x1)
    for j=1:length(y1)
        for k=1:length(z1)
            %all possible wighting combinations
            weighted_wavelet1(:,i,j,k) = [(x1(i)*wavelet_C1(:,2)+y1(j)*wavelet_4(:,2)+z1(k)*wavelet_B14Z(:,2))];
            scaled_stat_wavelet1(:,i,j,k) = exstat1_scaled(:,2);
        end
    end
end

residual1 = weighted_wavelet1 - scaled_stat_wavelet1; %residual at all possible wighting combinations
mean_residual1 = mean(abs(residual1)); %mean residual at all weighting combinations
mean_residual1 = squeeze(mean_residual1); %convert from a 1x9x9x9 matrix to a 9x9x9 matrix

for i=1:length(mean_residual1)
    for j=1:length(mean_residual1)
        for k=1:length(mean_residual1)
            if i+j+k==10 %to ensure weighting sums to 1
                mean_resid1(i,j,k) = mean_residual1(i,j,k); %removes all weightings that don't sum to 1
            end
        end
    end
end
resid1 = nonzeros(mean_resid1); %removes zeros
minimum1=min(resid1); %the lowest residual between weighted stacked deterministic and statistical
[i1,j1,k1]=find(mean_residual1==minimum1); %finds position of lowest residual

fprintf('The optimal weighting of the wells (to 1 d.p in order of input) is: %.2f %.2f %.2f \n', i1/10, j1/10, k1/10)

%% - check the code works
x1=0.2;
y1=0.7;
z1=0.1;

weighted_wavelet1 = [exstat1(:,1),(x1*wavelet_C1(:,2)+y1*wavelet_4(:,2)+z1*wavelet_B14Z(:,2))];
diff1 = mean(abs(weighted_wavelet1(:,2) - exstat1_scaled(:,2)));
%%

figure;
hold on
grid on
%plot(wavelet_4(:,1),wavelet_4(:,2), 'k')
%plot(wavelet_C1(:,1),wavelet_C1(:,2), 'b')
%plot(wavelet_B14Z(:,1),wavelet_B14Z(:,2), 'r')
plot(exstat1_scaled(:,1),exstat1_scaled(:,2), '-ko','Markersize', 4, 'MarkerFaceColor', [0,0,0], 'LineWidth',1.5)
plot(weighted_wavelet1(:,1),weighted_wavelet1(:,2), '-bx', 'Markersize', 6, 'MarkerFaceColor', [0,0,1], 'LineWidth',1.5) %0.2 0.7 0.1 weighting
plot(ricker_scaled(:,1),ricker_scaled(:,2), '-s', 'Color', [0,0.5,0],'Markersize', 5, 'MarkerFaceColor', [0,0.5,0], 'LineWidth',1.5)
xlim([-128 128])
ylim([-9e6 7e6])
xlabel('Time (ms)')
ylabel('Scaled Amplitude')
legend('Scaled statistical wavelet', 'Weighted stacked deterministic wavelet', 'Scaled ricker wavelet', 'FontSize',8)
hold off

figure; %absolute residual in amplitude 
hold on
grid on
plot(stacked_wavelet(:,1),abs(stacked_wavelet(:,2)-exstat1_scaled(:,2)), '-rs','Markersize', 5, 'MarkerFaceColor', [1,0,0], 'LineWidth',1.5) %equal weighting
plot(weighted_wavelet1(:,1),abs(weighted_wavelet1(:,2)-exstat1_scaled(:,2)), '-bx', 'Markersize', 5, 'MarkerFaceColor', [0,0,1], 'LineWidth',1.5) %0.2 0.7 0.1 weighting
xlabel('Time (ms)')
ylabel('Absolute Residual Amplitude')
xlim([-128 128])
legend('Stacked deterministic wavelet', 'Weighted stacked deterministic wavelet', 'FontSize',8)

figure; %Stacked wavelet - equally weighted
hold on
grid on
%rectangle('Position',[-40 -14e6 80 22e6],'FaceColor',[.75 .75 .75],'EdgeColor','w') %starting coordinate, xlength, ylength
plot(wavelet_C1(:,1),wavelet_C1(:,2), 'b')
plot(wavelet_4(:,1),wavelet_4(:,2), 'g')
plot(wavelet_B14Z(:,1),wavelet_B14Z(:,2), 'r')
plot(stacked_wavelet(:,1),stacked_wavelet(:,2), 'k', 'LineWidth',1.5) %equal weighting
xlabel('Time (ms)')
ylabel('Amplitude')
legend('20/06a-C1 deterministic wavelet', '20/06-4 deterministic wavelet', '20/06a-B14Z deterministic wavelet', 'Stacked deterministic wavelet', 'FontSize',8)
xlim([-128 128])
ylim([-14e6 8e6])

