%% Multi-variate regression 
% Set data

cell_1_data = VmFilt_A2_L;
cell_2_data = VmFilt_A2_R;


BIN_SIZE = 0.050; % s
DT_EPHYS = ephys_SR * BIN_SIZE;
DT_YAW   = ball_SR * BIN_SIZE;

t_down = squeeze(mean(reshape(t_all{1}, [DT_EPHYS, length(t_all{1})/DT_EPHYS]), 1));

cell_1_down = squeeze(mean(reshape(cell_1_data, [ DT_EPHYS, length(cell_2_data)/DT_EPHYS ] ),1));
cell_2_down = squeeze(mean(reshape(cell_2_data, [ DT_EPHYS, length(cell_2_data)/DT_EPHYS ] ),1));

yaw_all_down = squeeze(mean(reshape(yaw_all{1}, [ DT_YAW, length( yaw_all{1} )/DT_YAW ]),1));
fwd_all_down = squeeze(mean(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));
fwd_all_std = squeeze(std(reshape(fwd_all{1}, [ DT_YAW, length( fwd_all{1} )/DT_YAW ]),1));

% Exclude trials when the fly is standing
cell_1_down_temp = cell_1_down;
cell_2_down_temp = cell_2_down;
yaw_all_down_temp = yaw_all_down;
fwd_all_down_temp = fwd_all_down;

cell_1_down = [];
cell_2_down = [];
yaw_all_down = [];
fwd_all_down = [];

FWD_STD_THRESHOLD = 0.01;

i=1;
while( i < length(fwd_all_down_temp) )
    
    cur_fwd_std = fwd_all_std(i);
    
    if( cur_fwd_std > FWD_STD_THRESHOLD )
        cell_1_down(end+1) = cell_1_down_temp(i);
        cell_2_down(end+1) = cell_2_down_temp(i);
        yaw_all_down(end+1) = yaw_all_down_temp(i);
        fwd_all_down(end+1) = fwd_all_down_temp(i);
    end
    
    i = i + 1;
end

%% Plot sanity check variables

figure;

subplot(2,2,1);
scatter(cell_1_down, yaw_all_down);
grid on;

subplot(2,2,2);
scatter(cell_1_down, fwd_all_down);
grid on;

subplot(2,2,3);
scatter(cell_2_down, yaw_all_down);
grid on;

subplot(2,2,4);
scatter(cell_2_down, fwd_all_down);
grid on;

%% Multi-variate regression using mvregress

 % X = n x p
 % Y = n x d
 % beta = p x d

%X = [cell_1_down', cell_2_down'];
%Y = [fwd_all_down', yaw_all_down' ];
X = [ cell_1_down', cell_2_down'];
Y = [ fwd_all_down' ];
[ beta, Sigma, E, CovB, logL ] = mvregress( X, Y );

figure;
plot(E)

% calculate_r_squared_adjusted( beta, Y, X );
% 1 = fwd, 2 = yaw
SStotal = [];
SSresid = [];
R2adjusted = [];

n = size(Y,1);
p = size(X,2);

for j = 1:size(Y,2)
    SSresid(j) = sum(squeeze(E(:,j)).^2);
    %SStotal(j) = (length(Y(:,j))-1) * var(Y(:,j));
    
    SStotal(j) = sum((Y(:,j) - mean(squeeze(Y(:,j)))).^2);
    
    %R2adjusted(j) = 1 - ((SSresid(j) / SStotal(j))*((n-1)/(n-p-1)));
    R2adjusted(j) = 1 - (SSresid(j) / SStotal(j));
end

figure;
bar(R2adjusted);

