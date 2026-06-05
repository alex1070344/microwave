clear

wn = 600;
minn = 100;
%% 計算橫軸
x = linspace(-4,4,81);
y = linspace(0,2,21);
w = linspace(600,1200,1+wn);
m = linspace(0,10,1+minn);

%% 計算各個參數的u(x)值
% 溫度(低 中 高)81點
templow = (x <= 0) .* (-x/4);
tempmed = (x>=-3 & x <= 0) .* (x/3+1) + (x > 0 & x <= 3) .* (-x/3+1);
temphigh = (x >= 0) .* (x/4);
templen = length(templow);
% 重量(輕 中 重)21點
wlight = (y <= 1) .* (1-y);
wmed = (y >= 0 & y <= 1) .* y + (y > 1 & y <= 2) .* (2-y);
wheavy = (y >= 1) .* (y-1);
wlen = length(wlight);
% 功率(低 中 高)601點
wattlow = (w <= 800) .* (4-w/200);
wattmed = (w <= 900) .* (w/300-2) + (w>900 & w<1200).*(4-w/300);
watthigh = (w >= 1000) .* (w/200-5);
wattlen = length(wattlow);
% 時間(短 中 長)101點
minshort = (m <= 5) .* (1-m/5);
minmed = (m <= 5) .* (m/5) + (m>5) .* (-m/5+2);
minlong = (m >= 5) .* (m/5-1);
minlen = length(minshort);

% 整合(繪圖)
namecol = {templow,tempmed,temphigh,wlight,wmed,wheavy,wattlow,wattmed,watthigh,minshort,minmed,minlong};
namehorizon = {x,y,w,m};
xlab = {'溫度(°C)','重量(kg)','功率(watt)','時間(min)'}; 
leg = {'low','median','high';'light','median','heavy';'low','median','high';'short','median','long'};
ind_name=1;
figure(1);
for fig = 1:4
    subplot(2,2,fig)
    plot(namehorizon{fig}, namecol{ind_name}, 'b', 'LineWidth', 2);
    hold on;
    plot(namehorizon{fig}, namecol{ind_name+1}, 'g', 'LineWidth',2);
    hold on;
    plot(namehorizon{fig}, namecol{ind_name+2}, 'r', 'LineWidth', 2);
    ind_name=ind_name+3;
    xlabel(sprintf('%s', xlab{fig}));
    ylabel('u(x)');
    title(sprintf('figure : %d',fig))
    legend(leg{fig,1},leg{fig,2},leg{fig,3})
    grid on;
end

%預分配
watt = zeros(length(x),length(y));
wattmat = zeros(length(x),length(y));
wattmod = zeros(length(x),length(y));
wattcen = zeros(length(x),length(y));

mins = zeros(length(x),length(y));
minmat = zeros(length(x),length(y));
minmod = zeros(length(x),length(y));
mincen = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        % alpha cut(取小)
        alphacut = [min(templow(i),wheavy(j)), min(templow(i),wmed(j)), min(templow(i),wlight(j));
            min(tempmed(i),wheavy(j)), min(tempmed(i),wmed(j)), min(tempmed(i),wlight(j));
            min(temphigh(i),wheavy(j)), min(temphigh(i),wmed(j)), min(temphigh(i),wlight(j))];

        % 功率cutted(watt)
        wattresultcut = cell(1,9);
        for cutn = 1:3
            wattresultcut{1,cutn} = min(alphacut(1,cutn),watthigh);
            wattresultcut{1,cutn+3} = min(alphacut(2,cutn),wattmed);
            wattresultcut{1,cutn+6} = min(alphacut(3,cutn),wattlow);
        end
        wattresult = zeros(1, length(w));
        for res = 1:length(wattresultcut)
            wattresult = max(wattresult,wattresultcut{res});
        end
        % 結果(center of gravity)
        sum_watt = sum(wattresult);
        watt(i,j) = sum(wattresult .* w) / sum_watt;

        % mean of maxima
        max_val = max(wattresult);
        max_indices = find(wattresult == max_val);% 取最大值區間
        wattmat(i,j) = sum(w(max_indices))/length(max_indices);
        
        % modified mean of maxima
        wattmod(i,j) = (w(max_indices(1))+w(max_indices(length(max_indices))))/2;
        
        % center avarage
        sums1 = 0; % 分子
        sums2 = 0; % 分母
        for J = 1:length(wattresultcut)
            max_val = max(wattresultcut{J});
            max_indices = find(wattresultcut{J} == max_val);% 取最大值區間
            max_ind = round(median(max_indices));% 找出中點
            sums1 = sums1 + w(max_ind)*wattresultcut{J}(max_ind);
            sums2 = sums2 + wattresultcut{J}(max_ind);
        end
        wattcen(i,j) = sums1/sums2;

        % 時間cutted(times)
        minresultcut = cell(1,9);
        for cutn = 1:3
            minresultcut{1,cutn} = min(alphacut(cutn,1),minlong);
            minresultcut{1,cutn+3} = min(alphacut(cutn,2),minmed);
            minresultcut{1,cutn+6} = min(alphacut(cutn,3),minshort);
        end

        minresult = zeros(size(minresultcut{1}));
        for res = 1:9
            minresult = max(minresult,minresultcut{res});
        end

        % 結果(center of gravity)
        sum_min = sum(minresult);
        mins(i,j) = sum(minresult .* m) / sum_min;

        % mean of maxima
        max_val = max(minresult);
        max_indices = find(minresult == max_val);% 取最大值區間
        minmat(i,j) =  sum(m(max_indices))/length(max_indices);
        
        % modified mean of maxima
        minmod(i,j) = (m(max_indices(1))+m(max_indices(length(max_indices))))/2;
        
        % center avarage
        sums1 = 0; % 分子
        sums2 = 0; % 分母
        for J = 1:length(minresultcut)
            max_val = max(minresultcut{J});
            max_indices = find(minresultcut{J} == max_val);% 取最大值區間
            max_ind = round(median(max_indices));% 找出中點
            sums1 = sums1 + m(max_ind)*minresultcut{J}(max_ind);
            sums2 = sums2 + minresultcut{J}(max_ind);
        end
        mincen(i,j) = sums1/sums2;
    end
end

% 繪圖
resultsname = {'Center of Gravity','Mean of Maxima','Modified Mean of Maxima','Center Avarage'};
%%分開繪圖
% 功率(w)
% wresult = {watt,wattmat,wattmod,wattcen};
% for n_plot = 1:length(resultsname)
%     figure(n_plot+1)
%     surf(y,x,wresult{n_plot});
%     xlabel('重量(kg)');
%     ylabel('溫度(°C)');
%     zlabel('功率(w)');
%     title(sprintf('功率關係圖 : %s'),resultsname{n_plot})
% end
% 
% % 時間(min)
% mresult = {mins,minmat,minmod,mincen};
% for n_plot = 1:length(resultsname)
%     figure(n_plot+5);
%     surf(y,x,mresult{n_plot});
%     xlabel('重量(kg)');
%     ylabel('溫度(°C)');
%     zlabel('時間(min)');
%     title(sprintf('時間關係圖 : %s'),resultsname{n_plot})
% end

%%合併繪圖
% 功率(w)
wresult = {watt,wattmat,wattmod,wattcen};
for n_plot = 1:length(resultsname)
    figure(2)
    subplot(2,2,n_plot)
    surf(y,x,wresult{n_plot});
    xlabel('重量(kg)');
    ylabel('溫度(°C)');
    zlabel('功率(w)');
    title(sprintf('功率關係圖 : %s'),resultsname{n_plot})
end

% 時間(min)
mresult = {mins,minmat,minmod,mincen};
for n_plot = 1:length(resultsname)
    figure(3);
    subplot(2,2,n_plot)
    surf(y,x,mresult{n_plot});
    xlabel('重量(kg)');
    ylabel('溫度(°C)');
    zlabel('時間(min)');
    title(sprintf('時間關係圖 : %s'),resultsname{n_plot})
end