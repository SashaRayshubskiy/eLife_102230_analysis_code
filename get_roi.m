function [roi_points] = get_roi()

hold on;

order_avg    = [ rgb('Blue'); rgb('Green'); rgb('Red'); rgb('Black'); rgb('Purple'); rgb('Brown'); rgb('Indigo'); rgb('DarkRed') ];
order_single = [ rgb('LightBlue'); rgb('LightGreen'); rgb('LightSalmon'); rgb('LightGray'); rgb('Violet'); rgb('Bisque'); rgb('Plum'); rgb('LightPink') ];
colorindex = 0;

[xv, yv] = (getline(gca, 'closed'));

currcolor_avg    = order_avg(1+mod(colorindex,size(order_single,1)),:);
plot(xv, yv, 'Linewidth', 1,'Color',currcolor_avg);
text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor_avg,'FontSize',12);

roi_points = [xv, yv];
end