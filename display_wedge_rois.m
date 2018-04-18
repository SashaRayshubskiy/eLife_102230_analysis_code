function display_wedge_rois( rois )

% Show wedge ROIs
colorindex = 0;
ac = get_analysis_constants;
order    = ac.order;

WEDGE_COUNT = 16;
for w = 1:WEDGE_COUNT
    xv = rois{w}{1};
    yv = rois{w}{2};
        
    currcolor    = order(1+mod(colorindex,size(order,1)),:);

    plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
    
    colorindex = colorindex + 1;
end

end

