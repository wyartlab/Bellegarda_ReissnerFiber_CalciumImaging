function rois = getROIcell

    rois = cell(0);
    n = 1;
    while(1)  
        [x,y] = ginput;
        if(length(x) == 0)
            break;
        end
        patch(x,y,'r','FaceAlpha',0.49);
        text(x(1),y(1),num2str(n),'Color','g');
        rois = [rois [ x y ] ];
        n = n+1;
    end
    
%     saveas(gcf,['ROI-capture-' datestr(clock,30)],'png');
%     delete(findobj('FaceAlpha',0.49));
%     delete(findobj('Color','g'));
    
end
    
    
        
