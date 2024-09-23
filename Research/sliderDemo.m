function sliderDemo
    clear; clc;

    f = figure(1);
    
    x = 0:0.1:2*pi;
    y = @(A) A*sin(x);
    A = 1;
    p = plot(x, y(A));

    axis tight
    axis([0 2*pi -10 10])
    set(gca, 'position', [0.1 0.25 0.85 0.7]);

    
    h = uicontrol(...
    'parent' , f,...
    'units' , 'normalized',...
    'style', 'slider',...
    'position', [0.05 0.05 0.9 0.05],...
    'min' , 1,...
    'max' , 10,...
    'value' , A,...
    'callback', @sliderCallback);

    hLstn = handle.listener(h,'ActionEvent',@sliderCallback); %#ok

    function sliderCallback(~,~)
        delete(p);
        p = plot(x, y(get(h,'value')));
        axis tight
        axis([0 2*pi -10 10])
    end

end