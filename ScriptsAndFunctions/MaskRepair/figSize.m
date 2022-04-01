function figSize(h,w,hfig)
% Set default figure size to 1/x wide by 1/y tall.
% Returns initial setting.

scrsz = get(0,'ScreenSize');
% figurePosBackup = get(0, 'DefaultFigurePosition');
if nargin < 3
%     set(0,'DefaultFigurePosition',[1, scrsz(4), scrsz(3)/w, scrsz(4)/h]);
    hfig = gcf;
end

set(hfig,'Position',[1, 1, scrsz(3)/w, scrsz(4)/h])
    