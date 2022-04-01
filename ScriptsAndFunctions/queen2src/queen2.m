classdef queen2 < handle  
  properties (Dependent)
      hAxes
      hLine
  end
  properties (Hidden, Access = protected)
      MainAxes
      version    
      statusBox
      hCol
      handles
      fontsize12
      args
      x
      y
      z
      u
      v
      w
      P
      rb_val
      dpdx
      dpdy
      xnodes
      ynodes
      znodes
  end
  properties (Constant, Hidden, Access = protected)
      ptrTypes = {'hand', 'arrow'}
      wsize = [800 600];
  end
  methods
    % Class Constructor
    function obj = queen2(varargin)
        
        os = getenv('OS');
        if ~isempty(regexp(os,'indow'))
            obj.fontsize12 = 10;
        else
            obj.fontsize12 = 12;
        end
        
        obj.version = '2.0';

        if nargin > 0               
            % Error checking
            errorcheck(obj);
        end

        if mod(nargin-2,2)
            error('Additional arguments must be Param/Value pairs');
        end

        % Create main GUI
        createGUI(obj);

        % Load default parameters
        args = varargin;  
        if isempty(args)
            obj.loadParams('default');
        end
        
    end
    
    % Handle set(app.hAxes) calls
    function out = get.hAxes(obj)
        out = obj.MainAxes;
    end
  end
  methods (Hidden, Access = protected)
%--------------------------------------------------------------------------
% Function for creating the GUI
      function createGUI(obj) 
          % Create main window
          hFig = figure(...
            'Name'                            , 'queen 2.0', ...
            'Numbertitle'                     , 'off', ...
            'Units'                           , 'pixels', ...
            'DockControls'                    , 'on', ...
            'Toolbar'                         , 'none', ...
            'Menubar'                         , 'none', ...
            'Color'                           , [0 0 0], ...
            'DoubleBuffer'                    , 'on', ...
            'Renderer'                        , 'OpenGL', ...
            'ResizeFcn'                       , @obj.figResizeFcn, ...
            'KeyPressFcn'                     , @obj.keypressFcn, ...
            'DeleteFcn'                       , @obj.figDeleteFcn, ...
            'Visible'                         , 'on', ...
            'HandleVisibility'                , 'callback', ...
            'WindowButtonMotionFcn'           , @obj.motionFcn, ...
            'Tag'                             , 'MainFigure', ...
            'defaultUIControlBackgroundColor' , 'w', ...
            'defaultUIControlFontName'        , 'Arial', ...
            'defaultUIControlFontSize'        , obj.fontsize12, ...
            'defaultUIPanelUnits'             , 'pixels', ...
            'defaultUIPanelPosition'          , [0 0 obj.wsize(1) 50], ...
            'defaultAxesUnits'                , 'pixels', ...              
            'defaultAxesXColor'               , 'w', ...
            'defaultAxesYColor'               , 'w', ...
            'defaultAxesZColor'               , 'w', ...
            'defaultAxesFontName'             , 'Arial', ...
            'defaultAxesFontSize'             , obj.fontsize12, ...
            'defaultAxesPosition'             , [385 180 440 330], ...
            'Position'                        , [0 0 obj.wsize] ...
          );
          cameratoolbar(hFig,'show','NoReset')
      
          % Create main axes
          obj.MainAxes = axes('Parent', hFig);
 
          obj.initializeMainAxes();
          
          % Create Processing Parameters menu
          hMenu0 = uimenu('Parent', hFig,  'Label','Processing Parameters');
          % Create View Sub-Menu
          uimenu('Parent', hMenu0, 'Label', 'View...',...
                 'Callback', @obj.viewParams);          
          % Create Edit Sub-Menu
          uimenu('Parent', hMenu0, 'Label', 'Edit...',...
                 'Callback', @obj.editParams);

          % Create Plot Options menu
          hMenu1 = uimenu('Parent', hFig,  'Label','Plot Options');             
          uimenu('Parent', hMenu1, 'Label', 'Clear Axes...',...
                 'Callback', @obj.resetAxes);
          uimenu('Parent', hMenu1, 'Label', 'Edit Title/Labels...',...
                 'Callback', @obj.labelCallback);
          hThemeMenu = uimenu( ...
                              'Parent'    , hMenu1,...
                              'Tag'       , 'ThemeSel', ...
                              'Label'     , 'Change Theme...' ...
                              );
          hTheme(1) = uimenu('Parent',hThemeMenu,'Label','Black');
          hTheme(2) = uimenu('Parent',hThemeMenu,'Label','Blue');
          hTheme(3) = uimenu('Parent',hThemeMenu,'Label','White');
          set(hTheme, {'Checked'} , {'on'; 'off'; 'off'}, ...
                      'Callback', @obj.changeThemeFcn);
          uimenu('Parent', hFig, 'Label', 'Help'             , ... 
                 'Callback', @obj.helpFcn, 'Separator', 'on');
          uimenu('Parent', hFig, 'Label', 'About'            ,...
                 'Callback', @obj.aboutFcn);
          % UI Context Menu for plot customization
          hContext = uicontextmenu('Parent', hFig, 'Tag', 'ContextMenu');
          hCPlotOptions = uimenu('Parent', hContext, 'Label', 'Plot Options');
          hContextMenuItems(1) = uimenu('Parent', hCPlotOptions, 'Label', 'Clear axes');

          set(hContextMenuItems, 'Callback', @obj.contextmenuFcn);

          % Folder selection text boxes and corresponding labels
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'HorizontalAlignment' , 'center', ...
              'FontWeight'          , 'Bold', ...
              'BackgroundColor'     , [.8 .8 .8], ...
              'String'              , 'Velocity Data Folder', ...
              'Position'            , [20 475 250 18] ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'velDataFolderTxt', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'Position'            , [20 450 250 24] ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'HorizontalAlignment' , 'center', ...
              'FontWeight'          , 'Bold', ...
              'BackgroundColor'     , [.8 .8 .8], ...
              'String'              , 'Interface Data Folder', ...
              'Position'            , [20 375 250 18] ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'ifaceDataFolderTxt', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'Position'            , [20 350 250 24] ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'HorizontalAlignment' , 'center', ...
              'FontWeight'          , 'Bold', ...
              'BackgroundColor'     , [.8 .8 .8], ...
              'String'              , 'Processing Parameters File', ...
              'Position'            , [20 305 250 18] ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'pParamFileTxt', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'Position'            , [20 280 250 24] ...
          );
          % Folder selection buttons
          uicontrol(...
              'Parent'            , hFig, ...
              'Style'             , 'pushbutton', ...
              'Tag'               , 'velDataFolder', ...
              'Position'          , [285 460 60 24], ...
              'String'            , 'Browse', ...
              'ToolTipString'     , ''  , ...
              'Callback'          , @obj.setDataFolder ...
          );
          uicontrol(...
              'Parent'            , hFig, ...
              'Style'             ,'pushbutton', ...
              'Tag'               ,'ifaceDataFolder' , ...
              'Position'          , [285 350 60 24], ...
              'String'            , 'Browse', ...
              'ToolTipString'     , '', ...
              'Callback'          , @obj.setDataFolder ...
          );
          uicontrol(...
              'Parent'            , hFig, ...
              'Style'             , 'pushbutton', ...
              'Tag'               , 'pParamFile', ...
              'Position'          , [285 280 60 24], ...
              'String'            , 'Browse', ...
              'ToolTipString'     , '', ...
              'Callback'          , @obj.setDataFolder ...
          );
          % Radio buttons
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'txt2d', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'Enable'              , 'inactive', ...
              'ForegroundColor'     , 'w', ...
              'BackgroundColor'     , 'k', ...
              'Position'            , [110 420 20 20], ...
              'String'              , '2D' ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'Style'               , 'radiobutton', ...
              'Tag'                 , '2d',...
              'Position'            , [90 420 20 20], ...
              'Value'               ,  1, ...
              'ToolTipString'       , '', ...
              'ForegroundColor'     , 'k', ...
              'BackgroundColor'     , 'k', ...
              'Callback'            , @obj.setDim ...
          );
          obj.rb_val = '2d';
          uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'txt3d', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'Enable'              , 'inactive', ...
              'ForegroundColor'     , 'w', ...
              'BackgroundColor'     , 'k', ...
              'Position'            , [170 420 20 20], ...
              'String'              , '3D' ...
          );
          uicontrol(...
              'Parent'              , hFig, ...
              'Style'               , 'radiobutton', ...
              'Tag'                 , '3d',...
              'Position'            , [150 420 20 20], ...
              'Value'               , 0, ...
              'ToolTipString'       , '', ...
              'ForegroundColor'     , 'k', ...
              'BackgroundColor'     , 'k', ...
              'Callback'            , @obj.setDim ...
              );
          % Stop button
          uicontrol(...
              'Parent'            , hFig, ...
              'Style'             , 'pushbutton', ...
              'Tag'               , 'Stop!',...
              'Position'          , [120 140 72 72], ...
              'String'            , 'Stop', ...
              'FontSize'          , 20, ...
              'ToolTipString'     , '', ...
              'BackgroundColor'   , 'r', ...
              'ForegroundColor'   , 'w', ...
              'Value'             , 0, ...
              'Callback'          , @obj.stop...
          );
          % Control panel
          hPanel = uipanel(...
              'Parent'            , hFig              , ...
              'BackgroundColor'   , [.5 .5 .5]        , ...
              'Tag'               , 'ControlPanelAxes', ...
              'BorderType'        , 'etchedin'...
          );
          % Create action buttons
          uicontrol(...
              'Parent'            , hPanel, ...
              'Style'             ,'pushbutton', ...
              'Tag'               ,'plot_vel', ...
              'Position'          ,[420 5 180 36], ...
              'String'            ,'Plot Velocity Fields', ...
              'ToolTipString'     ,'', ...
              'Callback'          , @obj.performAction ...
          );
          uicontrol(...
              'Parent'            , hPanel, ...
              'Style'             , 'pushbutton', ...
              'Tag'               , 'comp_press', ...
              'Position'          , [610 5 180 36], ...
              'String'            , 'Compute Pressure Fields' , ...
              'ToolTipString'     , '', ...
              'Callback'          , @obj.performAction ...
          );
          % Status message box
          obj.statusBox = uicontrol(...
              'Parent'              , hFig, ...
              'style'               , 'text', ...
              'Tag'                 , 'statusBox', ...
              'HorizontalAlignment' , 'left', ...
              'FontWeight'          , 'Bold', ...
              'FontSize'            , obj.fontsize12, ...
              'String'              , 'Status: ', ...
              'Position'            , [285 560 500 20] ...
          );        
          obj.handles = guihandles(hFig);
          % Make axes invisible from outside
          set(findobj(hFig, 'type', 'axes'),...
              'HandleVisibility', 'callback');         
          movegui(hFig, 'center');
          set(hFig, 'Visible' , 'on');
      end %createGUI
%--------------------------------------------------------------------------
% Function for loading the processing parameters
      function stop(obj,varargin)
            cmdWindow = com.mathworks.mde.cmdwin.CmdWin.getInstance();
            cmdWindow.grabFocus();

            %2) Wait for focus transfer to complete (up to 2 seconds)
            focustransferTimer = tic;
            while ~cmdWindow.isFocusOwner
                pause(0.1);  %Pause some small interval
                if (toc(focustransferTimer) > 2)
                    error('Error transferring focus for CTRL+C press.')
                end
            end

            %3) Use Java robot to execute a CTRL+C in the (now focused) command window.

            %3.1)  Setup a timer to relase CTRL + C in 1 second
            %  Try to reuse an existing timer if possible (this would be a holdover
            %  from a previous execution)
            t_all = timerfindall;
            releaseTimer = [];
            ix_timer = 1;
            while isempty(releaseTimer) && (ix_timer<= length(t_all))
                if isequal(t_all(ix_timer).TimerFcn, @releaseCtrl_C)
                    releaseTimer = t_all(ix_timer);
                end
                ix_timer = ix_timer+1;
            end
            if isempty(releaseTimer)
                releaseTimer = timer;
                releaseTimer.TimerFcn = @releaseCtrl_C;
            end
            releaseTimer.StartDelay = 1;
            start(releaseTimer);

            %3.2)  Press CTRL+C
            pressCtrl_C
      end
%--------------------------------------------------------------------------
% Function for loading the processing parameters
      function loadParams(obj,varargin)
        if strcmp(varargin{1},'default')
            filename = 'parameters2';
            obj.args = eval(filename);        
            o = obj.handles.pParamFileTxt;
            set(o,'String',[filename '.m']);
            set(o,'ToolTip',[pwd filesep filename '.m']);
        elseif strcmp(varargin{1},'user')
            o = obj.handles.pParamFileTxt;
            s = get(o,'String');
            filename = s(1:end-2);
            obj.args = eval(filename);
            clear args;
        end
        o = obj.handles.velDataFolderTxt;
        set(o,'String',obj.args.datafolder);
        set(o,'ToolTip',['Data files: ' obj.args.datafolder filesep obj.args.inroot '*.dat']);
        o = obj.handles.ifaceDataFolderTxt;
        if obj.args.blanking
          set(o,'String',obj.args.blankingfolder);
          set(o,'ToolTip',['Interface data files: ' obj.args.blankingfolder, ...
                           filesep obj.args.blankingroot '*.dat']);
        else
          set(o,'String','Select a file to enable.');
          set(o,'ForegroundColor','r');
        end        
      end
%--------------------------------------------------------------------------
% Function for intializing the plot axes
      function initializeMainAxes(obj,varargin)
          ax = obj.MainAxes;
          hcol= colorbar('Peer'       , ax, ...
                         'Units'      , 'Pixels', ...
                         'Color'      , 'w', ...
                         'YColor'     , 'w', ...
                         'Location'   , 'SouthOutside', ...
                         'FontSize'   , obj.fontsize12 ...
                         );
          cpos=get(hcol,'Position');
          cpos(1)=cpos(1)+cpos(3)/8;
          cpos(2)=100;
          cpos(3)=3*cpos(3)/4; % Halve the thickness
          cpos(4)=cpos(4)/2; % Halve the thickness
          set(hcol,'position',cpos)
          xlabel(hcol,'Pa','Color','w')
          colormap(hcol,customcolormap_preset('orange-white-purple'))
          obj.hCol = hcol;
          set(ax,'Color','none')
          xlabel (ax,'x (m)','Color','w')
          ylabel (ax,'y (m)','Color','w')
          zlabel (ax,'z (m)','Color','w')
          title(ax,'','Color','w')
          axis(ax,'square')
          box(ax,'off')
          hold(ax,'on')
          caxis(ax,[0 1]);
      end %initializeMainAxes
%--------------------------------------------------------------------------
% Called when the context menu on a graphics object is selected
      function contextmenuFcn(obj, hObj, varargin) %#ok<INUSL>
          val = get(hObj, 'Label');
          switch get(get(hObj, 'Parent'), 'Label')
              case 'Line Style'
                  set(gco, 'LineStyle', val);
              case 'Line Width'
                  set(gco, 'LineWidth', str2double(val));
              case 'Line Color'
                  set(gco, 'Color', str2num(val)); %#ok<ST2NM>
              case 'Marker'
                  set(gco, 'Marker', val);
              case 'Marker Size'
                  set(gco, 'MarkerSize', str2double(val));
          end
      end %contextmenuFcn
%--------------------------------------------------------------------------
% Checks for file information to be entered
      function ret = checkFileInfo(obj)
          ret = 1;
          box = findobj('Tag','velDataFolderTxt');
          vel = get(box,'ToolTip');
          box = findobj('Tag','ifaceDataFolderTxt');
          iface = get(box,'ToolTip');
          % Ensure that user has entered the necessary input data
          if ( isempty(vel) )
            % Flash error message if input data are missing
            obj.setMessage(obj.statusBox,'error','Please enter input data!',4)
            obj.setMessage(obj.statusBox,'default','',0)
            ret = 0;
          end
      end
%--------------------------------------------------------------------------
% Plots raw input data (velocity field) for verifying if it is in the 
% correct format
      function plotVelocity(obj)
        dim = 2;
        o = findobj('Tag','3d');
        if (get(o(1),'Val')==1) 
          dim = 3;
        end
        obj.resetAxes();
        obj.setMessage(obj.statusBox,'status',...
                       'Verifying input file...',...
                       .25);        
        caxis(obj.MainAxes,[0 1]);
        set(obj.hCol,'Visible','off');
        
        switch dim
            case 2
              for n = obj.args.first:obj.args.increment:obj.args.last
                    fnamev=[obj.args.datafolder filesep obj.args.inroot,...
                    num2str(n, obj.args.numformat) obj.args.fileextension];
                    vel = dlmread(fnamev,obj.args.separator,obj.args.numheaderlines,0);    
                if size(vel,2)~=4
                    errordlg('Wrong number of columns in input file.')
                    obj.setMessage(obj.statusBox,'default','',0)                
                    return
                end
                    a = obj.args;
                    obj.setMessage(obj.statusBox,'status',...
                    ['Plotting velocity field number ' num2str(n) ' of ' num2str(obj.args.last)],...
                    1);                               
                    vel = sortrows(vel,[1 -2]);                                                      
                    xnodes = numel(unique(vel(:,1)));                                                
                    ynodes = numel(unique(vel(:,2)));                                                
                    x = reshape(vel(:,1),ynodes,xnodes)*a.lengthcalib_axis;
                    y = reshape(vel(:,2),ynodes,xnodes)*a.lengthcalib_axis;
                    u = reshape(vel(:,3),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);
                    v = reshape(vel(:,4),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);
                    axes(obj.MainAxes);
                    cla(obj.hAxes)
                    quiver(x,y,u,v,'r')
                    axis(obj.hAxes,'tight')
                    pause(obj.args.plot_delay)
                    
              end              
            case 3
              for n = obj.args.first:obj.args.increment:obj.args.last
                    fnamev=[obj.args.datafolder filesep obj.args.inroot,...
                    num2str(n, obj.args.numformat) obj.args.fileextension];
                    vel = dlmread(fnamev,obj.args.separator,obj.args.numheaderlines,0);
                if size(vel,2)~=6
                    errordlg('Wrong number of columns in input file.')
                    obj.setMessage(obj.statusBox,'default','',0)                
                    return
                end                
                    a = obj.args;
                    obj.setMessage(obj.statusBox,'status',...
                    ['Plotting velocity field number ' num2str(n) ' of ' num2str(obj.args.last)],...
                    1);   
                    xnodes = numel(unique(vel(:,1)));
                    ynodes = numel(unique(vel(:,2)));
                    znodes = numel(unique(vel(:,3)));
                    X = reshape(vel(:,1),znodes,ynodes,xnodes)*a.lengthcalib_axis;
                    Y = reshape(vel(:,2),znodes,ynodes,xnodes)*a.lengthcalib_axis;
                    Z = reshape(vel(:,3),znodes,ynodes,xnodes)*a.lengthcalib_axis;
                    U = reshape(vel(:,4),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);
                    V = reshape(vel(:,5),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);
                    W = reshape(vel(:,6),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);
                    switch obj.args.gradientplane
                        case 'yz'
                            cla(obj.hAxes)
                            quiver(squeeze(Y(:,:,round(xnodes/2))),squeeze(Z(:,:,round(xnodes/2))),squeeze(V(:,:,round(xnodes/2))),squeeze(W(:,:,round(xnodes/2))),'r')
                            xlabel('y(m)')
                            ylabel('z(m)')
                        case 'xy'
                            cla(obj.hAxes)
                            quiver(squeeze(X(round(znodes/2),:,:)),squeeze(Y(round(znodes/2),:,:)),squeeze(U(round(znodes/2),:,:)),squeeze(V(round(znodes/2),:,:)),'r')
                            xlabel('x(m)')
                            ylabel('y(m)')
                        case 'iso'
                            cla(obj.hAxes)
                            quiver(squeeze(X(round(znodes/2),:,:)),squeeze(Y(round(znodes/2),:,:)),squeeze(U(round(znodes/2),:,:)),squeeze(V(round(znodes/2),:,:)),'r')
                            xlabel('x(m)')
                            ylabel('y(m)')
                    end
                    view(obj.hAxes,2)
                    axis(obj.hAxes,'tight')
                    pause(obj.args.plot_delay)
             end
        end
        obj.setMessage(obj.statusBox,'status','DONE!',4)
        obj.setMessage(obj.statusBox,'default','',0)        
      end
%--------------------------------------------------------------------------
% Called when any of the action buttons is clicked
      function performAction(obj, varargin)
        % Check for file input info       
        if ~obj.checkFileInfo();
            return
        end
        obj.loadParams('user');
        % Determine which action is to be performed
        buttonTag = get(varargin{1},'Tag');
        switch buttonTag
          case 'plot_vel'
            obj.plotVelocity()
            return
          case 'comp_press'
            obj.args.plotfigure = 1;
        end
        % 2d/3d case according to selected radio button
        dim = 2;
        o = findobj('Tag','3d');
        if (get(o(1),'Val')==1) 
          dim = 3;
        end
        % Call to process data using .p file
        try
            obj.queenFcn(dim);
        catch err
            errordlg(err.message);
            % Flash error message if input data are missing
                obj.setMessage(obj.statusBox,'error','Error executing script!',4)
            try
                obj.setMessage(obj.statusBox,'default','',0)
            catch
                disp 'Program ended...'
            end
        end
      end %performAction
%--------------------------------------------------------------------------
% Called when the radio buttons (2d/3d) are pressed
      function setDim(obj, varargin)
          
        % Determine which action is to be performed
        buttonTag = get(varargin{1},'Tag');
        
        if ~strcmp(obj.rb_val,buttonTag)
            obj.rb_val = buttonTag;
        else
            set(varargin{1},'Value',1)            
            return;
        end
        
        switch buttonTag
          case '2d'
            o = findobj('Tag','3d');
            set(o,'Value',0);
          case '3d'
            o = findobj('Tag','2d');
            set(o,'Value',0);
        end
        obj.resetAxes();
      end
%--------------------------------------------------------------------------
% Called when the browse buttons are pressed
      function setDataFolder(obj, varargin)
        buttonTag = get(varargin{1},'Tag');
        boxTag = [buttonTag 'Txt'];
        box = findobj('Tag',boxTag);

        switch buttonTag
          case 'velDataFolder'
            dirpath = uigetdir;
            % Check if user hits cancel
            if dirpath==0
              return;
            end
            dirname =  regexp(dirpath,filesep,'split');
            files = dir(fullfile(dirpath,'*.dat'));
            if size(files,1)==0
                obj.setMessage(obj.statusBox,'error','No data (.dat) files found!',4)
                obj.setMessage(obj.statusBox,'default','',0)
                set(box,'String','','ToolTip','');
                return;
            end
            
            filename = files(1,1).name;
            inroot = regexp(filename,'\.','split');
            inroot = inroot{1,1}(1,:);
            strfind(inroot,obj.args.inroot);
%             keyboard
            if (strfind(inroot,obj.args.inroot)==1)
              set(box,'String',dirname(end),'ToolTip',dirpath);
            else
              obj.setMessage(obj.statusBox,'error',...
                         'inroot parameter and file names don''t match!',...
                         4);
              obj.setMessage(obj.statusBox,'default','',0)
              set(box,'String','','ToolTip','');
            end
            
          case 'ifaceDataFolder'
            dirpath = uigetdir;
            % Check if user hits cancel
            if dirpath==0
              return;
            end
           
            obj.args.blanking = 1;
            
            o = findobj('Tag','pParamFileTxt');
            col = get(o,'ForegroundColor');
            dirname =  regexp(dirpath,filesep,'split');
%             keyboard
            if length(box)>1
              box = box(1);
              col = col{1};
            end
            set(...
                box,...
                'String'           ,dirname(end),...
                'ToolTip'          ,dirpath,...
                'ForegroundColor'  ,col...
                );
          case 'pParamFile'
            [filename, pathname] = uigetfile( ...
              {'*.m;*.fig;*.mat;*.slx;*.mdl',...
               'MATLAB Files (*.m,*.fig,*.mat,*.slx,*.mdl)';
               '*.m',  'Code files (*.m)'; ...
               '*.fig','Figures (*.fig)'; ...
               '*.mat','MAT-files (*.mat)'; ...
               '*.mdl;*.slx','Models (*.slx, *.mdl)'; ...
               '*.*',  'All Files (*.*)'}, ...
               'Pick a file');
            % Check if user hits cancel
            if filename==0
              return;
            end
            set(box,'String',filename,...
                'ToolTip',[pathname filename]);
            obj.loadParams('user');            
        end
      end %setDataFolder
%--------------------------------------------------------------------------
% Called when "Edit title/labels..." menu is selected
      function labelCallback(obj, varargin)

          answer = inputdlg({'Title:', 'X-Label:', 'Y-Label:', 'Z-Label:'}, ...
              'Enter labels', 1, {get(get(obj.MainAxes, 'Title'), 'String'), ...
              get(get(obj.MainAxes, 'XLabel'), 'String'), ...
              get(get(obj.MainAxes, 'YLabel'), 'String'), ...;
              get(get(obj.MainAxes, 'ZLabel'), 'String')});

          if ~isempty(answer)
              title(obj.MainAxes , answer{1});
              xlabel(obj.MainAxes, answer{2});
              ylabel(obj.MainAxes, answer{3});
              zlabel(obj.MainAxes, answer{3});
          end

      end %labelCallback
%--------------------------------------------------------------------------
% Called when "View Processing Parameters" menu is selected
      function viewParams(obj,varargin)
          obj.loadParams('user');
          c = {paramsCells(obj.args)};
          
          h = figure('Units','pixels',...
                     'Position',[500 500 300 600],...
                     'Menubar','none',...
                     'Name','Processing Parameters',...
                     'Numbertitle','off',...
                     'Resize','off');
          
          l = uicontrol('Parent'      , h, ...
                        'Style'       , 'List', ...
                        'Unit'        , 'Pixels', ...
                        'Position'    , [0 0 280 580],...
                        'Min'         , 0,...
                        'Max'         , 2,...
                        'FontSize'    , obj.fontsize12 ...
                      );
          set(l,'String',c{1});
          movegui(h, 'center');
      end
%--------------------------------------------------------------------------
% Called when "Edit Processing Parameters" menu is selected
      function editParams(obj,varargin)
          o = findobj('Tag','pParamFileTxt');
          filename = get(o(1),'ToolTip');
          edit(filename);          
      end
%--------------------------------------------------------------------------
% Called when "Help..." menu is selected
      function helpFcn(obj, varargin)  %#ok<INUSD>

          helpdlg({'Please consult the README2 file for additional information.'},...
                   'Help for queen');

      end %helpFcn
%--------------------------------------------------------------------------
% Called when change theme is selected
      function changeThemeFcn(obj, varargin)
          o = get(varargin{1},'Parent');
          c = get(o,'Children');
          set(c,{'Checked'},{'off';'off';'off'});
          set(varargin{1},'Checked','on');
          hFig = get(obj.MainAxes,'Parent');
          
          switch get(varargin{1},'Label')
            case 'Blue'
              set(hFig,'Color','b')
              set(obj.MainAxes,'XColor','w')
              set(obj.MainAxes,'YColor','w')
              set(obj.MainAxes,'ZColor','w')
              o = findobj('Tag','velDataFolderTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','ifaceDataFolderTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','pParamFileTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              s =  get(get(obj.MainAxes, 'XLabel'), 'String');
              xlabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'YLabel'), 'String');
              ylabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'ZLabel'), 'String');
              zlabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'Title'), 'String');
              title (obj.MainAxes,s,'Color','w')
              o = findobj('Tag','txt2d');
              set(o,'BackgroundColor','b', 'ForegroundColor','w');
              o = findobj('Tag','txt3d');
              set(o,'BackgroundColor','b', 'ForegroundColor','w');
              o = findobj('Tag','statusBox');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','3d');
              set(o,'BackgroundColor','b', 'ForegroundColor','b');
              o = findobj('Tag','2d');
              set(o,'BackgroundColor','b', 'ForegroundColor','b');
              xlabel(obj.hCol,'Pa','Color','w');
            case 'Black'
              set(hFig,'Color','k')
              set(obj.MainAxes,'XColor','w')
              set(obj.MainAxes,'YColor','w')
              set(obj.MainAxes,'ZColor','w')
              o = findobj('Tag','velDataFolderTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','ifaceDataFolderTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','pParamFileTxt');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              s =  get(get(obj.MainAxes, 'XLabel'), 'String');
              xlabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'YLabel'), 'String');
              ylabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'ZLabel'), 'String');
              zlabel (obj.MainAxes,s,'Color','w')
              s =  get(get(obj.MainAxes, 'Title'), 'String');
              title (obj.MainAxes,s,'Color','w')              
              o = findobj('Tag','txt2d');
              set(o,'BackgroundColor','k', 'ForegroundColor','w');
              o = findobj('Tag','txt3d');
              set(o,'BackgroundColor','k', 'ForegroundColor','w');
              o = findobj('Tag','statusBox');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
               o = findobj('Tag','3d');
              set(o,'BackgroundColor','k', 'ForegroundColor','k');
              o = findobj('Tag','2d');
              set(o,'BackgroundColor','k', 'ForegroundColor','k');
              xlabel(obj.hCol,'Pa','Color','w');
            case 'White'
              set(hFig,'Color','w')
              set(obj.MainAxes,'XColor','k')
              set(obj.MainAxes,'YColor','k')
              set(obj.MainAxes,'ZColor','k')
              o = findobj('Tag','velDataFolderTxt');
              set(o,'BackgroundColor',[.2 .2 .2], 'ForegroundColor','w');
              o = findobj('Tag','ifaceDataFolderTxt');
              set(o,'BackgroundColor',[.2 .2 .2], 'ForegroundColor','w');
              o = findobj('Tag','pParamFileTxt');
              set(o,'BackgroundColor',[.2 .2 .2], 'ForegroundColor','w');
              s =  get(get(obj.MainAxes, 'XLabel'), 'String');
              xlabel (obj.MainAxes,s,'Color','k')
              s =  get(get(obj.MainAxes, 'YLabel'), 'String');
              ylabel (obj.MainAxes,s,'Color','k')
              s =  get(get(obj.MainAxes, 'ZLabel'), 'String');
              zlabel (obj.MainAxes,s,'Color','k')
              s =  get(get(obj.MainAxes, 'Title'), 'String');
              title (obj.MainAxes,s,'Color','k')
              o = findobj('Tag','txt2d');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','txt3d');
              set(o,'BackgroundColor','w', 'ForegroundColor','k');
              o = findobj('Tag','statusBox');
              set(o,'BackgroundColor',[.2 .2 .2], 'ForegroundColor','w');
              o = findobj('Tag','3d');
              set(o,'BackgroundColor','w', 'ForegroundColor','w');
              o = findobj('Tag','2d');
              set(o,'BackgroundColor','w', 'ForegroundColor','w');              
             xlabel(obj.hCol,'Pa','Color','k');
          end
          o = findobj('Tag','ifaceDataFolderTxt');
          s = get(o,'String');
          if strcmp(s,'Select a file to enable.')
            set(o,'ForegroundColor','r');
          end
      end
%--------------------------------------------------------------------------
% Called when "About..." menu is selected
      function aboutFcn(obj, varargin)

          helpdlg({
              'This program takes as input one or more text files with ' ... 
              '2D or 3D velocity field data and, optionally, files with coordinates '...
              'of a solid object in the flow, and computes the corresponding '...
              'pressure field. '...
              ' ', ...
              'Copyright John O. Dabiri 2013'}, ...
              sprintf('About queen %s', obj.version));

      end %aboutFcn      
%--------------------------------------------------------------------------
% Called when a key is pressed
      function keypressFcn(obj, varargin)
          k = get(obj.handles.MainFigure, 'CurrentKey');
          switch k
              case 'rightarrow'     % go to the next frame
              case 'leftarrow'      % go to the previous frame
              case 'uparrow'
              case 'downarrow'
          end
      end %keypressFcn
%--------------------------------------------------------------------------
% Called whenever the figure is resized to reposition the components
% "nicely"
      function figResizeFcn(obj, varargin)

      end %figResizeFcn
%--------------------------------------------------------------------------
% Called when the cursor is moved (useful to change the pointer icon)
      function motionFcn(obj, varargin)

      end %motionFcn
%--------------------------------------------------------------------------
% Called when figure is closed
      function figDeleteFcn(obj, varargin)
        delete(obj);
        cleanup;
      end %figDeleteFcn
%--------------------------------------------------------------------------
% Error checking
      function errorcheck(obj)

      end %errorcheck
%--------------------------------------------------------------------------
% Reset figure axes      
      function resetAxes(obj,varargin)
        % 2d/3d case according to selected radio button
        dim = 2;
        
        if isempty(varargin)
            % From this app instance
            o = findobj('Tag','3d');
            if (get(o(1),'Val')==1) 
              dim = 3;
            end
        else
            % From clear menu
            o = findobj('Tag','3d');
            if (get(o(1),'Val')==1) 
              dim = 3;
            end
        end
        cla(obj.hAxes)
        switch dim
            case 2
                view(obj.hAxes,0,90); 
                set(obj.hAxes,'XLim',[0 1],'YLim',[0 1],'ZLim',[-1 1])
                set(obj.hCol,'Visible','on');
                xlabel(obj.hCol,'Pa');
                set(obj.hAxes,'Position',[385 180 440 330])
                caxis(obj.hAxes,[0 1]);
                set(obj.hAxes,'Color','None')
            case 3
                view(obj.hAxes,3); 
                set(obj.hAxes,'XLim',[0 1],'YLim',[0 1],'ZLim',[0 1])
                set(obj.hCol,'Visible','off');
                set(obj.hAxes,'Position',[405 200 360 270])
                caxis(obj.hAxes,[0 1]);
                set(obj.hAxes,'Color','None')                  
        end
        title(obj.hAxes,'');
        pause(0.25)
      end     
%--------------------------------------------------------------------------
% Status messages
      function setMessage(obj,varargin)
        textBox = varargin{1};
        type = varargin{2};
        msg = varargin{3};
        delay = varargin{4};
        switch type
          case 'default'
            set(textBox               , ...
                'String'              , 'Status: ', ...
                'HorizontalAlignment' , 'left' ...
                );
          t = findobj('Tag','ThemeSel');
          if length(t)>1
             t = t(1);
          end
          c = get(t,'Children');
          try
            o = get(c,'Checked');
          catch
            o = get(c{1},'Checked');
          end
          w = strcmp(o,'on');

          obj.changeThemeFcn(c(w));
          case 'status'
            set(textBox               , ...
                'String'              , msg, ...
                'BackgroundColor'     , [0 .4 0], ...
                'HorizontalAlignment' , 'center', ...
                'ForegroundColor'     , 'w'...
                )
          case 'error'
            set(textBox               , ...
                'String'              , msg, ...
                'BackgroundColor'     , [1 0 0], ...
                'HorizontalAlignment' , 'center', ...
                'ForegroundColor'     , 'w'...
                )
        end
        pause(delay)
      end
%--------------------------------------------------------------------------
% Main processing function     
      function queenFcn(obj,varargin)
        dim  = varargin{1};
        for n = obj.args.first:obj.args.increment:obj.args.last
          obj.resetAxes();
          hcol = obj.hCol;
          obj.setMessage(obj.statusBox,'status',...
                             ['Loading velocity field number ' num2str(n) ' of ' num2str(obj.args.last)],...
                             .25);
          fnamev=[obj.args.datafolder filesep obj.args.inroot,...
                num2str(n, obj.args.numformat) obj.args.fileextension];
          fullvel(:,:,((n-obj.args.first)/obj.args.increment)+1) = dlmread(fnamev,obj.args.separator,obj.args.numheaderlines,0);
              
          switch dim
              case 2
                if size(fullvel,2)~=4
                   errordlg('Wrong number of columns in input file.')
                   obj.setMessage(obj.statusBox,'default','',0)                
                   return
                end
              case 3
                if size(fullvel,2)~=6
                   errordlg('Wrong number of columns in input file.')
                   obj.setMessage(obj.statusBox,'default','',0)                
                   return
                end
          end
          % Optional domain blanking
          if obj.args.blanking == 1    
             switch dim
                 case 2 
                    fnameb=[obj.args.blankingfolder filesep obj.args.blankingroot...                              
                             num2str(n, obj.args.numformat) obj.args.fileextension];
                    blank = dlmread(fnameb,obj.args.separator,obj.args.numheaderlines,0);
                    blankindex = inpolygon(squeeze(fullvel(:,1,((n-obj.args.first)/obj.args.increment)+1)),squeeze(fullvel(:,2,((n-obj.args.first)/obj.args.increment)+1)),blank(:,1),blank(:,2));         
                    fullvel(blankindex,3,((n-obj.args.first)/obj.args.increment)+1) = NaN;
                    fullvel(blankindex,4,((n-obj.args.first)/obj.args.increment)+1) = NaN;
                 case 3
                    fnameb=[obj.args.blankingfolder filesep obj.args.blankingroot...                              
                             num2str(n, obj.args.numformat) obj.args.fileextension];
                    blank = dlmread(fnameb,obj.args.separator,obj.args.numheaderlines,0);
                    blankindex = dsearchn(cat(2,squeeze(fullvel(:,1,((n-obj.args.first)/obj.args.increment)+1)),squeeze(fullvel(:,2,((n-obj.args.first)/obj.args.increment)+1)),squeeze(fullvel(:,3,((n-obj.args.first)/obj.args.increment)+1))),blank);         
                    fullvel(blankindex,4,((n-obj.args.first)/obj.args.increment)+1) = NaN;
                    fullvel(blankindex,5,((n-obj.args.first)/obj.args.increment)+1) = NaN; 
                    fullvel(blankindex,6,((n-obj.args.first)/obj.args.increment)+1) = NaN;
             end
          end
        end
        % Optional temporal smoothing filter
        fullvelsmooth = fullvel;
        if obj.args.smooth_t == 1
            switch dim
                case 2
                    for point = 1:length(fullvel(:,1,1))
                        obj.setMessage(obj.statusBox,'status',...
                               ['Temporal filter ' num2str((100*(point/length(fullvel(:,1,1)))),'%1.1f') '% complete'],...
                               0);
                        if sum(~isnan(squeeze(fullvel(point,3,:)))) > 1 && sum(~isnan(squeeze(fullvel(point,4,:)))) > 1
                            fullvelsmooth(point,3,find(~isnan(squeeze(fullvel(point,3,:))))) = feval(fit(find(~isnan(squeeze(fullvel(point,3,:)))),squeeze(fullvel(point,3,~isnan(squeeze(fullvel(point,3,:))))),'smoothingspline','SmoothingParam',0.05),find(~isnan(squeeze(fullvel(point,3,:)))));
                            fullvelsmooth(point,4,find(~isnan(squeeze(fullvel(point,4,:))))) = feval(fit(find(~isnan(squeeze(fullvel(point,4,:)))),squeeze(fullvel(point,4,~isnan(squeeze(fullvel(point,4,:))))),'smoothingspline','SmoothingParam',0.05),find(~isnan(squeeze(fullvel(point,4,:)))));
                        end
                    end
                case 3
                    for point = 1:length(fullvel(:,1,1))
                        obj.setMessage(obj.statusBox,'status',...
                               ['Temporal filter ' num2str((100*(point/length(fullvel(:,1,1)))),'%1.1f') '% complete'],...
                               0);
                        if sum(~isnan(squeeze(fullvel(point,4,:)))) > 1 && sum(~isnan(squeeze(fullvel(point,5,:)))) > 1 && sum(~isnan(squeeze(fullvel(point,6,:)))) > 1
                            fullvelsmooth(point,4,find(~isnan(squeeze(fullvel(point,4,:))))) = feval(fit(find(~isnan(squeeze(fullvel(point,4,:)))),squeeze(fullvel(point,4,~isnan(squeeze(fullvel(point,4,:))))),'smoothingspline','SmoothingParam',0.05),find(~isnan(squeeze(fullvel(point,4,:)))));
                            fullvelsmooth(point,5,find(~isnan(squeeze(fullvel(point,5,:))))) = feval(fit(find(~isnan(squeeze(fullvel(point,5,:)))),squeeze(fullvel(point,5,~isnan(squeeze(fullvel(point,5,:))))),'smoothingspline','SmoothingParam',0.05),find(~isnan(squeeze(fullvel(point,5,:)))));
                            fullvelsmooth(point,6,find(~isnan(squeeze(fullvel(point,6,:))))) = feval(fit(find(~isnan(squeeze(fullvel(point,6,:)))),squeeze(fullvel(point,6,~isnan(squeeze(fullvel(point,6,:))))),'smoothingspline','SmoothingParam',0.05),find(~isnan(squeeze(fullvel(point,6,:)))));
                        end
                    end
            end
        end
        
        for n = obj.args.first:obj.args.increment:obj.args.last-obj.args.increment
          obj.resetAxes();
          hcol = obj.hCol;
          switch dim
            case 2
              obj.setMessage(obj.statusBox,'status',...
                         ['Processing velocity field number ' num2str(n) ' of ' num2str(obj.args.last-obj.args.increment)],...
                         .25);                 
              [obj.x,obj.y,obj.xnodes,obj.ynodes,...
               obj.dpdx,obj.dpdy,obj.u,obj.v,obj.P] = bishop2(fullvelsmooth,n,dim,obj.args);

              % Plot figure
              if obj.args.plotfigure == 1
                set(obj.hAxes,'Color','w')                  
                contourf(obj.hAxes,obj.x,obj.y,obj.P,75,'LineColor','none');
                axis(obj.hAxes,'tight')
                title(['PRESSURE FIELD FOR VELOCITY FIELD #' num2str(n)])
                caxislimit = max(max(abs(obj.P)));
                caxis(obj.hAxes,[-caxislimit caxislimit]);
                colormap(customcolormap_preset('orange-white-purple'))
                pause(obj.args.plot_delay)
              end           
              % Write data to file
              if obj.args.export_data == 1
                switch obj.args.export_format
                    case 'ascii'
                        obj.setMessage(obj.statusBox,'status',...
                                       'Exporting data to ASCII file...',...
                                       1);
                        press = [reshape(obj.x,obj.xnodes*obj.ynodes,1) reshape(obj.y,obj.xnodes*obj.ynodes,1) ...
                        reshape(obj.u,obj.xnodes*obj.ynodes,1) reshape(obj.v,obj.xnodes*obj.ynodes,1) ...
                        reshape(obj.dpdx,obj.xnodes*obj.ynodes,1) reshape(obj.dpdy,obj.xnodes*obj.ynodes,1) ...
                        reshape(obj.P,obj.xnodes*obj.ynodes,1) reshape(abs(obj.P),obj.xnodes*obj.ynodes,1)]; 

                        dlmwrite([obj.args.datafolder filesep obj.args.outroot num2str(n,obj.args.numformat), ...
                                  obj.args.fileextension],press)
                    case 'tecplot'
                        obj.setMessage(obj.statusBox,'status',...
                               'Exporting data to Tecplot file...',...
                               1);
                        press = [reshape(obj.x,obj.xnodes*obj.ynodes,1) reshape(obj.y,obj.xnodes*obj.ynodes,1) ...
                                 reshape(obj.u,obj.xnodes*obj.ynodes,1) reshape(obj.v,obj.xnodes*obj.ynodes,1) ...
                                 reshape(obj.dpdx,obj.xnodes*obj.ynodes,1) reshape(obj.dpdy,obj.xnodes*obj.ynodes,1) ...
                                 reshape(obj.P,obj.xnodes*obj.ynodes,1) reshape(abs(obj.P),obj.xnodes*obj.ynodes,1)];
                        press(isnan(press)) = 0;    
                        if n == obj.args.first
                           pid = fopen([obj.args.datafolder filesep obj.args.outroot 'combined' obj.args.fileextension],'w');                            
                           fprintf(pid,'TITLE = pressure field\n');
                           fprintf(pid,['VARIABLES = x(m), y(m), u(m/s), v(m/s), ' ...
                                        'dp/dx(Pa/m), dp/dy(Pa/m), p(Pa), |p|(Pa)\n']);
                        else
                           pid = fopen([obj.args.datafolder filesep obj.args.outroot 'combined' obj.args.fileextension],'a');                         
                        end
                        fprintf(pid,['ZONE T= pressure_' num2str(n) ' I=' num2str(obj.xnodes)...
                                     ', J=' num2str(obj.ynodes) ', F=POINT\n']);
                        fprintf(pid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n', press');
                        fclose(pid);
                end
              end
            case 3
                  obj.setMessage(obj.statusBox,'status',...
                                 ['Processing velocity field number ' num2str(n,'%04d')],...
                                 .25);              
                  [obj.x,obj.y,obj.z,obj.xnodes,obj.ynodes,obj.znodes,obj.u,obj.v,obj.w,obj.P] = bishop2(fullvelsmooth,n,dim,obj.args,obj.statusBox);
              % Plot figure
              if obj.args.plotfigure == 1
                  p = patch(isosurface(obj.x,obj.y,obj.z,obj.P,mean(mean(mean(obj.P)))));
                  title(['MEAN PRESSURE ISOSURFACE FOR VELOCITY FIELD #' num2str(n)],'FontWeight','bold')
                  set(p,'FaceColor','green','EdgeColor','none');
                  axis(obj.hAxes,'tight')
                  daspect(obj.hAxes,[1,1,1])
                  camlight 
                  lighting gouraud
              end
              % Write data to file
              if obj.args.export_data == 1
                switch obj.args.export_format
                    case 'ascii'
                        obj.setMessage(obj.statusBox,'status',...
                                       'Exporting data to ASCII file...',...
                                       1);
                        press = [squeeze(reshape(obj.x,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.y,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.z,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.u,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.v,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.w,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.P,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(abs(obj.P),obj.xnodes*obj.ynodes*obj.znodes,1,1))];
                        dlmwrite([obj.args.datafolder filesep obj.args.outroot num2str(n,obj.args.numformat), ...
                                  obj.args.fileextension],press)
                    case 'tecplot'
                        obj.setMessage(obj.statusBox,'status',...
                               'Exporting data to Tecplot file...',...
                               1);
                        press = [squeeze(reshape(obj.x,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.y,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.z,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.u,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.v,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.w,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(obj.P,obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 squeeze(reshape(abs(obj.P),obj.xnodes*obj.ynodes*obj.znodes,1,1))...
                                 ];
                        press(isnan(press)) = 0; 
                        if n == obj.args.first
                           pid = fopen([obj.args.datafolder filesep obj.args.outroot 'combined' obj.args.fileextension],'w');
                           fprintf(pid,'TITLE = pressure field\n');
                           fprintf(pid,['VARIABLES =x(m), y(m), z(m), u(m/s), v(m/s), w(m/s), ' ...
                                        'p(Pa), |p|(Pa)\n']);
                        else
                           pid = fopen([obj.args.datafolder filesep obj.args.outroot 'combined' obj.args.fileextension],'a');                     
                        end
                        fprintf(pid,['ZONE T= pressure_' num2str(n) ' I=' num2str(obj.znodes)...
                                     ', J=' num2str(obj.ynodes)  ', K=' num2str(obj.xnodes) ', F=POINT\n']);
                        fprintf(pid,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f \n', press');   
                        fclose(pid);
                end
              end  
          end
        end
        obj.setMessage(obj.statusBox,'status','DONE!',4)
        obj.setMessage(obj.statusBox,'default','',0)        
      end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions
%--------------------------------------------------------------------------
% Housekeeping
function cleanup
  close all; clear all; clc
end
function c = paramsCells(s)
  f = fields(s);
  for i = 1:size(f,1)
    if ischar(getfield(s,f{i}))
      val = getfield(s,f{i});
    else
      val = num2str(getfield(s,f{i}));
    end   
    c{i,1} = [f{i} ' : ' val];
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pressCtrl_C
    import java.awt.Robot;
    import java.awt.event.*;
    SimKey=Robot;
    SimKey.keyPress(KeyEvent.VK_CONTROL);
    SimKey.keyPress(KeyEvent.VK_C);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function releaseCtrl_C(ignore1, ignore2)
    import java.awt.Robot;
    import java.awt.event.*;
    SimKey=Robot;
    SimKey.keyRelease(KeyEvent.VK_CONTROL);
    SimKey.keyRelease(KeyEvent.VK_C);
end