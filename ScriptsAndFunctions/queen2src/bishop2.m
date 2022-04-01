function varargout=bishop(fullvelsmooth,n,dim,varargin)
    a = varargin{1};
    switch dim
        case 2
            vel = squeeze(fullvelsmooth(:,:,((n-a.first)/a.increment)+1));                      % extract velocity field #1 from full data set
            vel2 = squeeze(fullvelsmooth(:,:,((n-a.first)/a.increment)+1+1));                   % extract velocity field #2 from full data set
            vel = sortrows(vel,[1 -2]);                                                         % sort rows for processing
            vel2 = sortrows(vel2,[1 -2]);                                                       % sort rows for processing
            xnodes = numel(unique(vel(:,1)));                                                   % determine number of velocity nodes in x-direction
            ynodes = numel(unique(vel(:,2)));                                                   % determine number of velocity nodes in y-direction
            x = reshape(vel(:,1),ynodes,xnodes)*a.lengthcalib_axis;                             % convert to meters
            y = reshape(vel(:,2),ynodes,xnodes)*a.lengthcalib_axis;                             % convert to meters
            u = reshape(vel(:,3),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);            % convert to m/s
            v = reshape(vel(:,4),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);            % convert to m/s
            u2 = reshape(vel2(:,3),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);          % convert to m/s
            v2 = reshape(vel2(:,4),ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);          % convert to m/s
            dx = abs(x(1,2)-x(1,1));                                                            % compute x-dimension grid spacing
            dy = abs(y(2,1)-y(1,1));                                                            % compute y-dimension grid spacing
            Dx = ones(size(x))*dx;                                                              % matrix form of x-dimension grid spacing
            Dy = ones(size(y))*dy;                                                              % matrix form of y-dimension grid spacing
            ua = (u+u2)/2;                                                                      % compute advection velocity
            va = (v+v2)/2;                                                                      % compute advection velocity
            xa = x + ua*a.deltaT;                                                               % compute advected fluid particle new positions    
            ya = y + va*a.deltaT;                                                               % compute advected fluid particle new positions   
            u2a = griddata(x,y,u2,xa,ya);                                                       % compute advected fluid particle velocities at new positions     
            v2a = griddata(x,y,v2,xa,ya);                                                       % compute advected fluid particle velocities at new positions     
            if a.viscous == 0                                                                   % compute dp/dx and dp/dy
                dpdx = -a.rho*((u2a-u)/a.deltaT);                                               % omit viscous term
                dpdy = -a.rho*((v2a-v)/a.deltaT);                                               % omit viscous term
            else
                dpdx = -a.rho*(((u2a-u)/a.deltaT) - a.nu*4*del2(ua,dx,dy));                     % include viscous term
                dpdy = -a.rho*(((v2a-v)/a.deltaT) - a.nu*4*del2(va,dx,dy));                     % include viscous term
            end
            x_min = min(min(x)) + a.nodecrop*dx;                                                % crop boundaries of domain
            x_max = max(max(x)) - a.nodecrop*dx;
            y_min = min(min(y)) + a.nodecrop*dy;
            y_max = max(max(y)) - a.nodecrop*dy;
            dpdx(y>y_max | y<y_min | x>x_max | x<x_min) = 0;
            dpdy(y>y_max | y<y_min | x>x_max | x<x_min) = 0;
            if a.smooth_dp == 1                                                                 % optional nearest neighbor filter
                dpdx = nanmean(cat(3,dpdx,circshift(dpdx,[-1 -1]),circshift(dpdx,[0 -1]),circshift(dpdx,[1 -1]),circshift(dpdx,[1 0]),circshift(dpdx,[1 1]),circshift(dpdx,[0 1]),circshift(dpdx,[-1 1]),circshift(dpdx,[-1 0])),3);
                dpdy = nanmean(cat(3,dpdy,circshift(dpdy,[-1 -1]),circshift(dpdy,[0 -1]),circshift(dpdy,[1 -1]),circshift(dpdy,[1 0]),circshift(dpdy,[1 1]),circshift(dpdy,[0 1]),circshift(dpdy,[-1 1]),circshift(dpdy,[-1 0])),3);
            end
            p_left = cumsum(dpdx.*Dx,2);                                                        % compute pressure in domain via two horizontal and two vertical integration paths
            p_top = -cumsum(dpdy.*Dy,1);
            p_bottom = flipud(cumsum(flipud(dpdy.*Dy),1));
            p_right = -fliplr(cumsum(fliplr(dpdx.*Dx),2));
            dpdx_ul = dpdx;                                                                     % initialize dp/dx matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)
            dpdx_ur = fliplr(dpdx);
            dpdx_ll = flipud(dpdx);
            dpdx_lr = fliplr(flipud(dpdx));
            dpdy_ul = dpdy;                                                                     % initialize dp/dy matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)                 
            dpdy_ur = fliplr(dpdy);                                                  
            dpdy_ll = flipud(dpdy);
            dpdy_lr = fliplr(flipud(dpdy));
            for i = 1:length(dpdx(1,:))
                dpdx_ul1(:,i) = circshift(dpdx_ul(:,i),-(i-1));                                 % shift diagonals vertically for summation
                dpdx_ur1(:,i) = circshift(dpdx_ur(:,i),-(i-1));
                dpdx_ll1(:,i) = circshift(dpdx_ll(:,i),-(i-1));
                dpdx_lr1(:,i) = circshift(dpdx_lr(:,i),-(i-1));
                dpdy_ul1(:,i) = circshift(dpdy_ul(:,i),-(i-1));                                 % shift diagonals vertically for summation
                dpdy_ur1(:,i) = circshift(dpdy_ur(:,i),-(i-1));
                dpdy_ll1(:,i) = circshift(dpdy_ll(:,i),-(i-1));
                dpdy_lr1(:,i) = circshift(dpdy_lr(:,i),-(i-1));
            end
            px_ul1 = cumsum(dpdx_ul1.*Dx,2);                                                    % summation of dp/dx for lower triangle of domain
            px_ur1 = -cumsum(dpdx_ur1.*Dx,2);
            px_ll1 = cumsum(dpdx_ll1.*Dx,2);
            px_lr1 = -cumsum(dpdx_lr1.*Dx,2);
            py_ul1 = -cumsum(dpdy_ul1.*Dy,2);                                                   % summation of dp/dy for lower triangle of domain
            py_ur1 = -cumsum(dpdy_ur1.*Dy,2);
            py_ll1 = cumsum(dpdy_ll1.*Dy,2);
            py_lr1 = cumsum(dpdy_lr1.*Dy,2);
            for i = 1:length(dpdx(1,:))
                px_ul1(:,i) = circshift(px_ul1(:,i),i-1);                                       % shift diagonals back to original position
                px_ur1(:,i) = circshift(px_ur1(:,i),i-1);
                px_ll1(:,i) = circshift(px_ll1(:,i),i-1);
                px_lr1(:,i) = circshift(px_lr1(:,i),i-1);
                py_ul1(:,i) = circshift(py_ul1(:,i),i-1);                                       % shift diagonals back to original position
                py_ur1(:,i) = circshift(py_ur1(:,i),i-1);
                py_ll1(:,i) = circshift(py_ll1(:,i),i-1);
                py_lr1(:,i) = circshift(py_lr1(:,i),i-1);
            end
            for i = 1:length(dpdx(:,1))
                dpdx_ul2(i,:) = circshift(dpdx_ul(i,:),[0 -(i-1)]);                             % shift diagonals horizontally for summation
                dpdx_ur2(i,:) = circshift(dpdx_ur(i,:),[0 -(i-1)]);
                dpdx_ll2(i,:) = circshift(dpdx_ll(i,:),[0 -(i-1)]);
                dpdx_lr2(i,:) = circshift(dpdx_lr(i,:),[0 -(i-1)]);
                dpdy_ul2(i,:) = circshift(dpdy_ul(i,:),[0 -(i-1)]);                             % shift diagonals horizontally for summation
                dpdy_ur2(i,:) = circshift(dpdy_ur(i,:),[0 -(i-1)]);
                dpdy_ll2(i,:) = circshift(dpdy_ll(i,:),[0 -(i-1)]);
                dpdy_lr2(i,:) = circshift(dpdy_lr(i,:),[0 -(i-1)]);
            end
            px_ul2 = cumsum(dpdx_ul2.*Dx,1);                                                    % summation of dp/dx for upper triangle of domain
            px_ur2 = -cumsum(dpdx_ur2.*Dx,1);
            px_ll2 = cumsum(dpdx_ll2.*Dx,1);
            px_lr2 = -cumsum(dpdx_lr2.*Dx,1);
            py_ul2 = -cumsum(dpdy_ul2.*Dy,1);                                                   % summation of dp/dy for upper triangle of domain
            py_ur2 = -cumsum(dpdy_ur2.*Dy,1);
            py_ll2 = cumsum(dpdy_ll2.*Dy,1);
            py_lr2 = cumsum(dpdy_lr2.*Dy,1);
            for i = 1:length(dpdx(:,1))
                px_ul2(i,:) = circshift(px_ul2(i,:),[0 i-1]);                                   % shift diagonals back to original position
                px_ur2(i,:) = circshift(px_ur2(i,:),[0 i-1]);
                px_ll2(i,:) = circshift(px_ll2(i,:),[0 i-1]);
                px_lr2(i,:) = circshift(px_lr2(i,:),[0 i-1]);
                py_ul2(i,:) = circshift(py_ul2(i,:),[0 i-1]);                                   % shift diagonals back to original position
                py_ur2(i,:) = circshift(py_ur2(i,:),[0 i-1]);
                py_ll2(i,:) = circshift(py_ll2(i,:),[0 i-1]);
                py_lr2(i,:) = circshift(py_lr2(i,:),[0 i-1]);
            end
            px_ull = tril(px_ul1,-1);                                                           % extract lower trangle of computed pressure fields
            px_url = tril(px_ur1,-1);
            px_lll = tril(px_ll1,-1);
            px_lrl = tril(px_lr1,-1);
            py_ull = tril(py_ul1,-1);
            py_url = tril(py_ur1,-1);
            py_lll = tril(py_ll1,-1);
            py_lrl = tril(py_lr1,-1);
            px_ulu = triu(px_ul2);                                                              % extract upper triangle of computed pressure fields
            px_uru = triu(px_ur2);
            px_llu = triu(px_ll2);
            px_lru = triu(px_lr2);
            py_ulu = triu(py_ul2);
            py_uru = triu(py_ur2);
            py_llu = triu(py_ll2);
            py_lru = triu(py_lr2);
            px_ul = px_ulu+px_ull;                                                              % combine upper and lower triangles of computed pressure fields
            px_ur = fliplr(px_uru+px_url);
            px_ll = flipud(px_llu+px_lll);
            px_lr = fliplr(flipud(px_lru+px_lrl));
            py_ul = py_ulu+py_ull;
            py_ur = fliplr(py_uru+py_url);
            py_ll = flipud(py_llu+py_lll);
            py_lr = fliplr(flipud(py_lru+py_lrl));
            p_ul = px_ul+circshift(py_ul,[0 -1]);                                               % stitch diagonal staircase integration paths together
            p_ur = px_ur+circshift(py_ur,[0 1]);
            p_ll = px_ll+circshift(py_ll,[0 -1]);
            p_lr = px_lr+circshift(py_lr,[0 1]);
            P = nanmedian(cat(3,p_left,p_top,p_bottom,p_right,p_ul,p_ur,p_ll,p_lr),3);          % select median pressure value
            if a.smooth_p == 1                                                                  % optional nearest neighbor filter
                P = nanmean(cat(3,P,circshift(P,[-1 -1]),circshift(P,[0 -1]),circshift(P,[1 -1]),circshift(P,[1 0]),circshift(P,[1 1]),circshift(P,[0 1]),circshift(P,[-1 1]),circshift(P,[-1 0])),3);
            end
            varargout{1} = x;
            varargout{2} = y;
            varargout{3} = xnodes;
            varargout{4} = ynodes;
            varargout{5} = dpdx;
            varargout{6} = dpdy;
            varargout{7} = u;
            varargout{8} = v;
            varargout{9} = P;
            
        case 3    
            textbox = varargin{end};
            vel = squeeze(fullvelsmooth(:,:,((n-a.first)/a.increment)+1));                      % extract velocity field #1 from full data set
            vel2 = squeeze(fullvelsmooth(:,:,((n-a.first)/a.increment)+1+1));                   % extract velocity field #2 from full data set
            vel = sortrows(vel,[1 -2 3]);                                                       % sort rows for processing
            vel2 = sortrows(vel2,[1 -2 3]);                                                     % sort rows for processing
            xnodes = numel(unique(vel(:,1)));                                                   % determine number of velocity nodes in x-direction
            ynodes = numel(unique(vel(:,2)));                                                   % determine number of velocity nodes in y-direction
            znodes = numel(unique(vel(:,3)));                                                   % determine number of velocity nodes in z-direction
            X = reshape(vel(:,1),znodes,ynodes,xnodes)*a.lengthcalib_axis;                      % convert to meters
            Y = reshape(vel(:,2),znodes,ynodes,xnodes)*a.lengthcalib_axis;                      % convert to meters
            Z = reshape(vel(:,3),znodes,ynodes,xnodes)*a.lengthcalib_axis;                      % convert to meters
            U = reshape(vel(:,4),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);     % convert to m/s
            V = reshape(vel(:,5),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);     % convert to m/s
            W = reshape(vel(:,6),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);     % convert to m/s
            U2 = reshape(vel2(:,4),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);   % convert to m/s
            V2 = reshape(vel2(:,5),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);   % convert to m/s
            W2 = reshape(vel2(:,6),znodes,ynodes,xnodes)*(a.lengthcalib_vel/a.timecalib_vel);   % convert to m/s
            dX = abs(X(1,1,2)-X(1,1,1));                                                        % compute x-dimension grid spacing
            dY = abs(Y(1,2,1)-Y(1,1,1));                                                        % compute y-dimension grid spacing
            dZ = abs(Z(2,1,1)-Z(1,1,1));                                                        % compute z-dimension grid spacing
            DX = ones(size(X))*dX;                                                              % matrix form of x-dimension grid spacing
            DY = ones(size(Y))*dY;                                                              % matrix form of y-dimension grid spacing
            DZ = ones(size(Z))*dZ;                                                              % matrix form of z-dimension grid spacing
            if strcmp(a.gradientplane, 'xy') == 1 | strcmp(a.gradientplane, 'iso') == 1         % check gradientplane selection
                for h = 1:znodes
                    msg = ['Computing pressure in x-y plane for z-node #' num2str(h) ' of ' num2str(znodes)];  
                    set(textbox               , ...
                        'String'              , msg, ...
                        'BackgroundColor'     , [0 .4 0], ...
                        'HorizontalAlignment' , 'center', ...
                        'ForegroundColor'     , 'w'...
                        )
                    pause(.005)
                    x = squeeze(X(h,:,:));                                                      % extract current x-y plane data from full data set
                    y = squeeze(Y(h,:,:));
                    u = squeeze(U(h,:,:));
                    v = squeeze(V(h,:,:));
                    u2 = squeeze(U2(h,:,:));
                    v2 = squeeze(V2(h,:,:));
                    Dx = squeeze(DX(h,:,:));
                    Dy = squeeze(DY(h,:,:));
                    ua = (u+u2)/2;                                                              % compute advection velocity
                    va = (v+v2)/2;                                                              % compute advection velocity
                    xa = x + ua*a.deltaT;                                                       % compute advected fluid particle new positions  
                    ya = y + va*a.deltaT;                                                       % compute advected fluid particle new positions  
                    u2a = griddata(x,y,u2,xa,ya);                                               % compute advected fluid particle velocities at new positions     
                    v2a = griddata(x,y,v2,xa,ya);                                               % compute advected fluid particle velocities at new positions     
                    if a.viscous == 0                                                           % compute dp/dx and dp/dy
                        dpdx = -a.rho*((u2a-u)/a.deltaT);                                       % omit viscous term
                        dpdy = -a.rho*((v2a-v)/a.deltaT);                                       % omit viscous term
                    else
                        dpdx = -a.rho*(((u2a-u)/a.deltaT) - a.nu*4*del2(ua,dX,dY));             % include viscous term
                        dpdy = -a.rho*(((v2a-v)/a.deltaT) - a.nu*4*del2(va,dX,dY));             % include viscous term
                    end
                    x_min = min(min(x)) + a.nodecrop*dX;                                        % crop boundaries of domain
                    x_max = max(max(x)) - a.nodecrop*dX;
                    y_min = min(min(y)) + a.nodecrop*dY;
                    y_max = max(max(y)) - a.nodecrop*dY;
                    dpdx(y>y_max | y<y_min | x>x_max | x<x_min) = 0;
                    dpdy(y>y_max | y<y_min | x>x_max | x<x_min) = 0;
                    if a.smooth_dp == 1                                                         % optional nearest neighbor filter
                        dpdx = nanmean(cat(3,dpdx,circshift(dpdx,[-1 -1]),circshift(dpdx,[0 -1]),circshift(dpdx,[1 -1]),circshift(dpdx,[1 0]),circshift(dpdx,[1 1]),circshift(dpdx,[0 1]),circshift(dpdx,[-1 1]),circshift(dpdx,[-1 0])),3);
                        dpdy = nanmean(cat(3,dpdy,circshift(dpdy,[-1 -1]),circshift(dpdy,[0 -1]),circshift(dpdy,[1 -1]),circshift(dpdy,[1 0]),circshift(dpdy,[1 1]),circshift(dpdy,[0 1]),circshift(dpdy,[-1 1]),circshift(dpdy,[-1 0])),3);
                    end
                    p_left_xy(h,:,:) = cumsum(dpdx.*Dx,2);                                      % compute pressure in x-y plane via two x and two y integration paths
                    p_top_xy(h,:,:) = -cumsum(dpdy.*Dy,1);
                    p_bottom_xy(h,:,:) = flipud(cumsum(flipud(dpdy.*Dy),1));
                    p_right_xy(h,:,:) = -fliplr(cumsum(fliplr(dpdx.*Dx),2));
                    dpdx_ul = dpdx;                                                             % initialize dp/dx matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)
                    dpdx_ur = fliplr(dpdx);
                    dpdx_ll = flipud(dpdx);
                    dpdx_lr = fliplr(flipud(dpdx));
                    dpdy_ul = dpdy;                                       
                    dpdy_ur = fliplr(dpdy);                                                     % initialize dp/dy matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)
                    dpdy_ll = flipud(dpdy);
                    dpdy_lr = fliplr(flipud(dpdy));
                    for i = 1:length(dpdx(1,:))
                        dpdx_ul1(:,i) = circshift(dpdx_ul(:,i),-(i-1));                         % shift diagonals vertically for summation
                        dpdx_ur1(:,i) = circshift(dpdx_ur(:,i),-(i-1));
                        dpdx_ll1(:,i) = circshift(dpdx_ll(:,i),-(i-1));
                        dpdx_lr1(:,i) = circshift(dpdx_lr(:,i),-(i-1));
                        dpdy_ul1(:,i) = circshift(dpdy_ul(:,i),-(i-1));                         % shift diagonals vertically for summation
                        dpdy_ur1(:,i) = circshift(dpdy_ur(:,i),-(i-1));
                        dpdy_ll1(:,i) = circshift(dpdy_ll(:,i),-(i-1));
                        dpdy_lr1(:,i) = circshift(dpdy_lr(:,i),-(i-1));
                    end
                    px_ul1 = cumsum(dpdx_ul1.*Dx,2);                                            % summation of dp/dx for lower triangle of domain
                    px_ur1 = -cumsum(dpdx_ur1.*Dx,2);
                    px_ll1 = cumsum(dpdx_ll1.*Dx,2);
                    px_lr1 = -cumsum(dpdx_lr1.*Dx,2);
                    py_ul1 = -cumsum(dpdy_ul1.*Dy,2);                                           % summation of dp/dy for lower triangle of domain
                    py_ur1 = -cumsum(dpdy_ur1.*Dy,2);
                    py_ll1 = cumsum(dpdy_ll1.*Dy,2);
                    py_lr1 = cumsum(dpdy_lr1.*Dy,2);
                    for i = 1:length(dpdx(1,:))
                        px_ul1(:,i) = circshift(px_ul1(:,i),i-1);                               % shift diagonals back to original position
                        px_ur1(:,i) = circshift(px_ur1(:,i),i-1);
                        px_ll1(:,i) = circshift(px_ll1(:,i),i-1);
                        px_lr1(:,i) = circshift(px_lr1(:,i),i-1);
                        py_ul1(:,i) = circshift(py_ul1(:,i),i-1);                               % shift diagonals back to original position
                        py_ur1(:,i) = circshift(py_ur1(:,i),i-1);
                        py_ll1(:,i) = circshift(py_ll1(:,i),i-1);
                        py_lr1(:,i) = circshift(py_lr1(:,i),i-1);
                    end
                    for i = 1:length(dpdx(:,1))
                        dpdx_ul2(i,:) = circshift(dpdx_ul(i,:),[0 -(i-1)]);                     % shift diagonals horizontally for summation
                        dpdx_ur2(i,:) = circshift(dpdx_ur(i,:),[0 -(i-1)]);
                        dpdx_ll2(i,:) = circshift(dpdx_ll(i,:),[0 -(i-1)]);
                        dpdx_lr2(i,:) = circshift(dpdx_lr(i,:),[0 -(i-1)]);
                        dpdy_ul2(i,:) = circshift(dpdy_ul(i,:),[0 -(i-1)]);                     % shift diagonals horizontally for summation
                        dpdy_ur2(i,:) = circshift(dpdy_ur(i,:),[0 -(i-1)]);
                        dpdy_ll2(i,:) = circshift(dpdy_ll(i,:),[0 -(i-1)]);
                        dpdy_lr2(i,:) = circshift(dpdy_lr(i,:),[0 -(i-1)]);
                    end
                    px_ul2 = cumsum(dpdx_ul2.*Dx,1);                                            % summation of dp/dx for upper triangle of domain
                    px_ur2 = -cumsum(dpdx_ur2.*Dx,1);
                    px_ll2 = cumsum(dpdx_ll2.*Dx,1);
                    px_lr2 = -cumsum(dpdx_lr2.*Dx,1);
                    py_ul2 = -cumsum(dpdy_ul2.*Dy,1);                                           % summation of dp/dy for upper triangle of domain
                    py_ur2 = -cumsum(dpdy_ur2.*Dy,1);
                    py_ll2 = cumsum(dpdy_ll2.*Dy,1);
                    py_lr2 = cumsum(dpdy_lr2.*Dy,1);
                    for i = 1:length(dpdx(:,1))
                        px_ul2(i,:) = circshift(px_ul2(i,:),[0 i-1]);                           % shift diagonals back to original position
                        px_ur2(i,:) = circshift(px_ur2(i,:),[0 i-1]);
                        px_ll2(i,:) = circshift(px_ll2(i,:),[0 i-1]);
                        px_lr2(i,:) = circshift(px_lr2(i,:),[0 i-1]);
                        py_ul2(i,:) = circshift(py_ul2(i,:),[0 i-1]);                           % shift diagonals back to original position
                        py_ur2(i,:) = circshift(py_ur2(i,:),[0 i-1]);
                        py_ll2(i,:) = circshift(py_ll2(i,:),[0 i-1]);
                        py_lr2(i,:) = circshift(py_lr2(i,:),[0 i-1]);
                    end
                    px_ull = tril(px_ul1,-1);                                                   % extract lower trangle of computed pressure fields
                    px_url = tril(px_ur1,-1);
                    px_lll = tril(px_ll1,-1);
                    px_lrl = tril(px_lr1,-1);
                    py_ull = tril(py_ul1,-1);
                    py_url = tril(py_ur1,-1);
                    py_lll = tril(py_ll1,-1);
                    py_lrl = tril(py_lr1,-1);
                    px_ulu = triu(px_ul2);                                                      % extract upper trangle of computed pressure fields
                    px_uru = triu(px_ur2);
                    px_llu = triu(px_ll2);
                    px_lru = triu(px_lr2);
                    py_ulu = triu(py_ul2);
                    py_uru = triu(py_ur2);
                    py_llu = triu(py_ll2);
                    py_lru = triu(py_lr2);
                    px_ul = px_ulu+px_ull;                                                      % combine upper and lower triangles of computed pressure fields
                    px_ur = fliplr(px_uru+px_url);
                    px_ll = flipud(px_llu+px_lll);
                    px_lr = fliplr(flipud(px_lru+px_lrl));
                    py_ul = py_ulu+py_ull;
                    py_ur = fliplr(py_uru+py_url);
                    py_ll = flipud(py_llu+py_lll);
                    py_lr = fliplr(flipud(py_lru+py_lrl));
                    p_ul_xy(h,:,:) = px_ul+circshift(py_ul,[0 -1]);                             % stitch diagonal staircase integration paths together
                    p_ur_xy(h,:,:) = px_ur+circshift(py_ur,[0 1]);
                    p_ll_xy(h,:,:) = px_ll+circshift(py_ll,[0 -1]);
                    p_lr_xy(h,:,:) = px_lr+circshift(py_lr,[0 1]);
                end
            end 
            if strcmp(a.gradientplane, 'yz') == 1 | strcmp(a.gradientplane, 'iso') == 1         % check gradientplane selection
                for j = 1:xnodes
                    msg = ['Computing pressure in y-z plane for x-node #' num2str(j) ' of ' num2str(xnodes)];  
                    set(textbox               , ...
                        'String'              , msg, ...
                        'BackgroundColor'     , [0 .4 0], ...
                        'HorizontalAlignment' , 'center', ...
                        'ForegroundColor'     , 'w'...
                        )
                    pause(.005)
                    y = squeeze(Y(:,:,j));                                                      % extract current y-z plane data from full data set
                    z = squeeze(Z(:,:,j));
                    v = squeeze(V(:,:,j));
                    w = squeeze(W(:,:,j));
                    v2 = squeeze(V2(:,:,j));
                    w2 = squeeze(W2(:,:,j));
                    Dy = squeeze(DY(:,:,j));
                    Dz = squeeze(DZ(:,:,j));
                    va = (v+v2)/2;                                                              % compute advection velocity
                    wa = (w+v2)/2;                                                              % compute advection velocity
                    ya = y + va*a.deltaT;                                                       % compute advected fluid particle new positions
                    za = z + wa*a.deltaT;                                                       % compute advected fluid particle new positions
                    v2a = griddata(y,z,v2,ya,za);                                               % compute advected fluid particle velocities at new positions 
                    w2a = griddata(y,z,w2,ya,za);                                               % compute advected fluid particle velocities at new positions 
                    if a.viscous == 0                                                           % compute dp/dx and dp/dy
                        dpdy = -a.rho*((v2a-v)/a.deltaT);                                       % omit viscous term
                        dpdz = -a.rho*((w2a-w)/a.deltaT);                                       % omit viscous term
                    else
                        dpdy = -a.rho*(((v2a-v)/a.deltaT) - a.nu*4*del2(va,dY,dZ));             % include viscous term
                        dpdz = -a.rho*(((w2a-w)/a.deltaT) - a.nu*4*del2(wa,dY,dZ));             % include viscous term
                    end
                    y_min = min(min(y)) + a.nodecrop*dY;                                        % crop boundaries of domain
                    y_max = max(max(y)) - a.nodecrop*dY;
                    z_min = min(min(z)) + a.nodecrop*dZ;
                    z_max = max(max(z)) - a.nodecrop*dZ;
                    dpdy(z>z_max | z<z_min | y>y_max | y<y_min) = 0;
                    dpdz(z>z_max | z<z_min | y>y_max | y<y_min) = 0;
                    if a.smooth_dp == 1                                                         % optional nearest neighbor filter
                        dpdy = nanmean(cat(3,dpdy,circshift(dpdy,[-1 -1]),circshift(dpdy,[0 -1]),circshift(dpdy,[1 -1]),circshift(dpdy,[1 0]),circshift(dpdy,[1 1]),circshift(dpdy,[0 1]),circshift(dpdy,[-1 1]),circshift(dpdy,[-1 0])),3);
                        dpdz = nanmean(cat(3,dpdz,circshift(dpdz,[-1 -1]),circshift(dpdz,[0 -1]),circshift(dpdz,[1 -1]),circshift(dpdz,[1 0]),circshift(dpdz,[1 1]),circshift(dpdz,[0 1]),circshift(dpdz,[-1 1]),circshift(dpdz,[-1 0])),3);
                    end
                    p_back_yz(:,:,j) = cumsum(dpdz.*Dz,1);                                      % compute pressure in y-z plane via two y and two z integration paths
                    p_top_yz(:,:,j) = -cumsum(dpdy.*Dy,2);
                    p_bottom_yz(:,:,j) = fliplr(cumsum(fliplr(dpdy.*Dy),2));
                    p_front_yz(:,:,j) = -flipud(cumsum(flipud(dpdz.*Dz),1));
                    dpdy_ul = dpdy;                                                             % initialize dp/dy matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)
                    dpdy_ur = fliplr(dpdy);
                    dpdy_ll = flipud(dpdy);
                    dpdy_lr = fliplr(flipud(dpdy));
                    dpdz_ul = dpdz;                                       
                    dpdz_ur = fliplr(dpdz);                                                     % initialize dp/dz matrices to integrate in domain from upper left (ul), upper right (ur), lower left (ll), and lower right (lr)
                    dpdz_ll = flipud(dpdz);
                    dpdz_lr = fliplr(flipud(dpdz));
                    dpdy_ul1 = [];                                                              % clear derived dp/dy quantities due to use in x-y loop
                    dpdy_ur1 = [];  
                    dpdy_ll1 = [];  
                    dpdy_lr1 = [];  
                    for i = 1:length(dpdy(1,:))
                        dpdy_ul1(:,i) = circshift(dpdy_ul(:,i),-(i-1));                         % shift diagonals vertically for summation
                        dpdy_ur1(:,i) = circshift(dpdy_ur(:,i),-(i-1));
                        dpdy_ll1(:,i) = circshift(dpdy_ll(:,i),-(i-1));
                        dpdy_lr1(:,i) = circshift(dpdy_lr(:,i),-(i-1));
                        dpdz_ul1(:,i) = circshift(dpdz_ul(:,i),-(i-1));                         % shift diagonals vertically for summation
                        dpdz_ur1(:,i) = circshift(dpdz_ur(:,i),-(i-1));
                        dpdz_ll1(:,i) = circshift(dpdz_ll(:,i),-(i-1));
                        dpdz_lr1(:,i) = circshift(dpdz_lr(:,i),-(i-1));
                    end
                    py_ul1 = -cumsum(dpdy_ul1.*Dy,2);                                           % summation of dp/dy for upper triangle of domain
                    py_ur1 = cumsum(dpdy_ur1.*Dy,2);
                    py_ll1 = -cumsum(dpdy_ll1.*Dy,2);
                    py_lr1 = cumsum(dpdy_lr1.*Dy,2);
                    pz_ul1 = cumsum(dpdz_ul1.*Dz,2);                                            % summation of dp/dz for upper triangle of domain
                    pz_ur1 = cumsum(dpdz_ur1.*Dz,2);
                    pz_ll1 = -cumsum(dpdz_ll1.*Dz,2);
                    pz_lr1 = -cumsum(dpdz_lr1.*Dz,2);
                    for i = 1:length(dpdy(1,:))
                        py_ul1(:,i) = circshift(py_ul1(:,i),i-1);                               % shift diagonals back to original position
                        py_ur1(:,i) = circshift(py_ur1(:,i),i-1);
                        py_ll1(:,i) = circshift(py_ll1(:,i),i-1);
                        py_lr1(:,i) = circshift(py_lr1(:,i),i-1);
                        pz_ul1(:,i) = circshift(pz_ul1(:,i),i-1);                               % shift diagonals back to original position
                        pz_ur1(:,i) = circshift(pz_ur1(:,i),i-1);
                        pz_ll1(:,i) = circshift(pz_ll1(:,i),i-1);
                        pz_lr1(:,i) = circshift(pz_lr1(:,i),i-1);
                    end
                    dpdy_ul2 = [];                                                              % clear derived dp/dy quantities due to use in x-y loop
                    dpdy_ur2 = [];  
                    dpdy_ll2 = [];  
                    dpdy_lr2 = [];  
                    for i = 1:length(dpdy(:,1))
                        dpdy_ul2(i,:) = circshift(dpdy_ul(i,:),[0 -(i-1)]);                     % shift diagonals horizontally for summation
                        dpdy_ur2(i,:) = circshift(dpdy_ur(i,:),[0 -(i-1)]);
                        dpdy_ll2(i,:) = circshift(dpdy_ll(i,:),[0 -(i-1)]);
                        dpdy_lr2(i,:) = circshift(dpdy_lr(i,:),[0 -(i-1)]);
                        dpdz_ul2(i,:) = circshift(dpdz_ul(i,:),[0 -(i-1)]);                     % shift diagonals horizontally for summation
                        dpdz_ur2(i,:) = circshift(dpdz_ur(i,:),[0 -(i-1)]);
                        dpdz_ll2(i,:) = circshift(dpdz_ll(i,:),[0 -(i-1)]);
                        dpdz_lr2(i,:) = circshift(dpdz_lr(i,:),[0 -(i-1)]);
                    end
                    py_ul2 = -cumsum(dpdy_ul2.*Dy,1);                                           % summation of dp/dy for lower triangle of domain
                    py_ur2 = cumsum(dpdy_ur2.*Dy,1);
                    py_ll2 = -cumsum(dpdy_ll2.*Dy,1);
                    py_lr2 = cumsum(dpdy_lr2.*Dy,1);
                    pz_ul2 = cumsum(dpdz_ul2.*Dz,1);                                            % summation of dp/dz for lower triangle of domain
                    pz_ur2 = cumsum(dpdz_ur2.*Dz,1);
                    pz_ll2 = -cumsum(dpdz_ll2.*Dz,1);
                    pz_lr2 = -cumsum(dpdz_lr2.*Dz,1);
                    for i = 1:length(dpdy(:,1))
                        py_ul2(i,:) = circshift(py_ul2(i,:),[0 i-1]);                           % shift diagonals back to original position
                        py_ur2(i,:) = circshift(py_ur2(i,:),[0 i-1]);
                        py_ll2(i,:) = circshift(py_ll2(i,:),[0 i-1]);
                        py_lr2(i,:) = circshift(py_lr2(i,:),[0 i-1]);
                        pz_ul2(i,:) = circshift(pz_ul2(i,:),[0 i-1]);                           % shift diagonals back to original position
                        pz_ur2(i,:) = circshift(pz_ur2(i,:),[0 i-1]);
                        pz_ll2(i,:) = circshift(pz_ll2(i,:),[0 i-1]);
                        pz_lr2(i,:) = circshift(pz_lr2(i,:),[0 i-1]);
                    end
                    py_ull = tril(py_ul1,-1);                                                   % extract lower trangle of computed pressure fields
                    py_url = tril(py_ur1,-1);
                    py_lll = tril(py_ll1,-1);
                    py_lrl = tril(py_lr1,-1);
                    pz_ull = tril(pz_ul1,-1);
                    pz_url = tril(pz_ur1,-1);
                    pz_lll = tril(pz_ll1,-1);
                    pz_lrl = tril(pz_lr1,-1);
                    py_ulu = triu(py_ul2);                                                      % extract upper trangle of computed pressure fields
                    py_uru = triu(py_ur2);
                    py_llu = triu(py_ll2);
                    py_lru = triu(py_lr2);
                    pz_ulu = triu(pz_ul2);
                    pz_uru = triu(pz_ur2);
                    pz_llu = triu(pz_ll2);
                    pz_lru = triu(pz_lr2);
                    py_ul = py_ulu+py_ull;                                                      % combine upper and lower triangles of computed pressure fields
                    py_ur = fliplr(py_uru+py_url);
                    py_ll = flipud(py_llu+py_lll);
                    py_lr = fliplr(flipud(py_lru+py_lrl));
                    pz_ul = pz_ulu+pz_ull;
                    pz_ur = fliplr(pz_uru+pz_url);
                    pz_ll = flipud(pz_llu+pz_lll);
                    pz_lr = fliplr(flipud(pz_lru+pz_lrl));
                    p_ul_yz(:,:,j) = py_ul+circshift(pz_ul,[0 -1]);                             % stitch diagonal staircase integration paths together
                    p_ur_yz(:,:,j) = py_ur+circshift(pz_ur,[0 1]);
                    p_ll_yz(:,:,j) = py_ll+circshift(pz_ll,[0 -1]);
                    p_lr_yz(:,:,j) = py_lr+circshift(pz_lr,[0 1]);
                end
            end
            if strcmp(a.gradientplane, 'iso') == 1                                              % use all integration paths if velocity gradients are similar in x-y and y-z planes
                P = nanmedian(cat(4,p_left_xy,p_top_xy,p_bottom_xy,p_right_xy,p_front_yz,p_back_yz,p_top_yz,p_bottom_yz,p_ul_xy,p_ur_xy,p_ll_xy,p_lr_xy,p_ul_yz,p_ur_yz,p_ll_yz,p_lr_yz),4);
            elseif strcmp(a.gradientplane, 'xy') == 1                                           % use x-y plane integration paths if velocity gradients are greater in x-y plane than in y-z plane
                P = nanmedian(cat(4,p_left_xy,p_top_xy,p_bottom_xy,p_right_xy,p_ul_xy,p_ur_xy,p_ll_xy,p_lr_xy),4);   
            elseif strcmp(a.gradientplane, 'yz') == 1                                           % use y-z plane integration paths if velocity gradients are greater in y-z plane than in x-y plane
                P = nanmedian(cat(4,p_front_yz,p_back_yz,p_top_yz,p_bottom_yz,p_ul_yz,p_ur_yz,p_ll_yz,p_lr_yz),4);
            end
  
            if a.smooth_p == 1                                                                  % optional nearest neighbor filter
                for h = 1:znodes                                                                % smooth x-y planes
                    P(h,:,:) = nanmean(cat(3,squeeze(P(h,:,:)),circshift(squeeze(P(h,:,:)),[-1 -1]),circshift(squeeze(P(h,:,:)),[0 -1]),circshift(squeeze(P(h,:,:)),[1 -1]),circshift(squeeze(P(h,:,:)),[1 0]),circshift(squeeze(P(h,:,:)),[1 1]),circshift(squeeze(P(h,:,:)),[0 1]),circshift(squeeze(P(h,:,:)),[-1 1]),circshift(squeeze(P(h,:,:)),[-1 0])),3);
                end
                for j = 1:xnodes                                                                % smooth y-z planes
                    P(:,:,j) = nanmean(cat(3,squeeze(P(:,:,j)),circshift(squeeze(P(:,:,j)),[-1 -1]),circshift(squeeze(P(:,:,j)),[0 -1]),circshift(squeeze(P(:,:,j)),[1 -1]),circshift(squeeze(P(:,:,j)),[1 0]),circshift(squeeze(P(:,:,j)),[1 1]),circshift(squeeze(P(:,:,j)),[0 1]),circshift(squeeze(P(:,:,j)),[-1 1]),circshift(squeeze(P(:,:,j)),[-1 0])),3);
                end
            end
            varargout{1} = X;
            varargout{2} = Y;
            varargout{3} = Z;
            varargout{4} = xnodes;
            varargout{5} = ynodes;
            varargout{6} = znodes;
            varargout{7} = U;
            varargout{8} = V;
            varargout{9} = W;           
            varargout{10} = P;
    end
end