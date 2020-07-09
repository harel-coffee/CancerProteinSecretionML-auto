function [] = grouped_scatter(plot_data,group_data,dims,cmap,marker_size)
%GROUPED_SCATTER    Scatter plot, grouping points by color and lines.
%
% GROUPED_SCATTER(PLOT_DATA,GROUP_DATA,DIMS,CMAP) plots points given a
% matrix PLOT_DATA of coordinates. The points are grouped according
% to GROUP_DATA by color and are connected by lines drawn to the center of
% the points in the group.
%
%--------------------------------- INPUTS ---------------------------------
%
% plot_data   NxDIMS matrix of point coordinates.
%
% group_data  Nx1 vector or cell array indicating which group each point
%             belongs to.
%
% dims        (2 or 3) Number of dimensions to plot in. 
%             [DEFAULT = 2]
%
% cmap        Colormap to use for group coloring, supplied as a string or
%             an Mx3 matrix of rgb values, where M is the number of groups.
%             [DEFAULT = 'lines']
%
% marker_size   Scalar or vector of marker sizes. If a vector of marker 
%               sizes is provided, they will be scaled to span a range 
%               specified in ADDITIONAL SETTINGS below. 
%               [DEFAULT = 30]
%

%******* ADDITIONAL SETTINGS *******
draw_lines = false;
label_points = false;
label_align = 'left';
marker_symbol = 'o';
marker_size_range = [10,300];
marker_face_alpha = 0.6;  % 0 (transparent) to 1 (opaque)
marker_edge_alpha = 1;  % 0 (transparent) to 1 (opaque)
marker_outline = 'none';  % specify 'none', 'group', or a color (e.g, 'r', or [1 0 0])
marker_fill = true;
%***********************************

% handle input arguments
if nargin < 5
    marker_size = 30;
end
if nargin < 4 || isempty(cmap)
    cmap = 'lines';
end
if nargin < 3 || isempty(dims)
    dims = 2;
end
if nargin < 2 || isempty(group_data)
    group_data = repmat({'n/a'},size(plot_data,1),1);
end

% verify marker_size input is correct dimensions
if length(marker_size) == 1
    marker_size = repmat(marker_size,size(plot_data,1),1);
elseif length(marker_size) ~= size(plot_data,1)
    error('MARKER_SIZE must be a scalar or a vector of length equal to the number of rows in PLOT_DATA.');
else
    marker_size = (marker_size - min(marker_size))./range(marker_size).*range(marker_size_range) + marker_size_range(1);
end

% verify group_data input is correct dimensions
if length(group_data) ~= size(plot_data,1)
    error('Dimensions of PLOT_DATA and GROUP_DATA are inconsistent.');
end

% try to convert GROUP_DATA to cell array of strings
if ~iscell(group_data)
    group_data = arrayfun(@num2str,group_data,'UniformOutput',false);
end 

group_names = unique(group_data);  % get unique list of group names or IDs
if ischar(cmap)
    eval(['c = ',cmap,'(',num2str(length(unique(group_names))),');']);
else
    c = cmap;
end

if dims == 2
    
    % 2-D scatter plot
    text_coords = [];
    for i = 1:length(group_names)
        
        % determine which points belong to the current group
        ind = ismember(group_data,group_names{i});
        
        % determine the center of the points
        cent_x = mean(plot_data(ind,1));
        cent_y = mean(plot_data(ind,2));
        
        % zip the center into the points
        x = [plot_data(ind,1),repmat(cent_x,sum(ind),1)]'; x = x(:);
        y = [plot_data(ind,2),repmat(cent_y,sum(ind),1)]'; y = y(:);
        
        % plot data, coloring points by group
        if ( marker_fill )
            h = scatter(plot_data(ind,1),plot_data(ind,2),marker_size(ind),c(i,:),marker_symbol,'filled');
        else
            h = scatter(plot_data(ind,1),plot_data(ind,2),marker_size(ind),c(i,:),marker_symbol);
        end
        
        % adjust marker settings
        if strcmpi(marker_outline,'none')
            set(h,'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        elseif strcmpi(marker_outline,'group')
            set(h,'MarkerEdgeColor',c(i,:),'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        else
            set(h,'MarkerEdgeColor',marker_outline,'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        end
        hold on

        % add lines to connect point groups
        if ( draw_lines )
            plot(x,y,'-','Color',c(i,:));
        end
        
        % save coordinates of group center to add a text label to afterward
        text_coords = [text_coords;[cent_x,cent_y]];
        
    end
    
    % add labels at the end, so text is on top of points and lines
    if ( label_points )
        text(text_coords(:,1),text_coords(:,2),group_names,'HorizontalAlignment',label_align);
    end

elseif dims == 3
    
    % 3-D scatter plot
    text_coords = [];
    for i = 1:length(group_names)
        
        % determine which points belong to the current group
        ind = ismember(group_data,group_names{i});
        
        % determine the center of the points
        cent_x = mean(plot_data(ind,1));
        cent_y = mean(plot_data(ind,2));
        cent_z = mean(plot_data(ind,3));
        
        % zip the center into the points
        x = [plot_data(ind,1),repmat(cent_x,sum(ind),1)]'; x = x(:);
        y = [plot_data(ind,2),repmat(cent_y,sum(ind),1)]'; y = y(:);
        z = [plot_data(ind,3),repmat(cent_z,sum(ind),1)]'; z = z(:);
        
        % plot data, coloring points by group
        if ( marker_fill )
            h = scatter3(plot_data(ind,1),plot_data(ind,2),plot_data(ind,3),marker_size(ind),c(i,:),'filled');
        else
            h = scatter3(plot_data(ind,1),plot_data(ind,2),plot_data(ind,3),marker_size(ind),c(i,:));
        end
        
        if strcmpi(marker_outline,'none')
            set(h,'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        elseif strcmpi(marker_outline,'group')
            set(h,'MarkerEdgeColor',c(i,:),'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        else
            set(h,'MarkerEdgeColor',marker_outline,'MarkerEdgeAlpha',marker_edge_alpha,'MarkerFaceAlpha',marker_face_alpha);
        end
        hold on
        
        % add lines to connect point groups
        if ( draw_lines )
            plot3(x,y,z,'-','Color',c(i,:)); grid on;
        end
        
        % save coordinates of group center to add a text label to afterward
        text_coords = [text_coords;[cent_x,cent_y,cent_z]];
        
    end
    
    % add labels at the end, so text is on top of points and lines
    if ( label_points )
        text(text_coords(:,1),text_coords(:,2),text_coords(:,3),group_names,'HorizontalAlignment',label_align);
    end
    
    grid on
    
else
    error('Invalid number of dimensions. DIMS must be 2 or 3.');
end





