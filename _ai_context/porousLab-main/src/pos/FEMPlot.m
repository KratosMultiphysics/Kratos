%% FEMPlot Class
% This class provides methods to visualize graphical results from Finite 
% Element Method (FEM) analysis. It includes functionalities to plot 
% meshes, evaluate fields along segments, and compute bounding boxes for 
% models.
%
%% Methods
% * *getBoundingBox*: Computes the bounding box of the model with a 
%                     tolerance.
% * *plotMesh*: Plots the continuum elements of the mesh. It also
%               configures the figure wth default font, size and colormap.
% * *plotFieldAlongSegment*: Plots the specified field along a segment
%                            defined by two points.
% *getElementsPatches*: It combines faces, vertices and vertex data from
%                       all the elements in the model, returning combined
%                       faces connectivity, vertices coordinates and
%                       combined vertex data.
%
%% Authors
% * Danilo Cavalcanti (dborges@cimne.upc.edu)
%
%% Class Definition
classdef FEMPlot < handle
    %% Public properties
    properties
        model = []; % Handle to an object of the model class
    end
    %% Constant properties
    properties (Constant, Access = private)
        DEFAULT_FONT_NAME           = 'Times';        % Default font name for plots
        DEFAULT_FONT_SIZE           = 16;             % Default font size for plots
        DEFAULT_COLOR               = '-b';           % Default color and line style for plots
        DEFAULT_LINE_WIDTH          = 1.5;            % Default line width for plots
        DEFAULT_FACE_COLOR          = 'interp';       % Default face color scheme for plotting
        DEFAULT_MARKER_TYPE         = 'None';         % Default type of marker
        DEFAULT_MARKER_FACE_COLOR   = 'k';            % Default color of the face of the marker
        DEFAULT_EDGES_THICKNESS     = 0.7;            % Default thickness of the lines
        DEFAULT_COLORMAP_TYPE       = 'jet';          % Default colormap scheme
    end
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function this = FEMPlot(model)
            if (nargin > 0)
                this.model = model;
            end
        end
    end
    
    %% Public methods
    methods
        %------------------------------------------------------------------
        % Computes the bounding box of the model with a tolerance.
        % Note:
        %   - The tolerance is calculated as 10% of the absolute difference 
        %     between the maximum and minimum x-coordinates of the model's nodes.
        %
        function bbox = getBoundingBox(this)
            minX = Inf; minY = Inf; maxX = -Inf; maxY = -Inf;
            for el = 1:this.model.nelem
                res = this.model.element(el).type.result;
                minX = min(minX,min(res.vertices(:,1)));
                minY = min(minY,min(res.vertices(:,2)));
                maxX = max(maxX,max(res.vertices(:,1)));
                maxY = max(maxY,max(res.vertices(:,2)));
            end
            tolx = 0.10*abs(max(this.model.NODE(:,1)) - min(this.model.NODE(:,1)));
            toly = 0.10*abs(max(this.model.NODE(:,1)) - min(this.model.NODE(:,1)));
            bbox(1) = minX - tolx;
            bbox(2) = minY - toly;
            bbox(3) = maxX + tolx;
            bbox(4) = maxY + toly;
        end

        %------------------------------------------------------------------
        % Plots the continuum elements of the mesh
        function plotMesh(this, ax)
            % Verifica se o eixo foi passado
            if nargin < 2 || isempty(ax)
                figure;         % Cria nova figura
                ax = gca;       % Usa o eixo atual
            else
                axes(ax);       % Define o eixo alvo
                cla(ax);        % Limpa o conteúdo
            end
        
            hold(ax, 'on');
        
            % Get patches
            [Faces, Vertices, VertexData] = this.getElementsPatches();
        
            % Plot the continuum elements
            patch(ax, 'Faces', Faces, ...
                      'Vertices', Vertices, ...
                      'FaceVertexCData', VertexData, ...
                      'FaceColor', this.DEFAULT_FACE_COLOR, ...  
                      'LineWidth', this.DEFAULT_EDGES_THICKNESS, ...
                      'LineStyle', '-', ...
                      'Marker', this.DEFAULT_MARKER_TYPE);
        
            % Set the colormap
            colormap(ax, this.DEFAULT_COLORMAP_TYPE);
        
            % Get the bounding box
            bbox = this.getBoundingBox();
            xlim(ax, [bbox(1) bbox(3)]);
            ylim(ax, [bbox(2) bbox(4)]);
        
            % Formatting the figure
            box(ax, 'on');
            grid(ax, 'on');
            axis(ax, 'equal');
            set(ax, 'FontName', this.DEFAULT_FONT_NAME);
            set(ax, 'FontSize', this.DEFAULT_FONT_SIZE);
        end

        %------------------------------------------------------------------
        % Plots the field along a segment defined by two points Xi and Xf
        function plotFieldAlongSegment(this, field, Xi, Xf, npts, axisPlot, ax)

            % Coordinates of the points where the given field is going
            % to be evaluated
            X = [linspace(Xi(1), Xf(1),npts);linspace(Xi(2), Xf(2),npts)]';

            % Initialize vector at these points
            F = zeros(size(X,1),1);   % Field vector
            S = zeros(size(X,1),1);   % Longitudinal coordinate vector

            % Calculate the pressure field in the points of the segment
            for i = 1:npts
                
                % Longitudinal coordinate of the point X(i) wrt to Xi
                S(i) = sqrt((X(i,1) - Xi(1))*(X(i,1) - Xi(1)) + (X(i,2) - Xi(2))*(X(i,2) - Xi(2)));
                
                % Find in each element this point is inside
                elem = findElementInMesh(this.model.NODE, this.model.ELEM, X(i,:));
                
                % Calculate the given field in the point X using the
                % shape function of the element
                F(i) = this.model.evaluateField(field, elem, X(i,:));
                
            end

            % Initialize, plot and configure the figure
            cla(ax);
            hold(ax, 'on');
            box(ax, 'on');
            grid(ax, 'on');
            
            % Plota conforme o eixo escolhido
            if strcmp(axisPlot, 'y')
                plot(ax, F, S, this.DEFAULT_COLOR, 'LineWidth', this.DEFAULT_LINE_WIDTH);
                xlabel(ax, field);
                ylabel(ax, 'Longitudinal distance (m)');
            elseif strcmp(axisPlot, 'x')
                plot(ax, S, F, this.DEFAULT_COLOR, 'LineWidth', this.DEFAULT_LINE_WIDTH);
                ylabel(ax, field);
                xlabel(ax, 'Longitudinal distance (m)');
            end
            
            % Aplica formatação
            set(ax, 'FontName', this.DEFAULT_FONT_NAME);
            set(ax, 'FontSize', this.DEFAULT_FONT_SIZE);

        end

        %------------------------------------------------------------------
        % Plots the field along a discontinuity
        function plotFieldAlongDiscontinuity(this, S, F, field, axisPlot, ax)

            % Initialize, plot and configure the figure
            cla(ax);
            hold(ax, 'on');
            box(ax, 'on');
            grid(ax, 'on');

            % Set the limits 
            fmin = min(F);
            fmax = max(F);
            pad = max(0.02*(fmax-fmin),max(0.001*mean(F),0.01));
            
            % Plota conforme o eixo escolhido
            if strcmp(axisPlot, 'y')
                plot(ax, F, S, this.DEFAULT_COLOR, 'LineWidth', this.DEFAULT_LINE_WIDTH);
                xlabel(ax, field);
                ylabel(ax, 'Longitudinal distance (m)');
                xlim(ax, [fmin-pad, fmax+pad]);
            elseif strcmp(axisPlot, 'x')
                plot(ax, S, F, this.DEFAULT_COLOR, 'LineWidth', this.DEFAULT_LINE_WIDTH);
                ylabel(ax, field);
                xlabel(ax, 'Longitudinal distance (m)');
                ylim(ax, [fmin-pad, fmax+pad]);
            end
            
            % Aplica formatação
            set(ax, 'FontName', this.DEFAULT_FONT_NAME);
            set(ax, 'FontSize', this.DEFAULT_FONT_SIZE);         

        end

        %------------------------------------------------------------------
        % Combines faces, vertices, and vertex data from all elements in the model.
        function [allFaces, allVertices, allVertexData] = getElementsPatches(this)
            
            % Initialize combined matrices
            allFaces      = {};      % Combined faces connectivity
            allVertices   = [];      % Combined vertices coordinates
            allVertexData = [];      % Combined vertex data
            vertexOffset  = 0;       % Offset for face indices due to combined vertices
            
            % Loop through all elements to build combined matrices
            for el = 1:this.model.nelem
                % Get the current element and its results
                element = this.model.element(el).type;
                res = element.result;
                
                % Adjust face indices for combined vertices
                faces = res.faces + vertexOffset;

                % Store the faces as a cell array to handle different numbers of faces
                allFaces{end+1} = faces;
                
                % Append to combined matrices
                allVertices   = [allVertices; res.vertices];
                allVertexData = [allVertexData; res.vertexData];
                
                % Update vertex offset for the next element
                vertexOffset = vertexOffset + size(res.vertices, 1);
            end

            % After processing, you can convert allFaces into a matrix if needed
            % by padding with NaNs or zeros to ensure uniform row sizes
            maxFaces = max(cellfun(@numel, allFaces));  % Find the maximum number of faces
            paddedFaces = NaN(numel(allFaces), maxFaces);  % Pre-allocate matrix with NaNs
            
            for i = 1:numel(allFaces)
                paddedFaces(i, 1:numel(allFaces{i})) = allFaces{i};  % Pad each row
            end
            
            allFaces = paddedFaces;  % Use the padded faces matrix
        end
    end
end