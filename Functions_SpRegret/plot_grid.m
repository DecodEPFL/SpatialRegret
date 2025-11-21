function h = plot_grid(G, shape, color, graph_are_undirected, show_labels)
    if nargin<3, color = 'r';          end
    if nargin<4, graph_are_undirected = false; end
    if nargin<5, show_labels = true;   end

    assert(numel(G)==prod(shape)^2, ...
      "Error. Shape does not match with the size of G")
    m = shape(1);  n = shape(2);

    G_tilde = bin(G - eye(size(G,1)));
    graphObj = digraph(G_tilde);
    if graph_are_undirected
      graphObj = graph(graphObj.Edges.EndNodes(:,1), ...
                       graphObj.Edges.EndNodes(:,2));
      if ismultigraph(graphObj)
        graphObj = simplify(graphObj);
      end
    end

    [Y,X] = meshgrid(m:-1:1, 1:n);
    coords = [X(:), Y(:)];

    % build labels
    N = numnodes(graphObj);
    if show_labels
        labels = string(1:N);
    else
        labels = strings(N,1);
    end

    h = plot(graphObj, ...
        'XData', coords(:,1), 'YData', coords(:,2), ...
        'NodeColor','w', 'NodeLabel', labels, ...
        'MarkerSize',6, 'LineWidth',3, 'EdgeColor',color);

    axis equal
    set(gca, 'XTick',[], 'YTick',[], 'XColor','none','YColor','none')
    box on

    hold on
    size_circles = 0.1;
    for k = 1:numel(h.XData)
      rectangle('Position',...
        [h.XData(k)-size_circles/2, h.YData(k)-size_circles/2, ...
         size_circles, size_circles], ...
        'Curvature',[1,1], 'EdgeColor','k','LineWidth',1);
    end
    hold off
end
