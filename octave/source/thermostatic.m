 node_y = [1;1.2;1.5;1.8;2]*ones(1,11);
 node_x = ones(5,1)*[1,1.05,1.1,1.2, ...
           1.3,1.5,1.7,1.8,1.9,1.95,2];
 nodes = [node_x(:), node_y(:)];

 [h,w] = size (node_x);
 elems = [];
 for idx = 1:w-1
   widx = (idx-1)*h;
   elems = [elems; ...
     widx+[(1:h-1);(2:h);h+(1:h-1)]'; ...
     widx+[(2:h);h+(2:h);h+(1:h-1)]' ];
 endfor

 E = size (elems,1); # No. of simplices
 N = size (nodes,1); # No. of vertices
 D = size (elems,2); # dimensions+1 

 ## Element conductivity
  conductivity = [1*ones(1,16), ...
         2*ones(1,48), 1*ones(1,16)];

  ## Connectivity matrix
  C = sparse ((1:D*E), reshape (elems', ...
         D*E, 1), 1, D*E, N);

  ## Calculate system matrix
  Siidx = floor ([0:D*E-1]'/D) * D * ...
         ones(1,D) + ones(D*E,1)*(1:D) ;
  Sjidx = [1:D*E]'*ones (1,D);
  Sdata = zeros (D*E,D);
  dfact = factorial (D-1);
  for j = 1:E
     a = inv ([ones(D,1), ...
         nodes(elems(j,:), :)]);
     const = conductivity(j) * 2 / ...
         dfact / abs (det (a));
     Sdata(D*(j-1)+(1:D),:) = const * ...
         a(2:D,:)' * a(2:D,:);
  endfor
  ## Element-wise system matrix
  SE = sparse(Siidx,Sjidx,Sdata);
  ## Global system matrix
  S = C'* SE *C; 
  
  
   ## Dirichlet boundary conditions
  D_nodes = [1:5, 51:55];
  D_value = [10*ones(1,5), 20*ones(1,5)];

  V = zeros (N,1);
  V(D_nodes) = D_value;
  idx = 1:N; # vertices without Dirichlet
             # boundary condns
  idx(D_nodes) = [];

  ## Neumann boundary conditions.  Note that
  ## N_value must be normalized by the
  ## boundary length and element conductivity
  N_nodes = [];
  N_value = [];

  Q = zeros (N,1);
  Q(N_nodes) = N_value;

  V(idx) = S(idx,idx) \ ( Q(idx) - ...
            S(idx,D_nodes) * V(D_nodes)); 
            
  elemx = elems(:,[1,2,3,1])';
  xelems = reshape (nodes(elemx, 1), 4, E);
  yelems = reshape (nodes(elemx, 2), 4, E);
  velems = reshape (V(elemx), 4, E);
  plot3 (xelems,yelems,velems,"k");
  print "grid.eps"; 
  
  