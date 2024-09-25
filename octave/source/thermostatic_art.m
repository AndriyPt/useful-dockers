clear *

h = 0.1

function xy_new = Deform(xy)
xy_new = [xy(:,1).*cos(xy(:,2)), xy(:,1).*sin(xy(:,2))];
endfunction

function u = f_u_exact(xy)
u = exp(-2*xy(:,2));
endfunction

function u = f_DDu_exact(xy)
u = -4*(1+xy(:,1).^2).*exp(-2*xy(:,2));
endfunction

function a = f_a(xy)
a = 1 + xy(:,1).^2;
endfunction

FEMmesh = CreateMeshTriangle('Test',[1,0,-1;2,0,-1;2,pi/2,-2;1,pi/2,-1],h^2);
FEMmesh = MeshDeform(FEMmesh,'Deform');

figure(1); 
FEMtrimesh(FEMmesh);
FEMmesh = MeshUpgrade(FEMmesh,'cubic');

u = BVP2Dsym(FEMmesh,'f_a',0,'f_DDu_exact','f_u_exact',0,0);

figure(2); 
FEMtrimesh(FEMmesh,u);
xlabel('x'); ylabel('y'); zlabel('u'); 
view([-150,30]);
u_exact = f_u_exact(FEMmesh.nodes);
L2Error = sqrt(FEMIntegrate(FEMmesh,(u-u_exact).^2))

