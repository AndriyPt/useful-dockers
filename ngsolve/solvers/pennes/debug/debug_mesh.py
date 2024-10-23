from netgen.csg import *

brick = OrthoBrick(Pnt(-2,-2,-2),Pnt(2,2,2)).bc('outer')
sphere = Sphere(Pnt(0,0,0),1).bc('sphere')

# halfsphere = sphere * Plane(Pnt(0,0,0),Vec(1,0,0)).bc('plane')
# halfsphere = Plane(Pnt(0,0,0),Vec(1,0,0)).bc('plane')
halfsphere = sphere
box = brick-sphere
geo = CSGeometry()
# geo.Add(box.col([1,0,0]).transp())
# geo.Add(halfsphere.col([0,0,1]),bcmod=[(box,"halfsphere")])
geo.Add(halfsphere)
geo.Draw()

ngmesh = geo.GenerateMesh()
mesh = Mesh(ngmesh)
mesh.GetBoundaries()
