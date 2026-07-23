using GridapGmsh: gmsh

gmsh.initialize()

gmsh.model.add("Geo")

lc=0.05

p1 = gmsh.model.geo.addPoint(1,0,0,lc)
p2 = gmsh.model.geo.addPoint(0,1,0,lc)
p3 = gmsh.model.geo.addPoint(-1,0,0,lc)
p4 = gmsh.model.geo.addPoint(0,-1,0,lc)
center = gmsh.model.geo.addPoint(0,0,0,lc)

l1 = gmsh.model.geo.addCircleArc(p1,center,p2)
l2 = gmsh.model.geo.addCircleArc(p2,center,p3)
l3 = gmsh.model.geo.addCircleArc(p3,center,p4)
l4 = gmsh.model.geo.addCircleArc(p4,center,p1)

circle = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4])
disk = gmsh.model.geo.addPlaneSurface([circle])

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1,[l1,l2,l3,l4],-1,"boundary")
gmsh.model.addPhysicalGroup(2,[disk],-1,"domain")

gmsh.model.mesh.generate(2)

gmsh.write(joinpath(@__DIR__, "geo-circle.msh"))

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()







