using GridapGmsh: gmsh

gmsh.initialize()

gmsh.model.add("Geo_rectangle")

lc=0.05

p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(pi, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(pi, pi, 0, lc)
p4 = gmsh.model.geo.addPoint(0, pi, 0, lc)

l1 = gmsh.model.geo.addLine(p1, p2) 
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

rectangle = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
rectangle_area = gmsh.model.geo.addPlaneSurface([rectangle])

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], -1, "boundary")
gmsh.model.addPhysicalGroup(2, [rectangle_area], -1, "domain")

gmsh.model.mesh.generate(2)

gmsh.write(joinpath(@__DIR__, "geo-rectangle.msh"))

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()
