using GridapGmsh: gmsh

gmsh.initialize()

gmsh.model.add("Geo_ellipse")

lc=0.05

p1 = gmsh.model.geo.addPoint(2, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(0, 1, 0, lc)
p3 = gmsh.model.geo.addPoint(-2, 0, 0, lc)
p4 = gmsh.model.geo.addPoint(0, -1, 0, lc)
center = gmsh.model.geo.addPoint(0, 0, 0, lc)

l1 = gmsh.model.geo.addEllipseArc(p1, center, p1, p2)
l2 = gmsh.model.geo.addEllipseArc(p2, center, p1, p3)
l3 = gmsh.model.geo.addEllipseArc(p3, center, p1, p4)
l4 = gmsh.model.geo.addEllipseArc(p4, center, p1, p1)

ellipse = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
ellipse_area = gmsh.model.geo.addPlaneSurface([ellipse])

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l1, l2, l3, l4], -1, "boundary")
gmsh.model.addPhysicalGroup(2, [ellipse_area], -1, "domain")

gmsh.model.mesh.generate(2)

gmsh.write(joinpath(@__DIR__, "geo-ellipse.msh"))

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()
