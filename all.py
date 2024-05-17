import cone, ellipsoid, hyperbolicparaboloid, hyperboloid, paraboliccylinder


for w in [1.6, 2, 2.8]:
    cone.generate(0.5, 0.5, width_mm=w)
    cone.generate(0.5, 1, width_mm=w)
    cone.generate(2, 1, width_mm=w)

    ellipsoid.generate(width_mm=w)
    hyperbolicparaboloid.generate(width_mm=w)
    hyperboloid.generate(width_mm=w)
    paraboliccylinder.generate(width_mm=w)
