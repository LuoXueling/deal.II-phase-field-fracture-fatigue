// Monopile OWT -- unstructured C3D4 mesh source geometry.
//
// The Abaqus BaseSolidRevolve route cannot produce an unstructured mesh: the
// swept solid carries a parametric (theta, z) grid that the free tet mesher
// inherits, so the surface comes out as continuous circumferential rings and
// axial columns regardless of seeding. Perturbing node positions does not help
// either -- jitter moves nodes but cannot change connectivity, so the rings
// survive (measured: 26 distinct surface z-levels vs 1648 for this mesh).
//
// Gmsh builds the same solid in the OpenCASCADE kernel and triangulates it with
// Frontal-Delaunay, which has no sweep parameterization to inherit.
SetFactory("OpenCASCADE");

z_mud = -30.0; z_mwl = 0.0; z_hub = 87.0;
r_p_out = 2.75; r_p_in = 2.69;   // pile:  D=5.5, t=0.060
r_t_out = 1.95; r_t_in = 1.89;   // tower: D=3.9, t=0.060

// The EXACT closed profile from the Abaqus sketch, in the (r,z) half-plane,
// revolved 2*pi about z. One solid -- no booleans, which self-intersected at
// the pile/tower junction and produced "PLC Error: segment and facet intersect".
Point(1) = {r_p_in,  0, z_mud};
Point(2) = {r_p_in,  0, z_mwl};
Point(3) = {r_t_in,  0, z_hub};
Point(4) = {r_t_out, 0, z_hub};
Point(5) = {r_p_out, 0, z_mwl};
Point(6) = {r_p_out, 0, z_mud};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4};
Line(4) = {4,5}; Line(5) = {5,6}; Line(6) = {6,1};
Curve Loop(1) = {1,2,3,4,5,6};
Plane Surface(1) = {1};
Extrude { {0,0,1}, {0,0,0}, 2*Pi } { Surface{1}; }

// Refinement is wanted ONLY on the load-facing (-X) half of the two loaded
// ends. Everywhere else -- the whole leeward half, and the mid-tower on both
// sides -- is held at the coarse size.
//
// Two multiplied factors:
//
//  AXIAL   size = cap + (mid-cap)*(d/2)^0.20 , capped at mid
//          d = distance to the nearer loaded end. cap = 0.005 -> 12 elements
//          through the 0.060 wall at the loaded face.
//          The ramp MUST grow fast: 0.005 tets cost ~1.4e6 elements per
//          centimetre of reach (the -X half of both caps is ~1 m^3 per metre,
//          and a 0.005 tet is 1.5e-8 m^3). Exponent 0.20 over a 2 m reach
//          gives 0.005 at the face, 0.18 by 1 cm, 0.28 by 10 cm, coarse by 2 m.
//
//  SECTOR  the axial refinement is confined to a ~5 m ARC centred on the -X
//          apex, not the whole load-facing half. At the pile (r=2.75) the
//          load-facing half-circumference is 8.64 m, so 5 m is +/-52 deg about
//          -X rather than the full +/-90 deg.
//          MathEval has no ">" operator and no Step(), so the 0/1 switch is
//          built from a normalised ratio: with u = -x/r (u=+1 at the -X apex,
//          0 at the seam, -1 on the leeward apex), the refined sector is
//          u > cos(halfangle). off = Max(0, sign(cos_ha - u)) evaluates to 1
//          outside the sector and 0 inside, and the term (mid-cap)*off lifts
//          everything outside the arc straight to h_mid.
//
// mid is capped at 0.50: the wall is only 0.060 thick, and a target much
// larger than the wall makes the surface mesh self-intersect -- at mid=1.0
// gmsh reports a PLC error and emits ZERO volume elements.
// ---- knobs -----------------------------------------------------------
// REFINE = 1 : graded mesh, clustered on the load-facing half of both caps
// REFINE = 0 : UNIFORM mesh at h_mid everywhere -- the verification mesh.
//              Same geometry, same surface ids, same loads; only the size
//              field differs, so results can be compared directly.
// Override from the command line without editing this file:
//   gmsh monopile_gmsh.geo -3 -setnumber REFINE 0 -format inp -o uniform.inp
DefineConstant[ REFINE = 1 ];
DefineConstant[ h_cap  = 0.01  ];   // size at the loaded faces (refined only)
DefineConstant[ h_mid  = 0.50  ];   // size everywhere else / uniform size
// A band of thickness `plateau` measured in from each loaded face is held at
// h_cap, so there really is a layer of ~h_cap elements there rather than the
// size starting to grow the instant it leaves the surface. The ramp to h_mid
// then runs from the end of the plateau out to `reach`.
DefineConstant[ plateau = 0.02 ];   // thickness held at h_cap (m)
// reach/pw set the AXIAL extent of the refinement, measured from each loaded
// end face. A short reach (0.05) makes the decay a CLIFF -- everything past
// 5 cm is already at h_mid, so the mesh reads as refined-then-abrupt rather
// than smoothly decaying. 1.0 m with a linear ramp decays smoothly:
//   0.005 at the face, 0.018 at 2 cm, 0.070 at 10 cm, 0.13 at 20 cm,
//   0.33 at 50 cm, 0.49 at 75 cm, 0.65 at 1 m.
// The longer reach is affordable because arc_len below confines the fine zone
// to a narrow sector instead of the whole load-facing half.
DefineConstant[ reach  = 1.0   ];   // axial extent of the refinement (m)
DefineConstant[ pw     = 1.0   ];   // ramp shape after the plateau (1.0 = linear)
DefineConstant[ arc_len = 5.0  ];   // circumferential extent of refinement (m)

// The two loaded ends are treated DIFFERENTLY.
//
// MUDLINE (z=-30): the crack-resolution end. Fine (h_cap), confined to the
//   load-facing arc, graded out over `reach` -- the knobs above.
//
// HUB (z=+87): only needs to resolve the applied RNA traction, not a crack.
//   Refining one side of it is wrong here: the RNA weight and thrust act over
//   the WHOLE annulus, so the refinement has to be all-round, and it only has
//   to be modestly finer than h_mid. At h_mid=0.65 the hub ring (12.25 m
//   circumference) carries just 19 elements around, too few to represent the
//   traction; h_hub=0.40 gives ~31.
//   Converged: a sweep of 0.40/0.30/0.20/0.15 moved the tip displacement by
//   only 0.2%/0.3%/0.02%, and the reaction force not at all (-2.189 MN at every
//   level), so 0.40 already resolves the traction. Refining further just costs
//   elements -- the residual error against the hex reference is set by h_mid,
//   not by the hub.
DefineConstant[ h_hub     = 0.40 ];  // size at the hub cap (all round)
DefineConstant[ hub_reach = 1.0  ];  // axial extent of the hub refinement (m)

If (REFINE)
  // Half-angle of the refined arc, taken at the pile radius (the widest part,
  // so a 5 m arc there is >=5 m everywhere the pile is loaded).
  cos_ha = Cos(arc_len / (2.0 * 2.75));
  // MUDLINE term: fine, arc-limited, distance measured from z=-30 only.
  //   term 1: h_cap
  //   term 2: axial ramp, shifted by the plateau so the size stays exactly at
  //           h_cap for the first `plateau` metres and only then climbs.
  //   term 3: sector gate -- 0 inside the arc, 1 outside, so outside jumps to
  //           h_mid. u = -x/r is the cosine of the angle from the -X apex.
  Field[1] = MathEval;
  Field[1].F = Sprintf("Min(%g + %g*( (Max(0, Abs(z+30) - %g)/%g)^%g ) + %g*Max(0, ((%g) - (-x/Sqrt(x*x+y*y+1e-12)))/Sqrt(((%g) - (-x/Sqrt(x*x+y*y+1e-12)))^2 + 1e-12)), %g)",
                       h_cap, h_mid - h_cap, plateau, reach - plateau, pw,
                       h_mid - h_cap, cos_ha, cos_ha, h_mid);

  // HUB term: modest, ALL ROUND (no sector gate), distance from z=+87 only.
  Field[2] = MathEval;
  Field[2].F = Sprintf("Min(%g + %g*(Abs(z-87)/%g), %g)",
                       h_hub, h_mid - h_hub, hub_reach, h_mid);

  // The mesh takes the finer of the two requests everywhere.
  Field[3] = Min;
  Field[3].FieldsList = {1, 2};
  BG = 3;
Else
  // Uniform. Note h_mid is the SAME value the refined mesh coarsens to, so
  // the two meshes agree away from the caps and differ only where the
  // refinement acts -- which is what makes them comparable.
  Field[1] = MathEval;
  Field[1].F = Sprintf("%g", h_mid);
  BG = 1;
EndIf
Background Field = BG;

// Take the size ONLY from the field above -- otherwise gmsh adds its own
// curvature/point-based refinement and the grading is no longer what we asked.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.Algorithm = 6;      // 2D Frontal-Delaunay
Mesh.Algorithm3D = 1;    // 3D Delaunay
Mesh.ElementOrder = 1;   // linear tets -> C3D4
Mesh.Optimize = 1;
