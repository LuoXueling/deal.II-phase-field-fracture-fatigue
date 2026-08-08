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
// HALF MODEL. The profile lies in the +X half-plane, so revolving by Pi
// (not 2*Pi) sweeps it through +Y only and leaves the y=0 plane as two flat
// exposed faces -- the sweep start (+X side) and the sweep end (-X side).
//
// This is legitimate because EVERY ingredient of the problem is symmetric
// about y=0: the geometry is a surface of revolution, the loads in
// monopile_*_boundary.txt all act along +X ("triangularneumann 1 0 0"), and
// the refined sector is centred on the -X apex, which lies ON the y=0 plane.
// So the solution must be symmetric, and half the domain carries it.
//
// The cut plane needs a SYMMETRY (roller) condition: u_y = 0, with u_x and
// u_z free. In monopile_*_boundary.txt that is a single line per exposed
// face id, constraining component 1 only --
//     <id> dirichlet 1 0.0
// which AbstractField::setup_dirichlet_boundary_condition applies through
// fields.component_masks, so it fixes y alone and leaves x/z unconstrained.
// Do NOT use "<id> dirichlet 0 0.0 / 1 0.0 / 2 0.0" there: that is an
// encastre, not a roller, and would clamp the crack face shut.
//
// Halving the domain frees up budget, so arc_len below is doubled to keep
// the same absolute refined width in metres as the full model had.
Extrude { {0,0,1}, {0,0,0}, Pi } { Surface{1}; }

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
// h_cap is a TARGET, not a bound: gmsh's Delaunay refinement lands ABOVE the
// requested size, so h_cap must be pre-divided by an inflation factor. That
// factor depends on the GRADING and on whether netgen optimisation runs -- NOT
// on arc_len or on h_cap itself. Measured history:
//   pw=0.5, Delaunay only        -> ~1.32x  (1.27, 1.314, 1.319, 1.323)
//   pw=1.5, Delaunay only        -> ~1.645x (1.641/1.642/1.643 at arc 1.2 for
//                                   h_cap 0.0038/0.0030/0.0026; 1.648/1.647
//                                   at arc 0.5 -- arc- and h_cap-independent)
//   pw=2.5, + OptimizeNetgen     -> ~1.575x (netgen redistributes nodes after
//                                   the size field is applied, so it shifts
//                                   the factor DOWN from the 1.645 above)
//
// The two earlier factors were measured on tets whose centroid lies inside the
// PLATEAU, where the field is exactly h_cap. h_cap below is instead calibrated
// directly against the quantity that matters -- the median MAX EDGE of tets in
// the refined CORE (inside arc_len, inside the plateau) -- because that is the
// length the phase field actually has to resolve:
//   h_cap = 0.00300 -> core 0.00473  (94.5% of 0.005, 4.23 per l_phi), 276325 tets
//   h_cap = 0.00317 -> core 0.00499  (99.8% of 0.005, 4.01 per l_phi), 241218 tets
// Note the MEAN edge of a tet runs ~18% below its max edge (0.00400 vs 0.00473
// at h_cap=0.0030), so "element size" is ambiguous by that much -- these are
// all max-edge figures.
//
// NOTE: parameters/Monopile*.prm asks for cells <= 0.2*l_phi = 0.004, i.e. 5
// per l_phi. The value below targets 0.005 (4.0 per l_phi) and so is one step
// coarser than that criterion; h_cap = 0.0025 would satisfy it, at roughly
// twice the element count.
DefineConstant[ h_cap  = 0.00317 ]; // request; core max-edge comes out ~0.00499
DefineConstant[ h_mid  = 0.65  ];   // size everywhere else / uniform size
// A band of thickness `plateau` measured in from each loaded face is held at
// h_cap, so there really is a layer of ~h_cap elements there rather than the
// size starting to grow the instant it leaves the surface. The ramp to h_mid
// then runs from the end of the plateau out to `reach`.
DefineConstant[ plateau = 0.025 ];  // thickness held at h_cap (m)
// reach/pw set the AXIAL extent of the refinement, measured from each loaded
// end face. A short reach (0.05) makes the decay a CLIFF -- everything past
// 5 cm is already at h_mid, so the mesh reads as refined-then-abrupt rather
// than smoothly decaying.
//
// pw is the ramp EXPONENT, and it must be >= 1. At pw < 1 the ramp leaves the
// plateau with INFINITE slope -- pw=0.5 grew the requested size 4.6x within
// half a millimetre of the plateau edge, which gmsh cannot honour, so it
// emitted tets that kept the fine in-plane size inherited from the surface
// but stretched axially to meet the coarse request above. That is the
// "refined on the bottom face but elongated axially" failure.
//
// pw must also RISE with h_cap: the ramp spans a size ratio of h_mid/h_cap
// over a fixed `reach`, so making h_cap finer steepens the gradient per
// element and needs a flatter exponent to compensate. Measured at
// h_cap=0.0030, arc_len=0.5 (p99 aspect / p99 max-min edge ratio / size jump
// across the plateau edge / tets):
//   pw=1.5  1.656 / 3.560 / 23.0% / 178k
//   pw=2.0  1.544 / 3.019 /  9.0% / 213k
//   pw=2.5  1.542 / 2.961 /  2.8% / 258k   <- knee: near-field saturates here
//   pw=3.0  1.519 / 2.894 /  0.9% / 310k
//   pw=4.0  1.539 / 2.887 /  0.3% / 414k
// Past 2.5 the near-field metrics are flat while the count keeps climbing
// (pw=4 costs 60% more than pw=2.5 for 0.2% quality). The FAR field does keep
// improving though -- median axial dz at 8 cm depth is 0.0157 at pw=2.5 vs
// 0.0093 at pw=3.0 and 0.0049 at pw=4.0 (plateau elements are 0.0035). So if
// damage ever reaches ~10 cm in from the face, pw=3.0 is worth the +20%.
// (`max` aspect ratios are NOT a useful signal here: they track single random
// slivers in an unstructured Delaunay mesh and are non-monotone in pw.)
DefineConstant[ reach  = 0.3   ];   // axial extent of the refinement (m)
DefineConstant[ pw     = 2.5   ];   // ramp exponent after the plateau (>= 1)
// arc_len pays for h_cap: halving the element size costs 8x per unit volume,
// so the fine zone is narrowed to keep the total affordable. At 0.5 m the
// refined sector is only +/-5.2 deg about the -X apex (25 l_phi wide, ~100
// elements across at h=0.005), and the whole mesh is 179511 tets -- only 25%
// more than the 143911-tet mesh it replaces, which had 2x coarser elements.
// WARNING: the sector gate is a 0/1 switch, not a ramp, so the size jumps
// straight to h_mid at the arc edge. If damage nucleates outside +/-5.2 deg
// it grows into coarse elements immediately. Widen arc_len if the phi field
// in the output reaches the edge of the refined sector.
// Doubled from 0.5 when the model was halved: the half model costs half as
// much per unit refined width, so 1.0 m here buys the same absolute coverage
// budget the full model spent on 0.5 m. Note the refined sector STRADDLES the
// y=0 cut plane (it is centred on the -X apex, which lies in that plane), so
// in the half model only HALF of this arc is actually present -- 1.0 m of arc
// means 0.5 m of mesh, +/-10.4 deg measured from -X but existing only on the
// +Y side.
DefineConstant[ arc_len = 1.0  ];   // circumferential extent of refinement (m)
// arc_ramp / pw_arc are the CIRCUMFERENTIAL analogue of reach / pw. arc_len
// is the fully-refined core; outside it the size ramps to h_mid over a band
// arc_ramp metres wide (measured on each side, at the pile radius), with
// exponent pw_arc. Set arc_ramp = 0 to recover the old hard sector gate --
// but see the note on Field[1]: that step is what produced the stray
// high-aspect elements at the arc edge.
// NOTE the total refined footprint is now arc_len + 2*arc_ramp, so raising
// arc_ramp costs elements over the FULL refined depth, not just near the cap.
//
// Set EQUAL to arc_len: the transition band on each side is as long as the
// fully refined core itself, so the size decays over the same angular scale
// it was held constant over, rather than being rushed. Total footprint is
// arc_len + 2*arc_ramp = 3 m (+/-31.2 deg), of which the core is the middle
// third. Measured cost, all with the core guard active and 0 degenerate cells:
//   arc_ramp   tets     core med edge   per l_phi
//     0.3    234218        0.00473        4.23
//     0.5    245242        0.00472        4.24
//     1.0    275710        0.00473        4.23
// The core is untouched by the choice -- Field[4] pins it -- so the extra
// elements all land in the transition, which is the point. 1.0 costs +18%
// over 0.3.
DefineConstant[ arc_ramp = 1.0 ];   // circumferential transition width (m)
DefineConstant[ pw_arc   = 2.5 ];   // circumferential ramp exponent (>= 1)

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
  //
  // The sector term used to be a 0/1 STEP built from sign(): 0 inside the arc,
  // 1 outside, so the requested size jumped h_cap -> h_mid (0.0030 -> 0.65, a
  // 217x discontinuity) across zero distance at the arc edge. That is exactly
  // the pathology `pw` fixes axially, left un-fixed circumferentially, and it
  // is where the stray max-aspect-ratio elements of 14-22 came from while the
  // p99 stayed at 1.5 -- a few tets straddling the edge, stretched to bridge a
  // gradient no mesher can honour.
  //
  // It is now a RAMP over a band `arc_ramp` metres wide, in the same shape as
  // the axial one. In cosine space:
  //     cos_ha       = cos at the edge of the fully refined arc
  //     cos_ha_outer = cos at the outer edge of the transition band
  // and s = (cos_ha - u)/(cos_ha - cos_ha_outer) runs 0 -> 1 across the band
  // (negative inside the arc, >1 beyond it), clamped by Max(0,..)/Min(1,..).
  // MathEval has no conditionals, hence the clamp rather than a branch.
  cos_ha_outer = Cos((arc_len + 2.0 * arc_ramp) / (2.0 * 2.75));
  Field[1] = MathEval;
  Field[1].F = Sprintf("Min(%g + %g*( (Max(0, Abs(z+30) - %g)/%g)^%g ) + %g*( Min(1, Max(0, ((%g) - (-x/Sqrt(x*x+y*y+1e-12)))/(%g)))^%g ), %g)",
                       h_cap, h_mid - h_cap, plateau, reach - plateau, pw,
                       h_mid - h_cap, cos_ha, cos_ha - cos_ha_outer, pw_arc,
                       h_mid);

  // GUARD: the refined region must stay at h_cap no matter what the ramps do.
  // The two ramp terms above are ADDITIVE, so a point that is inside the arc
  // but slightly past the axial plateau -- or vice versa -- picks up a
  // contribution from the other term and drifts coarser than h_cap. This field
  // pins the core box to exactly h_cap and is combined with Min below, so it
  // can only ever make the mesh FINER than Field[1], never coarser.
  //
  // The box is the fully refined arc (|angle| <= half of arc_len) over the
  // axial plateau. Written with the same clamped-ratio trick: g = 0 inside the
  // box and 1 outside, so the value is h_cap inside and h_mid outside (where
  // it is inert, because Field[1] is already <= h_mid there).
  Field[4] = MathEval;
  Field[4].F = Sprintf("%g + %g*Min(1, Max(0, Abs(z+30) - %g)/1e-6 + Max(0, ((%g) - (-x/Sqrt(x*x+y*y+1e-12))))/1e-6)",
                       h_cap, h_mid - h_cap, plateau, cos_ha);

  // HUB term: modest, ALL ROUND (no sector gate), distance from z=+87 only.
  Field[2] = MathEval;
  Field[2].F = Sprintf("Min(%g + %g*(Abs(z-87)/%g), %g)",
                       h_hub, h_mid - h_hub, hub_reach, h_mid);

  // The mesh takes the finer of the requests everywhere. Field[4] is the
  // refined-core guard, so including it here is what makes it a floor: Min
  // can only pull the size DOWN to h_cap inside the core, never up.
  Field[3] = Min;
  Field[3].FieldsList = {1, 2, 4};
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
// MANDATORY, not a quality nicety. Without it the 3D Delaunay pass closes a
// handful of slivers using ONLY mudline-cap surface nodes -- all four vertices
// at z=-30, so the tet is exactly FLAT. deal.II rejects those outright:
//   GridTools::invert_cells_with_negative_measure ... line 992
//   The violated condition was: GridTools::cell_measure(...) > 0
// and every rank aborts before the first solve. At h_cap=0.0030/arc_len=0.5
// there were 5 such cells in 259112 (0.002%) -- rare enough to survive casual
// inspection, fatal on load. Netgen optimisation removes them and is better
// on every other axis at the same time (measured, same .geo, same knobs):
//                      tets   zero-vol  elem size  per l_phi  p99 AR  max AR
//   Delaunay only    259112      5       0.00493     4.05      1.546   14.4
//   + OptimizeNetgen 215217      0       0.00473     4.23      1.383    3.68
// i.e. 17% FEWER elements that are also FINER and better shaped -- the
// optimiser deletes wasted slivers rather than trading resolution away.
// (Algorithm3D=10 (HXT) also yields zero degenerate cells and is cheaper
// still at 172141, but lands at 0.00521 -- 4% coarser than the 0.005 target,
// only 3.84 per l_phi -- so it gives up resolution exactly where it is wanted.)
// NOTE: netgen redistributes nodes AFTER the size field is applied, so it
// changes the h_cap inflation factor from 1.645 to ~1.576. h_cap=0.0030 still
// lands within 5% of 0.005, so the value above needs no retune, but re-measure
// the factor before trusting it at a different h_cap.
Mesh.OptimizeNetgen = 1;
