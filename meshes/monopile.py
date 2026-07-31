# -*- coding: mbcs -*-
#
# Monopile OWT. The mesh is no longer built inside Abaqus: it is imported from
# gmsh (monopile_gmsh.geo), which produces an UNSTRUCTURED C3D4 mesh graded on
# the load-facing side of both loaded ends.
#
# Usage:
#   gmsh monopile_gmsh.geo -3 -format inp -o monopile_gmsh.inp
#   abq2022 cae noGUI=monopile.py [-- <mesh>.inp]
#
# Sections 2 (parameters), 3 (helpers) and 8 (assembly & loads) are unchanged
# from the original hex script. Sections 4-7 -- geometry, partitioning, surface
# picking and hex seeding -- are replaced by the import plus a surface rebuild,
# because an orphan mesh has no geometric faces for findAt to pick.
from abaqus import *
from abaqusConstants import *
from caeModules import *
import regionToolset
import mesh
import assembly
import math
import os
import sys

# -------------------------------------------------------------------
# 1. INITIALIZATION
# -------------------------------------------------------------------
INP_IN = 'monopile_gmsh.inp'          # produced by monopile_gmsh.geo
for _a in sys.argv[1:]:
    if _a.lower().endswith('.inp'):
        INP_IN = _a
        break
if not os.path.exists(INP_IN):
    raise IOError("%s not found -- run:\n"
                  "  gmsh monopile_gmsh.geo -3 -format inp -o %s" % (INP_IN, INP_IN))

model_name = 'OWT_' + os.path.splitext(os.path.basename(INP_IN))[0]
if model_name in mdb.models:
    del mdb.models[model_name]
my_model = mdb.ModelFromInputFile(name=model_name, inputFileName=INP_IN)

# -------------------------------------------------------------------
# 2. PARAMETERS
# -------------------------------------------------------------------
d_water    = 30.0
D_pile     = 5.5
t_pile     = 0.060
z_mudline  = -d_water    # -30.0
z_mwl      = 0.0
z_hub      = 87.0
D_top      = 3.9
D_base     = 5.5
t_tower    = 0.060

# Physics & Loads
rho_steel = 7850.0   
E_steel   = 210.0e9  
nu_steel  = 0.3
g         = 9.81     
V_avg_ref = 8.8      
H_s       = 5.0      
T_wave    = 10.0     
rho_air   = 1.225    
rho_water = 1025.0   
C_D_tower = 0.7
C_D_pile  = 0.7
C_M_pile  = 2.0
alpha     = 0.12     
z_ref     = 10.0     
M_RNA     = 350000.0 
D_rotor   = 126.0    
R_rotor   = D_rotor / 2.0
A_rotor   = math.pi * R_rotor**2
V_rated   = 11.4     

# (The old hex seeding sizes -- max/mid/min_size -- are gone: the mesh now
# comes from gmsh, so element sizing lives in monopile_gmsh.geo.)

# -------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# -------------------------------------------------------------------
def get_wind_speed(z):
    if z <= 0: return 0.0
    return V_avg_ref * (z / z_ref)**alpha

def get_thrust_coeff(V_hub):
    if V_hub <= V_rated: return min(7.0 / V_rated, 1.0) 
    else: return max((7.0 * V_rated**2) / (V_hub**3), 0.005)

def calc_tower_traction_stress(z_top):
    V_z = get_wind_speed(z_top) 
    return (rho_air * C_D_tower / math.pi) * (V_z**2)

def calc_wave_traction_stress():
    u_surf = (math.pi * H_s) / T_wave
    u_dot_surf = (2 * math.pi**2 * H_s) / (T_wave**2)
    drag = 0.5 * rho_water * C_D_pile * D_pile * (u_surf**2)
    inertia = rho_water * C_M_pile * (math.pi * D_pile**2 / 4.0) * u_dot_surf
    f_surf = drag + inertia
    return (2.0 * f_surf) / (math.pi * D_pile)

def calc_rna_tractions():
    V_hub = get_wind_speed(z_hub)
    C_T = get_thrust_coeff(V_hub)
    F_thrust = 0.5 * rho_air * A_rotor * C_T * (V_hub**2)
    r_ext = D_top / 2.0
    r_int = r_ext - t_tower
    A_top = math.pi * (r_ext**2 - r_int**2)
    return F_thrust / A_top, (M_RNA * g) / A_top

def get_outer_radius(z):
    if z <= 0.0:
        return D_pile / 2.0
    else:
        return (D_base - (z / z_hub) * (D_base - D_top)) / 2.0

wave_tr_val = calc_wave_traction_stress()
rna_tr_thrust, rna_tr_weight = calc_rna_tractions()

# -------------------------------------------------------------------
# 4. IMPORTED MESH  (replaces the old GEOMETRY / PARTITIONING / SEEDING)
# -------------------------------------------------------------------
z_struct = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
tower_levels = z_struct + [z_hub]

def get_inner_radius(z):
    return get_outer_radius(z) - (t_pile if z <= 0.0 else t_tower)

p = my_model.parts[my_model.parts.keys()[0]]
print("imported: %d nodes, %d elements" % (len(p.nodes), len(p.elements)))

# gmsh also writes the 1D/2D entities (T3D2 truss, CPS3 shell). Only the C3D4
# tets are wanted -- the rest would add spurious stiffness and mass.
junk = [e for e in p.elements if e.type != C3D4]
if len(junk) > 0:
    p.SetFromElementLabels(name='_junk', elementLabels=tuple([e.label for e in junk]))
    p.deleteElement(elements=p.sets['_junk'], deleteUnreferencedNodes=ON)
    print("removed %d non-solid elements" % len(junk))
print("after cleanup: %d nodes, %d elements" % (len(p.nodes), len(p.elements)))
p.Set(elements=p.elements, name='Set-All')

# -------------------------------------------------------------------
# 5. SURFACES  (rebuilt from element faces -- an orphan mesh has no faces
#               for findAt to pick)
# -------------------------------------------------------------------
# A tet face shared by two elements is interior; a face owned by exactly one
# element is on the boundary. Abaqus C3D4 face numbering:
FACE_NODES = {1: (0, 1, 2), 2: (0, 3, 1), 3: (1, 3, 2), 4: (0, 2, 3)}

coord = {}
for n in p.nodes:
    coord[n.label] = n.coordinates

face_count = {}
face_owner = {}
for el in p.elements:
    labels = [p.nodes[c].label for c in el.connectivity]
    for fid, idx in FACE_NODES.items():
        key = tuple(sorted([labels[i] for i in idx]))
        face_count[key] = face_count.get(key, 0) + 1
        face_owner[key] = (el.label, fid)
exterior = [(k, face_owner[k]) for k, c in face_count.items() if c == 1]
print("exterior faces: %d" % len(exterior))

def face_centroid(key):
    pts = [coord[l] for l in key]
    return (sum([q[0] for q in pts]) / 3.0,
            sum([q[1] for q in pts]) / 3.0,
            sum([q[2] for q in pts]) / 3.0)

def face_area(key):
    pts = [coord[l] for l in key]
    a1 = [pts[1][i] - pts[0][i] for i in range(3)]
    b1 = [pts[2][i] - pts[0][i] for i in range(3)]
    cr = (a1[1]*b1[2] - a1[2]*b1[1], a1[2]*b1[0] - a1[0]*b1[2], a1[0]*b1[1] - a1[1]*b1[0])
    return 0.5 * math.sqrt(sum([q*q for q in cr]))

def on_outer_wall(key):
    """Per-NODE vote against the wall at that node's own z.

    Not a centroid test: on a curved, tapered wall the face centroid sits
    inboard of BOTH surfaces by the chord sagitta, which here exceeds half the
    0.060 wall, so a centroid test rejects nearly every genuine outer face.
    """
    votes = 0
    for l in key:
        x, y, z = coord[l]
        r = math.hypot(x, y)
        if abs(r - get_outer_radius(z)) <= abs(r - get_inner_radius(z)):
            votes += 1
    return votes >= 2

def face_is_flat_at(key, z_plane):
    """True only if EVERY node of the face lies on the z_plane end cap.

    A centroid-distance test does not work. With a 0.05 m tolerance and coarse
    elements nothing but the cap is close enough, but once the cap is refined
    the wall faces just inside the end plane also fall within the tolerance and
    get swallowed into Surf-1/Surf-2 -- which inflated those annuli by 57% and,
    since they carry the RNA tractions, applied 14-17% too much total load.
    """
    for l in key:
        if abs(coord[l][2] - z_plane) > 1.0e-6:
            return False
    return True

# Bucket every exterior face into the surface it belongs to, using the same
# geometric criteria the original findAt pick points encoded: the loaded
# surfaces are on the OUTER wall, on the -X (load-facing) half, in their z band.
buckets = {}
for key, owner in exterior:
    cx, cy, cz = face_centroid(key)
    if face_is_flat_at(key, z_mudline):
        buckets.setdefault('Surf-1', []).append(owner)      # mudline annulus
        continue
    if face_is_flat_at(key, z_hub):
        buckets.setdefault('Surf-2', []).append(owner)      # hub annulus
        continue
    if not on_outer_wall(key) or cx >= 0.0:
        continue                                            # inner wall / leeward
    if cz <= z_mwl:
        buckets.setdefault('Surf-3', []).append(owner)      # windward pile
        continue
    for i in range(len(tower_levels) - 1):
        if tower_levels[i] <= cz < tower_levels[i + 1]:
            buckets.setdefault('Surf-%d' % (4 + i), []).append(owner)
            break

for name in sorted(buckets, key=lambda s: int(s.split('-')[1])):
    by_face = {}
    for elabel, fid in buckets[name]:
        by_face.setdefault(fid, []).append(elabel)
    kwargs = {'name': name}
    for fid, labels in by_face.items():
        kwargs['face%dElements' % fid] = p.elements.sequenceFromLabels(labels)
    p.Surface(**kwargs)

print("surfaces built: %d -> %s" % (len(p.surfaces),
      sorted(p.surfaces.keys(), key=lambda s: int(s.split('-')[1]))))

mud_nodes = p.nodes.getByBoundingBox(zMin=z_mudline - 0.01, zMax=z_mudline + 0.01)
p.Set(nodes=mud_nodes, name='Set-Mudline')

# --- ASSERT the loaded surfaces cover their bands, before any load is applied
def analytic_band_area(z0, z1):
    r0, r1 = get_outer_radius(z0), get_outer_radius(z1)
    sl = math.sqrt((r1 - r0)**2 + (z1 - z0)**2)
    return 0.5 * math.pi * (r0 + r1) * sl      # the -X half of the frustum

area_of = {}
owner_to_surf = {}
for name, owners in buckets.items():
    for o in owners:
        owner_to_surf[o] = name
for key, owner in exterior:
    nm = owner_to_surf.get(owner)
    if nm:
        area_of[nm] = area_of.get(nm, 0.0) + face_area(key)

expect = {'Surf-3': 0.5 * math.pi * D_pile * abs(z_mwl - z_mudline)}
for i in range(len(tower_levels) - 1):
    expect['Surf-%d' % (4 + i)] = analytic_band_area(tower_levels[i], tower_levels[i+1])
# Surf-1 / Surf-2 MUST be checked too: they carry the RNA weight and thrust, by
# far the largest tractions here, so an over-sized annulus silently rescales the
# whole load case.
expect['Surf-1'] = math.pi * ((D_pile/2.0)**2 - (D_pile/2.0 - t_pile)**2)
expect['Surf-2'] = math.pi * ((D_top/2.0)**2 - (D_top/2.0 - t_tower)**2)

print("--- surface coverage (actual vs analytic) ---")
bad = 0
for nm in sorted(expect, key=lambda s: int(s.split('-')[1])):
    act = area_of.get(nm, 0.0)
    exp = expect[nm]
    ratio = act / exp if exp else 0.0
    # 5% band: flat triangles chord a curved surface, so the meshed area lands a
    # shade under analytic. A missing or doubled band shows up near 0 or 2.
    ok = abs(ratio - 1.0) < 0.05
    if not ok: bad += 1
    print("  %-8s actual=%9.3f expected=%9.3f ratio=%.3f %s"
          % (nm, act, exp, ratio, "OK" if ok else "<<< MISMATCH"))
if bad:
    raise ValueError("%d loaded surface(s) mis-built -- loads would be wrong." % bad)
print("all loaded surfaces verified")

# -------------------------------------------------------------------
# 7. MATERIAL & SECTION
# -------------------------------------------------------------------
my_model.HomogeneousSolidSection(name='SteelSection', material='Steel', thickness=None)
mat = my_model.Material(name='Steel')
mat.Density(table=((rho_steel, ), ))
mat.Elastic(table=((E_steel, nu_steel), ))
p.SectionAssignment(region=p.sets['Set-All'], sectionName='SteelSection',
                    thicknessAssignment=FROM_SECTION)

# -------------------------------------------------------------------
# 8. ASSEMBLY & LOADS
# -------------------------------------------------------------------
a = my_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)

# ModelFromInputFile already instances the imported part. Adding a second
# instance leaves the FIRST one in the assembly unrestrained -- the BC below
# pins only OWT-1, so the duplicate body floats free. That shows up as six
# rigid body modes and the job dies with "numerical singularity" /
# "too many attempts".
for _nm in list(a.instances.keys()):
    del a.instances[_nm]
inst = a.Instance(name='OWT-1', part=p, dependent=ON)

my_model.StaticStep(name='Step-Load', previous='Initial', nlgeom=ON)
my_model.Gravity(name='Load-SelfWeight', createStepName='Step-Load', comp3=-g, distributionType=UNIFORM, field='')

region_bot = a.instances['OWT-1'].sets['Set-Mudline']
my_model.DisplacementBC(name='BC-Pinned', createStepName='Initial', region=region_bot, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

surf_2 = a.instances['OWT-1'].surfaces['Surf-2']
my_model.SurfaceTraction(name='Load-RNA-Weight', createStepName='Step-Load', region=surf_2, distributionType=UNIFORM, directionVector=((0,0,0),(0,0,1)), magnitude=-rna_tr_weight, traction=GENERAL)
my_model.SurfaceTraction(name='Load-RNA-Thrust', createStepName='Step-Load', region=surf_2, distributionType=UNIFORM, directionVector=((0,0,0),(1,0,0)), magnitude=rna_tr_thrust, traction=GENERAL)

surf_3 = a.instances['OWT-1'].surfaces['Surf-3']
my_model.SurfaceTraction(name='Load-Wave', createStepName='Step-Load', region=surf_3, distributionType=UNIFORM, directionVector=((0,0,0),(1,0,0)), magnitude=wave_tr_val, traction=GENERAL)

surf_id = 4
for i in range(len(tower_levels)-1):
    z_t = tower_levels[i+1]
    s_name = 'Surf-%d' % surf_id
    if s_name in a.instances['OWT-1'].surfaces.keys():
        seg_tr_val = calc_tower_traction_stress(z_t)
        surf_inst = a.instances['OWT-1'].surfaces[s_name]
        l_name = 'Load-Wind-Seg-%d' % surf_id
        my_model.SurfaceTraction(name=l_name, createStepName='Step-Load', region=surf_inst, distributionType=UNIFORM, directionVector=((0,0,0),(1,0,0)), magnitude=seg_tr_val, traction=GENERAL)
    surf_id += 1

print("Model Loaded: unstructured C3D4 mesh imported from gmsh.")