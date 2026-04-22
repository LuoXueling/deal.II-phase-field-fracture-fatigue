# -*- coding: mbcs -*-
from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import assembly
import math

# -------------------------------------------------------------------
# 1. INITIALIZATION
# -------------------------------------------------------------------
model_name = 'OWT_Final_Fixed_Bias_InnerOuter'
if model_name in mdb.models:
    del mdb.models[model_name]
my_model = mdb.Model(name=model_name)

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

# Meshing
max_size = 0.5
mid_size = 0.1
min_size = 0.01

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
# 4. GEOMETRY
# -------------------------------------------------------------------
p = my_model.Part(name='OWT_NoEmbed', dimensionality=THREE_D, type=DEFORMABLE_BODY)
d = p.datums
plane_yz = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
axis_z = p.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
t = p.MakeSketchTransform(sketchPlane=d[plane_yz.id], sketchUpEdge=d[axis_z.id], 
                          sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s = my_model.ConstrainedSketch(name='profile', sheetSize=200.0, transform=t)
s.ConstructionLine(point1=(0.0, -200.0), point2=(0.0, 200.0))

r_p_out    = D_pile / 2.0
r_p_in     = r_p_out - t_pile
r_tw_b_in  = (D_base / 2.0) - t_tower
r_tw_t_out = D_top / 2.0
r_tw_t_in  = r_tw_t_out - t_tower

s.Line(point1=(r_p_in, z_mudline), point2=(r_p_in, z_mwl))       
s.Line(point1=(r_tw_b_in, z_mwl), point2=(r_tw_t_in, z_hub))     
s.Line(point1=(r_tw_t_in, z_hub), point2=(r_tw_t_out, z_hub))    
s.Line(point1=(r_tw_t_out, z_hub), point2=(r_p_out, z_mwl))      
s.Line(point1=(r_p_out, z_mwl), point2=(r_p_out, z_mudline))     
s.Line(point1=(r_p_out, z_mudline), point2=(r_p_in, z_mudline))  

p.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
del s

# -------------------------------------------------------------------
# 5. PARTITIONING
# -------------------------------------------------------------------
# Structural Levels
z_struct = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
# Mesh Refinement Partitions (2m from Top and Bottom)
z_mesh = [-28.0, 85.0]

all_z_cuts = sorted(z_struct + z_mesh)

for z_lev in all_z_cuts:
    pt = p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=z_lev)
    p.PartitionCellByDatumPlane(datumPlane=p.datums[pt.id], cells=p.cells)

# Circumferential Partition (X=0)
plane_x0 = p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
p.PartitionCellByDatumPlane(datumPlane=p.datums[plane_x0.id], cells=p.cells)

# -------------------------------------------------------------------
# 6. SURFACES
# -------------------------------------------------------------------
# Surf-1: Bottom Face
faces_bot = p.faces.getByBoundingBox(zMin=z_mudline-0.1, zMax=z_mudline+0.1)
p.Surface(side1Faces=faces_bot, name='Surf-1')
p.Set(faces=faces_bot, name='Set-Mudline')

# Surf-2: Top Face
faces_top = p.faces.getByBoundingBox(zMin=z_hub-0.1, zMax=z_hub+0.1)
p.Surface(side1Faces=faces_top, name='Surf-2')

# Surf-3: Windward Pile (Split at -28.0)
pt_pile_bot = (-r_p_out, 0.0, -29.0)
pt_pile_main = (-r_p_out, 0.0, -14.0)
f_p_w = p.faces.findAt( (pt_pile_bot,), (pt_pile_main,) )
p.Surface(side1Faces=f_p_w, name='Surf-3')

# Surf-4 to Surf-12: Windward Tower Segments
tower_levels = z_struct + [z_hub] 
surf_counter = 4

for i in range(len(tower_levels)-1):
    z_b = tower_levels[i]
    z_t = tower_levels[i+1]
    
    pick_points = []
    
    # Check if split at 85.0
    if abs(z_b - 80.0) < 0.1 and abs(z_t - 87.0) < 0.1:
        # Split: 80-85 and 85-87
        z_pick1 = 82.5; r_pick1 = get_outer_radius(z_pick1)
        pick_points.append( (-r_pick1, 0.0, z_pick1) )
        z_pick2 = 86.0; r_pick2 = get_outer_radius(z_pick2)
        pick_points.append( (-r_pick2, 0.0, z_pick2) )
    else:
        # Standard
        z_pick = (z_b + z_t) / 2.0; r_pick = get_outer_radius(z_pick)
        pick_points.append( (-r_pick, 0.0, z_pick) )
        
    args = []
    for pt in pick_points: args.append((pt,))
    f_seg = p.faces.findAt(*args)
    p.Surface(side1Faces=f_seg, name='Surf-%d' % surf_counter)
    surf_counter += 1

# -------------------------------------------------------------------
# 7. MESHING & SEEDING (FIXED BIAS KEYWORDS + INNER/OUTER ARCS)
# -------------------------------------------------------------------
my_model.HomogeneousSolidSection(name='SteelSection', material='Steel', thickness=None)
mat = my_model.Material(name='Steel')
mat.Density(table=((rho_steel, ), ))
mat.Elastic(table=((E_steel, nu_steel), ))

p.Set(cells=p.cells, name='Set-All')
p.SectionAssignment(region=p.sets['Set-All'], sectionName='SteelSection', thicknessAssignment=FROM_SECTION)

p.setMeshControls(regions=p.cells, elemShape=HEX, technique=STRUCTURED)
elemType = mesh.ElemType(elemCode=C3D8, elemLibrary=STANDARD)
p.setElementType(regions=p.sets['Set-All'], elemTypes=(elemType,))

# --- SEEDING ---

# 1. Thickness Seeding (0.02 Uniform)
all_cuts_mesh = [z_mudline, -28.0] + z_struct + [85.0, z_hub]
for z in all_cuts_mesh:
    # Select horizontal edges on X=0 plane
    thick_edges = p.edges.getByBoundingBox(xMin=-0.01, xMax=0.01, zMin=z-0.01, zMax=z+0.01)
    real_thick = [e for e in thick_edges if e.getSize() < 0.2]
    if len(real_thick) > 0:
        # Positional arg for edges
        p.seedEdgeBySize(real_thick, size=min_size, constraint=FIXED)

# 2. Axial Bias (Top/Bottom 2m)
# Fix: Using 'end2Edges' to flip direction if 'end1Edges' was backwards

def split_inner_outer(edges, z_mid):
    inner = []
    outer = []
    r_outer_expected = get_outer_radius(z_mid)
    for e in edges:
        pt = e.pointOn[0]
        r_actual = abs(pt[1])
        if r_actual > (r_outer_expected - 0.03): outer.append(e)
        else: inner.append(e)
    return inner, outer

# Bottom 2m (-30 to -28)
edges_bot_2m = p.edges.getByBoundingBox(xMin=-0.01, xMax=0.01, zMin=-30.1, zMax=-27.9)
real_bot_2m = [e for e in edges_bot_2m if e.getSize() > 1.0]

if len(real_bot_2m) > 0:
    inner_bot, outer_bot = split_inner_outer(real_bot_2m, -29.0)
    
    # Inner: Use end1Edges
    if len(inner_bot) > 0:
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=inner_bot, 
                         minSize=min_size, maxSize=max_size, constraint=FIXED)
                         
    # Outer: Use end2Edges (Flipped vs Inner)
    if len(outer_bot) > 0:
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=outer_bot, 
                         minSize=min_size, maxSize=max_size, constraint=FIXED)

# Top 2m (85 to 87)
edges_top_2m = p.edges.getByBoundingBox(xMin=-0.01, xMax=0.01, zMin=84.9, zMax=87.1)
real_top_2m = [e for e in edges_top_2m if e.getSize() > 1.0]

if len(real_top_2m) > 0:
    inner_top, outer_top = split_inner_outer(real_top_2m, 86.0)
    
    # Inner: Use end2Edges
    if len(inner_top) > 0:
        p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=inner_top, 
                         minSize=min_size, maxSize=max_size, constraint=FIXED)
                         
    # Outer: Use end1Edges (Flipped vs Inner)
    if len(outer_top) > 0:
        p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=outer_top, 
                         minSize=min_size, maxSize=max_size, constraint=FIXED)

# 3. Circumferential Double Bias
# UPDATED: Select both Outer AND Inner arcs

# Bottom Arc (-30)
r_bot = D_pile / 2.0
r_bot_in = r_bot - t_pile # Inner radius
# Select outer (-r_bot) and inner (-r_bot_in) points
arcs_bot = p.edges.findAt( ((-r_bot, 0.0, z_mudline), ), 
                           ((-r_bot_in, 0.0, z_mudline), ) )

if len(arcs_bot) > 0:
    p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=arcs_bot, 
                     minSize=min_size, maxSize=mid_size, constraint=FINER)

# # Top Arc (87)
# r_top_exact = get_outer_radius(z_hub)
# r_top_in = r_top_exact - t_tower # Inner radius
# # Select outer and inner points
# arcs_top = p.edges.findAt( ((-r_top_exact, 0.0, z_hub), ), 
#                            ((-r_top_in, 0.0, z_hub), ) )

# if len(arcs_top) > 0:
#     p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=arcs_top, 
#                      minSize=min_size, maxSize=mid_size, constraint=FINER)

# Global Seed
p.seedPart(size=max_size, deviationFactor=0.1, minSizeFactor=0.1)

p.generateMesh()

# -------------------------------------------------------------------
# 8. ASSEMBLY & LOADS
# -------------------------------------------------------------------
a = my_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
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

print("Model Loaded: Fixed Bias Direction & Added Inner Arcs to Selection.")