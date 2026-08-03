from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from abaqusConstants import *
from regionToolset import Region
import math
import csv
import os, sys

# ============================================================
# USER INPUT
# ============================================================
# load case name
loadCaseName = 'TensionBiaxial'       # 'TensionUniaxial', 'TensionBiaxial', 'CompressiveUniaxial', 'CompressiveBiaxial', 'SimpleShear', 'PureShear', 'PureBending', 'ParabolicBending'


MODEL = 'Model-1'
PART  = 'RVE2D_FULL'

# ---- choose one translation zeta for this run -------------
zeta1 = 0.15   # shift in x
zeta2 = -0.15   # shift in y

# output folder for exported nodal csv
OUTDIR = r'OUTPUT_CSV'
ENERGYDIR = r'OUTPUT_ENERGY'

# Unit cell
Lx = 1.0
Ly = 1.0
R  = 0.30 * Lx
mesh_size = 0.05 * Lx

# Repeats
NX = 4
NY = 4

materialID = 'NEl1'
Rho_Inc = 1.12e-06
Rho_Matrix = 1.12e-06

C10_Inc = 400/100.0
C10_Matrix = 400.0

D1_Inc = 0.0011429*100.0
D1_Matrix = 0.0011429

E_m  = 70000.0
nu_m = 0.33
E_i  = 3500.0
nu_i = 0.33
rho_m = 1.12e-06
rho_i = rho_m

STATE_2D = 'plane_strain'
TOL = 1e-6

elongation = 0.5
gamma = 0.3 * elongation 
kappa = 0.15 * elongation
kappa_parabolic = 5 * kappa
# job name includes load case + zeta
JOB = '%s_zx_%0.3f_zy_%0.3f' % (loadCaseName, zeta1, zeta2)
JOB = JOB.replace('.', 'p').replace('-', 'm')

# ============================================================
# MODEL
# ============================================================
if MODEL not in mdb.models:
    mdb.Model(name=MODEL)

model = mdb.models[MODEL]

if JOB in mdb.jobs:
    del mdb.jobs[JOB]

Lx_tot = NX * Lx
Ly_tot = NY * Ly

# ============================================================
# GEOMETRY: build whole domain directly
# ============================================================
if PART in model.parts:
    del model.parts[PART]

sk = model.ConstrainedSketch(name='__profile__', sheetSize=10.0)
sk.rectangle(point1=(0.0, 0.0), point2=(Lx_tot, Ly_tot))

p = model.Part(name=PART, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p.BaseShell(sketch=sk)
del model.sketches['__profile__']

f0 = p.faces[0]
tr = p.MakeSketchTransform(sketchPlane=f0, sketchPlaneSide=SIDE1,
                           sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
sk2 = model.ConstrainedSketch(name='__profile__', sheetSize=10.0, transform=tr)
p.projectReferencesOntoSketch(sketch=sk2, filter=COPLANAR_EDGES)

# ------------------------------------------------------------
# Draw all shifted inclusions, including periodic-image circles
# so circles crossing specimen boundaries are represented.
# ------------------------------------------------------------
def circle_intersects_domain(cx, cy, R, xmin, xmax, ymin, ymax):
    return not (cx + R < xmin or cx - R > xmax or cy + R < ymin or cy - R > ymax)

# loop over real cells plus one image layer around them
for i in range(-1, NX + 1):
    for j in range(-1, NY + 1):
        cx = (i + 0.5) * Lx + zeta1
        cy = (j + 0.5) * Ly + zeta2
        if circle_intersects_domain(cx, cy, R, 0.0, Lx_tot, 0.0, Ly_tot):
            sk2.CircleByCenterPerimeter(center=(cx, cy), point1=(cx + R, cy))

p.PartitionFaceBySketch(faces=(f0,), sketch=sk2)
del model.sketches['__profile__']

# ============================================================
# MATERIALS + SECTIONS
# ============================================================
if 'Matrix' in model.materials.keys():
    del model.materials['Matrix']

if 'Inclusion' in model.materials.keys():
    del model.materials['Inclusion']

model.Material(name='Matrix')
model.Material(name='Inclusion')

if materialID in ('HEl1', 'NEl1'):
    model.materials['Matrix'].Density(table=((rho_m,),))
    model.materials['Matrix'].Elastic(table=((E_m, nu_m),))
    model.materials['Inclusion'].Density(table=((rho_i,),))
    model.materials['Inclusion'].Elastic(table=((E_i, nu_i),))
elif materialID in ('HNe1', 'NNe1'):
    model.materials['Matrix'].Density(table=((Rho_Matrix,),))
    model.materials['Matrix'].Hyperelastic(
        materialType=ISOTROPIC,
        type=NEO_HOOKE,
        table=((C10_Matrix, D1_Matrix),),
        testData=OFF,
        volumetricResponse=VOLUMETRIC_DATA
    )
    model.materials['Inclusion'].Density(table=((Rho_Inc,),))
    model.materials['Inclusion'].Hyperelastic(
        materialType=ISOTROPIC,
        type=NEO_HOOKE,
        table=((C10_Inc, D1_Inc),),
        testData=OFF,
        volumetricResponse=VOLUMETRIC_DATA
    )
elif materialID in ('HMo1', 'NMo1'):
    model.materials['Matrix'].Density(table=((Rho_Matrix,),))
    model.materials['Matrix'].Hyperelastic(
        materialType=ISOTROPIC,
        type=MOONEY_RIVLIN,
        table=((C10_Matrix, C01_Matrix, D1_Matrix),),
        testData=OFF,
        volumetricResponse=VOLUMETRIC_DATA
    )
    model.materials['Inclusion'].Density(table=((Rho_Inc,),))
    model.materials['Inclusion'].Hyperelastic(
        materialType=ISOTROPIC,
        type=MOONEY_RIVLIN,
        table=((C10_Inc, C01_Inc, D1_Inc),),
        testData=OFF,
        volumetricResponse=VOLUMETRIC_DATA
    )

if 'Sec-Matrix' in model.sections.keys():
    del model.sections['Sec-Matrix']

if 'Sec-Inclusion' in model.sections.keys():
    del model.sections['Sec-Inclusion']

model.HomogeneousSolidSection(name='Sec-Matrix', material='Matrix', thickness=1.0)
model.HomogeneousSolidSection(name='Sec-Inclusion', material='Inclusion', thickness=1.0)

# default all faces = matrix
p.Set(name='SET_ALL', faces=p.faces[:])
p.SectionAssignment(region=p.sets['SET_ALL'], sectionName='Sec-Matrix')

# faces whose centroid lies inside any shifted circle -> inclusion
# default all faces = matrix
if 'SET_ALL' in p.sets.keys():
    del p.sets['SET_ALL']

p.Set(name='SET_ALL', faces=p.faces[:])
p.SectionAssignment(region=p.sets['SET_ALL'], sectionName='Sec-Matrix')

# ------------------------------------------------------------
# inclusion faces: classify every face by its centroid
# and build the set using FaceArray indexing
# ------------------------------------------------------------
inc_idx = []

for k, face in enumerate(p.faces):
    c = face.getCentroid()[0]
    x = c[0]
    y = c[1]
    inside_any = False
    for i in range(-1, NX + 1):
        for j in range(-1, NY + 1):
            cx = (i + 0.5) * Lx + zeta1
            cy = (j + 0.5) * Ly + zeta2
            if (x - cx)**2 + (y - cy)**2 <= R**2:
                inside_any = True
                break
        if inside_any:
            break
    if inside_any:
        inc_idx.append(k)

print('Number of faces in part =', len(p.faces))
print('Number of inclusion faces found =', len(inc_idx))

if len(inc_idx) == 0:
    raise RuntimeError('No inclusion faces found. Check zeta and geometry.')

# build FaceArray by concatenating slices from p.faces
inc_faces = p.faces[inc_idx[0]:inc_idx[0]+1]
for k in inc_idx[1:]:
    inc_faces = inc_faces + p.faces[k:k+1]

if 'SET_INC' in p.sets.keys():
    del p.sets['SET_INC']

p.Set(name='SET_INC', faces=inc_faces)
p.SectionAssignment(region=p.sets['SET_INC'], sectionName='Sec-Inclusion')

# ============================================================
# ASSEMBLY + STEP
# ============================================================
a = model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)

for nm in list(a.instances.keys()):
    del a.instances[nm]

INST = PART + '-1'
inst = a.Instance(name=INST, part=p, dependent=ON)

STEP = 'Step-1'

if materialID in ('HEl1', 'NEl1'):
    model.StaticStep(name=STEP, previous='Initial', nlgeom=OFF,
                     maxNumInc=2, initialInc=1.0, minInc=1.0, maxInc=1.0)
else:
    model.StaticStep(name=STEP, previous='Initial', nlgeom=ON,
                     maxNumInc=10, initialInc=0.1, minInc=0.1, maxInc=0.1)

# ============================================================
# MESH
# ============================================================
p.seedPart(size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
p.setMeshControls(regions=p.faces[:], technique=FREE, elemShape=TRI)

if STATE_2D == 'plane_strain':
    p.setElementType(regions=(p.faces[:],),
                     elemTypes=(ElemType(elemCode=CPE6, elemLibrary=STANDARD),))
else:
    p.setElementType(regions=(p.faces[:],),
                     elemTypes=(ElemType(elemCode=CPS6, elemLibrary=STANDARD),))

p.generateMesh()

# regenerate assembly instance after meshing
a.regenerate()
inst = a.instances[INST]

# ============================================================
# BCs
# ============================================================

if loadCaseName == 'TensionUniaxial':
    top_u2 = elongation*Ly_tot
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    corner00 = None
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
        if abs(x - 0.0) < TOL and abs(y - 0.0) < TOL:
            corner00 = nd.label
    if len(top_labs) == 0 or len(bottom_labs) == 0:
        raise RuntimeError("TOP/BOTTOM empty. Increase TOL.")
    if corner00 is None:
        raise RuntimeError("Corner (0,0) not found. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    a.Set(name='CORNER00', nodes=inst.nodes.sequenceFromLabels((corner00,)))
    model.DisplacementBC(
        name='BC_BOTTOM_U2',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=UNSET, u2=0.0, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_TOP_U2',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=UNSET, u2=top_u2, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_FIX_CORNER_U1',
        createStepName='Initial',
        region=a.sets['CORNER00'],
        u1=0.0, u2=UNSET, ur3=UNSET
    )
elif loadCaseName == 'TensionBiaxial':
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
    print("Boundary counts:",
          "LEFT", len(left_labs),
          "RIGHT", len(right_labs),
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(left_labs) == 0 or len(right_labs) == 0:
        raise RuntimeError("LEFT/RIGHT empty. Increase TOL.")
    if len(bottom_labs) == 0 or len(top_labs) == 0:
        raise RuntimeError("BOTTOM/TOP empty. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    Ux = 0.5 *elongation*Lx_tot
    Uy = 0.5 *elongation*Ly_tot
    model.DisplacementBC(
        name='BC_LEFT_U1',
        createStepName=STEP,
        region=a.sets['LEFT'],
        u1=0.0, u2=UNSET, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_RIGHT_U1',
        createStepName=STEP,
        region=a.sets['RIGHT'],
        u1=Ux, u2=UNSET, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_BOTTOM_U2',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=UNSET, u2=0.0, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_TOP_U2',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=UNSET, u2=Uy, ur3=UNSET
    )
elif loadCaseName == 'CompressiveUniaxial':
    top_labs, bottom_labs = [], []
    corner00 = None
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
        if abs(x - 0.0) < TOL and abs(y - 0.0) < TOL:
            corner00 = nd.label
    print("Boundary counts:",
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(top_labs) == 0 or len(bottom_labs) == 0:
        raise RuntimeError("TOP/BOTTOM empty. Increase TOL.")
    if corner00 is None:
        raise RuntimeError("Corner (0,0) not found. Increase TOL.")
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    a.Set(name='CORNER00', nodes=inst.nodes.sequenceFromLabels((corner00,)))
    Uy = 0.3*elongation*Lx_tot   # negative = compression downward
    model.DisplacementBC(
        name='BC_BOTTOM_U2',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=UNSET, u2=0.0, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_TOP_U2',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=UNSET, u2=-Uy, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_FIX_CORNER_U1',
        createStepName='Initial',
        region=a.sets['CORNER00'],
        u1=0.0, u2=UNSET, ur3=UNSET
    )
elif loadCaseName == 'CompressiveBiaxial':
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
    print("Boundary counts:",
          "LEFT", len(left_labs),
          "RIGHT", len(right_labs),
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(left_labs) == 0 or len(right_labs) == 0:
        raise RuntimeError("LEFT/RIGHT empty. Increase TOL.")
    if len(bottom_labs) == 0 or len(top_labs) == 0:
        raise RuntimeError("BOTTOM/TOP empty. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    Ux = 0.25 * elongation*Lx_tot
    Uy = 0.25 * elongation*Ly_tot
    model.DisplacementBC(
        name='BC_LEFT_U1',
        createStepName=STEP,
        region=a.sets['LEFT'],
        u1=0.0, u2=UNSET, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_RIGHT_U1',
        createStepName=STEP,
        region=a.sets['RIGHT'],
        u1=-Ux, u2=UNSET, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_BOTTOM_U2',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=UNSET, u2=0.0, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_TOP_U2',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=UNSET, u2=-Uy, ur3=UNSET
    )
elif loadCaseName == 'SimpleShear':
    bottom_labs, top_labs = [], []
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
    print("Boundary counts:",
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(bottom_labs) == 0 or len(top_labs) == 0:
        raise RuntimeError("BOTTOM/TOP empty. Increase TOL.")
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    Ux = gamma * Ly_tot
    model.DisplacementBC(
        name='BC_BOTTOM_FIX',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=0.0, u2=0.0, ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_TOP_SHEAR',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=Ux, u2=UNSET, ur3=UNSET
    )
elif loadCaseName == 'PureShear':
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
    print("Boundary counts:",
          "LEFT", len(left_labs),
          "RIGHT", len(right_labs),
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(left_labs) == 0 or len(right_labs) == 0:
        raise RuntimeError("LEFT/RIGHT empty. Increase TOL.")
    if len(bottom_labs) == 0 or len(top_labs) == 0:
        raise RuntimeError("BOTTOM/TOP empty. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    Uy = gamma * Ly_tot
    Ux = gamma * Lx_tot
    model.DisplacementBC(
        name='BC_LEFTu2',
        createStepName=STEP,
        region=a.sets['LEFT'],
        u1=UNSET,
        u2=0.0,
        ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_Rightu2',
        createStepName=STEP,
        region=a.sets['RIGHT'],
        u1=UNSET,
        u2=Uy,
        ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_Bottomu1',
        createStepName=STEP,
        region=a.sets['BOTTOM'],
        u1=0.0,
        u2=UNSET,
        ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_Topu1',
        createStepName=STEP,
        region=a.sets['TOP'],
        u1=Ux ,
        u2=UNSET,
        ur3=UNSET
    )
elif loadCaseName == 'PureBending':
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    corner00 = None
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
        if abs(x - 0.0) < TOL and abs(y - 0.0) < TOL:
            corner00 = nd.label
    if len(left_labs) == 0 or len(right_labs) == 0:
        raise RuntimeError("LEFT/RIGHT empty. Increase TOL.")
    if corner00 is None:
        raise RuntimeError("Corner (0,0) not found. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    a.Set(name='CORNER00', nodes=inst.nodes.sequenceFromLabels((corner00,)))
    # curvature-like parameter
    # normalized Y coordinate in [-1, +1]
    bendFieldName = 'PUREBEND_Y_NORM'
    model.ExpressionField(
        name=bendFieldName,
        expression='((Y - %g)/%g)' % (0.5*Ly_tot, 0.5*Ly_tot),
        localCsys=None,
        description='Linear profile across Y for pure bending'
    )
    # amplitude for horizontal displacement
    # u1 ~ +/- const * (Y - Ly/2)
    U1_bend_amp =  kappa * Lx_tot 
    model.DisplacementBC(
        name='BC_LEFT_BEND_U1',
        createStepName=STEP,
        region=a.sets['LEFT'],
        distributionType=FIELD,
        fieldName=bendFieldName,
        u1=-U1_bend_amp,
        u2=UNSET,
        ur3=UNSET
    )
    model.DisplacementBC(
        name='BC_RIGHT_BEND_U1',
        createStepName=STEP,
        region=a.sets['RIGHT'],
        distributionType=FIELD,
        fieldName=bendFieldName,
        u1=+U1_bend_amp,
        u2=UNSET,
        ur3=UNSET
    )
    # minimal constraint to remove rigid vertical motion
    model.DisplacementBC(
        name='BC_FIX_CORNER_U2',
        createStepName='Initial',
        region=a.sets['CORNER00'],
        u1=UNSET,
        u2=0.0,
        ur3=UNSET
    )
elif loadCaseName == 'ParabolicBending':
    bend_amp = 0.40
    left_labs, right_labs, bottom_labs, top_labs = [], [], [], []
    corner00 = None
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - 0.0) < TOL:
            left_labs.append(nd.label)
        if abs(x - Lx_tot) < TOL:
            right_labs.append(nd.label)
        if abs(y - 0.0) < TOL:
            bottom_labs.append(nd.label)
        if abs(y - Ly_tot) < TOL:
            top_labs.append(nd.label)
        if abs(x - 0.0) < TOL and abs(y - 0.0) < TOL:
            corner00 = nd.label
    print("Boundary counts:",
          "LEFT", len(left_labs),
          "RIGHT", len(right_labs),
          "BOTTOM", len(bottom_labs),
          "TOP", len(top_labs))
    if len(left_labs) == 0 or len(right_labs) == 0:
        raise RuntimeError("LEFT/RIGHT empty. Increase TOL.")
    if corner00 is None:
        raise RuntimeError("Corner (0,0) not found. Increase TOL.")
    a.Set(name='LEFT',   nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
    a.Set(name='RIGHT',  nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
    a.Set(name='BOTTOM', nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
    a.Set(name='TOP',    nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
    a.Set(name='CORNER00', nodes=inst.nodes.sequenceFromLabels((corner00,)))
    # Parabolic shape in y:
    # f(y) = 4*(y/L)*(1-y/L)
    # max at mid-height, zero at top and bottom
    expr_parab = '(Y/L)*(1.0-Y/L)'
    bend_amp = kappa_parabolic * Lx_tot
    model.ExpressionField(
        name='PARAB_Y_FIELD',
        localCsys=None,
        expression=expr_parab.replace('L', str(Ly_tot))
    )
    model.DisplacementBC(
        name='BC_LEFT_PARABOLIC',
        createStepName=STEP,
        region=a.sets['LEFT'],
        u1=-bend_amp,
        u2=UNSET,
        ur3=UNSET,
        distributionType=FIELD,
        fieldName='PARAB_Y_FIELD'
    )
    model.DisplacementBC(
        name='BC_RIGHT_PARABOLIC',
        createStepName=STEP,
        region=a.sets['RIGHT'],
        u1=+bend_amp,
        u2=UNSET,
        ur3=UNSET,
        distributionType=FIELD,
        fieldName='PARAB_Y_FIELD'
    )
    model.DisplacementBC(
        name='BC_FIX_CORNER_U2',
        createStepName='Initial',
        region=a.sets['CORNER00'],
        u1=UNSET,
        u2=0.0,
        ur3=UNSET
    )           
else:
    raise RuntimeError('Unknown loadCaseName: %s' % loadCaseName)

model.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'LE', 'U', 'RF', 'COORD')
)

# ============================================================
# JOB
# ============================================================
mdb.Job(
    name=JOB, model=MODEL, type=ANALYSIS,
    memory=90, memoryUnits=PERCENTAGE,
    numCpus=1, numGPUs=0,
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE
)

mdb.jobs[JOB].submit(consistencyChecking=OFF)
mdb.jobs[JOB].waitForCompletion()



# ============================================================
# EXPORT ALL FRAMES: NODAL COORDINATES + DISPLACEMENTS
# ============================================================
from odbAccess import openOdb

if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)

caseOutDir = os.path.join(OUTDIR, loadCaseName)
if not os.path.isdir(caseOutDir):
    os.makedirs(caseOutDir)

if not os.path.isdir(ENERGYDIR):
    os.makedirs(ENERGYDIR)

caseEnergyDir = os.path.join(ENERGYDIR, loadCaseName)
if not os.path.isdir(caseEnergyDir):
    os.makedirs(caseEnergyDir)

odbPath = JOB + '.odb'
odb = openOdb(path=odbPath)

step = odb.steps[STEP]
instNameODB = INST.upper()
odbInst = odb.rootAssembly.instances[instNameODB]

# ============================================================
# EXPORT ALLSE HISTORY FOR THIS TRANSLATION
# ============================================================
histRegionKey = None

for key in step.historyRegions.keys():
    if 'Assembly ASSEMBLY' in key or key == 'Assembly ASSEMBLY':
        histRegionKey = key
        break

if histRegionKey is None:
    # fallback: take first history region that contains ALLSE
    for key, reg in step.historyRegions.items():
        if 'ALLSE' in reg.historyOutputs.keys():
            histRegionKey = key
            break

if histRegionKey is None:
    print('WARNING: No history region containing ALLSE was found.')
else:
    allseData = step.historyRegions[histRegionKey].historyOutputs['ALLSE'].data
    energyCsvName = '%s_RUC_%dx%d_ALLSE.csv' % (JOB, NX, NY)
    energyCsvPath = os.path.join(caseEnergyDir, energyCsvName)
    with open(energyCsvPath, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['LoadCase', 'NX', 'NY', 'zeta1', 'zeta2', 'FrameIndex', 'StepTime', 'ALLSE'])
        for i, pair in enumerate(allseData):
            stepTime = pair[0]
            allseVal = pair[1]
            writer.writerow([loadCaseName, NX, NY, zeta1, zeta2, i, stepTime, allseVal])
    print('Saved energy file:', energyCsvPath)

nFrames = len(step.frames)
print('Number of frames to export =', nFrames)

for iFrame, frame in enumerate(step.frames):
    uField = frame.fieldOutputs['U']
    uSub = uField.getSubset(region=odbInst, position=NODAL)
    # build displacement dictionary by node label
    uDict = {}
    for v in uSub.values:
        uDict[v.nodeLabel] = v.data
    csvName = '%s_f%04d_Nodal.csv' % (JOB, iFrame)
    csvPath = os.path.join(caseOutDir, csvName)
    with open(csvPath, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['FrameID', 'StepTime', 'NodeLabel', 'X', 'Y', 'U1', 'U2'])
        for nd in odbInst.nodes:
            lab = nd.label
            x = nd.coordinates[0]
            y = nd.coordinates[1]
            if lab in uDict:
                u1 = uDict[lab][0]
                u2 = uDict[lab][1]
            else:
                u1 = 0.0
                u2 = 0.0
            writer.writerow([iFrame, frame.frameValue, lab, x, y, u1, u2])
    print('Saved:', csvPath)




# ============================================================
# EXPORT SIDE REACTION FORCES (8 components per frame)
# ============================================================
# Pull the node sets from the ODB (assembly-level)
odbA = odb.rootAssembly
sideSetNames = ['LEFT', 'RIGHT', 'BOTTOM', 'TOP']
sideSets = {nm: odbA.nodeSets[nm] for nm in sideSetNames}

rfCsvName = '%s_RUC_%dx%d_RF_sides.csv' % (JOB, NX, NY)
rfCsvPath = os.path.join(caseEnergyDir, rfCsvName)

with open(rfCsvPath, 'w') as f:
    writer = csv.writer(f)
    writer.writerow([
        'LoadCase','NX','NY','zeta1','zeta2','FrameIndex','StepTime',
        'RF_LEFT_x','RF_LEFT_y',
        'RF_RIGHT_x','RF_RIGHT_y',
        'RF_BOTTOM_x','RF_BOTTOM_y',
        'RF_TOP_x','RF_TOP_y',
    ])
    for iFrame, frame in enumerate(step.frames):
        rfField = frame.fieldOutputs['RF']
        row = [loadCaseName, NX, NY, zeta1, zeta2, iFrame, frame.frameValue]
        for nm in sideSetNames:
            rfSub = rfField.getSubset(region=sideSets[nm], position=NODAL)
            Fx, Fy = 0.0, 0.0
            for v in rfSub.values:
                Fx += v.data[0]
                Fy += v.data[1]
            row.extend([Fx, Fy])
        writer.writerow(row)

print('Saved side reaction file:', rfCsvPath)

odb.close()
print('Finished exporting all frames for job:', JOB)