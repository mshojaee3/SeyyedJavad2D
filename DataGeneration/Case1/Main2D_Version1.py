# -*- coding: utf-8 -*-
# ============================================================
#  Main2D_LE.py
#  Linear-elastic 2D RUC array with translated inclusions.
#  ONE static step, ONE increment, 1% nominal strain.
#  Exports ONLY the final frame -> one nodal CSV per translation.
#
#  Boundary conditions are generated from a second-order
#  polynomial displacement field
#
#       u_i(X) = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k ,   Xc = X - X_center
#
#  imposed on the complete outer boundary (KUBC-type), which
#  covers uniaxial / biaxial / shear / bending / shear-gradient
#  and any combination of them in a single implementation.
#  A few "free" (mixed-BC, specimen-like) cases are also kept.
#
#  Run:  abaqus cae noGUI=Main2D_LE.py
# ============================================================

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
# USER INPUT  (this block is overwritten by the MATLAB driver)
# ============================================================
loadCaseName = 'UniaxialX'

MODEL = 'Model-1'
PART  = 'RVE2D_FULL'

# translation of the inclusion pattern inside each unit cell
zeta1 = 0.0
zeta2 = 0.0

# output folders
OUTDIR    = r'OUTPUT_CSV'
ENERGYDIR = r'OUTPUT_ENERGY'



# ------------------------------------------------------------
# FIXED TOTAL SPECIMEN SIZE
# ------------------------------------------------------------
Lx_tot = 1.0          # Total specimen length
Ly_tot = 1.0          # Total specimen height

# ------------------------------------------------------------
# Number of repeated unit cells
# ------------------------------------------------------------
NX = 1
NY = 1


# geometry / mesh fractions (fractions of Lx)
Rfrac    = 0.30
meshFrac = 0.25

# loading amplitude: 1 percent
strain0 = 0.01

# linear elastic phases
E_m  = 70000.0
nu_m = 0.33
E_i  = 3500.0
nu_i = 0.33
rho_m = 1.12e-06
rho_i = 1.12e-06

STATE_2D = 'plane_strain'
TOL = 1e-6

# ============================================================
# DERIVED GEOMETRY
# ============================================================

# Unit-cell dimensions
Lx = Lx_tot / float(NX)
Ly = Ly_tot / float(NY)

# Inclusion radius
R = Rfrac * min(Lx, Ly)

# Mesh size
mesh_size = meshFrac * min(Lx, Ly)

# Plate centre
Xc = 0.5 * Lx_tot
Yc = 0.5 * Ly_tot

Lx_tot = NX * Lx
Ly_tot = NY * Ly
Xc     = 0.5 * Lx_tot
Yc     = 0.5 * Ly_tot

JOB = '%s_zx_%0.3f_zy_%0.3f' % (loadCaseName, zeta1, zeta2)
JOB = JOB.replace('.', 'p').replace('-', 'm')

STEP = 'Step-1'

# ============================================================
# LOAD CASE LIBRARY
# ============================================================
# e   : nominal strain amplitude (1%)
# kap : curvature-type amplitude, scaled so that the maximum
#       gradient-induced strain on the boundary is also ~ e
def zerosG():
    return [[[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]]]


def buildLoadCase(name):
    e   = strain0
    Lmx = max(Lx_tot, Ly_tot)
    kap = 2.0 * e / Lmx
    H = [[0.0, 0.0], [0.0, 0.0]]
    G = zerosG()
    # ---------- pure affine (H only) ----------
    if name == 'UniaxialX':
        H = [[e, 0.0], [0.0, 0.0]]
    elif name == 'UniaxialY':
        H = [[0.0, 0.0], [0.0, e]]
    elif name == 'CompressionX':
        H = [[-e, 0.0], [0.0, 0.0]]
    elif name == 'CompressionY':
        H = [[0.0, 0.0], [0.0, -e]]
    elif name == 'BiaxialTension':
        H = [[e, 0.0], [0.0, e]]
    elif name == 'BiaxialCompression':
        H = [[-e, 0.0], [0.0, -e]]
    elif name == 'BiaxialUnequal':
        H = [[e, 0.0], [0.0, -0.5 * e]]
    elif name == 'PureShearNormal':          # isochoric, normal form
        H = [[e, 0.0], [0.0, -e]]
    elif name == 'SimpleShearXY':
        H = [[0.0, e], [0.0, 0.0]]
    elif name == 'SimpleShearYX':
        H = [[0.0, 0.0], [e, 0.0]]
    elif name == 'PureShearSym':             # symmetric shear, no rotation
        H = [[0.0, 0.5 * e], [0.5 * e, 0.0]]
    elif name == 'RotationShear':            # antisymmetric part only
        H = [[0.0, 0.5 * e], [-0.5 * e, 0.0]]
    elif name == 'CombinedTensionShear':
        H = [[e, 0.5 * e], [0.0, 0.0]]
    elif name == 'CombinedBiaxialShear':
        H = [[e, 0.5 * e], [0.5 * e, -0.5 * e]]
    # ---------- pure gradient (G only) ----------
    elif name == 'PureBendingX':             # beam bending about z, axis x
        G[0][0][1] = -kap
        G[0][1][0] = -kap
        G[1][0][0] = kap
    elif name == 'PureBendingY':             # bending with axis y
        G[1][1][0] = -kap
        G[1][0][1] = -kap
        G[0][1][1] = kap
    elif name == 'ShearGradientX':           # u1 = 1/2 kap X2^2 (flexure-like)
        G[0][1][1] = kap
    elif name == 'ShearGradientY':           # u2 = 1/2 kap X1^2 (parabolic)
        G[1][0][0] = kap
    elif name == 'StretchGradientX':         # u1 = 1/2 kap X1^2
        G[0][0][0] = kap
    elif name == 'StretchGradientY':         # u2 = 1/2 kap X2^2
        G[1][1][1] = kap
    elif name == 'DilatationGradient':
        G[0][0][0] = kap
        G[1][1][1] = kap
    # ---------- combined affine + gradient ----------
    elif name == 'TensionPlusBending':
        H = [[e, 0.0], [0.0, 0.0]]
        G[0][0][1] = -kap
        G[0][1][0] = -kap
        G[1][0][0] = kap
    elif name == 'ShearPlusBending':
        H = [[0.0, e], [0.0, 0.0]]
        G[0][0][1] = -kap
        G[0][1][0] = -kap
        G[1][0][0] = kap
    elif name == 'BiaxialPlusShearGradient':
        H = [[e, 0.0], [0.0, -0.5 * e]]
        G[0][1][1] = kap
        G[1][0][0] = kap
    elif name == 'GeneralMixed':
        H = [[e, 0.5 * e], [0.25 * e, -0.5 * e]]
        G[0][0][1] = -kap
        G[0][1][0] = -kap
        G[1][0][0] = kap
        G[0][1][1] = 0.5 * kap
    # ---------- mixed / free specimen-like cases ----------
    elif name == 'FreeUniaxialX':
        return {'kind': 'mixed', 'sub': 'ux', 'val': e * Lx_tot}
    elif name == 'FreeUniaxialY':
        return {'kind': 'mixed', 'sub': 'uy', 'val': e * Ly_tot}
    elif name == 'FreeCompressionY':
        return {'kind': 'mixed', 'sub': 'uy', 'val': -e * Ly_tot}
    elif name == 'FreeSimpleShear':
        return {'kind': 'mixed', 'sub': 'shear', 'val': e * Ly_tot}
    else:
        raise RuntimeError('Unknown loadCaseName: %s' % name)
    return {'kind': 'affine', 'H': H, 'G': G}


LC = buildLoadCase(loadCaseName)

# ============================================================
# MODEL
# ============================================================
if MODEL not in mdb.models:
    mdb.Model(name=MODEL)

model = mdb.models[MODEL]

if JOB in mdb.jobs:
    del mdb.jobs[JOB]

# ============================================================
# GEOMETRY
# ============================================================
if PART in model.parts:
    del model.parts[PART]

sk = model.ConstrainedSketch(name='__profile__', sheetSize=10.0 * max(Lx_tot, Ly_tot))
sk.rectangle(point1=(0.0, 0.0), point2=(Lx_tot, Ly_tot))

p = model.Part(name=PART, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p.BaseShell(sketch=sk)
del model.sketches['__profile__']

f0 = p.faces[0]
tr = p.MakeSketchTransform(sketchPlane=f0, sketchPlaneSide=SIDE1,
                           sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0))
sk2 = model.ConstrainedSketch(name='__profile__', sheetSize=10.0 * max(Lx_tot, Ly_tot),
                              transform=tr)
p.projectReferencesOntoSketch(sketch=sk2, filter=COPLANAR_EDGES)


def circle_intersects_domain(cx, cy, Rc, xmin, xmax, ymin, ymax):
    return not (cx + Rc < xmin or cx - Rc > xmax or cy + Rc < ymin or cy - Rc > ymax)


for i in range(-1, NX + 1):
    for j in range(-1, NY + 1):
        cx = (i + 0.5) * Lx + zeta1
        cy = (j + 0.5) * Ly + zeta2
        if circle_intersects_domain(cx, cy, R, 0.0, Lx_tot, 0.0, Ly_tot):
            sk2.CircleByCenterPerimeter(center=(cx, cy), point1=(cx + R, cy))

p.PartitionFaceBySketch(faces=(f0,), sketch=sk2)
del model.sketches['__profile__']

# ============================================================
# MATERIALS + SECTIONS  (linear elastic only)
# ============================================================
for m in ('Matrix', 'Inclusion'):
    if m in model.materials.keys():
        del model.materials[m]

model.Material(name='Matrix')
model.materials['Matrix'].Density(table=((rho_m,),))
model.materials['Matrix'].Elastic(table=((E_m, nu_m),))

model.Material(name='Inclusion')
model.materials['Inclusion'].Density(table=((rho_i,),))
model.materials['Inclusion'].Elastic(table=((E_i, nu_i),))

for s in ('Sec-Matrix', 'Sec-Inclusion'):
    if s in model.sections.keys():
        del model.sections[s]

model.HomogeneousSolidSection(name='Sec-Matrix', material='Matrix', thickness=1.0)
model.HomogeneousSolidSection(name='Sec-Inclusion', material='Inclusion', thickness=1.0)

if 'SET_ALL' in p.sets.keys():
    del p.sets['SET_ALL']

p.Set(name='SET_ALL', faces=p.faces[:])
p.SectionAssignment(region=p.sets['SET_ALL'], sectionName='Sec-Matrix')

# inclusion faces by centroid test
inc_idx = []
for k, face in enumerate(p.faces):
    c = face.getCentroid()[0]
    x, y = c[0], c[1]
    inside_any = False
    for i in range(-1, NX + 1):
        for j in range(-1, NY + 1):
            cx = (i + 0.5) * Lx + zeta1
            cy = (j + 0.5) * Ly + zeta2
            if abs(face.getSize() - math.pi*R**2) < 1e-6:
                inside_any = True
                break
        if inside_any:
            break
    if inside_any:
        inc_idx.append(k)

print('Number of faces in part        =', len(p.faces))
print('Number of inclusion faces found =', len(inc_idx))

if len(inc_idx) == 0:
    raise RuntimeError('No inclusion faces found. Check zeta and geometry.')

inc_faces = p.faces[inc_idx[0]:inc_idx[0] + 1]
for k in inc_idx[1:]:
    inc_faces = inc_faces + p.faces[k:k + 1]

if 'SET_INC' in p.sets.keys():
    del p.sets['SET_INC']

p.Set(name='SET_INC', faces=inc_faces)
p.SectionAssignment(region=p.sets['SET_INC'], sectionName='Sec-Inclusion')

# ============================================================
# ASSEMBLY + STEP  (single increment, geometrically linear)
# ============================================================
a = model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)

for nm in list(a.instances.keys()):
    del a.instances[nm]

INST = PART + '-1'
inst = a.Instance(name=INST, part=p, dependent=ON)

model.StaticStep(name=STEP, previous='Initial', nlgeom=OFF,
                 timePeriod=1.0, maxNumInc=1,
                 initialInc=1.0, minInc=1.0, maxInc=1.0)

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
a.regenerate()
inst = a.instances[INST]

# ============================================================
# BOUNDARY NODE SETS
# ============================================================
left_labs, right_labs, bottom_labs, top_labs, bnd_labs = [], [], [], [], []
corner00 = None

for nd in inst.nodes:
    x, y, z = nd.coordinates
    onL = abs(x - 0.0) < TOL
    onR = abs(x - Lx_tot) < TOL
    onB = abs(y - 0.0) < TOL
    onT = abs(y - Ly_tot) < TOL
    if onL:
        left_labs.append(nd.label)
    if onR:
        right_labs.append(nd.label)
    if onB:
        bottom_labs.append(nd.label)
    if onT:
        top_labs.append(nd.label)
    if onL or onR or onB or onT:
        bnd_labs.append(nd.label)
    if onL and onB:
        corner00 = nd.label

print('Boundary counts: LEFT %d RIGHT %d BOTTOM %d TOP %d ALL %d'
      % (len(left_labs), len(right_labs), len(bottom_labs), len(top_labs), len(bnd_labs)))

if min(len(left_labs), len(right_labs), len(bottom_labs), len(top_labs)) == 0:
    raise RuntimeError('One boundary set is empty. Increase TOL.')

if corner00 is None:
    raise RuntimeError('Corner (0,0) not found. Increase TOL.')

a.Set(name='LEFT',     nodes=inst.nodes.sequenceFromLabels(tuple(left_labs)))
a.Set(name='RIGHT',    nodes=inst.nodes.sequenceFromLabels(tuple(right_labs)))
a.Set(name='BOTTOM',   nodes=inst.nodes.sequenceFromLabels(tuple(bottom_labs)))
a.Set(name='TOP',      nodes=inst.nodes.sequenceFromLabels(tuple(top_labs)))
a.Set(name='BND_ALL',  nodes=inst.nodes.sequenceFromLabels(tuple(sorted(set(bnd_labs)))))
a.Set(name='CORNER00', nodes=inst.nodes.sequenceFromLabels((corner00,)))

# ============================================================
# BOUNDARY CONDITIONS
# ============================================================
def polyExpression(i, H, G):
    """String for u_i(X) = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k."""
    cX  = H[i][0]
    cY  = H[i][1]
    cXX = 0.5 * G[i][0][0]
    cXY = 0.5 * (G[i][0][1] + G[i][1][0])
    cYY = 0.5 * G[i][1][1]
    t = []
    if abs(cX) > 0.0:
        t.append('(%.12g)*(X-(%.12g))' % (cX, Xc))
    if abs(cY) > 0.0:
        t.append('(%.12g)*(Y-(%.12g))' % (cY, Yc))
    if abs(cXX) > 0.0:
        t.append('(%.12g)*(X-(%.12g))*(X-(%.12g))' % (cXX, Xc, Xc))
    if abs(cXY) > 0.0:
        t.append('(%.12g)*(X-(%.12g))*(Y-(%.12g))' % (cXY, Xc, Yc))
    if abs(cYY) > 0.0:
        t.append('(%.12g)*(Y-(%.12g))*(Y-(%.12g))' % (cYY, Yc, Yc))
    # empty string signals an identically zero component: an Abaqus
    # ExpressionField must contain at least one of X, Y, Z
    if len(t) == 0:
        return ''
    return ' + '.join(t)


if LC['kind'] == 'affine':
    H, G = LC['H'], LC['G']
    ex = polyExpression(0, H, G)
    ey = polyExpression(1, H, G)
    print('u1 field = ' + (ex if ex else '0 (constant)'))
    print('u2 field = ' + (ey if ey else '0 (constant)'))
    if ex:
        model.ExpressionField(name='UX_FIELD', localCsys=None, expression=ex)
        model.DisplacementBC(name='BC_BND_U1', createStepName=STEP,
                             region=a.sets['BND_ALL'],
                             u1=1.0, u2=UNSET, ur3=UNSET,
                             distributionType=FIELD, fieldName='UX_FIELD')
    else:
        model.DisplacementBC(name='BC_BND_U1', createStepName=STEP,
                             region=a.sets['BND_ALL'],
                             u1=0.0, u2=UNSET, ur3=UNSET)
    if ey:
        model.ExpressionField(name='UY_FIELD', localCsys=None, expression=ey)
        model.DisplacementBC(name='BC_BND_U2', createStepName=STEP,
                             region=a.sets['BND_ALL'],
                             u1=UNSET, u2=1.0, ur3=UNSET,
                             distributionType=FIELD, fieldName='UY_FIELD')
    else:
        model.DisplacementBC(name='BC_BND_U2', createStepName=STEP,
                             region=a.sets['BND_ALL'],
                             u1=UNSET, u2=0.0, ur3=UNSET)
else:  # mixed / free specimen-like
    sub = LC['sub']
    val = LC['val']
    if sub == 'ux':
        model.DisplacementBC(name='BC_LEFT', createStepName=STEP,
                             region=a.sets['LEFT'], u1=0.0, u2=UNSET, ur3=UNSET)
        model.DisplacementBC(name='BC_RIGHT', createStepName=STEP,
                             region=a.sets['RIGHT'], u1=val, u2=UNSET, ur3=UNSET)
        model.DisplacementBC(name='BC_PIN', createStepName='Initial',
                             region=a.sets['CORNER00'], u1=UNSET, u2=0.0, ur3=UNSET)
    elif sub == 'uy':
        model.DisplacementBC(name='BC_BOTTOM', createStepName=STEP,
                             region=a.sets['BOTTOM'], u1=UNSET, u2=0.0, ur3=UNSET)
        model.DisplacementBC(name='BC_TOP', createStepName=STEP,
                             region=a.sets['TOP'], u1=UNSET, u2=val, ur3=UNSET)
        model.DisplacementBC(name='BC_PIN', createStepName='Initial',
                             region=a.sets['CORNER00'], u1=0.0, u2=UNSET, ur3=UNSET)
    elif sub == 'shear':
        model.DisplacementBC(name='BC_BOTTOM', createStepName=STEP,
                             region=a.sets['BOTTOM'], u1=0.0, u2=0.0, ur3=UNSET)
        model.DisplacementBC(name='BC_TOP', createStepName=STEP,
                             region=a.sets['TOP'], u1=val, u2=0.0, ur3=UNSET)
    else:
        raise RuntimeError('Unknown mixed sub-case: %s' % sub)

model.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'E', 'U', 'RF', 'COORD', 'IVOL'), frequency=LAST_INCREMENT
)

# ============================================================
# JOB
# ============================================================
mdb.Job(name=JOB, model=MODEL, type=ANALYSIS,
        memory=90, memoryUnits=PERCENTAGE,
        numCpus=1, numGPUs=0,
        explicitPrecision=SINGLE, nodalOutputPrecision=FULL)

mdb.jobs[JOB].submit(consistencyChecking=OFF)
mdb.jobs[JOB].waitForCompletion()

# ============================================================
# EXPORT: FINAL FRAME ONLY  -> one nodal CSV per translation
# ============================================================
from odbAccess import openOdb

caseOutDir = os.path.join(OUTDIR, loadCaseName)
caseEnergyDir = os.path.join(ENERGYDIR, loadCaseName)
for d in (OUTDIR, caseOutDir, ENERGYDIR, caseEnergyDir):
    if not os.path.isdir(d):
        os.makedirs(d)

odb = openOdb(path=JOB + '.odb')
step = odb.steps[STEP]
odbInst = odb.rootAssembly.instances[INST.upper()]

frame = step.frames[-1]
print('Exporting final frame, frameValue =', frame.frameValue)

def getData(v):
    """Double precision data if available, single otherwise."""
    try:
        return v.dataDouble
    except:
        return v.data


# fixed component schema so the CSV header never changes
S_COMPS = ['S11', 'S22', 'S33', 'S12']
E_COMPS = ['E11', 'E22', 'E33', 'E12']

# ---- nodal displacement ------------------------------------
uSub = frame.fieldOutputs['U'].getSubset(region=odbInst, position=NODAL)
uDict = {}
for v in uSub.values:
    uDict[v.nodeLabel] = getData(v)


# ---- nodal stress / strain (element-nodal, averaged per node)
def nodalAveragedTensor(fieldName, wantedComps):
    """Extrapolate to nodes and average contributions of all elements
    sharing a node. Returns {nodeLabel: [c1, c2, ...]} in wantedComps order."""
    fld = frame.fieldOutputs[fieldName].getSubset(region=odbInst,
                                                  position=ELEMENT_NODAL)
    labels = list(fld.componentLabels)
    # index of each wanted component inside the ODB data vector
    idx = []
    for c in wantedComps:
        idx.append(labels.index(c) if c in labels else -1)
    acc, cnt = {}, {}
    for v in fld.values:
        lab = v.nodeLabel
        d = getData(v)
        if lab not in acc:
            acc[lab] = [0.0] * len(wantedComps)
            cnt[lab] = 0
        for q, ii in enumerate(idx):
            if ii >= 0:
                acc[lab][q] += d[ii]
        cnt[lab] += 1
    for lab in acc.keys():
        n = float(cnt[lab])
        acc[lab] = [val / n for val in acc[lab]]
    return acc


sDict = nodalAveragedTensor('S', S_COMPS)
eDict = nodalAveragedTensor('E', E_COMPS)

csvPath = os.path.join(caseOutDir, '%s_Nodal.csv' % JOB)
fcsv = open(csvPath, 'w')
w = csv.writer(fcsv)
w.writerow(['LoadCase', 'zeta1', 'zeta2', 'StepTime', 'NodeLabel', 'X', 'Y',
            'U1', 'U2'] + S_COMPS + E_COMPS)
for nd in odbInst.nodes:
    lab = nd.label
    x = nd.coordinates[0]
    y = nd.coordinates[1]
    if lab in uDict:
        u1, u2 = uDict[lab][0], uDict[lab][1]
    else:
        u1, u2 = 0.0, 0.0
    srow = sDict.get(lab, [0.0] * len(S_COMPS))
    erow = eDict.get(lab, [0.0] * len(E_COMPS))
    w.writerow([loadCaseName, zeta1, zeta2, frame.frameValue, lab, x, y, u1, u2]
               + list(srow) + list(erow))

fcsv.close()
print('Saved nodal file:', csvPath)

# ------------------------------------------------------------
# VOLUME AVERAGES over the whole domain
#   <A> = sum_gp (A_gp * IVOL_gp) / sum_gp IVOL_gp
# ------------------------------------------------------------
def volumeAverage(fieldName, wantedComps):
    fld = frame.fieldOutputs[fieldName].getSubset(region=odbInst,
                                                  position=INTEGRATION_POINT)
    vol = frame.fieldOutputs['IVOL'].getSubset(region=odbInst,
                                               position=INTEGRATION_POINT)
    volMap = {}
    for v in vol.values:
        d = getData(v)
        volMap[(v.elementLabel, v.integrationPoint)] = d if isinstance(d, float) else d[0]
        
        
    labels = list(fld.componentLabels)
    idx = []
    for c in wantedComps:
        idx.append(labels.index(c) if c in labels else -1)
    acc = [0.0] * len(wantedComps)
    Vtot = 0.0
    for v in fld.values:
        dV = volMap.get((v.elementLabel, v.integrationPoint), 0.0)
        if dV == 0.0:
            continue
        d = getData(v)
        for q, ii in enumerate(idx):
            if ii >= 0:
                acc[q] += d[ii] * dV
        Vtot += dV
    if Vtot <= 0.0:
        print('WARNING: total IVOL is zero for field %s' % fieldName)
        return [float('nan')] * len(wantedComps), 0.0
    return [val / Vtot for val in acc], Vtot


Savg, Vtot = volumeAverage('S', S_COMPS)
Eavg, _    = volumeAverage('E', E_COMPS)

# mean nodal displacement (unweighted arithmetic mean over all nodes)
nU = 0
u1sum, u2sum = 0.0, 0.0
for lab in uDict.keys():
    u1sum += uDict[lab][0]
    u2sum += uDict[lab][1]
    nU += 1

Uavg = [u1sum / nU, u2sum / nU] if nU > 0 else [float('nan'), float('nan')]

print('Volume  = %.6e' % Vtot)
print('<S>     =', Savg)
print('<E>     =', Eavg)
print('<U>     =', Uavg)

# ------------------------------------------------------------
# strain energy (single value) + side reaction forces
# ------------------------------------------------------------
allseVal = float('nan')
histKey = None
for key, reg in step.historyRegions.items():
    if 'ALLSE' in reg.historyOutputs.keys():
        histKey = key
        break

if histKey is not None:
    dataALLSE = step.historyRegions[histKey].historyOutputs['ALLSE'].data
    allseVal = dataALLSE[-1][1]
else:
    print('WARNING: ALLSE history not found.')

odbA = odb.rootAssembly
sideSetNames = ['LEFT', 'RIGHT', 'BOTTOM', 'TOP']
rfField = frame.fieldOutputs['RF']
rfRow = []
for nm in sideSetNames:
    rfSub = rfField.getSubset(region=odbA.nodeSets[nm], position=NODAL)
    Fx, Fy = 0.0, 0.0
    for v in rfSub.values:
        Fx += v.dataDouble[0]
        Fy += v.dataDouble[1]
    rfRow.extend([Fx, Fy])

sumPath = os.path.join(caseEnergyDir, '%s_Summary.csv' % JOB)
fsum = open(sumPath, 'w')
w = csv.writer(fsum)
w.writerow(['LoadCase', 'NX', 'NY', 'zeta1', 'zeta2', 'strain0', 'StepTime', 'ALLSE',
            'RF_LEFT_x', 'RF_LEFT_y', 'RF_RIGHT_x', 'RF_RIGHT_y',
            'RF_BOTTOM_x', 'RF_BOTTOM_y', 'RF_TOP_x', 'RF_TOP_y',
            'Volume',
            'Savg_11', 'Savg_22', 'Savg_33', 'Savg_12',
            'Eavg_11', 'Eavg_22', 'Eavg_33', 'Eavg_12',
            'Uavg_1', 'Uavg_2'])
w.writerow([loadCaseName, NX, NY, zeta1, zeta2, strain0, frame.frameValue, allseVal]
           + rfRow + [Vtot] + Savg + Eavg + Uavg)
fsum.close()
print('Saved summary file:', sumPath)

odb.close()
print('Finished job:', JOB)
