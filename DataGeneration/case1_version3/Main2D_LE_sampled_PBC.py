# -*- coding: utf-8 -*-
# ============================================================
#  Main2D_LE_sampled_PBC.py
#  Linear-elastic 2D RUC array with translated inclusions.
#  ONE static step, ONE increment, 1% nominal strain.
#  Exports ONLY the final frame -> one nodal CSV per translation.
#
#  BOUNDARY CONDITIONS -- PERIODIC (PBC), replacing the previous
#  full-boundary KUBC (kinematic uniform BC / quadratic expression
#  field). Only the affine (first-order) macroscopic strain H_ij is
#  applied; the gradient part G is NOT enforced under PBC (a warning
#  is printed if a load case that specifies nonzero G is run here).
#
#  Formulation: fluctuation u* = u - H_ij Xc_j is enforced periodic
#  node-by-node between opposite RVE boundaries via *Equation
#  constraints, and the macroscopic strain H is driven through two
#  dummy reference points (RP_X, RP_Y):
#
#       u_i(RIGHT)  - u_i(LEFT)   = H_i1 * Lx_tot   (= RP_X u_i)
#       u_i(TOP)    - u_i(BOTTOM) = H_i2 * Ly_tot   (= RP_Y u_i)
#
#  with the four corners tied by the usual parallelogram closure and
#  one corner (CORNER00) pinned to remove rigid-body translation.
#  This requires opposite RVE edges to mesh IDENTICALLY (same node
#  y/x-coordinates on LEFT/RIGHT and BOTTOM/TOP) -- see the MESH
#  section below, which seeds matching boundary sub-edges by number
#  (constraint=FIXED) before the free interior mesh is generated.
#
#  The "mixed" (Free*, specimen-like) load cases are unchanged and do
#  not use PBC.
#
#  Run:  abaqus cae noGUI=Main2D_LE_sampled_PBC.py
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
loadCaseName = 'StretchGradientX'

# ------------------------------------------------------------
# Sampled-case coefficients (point-set library, e.g. Run_LE_Sampled
# _Ensemble.m + sample_points_d09*.txt). Used ONLY when loadCaseName
# starts with SAMPLED_PREFIX below; ignored for all named cases in
# the library further down. Values are the already-scaled physical
# components (H in strain units, G in 1/length units), computed in
# MATLAB from the raw unit-sphere sample point.
# ------------------------------------------------------------
H11   = 0.0
H22   = 0.0
H12   = 0.0
G1_11 = 0.0
G1_22 = 0.0
G1_12 = 0.0
G2_11 = 0.0
G2_22 = 0.0
G2_12 = 0.0

# ------------------------------------------------------------
# Set True to also write the per-node CSV (X,Y,U,S,E at every node).
# NOT needed for Fbar/G (boundary-only least-squares fit + surface
# integral, computed below directly from the ODB) or for Pbar/Qbar
# (IVOL-weighted volumeAverage/volumeAverageAndMoment, also computed
# directly from the ODB at INTEGRATION_POINT). Turn this on only if
# you separately want the raw per-node fields for plotting or another
# diagnostic -- leaving it False skips a large per-node write and
# substantially speeds up a full sweep.
# ------------------------------------------------------------
EXPORT_NODAL_CSV = False

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

# geometry / mesh fractions
Rfrac    = 6/19
meshFrac = 0.05

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

# Unit-cell dimensions (total size is fixed, number of cells changes)
Lx = Lx_tot / float(NX)
Ly = Ly_tot / float(NY)

# Inclusion radius
R = Rfrac * min(Lx, Ly)

# Mesh size
mesh_size = meshFrac * min(Lx, Ly)

# Plate centre
Xc = 0.5 * Lx_tot
Yc = 0.5 * Ly_tot

JOB = '%s_zx_%0.3f_zy_%0.3f' % (loadCaseName, zeta1, zeta2)
JOB = JOB.replace('.', 'p').replace('-', 'm')

STEP = 'Step-1'

# ============================================================
# LOAD CASE LIBRARY
# ============================================================
def zerosG():
    return [[[0.0, 0.0], [0.0, 0.0]], [[0.0, 0.0], [0.0, 0.0]]]


SAMPLED_PREFIX = 'NC_k'   # must match CASE_PREFIX in the MATLAB driver


def buildLoadCase(name):
    e   = strain0
    Lmx = max(Lx_tot, Ly_tot)
    kap = 2.0 * e / Lmx
    H = [[0.0, 0.0], [0.0, 0.0]]
    G = zerosG()
    # ---------- sampled case (point-set library) ----------
    # H symmetric (3 dof) + G symmetric in its last two indices (6 dof),
    # i.e. the same 'reduced9' convention as pann_loadcases.py. The
    # component values themselves already carry the strain0 / bending
    # scaling; buildLoadCase only assembles them into full tensors.
    if name.startswith(SAMPLED_PREFIX):
        H[0][0] = H11
        H[1][1] = H22
        H[0][1] = H12
        H[1][0] = H12

        G[0][0][0] = G1_11
        G[0][1][1] = G1_22
        G[0][0][1] = G1_12
        G[0][1][0] = G1_12

        G[1][0][0] = G2_11
        G[1][1][1] = G2_22
        G[1][0][1] = G2_12
        G[1][1][0] = G2_12

        return {'kind': 'affine', 'H': H, 'G': G}
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
    elif name == 'PureShearNormal':
        H = [[e, 0.0], [0.0, -e]]
    elif name == 'SimpleShearXY':
        H = [[0.0, e], [0.0, 0.0]]
    elif name == 'SimpleShearYX':
        H = [[0.0, 0.0], [e, 0.0]]
    elif name == 'PureShearSym':
        H = [[0.0, 0.5 * e], [0.5 * e, 0.0]]
    elif name == 'RotationShear':
        H = [[0.0, 0.5 * e], [-0.5 * e, 0.0]]
    elif name == 'CombinedTensionShear':
        H = [[e, 0.5 * e], [0.0, 0.0]]
    elif name == 'CombinedBiaxialShear':
        H = [[e, 0.5 * e], [0.5 * e, -0.5 * e]]
    # ---------- pure gradient (G only) ----------
    elif name == 'PureBendingX':
        G[0][0][1] = -kap
        G[0][1][0] = -kap
        G[1][0][0] = kap
    elif name == 'PureBendingY':
        G[1][1][0] = -kap
        G[1][0][1] = -kap
        G[0][1][1] = kap
    elif name == 'ShearGradientX':
        G[0][1][1] = kap
    elif name == 'ShearGradientY':
        G[1][0][0] = kap
    elif name == 'StretchGradientX':
        G[0][0][0] = kap
    elif name == 'StretchGradientY':
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


# draw periodic translated pattern
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

# ------------------------------------------------------------
# inclusion faces by area test
# Full circles have area pi*R^2.
# This matches the interior inclusions associated with the NX*NY cells.
# ------------------------------------------------------------
inc_idx = []
target_area = math.pi * R**2
area_tol = max(1.0e-8, 1.0e-4 * target_area)

for k, face in enumerate(p.faces):
    A = face.getSize()
    if abs(A - target_area) < area_tol:
        inc_idx.append(k)

print('Number of faces in part         =', len(p.faces))
print('Number of inclusion faces found =', len(inc_idx))
print('Expected full inclusion faces   =', NX * NY)

if len(inc_idx) == 0:
    raise RuntimeError('No inclusion faces found. Check NX, NY, Rfrac, zeta and geometry.')

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
USE_PBC = (LC['kind'] == 'affine')


def _edgeEndpoints(edge):
    vs = edge.getVertices()
    if len(vs) != 2:
        return None
    p0 = p.vertices[vs[0]].pointOn[0]
    p1 = p.vertices[vs[1]].pointOn[0]
    return p0, p1


def _isStraight(edge, p0, p1):
    chord = math.sqrt((p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2)
    if chord < 1.0e-12:
        return False
    L = edge.getSize()
    return abs(L - chord) < 1.0e-6 * chord


def _classifyBoundaryEdges():
    """Straight sub-edges lying exactly on each of the four RVE
    boundary lines (there can be more than one per side where an
    inclusion circle intersects the boundary). Returns four lists of
    (edge, coordMin, coordMax), sorted by coordMin."""
    left, right, bottom, top = [], [], [], []
    for e in p.edges:
        ends = _edgeEndpoints(e)
        if ends is None:
            continue
        p0, p1 = ends
        if not _isStraight(e, p0, p1):
            continue
        x0, y0 = p0[0], p0[1]
        x1, y1 = p1[0], p1[1]
        if abs(x0) < TOL and abs(x1) < TOL:
            left.append((e, min(y0, y1), max(y0, y1)))
        elif abs(x0 - Lx_tot) < TOL and abs(x1 - Lx_tot) < TOL:
            right.append((e, min(y0, y1), max(y0, y1)))
        elif abs(y0) < TOL and abs(y1) < TOL:
            bottom.append((e, min(x0, x1), max(x0, x1)))
        elif abs(y0 - Ly_tot) < TOL and abs(y1 - Ly_tot) < TOL:
            top.append((e, min(x0, x1), max(x0, x1)))
    left.sort(key=lambda t: t[1])
    right.sort(key=lambda t: t[1])
    bottom.sort(key=lambda t: t[1])
    top.sort(key=lambda t: t[1])
    return left, right, bottom, top


if USE_PBC:
    leftE, rightE, bottomE, topE = _classifyBoundaryEdges()

    if len(leftE) != len(rightE):
        raise RuntimeError('PBC mesh: LEFT has %d boundary sub-edges but '
                           'RIGHT has %d -- the inclusion pattern is not '
                           'mirrored between the two sides.'
                           % (len(leftE), len(rightE)))
    if len(bottomE) != len(topE):
        raise RuntimeError('PBC mesh: BOTTOM has %d boundary sub-edges but '
                           'TOP has %d -- the inclusion pattern is not '
                           'mirrored between the two sides.'
                           % (len(bottomE), len(topE)))

    for (eL, yL0, yL1), (eR, yR0, yR1) in zip(leftE, rightE):
        if abs(yL0 - yR0) > TOL or abs(yL1 - yR1) > TOL:
            raise RuntimeError('PBC mesh: LEFT sub-edge [%.6g, %.6g] has no '
                               'matching RIGHT sub-edge (closest RIGHT '
                               'candidate [%.6g, %.6g]).'
                               % (yL0, yL1, yR0, yR1))
        n = max(1, int(round((yL1 - yL0) / mesh_size)))
        p.seedEdgeByNumber(edges=(eL,), number=n, constraint=FIXED)
        p.seedEdgeByNumber(edges=(eR,), number=n, constraint=FIXED)

    for (eB, xB0, xB1), (eT, xT0, xT1) in zip(bottomE, topE):
        if abs(xB0 - xT0) > TOL or abs(xB1 - xT1) > TOL:
            raise RuntimeError('PBC mesh: BOTTOM sub-edge [%.6g, %.6g] has no '
                               'matching TOP sub-edge (closest TOP '
                               'candidate [%.6g, %.6g]).'
                               % (xB0, xB1, xT0, xT1))
        n = max(1, int(round((xB1 - xB0) / mesh_size)))
        p.seedEdgeByNumber(edges=(eB,), number=n, constraint=FIXED)
        p.seedEdgeByNumber(edges=(eT,), number=n, constraint=FIXED)

    print('PBC mesh: matched %d LEFT/RIGHT sub-edge pair(s), %d BOTTOM/TOP '
          'sub-edge pair(s) with equal, FIXED element counts.'
          % (len(leftE), len(bottomE)))

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
if LC['kind'] == 'affine':
    # -------- PERIODIC BOUNDARY CONDITIONS (PBC) --------
    # Only H (affine part) is applied; G is intentionally ignored here.
    H = LC['H']
    Gnz = any(abs(LC['G'][i][j][k]) > 0.0
             for i in range(2) for j in range(2) for k in range(2))
    if Gnz:
        print('WARNING: loadCaseName=%s specifies nonzero G, but PBC only '
              'drives the affine part H through the reference points -- '
              'the gradient part G is NOT applied under PBC and is '
              'silently dropped.' % loadCaseName)

    # ---- locate the other three corners (corner00 already found) ----
    corner10 = corner01 = corner11 = None
    for nd in inst.nodes:
        x, y, z = nd.coordinates
        if abs(x - Lx_tot) < TOL and abs(y - 0.0) < TOL:
            corner10 = nd.label
        elif abs(x - 0.0) < TOL and abs(y - Ly_tot) < TOL:
            corner01 = nd.label
        elif abs(x - Lx_tot) < TOL and abs(y - Ly_tot) < TOL:
            corner11 = nd.label
    if None in (corner10, corner01, corner11):
        raise RuntimeError('PBC: could not find all four corner nodes.')

    a.Set(name='CORNER10', nodes=inst.nodes.sequenceFromLabels((corner10,)))
    a.Set(name='CORNER01', nodes=inst.nodes.sequenceFromLabels((corner01,)))
    a.Set(name='CORNER11', nodes=inst.nodes.sequenceFromLabels((corner11,)))
    cornerLabs = set([corner00, corner10, corner01, corner11])

    # ---- edge-interior node lists (corners excluded), paired by coord ----
    leftOnly = sorted([l for l in left_labs if l not in cornerLabs],
                      key=lambda l: inst.nodes.getFromLabel(l).coordinates[1])
    rightOnly = sorted([l for l in right_labs if l not in cornerLabs],
                       key=lambda l: inst.nodes.getFromLabel(l).coordinates[1])
    bottomOnly = sorted([l for l in bottom_labs if l not in cornerLabs],
                        key=lambda l: inst.nodes.getFromLabel(l).coordinates[0])
    topOnly = sorted([l for l in top_labs if l not in cornerLabs],
                     key=lambda l: inst.nodes.getFromLabel(l).coordinates[0])

    if len(leftOnly) != len(rightOnly):
        raise RuntimeError('PBC: LEFT has %d non-corner boundary nodes but '
                           'RIGHT has %d -- mesh is not periodic.'
                           % (len(leftOnly), len(rightOnly)))
    if len(bottomOnly) != len(topOnly):
        raise RuntimeError('PBC: BOTTOM has %d non-corner boundary nodes but '
                           'TOP has %d -- mesh is not periodic.'
                           % (len(bottomOnly), len(topOnly)))

    for lab_l, lab_r in zip(leftOnly, rightOnly):
        yl = inst.nodes.getFromLabel(lab_l).coordinates[1]
        yr = inst.nodes.getFromLabel(lab_r).coordinates[1]
        if abs(yl - yr) > TOL:
            raise RuntimeError('PBC: LEFT node %d (y=%.9g) has no matching '
                               'RIGHT node at the same y (closest RIGHT '
                               'candidate y=%.9g). Check mesh periodicity.'
                               % (lab_l, yl, yr))
    for lab_b, lab_t in zip(bottomOnly, topOnly):
        xb = inst.nodes.getFromLabel(lab_b).coordinates[0]
        xt = inst.nodes.getFromLabel(lab_t).coordinates[0]
        if abs(xb - xt) > TOL:
            raise RuntimeError('PBC: BOTTOM node %d (x=%.9g) has no matching '
                               'TOP node at the same x (closest TOP '
                               'candidate x=%.9g). Check mesh periodicity.'
                               % (lab_b, xb, xt))

    print('PBC: paired %d LEFT/RIGHT node(s), %d BOTTOM/TOP node(s), '
          '4 corners.' % (len(leftOnly), len(bottomOnly)))

    # ---- dummy reference points carrying the macroscopic strain H ----
    rpX_id = a.ReferencePoint(point=(Lx_tot + 0.25 * Lx_tot,
                                     -0.25 * Ly_tot, 0.0)).id
    rpY_id = a.ReferencePoint(point=(-0.25 * Lx_tot,
                                     Ly_tot + 0.25 * Ly_tot, 0.0)).id
    a.Set(name='RP_X', referencePoints=(a.referencePoints[rpX_id],))
    a.Set(name='RP_Y', referencePoints=(a.referencePoints[rpY_id],))

    # ---- periodicity equations: u(+side) - u(-side) - RP = 0 ----
    for k, (lab_l, lab_r) in enumerate(zip(leftOnly, rightOnly)):
        nL, nR = 'PBC_L_%d' % k, 'PBC_R_%d' % k
        a.Set(name=nL, nodes=inst.nodes.sequenceFromLabels((lab_l,)))
        a.Set(name=nR, nodes=inst.nodes.sequenceFromLabels((lab_r,)))
        model.Equation(name='EQ_LR_U1_%d' % k,
                       terms=((1.0, nR, 1), (-1.0, nL, 1), (-1.0, 'RP_X', 1)))
        model.Equation(name='EQ_LR_U2_%d' % k,
                       terms=((1.0, nR, 2), (-1.0, nL, 2), (-1.0, 'RP_X', 2)))

    for k, (lab_b, lab_t) in enumerate(zip(bottomOnly, topOnly)):
        nB, nT = 'PBC_B_%d' % k, 'PBC_T_%d' % k
        a.Set(name=nB, nodes=inst.nodes.sequenceFromLabels((lab_b,)))
        a.Set(name=nT, nodes=inst.nodes.sequenceFromLabels((lab_t,)))
        model.Equation(name='EQ_BT_U1_%d' % k,
                       terms=((1.0, nT, 1), (-1.0, nB, 1), (-1.0, 'RP_Y', 1)))
        model.Equation(name='EQ_BT_U2_%d' % k,
                       terms=((1.0, nT, 2), (-1.0, nB, 2), (-1.0, 'RP_Y', 2)))

    # ---- corner closure (parallelogram) ----
    model.Equation(name='EQ_C10_U1', terms=((1.0, 'CORNER10', 1),
                   (-1.0, 'CORNER00', 1), (-1.0, 'RP_X', 1)))
    model.Equation(name='EQ_C10_U2', terms=((1.0, 'CORNER10', 2),
                   (-1.0, 'CORNER00', 2), (-1.0, 'RP_X', 2)))
    model.Equation(name='EQ_C01_U1', terms=((1.0, 'CORNER01', 1),
                   (-1.0, 'CORNER00', 1), (-1.0, 'RP_Y', 1)))
    model.Equation(name='EQ_C01_U2', terms=((1.0, 'CORNER01', 2),
                   (-1.0, 'CORNER00', 2), (-1.0, 'RP_Y', 2)))
    model.Equation(name='EQ_C11_U1', terms=((1.0, 'CORNER11', 1),
                   (-1.0, 'CORNER00', 1), (-1.0, 'RP_X', 1), (-1.0, 'RP_Y', 1)))
    model.Equation(name='EQ_C11_U2', terms=((1.0, 'CORNER11', 2),
                   (-1.0, 'CORNER00', 2), (-1.0, 'RP_X', 2), (-1.0, 'RP_Y', 2)))

    # ---- remove rigid-body translation ----
    model.DisplacementBC(name='BC_PIN', createStepName='Initial',
                         region=a.sets['CORNER00'], u1=0.0, u2=0.0, ur3=UNSET)

    # ---- drive the macroscopic strain H through the reference points ----
    #      u_i(RIGHT)-u_i(LEFT)   = H_i1 * Lx_tot  -> RP_X u_i
    #      u_i(TOP)  -u_i(BOTTOM) = H_i2 * Ly_tot  -> RP_Y u_i
    model.DisplacementBC(name='BC_RPX', createStepName=STEP,
                         region=a.sets['RP_X'],
                         u1=H[0][0] * Lx_tot, u2=H[1][0] * Lx_tot, ur3=UNSET)
    model.DisplacementBC(name='BC_RPY', createStepName=STEP,
                         region=a.sets['RP_Y'],
                         u1=H[0][1] * Ly_tot, u2=H[1][1] * Ly_tot, ur3=UNSET)
else:
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
    try:
        return v.dataDouble
    except:
        return v.data


S_COMPS = ['S11', 'S22', 'S33', 'S12']
E_COMPS = ['E11', 'E22', 'E33', 'E12']

uSub = frame.fieldOutputs['U'].getSubset(region=odbInst, position=NODAL)
uDict = {}
for v in uSub.values:
    uDict[v.nodeLabel] = getData(v)


def nodalAveragedTensor(fieldName, wantedComps):
    fld = frame.fieldOutputs[fieldName].getSubset(region=odbInst,
                                                  position=ELEMENT_NODAL)
    labels = list(fld.componentLabels)
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


if EXPORT_NODAL_CSV:
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
else:
    print('EXPORT_NODAL_CSV = False: skipped per-node CSV export.')

# ------------------------------------------------------------
# VOLUME AVERAGES + FIRST MOMENTS over the whole domain
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


def volumeAverageAndMoment(fieldName, wantedComps):
    """IVOL-weighted volume average AND first spatial moment (against
    the domain-centred coordinates Xc, Yc) of an integration-point
    field, computed in one pass using the COORD/IVOL outputs already
    requested. Returns (avg, mom1, mom2, Vtot):

        avg[i]  = (1/V) * sum_ip  d_i * dV                 0th moment
        mom1[i] = (1/V) * sum_ip  d_i * (X-Xc) * dV         1st moment, X
        mom2[i] = (1/V) * sum_ip  d_i * (Y-Yc) * dV         1st moment, Y

    Used below for the double stress Q_ikl = <P_ik Xc_l>_V (P is
    approximated by the Cauchy stress S -- exact to leading order in
    strain0 for this small-strain analysis).
    """
    fld = frame.fieldOutputs[fieldName].getSubset(region=odbInst,
                                                  position=INTEGRATION_POINT)
    coordFld = frame.fieldOutputs['COORD'].getSubset(region=odbInst,
                                                      position=INTEGRATION_POINT)
    vol = frame.fieldOutputs['IVOL'].getSubset(region=odbInst,
                                               position=INTEGRATION_POINT)

    volMap = {}
    for v in vol.values:
        d = getData(v)
        volMap[(v.elementLabel, v.integrationPoint)] = d if isinstance(d, float) else d[0]

    coordMap = {}
    for v in coordFld.values:
        d = getData(v)
        coordMap[(v.elementLabel, v.integrationPoint)] = (d[0], d[1])

    labels = list(fld.componentLabels)
    idx = []
    for c in wantedComps:
        idx.append(labels.index(c) if c in labels else -1)

    nC = len(wantedComps)
    acc0 = [0.0] * nC
    acc1 = [0.0] * nC
    acc2 = [0.0] * nC
    Vtot = 0.0

    for v in fld.values:
        key = (v.elementLabel, v.integrationPoint)
        dV = volMap.get(key, 0.0)
        if dV == 0.0:
            continue
        xy = coordMap.get(key)
        if xy is None:
            continue
        xc1 = xy[0] - Xc
        xc2 = xy[1] - Yc
        d = getData(v)
        for q, ii in enumerate(idx):
            if ii >= 0:
                val = d[ii]
                acc0[q] += val * dV
                acc1[q] += val * xc1 * dV
                acc2[q] += val * xc2 * dV
        Vtot += dV

    if Vtot <= 0.0:
        print('WARNING: total IVOL is zero for field %s' % fieldName)
        nanList = [float('nan')] * nC
        return nanList, nanList, nanList, 0.0

    avg  = [a / Vtot for a in acc0]
    mom1 = [a / Vtot for a in acc1]
    mom2 = [a / Vtot for a in acc2]
    return avg, mom1, mom2, Vtot


Savg, Q1_S, Q2_S, Vtot = volumeAverageAndMoment('S', S_COMPS)
Eavg, _                = volumeAverage('E', E_COMPS)

# ------------------------------------------------------------
# Uavg -- AREA(VOLUME)-WEIGHTED, not a plain node-count mean.
#
# The old version:
#     nU = 0; u1sum=u2sum=0.0
#     for lab in uDict.keys(): u1sum += uDict[lab][0]; ... ; nU += 1
#     Uavg = [u1sum/nU, u2sum/nU]
#
# treats every FE NODE as equally weighted. Since the mesh is denser
# near the circular inclusion, this is a biased estimate of the true
# (1/Omega) * int u dV, and the bias tracks the inclusion position
# (zeta) through local mesh density -- exactly the drift seen in
# G111/G122/G211/G222 (the only G components whose formula subtracts
# Uavg) while G112/G121/G212/G221 (no Uavg term) stay exactly constant
# across translations.
#
# Fix: weight each element's nodal-averaged U by that element's own
# IVOL (already summed correctly for Savg/Eavg). For linear elements
# (CPE4/CPE3), the mean of an element's corner-node U values equals
# the exact area-average of the bilinear/linear interpolant over that
# element, so this reduces to a proper area-weighted average with no
# extra shape-function work needed.
# ------------------------------------------------------------
volField = frame.fieldOutputs['IVOL'].getSubset(region=odbInst,
                                                 position=INTEGRATION_POINT)
elemVol = {}
for v in volField.values:
    d = getData(v)
    dv = d if isinstance(d, float) else d[0]
    elemVol[v.elementLabel] = elemVol.get(v.elementLabel, 0.0) + dv

u1acc = 0.0
u2acc = 0.0
Vtot_u = 0.0
for elem in odbInst.elements:
    eVol = elemVol.get(elem.label, 0.0)
    if eVol <= 0.0:
        continue
    labs = elem.connectivity
    n = len(labs)
    u1sum = u2sum = 0.0
    nFound = 0
    for lab in labs:
        if lab in uDict:
            u1sum += uDict[lab][0]
            u2sum += uDict[lab][1]
            nFound += 1
    if nFound == 0:
        continue
    u1acc  += eVol * (u1sum / nFound)
    u2acc  += eVol * (u2sum / nFound)
    Vtot_u += eVol

Uavg = [u1acc / Vtot_u, u2acc / Vtot_u] if Vtot_u > 0.0 else [float('nan'), float('nan')]

print('Volume  = %.6e' % Vtot)
print('<S>     =', Savg)
print('<E>     =', Eavg)
print('<U>     =', Uavg)

# ------------------------------------------------------------
# DOUBLE STRESS  Q_ikl = <P_ik Xc_l>_V   (i, k, l in {1, 2})
#
# P (first Piola stress) is approximated by S (Cauchy stress): exact
# to leading order in strain0 for this small-strain analysis (P and S
# differ only by O(strain), ~1% here). S_COMPS = ['S11','S22','S33',
# 'S12'], so Q1_S/Q2_S (moments against Xc1, Xc2) index the same way.
# S12 = S21 (stress is symmetric, an exact property, not an
# approximation), so Q121 = Q211 and Q122 = Q212 exactly.
# ------------------------------------------------------------
Q111 = Q1_S[0]; Q112 = Q2_S[0]          # from S11
Q221 = Q1_S[1]; Q222 = Q2_S[1]          # from S22
Q121 = Q1_S[3]; Q122 = Q2_S[3]          # from S12
Q211 = Q1_S[3]; Q212 = Q2_S[3]          # from S21 = S12

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
        try:
            Fx += v.dataDouble[0]
            Fy += v.dataDouble[1]
        except:
            Fx += v.data[0]
            Fy += v.data[1]
    rfRow.extend([Fx, Fy])

# ------------------------------------------------------------
# PBC diagnostics: applied H, and RP_X/RP_Y reaction forces (the
# macroscopic force conjugate to H under periodic BC). NaN for
# non-PBC ('mixed') load cases.
# ------------------------------------------------------------
if LC['kind'] == 'affine':
    Happl_11, Happl_12, Happl_21, Happl_22 = (LC['H'][0][0], LC['H'][0][1],
                                              LC['H'][1][0], LC['H'][1][1])
else:
    Happl_11 = Happl_12 = Happl_21 = Happl_22 = float('nan')

RF_RPX = [float('nan'), float('nan')]
RF_RPY = [float('nan'), float('nan')]
if 'RP_X' in odbA.nodeSets.keys() and 'RP_Y' in odbA.nodeSets.keys():
    for v in rfField.getSubset(region=odbA.nodeSets['RP_X'], position=NODAL).values:
        d = getData(v)
        RF_RPX = [d[0], d[1]]
    for v in rfField.getSubset(region=odbA.nodeSets['RP_Y'], position=NODAL).values:
        d = getData(v)
        RF_RPY = [d[0], d[1]]

# ============================================================
# CALCULATED Fbar_ik and G_ikl -- BOUNDARY DISPLACEMENT ONLY.
#
# KUBC prescribes u_i(X) = H_ij Xc_j + 1/2 G_ijk Xc_j Xc_k EXACTLY at
# every boundary node (Xc = X - X_center). Two independent ways to
# recover H, G from that boundary trace are computed and both written
# out, so a mismatch between them is itself a diagnostic:
#
#   METHOD A -- LEAST-SQUARES FIT  (primary, machine precision)
#   Fit all boundary nodal (X, Y, u1, u2) to the quadratic basis
#       [1, Xc1, Xc2, 1/2 Xc1^2, Xc1 Xc2, 1/2 Xc2^2]
#       coeffs = [c_i, H_i1, H_i2, G_i11, G_i12, G_i22]
#   Because the BC is exact at every boundary node, this recovers the
#   applied H, G to round-off. The fitted constant c_i carries no
#   physical meaning here (the ansatz has none) and should be ~0 --
#   reported as fitResidMax / c1_fit / c2_fit diagnostics.
#
#   METHOD B -- CORRECTED SURFACE INTEGRAL  (independent cross-check)
#   S_ikm = oint u_i Xc_m n_k dS  (trapezoidal edge quadrature).
#   Product rule + divergence theorem give
#       S_ikm/V = G_ikm * (L_m^2/12) + delta_km * A_i
#   where A_i is the volume average of the PRESCRIBED POLYNOMIAL --
#   itself boundary-determined in closed form (NOT the FE field's own
#   <u>_V, which would leak interior/microstructure/translation
#   dependence into exactly the delta_km terms -- G111/G122/G211/G222
#   -- while leaving the others untouched; that was the earlier bug):
#       A_i = ( S_i11 + S_i22 ) / (4 V)
#   hence
#       G_ikm = ( S_ikm/V - delta_km A_i ) / (L_m^2/12)
#   This carries O(h^2) trapezoidal edge-quadrature error (~1e-5..1e-6
#   on a typical mesh) and is kept purely as a cross-check on Method A.
#
# Q (double stress, computed above via volumeAverageAndMoment) has NO
# such boundary shortcut -- it is the genuinely emergent, micro-
# structure-dependent quantity this whole pipeline exists to measure.
# ============================================================

def solveGauss(A, b):
    """Solve A x = b by Gaussian elimination with partial pivoting.
    Pure Python, no numpy dependency. A is n x n (list of lists), b is
    length n. Returns list x."""
    n = len(b)
    M = [list(A[r]) + [b[r]] for r in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1.0e-300:
            raise RuntimeError('Singular normal matrix at column %d -- '
                               'not enough distinct boundary nodes?' % col)
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for r in range(col + 1, n):
            f = M[r][col] / pv
            if f == 0.0:
                continue
            for c in range(col, n + 1):
                M[r][c] -= f * M[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        s = M[r][n]
        for c in range(r + 1, n):
            s -= M[r][c] * x[c]
        x[r] = s / M[r][r]
    return x


# ---- collect ALL boundary nodes (dict de-duplicates shared corners) ----
odbA = odb.rootAssembly
bndNodes = {}
for setName in ('LEFT', 'RIGHT', 'BOTTOM', 'TOP'):
    for arr in odbA.nodeSets[setName].nodes:
        for nd in arr:
            lab = nd.label
            if lab in uDict:
                u1, u2 = uDict[lab]
                bndNodes[lab] = (nd.coordinates[0], nd.coordinates[1], u1, u2)

nBnd = len(bndNodes)
if nBnd < 6:
    raise RuntimeError('Only %d boundary nodes found -- cannot fit the '
                       'quadratic ansatz (need >= 6).' % nBnd)

# ---------------- METHOD A: least-squares fit ----------------
nB = 6
ATA = [[0.0] * nB for _ in range(nB)]
ATb = [[0.0] * nB, [0.0] * nB]

for (x, y, u1, u2) in bndNodes.values():
    p = x - Xc
    q = y - Yc
    phi = [1.0, p, q, 0.5 * p * p, p * q, 0.5 * q * q]
    for a in range(nB):
        for c in range(nB):
            ATA[a][c] += phi[a] * phi[c]
        ATb[0][a] += phi[a] * u1
        ATb[1][a] += phi[a] * u2

fit = [solveGauss(ATA, ATb[0]), solveGauss(ATA, ATb[1])]

c1_fit, c2_fit = fit[0][0], fit[1][0]

Fbar_11 = 1.0 + fit[0][1]
Fbar_12 = fit[0][2]
Fbar_21 = fit[1][1]
Fbar_22 = 1.0 + fit[1][2]

G111, G112, G122 = fit[0][3], fit[0][4], fit[0][5]
G121 = G112
G211, G212, G222 = fit[1][3], fit[1][4], fit[1][5]
G221 = G212

fitResidMax = 0.0
for (x, y, u1, u2) in bndNodes.values():
    p = x - Xc
    q = y - Yc
    phi = [1.0, p, q, 0.5 * p * p, p * q, 0.5 * q * q]
    for i, uu in enumerate((u1, u2)):
        pred = 0.0
        for a in range(nB):
            pred += phi[a] * fit[i][a]
        d = abs(pred - uu)
        if d > fitResidMax:
            fitResidMax = d

# ---------------- METHOD B: corrected surface integral (cross-check) ----
def sideNodes(nm):
    pts = []
    for arr in odbA.nodeSets[nm].nodes:
        for nd in arr:
            pts.append((nd.label, nd.coordinates[0], nd.coordinates[1]))
    return pts


def edgeIntegrals(pts, axisIdx, uDictLocal):
    """Trapezoidal line integrals of u1, u2, u1*Xc1, u1*Xc2, u2*Xc1,
    u2*Xc2 along one boundary edge. axisIdx = 0 for a horizontal edge
    (BOTTOM/TOP, varies in X), 1 for a vertical edge (LEFT/RIGHT,
    varies in Y)."""
    pts = sorted(pts, key=lambda p: p[axisIdx + 1])
    n = len(pts)
    Iu1 = Iu2 = Iu1x1 = Iu1x2 = Iu2x1 = Iu2x2 = 0.0
    if n < 2:
        return (Iu1, Iu2, Iu1x1, Iu1x2, Iu2x1, Iu2x2)

    def sample(p):
        lab, x, y = p
        u1, u2 = uDictLocal.get(lab, (0.0, 0.0))
        xc1 = x - Xc
        xc2 = y - Yc
        return u1, u2, u1 * xc1, u1 * xc2, u2 * xc1, u2 * xc2

    prev = sample(pts[0])
    for i in range(1, n):
        ds = pts[i][axisIdx + 1] - pts[i - 1][axisIdx + 1]
        cur = sample(pts[i])
        Iu1   += 0.5 * (prev[0] + cur[0]) * ds
        Iu2   += 0.5 * (prev[1] + cur[1]) * ds
        Iu1x1 += 0.5 * (prev[2] + cur[2]) * ds
        Iu1x2 += 0.5 * (prev[3] + cur[3]) * ds
        Iu2x1 += 0.5 * (prev[4] + cur[4]) * ds
        Iu2x2 += 0.5 * (prev[5] + cur[5]) * ds
        prev = cur
    return (Iu1, Iu2, Iu1x1, Iu1x2, Iu2x1, Iu2x2)


L_Iu1, L_Iu2, L_Iu1x1, L_Iu1x2, L_Iu2x1, L_Iu2x2 = edgeIntegrals(sideNodes('LEFT'),   1, uDict)
R_Iu1, R_Iu2, R_Iu1x1, R_Iu1x2, R_Iu2x1, R_Iu2x2 = edgeIntegrals(sideNodes('RIGHT'),  1, uDict)
B_Iu1, B_Iu2, B_Iu1x1, B_Iu1x2, B_Iu2x1, B_Iu2x2 = edgeIntegrals(sideNodes('BOTTOM'), 0, uDict)
T_Iu1, T_Iu2, T_Iu1x1, T_Iu1x2, T_Iu2x1, T_Iu2x2 = edgeIntegrals(sideNodes('TOP'),    0, uDict)

Fbar_11_int = 1.0 + (-L_Iu1 + R_Iu1) / Vtot
Fbar_21_int = 0.0 + (-L_Iu2 + R_Iu2) / Vtot
Fbar_12_int = 0.0 + (-B_Iu1 + T_Iu1) / Vtot
Fbar_22_int = 1.0 + (-B_Iu2 + T_Iu2) / Vtot

varX = Lx_tot ** 2 / 12.0
varY = Ly_tot ** 2 / 12.0

S111 = (-L_Iu1x1 + R_Iu1x1) / Vtot
S112 = (-L_Iu1x2 + R_Iu1x2) / Vtot
S121 = (-B_Iu1x1 + T_Iu1x1) / Vtot
S122 = (-B_Iu1x2 + T_Iu1x2) / Vtot
S211 = (-L_Iu2x1 + R_Iu2x1) / Vtot
S212 = (-L_Iu2x2 + R_Iu2x2) / Vtot
S221 = (-B_Iu2x1 + T_Iu2x1) / Vtot
S222 = (-B_Iu2x2 + T_Iu2x2) / Vtot

# closed-form mean of the PRESCRIBED polynomial -- boundary data only,
# NOT the FE field's own <u>_V (that was the earlier bug).
A1 = 0.25 * (S111 + S122)
A2 = 0.25 * (S211 + S222)

G111_int = (S111 - A1) / varX
G112_int = (S112 - 0.0) / varY
G121_int = (S121 - 0.0) / varX
G122_int = (S122 - A1) / varY
G211_int = (S211 - A2) / varX
G212_int = (S212 - 0.0) / varY
G221_int = (S221 - 0.0) / varX
G222_int = (S222 - A2) / varY

print('Fbar (LSQ fit)  = [[%.10f, %.10f], [%.10f, %.10f]]'
      % (Fbar_11, Fbar_12, Fbar_21, Fbar_22))
print('Fbar (integral) = [[%.10f, %.10f], [%.10f, %.10f]]'
      % (Fbar_11_int, Fbar_12_int, Fbar_21_int, Fbar_22_int))
print('fitResidMax = %.3e   c1_fit = %.3e   c2_fit = %.3e'
      % (fitResidMax, c1_fit, c2_fit))

# ============================================================
# SUMMARY CSV
# ============================================================
sumPath = os.path.join(caseEnergyDir, '%s_Summary.csv' % JOB)
fsum = open(sumPath, 'w')
w = csv.writer(fsum)
w.writerow(['LoadCase', 'NX', 'NY', 'zeta1', 'zeta2', 'strain0', 'StepTime', 'ALLSE',
            'RF_LEFT_x', 'RF_LEFT_y', 'RF_RIGHT_x', 'RF_RIGHT_y',
            'RF_BOTTOM_x', 'RF_BOTTOM_y', 'RF_TOP_x', 'RF_TOP_y',
            'Volume', 'nBndNodes',
            'Savg_11', 'Savg_22', 'Savg_33', 'Savg_12',
            'Eavg_11', 'Eavg_22', 'Eavg_33', 'Eavg_12',
            'Uavg_1', 'Uavg_2',
            'Fbar_11', 'Fbar_12', 'Fbar_21', 'Fbar_22',
            'G111', 'G112', 'G121', 'G122', 'G211', 'G212', 'G221', 'G222',
            'Fbar_11_int', 'Fbar_12_int', 'Fbar_21_int', 'Fbar_22_int',
            'G111_int', 'G112_int', 'G121_int', 'G122_int',
            'G211_int', 'G212_int', 'G221_int', 'G222_int',
            'Q111', 'Q112', 'Q121', 'Q122', 'Q211', 'Q212', 'Q221', 'Q222',
            'fitResidMax', 'c1_fit', 'c2_fit',
            'Happl_11', 'Happl_12', 'Happl_21', 'Happl_22',
            'RF_RPX_1', 'RF_RPX_2', 'RF_RPY_1', 'RF_RPY_2'])
w.writerow([loadCaseName, NX, NY, zeta1, zeta2, strain0, frame.frameValue, allseVal]
           + rfRow + [Vtot, nBnd] + Savg + Eavg + Uavg
           + [Fbar_11, Fbar_12, Fbar_21, Fbar_22]
           + [G111, G112, G121, G122, G211, G212, G221, G222]
           + [Fbar_11_int, Fbar_12_int, Fbar_21_int, Fbar_22_int]
           + [G111_int, G112_int, G121_int, G122_int,
              G211_int, G212_int, G221_int, G222_int]
           + [Q111, Q112, Q121, Q122, Q211, Q212, Q221, Q222]
           + [fitResidMax, c1_fit, c2_fit]
           + [Happl_11, Happl_12, Happl_21, Happl_22]
           + RF_RPX + RF_RPY)
fsum.close()
print('Saved summary file:', sumPath)

odb.close()
print('Finished job:', JOB)
