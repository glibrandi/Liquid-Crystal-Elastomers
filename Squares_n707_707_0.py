########################################################
# Author: Gabriele Librandi, PhD 
# Materials Science & Mechanical Engineering 
# Harvard University
# written in Oct 2018
########################################################
# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy as np

##-----------------------------------------------
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
##-----------------------------------------------
Mdb()

## Geometry
L = 75.0 #length
h = 100.0 #height 
t = 15.0 #thickness
## Seeding size
S = 3.76 #t/4
TOL = S/5
TOLL = 0.1E-004
## Material
RHO = 1 #[Kg/m^3]
N_modes = 20
# E = 1.125e6 #[Pa]
# G = 0.375e6 #[Pa]
# K = 971e6 #[Pa]
# C10 = G/2 #[Pa]
# D1 = 2/K #[1/Pa]
# PR = 0.4998
# alpha = 80e-6
 
###########################
## Unit cell geometry #####
###########################
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-L, t/2), 
    point2=(L, -t/2))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseSolidExtrude(depth=h, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-2', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-3', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-4', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-2', ), 
    vector=(0.0, L, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-3', ), 
    vector=(0.0, L/2, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(0.0, L/2, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=90.0, axisDirection=(0.0, 0.0, 
    h), axisPoint=(-L, -t/2, 0.0), instanceList=('Part-1-3', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-3', ), 
    vector=((L*0.5)+t, -(L-t)*0.5, 0.0))
mdb.models['Model-1'].rootAssembly.rotate(angle=-90.0, axisDirection=(0.0, 0.0, 
    h), axisPoint=(+L, -t/2, 0.0), instanceList=('Part-1-4', ))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(-((L*0.5)+t), -(L-t)*0.5, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-3', ), 
    vector=((L-t)*0.5, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(-(L-t)*0.5, 0.0, 0.0))

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-4']), name='Part-2', 
    originalInstances=SUPPRESS)

partname = 'Part-2'

###########################
## Partitioning ###########
###########################

mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-45.0, 67.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-30.0, 
    82.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((30.0, 82.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((45.0, 
    67.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-45.0, -7.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-30.0, 
    7.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-45.0, 82.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-30.0, 
    67.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((30.0, 67.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((45.0, 
    82.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((30.0, 7.5, 100.0), )
    , point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((45.0, -7.5, 
    100.0), ))
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((-30.0, 27.5, 
    33.333333), )), datumPlane=mdb.models['Model-1'].parts['Part-2'].datums[2])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((10.0, 67.5, 
    66.666667), ), ((-55.0, 67.5, 66.666667), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[3])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((-30.0, 92.5, 
    33.333333), ), ((-55.0, 7.5, 33.333333), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[5])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((-45.0, -17.5, 
    33.333333), ), ((-55.0, 7.5, 33.333333), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[4])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((45.0, 92.5, 
    33.333333), ), ((55.0, -7.5, 33.333333), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[6])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((45.0, -17.5, 
    66.666667), ), ((55.0, 7.5, 66.666667), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[7])

mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-45.0, 112.5, 
    100.0), ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((
    -30.0, 112.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((30.0, 112.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((45.0, 
    112.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-75.0, 82.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-75.0, 
    67.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].DatumPlaneByTwoPoint(point1=
    mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-75.0, 7.5, 100.0), 
    ), point2=mdb.models['Model-1'].parts['Part-2'].vertices.findAt((-75.0, 
    -7.5, 100.0), ))
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((10.0, -7.5, 
    66.666667), ), ((10.0, 82.5, 33.333333), ), ((-55.0, 7.5, 33.333333), ), ((
    -55.0, 82.5, 33.333333), ), ((-30.0, -17.5, 66.666667), ), ((-45.0, 92.5, 
    66.666667), ), ((-45.0, 27.5, 66.666667), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[14])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((30.0, 47.5, 
    33.333333), ), ((10.0, -7.5, 66.666667), ), ((55.0, 82.5, 66.666667), ), ((
    10.0, 82.5, 33.333333), ), ((45.0, -17.5, 66.666667), ), ((45.0, 92.5, 
    33.333333), ), ((55.0, 7.5, 66.666667), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[15])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((35.0, 90.0, 0.0), ), (
    (40.0, 47.5, 0.0), ), ((-40.0, 47.5, 100.0), ), ((-35.0, 90.0, 100.0), ), (
    (35.0, 47.5, 0.0), ), ((55.0, 67.5, 33.333333), ), ((10.0, 82.5, 
    33.333333), ), ((-55.0, 82.5, 33.333333), ), ((45.0, 92.5, 33.333333), ), (
    (-45.0, 92.5, 66.666667), ), ((-35.0, 47.5, 100.0), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[16])
mdb.models['Model-1'].parts['Part-2'].PartitionCellByDatumPlane(cells=
    mdb.models['Model-1'].parts['Part-2'].cells.findAt(((35.0, -15.0, 100.0), 
    ), ((40.0, 47.5, 0.0), ), ((-40.0, 47.5, 100.0), ), ((-40.0, -15.0, 100.0), 
    ), ((35.0, 47.5, 0.0), ), ((10.0, -7.5, 66.666667), ), ((-55.0, -7.5, 
    66.666667), ), ((45.0, -17.5, 66.666667), ), ((-30.0, -17.5, 66.666667), ), 
    ((55.0, -7.5, 33.333333), ), ((-35.0, 47.5, 100.0), ), ), datumPlane=
    mdb.models['Model-1'].parts['Part-2'].datums[17])

###########################
## Mesh ###################
###########################

mdb.models['Model-1'].parts[partname].seedPart(deviationFactor=0.01, 
    minSizeFactor=0.01, size=S)
#mdb.models['Model-1'].parts[partName].setMeshControls(algorithm=ADVANCING_FRONT, regions=mdb.models['Model-1'].parts[partName].cells)

mdb.models['Model-1'].parts[partname].setMeshControls(algorithm=ADVANCING_FRONT, 
    regions=mdb.models['Model-1'].parts[partname].cells, technique=SWEEP)
mdb.models['Model-1'].parts[partname].setElementType(elemTypes=(ElemType(
    elemCode=C3D8, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
    distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts[partname].cells, ))

mdb.models['Model-1'].parts[partname].generateMesh()

###########################
## Set definitions ########
###########################

mdb.models['Model-1'].parts[partname].Set(name='bottom', nodes= 
    mdb.models['Model-1'].parts[partname].nodes.getByBoundingBox( -10000, -10000, -TOL, +10000, +10000, +TOL))


mdb.models['Model-1'].parts[partname].Set(name='left', nodes= 
    mdb.models['Model-1'].parts[partname].nodes.getByBoundingBox( -L-TOL, -1000, 0.0+TOL, -L+TOL, +1000, +1000))

mdb.models['Model-1'].parts[partname].Set(name='right', nodes= 
    mdb.models['Model-1'].parts[partname].nodes.getByBoundingBox( +L-TOL, -1000, 0.0+TOL, +L+TOL, +1000, +1000))

mdb.models['Model-1'].parts[partname].Set(name='down', nodes= 
    mdb.models['Model-1'].parts[partname].nodes.getByBoundingBox( -1000, -L/2-TOL, 0.0+TOL, +1000, -L/2+TOL, +1000))

mdb.models['Model-1'].parts[partname].Set(name='up', nodes= 
    mdb.models['Model-1'].parts[partname].nodes.getByBoundingBox( -1000, +(L+L/2)-TOL, 0.0+TOL, +1000, +(L+L/2)+TOL, +1000))

mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_1', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[24:25])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_2', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[31:32])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_3', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[9:10])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_4', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[4:5])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_1_UP', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1242:1243])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_1_BOTTOM', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1184:1185])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_1_LEFT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[620:621])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_1_RIGHT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[754:755])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_2_UP', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[974:975])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_2_RIGHT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[888:889])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_2_BOTTOM', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1136:1137])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_2_LEFT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[766:767])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_3_UP', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1172:1173])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_3_RIGHT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[472:473])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_3_BOTTOM', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1326:1327])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_3_LEFT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[232:233])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_4_UP', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[1148:1149])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_4_RIGHT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[212:213])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_4_BOTTOM', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[992:993])
mdb.models['Model-1'].parts['Part-2'].Set(name='vertex_4_LEFT', nodes=
    mdb.models['Model-1'].parts['Part-2'].nodes[460:461])

###########################
## Material ###############
###########################

mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Density(table=((RHO, ), ))
mdb.models['Model-1'].materials['Material-1'].UserMaterial(mechanicalConstants=
    (1.0, 0.49, cos(45*pi/180), sin(45*pi/180), 0.0))
mdb.models['Model-1'].HomogeneousSolidSection(material='Material-1', name=
    'Section-1', thickness=None)

###########################
## Coupling nodes for PBC #
###########################

# SETS OF PERIODIC NODE PAIRS
# UP AND DOWN 
UpDown_Index = []
for i in mdb.models['Model-1'].parts['Part-2'].sets['up'].nodes:
    TopCoordinate = i.coordinates
    for j in mdb.models['Model-1'].parts['Part-2'].sets['down'].nodes:
        DownCoordinate = j.coordinates
        if (fabs(TopCoordinate[0] - DownCoordinate[0]) + fabs(TopCoordinate[2] - DownCoordinate[2]) < TOLL):
            mdb.models['Model-1'].parts['Part-2'].SetFromNodeLabels(name='Up_Node_Pair_' + str(i.label), nodeLabels=(i.label,))
            mdb.models['Model-1'].parts['Part-2'].SetFromNodeLabels(name='Down_Node_Pair_' + str(i.label), nodeLabels=(j.label,))
            UpDown_Index.append(i)
            break

# LEFT AND RIGHT
RightLeft_Index = []
for i in mdb.models['Model-1'].parts['Part-2'].sets['right'].nodes:
    RightCoordinate = i.coordinates
    for j in mdb.models['Model-1'].parts['Part-2'].sets['left'].nodes:
        LeftCoordinate = j.coordinates
        if (fabs(RightCoordinate[1] - LeftCoordinate[1]) + fabs(RightCoordinate[2] - LeftCoordinate[2]) < TOLL):
            mdb.models['Model-1'].parts['Part-2'].SetFromNodeLabels(name='Right_Node_Pair_' + str(i.label), nodeLabels=(i.label,))
            mdb.models['Model-1'].parts['Part-2'].SetFromNodeLabels(name='Left_Node_Pair_' + str(i.label), nodeLabels=(j.label,))
            RightLeft_Index.append(i)
            break

Part_VP_A = mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Virtual_point_A', type=DEFORMABLE_BODY)
Part_VP_A.ReferencePoint(point=(0.0, 0.0, 0.0))

Part_VP_B = mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Virtual_point_B', type=DEFORMABLE_BODY)
Part_VP_B.ReferencePoint(point=(0.0, 0.0, 0.0))

# VIRTUAL POINTS INSTANTS 
InstVP_X = mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='VirtPointInst_X', part=Part_VP_A)
InstVP_Y = mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='VirtPointInst_Y', part=Part_VP_B)

# CREATES SETS FOR VIRTUAL NODES
Part_VP_A.Set(name='V_point_A', referencePoints=(Part_VP_A.referencePoints[1], ))
Part_VP_B.Set(name='V_point_B', referencePoints=(Part_VP_B.referencePoints[1], ))

# RIGHT AND LEFT EDGES
for i in RightLeft_Index:
    
    # COEFFICIENTS PREPARATION
    InDependCoord=mdb.models['Model-1'].parts['Part-2'].sets['Left_Node_Pair_' + str(i.label)].nodes[0].coordinates
    DependCoord=mdb.models['Model-1'].parts['Part-2'].sets['Right_Node_Pair_' + str(i.label)].nodes[0].coordinates
    
    coeff1=-fabs(DependCoord[0]-InDependCoord[0])
    coeff2=-fabs(DependCoord[1]-InDependCoord[1])
    coeff3=0
 
    # X-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='LR_Const_at_X_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Right_Node_Pair_' + str(i.label), 1), 
        (-1.0, 'Part-2-1.Left_Node_Pair_' + str(i.label), 1), 
        (coeff1, 'VirtPointInst_X.V_point_A', 1), 
        (coeff2, 'VirtPointInst_Y.V_point_B', 1)))
        
    # Y-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='LR_Const_at_Y_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Right_Node_Pair_' + str(i.label), 2), 
        (-1.0, 'Part-2-1.Left_Node_Pair_' + str(i.label), 2), 
        (coeff1, 'VirtPointInst_X.V_point_A', 2), 
        (coeff2, 'VirtPointInst_Y.V_point_B', 2)))

    # Z-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='LR_Const_at_Z_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Right_Node_Pair_' + str(i.label), 3), 
        (-1.0, 'Part-2-1.Left_Node_Pair_' + str(i.label), 3))) 

# UP AND DOWN EDGES
for i in UpDown_Index:
    
    # COEFFICIENTS PREPARATION
    InDependCoord=mdb.models['Model-1'].parts['Part-2'].sets['Down_Node_Pair_' + str(i.label)].nodes[0].coordinates
    DependCoord=mdb.models['Model-1'].parts['Part-2'].sets['Up_Node_Pair_' + str(i.label)].nodes[0].coordinates
    
    coeff1=-fabs(DependCoord[0]-InDependCoord[0])
    coeff2=-fabs(DependCoord[1]-InDependCoord[1])
    coeff3=0
   
    # X-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_X_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Up_Node_Pair_' + str(i.label), 1), 
        (-1.0, 'Part-2-1.Down_Node_Pair_' + str(i.label), 1), 
        (coeff1, 'VirtPointInst_X.V_point_A', 1), 
        (coeff2, 'VirtPointInst_Y.V_point_B', 1)))
        
    # Y-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_Y_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Up_Node_Pair_' + str(i.label), 2), 
        (-1.0, 'Part-2-1.Down_Node_Pair_' + str(i.label), 2), 
        (coeff1, 'VirtPointInst_X.V_point_A', 2), 
        (coeff2, 'VirtPointInst_Y.V_point_B', 2)))

    # Z-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_Z_' + str(i.label), terms=(
        ( 1.0, 'Part-2-1.Up_Node_Pair_' + str(i.label), 3), 
        (-1.0, 'Part-2-1.Down_Node_Pair_' + str(i.label), 3))) 

###########################
## Customize analysis #####
###########################

mdb.models['Model-1'].StaticStep(initialInc=0.1, maxInc=1.0, maxNumInc=10000000000000, 
    minInc=10e-40, name='temp_applied', nlgeom=ON, previous='Initial')

###########################
## Section assignments ####
###########################

mdb.models['Model-1'].parts[partname].Set(cells=
    mdb.models['Model-1'].parts[partname].cells, name='Set-4')

mdb.models['Model-1'].parts[partname].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts[partname].sets['Set-4'], sectionName=
    'Section-1', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

# X VIRTUAL POINT 
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'temp_applied', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='VP_along_x', region=InstVP_X.sets['V_point_A'], 
    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
# Y VIRTUAL POINT 
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'temp_applied', distributionType=UNIFORM, fieldName='', fixed=OFF,
    localCsys=None, name='VP_along_y', region=InstVP_Y.sets['V_point_B'], 
    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# FIX BOTTOM ### ===================================
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'temp_applied', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='BC-1', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-2-1'].sets['bottom'], 
    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Squares_n707_707_0', nodalOutputPrecision=SINGLE, 
    numCpus=2, numDomains=2, numGPUs=0, queue=None, scratch='', type=ANALYSIS, 
    userSubroutine='', waitHours=0, waitMinutes=0)

# mdb.jobs['Square_UNIT_CELL_1_LCE_n_001'].submit(consistencyChecking=OFF)
#     # mdb.jobs[CurJobName].submit(consistencyChecking=OFF, datacheckJob=True)
# mdb.jobs['Square_UNIT_CELL_1_LCE_n_001'].waitForCompletion()



