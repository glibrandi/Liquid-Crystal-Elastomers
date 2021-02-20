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
S = t/2
TOL = S/5
TOLL = 0.1E-004
N_modes = 60
RHO = 1
## Material
# RHO = 1070 #[Kg/m^3]
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
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(L, 0.0), 
    point2=(0, L))
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(L-t/2, +t/2), 
    point2=(+t/2, +L-t/2))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-1'].BaseSolidExtrude(depth=h, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-2', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-2', ), 
    vector=(L, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-3', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-3', ), 
    vector=(0.0, -L, 0.0))
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-1-4', 
    part=mdb.models['Model-1'].parts['Part-1'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-1-4', ), 
    vector=(L, -L, 0.0))

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-1-4']), name='Part-2', 
    originalInstances=SUPPRESS)

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-1', 
    part=mdb.models['Model-1'].parts['Part-2'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-2-2', 
    part=mdb.models['Model-1'].parts['Part-2'])
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-2-2', ), 
    vector=(0.0, -2*L, 0.0))

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-2-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-2-2']), name='Part-3', 
    originalInstances=SUPPRESS)

for i in range(1,9):

    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-3-'+str(i)+'', 
        part=mdb.models['Model-1'].parts['Part-3'])
    mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-3-'+str(i)+'', ), 
        vector=(+(i-1)*(2*L), 0.0, 0.0))

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-3-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-3-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-3-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-3-4'],
    mdb.models['Model-1'].rootAssembly.instances['Part-3-5'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-3-6'],
    mdb.models['Model-1'].rootAssembly.instances['Part-3-7'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-3-8']), name='Part-5', 
    originalInstances=SUPPRESS)

for i in range(1,5):
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-5-'+str(i)+'', 
        part=mdb.models['Model-1'].parts['Part-5'])

L_tot = 1200.0
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-5-3', ), 
    vector=(L_tot, 0.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-5-2', ), 
    vector=(0.0, 300.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-5-4', ), 
    vector=(L_tot, 300.0, 0.0))

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-5-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-5-2'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-5-3'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-5-4']), name='Part-8', 
    originalInstances=SUPPRESS)

mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(67.5, -75.0), 
    point2=(1582.5, 75.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((67.5, 0.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((825.0, 75.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((1582.5, 0.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((825.0, -75.0))
mdb.models['Model-1'].sketches['__profile__'].offset(distance=2000.0, 
    objectList=(mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((
    67.5, 0.0), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((825.0, 
    75.0), ), mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((
    1582.5, 0.0), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((825.0, 
    -75.0), )), side=RIGHT)
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-6', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-6'].BaseSolidExtrude(depth=100.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-6-1', 
    part=mdb.models['Model-1'].parts['Part-6'])

mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
    mdb.models['Model-1'].rootAssembly.instances['Part-6-1'], ), 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Part-8-1'], 
    name='Part-12', originalInstances=DELETE)

mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-12-1', 
    part=mdb.models['Model-1'].parts['Part-12'])
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-12-2', 
    part=mdb.models['Model-1'].parts['Part-12'])


mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-12-1', ), 
    vector=(-67.5, 75.0, 0.0))
mdb.models['Model-1'].rootAssembly.translate(instanceList=('Part-12-2', ), 
    vector=(7.5, 75.0, 0.0))
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanMerge(domain=GEOMETRY, 
    instances=(mdb.models['Model-1'].rootAssembly.instances['Part-12-1'], 
    mdb.models['Model-1'].rootAssembly.instances['Part-12-2']), name='Part-9', 
    originalInstances=SUPPRESS)
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(33.75, -22.5))
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(1515.0, 150.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 75.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((757.5, 150.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((1515.0, 75.0))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((757.5, 0.0))
mdb.models['Model-1'].sketches['__profile__'].offset(distance=2000.0, 
    objectList=(mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((
    0.0, 75.0), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((757.5, 
    150.0), ), mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((
    1515.0, 75.0), ), 
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((757.5, 0.0), 
    )), side=RIGHT)
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Part-10', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Part-10'].BaseSolidExtrude(depth=100.0, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Part-10-1', 
    part=mdb.models['Model-1'].parts['Part-10'])
mdb.models['Model-1'].rootAssembly.InstanceFromBooleanCut(cuttingInstances=(
    mdb.models['Model-1'].rootAssembly.instances['Part-10-1'], ), 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Part-9-1'], 
    name='Part-4', originalInstances=DELETE)

###########################
## Materials ##############
###########################
mdb.models['Model-1'].Material(name='Subroutine_1')
mdb.models['Model-1'].materials['Subroutine_1'].Density(table=((RHO, ), ))
mdb.models['Model-1'].materials['Subroutine_1'].Depvar(n=5)
mdb.models['Model-1'].materials['Subroutine_1'].UserMaterial(mechanicalConstants=
    (1.0, 0.49, 0.0, 1.515E+03, 0.0))
mdb.models['Model-1'].HomogeneousSolidSection(material='Subroutine_1', name=
    'Cube_section_1', thickness=None)

###########################
## Mesh ###################
###########################

mdb.models['Model-1'].parts['Part-4'].seedPart(deviationFactor=0.01, 
    minSizeFactor=0.01, size=S)

mdb.models['Model-1'].parts['Part-4'].setMeshControls(algorithm=MEDIAL_AXIS, 
    regions=mdb.models['Model-1'].parts['Part-4'].cells.findAt(((102.5, 7.5, 
    66.666667), )))
mdb.models['Model-1'].parts['Part-4'].setElementType(elemTypes=(ElemType(
    elemCode=C3D8, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
    distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
    ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Part-4'].cells.findAt(((102.5, 7.5, 
    66.666667), )), ))

mdb.models['Model-1'].parts['Part-4'].generateMesh()

###########################
## Set definitions ########
###########################

mdb.models['Model-1'].parts['Part-4'].Set(name='bottom', nodes= 
    mdb.models['Model-1'].parts['Part-4'].nodes.getByBoundingBox( -1000000, -1000000, -TOL, +1000000, +1000000, +TOL))
mdb.models['Model-1'].parts['Part-4'].Set(name='all', nodes= 
    mdb.models['Model-1'].parts['Part-4'].nodes.getByBoundingBox( -1000000, -1000000, -1000000, +1000000, +1000000, +1000000))
mdb.models['Model-1'].parts['Part-4'].Set(name='central_nodes', nodes= 
    mdb.models['Model-1'].parts['Part-4'].nodes.getByBoundingBox( -1000000, 37.5-TOL, 100-TOL, +1000000, 37.5+TOL,  100+TOL))
mdb.models['Model-1'].parts['Part-4'].Set(name='up', nodes= 
    mdb.models['Model-1'].parts['Part-4'].nodes.getByBoundingBox( -1000000, 150-TOL, -100, +1000000, 150+TOL,  100+TOL))
mdb.models['Model-1'].parts['Part-4'].Set(name='down', nodes= 
    mdb.models['Model-1'].parts['Part-4'].nodes.getByBoundingBox( -1000000, 0.0-TOL, -100, +1000000, 0.0+TOL,  100+TOL))

import numpy as np
import random

odb=openOdb(path='Squares_PBC_FINITE_SIZE_LCE_Nem_director_gradient_eigenfreq_20holes.odb')
pbpPartNodes=odb.rootAssembly.nodeSets[' ALL NODES']

#CREATE A MATRIX TO SAVE THE NEW COORDINATES OF ALL NODES
NewCoord=np.zeros((len(mdb.models['Model-1'].parts['Part-4'].nodes), 3))

#SELECT THE WANTED BUCKLING MODE AND IMPERFECTION WEIGHT
ImpFrames  = []
ImpWeights = []
for CImp in range(len(ImpFrames)):
    cframe = ImpFrames[CImp]
    firstFrame = odb.steps['Frequency_1'].frames[cframe]
    displacement = firstFrame.fieldOutputs['U']
    pbpDispField = displacement.getSubset(region=pbpPartNodes)
    pbpDisp = pbpDispField.values
# Imperfection Using Buckling Analysis Results
#---------------------------------------------------------------
    ind=0;
    IMP = ImpWeights[CImp]
#CALCULATE THE MODIFIED COORDINATES
    for i in mdb.models['Model-1'].parts['Part-4'].nodes:
        NewCoord[ind][0]=i.coordinates[0]+IMP*pbpDisp[ind].data[0]
        NewCoord[ind][1]=i.coordinates[1]+IMP*pbpDisp[ind].data[1]
        NewCoord[ind][2]=i.coordinates[2]+IMP*pbpDisp[ind].data[2]
        ind=ind+1

#SET THE NEW COORDINATES
    mdb.models['Model-1'].parts['Part-4'].editNode(
    nodes=mdb.models['Model-1'].parts['Part-4'].nodes,
    coordinates=NewCoord)

mdb.models['Model-1'].rootAssembly.regenerate()

###########################
## Coupling nodes for PBC #
###########################

# SETS OF PERIODIC NODE PAIRS
# UP AND DOWN 
UpDown_Index = []
for i in mdb.models['Model-1'].parts['Part-4'].sets['up'].nodes:
    TopCoordinate = i.coordinates
    for j in mdb.models['Model-1'].parts['Part-4'].sets['down'].nodes:
        DownCoordinate = j.coordinates
        if (fabs(TopCoordinate[0] - DownCoordinate[0]) + fabs(TopCoordinate[2] - DownCoordinate[2]) < TOLL):
            mdb.models['Model-1'].parts['Part-4'].SetFromNodeLabels(name='Up_Node_Pair_' + str(i.label), nodeLabels=(i.label,))
            mdb.models['Model-1'].parts['Part-4'].SetFromNodeLabels(name='Down_Node_Pair_' + str(i.label), nodeLabels=(j.label,))
            UpDown_Index.append(i)
            break


Part_VP_B = mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Virtual_point_B', type=DEFORMABLE_BODY)
Part_VP_B.ReferencePoint(point=(0.0, 0.0, 0.0))

# VIRTUAL POINTS INSTANTS 
InstVP_Y = mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='VirtPointInst_Y', part=Part_VP_B)

# CREATES SETS FOR VIRTUAL NODES
Part_VP_B.Set(name='V_point_B', referencePoints=(Part_VP_B.referencePoints[1], ))

# UP AND DOWN EDGES
for i in UpDown_Index:
    
    # COEFFICIENTS PREPARATION
    InDependCoord=mdb.models['Model-1'].parts['Part-4'].sets['Down_Node_Pair_' + str(i.label)].nodes[0].coordinates
    DependCoord=mdb.models['Model-1'].parts['Part-4'].sets['Up_Node_Pair_' + str(i.label)].nodes[0].coordinates
    
    coeff1=-fabs(DependCoord[0]-InDependCoord[0])
    coeff2=-fabs(DependCoord[1]-InDependCoord[1])
    coeff3=0
   
    # X-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_X_' + str(i.label), terms=(
        ( 1.0, 'Part-4-1.Up_Node_Pair_' + str(i.label), 1), 
        (-1.0, 'Part-4-1.Down_Node_Pair_' + str(i.label), 1),  
        (coeff2, 'VirtPointInst_Y.V_point_B', 1)))
        
    # Y-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_Y_' + str(i.label), terms=(
        ( 1.0, 'Part-4-1.Up_Node_Pair_' + str(i.label), 2), 
        (-1.0, 'Part-4-1.Down_Node_Pair_' + str(i.label), 2),  
        (coeff2, 'VirtPointInst_Y.V_point_B', 2)))

    # Z-COORDINATE OF DEPENDENT SET
    mdb.models['Model-1'].Equation(name='UD_Const_at_Z_' + str(i.label), terms=(
        ( 1.0, 'Part-4-1.Up_Node_Pair_' + str(i.label), 3), 
        (-1.0, 'Part-4-1.Down_Node_Pair_' + str(i.label), 3))) 


###########################
## Customize analysis #####
###########################

mdb.models['Model-1'].StaticStep(initialInc=0.000001, maxInc=0.001, maxNumInc=10000000000000000, 
    minInc=1e-40, name='temp_applied', nlgeom=ON, previous='Initial')

mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(numIntervals=
    100)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(timeMarks=
    OFF)

###########################
## Section assignments ####
###########################

mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].parts['Part-4'].Set(cells=
    mdb.models['Model-1'].parts['Part-4'].cells.findAt(((+TOL, +TOL, 
    h-TOL), )), name='Set-1')

mdb.models['Model-1'].parts['Part-4'].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=MIDDLE_SURFACE, region=
    mdb.models['Model-1'].parts['Part-4'].sets['Set-1'], sectionName=
    'Cube_section_1', thicknessAssignment=FROM_SECTION)

mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)

# Y VIRTUAL POINT 
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'temp_applied', distributionType=UNIFORM, fieldName='', fixed=OFF,
    localCsys=None, name='VP_along_y_1', region=InstVP_Y.sets['V_point_B'], 
    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# FIX BOTTOM ### ===================================
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName=
    'temp_applied', distributionType=UNIFORM, fieldName='', fixed=OFF, 
    localCsys=None, name='fix', region=
    mdb.models['Model-1'].rootAssembly.instances['Part-4-1'].sets['bottom']
    , u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)




